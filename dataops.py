from __future__ import print_function, division
import os
import sys
import logging
import __main__
import argparse
import tempfile
import shutil
import gzip
import string
import urllib
import atexit
import urllib
import csv
import datetime
from time import sleep
from Bio import SeqIO
from itertools import cycle
from glob import glob
from subprocess import Popen, PIPE
from HTSeq import FastqReader
from numpy import std
from collections import OrderedDict


log = logging.getLogger(__main__.__name__)

class Error(Exception):
	pass

class DataWrapper:
    """Object used to create a temporary folder where data and temporary files
    will be created.

    Provided files are usually in gzip format. This wrapper also ensures they
    are decompressed in place.

    Needs a ConfigParser object containing a data entry with information about
    each sample used.
    """

    def __init__(self, config):
        self.dir = tempfile.mkdtemp(prefix="genoqual_")
        atexit.register(lambda: self.cleanup)

        self.config = config
        self.data = []

    def cleanup(self):
        if os.path.isdir(self.dir):
            log.debug("Removing temporary folder and all its contents")
            shutil.rmtree(self.dir)

    def _get_files(self, prefix):
        """Get all files that match the given pattern in the input folder
        """
        files = glob(os.path.join(self.config.params["input_path"],
                                  "{0}*".format(prefix)))

        # Check that we have files with valid extensions
        valid = {".gz", ".gzip", ".fastq"}
        valid_files = []
        for file in files:
            if os.path.splitext(file)[-1] in valid:
                valid_files.append(file)
            else:
                log.warning("Unexpected file extension in file '%s'. "
                            "Ignoring file.", file)
        return tuple(sorted(valid_files))

    def prepare(self):
        log.debug("Making sure final output folder is clean")

        final_outdir = self.config.params["output_path"]
        if os.path.isdir(final_outdir):
            if not glob(os.path.join(final_outdir, "*")):
                shutil.rmtree(final_outdir)
            else:
                raise Error("Output folder '{0}' exists and is not empty. "
                            "Remove it before re-launching the analysis."
                            .format(final_outdir))

        log.debug("Looking for files matching identifiers in the config file")

        for entry in self.config.config:
            prefix = entry["prefix"]
            files = self._get_files(prefix)
            entry["original_files"] = files

            # Check that we have paired reads
            count = len(files)

            if not count:
                raise Error("No files found matching the prefix '{0}'"
                            .format(prefix))

            if count == 1:
                log.info("Assuming data as single-end on Sample ID '%s'.",
                         entry["Sample_ID"])

            elif count == 2 and self.config.params["meta"] == False:
                log.info("Assuming data as paired-end on Sample ID '%s'.",
                         entry["Sample_ID"])
            elif count == 4:
                log.info("Assuming data as NextSeq single-end on Sample ID '%s'.",
                         entry["Sample_ID"])
                self.config.params["nextseq"] = True
            elif count == 8:
                log.info("Assuming data as NextSeq paired-end on Sample ID '%s'.",
                         entry["Sample_ID"])
                self.config.params["nextseq"] = True

            elif self.config.params["meta"] == True:
                if count == 3:
                    log.info("All three files for metagenomics analysis were found for Sample ID '%s'.",
                             entry["Sample_ID"])
                elif count == 2 and '_I1_001.fastq' in ' '.join(files):
                    log.info("Two files for metagenomics analysis were found for Sample ID '%s'. Assuming R1 and I1.",
                             entry["Sample_ID"])
                else:
                    raise Error("Unexpected number of files for metagenomics analysis matching the prefix "
                                "'{0}'. Matching files were: '{1}'"
                                .format(prefix, files))
            else:
                raise Error("Unexpected number of files matching the prefix "
                            "'{0}'. Matching files were: '{1}'"
                            .format(prefix, files))

            fasta_files = []
            for file in files:
                # Decompress each gzip file
                fname = os.path.basename(file)
                name, ext = os.path.splitext(fname)

                if ext in [".gz", ".gzip"]:
                    out = os.path.join(self.dir, name)

                else:
                    out = os.path.join(self.dir, fname)

                # Keep track of the files in the temporary folder
                fasta_files.append(out)

            entry["files"] = tuple(fasta_files)

            self.data.append(entry)

    def run(self):
        for entry in self.data:
            for file, out in zip(entry["original_files"], entry["files"]):
                # Decompress each gzip file
                fname = os.path.basename(file)
                name, ext = os.path.splitext(fname)

                if ext in [".gz", ".gzip"]:
                    log.info("Decompressing '%s' into temporary folder", fname)
                    os.popen('zcat {0} > {1}'.format(file, out))
                # with gzip.open(file) as infile:
                # with open(out, 'w') as outfile:
                # while True:
                ## Read/Write 1 MB at a time
                # data = infile.read(1024 ** 2)
                # if not data:
                # break
                # outfile.write(data)
                else:
                    log.info("Copying '%s' to temporary folder", fname)
                    # This should only be fastq files
                    shutil.copyfile(file, out)
            if self.config.params["nextseq"] == True:
                log.info("Joining lanes R1 for sample {0}".format(entry["Sample_ID"]))
                self.join_lanes(entry["Sample_ID"], "R1")
                if len(entry["files"]) == 8:
                    log.info("Joining lanes R2 for sample {0}".format(entry["Sample_ID"]))
                    self.join_lanes(entry["Sample_ID"], "R2")

    def join_lanes(self, sampleID, direction):
        if not os.path.exists(os.path.join(self.dir, 'Merged_lanes')):
            os.makedirs(os.path.join(self.dir, 'Merged_lanes'))
        os.system("cd {0}; cat {1}_*_{2}_* > {0}/Merged_lanes/{1}_{2}.fastq ".format(self.dir, sampleID, direction))


class DataCollector:
    """Reaps all objects used in the pipeline and formats data in a friendly
    memory representation for text/html output
    """

    def __init__(self, pipeline):
        self.pipe = pipeline
        self.output = []
        self.demux = {}
        self.links = set([])
        self.stats_all_lanes = {}
        self.collect()

    def collect(self):
        columns = (
            "Sample ID",
            "Filename",
            "Q20",
            "Q30",
            "%Q30",
            "Total bases",
            "Total reads",
            "FastQC",
            "Primary organism",
            "Contaminants",
            "Total mapped bases",
            "Mean coverage",
            "Theoretical coverage Q20",
            "Theoretical coverage Q30",
            "Total mapped reads",
            "Total reads (used)",
            "Qualimap",
            "Demux read mean",
            "Demux read stdev",
            "Demux reads failed",
            "Demux",
            "Merged reads",
            "Merged reads failed",
        )

        self.output.append(columns)
        # These columns (index number) should be links in HTML output
        self.links.update([7, 9, 16, 20])  # , 20])

        sum_q20 = 0
        sum_q30 = 0
        sum_bases = 0
        sum_reads = 0
        sum_mapped_bases = 0
        sum_mapped_reads = 0
        sum_total_reads_used = 0
        multiqc_file = self.pipe.data.config.params["output_url"] + '/multiqc_report.html'

        for entry in self.pipe.data.config.config:
            sample_id = entry["Sample_ID"]
            self.stats_all_lanes[sample_id] = OrderedDict()
            self.stats_all_lanes[sample_id]['q20'] = 0
            self.stats_all_lanes[sample_id]['q30'] = 0
            self.stats_all_lanes[sample_id]['percent_q30'] = 0
            self.stats_all_lanes[sample_id]['total_bases'] = 0
            self.stats_all_lanes[sample_id]['total_reads'] = 0
            self.stats_all_lanes[sample_id]['fastqc'] = '...'
            self.stats_all_lanes[sample_id]['primary'] = ''
            self.stats_all_lanes[sample_id]['contam'] = ''
            self.stats_all_lanes[sample_id]['mapped_bases'] = ''
            self.stats_all_lanes[sample_id]['mean_cov'] = ''
            self.stats_all_lanes[sample_id]['cov_q20'] = ''
            self.stats_all_lanes[sample_id]['cov_q30'] = ''
            self.stats_all_lanes[sample_id]['mapped_reads'] = ''
            self.stats_all_lanes[sample_id]['used_reads'] = ''
            self.stats_all_lanes[sample_id]['qualimap'] = ''
            self.stats_all_lanes[sample_id]['demux_mean'] = ''
            self.stats_all_lanes[sample_id]['demux_stdev'] = ''
            self.stats_all_lanes[sample_id]['demux_failed'] = ''
            self.stats_all_lanes[sample_id]['demux'] = ''
            self.stats_all_lanes[sample_id]['merged_reads'] = ''
            self.stats_all_lanes[sample_id]['merged_failed'] = ''

            for filename in entry["files"]:
                if filename.endswith("I1_001.fastq"):
                    continue
                fname = os.path.basename(filename)

                quality = self.pipe.qual.data[filename]
                q20 = str(quality[20])
                q30 = str(quality[30])
                total_bases = str(quality["total_bases"])
                total_reads = str(quality["total_reads"])
                sum_q20 += int(q20)
                sum_q30 += int(q30)
                sum_bases += int(total_bases)
                sum_reads += int(total_reads)

                self.stats_all_lanes[sample_id]['q20'] += int(q20)
                self.stats_all_lanes[sample_id]['q30'] += int(q30)
                self.stats_all_lanes[sample_id]['total_bases'] += int(total_bases)
                self.stats_all_lanes[sample_id]['total_reads'] += int(total_reads)

                percent_q30 = str(round((quality[30] / quality["total_bases"]) * 100, 2))

                fastqc_file = self.pipe.qc.data[filename].replace(
                    self.pipe.data.dir,
                    self.pipe.data.config.params["output_url"],
                )

                # jcostaDamage
                if filename in self.pipe.tax.data:

                    try:
                        cont_file, primary_org = self.pipe.tax.data[filename]
                    except (KeyError, TypeError):
                        log.debug("Filename '%s' has no contaminants "
                                  "information, ignoring.", filename)
                        contaminants_file = primary_org = ""
                    else:
                        if primary_org == "=":
                            contaminants_file = primary_org
                        else:
                            contaminants_file = cont_file.replace(
                                self.pipe.data.dir,
                                self.pipe.data.config.params["output_url"],
                            )
                            self.stats_all_lanes[sample_id]['primary'] = primary_org
                            if os.path.basename(contaminants_file) != 'skipped':
                                self.stats_all_lanes[sample_id]['contam'] = '<a href={0}>Report</a>'.format(
                                    contaminants_file)
                            else:
                                self.stats_all_lanes[sample_id]['contam'] = 'skipped'


                else:
                    contaminants_file = '='
                    primary_org = '='

                try:
                    ref = entry["Reference"]
                except KeyError:
                    ref = None

                if ref and filename in self.pipe.cov.data:
                    cov = self.pipe.cov.data[filename]
                    total_mapped_bases = cov["total_mapped_bases"]
                    self.stats_all_lanes[sample_id]['mapped_bases'] = total_mapped_bases
                    mean_coverage = cov["mean_coverage"]
                    try:
                        self.stats_all_lanes[sample_id]['mean_cov'] = "{0:.2f}".format(float(mean_coverage))
                    except ValueError:
                        self.stats_all_lanes[sample_id]['mean_cov'] = mean_coverage
                    try:
                        theoretical_cov_Q20 = "{0:.2f}".format(quality[20] / cov["ref_length"])
                    except ValueError:
                        theoretical_cov_Q20 = 'ERROR'

                    try:
                        theoretical_cov_Q30 = "{0:.2f}".format(quality[30] / cov["ref_length"])
                    except ValueError:
                        theoretical_cov_Q30 = 'ERROR'

                    total_mapped_reads = cov["total_mapped_reads"]
                    self.stats_all_lanes[sample_id]['mapped_reads'] = total_mapped_reads
                    total_reads_used = cov["total_reads_used"]
                    self.stats_all_lanes[sample_id]['used_reads'] = total_reads_used
                    qualimap_file = cov["qualimap"]
                    self.stats_all_lanes[sample_id]['qualimap'] = qualimap_file  # get it as it is (tmp path)

                    if not qualimap_file in ("=", "ERROR"):
                        sum_mapped_bases += int(total_mapped_bases)
                        sum_mapped_reads += int(total_mapped_reads.split()[0])
                        sum_total_reads_used += int(total_reads_used)

                else:
                    total_mapped_bases = total_mapped_reads = ""
                    total_reads_used = qualimap_file = mean_coverage = ""
                    theoretical_cov_Q20 = ''
                    theoretical_cov_Q30 = ''

                if qualimap_file != "=":
                    qualimap_file = qualimap_file.replace(
                        self.pipe.data.dir,
                        self.pipe.data.config.params["output_url"],
                    )
                    if qualimap_file != '':
                        self.stats_all_lanes[sample_id]['qualimap'] = "<a href=\"{0}\">Report</a>".format(
                            qualimap_file)  # convert it to proper address

                if filename in self.pipe.dmux.data:
                    if self.pipe.dmux.data[filename] == "=":
                        demux_mean = "="
                        demux_stdev = "="
                        demux_read_failed = "="
                        demux_url = "="
                    else:
                        samples, demux = self.pipe.dmux.data[filename]

                        dmux_read = []
                        demux_mean = 0
                        demux_stdev = 0
                        demux_read_failed = 0

                        total = len(demux)

                        log.debug("Demux samples data is %s", samples)
                        log.debug("Demux data is %s", demux)
                        counts = demux["options"]["read_counts"]

                        for key1 in counts:
                            if key1 == "rejected":
                                demux_read_failed += counts["rejected"]
                                continue

                            for key2 in counts[key1]:
                                dmux_read.append(counts[key1][key2])

                        demux_read_failed = str(demux_read_failed)
                        demux_mean = "{0:.2f}".format(sum(dmux_read) / total)
                        demux_stdev = "{0:.2f}".format(std(dmux_read, ddof=1))

                        demux_file = os.path.join(
                            demux["options"]["samplepath"],
                            "demux.html"
                        )

                        self.collect_demux_merge(filename, demux_file)

                        demux_url = demux_file.replace(
                            self.pipe.data.dir,
                            self.pipe.data.config.params["output_url"],
                        )
                else:
                    demux_mean = ""
                    demux_stdev = ""
                    demux_read_failed = ""
                    demux_url = ""

                if filename in self.pipe.merge.data:
                    if self.pipe.merge.data[filename] == "=":   # is muda para == ?
                        merged_reads = "="
                        merged_reads_failed = "="
                    else:
                        merged_reads = 0
                        merged_reads_failed = 0

                        for val in self.pipe.merge.data[filename].itervalues():
                            merged_reads += val["reads_combined"]
                            merged_reads_failed += val["reads_notcombined"]

                        merged_reads = str(merged_reads)
                        merged_reads_failed = str(merged_reads_failed)
                        self.stats_all_lanes[sample_id]['merged_reads'] = merged_reads
                        self.stats_all_lanes[sample_id]['merged_failed'] = merged_reads_failed
                else:
                    merged_reads = ""
                    merged_reads_failed = ""

                self.output.append([
                    sample_id,
                    fname,
                    q20,
                    q30,
                    percent_q30,
                    total_bases,
                    total_reads,
                    fastqc_file,
                    primary_org,
                    contaminants_file,
                    total_mapped_bases,
                    mean_coverage,
                    theoretical_cov_Q20,
                    theoretical_cov_Q30,
                    total_mapped_reads,
                    total_reads_used,
                    qualimap_file,
                    demux_mean,
                    demux_stdev,
                    demux_read_failed,
                    demux_url,
                    merged_reads,
                    merged_reads_failed,
                ])

            if qualimap_file != "":
                self.stats_all_lanes[sample_id]['cov_q20'] = "{0:.2f}".format(
                    float(self.stats_all_lanes[sample_id]['q20']) / float(cov["ref_length"]))
                self.stats_all_lanes[sample_id]['cov_q30'] = "{0:.2f}".format(
                    float(self.stats_all_lanes[sample_id]['q30']) / float(cov["ref_length"]))

        self.output.append(['Total',
                            '',
                            str(sum_q20),
                            str(sum_q30),
                            "{0:.2f}".format(100 * sum_q30 / float(sum_bases)),
                            str(sum_bases), str(sum_reads),
                            multiqc_file, '', '',
                            str(sum_mapped_bases),
                            '', '', '',

                            str(sum_mapped_reads),
                            str(sum_total_reads_used),
                            ''])

    def collect_demux_merge(self, filename, demux_file):
        columns = (
            "Demux ID",
            "Seq.F",
            "Seq.R",
            "Demuxed reads",
            "Demux failed reads",
            "Merged reads",
            "Merge failed reads",
        )

        data = []

        data.append(columns)

        files, demux = self.pipe.dmux.data[filename]

        log.debug("Preparing demux data %s for file %s", demux, filename)
        samples = demux["options"]["samples"]
        stats = demux["options"]["read_counts"]

        for i, sample_id in enumerate(sorted(samples)):
            seqF, seqR = samples[sample_id]

            # Demuxed files
            demuxed_reads = str(stats[seqF][seqR])
            if i == 0:
                # First line of the table
                demux_failed_reads = str(stats["rejected"])
            else:
                demux_failed_reads = "="

            # Merged files
            if filename in self.pipe.merge.data:
                vals = self.pipe.merge.data[filename][sample_id]
                reads_combined = str(vals["reads_combined"])
                reads_notcombined = str(vals["reads_notcombined"])
            else:
                reads_combined = ""
                reads_notcombined = ""

            values = (sample_id, seqF, seqR, demuxed_reads, demux_failed_reads,
                      reads_combined, reads_notcombined)
            log.debug("Demux file values %s", values)
            data.append(values)

        self.demux[filename] = (demux_file, data)





