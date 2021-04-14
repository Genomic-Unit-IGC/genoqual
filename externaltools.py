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
import subprocess
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

class Process:
	def __init__(self, cmd):
		self.cmd = cmd

	def run(self, stdout=PIPE, stderr=PIPE, stdin=None, trust_exitcode=True,
	        env=None, piped=False):
		"""Launch or specified command in a subprocess

		:arg trust_exitcode: raises IOError if command exits with non-zero code
		:arg env: dictionary with shell environment to use
		:arg piped: don't wait for command to finish to allow piping output to
		            subsequent commands

		If piped is True trust_exitcode has no effect.
		"""
		log.debug("Launching subprocess %r", self.cmd)

		if isinstance(stdin, str):
			_in = PIPE
			input = stdin

			if piped:
				raise TypeError("Cannot mix stdin stdin='string' with "
				                "piped=True")
		else:
			_in = stdin
			input = None

		try:
			self.p = Popen(self.cmd, stdout=stdout, stderr=stderr, stdin=_in,
			               env=env)
		except OSError:
			log.error("Command %r failed to execute. Is the file in the path?",
			          self.cmd)
			raise

		# If stdin is a pipe from another process, close it to allow the first
		# process to receive SIGPIPE if the second process finishes before.
		try:
			name = stdin.name
		except AttributeError:
			pass
		else:
			if name == "<fdopen>":
				stdin.close()

		# Don't wait for process to finish if it's going to be piped with
		# another command.
		if piped:
			return

		out, err = self.p.communicate(input)
		self.finish()

		if trust_exitcode and self.p.returncode != 0:
			raise IOError("Command {0!r} failed with exit code {1}.\n"
			              "Stdout:\n{2}\n Stderr:\n{3}\n"
			              .format(self.cmd, self.p.returncode, out, err))

		log.debug("Subprocess finished with output:\n"
		          "Stdout:\n%s\nStderr:\n%s\n", out, err)

		return out, err

	def finish(self):
		# Wait for processes to finish to avoid zombies
		self.p.wait()

class BaseTool:
	"""Base skeleton for external tools or pipeline components.
	Implements a check phase to validate requirements prior to execution
	"""
	def __init__(self, datacfg):
		self.datacfg = datacfg
		self.data = {}

	def prepare(self):
		"""Perform pre-flight checks to ensure all necessary resources are
		available and settings are correctly defined
		"""
		raise NotImplementedError("Subclasses need to override this method")

	def run(self):
		raise NotImplementedError("Subclasses need to override this method")


class QualityStats(BaseTool):
	"""Perform calculation of q20 and q30 statistics
	"""
	def prepare(self):
		log.info("Preparing files for q20 and q30 analysis")
		for entry in self.datacfg.data:
			for file in entry["files"]:
				if file.endswith("I1_001.fastq"):
					continue # Skips the index file
				self.data[file] = None

	def run(self):
		self._compute_q(20, 30)

	def _compute_q(self, *args):
		for infile in self.data:
			filename = os.path.basename(infile)

			log.info("Calculating quality scores for file '%s'", filename)
			result = {x: 0 for x in args}
			total_bases = 0
			total_reads = 0
			with open(infile) as fh:
				for read in FastqReader(fh):
					total_bases += read.qual.size
					total_reads += 1
					for val in result:
						result[val] += read.qual[read.qual >= val].size

			# Keep result of calculation
			result["total_bases"] = total_bases
			result["total_reads"] = total_reads
			log.debug("Quality scores for file '%s': %r", filename, result)

			self.data[infile] = result


class MultiQC(BaseTool):
	"""Run MultiQC and write a summary report in the main results folder
	"""
	def __init__(self, *args, **kwargs):
		super(MultiQC, self).__init__(*args, **kwargs)

		self.destdir = self.datacfg.dir

	def prepare(self):
		log.info("Performing pre-flight checks on MultiQC")
		cmd = ("multiqc", "--version")
		p = Process(cmd)
		out, err = p.run()

		if err and not out.startswith(b"multiqc, version"):
			raise IOError("Unexpected output from MultiQC wrapper.\n"
			              "Stdout:\n{0}\n Stderr:\n{1}\n"
			              .format(out, err))

		log.info("MultiQC pre-flight checks finished successfully")

	def run(self):
		log.info("Running MultiQC")
		try:
			cmd = [
			        "multiqc",
			        "--no-data-dir", "-o", self.destdir, self.destdir
			]
			p = Process(cmd)
			out, err = p.run()

		except:
			raise IOError("MultiQC error.\n"
			              "Stdout:\n{0}\n Stderr:\n{1}\n"
			              .format(out, err))


class Qiime(BaseTool):

    def __init__(self, *args, **kwargs):
        super(Qiime, self).__init__(*args, **kwargs)

        self.destdir = os.path.join(self.datacfg.dir, "Qiime")
        self.infile = None
        self.infile_rev = None
        self.infile_fwd = None
        self.indexfile = None
        self.metafile = None
        self.users = []
        self.inputdir = None  # Useful only in the merged R1+R2 case
        self.original_index = None  #

    def prepare(self):
        log.info("Performing pre-flight checks on Qiime")
        cmd = ("qiime", "--version")
        p = Process(cmd)
        out, err = p.run()

        if err or not out.startswith(b"q2cli version"):
            raise IOError("Unexpected output from Qiime wrapper.\n"
                          "Stdout:\n{0}\n Stderr:\n{1}\n"
                          .format(out, err))

        cmd = ("biom")
        p = Process(cmd)
        out, err = p.run()

        if err or not out.startswith(b"Usage:"):
            raise IOError("Unexpected output from Qiime wrapper.\n"
                          "Stdout:\n{0}\n Stderr:\n{1}\n"
                          .format(out, err))

        for entry in self.datacfg.data:
            if not entry["Merge_reads"]:
                log.info('QIIME on R1')
                for file in entry["files"]:
                    if file.endswith("I1_001.fastq"):
                        self.indexfile = file
                        self.original_index = file
                    elif file.endswith("R1_001.fastq"):
                        self.infile = file
                        self.infile_fwd = file
                    elif file.endswith("R2_001.fastq"):
                        self.infile_rev = file
            else:
                log.info('QIIME on merged')
                self.inputdir = os.path.join(self.datacfg.dir, "Merged_reads", entry["Sample_ID"])
                for file in entry["files"]:
                    if file.endswith("R1_001.fastq"):
                        self.infile = self.inputdir + '/' + os.path.basename(file + ".extendedFrags.fastq")
                        self.infile_fwd = file
                        log.debug("Original R1 file should be:")
                        log.debug(self.infile_fwd)
                        log.debug("Merged file should be:")
                        log.debug(self.infile)
                    elif file.endswith("I1_001.fastq"):
                        self.original_index = file
                        self.indexfile = self.inputdir + '/' + os.path.basename(file.replace(".fastq", "_fixed.fastq"))
                        log.debug("Index file should be:")
                        log.debug(self.indexfile)
                    elif file.endswith("R2_001.fastq"):
                        self.infile_rev = file
                        log.debug("Original R2 file should be:")
                        log.debug(self.infile_rev)

        if not os.path.isfile(self.datacfg.config.params["metadata"]):
            raise Error("Metadata file not found at '{0}'"
                        .format(self.datacfg.config.params["metadata"]))
        else:
            self.metafile = self.datacfg.config.params["metadata"]
            log.info("Qiime pre-flight checks finished successfully")

       # out, err = self.check_metafile(self.metafile, os.path.join(self.datacfg.dir, "validate_metadata"))    # Metadata validation isn't necessary
       # if err or not out.startswith(b"No errors"):
       #     log.error("Error while validating metadata file.")
       #     log_metadata = open(glob(self.datacfg.dir + "/validate_metadata/*.log")[0]).read()
       #     log.error(log_metadata)
       #     sys.exit(1)

       #     log.info("Metadata file validated")
       #     shutil.rmtree(self.datacfg.dir + "/validate_metadata/")
       # log.info("Qiime pre-flight checks finished successfully")

    def run(self):
        log.info("Entering QIIME")
        if not os.path.isdir(self.destdir):
            os.mkdir(self.destdir)

        # Start the Qiime
        log.info("Writing parameters file")
        with open('%s/params.txt' % self.destdir, 'w') as params:
            params.write("pick_otus:enable_rev_strand_match\tTrue\n")
            params.write("biom-summarize-table:qualitative\tTrue\n")
            params.write("make_emperor:ignore_missing_samples\tTrue\n")
            params.write("beta_diversity_through_plots:ignore_missing_samples\tTrue\n")
            params.write(
                "beta_diversity:metrics\tbray_curtis,euclidean,unweighted_unifrac,weighted_unifrac,binary_jaccard,abund_jaccard\n")
            params.write("alpha_diversity:metrics\tobserved_species,chao1,shannon,PD_whole_tree\n")
        paramfile = self.destdir + '/params.txt'

        if self.inputdir != None:
            log.info("Correcting barcodes in order to match the merged fastq file...")
            self.fix_header_pairs()

        if self.indexfile != None and self.infile != None and self.metafile != None:
            log.debug(self.indexfile)
            log.debug(self.infile)
            log.debug(self.metafile)
            log.debug("Checking metafile for whitespaces...")
            self.metafile_correct(self.metafile)

            # Splitting libraries...
            try:
                out, err = self.split_libraries(self.infile, self.indexfile, self.metafile,
                                                self.destdir + "/temp_slout")
            except:
                log.info('Splitting failed. Trying to split without reversing barcodes...')
                shutil.rmtree(self.destdir + "/temp_slout")
                out, err = self.split_libraries(self.infile, self.indexfile, self.metafile,
                                                self.destdir + "/temp_slout", rev=0)

            # Reads separation by user...
            log.info("Starting separation by user")
            cmd = [
                "split_sequence_file_on_sample_ids.py",
                "-o", self.destdir + "/temp_slout/out_fasta",
                "-i", self.destdir + "/temp_slout/seqs.fna",
                "--file_type", "fasta"]
            p = Process(cmd)
            out, err = p.run()

            fastapath = self.destdir + "/temp_slout/out_fasta/"
            fastafiles = glob(fastapath + "*.fasta")
            for f in fastafiles:
                user = os.path.basename(f).split('.')[0]
                if user not in self.users:
                    self.users.append(user)
            for user in self.users:
                user_fastapath = self.destdir + '/' + user + '/slout'
                user_fastqpath = self.destdir + '/' + user + '/fastq'
                os.makedirs(user_fastapath)
                os.makedirs(user_fastqpath)
                os.system('cat %s/%s*.fasta > %s/seqs.fna' % (fastapath, user, user_fastapath))
                log.info('QIIME analysis for user %s' % user)

                # OTUs picking...
                try:
                    out, err = self.pick_otus(self.destdir + "/%s/slout/seqs.fna" % user,
                                              self.destdir + "/%s/otus" % user, paramfile)
                except:
                    log.info('OTU picking failed for user %s, with error:' % user)
                    log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))
                    continue

                # Biom summarizing...
                try:
                    out, err = self.biom_summarize(
                        self.destdir + "/%s/otus/otu_table_mc2_w_tax_no_pynast_failures.biom" % user,
                        self.destdir + "/%s/counts.txt" % user)
                except:
                    log.info('BIOM summarize failed for user %s, with error:' % user)
                    log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))
                    continue

                # Finding minimum count for rarefaction...(min:1000)
                min_count = 1000
                if os.path.exists(self.destdir + "/%s/counts.txt" % user):
                    c = open(self.destdir + "/%s/counts.txt" % user).readlines()
                    for line in c:
                        if line.startswith('%s.' % user):
                            min_count_tmp = int(float(line.split()[1]))
                            if min_count_tmp > min_count:
                                min_count = min_count_tmp
                                break

                # Core analyses...
                try:
                    out, err = self.core_analyses(
                        self.destdir + "/%s/otus/otu_table_mc2_w_tax_no_pynast_failures.biom" % user,
                        self.destdir + "/%s/cdout" % user, self.metafile,
                        self.destdir + "/%s/otus/rep_set.tre" % user, str(min_count), paramfile)
                except:
                    log.info('Core analyses failed for user %s, with error:' % user)
                    log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))
                    continue

            # Splitting again, without any quality constraint, to give the user the fastq files to repeat
            # the analysis.

            try:
                out, err = self.split_libraries(self.infile_fwd, self.indexfile, self.metafile,
                                                self.destdir + "/temp_slout_unfiltered", n=999, q=0, fastq=1)
            except:
                log.info('Splitting failed. Trying to split without reversing barcodes...')
                shutil.rmtree(self.destdir + "/temp_slout_unfiltered")
                out, err = self.split_libraries(self.infile_fwd, self.indexfile, self.metafile,
                                                self.destdir + "/temp_slout_unfiltered", n=999, q=0, rev=0, fastq=1)

            if self.infile_rev != None:  # i.e. input is paired-end
                log.info('Pairing reads after splitting into fastq format...')
                self.fix_reads_pairs()
            log.info('Splitting fastq files for users...')
            self.split_fastq()
            shutil.rmtree(self.destdir + "/temp_slout")
            shutil.copy(self.destdir + "/temp_slout_unfiltered/split_library_log.txt",
                        self.destdir + "/demultiplex_log.txt")
            shutil.rmtree(self.destdir + "/temp_slout_unfiltered")

        else:
            log.error("Could not find either the reads or the index file:"
                      "Reads:\n{0}\n Index:\n{1}\n".format(self.infile, self.indexfile))
            raise

    def split_libraries(self, infile, index, meta, dest, n=0, q=19, rev=1, fastq=0):
        cmd = [
            "split_libraries_fastq.py",
            "-o", dest,
            "-i", infile,
            "-b", index,
            "-m", meta,
            "-n", str(n),
            "-q", str(q)]

        if rev == 1:
            cmd = cmd + ["--rev_comp_mapping_barcodes", "--rev_comp_barcode"]
        if fastq == 1:
            cmd = cmd + ["--store_demultiplexed_fastq"]

        p = Process(cmd)
        out, err = p.run()
        return out, err

    def pick_otus(self, i, o, paramfile):
        log.info('Entered otu pick')
        cmd = [
            "pick_open_reference_otus.py",
            "-i", i,
            "-o", o,
            "-p", paramfile]
        p = Process(cmd)
        out, err = p.run()
        return out, err

    def check_metafile(self, i, outdir):
        log.info('Checking metafile...')
        cmd = [
            "validate_mapping_file.py",
            "-m", i,
            "-o", outdir]
        p = Process(cmd)
        out, err = p.run()
        return out, err

    def biom_summarize(self, i, o):
        cmd = [
            "biom", "summarize-table",
            "-i", i,
            "-o", o]
        p = Process(cmd)
        out, err = p.run()
        return out, err

    def core_analyses(self, i, o, meta, tree, min_count, params):
        cmd = [
            "core_diversity_analyses.py",
            "-i", i,
            "-o", o,
            "-m", meta,
            "-t", tree,
            "-e", min_count,
            "-p", params]
        p = Process(cmd)
        out, err = p.run()
        return out, err

    def split_fastq(self):
        users = []
        records_f = SeqIO.parse(open(self.destdir + "/temp_slout_unfiltered/seqs.fastq", "rU"), "fastq")
        if os.path.exists(self.destdir + "/temp_slout_unfiltered/seqs_rev.fastq"):
            records_r = SeqIO.parse(open(self.destdir + "/temp_slout_unfiltered/seqs_rev.fastq", "rU"), "fastq")
            for (forward, reverse) in izip(records_f, records_r):
                forward.description = ' '.join(forward.description.split()[1:3])
                user = forward.id.split('.')[0]
                if not user in users:
                    users.append(user)
                sample = '_'.join(forward.id.split('_')[0].split('.')[1:])
                forward.id = forward.description
                forward.description = ''
                assert forward.id.split()[0] == reverse.description.split()[0]
                user_fastqpath = self.destdir + '/%s/fastq' % user
                handleout_f = open(user_fastqpath + "/%s_L001_R1_001.fastq" % sample, "ab")
                handleout_r = open(user_fastqpath + "/%s_L001_R2_001.fastq" % sample, "ab")
                SeqIO.write(forward, handleout_f, 'fastq')
                SeqIO.write(reverse, handleout_r, 'fastq')
                handleout_f.close()
                handleout_r.close()

        else:
            log.debug("Found seqs.fastq only")
            for forward in records_f:
                forward.description = ' '.join(forward.description.split()[1:3])
                user = forward.id.split('.')[0]
                if not user in users:
                    users.append(user)
                sample = '_'.join(forward.id.split('_')[0].split('.')[1:])
                forward.id = forward.description
                forward.description = ''
                user_fastqpath = self.destdir + '/%s/fastq' % user
                handleout_f = open(user_fastqpath + "/%s_L001_R1_001.fastq" % sample, "ab")
                SeqIO.write(forward, handleout_f, 'fastq')
                handleout_f.close()

        counts = open(self.destdir + "/temp_slout_unfiltered/split_library_log.txt").readlines()
        for l in counts:
            user = l.split('.')[0]
            if user in users:
                fileout = open(self.destdir + '/%s/fastq/counts.txt' % user, 'a')
                fileout.write(l)
                fileout.close()

        return

    def fix_header_pairs(self):
        handle = open(self.infile, "rU")
        headers = []
        for record in SeqIO.parse(handle, "fastq"):
            headers.append(record.description)
        handle.close()
        headers = set(headers)
        handle = open(self.original_index, "rU")
        handleout = open(self.indexfile, "wb")
        for record in SeqIO.parse(handle, "fastq"):
            if record.description in headers:
                SeqIO.write(record, handleout, 'fastq')
        handle.close()
        handleout.close()
        return

    def fix_reads_pairs(self):
        handle = open(self.destdir + "/temp_slout_unfiltered/seqs.fastq", "rU")
        headers = []
        for record in SeqIO.parse(handle, "fastq"):
            headers.append(record.description.split()[1])
        handle.close()
        headers = set(headers)
        handle = open(self.infile_rev, "rU")
        handleout = open(self.destdir + "/temp_slout_unfiltered/seqs_rev.fastq", "wb")
        for record in SeqIO.parse(handle, "fastq"):
            if record.description.split()[0] in headers:
                SeqIO.write(record, handleout, 'fastq')
        handle.close()
        handleout.close()
        return

    def metafile_correct(self, filename):
        # Replaces whitespace with dot
        filedata = None
        with open(filename, 'r') as file:
            filedata = file.read()
        # Replace the target string
        filedata = filedata.replace(' ', '.')
        filedata = filedata.replace('_', '.')
        # Write the file out again
        with open(filename, 'w') as file:
            file.write(filedata)
        return


class QiimeOnControls(Qiime):
    def __init__(self, *args, **kwargs):
        super(QiimeOnControls, self).__init__(*args, **kwargs)
        self.results_path = self.datacfg.config.params["base_path"] + 'results/'
        self.input_path = self.datacfg.config.params["base_path"] + 'input/'
        self.current_run = self.datacfg.config.params["run_folder"]

    def run(self):
        log.info('Starting historical QIIME analysis of controls...')
        metaruns = []
        # Create a big single metadata file for all controls. Duplicate names are avoided by adding
        # the run ID to the samplename. Duplicate barcodes are not fixed.
        os.makedirs(self.destdir + '/Controls_analysis')
        control_metadata = open(self.destdir + '/Controls_analysis/GEU_metadata.csv', 'w')
        # control_metadata = open('/home/mtruglio/Desktop/GEU_metadata.csv', 'w')
        control_metadata.write('#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tControlType\tDescription\n')

        # Also, create a big single seqs.fna file from all GEU runs. The headers will be modified according to the
        # metadata, again in order to avoid duplicate names.
        # control_seqs = open('/home/mtruglio/Desktop/GEU_seqs.fna', 'wb')
        control_seqs = open(self.destdir + '/Controls_analysis/GEU_seqs.fna', 'wb')

        for item in os.listdir(self.results_path):
            if os.path.exists(os.path.join(self.results_path, item) + '/Qiime/GEU/slout/seqs.fna'):
                metaruns.append(item)
        metaruns += [self.current_run]
        for item in metaruns:
            metain = open(self.input_path + item + '/metadata.csv').readlines()

            for line in metain:
                if line.startswith('GEU'):
                    samplename = line.split()[0] + '.' + item.split('_')[1]  # make it unique
                    if ('Pos' in samplename) or ('pos' in samplename):
                        controltype = 'Pos'
                    else:
                        controltype = 'Neg'

                    if item == self.current_run:
                        controltype += '_current'

                    control_metadata.write(samplename + '\t' + '\t'.join(
                        line.split()[1:3]) + '\t' + controltype + '\t' + samplename + '\n')

            if item != self.current_run:
                records_fna = SeqIO.parse(
                    open(os.path.join(self.results_path, item) + '/Qiime/GEU/slout/seqs.fna', "rU"), "fasta")
            else:
                records_fna = SeqIO.parse(open(self.destdir + '/GEU/slout/seqs.fna', "rU"), "fasta")

            for r in records_fna:
                entry_n = r.id.split('_')[-1]
                r.id = r.id.split('_')[0] + '.%s' % item.split('_')[1] + '_%s ' % entry_n + ' '.join(
                    r.description.split()[1:])
                r.description = ''
                SeqIO.write(r, control_seqs, 'fasta')
        control_metadata.close()
        control_seqs.close()

        try:
            out, err = self.pick_otus(self.destdir + '/Controls_analysis/GEU_seqs.fna',
                                      self.destdir + '/Controls_analysis/otus', self.destdir + '/params.txt')
        except:
            log.info('OTU picking failed, with error:')
            log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))

        try:
            out, err = self.biom_summarize(
                self.destdir + "/Controls_analysis/otus/otu_table_mc2_w_tax_no_pynast_failures.biom",
                self.destdir + "/Controls_analysis/counts.txt")
        except:
            log.info('BIOM summarize failed, with error:')
            log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))

        counts = open(self.destdir + "/Controls_analysis/counts.txt").readlines()
        current_run_counts = []
        for line in counts:
            if line.startswith('GEU') and line.split()[0].split('.')[-1].strip(':') == self.current_run.split('_')[1]:
                c = int(float(line.split()[-1]))
                current_run_counts.append(c)
        min_count = min(current_run_counts)

        try:
            out, err = self.beta_diversity(
                self.destdir + "/Controls_analysis/otus/otu_table_mc2_w_tax_no_pynast_failures.biom",
                self.destdir + "/Controls_analysis/bdiv",
                self.destdir + "/Controls_analysis/otus/rep_set.tre",
                self.destdir + '/Controls_analysis/GEU_metadata.csv', min_count,
                self.destdir + '/params.txt')
        except:
            log.info('Beta diversity failed, with error:')
            log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))

        try:
            out, err = self.taxa_plots(
                self.destdir + "/Controls_analysis/otus/otu_table_mc2_w_tax_no_pynast_failures.biom",
                self.destdir + "/Controls_analysis/taxa_summary",
                self.destdir + '/Controls_analysis/GEU_metadata.csv',
                self.destdir + '/params.txt')
        except:
            log.info('Taxa plotting failed, with error:')
            log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))

    def beta_diversity(self, i, o, tree, meta, min_count, params):
        if min_count < 80:
            min_count = 80
        cmd = [
            "beta_diversity_through_plots.py",
            "-i", i,
            "-o", o,
            "-m", meta,
            "-t", tree,
            "-e", str(min_count),
            "-p", params]
        p = Process(cmd)
        out, err = p.run()
        return out, err

    def taxa_plots(self, i, o, meta, params):
        cmd = [
            "summarize_taxa_through_plots.py",
            "-i", i,
            "-o", o,
            "-m", meta,
            "-p", params]
        p = Process(cmd)
        out, err = p.run()
        return out, err


class FastQC(BaseTool):
    """Run FastQC and store the results in folders with names matching the
    input files
    """

    def __init__(self, *args, **kwargs):
        super(FastQC, self).__init__(*args, **kwargs)

        self.destdir = os.path.join(self.datacfg.dir, "FastQC")

    def prepare(self):
        log.info("Performing pre-flight checks on FastQC")
        cmd = ("fastqc", "-v")
        p = Process(cmd)
        out, err = p.run()

        if err or not out.startswith(b"FastQC v"):
            raise IOError("Unexpected output from FastQC wrapper.\n"
                          "Stdout:\n{0}\n Stderr:\n{1}\n"
                          .format(out, err))

        log.info("FastQC pre-flight checks finished successfully")

        log.info("Preparing files for FastQC analysis")
        for entry in self.datacfg.data:
            # Keep track of the directory where FastQC stored the output
            for file in entry["files"]:
                if file.endswith("I1_001.fastq"):
                    continue
                fname = os.path.basename(file)
                report = os.path.join(self.destdir, os.path.splitext(fname)[0] + "_fastqc.html")
                self.data[file] = report

    def run(self):
        if not os.path.isdir(self.destdir):
            os.mkdir(self.destdir)

        for infile in self.data:
            # jcostaDamage
            try:
                fname = os.path.basename(infile)
                log.info("Running FastQC on file '%s'", fname)
                cmd = [
                    "fastqc",
                    "-t", self.datacfg.config.params["args"].threads,
                    "-o", self.destdir,
                    infile
                ]
                p = Process(cmd)
                out, err = p.run()

                zip_file = os.path.dirname(self.data[infile]) + ".zip"
                log.debug("Removing FastQC zipped output '%s'", zip_file)
                os.remove(zip_file)
            except:
                continue


class TaxonomyClassifier(BaseTool):
    """Use a Seqtk and Blast to estimate the taxonomical distribution of the
    sampled reads
    """

    def __init__(self, *args, **kwargs):
        super(TaxonomyClassifier, self).__init__(*args, **kwargs)

        self.destdir = os.path.join(self.datacfg.dir, "Contaminants")

    def prepare(self):
        log.info("Performing pre-flight checks on Seqtk")
        cmd = ("seqtk",)
        p = Process(cmd)
        out, err = p.run(trust_exitcode=False)

        if b"seqtk <command> <arguments>" not in err:
            raise IOError("Unexpected output from Seqtk.\n"
                          "Stdout:\n{0}\n Stderr:\n{1}\n"
                          .format(out, err))

        log.info("Seqtk pre-flight checks finished successfully")

        log.info("Performing pre-flight checks on blastn")
        cmd = ("blastn", "-version")
        p = Process(cmd)
        out, err = p.run()

        if err or not out.startswith(b"blastn: "):
            raise IOError("Unexpected output from blastn.\n"
                          "Stdout:\n{0}\n Stderr:\n{1}\n"
                          .format(out, err))

        log.info("blastn pre-flight checks finished successfully")

    def run(self):
        log.info("Checking if all samples have a file with forward reads")
        for entry in self.datacfg.data:
            if entry["Contaminants"] and entry["skip_contaminants"] == False:
                log.debug("Contaminants requested for '%s'", entry["ID"])
                for infile in entry["files"]:
                    if not infile.endswith("L001_R1_001.fastq"):
                        # R2 should have the same results as R1,other lanes as well
                        # cont_file, primary_organism
                        self.data[infile] = (None, "=")
                        continue

                    self.data[infile] = None
            else:
                if entry["skip_contaminants"] == True:
                    log.info(
                        'Blast search will be skipped for sample {0} due to good alignment'.format(entry["Sample_ID"]))
                    for infile in entry["files"]:
                        if infile.endswith("L001_R1_001.fastq"):
                            self.data[infile] = ('skipped', os.path.splitext(entry["Reference"])[0])

                continue

        if not os.path.isdir(self.destdir):
            os.mkdir(self.destdir)

        for file in self.data:

            # jcostaDamage (except: at the end..daaaaa)
            try:
                # Skip R2 reads
                if self.data[file] is not None:
                    continue

                filename = os.path.basename(file)
                name, ext = os.path.splitext(filename)

                log.info("Running '%s' for taxonomic distribution analysis",
                         filename)

                outsample = os.path.join(self.destdir, name + ".sample")
                # jcostaDamage
                # SAMPLE
                sample = "1000"
                log.debug("Sampling %s sequences from '%s' and writing to '%s'",
                          sample, filename, outsample)

                cmd = ("seqtk", "sample", file, sample)
                p = Process(cmd)

                with open(outsample, 'w') as fh:
                    out, err = p.run(stdout=fh)

                outfasta = os.path.join(self.destdir, name + ".fasta")
                log.debug("Converting sample '%s' from FastQ to Fasta format and "
                          "writing to '%s'", outsample, outfasta)
                cmd = ("seqtk", "seq", "-A", outsample)
                p = Process(cmd)
                with open(outfasta, 'w') as fh:
                    out, err = p.run(stdout=fh)

                    outblast = os.path.join(self.destdir, name + ".txt")
                    log.info("BLASTing sampled sequences against nt database. "
                             "This may take a while...")
                    log.debug("Sampled file '%s' will be blasted against the nt "
                              "database located at '%s'",
                              outfasta, os.environ["BLASTDB"])

                cmd = (
                    "blastn",
                    "-task", "megablast",
                    "-outfmt", "6 sscinames sseqid",
                    "-ungapped",
                    "-max_target_seqs", "1",  # Keep only the best hit
                    "-max_hsps", "1",  # And only keep 1 HSP per hit
                    "-evalue", "1e-20",
                    "-num_threads", self.datacfg.config.params["args"].threads,
                    "-query", outfasta,
                    "-db", "nt",  # Look up at location set in BLASTDB env var
                )
                p = Process(cmd)

                with open(outblast, 'w') as fh:
                    out, err = p.run(stdout=fh, env=os.environ)

                log.info("Binning BLAST results per taxonomic group.")

                taxs = {}
                with open(outblast) as fh:
                    for line in fh:
                        # First line contains the Scientific taxonomy name
                        tax = line.rstrip().split("\t", 1)[0]
                        try:
                            taxs[tax] += 1
                        except KeyError:
                            taxs[tax] = 1

                outcounts = os.path.join(self.destdir, name + "_count.txt")
                log.debug("Writing binning result to '%s'", outcounts)

                with open(outcounts, 'w') as fh:
                    # Sort by species count in descending order
                    sortk = lambda x: x[1]
                    total_hits = 0
                    highest_hits = None
                    for sp, count in sorted(taxs.items(), key=sortk, reverse=True):
                        if highest_hits is None:
                            highest_hits = sp

                        fh.write("{0}\t{1}\n".format(count, sp))
                        total_hits += count

                    fh.write("{0}\t{1}{2}\n".format(
                        total_hits,
                        "Sequences that hit Blast nt database out of ",
                        sample))

                    # If BLAST fails to get even one hit
                    if highest_hits is None:
                        highest_hits = "No hits"
            except:
                continue

            self.data[file] = (outcounts, highest_hits)


class CoverageAndQuality(BaseTool):
    """Use Picard tools to count how many sequences align to a given reference
    """

    def __init__(self, *args, **kwargs):
        super(CoverageAndQuality, self).__init__(*args, **kwargs)
        self.input = {}
        self.destdir = os.path.join(self.datacfg.dir, "Quality")
        self.tmpdir = os.path.join(self.datacfg.dir, "Tmp")
        self.indexdir = os.path.join(self.datacfg.config.params["ref_path"],
                                     "index")
        self.tot_ref_length = 0
        self.nextseq = False

    def prepare(self):
        self.nextseq = self.datacfg.config.params['nextseq']
        log.info("Performing pre-flight checks on Bwa aligner")
        cmd = ("bwa")
        p = Process(cmd)
        out, err = p.run(trust_exitcode=False)

        if b"bwa <command> [options]" not in err:
            raise IOError("Unexpected output from Bwa.\n"
                          "Stdout:\n{0}\n Stderr:\n{1}\n"
                          .format(out, err))

        log.info("Performing pre-flight checks on SamTools")

        cmd = ("samtools",)
        p = Process(cmd)
        out, err = p.run(trust_exitcode=False)

        if b"Program: samtools" not in err:
            raise IOError("Unexpected output from Samtools.\n"
                          "Stdout:\n{0}\n Stderr:\n{1}\n"
                          .format(out, err))

        log.info("Performing pre-flight checks on qualimap")

        cmd = ("qualimap", "-h")
        p = Process(cmd)
        out, err = p.run()

        if b"QualiMap v." not in out:
            raise IOError("Unexpected output from Qualimap.\n"
                          "Stdout:\n{0}\n Stderr:\n{1}\n"
                          .format(out, err))

        log.info("Preparing output files and confirming references")
        for entry in self.datacfg.data:
            params = self.datacfg.config.params

            log.debug("Checking sample '%s' for reference", entry["Sample_ID"])
            try:
                ref = entry["Reference"]
            except KeyError:
                ref = None

            if not ref:
                log.debug("Got no reference in sample '%s'",
                          entry["Sample_ID"])
                continue

            reference = os.path.join(params["ref_path"], ref)
            index = os.path.join(self.indexdir, ref)

            if not os.path.isfile(reference):
                raise Error("Reference file '{0}' doesn't exist in '{1}'"
                            .format(ref, reference))

            # Reference length, for theoretical coverage
            self.tot_ref_length = 0
            for seq_record in SeqIO.parse(str(reference), "fasta"):
                output_line = '%s\t%i' % (seq_record.id, len(seq_record))
                self.tot_ref_length += len(seq_record)

            # Use the sample ID as folder location
            file = entry["ID"]
            fname = os.path.basename(file)
            folder = os.path.splitext(fname)[0]
            outputfolder = os.path.join(self.destdir, folder)
            samfile = os.path.join(self.tmpdir, entry["ID"] + ".sam")
            bamfile = os.path.join(self.tmpdir, entry["ID"] + ".bam")
            sampleid = entry["Sample_ID"]
            log.debug("Preparing to run %s against reference %s with index %s and length %s bp",
                      entry["files"], reference, index, self.tot_ref_length)
            log.debug("Using temporary files %s, %s and final files in %s",
                      samfile, bamfile, outputfolder)

            files = ()

            for file in entry["files"]:
                if not (file.endswith("R1_001.fastq") or file.endswith("R2_001.fastq")):
                    continue
                files = files + (file,)
            self.input[files] = (reference, index, samfile, bamfile,
                                 outputfolder, sampleid)

    def run(self):
        if not os.path.isdir(self.destdir):
            os.mkdir(self.destdir)

        if not os.path.isdir(self.indexdir):
            os.mkdir(self.indexdir)

        if not os.path.isdir(self.tmpdir):
            os.mkdir(self.tmpdir)

        for files in self.input:
            try:
                files_str = ", ".join(map(os.path.basename, files))
                if self.input[files] is None:
                    log.info("Skipping QC of %s - missing a reference",
                             files_str)
                    continue

                ref, index, sam, bam, folder, sampleid = self.input[files]
                ref_str = os.path.basename(ref)

                log.info("Running %s for quality control", files_str)

                index_file = index + ".bwt"
                if not os.path.isfile(index_file):
                    log.info("Indexing reference '%s'", ref_str)

                    cmd = ("bwa", "index", ref, "-p", index)
                    p = Process(cmd)
                    p.run()

                log.info("Aligning files %s against reference %s", files_str,
                         ref_str)

                if self.nextseq == True:
                    nextseq_files = ()
                    nextseq_files = nextseq_files + (
                    os.path.join(self.datacfg.dir, 'Merged_lanes/' + sampleid + "_R1.fastq"),)
                    if len(files) == 8:
                        nextseq_files = nextseq_files + (
                        os.path.join(self.datacfg.dir, 'Merged_lanes/' + sampleid + "_R2.fastq"),)
                    cmd = (
                              "bwa",
                              "mem",
                              "-t", self.datacfg.config.params["args"].threads,
                              index) + nextseq_files
                    p = Process(cmd)

                else:
                    cmd = (
                              "bwa",
                              "mem",
                              "-t", self.datacfg.config.params["args"].threads,
                              index) + files
                    p = Process(cmd)

                with open(sam, 'w') as fh:
                    p.run(stdout=fh)

                log.info("Converting formats SAM -> BAM and sorting")
                log.debug("SAM: '%s' -> BAM '%s'", sam, bam)
                os.popen("samtools view -@ {2} -Sb {0} > {1}".format(sam, bam.replace('.bam', '_unsorted.bam'),
                                                                     self.datacfg.config.params["args"].threads))
                os.popen("samtools sort -@ {2} -o {1} {0} ".format(bam.replace('.bam', '_unsorted.bam'), bam,
                                                                   self.datacfg.config.params["args"].threads))
                # cmd1 = ("samtools", "view", "-Sb", sam)
                # cmd2 = ("samtools", "sort", "-", os.path.splitext(bam)[0])
                # p1 = Process(cmd1)
                # p1.run()
                # p2 = Process(cmd2)
                # p2.run(stdin=p1.p.stdout)

                # Avoid zombies on the first process
                # p1.finish()

                log.debug("Removing SAM file")

                # Remove SAM file to reduce space
                os.remove(sam)
                os.remove(bam.replace('.bam', '_unsorted.bam'))

                log.info("Indexing BAM file")
                log.debug("Indexed file: '%s'", bam)

                cmd = ("samtools", "index", bam)
                p = Process(cmd)
                p.run()

                log.info("Running Qualimap on %s", files_str)
                cmd = (
                    "qualimap",
                    "bamqc",
                    "-nt", self.datacfg.config.params["args"].threads,
                    "-bam",
                    bam,
                    "-outdir", folder,
                )
                p = Process(cmd)
                # Before was: out, err = p.run()
                p.run()

                log.debug("Removing BAM files")

                # Remove BAM file and index
                for f in glob(bam + "*"):
                    os.remove(f)

                log.info("Parsing Qualimap output for number of mapped bases")

                stats = {}

                with open(os.path.join(folder, "genome_results.txt")) as fh:
                    for l in fh:
                        l = l.strip()

                        if l.startswith("number of reads ="):
                            stats["total_reads"] = l.split()[-1].replace(",", "")
                            continue

                        if l.startswith("number of mapped reads ="):
                            stats["mapped_reads"] = l.split()[-2].replace(",", "")
                            continue

                        if l.startswith("number of mapped bases ="):
                            stats["mapped_bases"] = l.split()[-2].replace(",", "")
                            continue

                        if l.startswith("mean coverageData ="):
                            stats["mean_coverage"] = l.split()[-1].replace("X", "")
                            continue

                perc_mapped_reads = (
                        int(stats["mapped_reads"]) / int(stats["total_reads"]) * 100
                )
                log.info(perc_mapped_reads)
                if perc_mapped_reads >= 80:
                    for entry in self.datacfg.data:
                        if entry["Sample_ID"] == os.path.splitext(os.path.basename(sam))[0] and entry["Contaminants"]:
                            log.info(
                                'More than 80% of the reads mapped to the given reference genome; contaminants analysis will be skipped.')
                            entry['skip_contaminants'] = True
                self.data[files[0]] = {
                    "total_mapped_reads": "{0} ({1:.2f}%)".format(
                        stats["mapped_reads"], perc_mapped_reads),
                    "total_mapped_bases": stats["mapped_bases"],
                    "total_reads_used": stats["total_reads"],
                    "mean_coverage": stats["mean_coverage"],
                    "ref_length": self.tot_ref_length,
                    "qualimap": os.path.join(folder, "qualimapReport.html"),
                }

                if len(files) == 2:
                    self.data[files[1]] = {
                        "total_mapped_reads": "=",
                        "total_mapped_bases": "=",
                        "total_reads_used": "=",
                        "mean_coverage": "=",
                        "ref_length": self.tot_ref_length,
                        "qualimap": "=",
                    }
            except:
                self.data[files[0]] = {
                    "total_mapped_reads": "ERROR",
                    "total_mapped_bases": "ERROR",
                    "total_reads_used": "ERROR",
                    "mean_coverage": "ERROR",
                    "ref_length": self.tot_ref_length,
                    "qualimap": "ERROR",
                }

                if len(files) == 2:
                    self.data[files[1]] = {
                        "total_mapped_reads": "=",
                        "total_mapped_bases": "=",
                        "total_reads_used": "=",
                        "mean_coverage": "=",
                        "ref_length": self.tot_ref_length,
                        "qualimap": "=",
                    }
                continue

        os.rmdir(self.tmpdir)


class Demultiplexer(BaseTool):
    """Use the first bases on each sequence and information contained in the
    demux_map file to split different sequences to each output.
    """

    def __init__(self, *args, **kwargs):
        super(Demultiplexer, self).__init__(*args, **kwargs)

        self.destdir = os.path.join(self.datacfg.dir, "Demultiplex")

    def _parse_demux_file(self, filename, filepath):
        data = {}

        samplepath = os.path.join(self.destdir, os.path.splitext(filename)[0])
        os.makedirs(samplepath)

        # Add 2 files for reads that don't match any tag
        rejected1 = os.path.join(samplepath, "rejected_1.fastq")
        rejected2 = os.path.join(samplepath, "rejected_2.fastq")
        data["rejected"] = (rejected1, rejected2)
        log.debug("Rejected demux files %s", data["rejected"])

        # Extra information useful during and after demuxing
        data["options"] = {}
        data["options"]["samplepath"] = samplepath
        stats = data["options"]["read_counts"] = {}
        stats["rejected"] = 0
        samples = data["options"]["samples"] = {}

        # NOTE Keep track of files for following pipeline steps (mergereads)
        files = []

        with open(filepath) as fh:
            header = None

            prev_tagsize = None

            for line_no, line in enumerate(csv.reader(fh)):
                if header is None:
                    header = line

                    expected = 7
                    cols = len(header)
                    if cols != expected:
                        raise Error("Demux file {0} has an unexpected {1} "
                                    "number of columns. Expected {2}".format(
                            filename, cols, expected))
                    continue

                sample_id = line[0]

                if not sample_id:
                    raise Error("Demux file {0} has invalid data. "
                                "Sample column cannot be empty - line {1}."
                                .format(filename, line_no + 1))

                sample1 = os.path.join(samplepath, sample_id + "_R1.fastq")
                sample2 = os.path.join(samplepath, sample_id + "_R2.fastq")
                log.debug("Demux files will be %s and %s", sample1, sample2)

                strandF, tagF, seqF, strandR, tagR, seqR = line[1:]

                # Check that indeed strands are in the correct order
                if not (strandF == "F" and strandR == "R"):
                    raise Error("Demux file {0} should have Forward "
                                "information on columns 2-4 and Reverse "
                                "information on columns 5-7".format(filename))

                if len(seqF) == len(seqR):
                    if prev_tagsize is None:
                        pass

                    elif len(seqF) != prev_tagsize:
                        log.debug("Tag previous tag had size %s while current "
                                  "has %s", prev_tagsize, len(seqF))
                        raise Error("Tags on file {0} have different lengths. "
                                    "current algorithm requires same size "
                                    "tags".format(filename))

                    prev_tagsize = len(seqF)

                else:
                    log.debug("Tag forward and tag reverse are %s and %s in "
                              "size", len(seqF), len(seqR))
                    raise Error("Forward and reverse tags on file {0} have "
                                "different lengths. Current algorithm "
                                "requires same size tags".format(filename))

                try:
                    data[seqF][seqR] = (sample1, sample2)
                    stats[seqF][seqR] = 0
                except KeyError:
                    data[seqF] = {seqR: (sample1, sample2)}
                    stats[seqF] = {seqR: 0}

                if sample_id in samples:
                    raise Error("Duplicate sample_id '{0}' on demux file {1} "
                                ". Sample_id column must be unique and "
                                "non-empty.".format(sample_id, filename))

                samples[sample_id] = (seqF, seqR)

                files.append((sample1, sample2))

            if prev_tagsize is None:
                raise Error("No tag found on file {0}".format(filename))

        log.info("Tag size to use in algorithm is %s", prev_tagsize)
        data["options"]["tagsize"] = prev_tagsize

        return files, data

    def _open_filehandles(self, demux):
        """Open files for demultiplexing. This method avoids having too many
        open files
        """
        for k1 in demux.keys():
            if k1 == "rejected":
                demux[k1] = (open(demux[k1][0], 'w'),
                             open(demux[k1][1], 'w'))

            elif k1 == "options":
                pass

            else:
                for k2 in demux[k1].keys():
                    demux[k1][k2] = (open(demux[k1][k2][0], 'w'),
                                     open(demux[k1][k2][1], 'w'))

    def _close_filehandles(self, demux):
        """Close files for demultiplexing. This method avoids having too many
        open files
        """
        for k1 in demux.keys():
            if k1 == "rejected":
                for fh in demux[k1]:
                    fh.close()
                demux[k1] = map(lambda x: x.name, demux[k1])

            elif k1 == "options":
                pass

            else:
                for k2 in demux[k1].keys():
                    for fh in demux[k1][k2]:
                        fh.close()
                    demux[k1][k2] = map(lambda x: x.name, demux[k1][k2])

    def _reject_seq(self, rejected, seq1, seq2):
        """Used to write Fastq sequences to the rejected/unclassified files
        """
        seq1.write_to_fastq_file(rejected[0])
        seq2.write_to_fastq_file(rejected[1])

    def prepare(self):
        for entry in self.datacfg.data:
            if not entry["Demux"]:
                continue

            demuxfile = entry["Demux"]
            demuxpath = os.path.join(self.datacfg.config.params["input_path"],
                                     demuxfile)
            log.debug("Using Demux at %s for sample %s", demuxpath,
                      entry["Sample_ID"])

            if not os.path.isfile(demuxpath):
                raise Error("Demux file {0} couldn't be found at {1}".format(
                    demuxfile, self.datacfg.config.params["input_path"]))

            log.info("Parsing demux file %s", demuxfile)

            for infile in entry["files"]:
                if not infile.endswith("R1_001.fastq"):
                    self.data[infile] = "="
                    continue

                demux_files, demux = self._parse_demux_file(demuxfile,
                                                            demuxpath)
                entry["demux_files"] = demux_files

                self.data[infile] = (entry["files"], demux)

    def run(self):
        for file in self.data:
            if self.data[file] == "=":
                continue

            files, demux = self.data[file]
            sample = os.path.basename(file).replace("_L001_R1_001.fastq", '')

            log.debug("Opening filehandles for demux of sample %s", sample)
            self._open_filehandles(demux)

            r1, r2 = files
            tagsize = demux["options"]["tagsize"]
            read_counts = demux["options"]["read_counts"]

            total_counts = {"classified": 0,
                            "unclassified": 0}

            with nested(open(r1), open(r2)) as (fh1, fh2):
                f1 = iter(FastqReader(fh1))
                f2 = iter(FastqReader(fh2))

                while True:
                    try:
                        seq1 = next(f1)
                        seq2 = next(f2)
                    except StopIteration:
                        break

                    key1 = seq1.seq[:tagsize]
                    key2 = seq2.seq[:tagsize]

                    match1 = get_by_distance_from_dict(demux, key1)

                    if match1 is None:
                        self._reject_seq(demux["rejected"], seq1, seq2)
                        total_counts["unclassified"] += 1
                        read_counts["rejected"] += 1
                        continue

                    match2 = get_by_distance_from_dict(demux[match1], key2)

                    if match2 is None:
                        self._reject_seq(demux["rejected"], seq1, seq2)
                        total_counts["unclassified"] += 1
                        read_counts["rejected"] += 1
                        continue

                    files = demux[match1][match2]

                    # Remove the tag sequence from final output
                    new_seq1 = seq1[tagsize:]
                    new_seq2 = seq2[tagsize:]
                    # Use the original name. [part] is appended when slicing
                    new_seq1.name = seq1.name
                    new_seq2.name = seq2.name

                    new_seq1.write_to_fastq_file(files[0])
                    new_seq2.write_to_fastq_file(files[1])

                    total_counts["classified"] += 1
                    read_counts[match1][match2] += 1

            log.debug("Results of demux are %s classified and %s unclassified",
                      total_counts["classified"], total_counts["unclassified"])
            demux["options"]["total_counts"] = total_counts

            log.debug("Closing filehandles for demux of sample %s", sample)
            self._close_filehandles(demux)


class MergeReads(BaseTool):
    """Use FLASH to fuse the paired reads into a single read.
    This process only makes sense if the dataset is 16S or if the sequenced
    fragments are expected to be smaller than ~500bp (300+300-overlap)
    """

    def __init__(self, *args, **kwargs):
        super(MergeReads, self).__init__(*args, **kwargs)
        self.nextseq = False
        self.destdir = os.path.join(self.datacfg.dir, "Merged_reads")

    def prepare(self):
        log.info("Performing pre-flight checks on Flash")
        self.nextseq = self.datacfg.config.params['nextseq']
        cmd = ("flash", "--version")
        p = Process(cmd)
        out, err = p.run(trust_exitcode=False)

        if b"FLASH v1.2.11" not in out:
            raise IOError("Unexpected output from Flash.\n"
                          "Stdout:\n{0}\n Stderr:\n{1}\n"
                          .format(out, err))

        log.info("Flash pre-flight checks finished successfully")

        for entry in self.datacfg.data:
            if not entry["Merge_reads"]:
                continue
            # if self.nextseq:
            #	pass
            # else:
            infile = None
            indexfile = None
            for file in entry["files"]:
                if file.endswith("L001_R1_001.fastq"):
                    infile = file
                elif file.endswith("I1_001.fastq"):
                    indexfile = file
                    print(indexfile)
                else:
                    self.data[file] = "="

            if not entry["Demux"]:
                if self.nextseq:
                    merged_lanes_dir = os.path.join(self.datacfg.dir, 'Merged_lanes/')
                    self.data[infile] = [tuple([merged_lanes_dir + entry["Sample_ID"] + "_R1.fastq",
                                                merged_lanes_dir + entry["Sample_ID"] + "_R2.fastq"])]
                else:
                    self.data[infile] = [entry["files"]]
                if indexfile != None:
                    self.data[infile][0] = tuple(x for x in self.data[infile][0] if x != indexfile)
            else:
                self.data[infile] = entry["demux_files"]

            print(infile, self.data[infile])

            final_files = []

            for f1, f2 in self.data[infile]:
                if not entry["Demux"] and not self.nextseq:
                    prefix = os.path.basename(f1).replace("_R1_001.fastq", "")
                else:
                    prefix = os.path.basename(f1).replace("_R1.fastq", "")
                outputdir = os.path.join(self.destdir, entry["Sample_ID"])
                final_files.append((f1, f2, outputdir, prefix))

            self.data[infile] = final_files

    def run(self):
        for file in self.data:
            if self.data[file] == "=":
                continue

            results = {}

            for f1, f2, outdir, sample_id in self.data[file]:
                if not os.path.isdir(outdir):
                    os.makedirs(outdir)

                logfile = os.path.join(outdir, sample_id + "_flash.log")

                log.info("Running Flash on sample files %s %s", f1, f2)
                cmd = (
                    "flash",
                    "-t", "1",
                    "-o", sample_id,
                    "-d", outdir,
                    "-M", "250",
                    f1, f2)
                p = Process(cmd)

                with open(logfile, 'w') as fh:
                    p.run(stdout=fh)

                log.info("Collecting FLASH output files")
                merged_reads = os.path.join(outdir,
                                            sample_id + ".extendedFrags.fastq")
                failed_f1 = os.path.join(outdir,
                                         sample_id + ".notCombined_1.fastq")
                failed_f2 = os.path.join(outdir,
                                         sample_id + ".notCombined_2.fastq")
                histogram = os.path.join(outdir, sample_id + ".histogram")
                hist = os.path.join(outdir, sample_id + ".hist")

                # Parse flash's logfile to retrieve counts
                total_reads_combined = 0
                total_reads_notcombined = 0

                with open(logfile) as fh:
                    for line in fh:
                        if "Combined pairs:" in line:
                            total_reads_combined = int(line.split()[-1])
                        elif "Uncombined pairs:" in line:
                            total_reads_notcombined = int(line.split()[-1])

                results[sample_id] = {
                    "merged_reads": merged_reads,
                    "not_combined_R1": failed_f1,
                    "not_combined_R2": failed_f2,
                    "hist": hist,
                    "histogram": histogram,
                    "logfile": logfile,
                    "reads_combined": total_reads_combined,
                    "reads_notcombined": total_reads_notcombined,
                }

            self.data[file] = results

