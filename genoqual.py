
from __future__ import print_function, division
import os
import sys
import logging
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

from externaltools import FastQC, QualityStats, CoverageAndQuality, TaxonomyClassifier, MultiQC, Demultiplexer, MergeReads, Qiime, QiimeOnControls
from dataops import DataWrapper, DataCollector


def hamming(str1, str2):
	"""Compute hamming distance between 2 strings

	"">>> hamming('abc', 'abc')""
	0
	">>> hamming('abc', 'abcd')""
	1
	"">>> hamming('abc', 'cbc')""
	1
	">>> hamming('abc', 'bca')""
	3
	"""
	return sum(map(str.__ne__, str1, str2))


def get_by_distance_from_dict(d, key, size_matters=True, distance=hamming,
                              threshold=1):
	"""Return all values of given dictionary that have a smaller or equal
	distance to the given threshold.

	If size_matters=True keys must have the same size or they will be treated
	as above threshold.

	Used for demultiplexing
	"""
	matches = []
	size = len(key)

	for k in d.iterkeys():
		# If key size is different, assume no possible match
		if size_matters and len(k) != size:
			continue

		if distance(k, key) <= threshold:
			matches.append(k)

	if len(matches) > 1:
		raise Error("Read with sequence {0} ambiguously matched more than "
		            "one tag {1}".format(key, matches))

	if matches:
		return matches[0]
	else:
		return None


class Error(Exception):
	pass


class ConfigParser:
	"""Reads the Config format (which is basically a CSV)
	"""
	def __init__(self, params):
		self.params = params
		self.configfile = self.params["config"]
		self.config = None

		# Ensure characters are valid to avoid problems afterwards
		self.transtable = str.maketrans(" .", "--")
		self.validchars = string.ascii_letters + string.digits + " .-_"

		# Columns expected to be present in the config file
		self._expected_columns = self.params["columns"]
		# Boolean columns (columns that should map to True/False
		self._booleans = self.params["booleans"]

		self.parse()

	def _booleanify(self, val):
		"""Try to convert given value to a boolean in a somewhat extensive way

		It will try to convert Yes/No 0/1 Y/N to a boolean. Other values will
		raise an exception.
		"""
		bools = {
		        "yes": True,
		        "no": False,
		        "y": True,
		        "n": False,
		        "1": True,
		        "0": False,
		        '': False,
		}

		try:
			return bools[val.lower()]
		except KeyError:
			raise Error("Value '{0}' is not a valid option for that column, "
			            "should be one of {1}".format(val, bools.keys()))

	def _sanitize(self, filename):
		fname = filename.translate(self.transtable)

		for char in fname:
			if char not in self.validchars:
				raise Error("Config file contains invalid characters in "
				            "column 'Sample ID'. Offending character is "
				            "{0!r} in '{1}'".format(char, filename))

		return fname

	def parse(self):
		"""Read all settings from given configfile
		"""
		if not os.path.isdir(self.params["input_path"]):
			raise Error("Folder '{0}' not found, did you specify the prefix "
			            "correctly?".format(self.params["input_path"]))

		if not os.path.isfile(self.configfile):
			raise Error("Config file not found at '{0}'"
			            .format(self.configfile))

		# Prevent parsing the file more than once if parse is called twice
		if self.config:
			log.info("Using cached config settings")
			return

		log.info("Parsing config file '%s'", os.path.basename(self.configfile))
		with open(self.configfile) as fh:
			for line in csv.reader(fh):
				# Ignore empty lines
				if not line:
					continue

				if self.config is None:
					self.config = []
					self._ids = line

					# Check that required columns are present in the file
					for col in self._expected_columns:
						if col not in self._ids:
							raise Error("Column '{0}' is missing from "
							            "'{0}'".format(col, self.configfile))

					continue

				data = {}

				for i, value in enumerate(zip(self._ids, line)):
					id, val = value
					data[id] = val

				# Keep an ID safe for use as filename prefix
				data["ID"] = self._sanitize(data["Sample_ID"])
				
				# Figure out what should be the filename from sample data
				data["prefix"] = "{0}_".format(data["ID"])

				# Create a variable that will control the execution of blastn for contaminants when necessary 
				data["skip_contaminants"] = False
				
				# Transform requested columns into a boolean
				for key, value in data.items():
					if key in self._booleans:
						data[key] = self._booleanify(value)

				log.debug("Parsed config line as %s", data)
				self.config.append(data)

		if self.config is None:
			raise Error("Provided config file is empty")


class Pipeline:
	"""Sequentially execute the steps of the pipeline, performing checks first
	and preparing all steps and only afterwards executing each component
	"""
	def __init__(self, config):
		self.config = config

		self._steps = []
		self.step("data", DataWrapper(self.config))		
		self.step("qc", FastQC(self.data))
		self.step("qual", QualityStats(self.data))
		self.step("cov", CoverageAndQuality(self.data))
		self.step("tax", TaxonomyClassifier(self.data))
		self.step("multiqc", MultiQC(self.data))
		self.step("dmux", Demultiplexer(self.data))
		self.step("merge", MergeReads(self.data))
		if self.config.params["args"].meta:
			self.step("meta", Qiime(self.data))
			self.step("metaControls", QiimeOnControls(self.data))

	def step(self, name, obj):
		setattr(self, name, obj)
		self._steps.append(obj)

	def prepare(self):    # This function is used to execute the prepare function in each of the appended tools and datawraper
		try:
			for step in self._steps:
				step.prepare()
		except:
			log.fatal("Cleaning up after fatal error during prepare phase...")
			self.cleanup()
			raise

	def run(self):     # This function is used to execute the run function in each of the appended tools and datawraper
		try:
			for step in self._steps:
				try:
					step.run()
				except:
					continue
		except:
			log.fatal("Cleaning up after fatal error during run phase...")
			self.cleanup()
			raise

	def cleanup(self):
		log.debug("Cleaning up Pipeline leftovers")
		self.data.cleanup()

	def finish(self):
		"""Export the results as HTML and move all result files to their final
		destination.
		"""
		log.debug("Removing .fastq files used in the analysis")
		for file in glob(os.path.join(self.data.dir, "*.fastq")):
			os.remove(file)
		
		if self.config.params["nextseq"] == True:
			log.debug("Removing .fastq files for NextSeq merged lanes")			
			shutil.rmtree(os.path.join(self.data.dir, "Merged_lanes"))
			
		log.debug("Preparing data for reporting")
		data = DataCollector(self)

		log.debug("Generating Tabular report")
		out = DataFormatter(data,self.config.params)
		out.write_tabular(os.path.join(self.data.dir,
		                               self.config.params["report_tab"]))

		#log.debug("Generating HTML report")
		#out.write_html(os.path.join(self.data.dir,
		#                            self.config.params["report_html"]))

		log.info("Moving all files to final destination")
		shutil.move(self.data.dir, self.config.params["output_path"])

		sleep(5)
		log.debug("Updating runs and samples db")
		#Runs the build_db.py module to build
		os.system("%s/build_db.py %s"%(os.path.dirname(os.path.realpath(__file__)), os.path.join(self.config.params["base_path"], self.config.params["output_folder"])))

		log.debug("Creating images for quality history")
		os.system("%s/graphs.py -o %s/tmp_images -d %s/runs_samples.db -r %s"%(os.path.dirname(os.path.realpath(__file__)), self.config.params["output_path"],
		                                                                       os.path.join(self.config.params["base_path"], self.config.params["output_folder"]), args.run_folder))

		log.debug("Generating HTML report")
		log.info('Joining paths: %s %s'%(self.config.params["output_path"], self.config.params["report_html"] ))
		out.write_html(os.path.join(self.config.params["output_path"],
		                            self.config.params["report_html"]))


		log.debug("Writing report to requested output location")
		extra_output = self.config.params["args"].output
		if extra_output:
			shutil.copy(self.config.params["output_file_html"], extra_output)

		log.debug("Correcting AFS permissions for POSIX strict clients")
		for root, _, files in os.walk(self.config.params["output_path"]):
			# Folder permissions
			os.chmod(root, 0o777)
			# File permissions
			for file in files:
				os.chmod(os.path.join(root, file), 0o666)

		log.info("Pipeline finished successfully")
		log.info("Results can be found at '%s'",
		         self.config.params["output_path"])


def parse_args():
	"""Initial operation arguments
	"""
	parser = argparse.ArgumentParser(description="Perform Quality Control on"
	                                 " Illumina reads")
	parser.add_argument('run_folder', #run_folder should go to the input folder
	                    help="Folder suffix where files with reads are "
	                    "located. E.g. 16 if folder name is Run_16")
	parser.add_argument('--base_path',
	                    help="Location of input/reference/results folders",
	                    default="/media/genomics/genomics2/Analysis")
	                    
	parser.add_argument('--output', '-o',
	                    help="Also store a copy of the HTML report in given "
	                    "location")
	parser.add_argument('--meta', '-m',
	                    help="Perform basic QIIME run",
	                    action='store_true')
	parser.add_argument('--blastdb',
	                    help="Location of blast databases (to be set as the "
	                    "env variable BLASTDB",
	                    default="/media/genomics/genomics2/BLAST/db")
	parser.add_argument('--url', '-u',
	                    help="Url to use in HTML output",
	                    default="geu-galaxy.igc.gulbenkian.pt/static/genoqual_results/") 
	parser.add_argument('--threads', '-t',
	                    help="Number of threads to use in subprocesses",
	                    default=os.environ.get("GALAXY_SLOTS", "8"))

	parser.add_argument('--verbose', '-v', action="count",
	                    help="Verbosity level. -vvv is the highest level")
	return parser.parse_args()


def prepare_settings(args):
	"""Creates a dictionary to organize and prepare the settings so that they can be read by the configparser
	"""
	d = {}

	d["base_path"] = args.base_path
	d["config_name"] = "QCconf.csv"
	d["columns"] = ("Sample_ID", "Description", "Contaminants", "Reference",
	                "Demux", "Merge_reads")
	d["booleans"] = {"Contaminants", "Merge_reads"}

	d["input_folder"] = "input"
	d["output_folder"] = "results"
	d["meta"] = args.meta
	d["nextseq"] = False
	d["ref_folder"] = "references"
	d["report_html"] = "report.html"
	d["report_tab"] = "report.txt"
	d["run_folder"] = f"Run_{args.run_folder}"

	d["input_path"] = os.path.join(d["base_path"], d["input_folder"],
	                               d["run_folder"])
	d["output_path"] = os.path.join(d["base_path"], d["output_folder"],
	                                d["run_folder"])
	d["output_url"] = (
	        "http://" +
	        urllib.parse.quote(urllib.parse.urljoin(args.url + d["run_folder"], ''), ':/')
#        urllib.quote(urlparse.urljoin(args.url, d["run_folder"]))
	)
	d["ref_path"] = os.path.join(d["base_path"], d["ref_folder"])
	d["config"] = os.path.join(d["input_path"], d["config_name"])
	if args.meta:
		d["metadata"] = os.path.join(d["input_path"], "metadata.csv")
	d["output_file_html"] = os.path.join(d["output_path"], d["report_html"])
	d["output_file_tab"] = os.path.join(d["output_path"], d["report_tab"])
	d["args"] = args

	# Also prepare environment settings that are required by some tools
	os.environ["BLASTDB"] = args.blastdb

	log.debug("Prepared settings are %s", d)

	return d


class DataFormatter:
    """Used to generate the HTML report that gets produced at the end of all
    steps
    """

    def __init__(self, data, config):
        self.data = data.output
        self.config = config
        self.stats_all_lanes = data.stats_all_lanes
        self.links = data.links
        self.demux = data.demux

    def write_tabular(self, filename):
        """Write a tab delimited version of the data
        """
        with open(filename, 'w') as fh:
            for line in self.data:
                fh.write("\t".join(line) + "\n")

        for filename in self.demux:
            path, values = self.demux[filename]
            file = os.path.splitext(path)[0] + ".txt"

            with open(file, 'w') as fh:
                for line in values:
                    fh.write("\t".join(line) + "\n")

    def encode_figures(self, output_path):
        encoded_figs = OrderedDict()
        for img in sorted(glob('%s/tmp_images/*.png' % output_path)):
            with open(img, "rb") as f:
                encoded_figs[os.path.basename(img)] = f.read().encode("base64")

        shutil.rmtree('%s/tmp_images/' % output_path)

        return encoded_figs

    def write_html(self, filename):
        """Write an html version of the data
        """
        self._html_file(filename, self.data, "Quality Control", self.links)

        for filename in self.demux:
            file, values = self.demux[filename]
            self._html_file(file, values, "Demultiplexing and merging", ())

    def _html_file(self, filename, values, title, links):

        encoded_figs = self.encode_figures(os.path.dirname(filename))
        row_color = cycle((
            "style='background-color:white;'",
            "style='background-color:#eefafc;'",
        ))

        i = datetime.datetime.now()
        with open(filename, 'w') as fh:
            fh.write("<!DOCTYPE html>\n")
            fh.write("<html>\n<head>\n<title>{0}</title>\n".format(title))
            fh.write(
                "<link rel='stylesheet' href='https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css' integrity='sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u' crossorigin='anonymous'>\n")
            fh.write("<script type='text/javascript' src='http://code.jquery.com/jquery-1.9.1.js'></script>\n")
            fh.write("<script>\n")
            fh.write("$(document).ready(function(){\n")
            fh.write("\t$('tbody > tr:not(\".header\")').hide();\n")
            fh.write("\t$('.header').click(function(){\n")
            fh.write("\t\t$(this).toggleClass('expand').nextUntil('tr.header').slideToggle(100);\n")
            fh.write("\t});\n")
            fh.write("})\n")
            fh.write("</script>\n")
            fh.write("<style>\n")
            fh.write("body {margin: 0; padding: 20px;}\n")
            fh.write("table {width: 100%; border-collapse: collapse;}\n")
            fh.write("table, th, td {border: 1px solid #cecece; vertical-align:middle !important}\n")
            fh.write("th {text-align:center; font-size:13px; vertical-align:middle}\n")
            fh.write("td {vertical-align: middle; text-align: center; font-size:12px}\n")
            fh.write("tr.header.expand{cursor:pointer;}\n")
            fh.write(
                ".header .sign:after{content: ''; background-image:url(http://icons.iconarchive.com/icons/icojam/blue-bits/16/math-minus-icon.png); background-size:10px 10px; display:inline-block; width:10px; height:10px; vertical-align:middle;}\n")
            fh.write(
                ".header.expand .sign:after{content:''; background-image:url(http://icons.iconarchive.com/icons/icojam/blue-bits/16/math-add-icon.png); background-size:10px 10px; display:inline-block; width:10px; height:10px;vertical-align:middle; }\n")
            fh.write("</style>\n")
            fh.write("</head>\n<body>\n")
            fh.write("<h2><b>GenoQual Report </b><small>v.1.3</small></h2>\n")
            fh.write("<p>Report generated on %s/%s/%s, %s:%s for run %s</p><br>" % (
            i.day, i.month, i.year, i.hour, i.minute, args.run_folder))
            fh.write("<div class='table-responsive'>\n")
            fh.write("<table class='table table-striped'>\n")
            for i, line in enumerate(values):
                if i == 0:
                    fh.write("<tr class='header' style='background-color:#006d7a;color: white;'>")
                elif i == len(values) - 1:
                    fh.write("<tr class='header' style='background-color:#fff7ea;'>")

                if self.config["nextseq"] == True:

                    if '_L001_R1_' in line[
                        1]:  # Before L001_R1 of a sample, we write the summary line for all the sample's lanes
                        sample = line[0].split('_')[0]
                        fh.write("<tr class='header expand' {0}>".format(next(row_color)))
                        fh.write("<td><span style='vertical-align:middle'>{0} </span><span class=\"sign\"></td>".format(
                            sample))
                        fh.write("<td></td>")  # empty field for filename
                        self.stats_all_lanes[sample]['percent_q30'] = round(
                            100 * float(self.stats_all_lanes[sample]['q30']) / self.stats_all_lanes[sample][
                                'total_bases'], 2)
                        for e in self.stats_all_lanes[sample]:
                            try:
                                value = int(self.stats_all_lanes[sample][e])
                            except ValueError:
                                fh.write("<td>{0}</td>".format(self.stats_all_lanes[sample][e]))
                            else:
                                fh.write("<td>{0:,}</td>".format(value))
                        fh.write("</tr>\n")
                        # Clear cov and tax fields for L001_R1, as they are now in the summary line.
                        line[8] = line[9] = '='
                        for n in range(10, 17):
                            line[n] = ''
                        for n in range(17, len(line)):
                            if line[n] != '':
                                line[n] = '='

                    if i != 0 and i != len(values) - 1:
                        fh.write("<tr {0}>".format(next(row_color)))
                else:
                    if i != 0 and i != len(values) - 1:
                        fh.write("<tr class='header' {0}>".format(next(row_color)))

                for j, elem in enumerate(line):
                    if i == 0:
                        fh.write("<th>")
                    else:
                        fh.write("<td>")
                    if elem == "Total":
                        fh.write('<b>')
                    if elem and not elem in ["=", "skipped"] and j in links and i != 0:
                        name = os.path.basename(elem)
                        if i == len(values) - 1:
                            fh.write("<a href='{0}'>MultiQC report</a>".format(elem))
                        else:
                            fh.write("<a href='{0}'>Report</a>".format(elem))

                    else:
                        try:
                            elem = int(elem)
                        except ValueError:
                            fh.write("{0}".format(elem))
                        else:
                            # Use a comma as thousands delimiter
                            fh.write("{0:,}".format(elem))
                    if i == 0:
                        fh.write("</th>")
                    else:
                        if elem == "Total":
                            fh.write('</b>')
                        fh.write("</td>")
                fh.write("</tr>\n")
            fh.write("</table>\n</div>\n")
            fh.write("<div class='table-responsive'>\n")
            fh.write("<table>\n")
            fh.write("<tr class='header'>")
            for c, img in enumerate(encoded_figs.keys()):
                if c < 4 or c == 7:
                    if (c % 2 == 0 and c != 0) or c == 7:
                        fh.write("</tr>\n<tr class='header'>")
                    fh.write('<td><img src="data:image/png;base64,%s" class="img-rounded" width="80%%"></td>\n' %
                             encoded_figs[img])

                    if c == 7:
                        fh.write(
                            '<td><p style="font-size:18px; font-weight:bold;padding-top: 5px;margin-bottom: 15px;">Historical unweighted unifrac PCoA plot for controls</p>\n'
                            '<iframe src="%s/Qiime/Controls_analysis/bdiv/unweighted_unifrac_emperor_pcoa_plot/index.html" style="margin:0;width:100%%;height:450px">'
                            'Alternative text for browsers that do not understand IFrames.</iframe></td></tr>' %
                            self.config["output_url"])
                if c >= 4 and c < 7:
                    fh.write("<tr class='header'>\n<td colspan=2>\n")
                    fh.write(
                        '<p style="text-align:center;"><img src="data:image/png;base64,%s" class="img-rounded" width="95%%"></p>\n' %
                        encoded_figs[img])
                    fh.write("</td></tr>\n")

            fh.write("</table>\n</div>\n")
            fh.write("</body>\n</html>\n")


def main(args):
	# Setting things up acording to the parser arguments
	settings = prepare_settings(args)
	config = ConfigParser(settings)
	# Pre-flight check on all steps
	pp = Pipeline(config)
	pp.prepare()
	# Actually compute stuff
	pp.run()
	# Generate HTML report with links to files and copy files to final location
	pp.finish()


if __name__ == "__main__":
	args = parse_args()
	if args.verbose == 1:
		level = logging.WARN
	elif args.verbose == 2:
		level = logging.INFO
	elif args.verbose >= 3:
		level = logging.DEBUG
	else:
		level = logging.ERROR

	logging.basicConfig(
	        format="%(asctime)s - %(levelname)s - %(message)s",
	        level=level,
	)
	log = logging.getLogger(__name__)


	try:
		main(args)
	except Error as e:
		log.fatal(e)
		sys.exit(1)
	except Exception as e:
		log.exception(e)
		sys.exit(1)
	

# vim: ai sts=4 et sw=4
