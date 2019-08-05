# Import standard python libraries.
import argparse
import sys
import datetime

# Import third-party python modules.
try:
	import pymysql as pms
except ImportError:
	print("Failed to import module. Please make sure pymysql is installed, "
		  "and try again")
	sys.exit(1)

# Import functions from other k_phamerate scripts.
from classes.MySQLConnectionHandler import MySQLConnectionHandler
from classes.MMseqsPhameratorHandler import MMseqsPhameratorHandler

# Set up argparse to interact with users at the command line interface.
script_description = """
This is a script intended to quickly, sensitively, and accurately sort all 
proteins in a Phamerator database into phams using MMseqs2. (Steinegger & 
Soding, 2017) 
"""
parser = argparse.ArgumentParser(description=script_description)
parser.add_argument("database_name", metavar="db", type=str, nargs=1,
					help="name of the Phamerator database to be phamerated")

# Parse command line arguments.
args = parser.parse_args()
database = args.database_name[0].split("=")[-1]

# Establish the database connection using the MySQLConnectionHandler object
mysql_handler = MySQLConnectionHandler()
mysql_handler.database = database
mysql_handler.open_connection()

# Script will only continue if connection was successfully established
if mysql_handler.connection_status() is True:
	pass
else:
	sys.exit(1)

# Record pham analysis start time.
start = datetime.datetime.now()

# Create and set up MMseqsPhameratorHandler object
pham_handler = MMseqsPhameratorHandler(mysql_handler)
pham_handler.filter_redundant = True

# Set MMseqs2 clustering parameters
pham_handler.threads = 1
pham_handler.verbosity = 0
pham_handler.single_step = True
pham_handler.max_seqs = 1000
pham_handler.min_seq_id = 0.4
pham_handler.coverage = 0.8
pham_handler.alignment_mode = 3
pham_handler.coverage_mode = 0
pham_handler.cluster_mode = 0

# Refresh the temp directory.  Default action is to delete the directory if
# it exists and then recreate it, because MMseqs2 will use files from
# previous runs if they still exist.  Default temp directory is /tmp/MMseqs2/
pham_handler.refresh_temp_dir()

# Get current pham data from the database.
pham_handler.get_existing_pham_data()

# Get unphamerated GeneIDs.
pham_handler.get_unphamerated_geneids()

# Write all GeneIDs and translations to fasta file or just the non-redundant
# translations if mysql_handler.filter_redundant is True
pham_handler.write_geneids_and_translations_to_fasta()

# Convert the fasta file to mmseqsdb internal format.
pham_handler.convert_fasta_to_mmseqsdb()

# Cluster the database using the indicated parameters.
pham_handler.cluster_database()

# Convert and parse the output back into a python dictionary object.
pham_handler.convert_and_parse_output()

# Put all GeneIDs back together again, if duplicates were filtered out.
# Nothing will happen at this stage if duplicates were not filtered out.
pham_handler.reintroduce_duplicates()

# Preserve pham names by finding best matches from the old phams.
pham_handler.conserve_pham_names()

# Put data back into the database.
pham_handler.reinsert_pham_data()

# Fix miscolored phams and orphams.  Phams are miscolored if they're white,
# orphams are miscolored if they are not white.
pham_handler.fix_miscolored_phams()

# Record end time
end = datetime.datetime.now()

# Print start and end times so the user knows how long it took.
print("Start: {}".format(start))
print("End:   {}".format(end))
