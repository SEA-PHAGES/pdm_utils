#!/usr/bin/env python
#PYTHON code for executing parallel conserved domain searches
#Charles Bowman, adapted from original code from Matt Bogel
#cab106@pitt.edu
#Updated 20161026 by Travis Mavrich
#Updated 20180515 by Christian Gauthier to account for differences in cdd database basenames
#NCBI Legacy Blast toolkit no longer supported by Biopython.
#Script has been updated to use rpsblast+ from NCBI BLAST+ toolkit.
#Be sure to change rpsblast_exe to fit your system
#Change SQL server at 82, 107, 126 if yours differ from localhost

# Import standard libraries
import os
import sys
import getpass

# Import 3rd party libraries
try:
	import MySQLdb as mdb
	import pp
except ImportError:
	print "Unable to load one or more 3rd party libraries (MySQLdb or pp)."
	print "Install libraries and try again."
	sys.exit(1)
	
# Retrieve command line arguments
try:
	database = sys.argv[1]
	rpsblast_db_dri = sys.argv[2]
except IndexError:
	print """
This is a python script to identify conserved domains in phage genes.
It requires two argument(s):

	1) name of MySQL database that will be updated (e.g. 'Actino_Draft').
	2) path to conserved domain database file basename (e.g. ~/Databases/cdd_new/).
	
Note: script assumes the CDD database basename is the same as the name of the file that houses it.
Note: script assumes path to NCBI rpsblast+ executable is '/usr/bin/rpsblast+'.
"""
	sys.exit(1)

#Verify the cdd folder exists
#Expand home directory in case "~" home variable was used
rpsblast_db_dir = os.path.expanduser(rpsblast_db_dir)

#Expand the path to make sure it is a complete directory path (in case user inputted path with './path/to/folder')
rpsblast_db_dir = os.path.abspath(rpsblast_db_dir)

if os.path.isdir(rpsblast_db_dir) == False:
    print "\n\nInvalid input for CDD folder.\n\n"
    sys.exit(1)

#Add the basename for all files that constitute the cdd
rpsblast_db = os.path.join(rpsblast_db_dir,rpsblast_db_dir.split('/')[-1])

def search(geneid, translation, database, username, password, cd_db):
	#IMPORT STUFF
	import Bio

	from Bio.Blast.Applications import NcbirpsblastCommandline
	from Bio.Blast import NCBIXML
	import MySQLdb as mdb

	#DEFINE STUFF - Change variables here for executable
	rpsblast_exe = "/usr/bin/rpsblast+"
	query_filename = "/tmp/" + geneid + ".txt"
	output_filename = "/tmp/" + geneid + "_rps_out.xml"
	E_VALUE_THRESH = 0.001	#Adjust the expectation cut-off here

	#WRITE STUFF
	f = open(query_filename,'w')
	f.write(">" + geneid + "\n" + translation)
	f.close()

	#Compile the rpsblast command that will be executed.
	#outfmt. sets the format of the cdd data. 5 = XML format
	rps_command = NcbirpsblastCommandline(cmd=rpsblast_exe, db=cd_db, query= query_filename, evalue=E_VALUE_THRESH,outfmt=5,out=output_filename)
	rps_command()
	output_handle = open(output_filename,"r")

	#PARSE STUFF
	for record in NCBIXML.parse(output_handle):
		if record.alignments:
			for align in record.alignments:
				for hsp in align.hsps:
					align.hit_def = align.hit_def.replace("\"", "\'")
					con=False #initialize this variable. In case connection can't be made in the try clause, the finally clause to close con won't fail.
					try:
						descList = align.hit_def.split(',')
						if len(descList) >= 3:
							DomainID, Name = descList[0], descList[1]
							description = ','.join(descList[2:])
						elif len(descList) == 2:
							DomainID, description = descList[0], descList[1]
							Name = None
						elif len(descList) == 1:
							description = descList[0]
							DomainID, Name = None
						try: DomainID, Name, description = DomainID.strip(), Name.strip(), description.strip()
						except: pass # if DomainID, Name or description are None, strip() raises an objection

						#Connect to mysql and post hit
						con = mdb.connect('localhost', username, password, database)
						cur = con.cursor()
						sqlQuery = """insert ignore into domain (hit_id, DomainID, Name, description) VALUES ("%s", "%s", "%s", "%s")""" % (align.hit_id, DomainID, Name, description)
						cur.execute(sqlQuery)
						sqlQuery = """insert into gene_domain (geneid, hit_id, expect, query_start, query_end) VALUES ("%s", "%s", %s, %s, %s)""" % (geneid, align.hit_id, float(hsp.expect), int(hsp.query_start), int(hsp.query_end))
						cur.execute(sqlQuery)
					except mdb.Error, e:
					  	if e[0] == 1062:
					  		#print "Error %d: %s" % (e.args[0],e.args[1])
					  		print "%s. This hit will be ignored." % e.args[1]
					  	else:
					  		sys.exit(1)
					finally:
						assert hsp.expect <= E_VALUE_THRESH
						if con:
							cur.execute('COMMIT')
							con.close()


	#Now that cdd data has been added, mark the cdd_status field to 1 so that it won't be re-checked in future cdd search
	con=False #initialize this variable. In case connection can't be made in the try clause, the finally clause to close con won't fail.
	try:
		#connect to sql, post cdd_status, and commit changes in this thread
		con = mdb.connect('localhost', username, password, database)
		cur = con.cursor()
		sqlQuery = "update gene set cdd_status = 1 where geneid = '" + geneid + "'"
		cur.execute(sqlQuery)
		cur.execute('COMMIT')
	except mdb.Error, e:
		print "Error %d: %s" % (e.args[0],e.args[1])
		sys.exit(1)
	finally:
		if con:
			con.close()

#GET STUFF
print "\n\n"
username = getpass.getpass(prompt='mySQL username:')
print "\n\n"
password = getpass.getpass(prompt='mySQL password:')
print "\n\n"
con=False #initialize this variable. In case connection can't be made in the try clause, the finally clause to close con won't fail.
try:
	print "\n\n\nFetching Genes."
	con = mdb.connect('localhost', username, password, database)
	cur = con.cursor()
	cur.execute("select GeneID, translation from gene where cdd_status < 1")
	tuples = cur.fetchall()

except mdb.Error, e:

	print "Unable to connect to the database."
	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)

finally:
	if con:
		con.close()

#IF THERE IS STUFF, PROCESS STUFF
if tuples:
	#Set up pp server
	job_server = pp.Server(secret="password")
	print "pp initialized, " + `job_server.get_ncpus()` + " CPUs in use"

	#make the jobs
	jobs = []
	for tuple in tuples:
		jobs.append(job_server.submit(search, (tuple[0], tuple[1], database, username, password, rpsblast_db),(),()))
	numgenes = len(tuples)
	print "%s searches submitted...\n\n\n" % numgenes

	#WAIT FOR STUFF TO BE DONE
	counter = 0
	for job in jobs:
		counter += 1
		print "%s/%s" %(counter,numgenes)
		result = job()
else:
	print "No genes to process."

#Close script.
print "\n\n\n\nCDD script completed."
