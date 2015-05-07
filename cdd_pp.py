#!/usr/bin/env python
#PYTHON code for executing parallel conserved domain searches
#Charles Bowman, adapted from original code from Matt Bogel
#cab106@pitt.edu

#BE SURE TO CHANGE rpsblast_exe and rpsblast_db to fit your system
#Change SQL server at 82, 107, 126 if yours differ from localhost

#Add column:
#ALTER TABLE `Mycobacteriophage_Draft`.`gene` ADD COLUMN `cdd_status` TINYINT(1) NOT NULL  AFTER `blast_status` ;
#SET SQL_SAFE_UPDATES = 0;
#update gene
#set cdd_status = 0;
 #Or 1 if already completed for current DB
#SET SQL_SAFE_UPDATES = 1;
#To use on legacy databases

#To reset STUFF for update of cdd
#Truncate table gene_domain;
#SET SQL_SAFE_UPDATES = 0;
#update gene
#set cdd_status = 0;
#SET SQL_SAFE_UPDATES = 1;
#SET SQL_SAFE_UPDATES = 0;
#delete from domain;
#SET SQL_SAFE_UPDATES = 1;


import os, sys
import MySQLdb as mdb
import pp
import getpass

try:
	database = sys.argv[1]
except:
	print "Incorrect Parameters - ./cdd_pp.py DATABASE"
	print "Be sure to point rpsblast_exe and rpsblast_db to the proper directory in the file!"
	sys.exit(1)

def search(geneid, translation, database, username, password):

	#IMPORT STUFF
	import Bio
	from Bio.Blast import NCBIStandalone
	from Bio.Blast import NCBIXML
	import MySQLdb as mdb

	#DEFINE STUFF - Change variables here for executable and CDD locations
	rpsblast_exe = "/home/cbowman/Applications/BLAST/bin/rpsblast"
	rpsblast_db = "/home/cbowman/Databases/CDD/Cdd"
	query_filename = "/tmp/" + geneid + ".txt"
	E_VALUE_THRESH = 0.001	#Adjust the expectation cut-off here

	#WRITE STUFF
	f = open(query_filename,'w')
	f.write(">" + geneid + "\n" + translation)
	f.close()

	output_handle, error_handle = NCBIStandalone.rpsblast(rpsblast_exe, rpsblast_db, query_filename, expectation=E_VALUE_THRESH)
	
	#PARSE STUFF
	for record in NCBIXML.parse(output_handle):
		if record.alignments:
			for align in record.alignments:
				for hsp in align.hsps:
					align.hit_def = align.hit_def.replace("\"", "\'")
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
					  		print "Error %d: %s" % (e.args[0],e.args[1])
					  		print "Ignoring duplicate hit"
					  	else:
					  		sys.exit(1)
					    
					finally:    
		
						assert hsp.expect <= E_VALUE_THRESH
						
						if con:    
							cur.execute('COMMIT')
							con.close()

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
print "Be sure to point rpsblast_exe and rpsblast_db to the proper directory in the file!"
username = getpass.getpass(prompt='mySQL username:')
password = getpass.getpass(prompt='mySQL password:')

try:
	print "Fetching Genes"
	con = mdb.connect('localhost', username, password, database)
	cur = con.cursor()
	cur.execute("select GeneID, translation from gene where cdd_status < 1")
	tuples = cur.fetchall()

except mdb.Error, e:
  
	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)
    
finally:    
	if con:    
		con.close()

#IF THERE IS STUFF, PROCESS STUFF
if tuples:
	#Set up pp server
	job_server = pp.Server(secret="butt")
	print "pp initialized, " + `job_server.get_ncpus()` + " CPUs in use"
	
	#make the jobs
	jobs = []
	for tuple in tuples:
		jobs.append(job_server.submit(search, (tuple[0], tuple[1], database, username, password),(),()))
	print "Searches Submitted, waiting for jobs to complete"
	numgenes = len(tuples)
	counter = 0
	
	print `numgenes` + " searches to perform, please be patient."
	
	#WAIT FOR STUFF TO BE DONE
	for job in jobs:
		counter = counter + 1
		print `counter` + " / " + `numgenes`
		result = job()
	print "done..."
else:
	print "No genes to process..."	
