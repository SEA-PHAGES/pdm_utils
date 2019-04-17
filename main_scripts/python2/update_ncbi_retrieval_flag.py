#!/usr/bin/env python
#PYTHON code for updating the NCBI retrieval flag from CSV
#Travis Mavrich
#Modified scripts written by Cbowman


import MySQLdb as mdb
import sys

#Get the parameters
try:
	database = sys.argv[1]
	infile = sys.argv[2]
except:
	print "Incorrect Parameters - ./update_ncbi_retrieval_flag.py DATABASE RETRIEVAL_FLAG_CSV"
	sys.exit(1)

try:

	#Get the ncbi retrieval information
	f = open(infile,'r')
	lines= f.read().splitlines()
	f.close()

	#Connect and update data
	con = mdb.connect('localhost', 'root', 'phage', database)
	cur = con.cursor()

	for line in lines:
		line = line.split(',')
		print "update phage set RetrieveRecord = '" + str(line[1]) + "' where PhageID = '" + str(line[0]) + "';"
		cur.execute("update phage set RetrieveRecord = '" + str(line[1]) + "' where PhageID = '" + str(line[0]) + "';")

	cur.execute("COMMIT")
	con.close()

except mdb.Error, e:

	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)
