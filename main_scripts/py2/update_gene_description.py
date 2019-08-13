#!/usr/bin/env python
#PYTHON code for updating gene description data from CSV
#Travis Mavrich
#Modified scripts written by Cbowman


import MySQLdb as mdb
import sys, os

#GET u/n and pw
username = getpass.getpass(prompt='mySQL username:')
password = getpass.getpass(prompt='mySQL password:')


#Get the parameters
try:
	database = sys.argv[1]
	infile = sys.argv[2]
except:
	print "Incorrect Parameters - ./update_description.py DATABASE GENE_DESCRIPTION_CSV"
	sys.exit(1)

#Get the gene information
try:

	f = open(infile,'r')
	lines= f.read().splitlines()
	f.close()

	#Connect
	con = mdb.connect('localhost', username, password, database)
	cur = con.cursor()

	for line in lines:
		line = line.split(',')
		print "update gene set Notes = '" + str(line[1]) + "' where GeneID = '" + str(line[0]) + "';"
		cur.execute("update gene set Notes = '" + str(line[1]) + "' where GeneID = '" + str(line[0]) + "';")

	cur.execute("COMMIT")
	con.close()

except mdb.Error, e:

	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)
