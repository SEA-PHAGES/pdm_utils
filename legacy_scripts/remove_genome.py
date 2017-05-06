#!/usr/bin/env python
#PYTHON code for deleting a list of phage genomes
#Travis Mavrich
#20160415
#Based on code previously written by C. Bowman.
#If a phage name is provided that is not present as a PhageID, nothing will be removed and no error will be thrown.


import MySQLdb as mdb
import sys, os, getpass

#Get the parameters
try:
	database = sys.argv[1]
	infile = sys.argv[2]
except:
	print "Incorrect Parameters - ./remove_genome.py DATABASE GENOME.CSV"
	sys.exit(1)

#Get username and password
username = getpass.getpass(prompt='mySQL username:')
password = getpass.getpass(prompt='mySQL password:')


#Iterate through a list of PhageIDs and delete them from the database
try:
	
	f = open(infile,'r')
	lines= f.read().splitlines()
	f.close()
	
	#Connect
	con = mdb.connect('localhost', username, password, database)
	cur = con.cursor()

	for line in lines:
		line = line.split(',')
		print "delete from phage where PhageID = '" + str(line[0]) + "';"
		cur.execute("delete from phage where PhageID = '" + str(line[0]) + "';")
	cur.execute("COMMIT")
	con.close()

except mdb.Error, e:
  
	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)

