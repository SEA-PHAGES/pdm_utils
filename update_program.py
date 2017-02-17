#!/usr/bin/env python
#PYTHON code for updating phage.Program data from CSV
#Travis Mavrich
#Modified scripts written by Cbowman
#20170217


import MySQLdb as mdb
import sys, os, datetime, getpass

#Get the parameters
try:
	database = sys.argv[1]
	infile = sys.argv[2]
except:
	print "Incorrect Parameters - ./update_program.py DATABASE PROGRAM_CSV"
	sys.exit(1)

#GET u/n and pw
username = getpass.getpass(prompt='mySQL username:')
password = getpass.getpass(prompt='mySQL password:')


#Update genome information
try:
	
	f = open(infile,'r')
	lines= f.read().splitlines()
	f.close()
	
	#Connect
	con = mdb.connect('localhost', username, password, database)
	cur = con.cursor()

	for line in lines:
		line = line.split(',')
		#After split, there should be two fields. 0 = phageID.  1 = Program.

		print "UPDATE phage SET Program = '" + str(line[1]) + "' where PhageID = '" + str(line[0]) + "';"
		cur.execute("UPDATE phage SET Program = '" + str(line[1]) + "' where PhageID = '" + str(line[0]) + "';")		

	cur.execute("COMMIT")
	con.close()

except mdb.Error, e:
  
	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)

