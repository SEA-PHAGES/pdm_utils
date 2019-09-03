#!/usr/bin/env python
#PYTHON code for updating cluster from CSV
#Cbowman


import MySQLdb as mdb
import sys, os, random, colorsys, datetime

#Get the parameters
try:
	database = sys.argv[1]
	infile = sys.argv[2]
except:
	print "Incorrect Parameters - ./update_cluster.py DATABASE CLUSTERCSV"
	sys.exit(1)

#Get the gene information
try:

	f = open(infile,'r')
	lines= f.read().splitlines()
	f.close()

	#Connect
	con = mdb.connect('localhost', 'root', 'phage', database)
	cur = con.cursor()

	for line in lines:
		line = line.split(',')
		print "update phage set cluster = '" + str(line[1]) + "' where PhageID = '" + str(line[0]) + "';"
		cur.execute("update phage set cluster = '" + str(line[1]) + "' where PhageID = '" + str(line[0]) + "';")
	cur.execute("COMMIT")
	con.close()

except mdb.Error, e:

	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)
