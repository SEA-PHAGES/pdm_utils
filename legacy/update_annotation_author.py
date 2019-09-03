#!/usr/bin/env python
#PYTHON code for updating the AnnotationAuthor field from CSV
#Travis Mavrich
#Modified scripts written by Cbowman


import MySQLdb as mdb
import sys

#Get the parameters
try:
	database = sys.argv[1]
	infile = sys.argv[2]
except:
	print "Incorrect Parameters - ./update_annotation_author.py DATABASE ANNOTATION_AUTHOR_CSV"
	sys.exit(1)

try:

	#Get the author information
	f = open(infile,'r')
	lines= f.read().splitlines()
	f.close()

	#Connect and update data
	con = mdb.connect('localhost', 'root', 'phage', database)
	cur = con.cursor()

	for line in lines:
		line = line.split(',')
		print "update phage set AnnotationAuthor = '" + str(line[1]) + "' where PhageID = '" + str(line[0]) + "';"
		cur.execute("update phage set AnnotationAuthor = '" + str(line[1]) + "' where PhageID = '" + str(line[0]) + "';")

	cur.execute("COMMIT")
	con.close()

except mdb.Error, e:

	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)
