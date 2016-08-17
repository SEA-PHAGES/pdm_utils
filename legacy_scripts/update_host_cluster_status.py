#!/usr/bin/env python
#PYTHON code for updating cluster, host, and status data from CSV
#Travis Mavrich
#Modified scripts written by Cbowman


import MySQLdb as mdb
import sys, os, random, colorsys, datetime, getpass

#Get the parameters
try:
	database = sys.argv[1]
	infile = sys.argv[2]
except:
	print "Incorrect Parameters - ./update_host_cluster_status.py DATABASE HOST_CLUSTER_STATUS_CSV"
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
		#After split, there should be four fields. 0 = phage_name.  1 = host. 2 = cluster. 3 = status.
		

		#If the genome status is "draft", then add this suffix to the phage name
		if str(line[3]).lower() == "draft":
			line[0] = str(line[0]) + "_Draft"


		#Host info exported from phagesdb contains genus, species, and strain info. Only the genus should be retained.
		#The default delimiter for split() is whitespace.
		hostinfo = str(line[1]).split()
		print "update phage set HostStrain = '" + hostinfo[0] + "' where name = '" + str(line[0]) + "';"
		cur.execute("update phage set HostStrain = '" + hostinfo[0] + "' where name = '" + str(line[0]) + "';")
		
		#If the phage is a singleton, Cluster should be set to NULL.
		if str(line[2]).lower() == "singleton":
			print "update phage set cluster = NULL" + " where name = '" + str(line[0]) + "';"
			cur.execute("update phage set cluster = NULL" + " where name = '" + str(line[0]) + "';")			
        else:
			print "update phage set cluster = '" + str(line[2]) + "' where name = '" + str(line[0]) + "';"
			cur.execute("update phage set cluster = '" + str(line[2]) + "' where name = '" + str(line[0]) + "';")


		print "update phage set status = '" + str(line[3]) + "' where name = '" + str(line[0]) + "';"
		cur.execute("update phage set status = '" + str(line[3]) + "' where name = '" + str(line[0]) + "';")

	cur.execute("COMMIT")
	con.close()

except mdb.Error, e:
  
	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)

