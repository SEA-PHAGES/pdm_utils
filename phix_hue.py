#!/usr/bin/env python
#PYTHON code for fixing phamily coloring errors
#Charles Bowman
#20150210
 

import MySQLdb as mdb
import sys, os, random, colorsys, datetime, getpass

#Get the command line parameters
try:
	database = sys.argv[1]
except:
	print "Incorrect Parameters - no database given"
	sys.exit(1)

#GET u/n and pw

username = getpass.getpass(prompt='mySQL username:')
password = getpass.getpass(prompt='mySQL password:')

#fix miscolored phams
try:
	#Connect
	con = mdb.connect('localhost', username, password, database)
	cur = con.cursor()
	print "Connected to Database - Phetching Phalsely Hued Phams"
	print "...\n"

	cur.execute("Select * from (select b.id, count(GeneID) as count, a.name, b.color from pham as a inner join pham_color as b on a.name = b.name group by a.name) as c where color = '#FFFFFF' and count > 1")
	tuples = cur.fetchall()

	print `len(tuples)` + " miscolored phamilies found."
	if len(tuples) > 0:
		print " phixing that phor you..."

	for tuple in tuples:
		pham_id, count, name, color = tuple

		h=s=v=0
		while h <= 0: h = random.random()
		while s < 0.5: s = random.random()
		while v < 0.8: v = random.random()
		rgb = colorsys.hsv_to_rgb(h,s,v)
		rgb = (rgb[0]*255,rgb[1]*255,rgb[2]*255)
		hexrgb = '#%02x%02x%02x' % rgb
		newColor = hexrgb

		cur.execute("update pham_color set color = %s where id = %s" , (newColor, pham_id))

except mdb.Error, e:
  
	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)
    
finally:
        
    cur.execute("COMMIT")
    if con:    
        con.close()

#build miscolored orphams
try:
	#Connect
	con = mdb.connect('localhost', username, password, database)
	cur = con.cursor()
	print "\nConnected to Database - Phetching Phalsely Phlagged Orphams"
	print "...\n"

	cur.execute("Select * from (select b.id, count(GeneID) as count, a.name, b.color from pham as a inner join pham_color as b on a.name = b.name group by a.name) as c where color != '#FFFFFF' and count = 1")
	tuples = cur.fetchall()

	print `len(tuples)` + " miscolored orphams found."
	if len(tuples) > 0:
		print " phixing that phor you..."

	for tuple in tuples:
		pham_id, count, name, color = tuple
		newColor = "#FFFFFF"

		cur.execute("update pham_color set color = %s where id = %s" , (newColor, pham_id))


except mdb.Error, e:
  
	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)
    
finally:
        
    cur.execute("COMMIT")
    if con:    
        con.close()

print "\nYour database has been phixed!"
