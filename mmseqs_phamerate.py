#!/usr/bin/env python
#Python2.7.12 code for executing k-mer based clustering of database proteins
#Initially written by Charles Bowman for using kClust and hhSuite to cluster
#database proteins - 20150126
#Modified by Christian Gauthier to use the much faster MMseqs2 - 20190304

import MySQLdb as mdb
import sys
import os
import random
import colorsys
import datetime
import getpass

#Get the command line parameters
try:
	database = sys.argv[1]
except:
	print "Usage: ./k_phamerate.py <name of MySQL database to be phamerated>"
	print "e.g. ./k_phamerate.py Actino_Draft"
	sys.exit(1)

#GET u/n and pw

username = getpass.getpass(prompt='mySQL username:')
password = getpass.getpass(prompt='mySQL password:')

starttime = datetime.datetime.now()

new_set = set()
pre_phams = {}
pre_colors = {}

#build pre-dictionary
try:
	#Connect
	con = mdb.connect('localhost', username, password, database)
	cur = con.cursor()
	print "Connected to Database - Fetching Phams, GeneIDs, and Translations"
	print "..."

	cur.execute("SELECT a.name, a.GeneID, b.color FROM (SELECT p.name, g.GeneID FROM gene g INNER JOIN pham p ON g.GeneID = p.GeneID) AS a INNER JOIN pham_color AS b ON a.name = b.name ORDER BY a.name ASC")
	tuples = cur.fetchall()

	print "{} old genes".format(len(tuples))
	for t in tuples:
		name, GeneID, color = t
		if name in pre_phams.keys():
			pre_phams[name] = pre_phams[name] | set([GeneID])
		else:
			pre_phams[name] = set([GeneID])
			pre_colors[name] = color

	cur.execute("SELECT geneid FROM gene WHERE geneid NOT IN (SELECT g.GeneID FROM gene g INNER JOIN pham p ON g.geneid = p.geneid)")
	tuples = cur.fetchall()

	print "{} new genes".format(len(tuples))

	for t in tuples:
		GeneID = t[0]
		new_set = new_set | set([GeneID])

except mdb.Error, e:

	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)

finally:
	if con:
		con.close()

###MMseqs2 Analysis Start
#Get the gene information
try:
	#Connect
	con = mdb.connect('localhost', username, password, database)
	cur = con.cursor()
	print "Connected to Database - Fetching GeneIDs and Translations"

	#Clear previous pham data
	cur.execute("TRUNCATE TABLE pham")
	cur.execute("TRUNCATE TABLE pham_color")
	cur.execute("COMMIT")

	#Get gene data
	cur.execute("SELECT GeneID, translation FROM gene")
	tuples = cur.fetchall()

	print "{} total genes".format(len(tuples))

	#Clear any existing temporary MMseqs2 directory
	clear_mmseqs_dir = "rm -r /tmp/MMseqs2/"
	os.system(clear_mmseqs_dir)

	#Build temporary MMseqs2 directory
	make_mmseqs_dir = "mkdir /tmp/MMseqs2/"
	os.system(make_mmseqs_dir)

	#Filter out duplicate protein sequences
	all_gs_ts = []
	for t in tuples:
		all_gs_ts.append([t[0], t[1]])
	translations_and_geneids = {}
	
	for i in all_gs_ts:
		geneids = translations_and_geneids.get(i[1], set())
		geneids.add(i[0])
		translations_and_geneids[i[1]] = geneids
		
	duplicate_groups = {}

	#Build the query file
	f = open('/tmp/MMseqs2/tempquery.txt','w')
	for prot_seq in translations_and_geneids.keys():
		selected_geneid = random.sample(translations_and_geneids[prot_seq], 1)[0]
		duplicate_groups[selected_geneid] = translations_and_geneids[prot_seq]
		f.write(">{}\n{}\n".format(selected_geneid, prot_seq.replace('-','M')))
	
	#Run without filtering duplicate protein sequences
	#for tuple in tuples:
	#	f.write(">" + tuple[0] + '\n')
	#	f.write(tuple[1].replace('-','M') + '\n')
except mdb.Error, e:
	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)

finally:
	f.close()
	if con:
		con.close()

#Create mmseqs database
createdb = ("mmseqs createdb /tmp/MMseqs2/tempquery.txt /tmp/MMseqs2/sequenceDB")
print "MMseqs2 command: " + createdb
os.system(createdb)

#Run MMseqs2 clustering using --cascaded option
clusterdb = ("mmseqs cluster /tmp/MMseqs2/sequenceDB /tmp/MMseqs2/clusterDB /tmp/MMseqs2 --remove-tmp-files --threads 4 -v 0 --min-seq-id 0.40 -c 0.80 --alignment-mode 3 --cov-mode 0 --cluster-mode 0")
print "MMseqs2 command: " + clusterdb
os.system(clusterdb)

#Convert to more readable format
convertdb = ("mmseqs createseqfiledb /tmp/MMseqs2/sequenceDB /tmp/MMseqs2/clusterDB /tmp/MMseqs2/mmseqs_output.txt")
print "MMseqs2 command: " + convertdb
os.system(convertdb)

phams_count = {}

#For this more readable format, the file is essentially a Pearson (FASTA) format where lines that start a new pham begin with a null bit ("\x00"), and the end of file (last line) is also a null bit
#Parse output
f = open("/tmp/MMseqs2/mmseqs_output.txt", "r")
line = f.readline()

post_phams = {}
number = 1
geneids = []

while line != "\x00":
	if line[0] == ">":
		geneids.append(line[1:].rstrip("\n").rstrip(" "))
	elif "\x00" in line:
		post_phams[number] = set(geneids)
		phams_count[number] = len(geneids)
		number += 1
		geneids = [line[2:].rstrip("\n").rstrip(" ")]
	else:
		pass
	line = f.readline()
post_phams[number] = set(geneids)

# re-introduce duplicates to their groups
for pham in post_phams.keys():
	geneids = post_phams[pham]
	new_geneids = []
	for gene in geneids:
		group = duplicate_groups[gene]
		for geneid in group:
			new_geneids.append(geneid)
	post_phams[pham] = set(new_geneids)

total = 0
for pham in post_phams.keys():
	length = len(post_phams[pham])
	total += length
	
print total

f.close()

#Do Phamily Conservationatory Measures
phams_final = {}
phams_colors = {}

post_phams_temp = post_phams.copy()


#f = open("/tmp/post_cons.txt", "w")
#g = open("/tmp/pre_cons.txt", "w")

outcount = 0
total = len(pre_phams)

#Iterate through old and new phams
for pre_key in pre_phams:

	outcount += 1
	print "Pham Name Conservation Search: " + `outcount` + " / " + `total`

	pre = pre_phams[pre_key]
	#g.write(`pre_key` + `pre` + "\n")

	if (pre_key in phams_final.keys()):
			continue;

	for post_key in post_phams_temp:

		post = post_phams_temp[post_key]

		if pre & post == set():
			continue

		#Case 1 + 5 (Identity and Subtraction)
		if pre == post:
			phams_final[pre_key] = post
			phams_colors[pre_key] = pre_colors[pre_key]
			post_phams.pop(post_key, None)
			#f.write(post_key + `post` + "1 + 5\n")
			#print "1 or 5"
			break

		#Case 2 and 4 (Addition and Join) - PHAM GREW
		elif post - pre != set():

			#Case 2 and 4 (Addition and Join)
			if post & new_set != set():

				#Case 4 - Join with new gene
				if (post - (post & new_set)) - pre != set():
					#print "4"
					#f.write(post_key + `post` + "4\n")
					break

				#Case 2 - Addition with new gene
				phams_final[pre_key] = post
				phams_colors[pre_key] = pre_colors[pre_key]
				post_phams.pop(post_key, None)
				#f.write(post_key + `post` + "2\n")
				#print `pre_key` + `phams_final[pre_key]`
				break

			#Case 4 - Join without new gene
			else:
				#f.write(post_key + `post` + "4\n")
				#print "4"
				break

		#Case 3 - split - PHAM SHRANK, BUT NOT BY REMOVAL
		elif pre - post != set():
			#f.write(post_key + `post` + "3\n")
			#print "3"
			break

phams_final[0] = "placeholder"
highest_pham = max(map(int, phams_final.keys())) + 1
del phams_final[0]
#Reassign data for split or joined phams

#f.close()
#g.close()

#f = open("/tmp/phams_final.txt", "w")

#for key in phams_final:
	#f.write(`key` + `phams_final[key]` + "\n")

#f.close()

#f = open("/tmp/post_phams.txt", "w")

for key in post_phams:

	#f.write(key + `post_phams[key]` + "\n")

	new_key = highest_pham
	highest_pham += 1

	phams_final[new_key] = post_phams[key]

	if len(post_phams[key]) > 1:
		h=s=v=0
		while h <= 0: h = random.random()
		while s < 0.5: s = random.random()
		while v < 0.8: v = random.random()
		rgb = colorsys.hsv_to_rgb(h,s,v)
		rgb = (rgb[0]*255,rgb[1]*255,rgb[2]*255)
		hexrgb = '#%02x%02x%02x' % rgb
		phams_colors[new_key] = hexrgb
	else:
		phams_colors[new_key] = '#FFFFFF'

#f.close()
#Insert Data
try:
	con = mdb.connect('localhost', username, password, database)
	print "Connected to Database - Inserting Data"
	cur = con.cursor()

	cur.execute("TRUNCATE TABLE pham")
	cur.execute("TRUNCATE TABLE pham_color")
	cur.execute("COMMIT")

	f = open("/tmp/phinal_pham_insert_log.txt", "w")

	#do the pham inserts
	for key in phams_final:

		for gene in phams_final[key]:
			#print "*****"
			#print "Inserting Gene " + key + " - " + phams_firstiter[key] + " pham - " + phams[key]
			#print ""

			f.write(`gene` + "\t" + `key` + "\n")

			cur.execute("INSERT INTO pham(geneid, name) VALUES (%s, %s)" , (gene, key))
	cur.execute("COMMIT")

	f.close()

	#do the color inserts
	print "Inserting Colors"
	for key in phams_colors:
		#print "Pham " + key + " - " + colors[key]
		cur.execute("INSERT IGNORE INTO pham_color(name, color) VALUES (%s, %s)" , (key, phams_colors[key]))
	cur.execute("COMMIT")

except mdb.Error, e:

	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)

finally:
	if con:
		con.close()

#Do housekeeping on color settings for new phams or orphams
try:
	#Connect
	con = mdb.connect('localhost', username, password, database)
	cur = con.cursor()
	print "Connected to Database - Phetching Phalsely Hued Phams"
	print "...\n"
	#Fixed SQL statement to comply with new MYSQL standard where all items selected that aren't functionally dependent on "group by" must be included in "group by"
	#Error Code: 1055.  Expression #1 of SELECT list is not in GROUP BY clause and contains nonaggregated column 'b.id' which is not functionally dependent on columns in GROUP BY clause; this is incompatible with sql_mode=full_group_by_only
	cur.execute("SELECT * FROM (SELECT b.id, COUNT(GeneID) AS count, a.name, b.color FROM pham AS a INNER JOIN pham_color AS b ON a.name = b.name GROUP BY a.name, b.id) AS c WHERE color = '#FFFFFF' AND count > 1")
	tuples = cur.fetchall()

	print "{} miscolored phamilies found.".format(len(tuples))
	if len(tuples) > 0:
		print " phixing that phor you..."

	for t in tuples:
		pham_id, count, name, color = t

		h=s=v=0
		while h <= 0: h = random.random()
		while s < 0.5: s = random.random()
		while v < 0.8: v = random.random()
		rgb = colorsys.hsv_to_rgb(h,s,v)
		rgb = (rgb[0]*255,rgb[1]*255,rgb[2]*255)
		hexrgb = '#%02x%02x%02x' % rgb
		newColor = hexrgb

		cur.execute("UPDATE pham_color SET color = {} WHERE id = {}".format(newColor, pham_id))

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
	#Fixed SQL statement to comply with new MYSQL standard where all items selected that aren't functionally dependent on "group by" must be included in "group by"
	#Error Code: 1055.  Expression #1 of SELECT list is not in GROUP BY clause and contains nonaggregated column 'b.id' which is not functionally dependent on columns in GROUP BY clause; this is incompatible with sql_mode=only_full_group_by
	cur.execute("SELECT * FROM (SELECT b.id, COUNT(GeneID) AS count, a.name, b.color FROM pham AS a INNER JOIN pham_color AS b ON a.name = b.name GROUP BY a.name, b.id) AS c WHERE color != '#FFFFFF' AND count = 1")
	tuples = cur.fetchall()

	print "{} miscolored orphams found.".format(len(tuples))
	if len(tuples) > 0:
		print " phixing that phor you..."

	for t in tuples:
		pham_id, count, name, color = t
		newColor = "#FFFFFF"

		cur.execute("UPDATE pham_color SET color = {} WHERE id = {}".format(newColor, pham_id))

except mdb.Error, e:

	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)

finally:
	cur.execute("COMMIT")
	if con:
		con.close()

endtime = datetime.datetime.now()

print starttime
print endtime
