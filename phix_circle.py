#!/usr/bin/env python
#Code for doing BLAST and Clustal alignments for phamily circles
#charles bowman
#cab106@pitt.edu
#20150325

#Suppress biopython depreciation warnings
import warnings
warnings.filterwarnings("ignore")

import os, sys, subprocess
import Bio
from itertools import *
from Bio import AlignIO
from Bio.Blast import NCBIStandalone
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline
import MySQLdb as mdb
import getpass

#Globals
blast_exe = "/home/cbowman/Applications/BLAST/bin/bl2seq"
query_file = "/home/cbowman/Desktop/Query.fasta"
subject_file = "/home/cbowman/Desktop/subject.fasta"
result_file = "/home/cbowman/Desktop/result.txt"
clustalw_infile = "/home/cbowman/Desktop/clustalw.fasta"

#define db vars
try:
	database = sys.argv[1]
	pham = sys.argv[2]
except:
	print "Incorrect Parameters!\n./phix_circle.py DATABASE PHAM"
	sys.exit(1)

#get SQL stuff
username = getpass.getpass(prompt='mySQL username:')
password = getpass.getpass(prompt='mySQL password:')

#connect to SQL
try:
	#Connect
	con = mdb.connect('localhost', username, password, database)
	cur = con.cursor()
	print "Connected to Database - Phetching Phamily data"
	print "...\n"

	#gene geneid and translation for genes in pham
	cur.execute("select g.geneid, g.translation from gene as g inner join pham as p on g.geneid = p.GeneID where p.name = " + `pham`)
	tuples = cur.fetchall()

except mdb.Error, e:
  
	print "Error %d: %s" % (e.args[0],e.args[1])
	sys.exit(1)
    
finally:

    if con:    
        con.close()

#process the data into some data structures
genes = []
translations = {}

for tuple in tuples:
	genes.append(tuple[0])
	translations[tuple[0]] = tuple[1]

#get the minimal combination of alignments
tasks = combinations(genes, 2)

#ADD CODE FOR not having to redo stuff here

deletes = [];
inserts = [];
#Do the work
for task in tasks:
	print task
	
	#Write the needed files
	f = open(query_file, 'w')
	f.write(">" + task[0] + "\n" + translations[task[0]])
	f.close()
	f = open(subject_file, 'w')
	f.write(">" + task[1] + "\n" + translations[task[1]])
	f.close()
	f = open(clustalw_infile, 'w')
	f.write(">" + task[0] + "\n" + translations[task[0]] + "\n")
	f.write(">" + task[1] + "\n" + translations[task[1]])
	f.close()

	#BLAST code
	bashCom = '%s -p blastp -i %s -D 1 -F F -j %s -o %s' % (blast_exe, query_file, subject_file, result_file)
	#print bashCom
	os.system(bashCom)
	f = open(result_file, 'r').readlines()
	result = f[3].split('\t')
	blast_e = result[-2]
	blast_bit = int(result[-1])

	#Clustal code
	cline = ClustalwCommandline("clustalw", infile = clustalw_infile)
	stdout, stderr = cline()
	alignment = AlignIO.read(clustalw_infile.replace('.fasta', '.aln'), "clustal")
	length = alignment.get_alignment_length()
	stars = alignment._star_info.count('*')
	clustal_score = float(stars)/length

	#Stow away the remove statement
	deletes.append('DELETE FROM scores_summary where query like "%s" and subject like "%s"' % (task[0], task[1]))

	#stow away the insert statement
	inserts.append('REPLACE INTO scores_summary (query, subject, blast_score, clustalw_score, blast_bit_score) VALUES ("%s", "%s", %s, %s, %s)' % (task[0], task[1], blast_e,clustal_score,blast_bit))

print "...\nRemoving old then inserting new data into database\n..."
#do Inserts
try:
    con = mdb.connect('localhost', username, password, database)
    con.autocommit(False)
    cur = con.cursor()

    cur.execute("START TRANSACTION")

    for delete in deletes:
    		cur.execute(delete)

    for insert in inserts:
        cur.execute(insert)

    cur.execute("COMMIT")
    con.autocommit(True)

except mdb.Error, e:
    print ("error inserting data")
    print insert
    print "Error %d: %s" % (e.args[0],e.args[1])
    cur.execute("ROLLBACK")
    cur.execute("SET autocommit = 1")
    con.close()

print "Done!  You can now draw your phamily circle!"