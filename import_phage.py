#Phage Import Script
#Charles Bowman and Damani Brown
#University of Pittsburgh
#20150216
#WORKS FOR .gb or .txt files

#!/usr/bin/env python
import sys, os, getpass
from Bio import SeqIO
import MySQLdb as mdb

#Get the command line parameters
try:
    database = sys.argv[1]
    phageListDir = sys.argv[2]
except:
    print "Incorrect Parameters - python import_phage.py DATABASE DIRECTORY"
    sys.exit(1)

#GET u/n and pw
username = getpass.getpass(prompt='mySQL username:')
password = getpass.getpass(prompt='mySQL password:')
joinCounter = 0

genecount = 0;
geneIdSet = set()

f = open("/tmp/phage_import_log.txt", "w")
print "Adding Phages:"

#Iterating over each file in the directory. 

files =  [X for X in os.listdir(phageListDir) if os.path.isfile(os.path.join(phageListDir,X))]

print files

for seqFile in files:
    locusCounter = 1
    seqFile = phageListDir + seqFile
    if(str(seqFile[-3:]) == ".gb" or str(seqFile[-3:]) == "gbf" or str(seqFile[-3:]) == "txt"):
        statements = []
        #The file is parsed to grab the header information
        for seq_record in SeqIO.parse(seqFile, "genbank"):
            
            #This block is essential, skip adding if errors
            try:
                #sequence
                phageSeq = seq_record.seq.upper()
                #length
                seqLength = len(phageSeq)
                #GC
                seqGC = 100 * (float(phageSeq.count('G')) + float(phageSeq.count('C'))) / float(seqLength)
                #name (ID)
                
                phageName = seq_record.annotations["organism"].split(' ')[-1]
                if phageName == "Unclassified.":
                    phageName = seq_record.annotations["organism"].split(' ')[-2]
                print phageName
                
            except:
                f.write("Problem adding " + `seqFile` + " in Header Block\n")
                continue

            #Nonessential Stuff
            try:
                phageHost = seq_record.annotations["organism"].split(' ')[0]
                if phageHost == phageName:
                	phageHost = ""
            except:
                phageHost = ""
            
            try:
                phageNotes = str(seq_record.annotations["comment"])
            except:
                phageNotes = ""

            try:
                #accession
                accessionNum = seq_record.annotations["accessions"][0]
            except:
                accessionNum = ""
            
            
            #print("Phage: " + phageName + "\nAccession Number: " + accessionNum)
            #print("The host of this phage is: " +phageHost)
            #print("Length " + str(len(phageSeq)))
            #print("Phage Information \n================== \n"
            #      + phageNotes + "\n")
            
            #Next each CDS feature is parsed from the file, and we grab and store that information

            statements.append("""INSERT INTO phage (PhageID, Accession, Name, HostStrain, Sequence, SequenceLength, GC) VALUES ("%s","%s","%s","%s","%s",%s,%s)""" % (phageName, accessionNum, phageName, phageHost, phageSeq, seqLength, seqGC))

            CDSCount = 0
            addCount = 0
            for feature in seq_record.features:
                if feature.type == "CDS":
                    CDSCount += 1

                    #This block is essential, skip adding if errors
                    try:

                        #GeneID
                        try:
                            geneID = feature.qualifiers["locus_tag"][0]
                        except:
                            geneID = phageName+"_"+str(locusCounter)

                        if(geneID in geneIdSet):
                            dup = True
                            while dup:
                                geneID = input("Duplicate GeneID for " + `geneID` + ". Enter new GeneID:")
                                if (geneID not in geneIDSet):
                                    dup = False
                        else:
                            geneIdSet | set(geneID)

                        #Name
                        if (geneID.split('_')[-1].isdigit()):
                            geneName = geneID.split('_')[-1]
                        else:
                            geneName = locusCounter                            

                        #Start/Stop
                        if str(feature.location)[:4] == "join":

                            temp = str(feature.location).strip("join{}").split(",")[1]
                            startCoord = int(temp.split(":")[0].split("[")[1])
                            stopCoord = int(temp.split(":")[1].split("]")[0])

                        else:

                            strStart = str(feature.location.start)
                            strStop = str(feature.location.end)

                            if (strStart.isdigit() and strStop.isdigit()):
                                startCoord = int(strStart)
                                stopCoord = int(strStop)
                            else:
                                f.write("Non-traditional Coord " + geneID + " " + strStart + " " + strStop + "\n")
                                locusCounter += 1 
                                continue

                        #Length (via Translation)
                        translation = feature.qualifiers["translation"][0]
                        geneLen = (len(translation) * 3) + 3  #Add 3 for the stop codon...

                        #Parse orientation
                        if feature.strand == 1:
                            orientation = "Forward"
                        elif feature.strand == -1:
                            orientation = "Reverse"
                        #ssRNA phages
                        elif feature.strand is None:
                            orientation = "Forward"
                        else:
                            locusCounter += 1
                            f.write("Error Adding Gene " + `locusCounter` + " of " + phageName + "\n")
                            print feature
                            continue

                        '''#codon
                        startCodon = ""
                        stopCodon = ""
                        if(orientation == "Forward"):
                            for x in range(0,3):
                                startCodon = startCodon + phageSeq[startCoord + x]
                            for x in range(2,-1,-1):
                                stopCodon = stopCodon + phageSeq[stopCoord - x]
                        else:
                            for x in range(0,3):
                                temp = phageSeq[startCoord + x]
                                if temp == 'C':
                                    temp = 'G'
                                elif temp == 'G':
                                    temp = 'C'
                                elif temp == 'A':
                                    temp = 'T'
                                elif temp == 'T':
                                    temp = 'A'
                                stopCodon = stopCodon + temp
                                
                            for x in range(2,-1,-1):
                                temp = phageSeq[stopCoord - x]
                                if temp == 'C':
                                    temp = 'G'
                                elif temp == 'G':
                                    temp = 'C'
                                elif temp == 'A':
                                    temp = 'T'
                                elif temp == 'T':
                                    temp = 'A'
                                startCodon = startCodon + temp'''

                    except:
                        locusCounter += 1
                        f.write( "Error Adding Gene " + `locusCounter` + " of " + phageName + "\n")
                        f.write( "Error: " + `sys.exc_info()[0]` + " at line " + `sys.exc_info()[-1].tb_lineno` + "\n")
                        continue
                          
                    typeID = feature.type

                    try:
                        if feature.qualifiers["product"][0].lower() == "hypothetical protein":
                            featureNote = ""
                        else:
                            featureNote = feature.qualifiers["product"][0]
                    except:
                        True

                    addCount+= 1
                    locusCounter += 1
                    genecount += 1


                    statements.append("""INSERT INTO gene (GeneID, PhageID, Start, Stop, Length, Name, TypeID, translation, Orientation, Notes) VALUES ("%s","%s",%s,%s,%s,"%s","%s","%s","%s","%s");""" % (geneID, phageName,startCoord, stopCoord, geneLen, geneName, typeID, translation, orientation[0], featureNote)) 
                    
                    #print(locusTag)
                    #print(typeID)
                    #print("Type of protein is: " + str(typeOfProtein))
                    #print("This is the length of the CDS: "+ str(lengthOfCDS))
                    #print("Coordinates "+ str(startCoord) +":"+ str(stopCoord))
                    #print("Start Sequence: " +startCodon)
                    #print("Stop Sequence: " +stopCodon)
                    #print("The orientation of the genome is: " + orientation + "\n")
                    
                else:
                    continue

            #this is your confirmation
            f.write( `phageName` + ": " + `addCount` + "/" + `CDSCount` + "\n")
          
        #SQL time
        try:
            con = mdb.connect('localhost', username, password, database)
            con.autocommit(False)
            cur = con.cursor()

            cur.execute("START TRANSACTION")

            for statement in statements:
                cur.execute(statement)

            cur.execute("COMMIT")
            con.autocommit(True)

        except:
            print "error inserting " + `phageName` + "\n"
            f.write("error inserting " + `phageName` + "\n")
            f.write( "Error: " + `sys.exc_info()[0]`+ ":" +  `sys.exc_info()[1]` + "at: " + `sys.exc_info()[2]` + "\n")
            cur.execute("ROLLBACK")
            cur.execute("SET autocommit = 1")
            

    else:
        True
f.close()
