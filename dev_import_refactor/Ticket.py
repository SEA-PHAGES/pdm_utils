"""Represents a structure to contain directions for how to parse and import
genomes into Phamerator."""






class ImportTicket:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.ticket_type = '' #Add, Replace, Remove, UPDATE
        self.primary_phage_id = '' #Genome that will be added or updated
        self.host = ''
        self.cluster = '' #Singleton should be reported as 'singleton'
        self.subcluster = ''
        self.status = ''
        self.annotation_author = ''
        self.description_field = ''
        self.accession = ''
        self.run_mode = ''
        self.secondary_phage_id = '' #Genomd that will be removed or replaced











	#0 = Import action (unchanged)
	#1 = New PhageID (unchanged)
	#2 = Host (unchanged)
	#3 = Cluster (unchanged)
	#4 = Status
	#5 = Feature field (unchanged)
	#6 = PhageID to be removed
	#7 = Accession
	#8 = Subcluster
	#9 = AnnotationAuthor
	#10 = Run Mode



	#Make sure "none" and "retrieve" indications are lowercase,
	#as well as "action", "status", "feature", and "author" fields are lowercase


    #Todo finish this function:
    def clean_import_ticket(self):


        self.ticket_type
        self.primary_phage_id
        self.host
        self.cluster
        self.subcluster
        self.status
        self.annotation_author
        self.description_field
        self.accession
        self.run_mode
        self.secondary_phage_id



	row[0] = row[0].lower()
	if row[1].lower() == "none":
		row[1] = row[1].lower()
	if (row[2].lower() == "none" or row[2].lower() == "retrieve"):
		row[2] = row[2].lower()
	if (row[3].lower() == "none" or row[3].lower() == "retrieve"):
		row[3] = row[3].lower()
	row[4] = row[4].lower()
	row[5] = row[5].lower()
	if row[6].lower() == "none":
		row[6] = row[6].lower()
	if (row[7].lower() == "none" or row[7].lower() == "retrieve"):
		row[7] = row[7].lower()
	if (row[8].lower() == "none" or row[8].lower() == "retrieve"):
		row[8] = row[8].lower()
	row[9] = row[9].lower()
	row[10] = row[10].lower()











    #Below may need to go to the 'prepare_import_tickets' script

	#If either the Host, Cluster, Subcluster or Accession data needs to be retrieved,
	#try to access the data in phagesdb before proceeding
	if (row[2] == "retrieve" or \
		row[3] == "retrieve" or \
		row[7] == "retrieve" or \
		row[8] == "retrieve"):
		try:

			phage_url = api_prefix + row[1] + api_suffix
			online_data_json = urllib.urlopen(phage_url)
			online_data_dict = json.loads(online_data_json.read())
		except:
			phage_url = ""
			online_data_json = ""
			online_data_dict = {}

			write_out(output_file,"\nError: unable to retrieve Host, Cluster, Subcluster, or Accession data for phage %s from phagesdb." %row[1])
			#Host
			if row[2] == "retrieve":
				row[2] = "none"
			#Cluster
			if row[3] == "retrieve":
				row[3] = "none"
			#Accession
			if row[7] == "retrieve":
				row[7] = "none"
			#Subcluster
			if row[8] == "retrieve":
				row[8] = "none"
			table_errors += 1


	#Make sure the requested action is permissible
	#This script currently only allows four actions: add, remove, replace, and update
	if row[0] in action_set:
		if row[0] == "add":
			add_total += 1
		elif row[0] == "remove":
			remove_total += 1
		elif row[0] == "replace":
			replace_total += 1
		elif row[0] == "update":
			update_total += 1
		else:
			pass
	else:
		write_out(output_file,"\nError: %s is not a permissible action." %row[0])
		table_errors += 1
