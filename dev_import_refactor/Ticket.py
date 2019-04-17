"""Represents a structure to contain directions for how to parse and import
genomes into Phamerator."""






class ImportTicket:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.type = '' #Add, Replace, Remove, UPDATE
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





	#Make sure "none" and "retrieve" indications are lowercase,
	#as well as "action", "status", "feature", and "author" fields are lowercase


    def check_case(self):

    	self.type = self.type.lower()

    	if self.primary_phage_id.lower() == "none":
    		self.primary_phage_id = self.primary_phage_id.lower()

    	if (self.host.lower() == "none" or self.host.lower() == "retrieve"):
    		self.host = self.host.lower()

    	if (self.cluster.lower() == "none" or self.cluster.lower() == "retrieve"):
    		self.cluster = self.cluster.lower()

    	if (self.subcluster.lower() == "none" or self.subcluster.lower() == "retrieve"):
    		self.subcluster = self.subcluster.lower()

        self.status = self.status.lower()

    	self.description_field = self.description_field.lower()

    	if (self.accession.lower() == "none" or self.accession.lower() == "retrieve"):
    		self.accession = self.accession.lower()

    	self.annotation_author = self.annotation_author.lower()


    	if self.secondary_phage_id.lower() == "none":
    		self.secondary_phage_id = self.secondary_phage_id.lower()

    	self.run_mode = self.run_mode.lower()









	# If either the Host, Cluster, Subcluster or Accession data needs to be retrieved,
	# try to access the data in phagesdb before proceeding

    def check_retrieve_status(self):

		#Host
		if self.host == "retrieve":
			self.host = "none"
		#Cluster
		if self.cluster == "retrieve":
			self.cluster = "none"
		#Accession
		if self.accession == "retrieve":
			self.accession = "none"
		#Subcluster
		if self.subcluster == "retrieve":
			self.subcluster = "none"

        #TODO Handle error better
		table_errors += 1





	# Make sure the requested action is permissible
	# This script currently only allows four actions:
    # add, remove, replace, and update.
    def check_type(self):

    	if not (self.type == "add" or \
            self.type == "remove" or \
            self.type == "replace" or \
            self.type == "update"):

            #TODO handle this error better
    		write_out(output_file,"\nError: %s is not a permissible action." %self.type)
    		table_errors += 1










	#Modify fields if needed

	#Modify Host if needed
    def check_host(self,value):

    	if self.host == "retrieve":

    			self.host = value

                #TODO handle error better
    			write_out(output_file,"\nError: unable to retrieve Host data for phage %s from phagesdb." %self.primary_phage_id)
    			self.host = "none"
    			table_errors += 1

    	if self.host != "none":

    		self.host = self.host.split(' ')[0] #Keep only the genus in the host data field and discard the rest

            #TODO handle this error better
            if self.host not in phageHost_set:
    			print "The host strain %s is not currently in the database." % self.host
    			table_errors +=  question("\nError: %s is not the correct host for %s." % (self.host,self.primary_phage_id)) #errors will be incremented if host was not correct








	# Modify Cluster and Subcluster if needed
    def check_cluster_subcluster(self):

    	#Check Subcluster data
    	if self.subcluster != "none":

            #TODO handle errors better
    		if self.subcluster not in phageSubcluster_set:


    			print "The Subcluster %s is not currently in the database." % self.subcluster
    			table_errors +=  question("\nError: %s is not the correct Subcluster for %s." % (self.subcluster,self.primary_phage_id))

    		if len(self.subcluster) > 5:
    			write_out(output_file,"\nError: phage %s Subcluster designation %s exceeds character limit." % (self.primary_phage_id,self.subcluster))
    			table_errors += 1


    	#Check Cluster data
    	if self.cluster != "none":
    		if self.cluster.lower() == "singleton":
    			self.cluster = self.cluster.lower()

    		if (self.cluster not in phageCluster_set and self.cluster != "singleton"):
    			print "The Cluster %s is not currently in the database." % self.cluster
    			table_errors +=  question("\nError: %s is not the correct Cluster for %s." % (self.cluster,self.primary_phage_id))

    		if (self.cluster != "singleton" and len(self.cluster) > 5):
    			write_out(output_file,"\nError: phage %s Cluster designation %s exceeds character limit." % (self.primary_phage_id,self.cluster))
    			table_errors += 1

    		#If Singleton of Unknown Cluster, there should be no Subcluster
    		if (self.cluster == "singleton" or self.cluster == "UNK"):
    			if self.subcluster != "none":
    				write_out(output_file,"\nError: phage %s Cluster and Subcluster discrepancy." % self.primary_phage_id)
    				table_errors += 1

    		#If not Singleton or Unknown or none, then Cluster should be part
    		#of Subcluster data and the remainder should be a digit
    		elif self.subcluster != "none":

    			if (self.subcluster[:len(self.cluster)] != self.cluster or \
    				self.subcluster[len(self.cluster):].isdigit() == False):

    				write_out(output_file,"\nError: phage %s Cluster and Subcluster discrepancy." % self.primary_phage_id)
    				table_errors += 1
    		else:
    			pass



	#Modify Status if needed
    def check_status(self):

    	if (self.status not in phageStatus_set and self.status != "none"):
    			print "The status %s is not currently in the database." % self.status
    			table_errors +=  question("\nError: %s is not the correct status for %s." % (self.status,self.primary_phage_id))
    	if len(self.status) > 5:
    		write_out(output_file,"\nError: the status %s exceeds character limit." % self.status)
    		table_errors += 1


	#Modify Description Qualifier if needed
    def check_description_field(self):

    	if (self.description_field not in description_set and self.description_field != "none"):
    		print self.description_field + " is an uncommon qualifier."
    		table_errors += question("\nError: %s is an incorrect qualifier." % self.description_field)


	#Modify Accession if needed
    def check_accession(self):

    	if self.accession == "retrieve":

    		#On phagesdb, phages should always have a Genbank Accession field. If no accession, it will be ""
    		try:
    			self.accession = online_data_dict['genbank_accession']
    			if self.accession != "":
    				self.accession = self.accession.strip() #Sometimes accession data from phagesdb have whitespace characters
    				self.accession = self.accession.split('.')[0] #Sometimes accession data from phagesdb have version suffix
    			else:
    				self.accession = "none"

    		except:
    			write_out(output_file,"\nError: unable to retrieve Accession data for phage %s from phagesdb." %self.primary_phage_id)
    			self.accession = "none"
    			table_errors += 1

    	elif self.accession.strip() == "":
    		self.accession = "none"



	#Modify AnnotationAuthor
    def check_annotation_author(self):

    	#Author should only be 'hatfull','gbk', or 'none'.
    	if self.annotation_author == "hatfull":
    		self.annotation_author = "1"
    	elif self.annotation_author == "gbk":
    		self.annotation_author = "0"
    	elif self.annotation_author == "none":
    		self.annotation_author = "none"
    	else:
    		self.annotation_author = "error"



	#Make sure run mode is permissible
    def check_run_mode(self):

    	if self.run_mode not in set(run_mode_options_dict.keys()):
    		write_out(output_file,"\nError: run mode is not valid for phage %s." %self.primary_phage_id)
    		table_errors += 1
    	elif self.run_mode == 'custom':

            #TODO make this variable accessible
    		run_mode_custom_total += 1
    	else:
    		pass





















	#Rules for how each field is populated differs depending on each specific action
    def check_ticket(self):


    	#Update
    	if self.type == "update":

    		#FirstPhageID
    		if self.primary_phage_id not in phageId_set:
    			write_out(output_file,"\nError: %s is not a valid PhageID in the database." %primary_phage_id)
    			table_errors += 1

    		#Host, Cluster, Status
    		if (self.host == "none" or \
    			self.cluster == "none" or \
    			self.status == "none"):

    			write_out(output_file,"\nError: %s does not have correctly populated HostStrain, Cluster, Subcluster, or Status fields." %self.primary_phage_id)
    			table_errors += 1

    		#Description
    		if self.description_field != "none":
    			write_out(output_file,"\nError: %s does not have correctly populated Description field." %self.primary_phage_id)
    			table_errors += 1

    		#SecondPhageID
    		if self.secondary_phage_id != "none":
    			write_out(output_file,"\nError: %s should not have a genome listed to be removed." %self.primary_phage_id)
    			table_errors += 1

    		#Accession = it will either be an accession or it will be "none"
    		#Subcluster = it will either be a Subcluster or it will be "none"

    		#Author
    		if self.annotation_author != '1' and self.annotation_author != '0':
    			write_out(output_file,"\nError: %s does not have correctly populated Author field." %self.primary_phage_id)
    			table_errors += 1

    		#Run Mode
    		if self.run_mode != 'none':
    			write_out(output_file,"\nError: %s does not have correctly populated Run Mode field." %self.primary_phage_id)
    			table_errors += 1


    	#Add
    	elif self.type == "add":

    		#FirstPhageID
    		if self.primary_phage_id in phageId_set:
    			write_out(output_file,"\nError: %s is already a PhageID in the database. This genome cannot be added to the database." %self.primary_phage_id)
    			table_errors += 1

    		#FirstPhageID, Host, Cluster, Status, Description
    		if (self.primary_phage_id == "none" or \
    			self.host == "none" or \
    			self.cluster == "none" or \
    			self.status == "none" or \
    			self.description_field == "none"):

    			write_out(output_file,"\nError: %s does not have correctly populated fields." %self.primary_phage_id)
    			table_errors += 1

    		#Status
    		if self.status == "final":
    			print self.primary_phage_id + " to be added is listed as Final status, but no Draft (or other) genome is listed to be removed."
    			table_errors +=  question("\nError: %s is not the correct status for %s." % (self.status,self.primary_phage_id))

    		#SecondPhageID
    		if self.secondary_phage_id != "none":
    			write_out(output_file,"\nError: %s to be added should not have a genome indicated for removal." %self.primary_phage_id)
    			table_errors += 1

    		#Accession = it will either be an accession or it will be "none"
    		#Subcluster = it will either be a Subcluster or it will be "none"

    		#Author
    		if self.annotation_author != '1' and self.annotation_author != '0':
    			write_out(output_file,"\nError: %s does not have correctly populated Author field." %self.primary_phage_id)
    			table_errors += 1

    		#Run Mode
    		if self.run_mode == 'none':
    			write_out(output_file,"\nError: %s does not have correctly populated Run Mode field." %self.primary_phage_id)
    			table_errors += 1


    	#Remove
    	elif self.type == "remove":

    		#FirstPhageID,Host, Cluster, Subcluster, Status, Description, Accession, Author, Run Mode
    		if (self.primary_phage_id != "none" or \
    			self.host != "none" or \
    			self.cluster != "none" or \
    			self.status != "none" or \
    			self.description_field != "none" or \
    			self.accession != "none" or \
    			self.subcluster != "none" or \
    			self.annotation_author != "none" or \
    			self.run_mode != "none"):


    			write_out(output_file,"\nError: %s to be removed does not have correctly populated fields." %self.secondary_phage_id)
    			table_errors += 1

    		#SecondPhageID
    		if self.secondary_phage_id not in phageId_set:
    			write_out(output_file,"\nError: %s is not a valid PhageID. This genome cannot be dropped from the database." %self.secondary_phage_id)
    			table_errors += 1


    	#Replace
    	elif self.type == "replace":

    		#FirstPhageID
    		if self.primary_phage_id == "none":
    			write_out(output_file,"\nError: %s is not a valid PhageID. %s genome cannot be replaced." % (self.primary_phage_id,self.secondary_phage_id))
    			table_errors += 1

    		#FirstPhageID. If replacing a genome, ensure that if the genome to
    		#be removed is not the same, that the new genome added has a unique name
    		if (self.primary_phage_id in phageId_set and self.primary_phage_id != self.secondary_phage_id):
    			write_out(output_file,"\nError: %s is already a PhageID in the database. This genome cannot be added to the database." %self.primary_phage_id)
    			table_errors += 1

    		#Host,Cluster,Status,Description
    		if (self.host == "none" or \
    			self.cluster == "none" or \
    			self.status == "none" or \
    			self.description_field == "none"):

    			write_out(output_file,"\nError: %s does not have correctly populated fields." %self.primary_phage_id)
    			table_errors += 1

    		#SecondPhageID
    		if self.secondary_phage_id not in phageId_set:
    			write_out(output_file,"\nError: %s is not a valid PhageID. This genome cannot be dropped from the database." %self.secondary_phage_id)
    			table_errors += 1

    		if self.primary_phage_id != self.secondary_phage_id:
    			print "%s to replace %s is not spelled the same." %(self.primary_phage_id,self.secondary_phage_id)
    			table_errors +=  question("\nError: Phage %s is not spelled the same as phage %s." % (self.primary_phage_id,self.secondary_phage_id))

    		#Accession = it will either be an accession or it will be "none"
    		#Subcluster = it will either be a Subcluster or it will be "none"

    		#Author
    		if self.annotation_author != '1' and self.annotation_author != '0':
    			write_out(output_file,"\nError: %s does not have correctly populated Author field." %self.primary_phage_id)
    			table_errors += 1

    		#Run Mode
    		if self.run_mode == 'none':
    			write_out(output_file,"\nError: %s does not have correctly populated Run Mode field." %self.primary_phage_id)
    			table_errors += 1



    	else:
    		pass








###
