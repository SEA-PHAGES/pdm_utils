"""Represents a structure to contain directions for how to parse and import
genomes into Phamerator."""

import Eval




class ImportTicket:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.type = '' # Add, Replace, Remove, UPDATE
        self.primary_phage_id = '' # Genome that will be added or updated
        self.host = ''
        self.cluster = '' # Singleton should be reported as 'singleton'
        self.subcluster = ''
        self.status = ''
        self.annotation_author = ''
        self.description_field = ''
        self.accession = ''
        self.run_mode = ''
        self.secondary_phage_id = '' # Genome that will be removed or replaced




        # Initialize calculated attributes
        self.match_strategy = '' # phage_id or filename
        self.evaluations = []





    def set_evaluation(self, type, message1, message2 = None):

        if type == "warning":
            eval_object = Eval.construct_warning(message1,message2)

        elif type == "error":
            eval_object = Eval.construct_error(message1)

        else:
            eval_object = Eval.EvalResult()

        self.evaluations.append(eval_object)





	# Make sure "none" and "retrieve" indications are lowercase, as well as
    # "action", "status", "feature", and "author" fields are lowercase.


    def check_case(self):

        self.type = self.type.lower()

        if self.primary_phage_id.lower() == "none":
            self.primary_phage_id = self.primary_phage_id.lower()

        if (self.host.lower() == "none" or \
            self.host.lower() == "retrieve"):

            self.host = self.host.lower()

        if (self.cluster.lower() == "none" or \
            self.cluster.lower() == "retrieve"):

            self.cluster = self.cluster.lower()

        if (self.subcluster.lower() == "none" or
            self.subcluster.lower() == "retrieve"):

            self.subcluster = self.subcluster.lower()

        self.status = self.status.lower()

        self.description_field = self.description_field.lower()

        if (self.accession.lower() == "none" or \
            self.accession.lower() == "retrieve"):

            self.accession = self.accession.lower()

            self.annotation_author = self.annotation_author.lower()

        if self.secondary_phage_id.lower() == "none":
            self.secondary_phage_id = self.secondary_phage_id.lower()

        self.run_mode = self.run_mode.lower()




    # If either the Host, Cluster, Subcluster or Accession data needs to be
    # retrieved, try to access the data in phagesdb before proceeding

    def check_retrieve_status(self):

        #Host
        if self.host == "retrieve":
            self.host = "none"

            self.set_evaluation("error", \
            "Host data was not retrieved from Phagesdb.")

        #Cluster
        if self.cluster == "retrieve":
            self.cluster = "none"

            self.set_evaluation("error", \
            "Cluster data was not retrieved from Phagesdb.")

        #Subcluster
        if self.subcluster == "retrieve":
            self.subcluster = "none"

            self.set_evaluation("error", \
            "Subcluster data was not retrieved from Phagesdb.")

        #Accession
        if self.accession == "retrieve":
            self.accession = "none"

            self.set_evaluation("error", \
            "Accession data was not retrieved from Phagesdb.")






    # Make sure the requested action is permissible
    # This script currently only allows four actions:
    # add, remove, replace, and update.
    def check_type(self):

        if not (self.type == "add" or \
            self.type == "remove" or \
            self.type == "replace" or \
            self.type == "update"):

            self.set_evaluation("error", \
            "The ticket type %s is not valid." % self.type)







    # Modify fields if needed

    # Modify Host if needed
    def check_host(self,value,host_set):


        #TODO this error is needed unless elsewhere in script an error is thrown if unable to retrieve data from phagesdb
        # if self.host == "retrieve":
        #
        #     self.host = value
        #
        #     write_out(output_file,"\nError: unable to retrieve Host data for phage %s from phagesdb." %self.primary_phage_id)
        #     self.host = "none"
        #     table_errors += 1

        if self.host != "none":

            # Keep only the genus in the host data field and discard the rest
            self.host = self.host.split(' ')[0]

            #TODO need to implement it elsewhere?
            #TODO define host_set
            if self.host not in host_set:

                self.set_evaluation("warning", \
                    "The host strain %s is not currently in the database." % \
                    self.host, \
                    "The host strain %s is not correct." % self.host)




	# Modify Cluster and Subcluster if needed
    def check_cluster_subcluster(self):

        # Check Subcluster data
        if self.subcluster != "none":

            if self.subcluster not in phageSubcluster_set:

                self.set_evaluation("warning", \
                    "The Subcluster %s is not currently in the database." % \
                    self.subcluster, \
                    "The Subcluster %s is not correct." % self.subcluster)


            if len(self.subcluster) > 5:

                self.set_evaluation("error", \
                    "The Subcluster designation %s exceeds the character limit." % \
                    self.subcluster)


        # Check Cluster data
        if self.cluster != "none":
            if self.cluster.lower() == "singleton":
                self.cluster = self.cluster.lower()

            if (self.cluster not in phageCluster_set and \
                self.cluster != "singleton"):


                self.set_evaluation("warning", \
                    "The Cluster %s is not currently in the database." % \
                    self.cluster, \
                    "The Cluster %s is not correct." % self.cluster)


            if (self.cluster != "singleton" and len(self.cluster) > 5):

                self.set_evaluation("error", \
                    "The Cluster designation %s exceeds the character limit." % \
                    self.cluster)


            # If Singleton of Unknown Cluster, there should be no Subcluster
            if (self.cluster == "singleton" or self.cluster == "UNK"):
                if self.subcluster != "none":

                    self.set_evaluation("error", \
                        "There is a discrepancy between the Cluster and Subcluster \
                        designations.")



            # If not Singleton or Unknown or none, then Cluster should be part
            # of Subcluster data and the remainder should be a digit
            elif self.subcluster != "none":

                if (self.subcluster[:len(self.cluster)] != self.cluster or \
                    self.subcluster[len(self.cluster):].isdigit() == False):

                    self.set_evaluation("error", \
                        "There is a discrepancy between the Cluster and Subcluster \
                        designations.")

            else:
                pass

	# Modify Status if needed
    def check_status(self):

        if (self.status not in phageStatus_set and self.status != "none"):

            self.set_evaluation("warning", \
                "The status %s is not currently in the database." % self.status, \
                "The status %s is not correct." % self.status)

        if len(self.status) > 5:

            self.set_evaluation("error", \
                "The status designation %s exceeds the character limit." % \
                self.status)


    #Modify Description Qualifier if needed
    def check_description_field(self):

        if (self.description_field not in description_set and self.description_field != "none"):

            self.set_evaluation("warning", \
                "The description field %s is not commonly used." % self.description_field, \
                "The description field %s is not correct." % self.description_field)




    # Modify Accession if needed
    def check_accession(self):

        if self.accession == "retrieve":

            #TODO handle this better
            # On phagesdb, phages should always have a Genbank Accession field. If no accession, it will be ""
            try:
                self.accession = online_data_dict['genbank_accession']

                if self.accession != "":

                    # Sometimes accession data from phagesdb have whitespace
                    # characters or the version suffix
                    self.accession = self.accession.strip()
                    self.accession = self.accession.split('.')[0]

                else:

                    self.accession = "none"

            #TODO handle error better.
            except:
                write_out(output_file,"\nError: unable to retrieve Accession data for phage %s from phagesdb." %self.primary_phage_id)
                self.accession = "none"
                table_errors += 1

        elif self.accession.strip() == "":
            self.accession = "none"



	# Modify AnnotationAuthor
    def check_annotation_author(self):

        # Author should only be 'hatfull','gbk', or 'none'.
        if self.annotation_author == "hatfull":
            self.annotation_author = "1"
        elif self.annotation_author == "gbk":
            self.annotation_author = "0"
        elif self.annotation_author == "none":
            self.annotation_author = "none"
        else:
            self.set_evaluation("error", \
                "The annotation author designation %s is not correct." % \
                self.annotation_author)

            self.annotation_author = "error"



	# Make sure run mode is permissible
    def check_run_mode(self,run_mode_options_dict):

        if self.run_mode not in set(run_mode_options_dict.keys()):

            self.set_evaluation("error", \
                "The run mode %s is not valid." % self.run_mode)

        elif self.run_mode == 'custom':

            #TODO make this variable accessible
            run_mode_custom_total += 1

        else:
            pass



    # Rules for how each field is populated differs depending on each specific action
    def check_ticket(self,phage_id_set):

        if self.type == "update":
            self.check_update_ticket(phage_id_set)

        elif self.type == "add":
            self.check_add_ticket(phage_id_set)

        elif self.type == "remove":
            self.check_remove_ticket(phage_id_set)

        elif self.type == "replace":
            self.check_replace_ticket(phage_id_set)

        else:
            pass



    def check_update_ticket(self,phage_id_set):

        if self.primary_phage_id not in phageId_set:
            self.set_evaluation("error", \
                "The %s is not a valid PhageID in the database." % self.primary_phage_id)

        if self.host == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Host field." % self.primary_phage_id)

        if self.cluster == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Cluster field." % self.primary_phage_id)

        if self.subcluster == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Subcluster field." % self.primary_phage_id)

        if self.status == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Status field." % self.primary_phage_id)

        if self.description_field != "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Description field." % self.primary_phage_id)

        if self.secondary_phage_id != "none":
            self.set_evaluation("error", \
                "Phage %s should not have a genome listed to be removed." % self.primary_phage_id)


        #Accession = it will either be an accession or it will be "none"
        #Subcluster = it will either be a Subcluster or it will be "none"

        if self.annotation_author != '1' and self.annotation_author != '0':
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Author field." % self.primary_phage_id)

        if self.run_mode != 'none':
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Run Mode field." % self.primary_phage_id)



    def check_add_ticket(self,phage_id_set):

        if self.primary_phage_id in phage_id_set:

            self.set_evaluation("error", \
                "Phage %s is already a PhageID in the database. \
                This genome cannot be added to the database." % self.primary_phage_id)

        if self.primary_phage_id == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Primary PhageID field." % self.primary_phage_id)

        if self.host == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Host field." % self.primary_phage_id)

        if self.cluster == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Cluster field." % self.primary_phage_id)

        if self.subcluster == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Subcluster field." % self.primary_phage_id)

        if self.status == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Status field." % self.primary_phage_id)

        if self.description_field == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Description field." % self.primary_phage_id)

        if self.status == "final":
            self.set_evaluation("warning", \
                "The phage %s to be added is listed as Final status, but no Draft (or other) genome is listed to be removed." % self.primary_phage_id,
                "The phage %s does not have the correct status" % self.primary_phage_id)

        if self.secondary_phage_id != "none":
            self.set_evaluation("error", \
                "Phage %s to be added should not have a genome indicated for removal." % self.primary_phage_id)

        #Accession = it will either be an accession or it will be "none"
        #Subcluster = it will either be a Subcluster or it will be "none"

        if self.annotation_author != '1' and self.annotation_author != '0':
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Author field." % self.primary_phage_id)

        if self.run_mode == 'none':
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Run Mode field." % self.primary_phage_id)




    def check_remove_ticket(self,phage_id_set):

        if self.primary_phage_id != "none":
            self.set_evaluation("error", \
                "Secondary Phage %s does not have correctly populated Primary PhageID field." % self.secondary_phage_id)

        if self.host != "none":
            self.set_evaluation("error", \
                "Secondary Phage %s does not have correctly populated Host field." % self.secondary_phage_id)

        if self.cluster != "none":
            self.set_evaluation("error", \
                "Secondary Phage %s does not have correctly populated Cluster field." % self.secondary_phage_id)

        if self.status != "none":
            self.set_evaluation("error", \
                "Secondary Phage %s does not have correctly populated Status field." % self.secondary_phage_id)

        if self.description_field != "none":
            self.set_evaluation("error", \
                "Secondary Phage %s does not have correctly populated Description Field field." % self.secondary_phage_id)

        if self.accession != "none":
            self.set_evaluation("error", \
                "Secondary Phage %s does not have correctly populated Accession field." % self.secondary_phage_id)

        if self.subcluster != "none":
            self.set_evaluation("error", \
                "Secondary Phage %s does not have correctly populated Subcluster field." % self.secondary_phage_id)

        if self.annotation_author != "none":
            self.set_evaluation("error", \
                "Secondary Phage %s does not have correctly populated AnnotationAuthor field." % self.secondary_phage_id)

        if self.run_mode != "none":
            self.set_evaluation("error", \
                "Secondary Phage %s does not have correctly populated Run Mode field." % self.secondary_phage_id)

        if self.secondary_phage_id not in phage_id_set:
            self.set_evaluation("error", \
                "Secondary Phage %s is not a valid PhageID. This genome cannot be dropped from the database." % self.secondary_phage_id)





    def check_replace_ticket(self,phage_id_set):

        if self.primary_phage_id == "none":

            self.set_evaluation("error", \
                "Primary Phage %s is not a valid PhageID. %s genome cannot be replaced." % (self.primary_phage_id,self.secondary_phage_id))


        #FirstPhageID. If replacing a genome, ensure that if the genome to
        #be removed is not the same, that the new genome added has a unique name
        if (self.primary_phage_id in phage_id_set and \
            self.primary_phage_id != self.secondary_phage_id):

            self.set_evaluation("error", \
                "Primary Phage %s is already a PhageID in the database. This genome cannot be added to the database." % self.primary_phage_id)

        if self.host == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Host field." % self.primary_phage_id)

        if self.cluster == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Cluster field." % self.primary_phage_id)

        if self.status == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Status field." % self.primary_phage_id)

        if self.description_field == "none":
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Description Field field." % self.primary_phage_id)

        if self.secondary_phage_id not in phage_id_set:
            self.set_evaluation("error", \
                "Secondary Phage %s is not a valid PhageID. This genome cannot be dropped from the database." % self.secondary_phage_id)

        if self.primary_phage_id != self.secondary_phage_id:
            self.set_evaluation("warning", \
                "The Primary Phage %s and Secondary phage %s are not spelled the same." % (self.primary_phage_id,self.secondary_phage_id),
                "The Primary Phage %s and Secondary phage %s are not spelled the same." % (self.primary_phage_id,self.secondary_phage_id))


        #Accession = it will either be an accession or it will be "none"
        #Subcluster = it will either be a Subcluster or it will be "none"

        if self.annotation_author != '1' and self.annotation_author != '0':
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Author field." % self.primary_phage_id)

        if self.run_mode == 'none':
            self.set_evaluation("error", \
                "Phage %s does not have correctly populated Run Mode field." % self.primary_phage_id)









    # Group of all functions to evaluate the structure of the ticket.
    def validate_ticket(self):

        self.check_retrieve_status()
        self.check_type()
        self.check_host()
        self.check_cluster_subcluster()
        self.check_status()
        self.check_description_field()
        self.check_accession()
        self.check_annotation_author()
        self.check_run_mode()
        self.check_ticket()






###
