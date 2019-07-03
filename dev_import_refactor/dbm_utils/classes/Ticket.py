"""Represents a structure to contain directions for how to parse and import
genomes into Phamerator."""

from classes import Eval




class ImportTicket:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.type = "" # Add, Replace, Remove, UPDATE
        self.primary_phage_id = "" # Genome that will be added or updated
        self.host = ""
        self.cluster = "" # Singleton should be reported as 'singleton'
        self.subcluster = ""
        self.status = ""
        self.annotation_author = ""
        self.description_field = ""
        self.accession = ""
        self.run_mode = ""
        self.secondary_phage_id = "" # Genome that will be removed or replaced




        # Initialize calculated attributes
        self.match_strategy = "" # phage_id or filename
        self.evaluations = []




    #
    # def set_evaluation(self, type, message1 = None, message2 = None):
    #
    #     if type == "warning":
    #         eval_object = Eval.construct_warning(message1, message2)
    #
    #     elif type == "error":
    #         eval_object = Eval.construct_error(message1)
    #
    #     else:
    #         eval_object = Eval.EvalResult()
    #
    #     self.evaluations.append(eval_object)





	# Make sure "none" and "retrieve" indications are lowercase, as well as
    # "action", "status", "feature", and "author" fields are lowercase.


    def set_case(self):

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

    def clear_retrieve_status(self):

        #Host
        if self.host == "retrieve":
            self.host = "none"
            result = "Host data was not retrieved from Phagesdb."
            status = "error"
        else:
            result = "Host data is populated."
            status = "correct"

        definition = "Confirm host data is populated."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)





        #Cluster
        if self.cluster == "retrieve":
            self.cluster = "none"
            result = "Cluster data was not retrieved from Phagesdb."
            status = "error"
        else:
            result = "Cluster data is populated."
            status = "correct"

        definition = "Confirm Cluster data is populated."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)





        #Subcluster
        if self.subcluster == "retrieve":
            self.subcluster = "none"
            result = "Subcluster data was not retrieved from Phagesdb.")
            status = "error"
        else:
            result = "Subcluster data is populated."
            status = "correct"

        definition = "Confirm Subcluster data is populated."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


        #Accession
        if self.accession == "retrieve":
            self.accession = "none"
            result = "Accession data was not retrieved from Phagesdb."
            status = "error"
        else:
            result = "Accession data is populated."
            status = "correct"

        definition = "Confirm Accession data is populated."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)





    # Make sure the requested action is permissible
    # This script currently only allows four actions:
    # add, remove, replace, and update.
    def check_type(self):

        if not (self.type == "add" or \
            self.type == "remove" or \
            self.type == "replace" or \
            self.type == "update"):

            result = "The ticket type %s is not valid." % self.type
            status = "error"
        else:
            result = "The ticket type is valid."
            status = "correct"

        definition = "Confirm validity of ticket type."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)





    # Modify fields if needed

    # Modify Host if needed
    def check_host(self, host_set):

        if self.host != "none":

            # Keep only the genus in the host data field and discard the rest
            self.host = self.host.split(' ')[0]

            # TODO need to implement it elsewhere?
            # TODO define host_set
            if self.host not in host_set:

                result = \
                    "The host strain %s is not currently in the database." % \
                    self.host
                status = "error"
            else:
                result = "The host strain is not new."
                status = "correct"

            definition = "Check whether host strain is new."
            eval = Eval.Eval(id = "TICKET", \
                            definition = definition, \
                            result = result, \
                            status = status)
            self.evaluations.append(eval)









    # Check Subcluster data
    def check_subcluster(self, phage_subcluster_set):

        if self.subcluster != "none":

            if self.subcluster not in phage_subcluster_set:

                result = \
                    "The Subcluster %s is not currently in the database." % \
                    self.subcluster
                status = "error"
            else:
                result = "The Subcluster is not new."
                status = "correct"

            definition = "Check whether Subcluster is new."
            eval = Eval.Eval(id = "TICKET", \
                            definition = definition, \
                            result = result, \
                            status = status)
            self.evaluations.append(eval)


            if len(self.subcluster) > 5:

                result = \
                    "The Subcluster designation %s exceeds the character limit." % \
                    self.subcluster
                status = "error"
            else:
                result = "The Subcluster designation length is valid."
                status = "correct"

            definition = "Check whether Subcluster designation length is valud."
            eval = Eval.Eval(id = "TICKET", \
                            definition = definition, \
                            result = result, \
                            status = status)
            self.evaluations.append(eval)






	# Modify Cluster and Subcluster if needed
    def check_cluster(self, phage_cluster_set):

        # Check Cluster data
        if self.cluster != "none":

            if self.cluster.lower() != "singleton":

                if self.cluster not in phage_cluster_set:

                    result = \
                        "The Cluster %s is not currently in the database." % \
                        self.cluster
                    status = "error"
                else:
                    result = "The Cluster designation is valid."
                    status = "correct"

                definition = "Check whether Cluster is new."
                eval = Eval.Eval(id = "TICKET", \
                                definition = definition, \
                                result = result, \
                                status = status)
                self.evaluations.append(eval)

                if len(self.cluster) > 5:
                    result = \
                        "The Cluster designation %s exceeds the character limit." % \
                        self.cluster
                    status = "error"
                else:
                    result = "The Cluster designation length is valid."
                    status = "correct"

                definition = "Check whether Cluster designation length is valid."
                eval = Eval.Eval(id = "TICKET", \
                                definition = definition, \
                                result = result, \
                                status = status)
                self.evaluations.append(eval)


            else:
                self.cluster = self.cluster.lower()





	# Compare Cluster and Subcluster if needed
    def check_cluster_subcluster(self):

        # If Singleton or Unknown Cluster, there should be no Subcluster
        if (self.cluster == "singleton" or \
            self.cluster == "UNK" or \
            self.cluster == "none"):

            if self.subcluster != "none":
                result = "There is a discrepancy between the " + \
                            "Cluster and Subcluster designations."
                status = "error"
            else:
                result = "The Cluster and Subcluster are valid."
                status = "correct"

        # If not Singleton or Unknown or none, then Cluster should be part
        # of Subcluster data and the remainder should be a digit
        elif self.subcluster != "none":

            if (self.subcluster[:len(self.cluster)] != self.cluster or \
                self.subcluster[len(self.cluster):].isdigit() == False):

                result = "There is a discrepancy between the " + \
                            "Cluster and Subcluster designations."
                status = "error"
            else:
                result = "The Cluster and Subcluster are valid."
                status = "correct"

        else:
            result = "The Cluster and Subcluster are valid."
            status = "correct"

        definition = "Check for discrepancy between Cluster and Subcluster."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)

        self.evaluations.append(eval)










	# Modify Status if needed
    def check_status(self, phage_status_set):

        if (self.status not in phage_status_set and self.status != "none"):

            result = \
                "The status %s is not currently in the database." % self.status
            status = "error"
        else:
            result = "The status is not new."
            status = "correct"

        definition = "Check whether status is new."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)

        if len(self.status) > 5:
            result = "The status designation %s exceeds the character limit." % \
                self.status
            status = "error"
        else:
            result = "The status designation length is valid."
            status = "correct"

        definition = "Check whether status designation length is valid."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


    # Modify Description Qualifier if needed
    def check_description_field(self, description_set):

        if (self.description_field not in description_set and \
            self.description_field != "none"):

            result = "The description field %s is not commonly used." % \
                    self.description_field
            status = "error"
        else:
            result = "The description field is common."
            status = "correct"

        definition = "Check whether the description field is commonly used."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)



    # Modify Accession if needed
    def check_accession(self):

        self.accession = self.accession.strip()
        self.accession = self.accession.split('.')[0]

        if self.accession == "":
            self.accession = "none"


    # TODO improve this method = create separate, general function to convert
    # author value same way as with strand.
    # TODO the logic for creating Eval needs to be improved, so that
    # a 'correct' Eval is also created.
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
            self.annotation_author = "error"
            result = "The annotation author designation %s is not correct." \
                    % self.annotation_author
            status = "error"
            definition = "Check whether annotation_author is valid."
            eval = Eval.Eval(id = "TICKET", \
                            definition = definition, \
                            result = result, \
                            status = status)
            self.evaluations.append(eval)


	# Make sure run mode is permissible
    def check_run_mode(self,run_mode_options_dict):

        if self.run_mode == "custom":

            # TODO make this variable accessible or return a value
            # run_mode_custom_total += 1

            pass

        elif self.run_mode not in set(run_mode_options_dict.keys()):

            result = "The run mode %s is not valid." % self.run_mode
            status = "error"
        else:
            result = "The run mode is valid."
            status = "correct"

        definition = "Check whether the run mode is valid."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)








    def check_update_ticket(self,phage_id_set):

        if self.primary_phage_id not in phage_id_set:
            result = "The %s is not a valid PhageID in the database." \
                    % self.primary_phage_id
            status = "error"
        else:
            result = "The PhageID is valid."
            status = "correct"

        definition = "Check whether the primary PhageID is valid."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)

        if self.host == "none":
            result = "Phage %s does not have correctly populated Host field." \
                    % self.primary_phage_id
            status = "error"
        else:
            result = "The Host is valid."
            status = "correct"

        definition = "Check whether the Host is valid."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)










        if self.cluster == "none":
            result = "Phage %s does not have correctly populated Cluster field." \
                    % self.primary_phage_id
            status = "error"
        else:
            result = "The Cluster field is valid."
            status = "correct"

        definition = "Check whether the Cluster field is valid."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)




        if self.status == "none":
            result = "Phage %s does not have correctly populated Status field." \
                    % self.primary_phage_id
            status = "error"
        else:
            result = "The Status field is valid."
            status = "correct"

        definition = "Check whether the Status field is valid."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


        if self.description_field != "none":
            result = "Phage %s does not have correctly populated Description field." \
                    % self.primary_phage_id
            status = "error"
        else:
            result = "The Description field is valid."
            status = "correct"

        definition = "Check whether the Description field is valid."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)






        if self.secondary_phage_id != "none":
            result = "Phage %s should not have a genome listed to be removed." \
                    % self.primary_phage_id
            status = "error"
        else:
            result = "The Secondary PhageID field is valid."
            status = "correct"

        definition = "Check whether the Secondary PhageID field is valid."
        eval = Eval.Eval(id = "TICKET", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)


        #Accession = it will either be an accession or it will be "none"
        #Subcluster = it will either be a Subcluster or it will be "none"

        if (self.annotation_author != "1" and self.annotation_author != "0"):
            result = "Phage %s does not have correctly populated Author field." \
                    % self.primary_phage_id
            status = "error"

        if self.run_mode != "none":
            result = "Phage %s does not have correctly populated Run Mode field." \
                    % self.primary_phage_id
            status = "error"



    def check_add_ticket(self,phage_id_set):

        if self.primary_phage_id in phage_id_set:

            result = "Phage %s is already a PhageID in the database." + \
                    "This genome cannot be added to the database." \
                    % self.primary_phage_id
            status = "error"

        if self.primary_phage_id == "none":
            result = "Phage %s does not have correctly populated Primary PhageID field." \
                    % self.primary_phage_id
            status = "error"

        if self.host == "none":
            result = "Phage %s does not have correctly populated Host field." \
                    % self.primary_phage_id
            status = "error"

        if self.cluster == "none":
            result = "Phage %s does not have correctly populated Cluster field." \
                    % self.primary_phage_id
            status = "error"

        if self.status == "none":
            result = "Phage %s does not have correctly populated Status field." \
                    % self.primary_phage_id
            status = "error"

        if self.status == "final":
            result = "The phage %s to be added is listed " + \
                    "as Final status, but no Draft (or other) genome " + \
                    " is listed to be removed."
            status = "error"

        if self.description_field == "none":
            result = "Phage %s does not have correctly populated Description field." \
                    % self.primary_phage_id
            status = "error"

        if self.secondary_phage_id != "none":
            result = "Phage %s to be added should not have a genome indicated for removal." \
                    % self.primary_phage_id
            status = "error"

        #Accession = it will either be an accession or it will be "none"
        #Subcluster = it will either be a Subcluster or it will be "none"

        if self.annotation_author != '1' and self.annotation_author != '0':
            result = "Phage %s does not have correctly populated Author field." \
                    % self.primary_phage_id
            status = "error"

        if self.run_mode == 'none':
            result = "Phage %s does not have correctly populated Run Mode field." \
                    % self.primary_phage_id
            status = "error"




    def check_remove_ticket(self,phage_id_set):
        error_msg4 = "Secondary Phage %s does not have correctly " + \
        "populated %s field."

        if self.primary_phage_id != "none":
            result = error_msg4 % (self.secondary_phage_id, "Primary PhageID")
            status = "error"

        if self.host != "none":
            result =  error_msg4 % (self.secondary_phage_id, "Host")
            status = "error"

        if self.cluster != "none":
            result =  error_msg4 % (self.secondary_phage_id, "Cluster")
            status = "error"

        if self.subcluster != "none":
            result =  error_msg4 % (self.secondary_phage_id, "Subcluster")
            status = "error"

        if self.status != "none":
            result =  error_msg4 % (self.secondary_phage_id, "Status")
            status = "error"

        if self.description_field != "none":
            result =  error_msg4 % (self.secondary_phage_id, "Description Field")
            status = "error"

        if self.accession != "none":
            result =  error_msg4 % (self.secondary_phage_id, "Accession")
            status = "error"

        if self.annotation_author != "none":
            result =  error_msg4 % (self.secondary_phage_id, "AnnotationAuthor")
            status = "error"

        if self.run_mode != "none":
            result =  error_msg4 % (self.secondary_phage_id, "Run Mode")
            status = "error"

        if self.secondary_phage_id not in phage_id_set:
            result = "Secondary Phage %s is not a valid PhageID. " + \
                    "This genome cannot be dropped from the database."
                    % self.secondary_phage_id
            status = "error"




    def check_replace_ticket(self,phage_id_set):

        if self.primary_phage_id == "none":

            result = "Primary Phage %s is not a valid PhageID. %s genome cannot be replaced."
                % (self.primary_phage_id,self.secondary_phage_id)
            status = "error"

        if self.primary_phage_id != self.secondary_phage_id:

            result = \
                "The Primary Phage %s and Secondary phage %s are not spelled the same."
                % (self.primary_phage_id,self.secondary_phage_id
            status = "error"

            #FirstPhageID. If replacing a genome, ensure that if the genome to
            #be removed is not the same, that the new genome added has a unique name
            if self.primary_phage_id in phage_id_set:
                result = "Primary Phage %s is already a PhageID " + \
                    "in the database. This genome cannot be added to the database."
                    % self.primary_phage_id
                status = "error"

        if self.host == "none":
            result = "Phage %s does not have correctly populated Host field." \
                % self.primary_phage_id
            status = "error"

        if self.cluster == "none":
            result = "Phage %s does not have correctly populated Cluster field." \
                % self.primary_phage_id
            status = "error"

        if self.status == "none":
            result = "Phage %s does not have correctly populated Status field." \
                % self.primary_phage_id
            status = "error"

        if self.description_field == "none":
            result = "Phage %s does not have correctly populated Description Field field." \
                % self.primary_phage_id
            status = "error"

        if self.secondary_phage_id not in phage_id_set:
            result = "Secondary Phage %s is not a valid PhageID. " + \
                "This genome cannot be dropped from the database."
                 % self.secondary_phage_id
            status = "error"



        #Accession = it will either be an accession or it will be "none"
        #Subcluster = it will either be a Subcluster or it will be "none"

        if self.annotation_author != '1' and self.annotation_author != '0':
            result = "Phage %s does not have correctly populated Author field." \
                % self.primary_phage_id
            status = "error"

        if self.run_mode == 'none':
            result = "Phage %s does not have correctly populated Run Mode field." \
                % self.primary_phage_id
            status = "error"







    # Rules for how each field is populated differs depending on
    # each specific action.
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
