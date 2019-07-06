"""Represents a structure to contain directions for how to parse and import
genomes into Phamerator."""

from classes import Eval
from functions import basic




class GenomeTicket:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.type = "" # Add, Replace, Remove, UPDATE
        self.run_mode = ""

        # Attribute used to populate Genome objects for any ticket type.
        self.primary_phage_id = ""

        # Attributes used to populate Genome objects for
        # 'update', 'add', and 'replace' ticket types.
        self.host = ""
        self.cluster = "" # Singleton should be reported as 'singleton'
        self.subcluster = ""
        self.status = ""
        self.annotation_author = ""
        self.description_field = ""
        self.accession = ""

        # Attribute used to populate Genome objects for 'replace' ticket types.
        self.secondary_phage_id = ""




        # Initialize calculated attributes
        self.evaluations = []



        # TODO this attribute may no longer be needed.
        self.match_strategy = "" # phage_id or filename









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





    # TODO this method might not be necessary. It could be implemented
    # in the Genome object.
    # TODO unit test.
    def set_annotation_author(self,value):
        """Convert author name listed in ticket to binary value if needed."""
        self.annotation_author = basic.convert_author(value)






    # If either the Host, Cluster, Subcluster or Accession data needs to be
    # retrieved, try to access the data in phagesdb before proceeding

    def clear_retrieve_status(self):

        #Host
        if self.host == "retrieve":
            self.host = "none"
            result1 = "Host data was not retrieved from Phagesdb."
            status1 = "error"
        else:
            result1 = "Host data is populated."
            status1 = "correct"

        definition1 = "Confirm host data is populated."
        eval1 = Eval.Eval(id = "TICKET", \
                        definition = definition1, \
                        result = result1, \
                        status = status1)
        self.evaluations.append(eval1)


        #Cluster
        if self.cluster == "retrieve":
            self.cluster = "none"
            result2 = "Cluster data was not retrieved from Phagesdb."
            status2 = "error"
        else:
            result2 = "Cluster data is populated."
            status2 = "correct"

        definition2 = "Confirm Cluster data is populated."
        eval2 = Eval.Eval(id = "TICKET", \
                        definition = definition2, \
                        result = result2, \
                        status = status2)
        self.evaluations.append(eval2)





        #Subcluster
        if self.subcluster == "retrieve":
            self.subcluster = "none"
            result3 = "Subcluster data was not retrieved from Phagesdb."
            status3 = "error"
        else:
            result3 = "Subcluster data is populated."
            status3 = "correct"

        definition3 = "Confirm Subcluster data is populated."
        eval3 = Eval.Eval(id = "TICKET", \
                        definition = definition3, \
                        result = result3, \
                        status = status3)
        self.evaluations.append(eval3)


        #Accession
        if self.accession == "retrieve":
            self.accession = "none"
            result4 = "Accession data was not retrieved from Phagesdb."
            status4 = "error"
        else:
            result4 = "Accession data is populated."
            status4 = "correct"

        definition4 = "Confirm Accession data is populated."
        eval4 = Eval.Eval(id = "TICKET", \
                        definition = definition4, \
                        result = result4, \
                        status = status4)
        self.evaluations.append(eval4)












    # Evaluations


    def check_type(self, set1):
        """Check if the type is populated with an
        expected value. Every ticket type is expected to have
        a non-empty value, so provide a set of possible values."""

        if self.type in set1:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if ticket type field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_primary_phage_id(self, set1, set2, expect):
        """Check if the primary_phage_id is populated with an
        expected value. Provide two sets of mutually exclusive values
        The first set should contain all current PhageIDs in the database.
        The second set should contain values representing NULL or 'empty'.
        The expect value indicates whether the field is
        expected to be within the first, second, or neither set."""

        output = basic.check_value_in_two_sets(\
                    self.primary_phage_id, set1, set2)
        if output == expect:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if primary_phage_id field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_secondary_phage_id(self, set1, set2, expect):
        """Check if the secondary_phage_id is populated with an
        expected value. Provide two sets of mutually exclusive values
        The first set should contain all current PhageIDs in the database.
        The second set should contain values representing NULL or 'empty'.
        The expect value indicates whether the field is
        expected to be within the first, second, or neither set."""

        output = basic.check_value_in_two_sets(\
                    self.secondary_phage_id, set1, set2)
        if output == expect:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if secondary_phage_id field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_host(self, set1, set2, expect):
        """Check if the host is populated with an
        expected value. Provide two sets of mutually exclusive values
        The first set should contain all current PhageIDs in the database.
        The second set should contain values representing NULL or 'empty'.
        The expect value indicates whether the field is
        expected to be within the first, second, or neither set."""

        output = basic.check_value_in_two_sets(\
                    self.host, set1, set2)
        if output == expect:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if host field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)
















    def check_subcluster(self, set1, set2, expect):
        """Check if the subcluster is populated with an
        expected value. Provide two sets of mutually exclusive values
        The first set should contain all current PhageIDs in the database.
        The second set should contain values representing NULL or 'empty'.
        The expect value indicates whether the field is
        expected to be within the first, second, or neither set."""

        output = basic.check_value_in_two_sets(\
                    self.subcluster, set1, set2)
        if output == expect:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if subcluster field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)



















    def check_cluster(self, set1, set2, expect):
        """Check if the cluster is populated with an
        expected value. Provide two sets of mutually exclusive values
        The first set should contain all current PhageIDs in the database.
        The second set should contain values representing NULL or 'empty'.
        The expect value indicates whether the field is
        expected to be within the first, second, or neither set."""

        output = basic.check_value_in_two_sets(\
                    self.cluster, set1, set2)
        if output == expect:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if cluster field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)















    # TODO move this to Genome method?
    def check_subcluster_structure(self):
        """Check whether the cluster field is structured appropriately."""

        if self.subcluster != "none":
            left, right = basic.split_string(self.subcluster)

            if (left.isalpha() == False or right.isdigit() == False):
                result = "Subcluster is not structured correctly."
                status = "error"
            else:
                result = "Subcluster is structured correctly."
                status = "correct"

        else:
            result = "Subcluster is empty."
            status = "not_evaluated"

        definition = "Check if subcluster field is structured correctly."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    # TODO move this to Genome method?
    def check_cluster_structure(self):
        """Check whether the cluster field is structured appropriately."""

        if self.cluster != "none":
            left, right = basic.split_string(self.cluster)

            if (right != "" or left.isalpha() == False):
                result = "Cluster is not structured correctly."
                status = "error"
            else:
                result = "Cluster is structured correctly."
                status = "correct"
        else:
            result = "Cluster is empty."
            status = "not_evaluated"


        definition = "Check if cluster field is structured correctly."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    # TODO move this to Genome method?
    def check_cluster_subcluster(self):
        """Check whether the cluster and subcluster fields are
        compatible."""

        output = basic.compare_cluster_subcluster(self.cluster, self.subcluster)
        if not output:
            result = "Cluster and Subcluster designations are not compatible."
            status = "error"
        else:
            result = "Cluster and Subcluster designations are compatible."
            status = "correct"

        definition = "Check for compatibility between Cluster and Subcluster."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)















    def check_status(self, set1, set2, expect):
        """Check if the status is populated with an
        expected value. Provide two sets of mutually exclusive values
        The first set should contain all current PhageIDs in the database.
        The second set should contain values representing NULL or 'empty'.
        The expect value indicates whether the field is
        expected to be within the first, second, or neither set."""

        output = basic.check_value_in_two_sets(\
                    self.status, set1, set2)
        if output == expect:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if status field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)












    def check_description_field(self, set1, set2, expect):
        """Check if the description_field is populated with an
        expected value. Provide two sets of mutually exclusive values
        The first set should contain all current PhageIDs in the database.
        The second set should contain values representing NULL or 'empty'.
        The expect value indicates whether the field is
        expected to be within the first, second, or neither set."""

        output = basic.check_value_in_two_sets(\
                    self.description_field, set1, set2)
        if output == expect:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if description_field field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)












    def check_accession(self, set1, set2, expect):
        """Check if the accession is populated with an
        expected value. Provide two sets of mutually exclusive values
        The first set should contain all current PhageIDs in the database.
        The second set should contain values representing NULL or 'empty'.
        The expect value indicates whether the field is
        expected to be within the first, second, or neither set."""

        output = basic.check_value_in_two_sets(\
                    self.accession, set1, set2)
        if output == expect:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if accession field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)













    def check_annotation_author(self, set1, set2, expect):
        """Check if the annotation_author is populated with an
        expected value. Provide two sets of mutually exclusive values
        The first set should contain all current PhageIDs in the database.
        The second set should contain values representing NULL or 'empty'.
        The expect value indicates whether the field is
        expected to be within the first, second, or neither set."""

        output = basic.check_value_in_two_sets(\
                    self.annotation_author, set1, set2)
        if output == expect:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if annotation_author field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)















    def check_run_mode(self, set1, set2, expect):
        """Check if the run_mode is populated with an
        expected value. Provide two sets of mutually exclusive values
        The first set should contain all current PhageIDs in the database.
        The second set should contain values representing NULL or 'empty'.
        The expect value indicates whether the field is
        expected to be within the first, second, or neither set."""

        output = basic.check_value_in_two_sets(\
                    self.run_mode, set1, set2)
        if output == expect:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if run_mode field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)










    def check_primary_secondary_phage_ids(self):
        """Check whether the primary_phage_id and secondary_phage_id
        fields are the same."""

        if self.primary_phage_id != self.secondary_phage_id:
            result = "The primary_phage_id and secondary_phage_id are different."
            status = "error"
        else:
            result = "The primary_phage_id and secondary_phage_id are the same."
            status = "correct"

        definition = "Check if primary_phage_id and secondary_phage_id " + \
                        "fields are compatible."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)












###
