"""Represents a structure to contain directions for how to parse and import
genomes into Phamerator."""

from classes import Eval
from functions import basic




class GenomeTicket:

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


    # TODO unit test.
    def set_host(self,value):

        # Keep only the genus in the host data field and discard the rest
        self.host = value.split(' ')[0]





    # TODO unit test.
    def set_accession(self,value):

        self.accession = value.strip()
        self.accession = value.split('.')[0]

        if self.accession == "":
            self.accession = "none"



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
    def check_type(self, value_set, expected = True):
        """Check whether the ticket type field is populated as expected.
        The value_set contains a list of possible values for
        the ticket type. The 'expected' parameter indicates whether
        the ticket type is expected to be an element of the set or not."""

        result, status = basic.check_value_in_set(
                            self.type,
                            value_set,
                            expected)

        definition = "Check if ticket type field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)






    # TODO HERE IN PROGRESS. TRYING TO FIGURE OUT HOW TO TEST PHAGE_ID TO
    # BE WITHIN  A CERTAIN SET AS WELL AS NOT IN EMPTY SET.
    # TODO in progress, replacing old method.
    def check_primary_phage_id(self,
                                filled_set,
                                filled_expected = False,
                                empty_expected = False):
        """Check if primary_phage_id meets expectations."""

        filled_correct = basic.check_value_in_set2(
                                            self.primary_phage_id,
                                            filled_set,
                                            filled_expected)


        empty_correct = basic.check_value_in_set2(
                                            self.primary_phage_id,
                                            set(["none"]),
                                            empty_expected)


        if empty_correct:
            result = "The field is populated correctly."
            status = "correct"
        elif (empty_correct and filled_correct):
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if primary_phage_id field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)

    # TODO HERE IN PROGRESS. TRYING TO FIGURE OUT HOW TO TEST PHAGE_ID TO
    # BE WITHIN  A CERTAIN SET AS WELL AS NOT IN EMPTY SET.
    # TODO in progress, replacing old method.



















    def check_primary_phage_id(self, value_set, expected = True):
        """Check whether the primary_phage_id field is populated as expected.
        The value_set contains a list of possible values for
        the primary_phage_id. The 'expected' parameter indicates whether
        the primary_phage_id is expected to be an element of the set or not."""

        result, status = basic.check_value_in_set(
                            self.primary_phage_id,
                            value_set,
                            expected)

        definition = "Check if primary_phage_id field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_secondary_phage_id(self, value_set, expected = True):
        """Check whether the secondary_phage_id field is populated as expected.
        The value_set contains a list of possible values for
        the secondary_phage_id. The 'expected' parameter indicates whether
        the secondary_phage_id is expected to be an element of the set or not."""

        result, status = basic.check_value_in_set(
                            self.secondary_phage_id,
                            value_set,
                            expected)

        definition = "Check if secondary_phage_id field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_host(self, value_set, expected = True):
        """Check whether the host field is populated as expected.
        The value_set contains a list of possible values for
        the host. The 'expected' parameter indicates whether
        the host is expected to be an element of the set or not."""

        result, status = basic.check_value_in_set(
                            self.host,
                            value_set,
                            expected)

        definition = "Check if host field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)




    def check_subcluster(self, value_set, expected = True):
        """Check whether the subcluster field is populated as expected.
        The value_set contains a list of possible values for
        the subcluster. The 'expected' parameter indicates whether
        the subcluster is expected to be an element of the set or not."""

        result, status = basic.check_value_in_set(
                            self.subcluster,
                            value_set,
                            expected)

        definition = "Check if subcluster field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)




    def check_cluster(self, value_set, expected = True):
        """Check whether the cluster field is populated as expected.
        The value_set contains a list of possible values for
        the cluster. The 'expected' parameter indicates whether
        the cluster is expected to be an element of the set or not."""

        result, status = basic.check_value_in_set(
                            self.cluster,
                            value_set,
                            expected)

        definition = "Check if cluster field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


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


    def check_status(self, value_set, expected = True):
        """Check whether the status field is populated as expected.
        The value_set contains a list of possible values for
        the status. The 'expected' parameter indicates whether
        the status is expected to be an element of the set or not."""

        result, status = basic.check_value_in_set(
                            self.status,
                            value_set,
                            expected)

        definition = "Check if status field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_description_field(self, value_set, expected = True):
        """Check whether the description_field field is populated as expected.
        The value_set contains a list of possible values for
        the description_field. The 'expected' parameter indicates whether
        the description_field is expected to be an element of the set or not."""

        result, status = basic.check_value_in_set(
                            self.description_field,
                            value_set,
                            expected)

        definition = "Check if description_field field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_accession(self, value_set, expected = True):
        """Check whether the accession field is populated as expected.
        The value_set contains a list of possible values for
        the accession. The 'expected' parameter indicates whether
        the accession is expected to be an element of the set or not."""

        result, status = basic.check_value_in_set(
                            self.accession,
                            value_set,
                            expected)

        definition = "Check if accession field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_annotation_author(self, value_set, expected = True):
        """Check whether the annotation_author field is populated as expected.
        The value_set contains a list of possible values for
        the annotation_author. The 'expected' parameter indicates whether
        the annotation_author is expected to be an element of the set or not."""

        result, status = basic.check_value_in_set(
                            self.annotation_author,
                            value_set,
                            expected)

        definition = "Check if annotation_author field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_run_mode(self, value_set, expected = True):
        """Check whether the run_mode field is populated as expected.
        The value_set contains a list of possible values for
        the run_mode. The 'expected' parameter indicates whether
        the run_mode is expected to be an element of the set or not."""

        result, status = basic.check_value_in_set(
                            self.run_mode,
                            value_set,
                            expected)

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
