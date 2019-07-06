"""Represents a structure to contain directions for how to parse and import
genomes into Phamerator."""

from classes import Eval
from functions import basic




class GenomeTicket:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.type = "" # Add, Replace, Remove, UPDATE

        # Attribute used to populate Genome objects for any ticket type.
        self.primary_phage_id = ""

        # Attribute used to evaluate all 'add' and 'replace' ticket types.
        self.run_mode = ""
        self.description_field = ""


        # Attributes used to populate Genome objects for
        # 'update', 'add', and 'replace' ticket types.
        self.host = ""
        self.cluster = ""
        self.subcluster = ""
        self.status = ""
        self.annotation_author = ""
        self.accession = ""

        # Attribute used to populate Genome objects for 'replace' ticket types.
        self.secondary_phage_id = ""




        # Initialize calculated attributes
        self.evaluations = []



        # TODO this attribute may no longer be needed.
        self.match_strategy = "" # phage_id or filename












    # When parsing data from an import ticket, some fields should be
    # lowercase regardless of the content, while other fields should
    # be lowercase only if they are set to 'none' or 'retrieve'.
    def set_type(self, value):
        """Set the ticket type."""
        self.type = value.lower()

    def set_description_field(self, value):
        """Set the description_field."""
        self.description_field = value.lower()

    def set_run_mode(self, value):
        """Set the run_mode."""
        self.run_mode = value.lower()






    # TODO implement in Genome object?
    def set_primary_phage_id(self, value):
        """Set the primary_phage_id."""
        self.primary_phage_id = basic.lower_case(value)

    # TODO implement in Genome object?
    def set_host(self, value):
        """Set the host."""
        self.host = basic.lower_case(value)

    # TODO implement in Genome object?
    def set_cluster(self, value):
        """Set the cluster."""
        self.cluster = basic.lower_case(value)

    # TODO implement in Genome object?
    def set_subcluster(self, value):
        """Set the subcluster."""
        self.subcluster = basic.lower_case(value)

    # TODO implement in Genome object?
    def set_status(self, value):
        """Set the status."""
        self.status = value.lower()

    # TODO implement in Genome object?
    def set_accession(self, value):
        """Set the accession."""
        self.accession = basic.lower_case(value)

    # TODO implement in Genome object?
    def set_annotation_author(self, value):
        """Set the annotation_author."""
        self.annotation_author = value.lower()

    # TODO implement in Genome object?
    def set_secondary_phage_id(self, value):
        """Set the secondary_phage_id."""
        self.secondary_phage_id = basic.lower_case(value)






    # TODO below this should probably be moved to the Genome object.
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
    # TODO above this should probably be moved to the Genome object.












    # Evaluations


    def check_type(self, set1):
        """Check if the type is populated with a value from
        a specific set of possible values."""

        if self.type in set1:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if ticket type field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_primary_phage_id(self, set1):
        """Check if the primary_phage_id is populated with an empty value.
        The provided set contains a list of possible empty values."""

        if self.primary_phage_id not in set1:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if primary_phage_id field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_secondary_phage_id(self, set1):
        """Check if the secondary_phage_id is populated with an empty value.
        The provided set contains a list of possible empty values."""

        if self.secondary_phage_id not in set1:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if secondary_phage_id field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_host(self, set1):
        """Check if the host is populated with an empty value.
        The provided set contains a list of possible empty values."""

        if self.host not in set1:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if host field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_subcluster(self, set1):
        """Check if the subcluster is populated with an empty value.
        The provided set contains a list of possible empty values."""

        if self.subcluster not in set1:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if subcluster field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_cluster(self, set1):
        """Check if the cluster is populated with an empty value.
        The provided set contains a list of possible empty values."""

        if self.cluster not in set1:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if cluster field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_status(self, set1):
        """Check if the status is populated with an empty value.
        The provided set contains a list of possible empty values."""

        if self.status not in set1:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if status field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_description_field(self, set1):
        """Check if the description_field is populated with an empty value.
        The provided set contains a list of possible empty values."""

        if self.description_field not in set1:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if description_field field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_accession(self, set1):
        """Check if the accession is populated with an empty value.
        The provided set contains a list of possible empty values."""

        if self.accession not in set1:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if accession field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_annotation_author(self, set1):
        """Check if the annotation_author is populated with an empty value.
        The provided set contains a list of possible empty values."""

        if self.annotation_author not in set1:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if annotation_author field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_run_mode(self, set1):
        """Check if the run_mode is populated with an empty value.
        The provided set contains a list of possible empty values."""

        if self.run_mode not in set1:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if run_mode field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)





###
