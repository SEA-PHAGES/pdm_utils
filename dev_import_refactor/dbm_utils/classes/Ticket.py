"""Represents a structure to contain directions for how to parse and import
genomes into Phamerator."""

from classes import Eval
from functions import basic




class GenomeTicket:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.id = "" # Unique identifier
        self.type = "" # Add, Replace, Remove, UPDATE

        # Attribute used to populate Genome objects for any ticket type.
        self.primary_phage_id = ""

        # Attribute used to evaluate all 'add' and 'replace' ticket types.
        self.run_mode = ""
        self.description_field = ""


        # Attributes used to populate Genome objects for
        # 'update', 'add', and 'replace' ticket types.
        self.host_genus = ""
        self.cluster = ""
        self.subcluster = ""
        self.annotation_status = ""
        self.annotation_author = ""
        self.annotation_qc = ""
        self.retrieve_record = ""
        self.accession = ""


        # TODO this attribute can probably be deleted
        # Attribute used to populate Genome objects for 'replace' ticket types.
        self.secondary_phage_id = ""




        self.evaluations = []
        self._parsed_fields = 0



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

    def set_primary_phage_id(self, value):
        """Set the primary_phage_id."""
        self.primary_phage_id = basic.lower_case(value)

    def set_host(self, value):
        """Set the host_genus."""
        self.host_genus = basic.lower_case(value)

    def set_cluster(self, value):
        """Set the cluster."""
        self.cluster = basic.lower_case(value)

    def set_subcluster(self, value):
        """Set the subcluster."""
        self.subcluster = basic.lower_case(value)

    def set_annotation_status(self, value):
        """Set the annotation_status."""
        self.annotation_status = value.lower()

    def set_accession(self, value):
        """Set the accession."""
        self.accession = basic.lower_case(value)

    def set_annotation_author(self, value):
        """Set the annotation_author."""
        self.annotation_author = value.lower()

    def set_annotation_qc(self, value):
        """Set the set_annotation_qc."""
        self.annotation_qc = basic.lower_case(value)

    def set_retrieve_record(self, value):
        """Set the set_retrieve_record."""
        self.retrieve_record = basic.lower_case(value)

    def set_secondary_phage_id(self, value):
        """Set the secondary_phage_id."""
        self.secondary_phage_id = basic.lower_case(value)






    # Evaluations


    def check_parsed_fields(self):
        """Check if the import table contained the correct amount of data."""

        if self._parsed_fields == 11:
            result = "Import table contains the correct amount of data."
            status = "correct"
        else:
            result = "Import table contains incorrect amount of data."
            status = "error"

        definition = "Check if import table contains the correct amount of data."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_type(self, value_set, expected = True):
        """Check if the type is populated with a value from
        a specific set of possible values."""

        output = basic.check_value_expected_in_set(
                    self.type, value_set, expected)

        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if ticket type field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_primary_phage_id(self, value_set, expected = True):
        """Check if the primary_phage_id is populated with an empty value.
        The provided set contains a list of possible empty values."""

        output = basic.check_value_expected_in_set(
                    self.primary_phage_id, value_set, expected)

        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if primary_phage_id field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_secondary_phage_id(self, value_set, expected = True):
        """Check if the secondary_phage_id is populated with an empty value.
        The provided set contains a list of possible empty values."""

        output = basic.check_value_expected_in_set(
                    self.secondary_phage_id, value_set, expected)

        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if secondary_phage_id field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_host_genus(self, value_set, expected = True):
        """Check if the host_genus is populated with an empty value.
        The provided set contains a list of possible empty values."""

        output = basic.check_value_expected_in_set(
                    self.host_genus, value_set, expected)

        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if host_genus field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_subcluster(self, value_set, expected = True):
        """Check if the subcluster is populated with an empty value.
        The provided set contains a list of possible empty values."""

        output = basic.check_value_expected_in_set(
                    self.subcluster, value_set, expected)

        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if subcluster field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_cluster(self, value_set, expected = True):
        """Check if the cluster is populated with an empty value.
        The provided set contains a list of possible empty values."""

        output = basic.check_value_expected_in_set(
                    self.cluster, value_set, expected)

        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if cluster field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_annotation_status(self, value_set, expected = True):
        """Check if the annotation_status is populated with an empty value.
        The provided set contains a list of possible empty values."""

        output = basic.check_value_expected_in_set(
                    self.annotation_status, value_set, expected)

        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if annotation_status field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_description_field(self, value_set, expected = True):
        """Check if the description_field is populated with an empty value.
        The provided set contains a list of possible empty values."""

        output = basic.check_value_expected_in_set(
                    self.description_field, value_set, expected)

        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if description_field field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_accession(self, value_set, expected = True):
        """Check if the accession is populated with an empty value.
        The provided set contains a list of possible empty values."""

        output = basic.check_value_expected_in_set(
                    self.accession, value_set, expected)

        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if accession field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_annotation_author(self, value_set, expected = True):
        """Check if the annotation_author is populated with an empty value.
        The provided set contains a list of possible empty values."""

        output = basic.check_value_expected_in_set(
                    self.annotation_author, value_set, expected)

        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if annotation_author field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_run_mode(self, value_set, expected = True):
        """Check if the run_mode is populated with an empty value.
        The provided set contains a list of possible empty values."""

        output = basic.check_value_expected_in_set(
                    self.run_mode, value_set, expected)

        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if run_mode field is correctly populated."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_duplicate_primary_phage_id(self, set_of_duplicates):
        """Check if the primary_phage_id is unique to this ticket by
        checking if it is found within a list of previously
        determined duplicate phage_ids."""

        if self.primary_phage_id in set_of_duplicates:
            result = "The primary_phage_id is not unique to this ticket."
            status = "error"
        else:
            result = "The primary_phage_id is unique to this ticket"
            status = "correct"

        definition = "Check if the primary_phage_id is unique to this ticket."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_duplicate_secondary_phage_id(self, set_of_duplicates):
        """Check if the primary_phage_id is unique to this ticket by
        checking if it is found within a list of previously
        determined duplicate phage_ids."""

        if self.secondary_phage_id in set_of_duplicates:
            result = "The secondary_phage_id is not unique to this ticket."
            status = "error"
        else:
            result = "The secondary_phage_id is unique to this ticket"
            status = "correct"

        definition = "Check if the secondary_phage_id is unique to this ticket."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)


    def check_duplicate_accession(self, set_of_duplicates):
        """Check if the accession is unique to this ticket by
        checking if it is found within a list of previously
        determined duplicate accessions."""

        if self.accession in set_of_duplicates:
            result = "The accession is not unique to this ticket."
            status = "error"
        else:
            result = "The accession is unique to this ticket"
            status = "correct"

        definition = "Check if the accession is unique to this ticket."
        eval = Eval.Eval("TICKET", definition, result, status)
        self.evaluations.append(eval)















###
