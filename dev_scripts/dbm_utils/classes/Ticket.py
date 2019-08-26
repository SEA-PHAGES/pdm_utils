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
        self.phage_id = ""

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




        self.evaluations = []
        self._parsed_fields = 0
        self._value_flag = False


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

    def set_phage_id(self, value):
        """Set the phage_id."""
        self.phage_id = basic.lower_case(value)

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

    def set_value_flag(self, value):
        """Sets the flag if any attributes contain the specified 'value'."""
        if value in vars(self).values():
            self._value_flag = True
        else:
            self._value_flag = False




    # Evaluations

    def check_parsed_fields(self, eval_id=None):
        """Check if the import table contained the correct amount of data."""

        if self._parsed_fields == 11:
            result = "Import table contains the correct amount of data."
            status = "correct"
        else:
            result = "Import table contains incorrect amount of data."
            status = "error"

        definition = "Check if import table contains the correct amount of data."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_type(self, value_set, expected=True, eval_id=None):
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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_phage_id(self, value_set, expected=True, eval_id=None):
        """Check if the phage_id is populated with an empty value.
        The provided set contains a list of possible empty values."""

        output = basic.check_value_expected_in_set(
                    self.phage_id, value_set, expected)

        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"

        definition = "Check if phage_id field is correctly populated."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_host_genus(self, value_set, expected=True, eval_id=None):
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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_subcluster(self, value_set, expected=True, eval_id=None):
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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_cluster(self, value_set, expected=True, eval_id=None):
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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_annotation_status(self, value_set, expected=True, eval_id=None):
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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_description_field(self, value_set, expected=True, eval_id=None):
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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_accession(self, value_set, expected=True, eval_id=None):
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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_annotation_author(self, value_set, expected=True, eval_id=None):
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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_run_mode(self, value_set, expected=True, eval_id=None):
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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_duplicate_phage_id(self, set_of_duplicates, eval_id=None):
        """Check if the phage_id is unique to this ticket by
        checking if it is found within a list of previously
        determined duplicate phage_ids."""

        if self.phage_id in set_of_duplicates:
            result = "The phage_id is not unique to this ticket."
            status = "error"
        else:
            result = "The phage_id is unique to this ticket"
            status = "correct"

        definition = "Check if the phage_id is unique to this ticket."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_duplicate_accession(self, set_of_duplicates, eval_id=None):
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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)


    def check_compatible_type_and_annotation_status(self, eval_id=None):
        """Check if the ticket type and annotation_status are compatible.

        If the ticket type is "add", then the annotation_status is not
        expected to be "final".
        If the ticket type is "replace", then the annotation_status is
        not expected to be "draft".
        """

        if (self.type == "add" and self.annotation_status == "final"):
            result = "The ticket type indicates that a genome" + \
                     "with 'final' annotation_status will be added," + \
                     "which is not expected."
            status = "error"
        elif (self.type == "replace" and self.annotation_status == "draft"):
            result = "The ticket type indicates that a genome" + \
                     "with 'draft' annotation_status will be replaced," + \
                     "which is not expected."
            status = "error"
        else:
            result = "The ticket type and annotation_status are expected."
            status = "correct"
        definition = "Check if the ticket type and annotation_status" \
                     + " are compatible."
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)












###
