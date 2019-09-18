"""Represents a structure to contain directions for how to parse and import
genomes into Phamerator."""

from pdm_utils.classes import eval
from pdm_utils.functions import basic




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
        self.eval_flags = {} # Dictionary of evaluation flags.

        # Attributes used to populate Genome objects for
        # 'update', 'add', and 'replace' ticket types.
        self.host_genus = ""
        self.cluster = ""
        self.subcluster = ""
        self.annotation_status = ""
        self.annotation_author = -1
        self.annotation_qc = -1
        self.retrieve_record = -1
        self.accession = ""


        self.evaluations = []
        self._value_flag = False








    # When parsing data from an import ticket, some fields should be
    # lowercase regardless of the content, while other fields should
    # be lowercase only if they are set to 'none' or 'retrieve'.
    def set_type(self, value):
        """Set the ticket type.

        :param value: Value to be set as the type.
        :type value: str
        """
        self.type = value.lower()

    def set_description_field(self, value):
        """Set the description_field.

        :param value: Value to be set as the description_field.
        :type value: str
        """
        self.description_field = value.lower()

    def set_run_mode(self, value):
        """Set the run_mode.

        :param value: Value to be set as the run_mode.
        :type value: str
        """
        self.run_mode = value.lower()

    def set_phage_id(self, value):
        """Set the phage_id.

        :param value: Value to be set as the phage_id.
        :type value: str
        """
        self.phage_id = basic.lower_case(value)

    def set_host(self, value):
        """Set the host_genus.

        :param value: Value to be set as the host_genus.
        :type value: str
        """
        self.host_genus = basic.lower_case(value)

    def set_cluster(self, value):
        """Set the cluster.

        :param value: Value to be set as the cluster.
        :type value: str
        """
        self.cluster = basic.lower_case(value)

    def set_subcluster(self, value):
        """Set the subcluster.

        :param value: Value to be set as the subcluster.
        :type value: str
        """
        self.subcluster = basic.lower_case(value)

    def set_annotation_status(self, value):
        """Set the annotation_status.

        :param value: Value to be set as the annotation_status.
        :type value: str
        """
        self.annotation_status = value.lower()

    def set_accession(self, value):
        """Set the accession.

        :param value: Value to be set as the accession.
        :type value: str
        """
        self.accession = basic.lower_case(value)

    def set_annotation_author(self, value):
        """Set the annotation_author.

        :param value: Value to be set as the annotation_author.
        :type value: str
        """
        self.annotation_author = basic.lower_case(value)

    def set_annotation_qc(self, value):
        """Set the set_annotation_qc.

        :param value: Value to be set as the annotation_qc.
        :type value: str
        """
        self.annotation_qc = basic.lower_case(value)

    def set_retrieve_record(self, value):
        """Set the set_retrieve_record.

        :param value: Value to be set as the retrieve_record.
        :type value: str
        """
        self.retrieve_record = basic.lower_case(value)

    def set_value_flag(self, value):
        """Sets the flag if any attributes contain the specified 'value'.

        :param value:
            Indicates the value that should be searched within
            the attributes.
        :type value: str
        """
        if value in vars(self).values():
            self._value_flag = True
        else:
            self._value_flag = False




    # Evaluations



    # TODO this is no longer needed.
    # def check_parsed_fields(self, eval_id=None):
    #     """Check if the import table contained the correct amount of data.
    #
    #     :param eval_id: Unique identifier for the evaluation.
    #     :type eval_id: str
    #     """
    #     if self._parsed_fields == 11:
    #         result = "Import table contains the correct amount of data."
    #         status = "correct"
    #     else:
    #         result = "Import table contains incorrect amount of data."
    #         status = "error"
    #
    #     definition = "Check if import table contains the correct amount of data."
    #     evl = eval.Eval(eval_id, definition, result, status)
    #     self.evaluations.append(evl)


    def check_type(self, check_set, expect=True, eval_id=None):
        """Check that the type is valid.

        :param check_set: Set of reference types.
        :type check_set: set
        :param expect:
            Indicates whether the type is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        output = basic.check_value_expected_in_set(
                    self.type, check_set, expect)
        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"
        definition = "Check if ticket type field is correctly populated."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_phage_id(self, check_set, expect=True, eval_id=None):
        """Check that the phage_id is valid.

        :param check_set: Set of reference phage_ids.
        :type check_set: set
        :param expect:
            Indicates whether the phage_id is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        output = basic.check_value_expected_in_set(
                    self.phage_id, check_set, expect)
        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"
        definition = "Check if phage_id field is correctly populated."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_host_genus(self, check_set, expect=True, eval_id=None):
        """Check that the host_genus is valid.

        :param check_set: Set of reference host_genus values.
        :type check_set: set
        :param expect:
            Indicates whether the host_genus is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        output = basic.check_value_expected_in_set(
                    self.host_genus, check_set, expect)
        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"
        definition = "Check if host_genus field is correctly populated."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_subcluster(self, check_set, expect=True, eval_id=None):
        """Check that the subcluster is valid.

        :param check_set: Set of reference subclusters.
        :type check_set: set
        :param expect:
            Indicates whether the subcluster is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        output = basic.check_value_expected_in_set(
                    self.subcluster, check_set, expect)
        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"
        definition = "Check if subcluster field is correctly populated."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_cluster(self, check_set, expect=True, eval_id=None):
        """Check that the cluster is valid.

        :param check_set: Set of reference clusters.
        :type check_set: set
        :param expect:
            Indicates whether the cluster is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        output = basic.check_value_expected_in_set(
                    self.cluster, check_set, expect)
        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"
        definition = "Check if cluster field is correctly populated."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_annotation_status(self, check_set, expect=True, eval_id=None):
        """Check that the annotation_status is valid.

        :param check_set: Set of reference annotation_status values.
        :type check_set: set
        :param expect:
            Indicates whether the annotation_status is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        output = basic.check_value_expected_in_set(
                    self.annotation_status, check_set, expect)
        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"
        definition = "Check if annotation_status field is correctly populated."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_description_field(self, check_set, expect=True, eval_id=None):
        """Check that the description_field is valid.

        :param check_set: Set of reference description_field values.
        :type check_set: set
        :param expect:
            Indicates whether the description_field is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        output = basic.check_value_expected_in_set(
                    self.description_field, check_set, expect)
        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"
        definition = "Check if description_field field is correctly populated."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_accession(self, check_set, expect=True, eval_id=None):
        """Check that the accession is valid.

        :param check_set: Set of reference accessions.
        :type check_set: set
        :param expect:
            Indicates whether the accession is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        output = basic.check_value_expected_in_set(
                    self.accession, check_set, expect)
        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"
        definition = "Check if accession field is correctly populated."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_annotation_author(self, check_set, expect=True, eval_id=None):
        """Check that the annotation_author is valid.

        :param check_set: Set of reference annotation_author values.
        :type check_set: set
        :param expect:
            Indicates whether the annotation_author is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        output = basic.check_value_expected_in_set(
                    self.annotation_author, check_set, expect)
        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"
        definition = "Check if annotation_author field is correctly populated."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_run_mode(self, check_set, expect=True, eval_id=None):
        """Check that the run_mode is valid.

        :param check_set: Set of reference run_mode values.
        :type check_set: set
        :param expect:
            Indicates whether the run_mode is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        output = basic.check_value_expected_in_set(
                    self.run_mode, check_set, expect)
        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"
        definition = "Check if run_mode field is correctly populated."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)

    def check_duplicate_id(self, set_of_duplicates, eval_id=None):
        """Check if the id is unique to this ticket by
        checking if it is found within a list of previously
        determined duplicate ids.

        :param set_of_duplicates: Set of reference duplicated values.
        :type set_of_duplicates: set
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        if self.id in set_of_duplicates:
            result = "The id is not unique to this ticket."
            status = "error"
        else:
            result = "The id is unique to this ticket"
            status = "correct"
        definition = "Check if the id is unique to this ticket."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)

    def check_duplicate_phage_id(self, set_of_duplicates, eval_id=None):
        """Check if the phage_id is unique to this ticket by
        checking if it is found within a list of previously
        determined duplicate phage_ids.

        :param set_of_duplicates: Set of reference duplicated values.
        :type set_of_duplicates: set
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        if self.phage_id in set_of_duplicates:
            result = "The phage_id is not unique to this ticket."
            status = "error"
        else:
            result = "The phage_id is unique to this ticket"
            status = "correct"
        definition = "Check if the phage_id is unique to this ticket."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_duplicate_accession(self, set_of_duplicates, eval_id=None):
        """Check if the accession is unique to this ticket by
        checking if it is found within a list of previously
        determined duplicate accessions.

        :param set_of_duplicates: Set of reference duplicated values.
        :type set_of_duplicates: set
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        if self.accession in set_of_duplicates:
            result = "The accession is not unique to this ticket."
            status = "error"
        else:
            result = "The accession is unique to this ticket"
            status = "correct"
        definition = "Check if the accession is unique to this ticket."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_compatible_type_and_annotation_status(self, eval_id=None):
        """Check if the ticket type and annotation_status are compatible.

        If the ticket type is 'add', then the annotation_status is not
        expected to be 'final'.
        If the ticket type is 'replace', then the annotation_status is
        not expected to be 'draft'.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
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
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_retrieve_record(self, check_set, expect=True, eval_id=None):
        """Check that the retrieve_record is valid.

        :param check_set: Set of reference accessions.
        :type check_set: set
        :param expect:
            Indicates whether retrieve_record is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        output = basic.check_value_expected_in_set(
                    self.retrieve_record, check_set, expect)
        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"
        definition = "Check if retrieve_record field is correctly populated."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)






###
