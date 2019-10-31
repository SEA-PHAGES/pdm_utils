"""Represents a structure to contain directions for how to parse and import
genomes into Phamerator."""

from pdm_utils.classes import eval
from pdm_utils.functions import basic




class GenomeTicket:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.id = None # Unique identifier
        self.type = "" # Add, Replace
        self.phage_id = ""
        self.run_mode = ""
        self.description_field = ""
        self.eval_flags = {} # Dictionary of evaluation flags.

        # Attributes used to populate Genome objects for
        # 'update', 'add', and 'replace' ticket types.
        self.data_retrieve = set() # Data that should be retrieved from PhagesDB.
        self.data_retain = set() # Data that should be retained from Phamerator.
        self.data_add = set() # Data to be added to genome from ticket.
        self.data_dict = {} # Original data from import table.

        # Used to check the structure of the ticket data.
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
    def check_attribute(self, attribute, check_set, expect=False, eval_id=None):
        """Check that the id is valid.

        :param attribute: Name of the ticket object attribute to evaluate.
        :type attribute: str
        :param check_set:
            Set of reference ids.
        :type check_set: set
        :param expect:
            Indicates whether the id is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id:
            Unique identifier for the evaluation.
        :type eval_id: str
        """
        try:
            test = True
            value1 = getattr(self, attribute)
        except:
            test = False
            value1 = None
        if test:
            value2 = basic.check_value_expected_in_set(
                        value1, check_set, expect)
            if value2:
                result = f"The {attribute} is valid."
                status = "correct"
            else:
                result = f"The {attribute} is not valid."
                status = "error"
        else:
            result = f"The {attribute} was not evaluated."
            status = "untested"
        definition = f"Check the {attribute} attribute."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_eval_flags(self, expect=True, eval_id=None):
        """Check that the eval_flags is valid.

        :param expect:
            Indicates whether the eval_flags is expected to contain data.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        keys = len(self.eval_flags)
        if keys > 0 and expect:
            output = True
        elif keys == 0 and not expect:
            output = True
        else:
            output = False

        if output:
            result = "The field is populated correctly."
            status = "correct"
        else:
            result = "The field is not populated correctly."
            status = "error"
        definition = "Check if eval_flags field is correctly populated."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_compatible_type_and_data_retain(self, eval_id=None):
        """Check if the ticket type and data_retain are compatible.

        If the ticket type is 'add', then the data_retain set is not
        expected to have any data.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        if (self.type == "add" and len(self.data_retain) > 0):
            result = ("The ticket type indicates that a genome "
                     "that will be added should also retain data, "
                     "from a genome in the database, which is not expected.")
            status = "error"
        else:
            result = "The ticket type and data_retain set are expected."
            status = "correct"
        definition = ("Check if the ticket type and data_retain "
                      "are compatible.")
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_valid_data_source(self, ref_set_attr, check_set, eval_id=None):
        """Check that the values in the specified attribute are valid.

        :param ref_set_attr:
            Name of the data_dict in the ticket to be evaluated
            (data_add, data_retain, data_retrieve)
        :type ref_set_attr: str
        :param check_set: Set of valid field names.
        :type check_set: set
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        if ref_set_attr == "data_add":
            ref_set = self.data_add
        elif ref_set_attr == "data_retain":
            ref_set = self.data_retain
        elif ref_set_attr == "data_retrieve":
            ref_set = self.data_retrieve
        else:
            ref_set = None
        if ref_set is not None:
            invalid_values = ref_set - check_set
            if len(invalid_values) == 0:
                result = "The field is populated correctly."
                status = "correct"
            else:
                result = (
                    "The field is not populated correctly. The following "
                    f"values are not permitted in '{ref_set_attr}': "
                    f"{list(invalid_values)}")
                status = "error"
        else:
            result = "Invalid field to be evaluated."
            status = "error"

        definition = f"Check if {ref_set_attr} field is correctly populated."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)
###
