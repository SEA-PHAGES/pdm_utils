"""Represents a structure to contain directions for how to parse and import
genomes into a MySQL database."""

from pdm_utils.classes import eval
from pdm_utils.functions import basic




class ImportTicket:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.id = None # Unique identifier
        self.type = "" # Add, Replace
        self.phage_id = ""
        self.eval_mode = ""
        self.description_field = ""
        self.eval_flags = {} # Dictionary of evaluation flags.

        # Attributes used to populate Genome objects for
        # 'update', 'add', and 'replace' ticket types.
        self.data_retrieve = set() # Data that should be retrieved from PhagesDB.
        self.data_retain = set() # Data that should be retained from the MySQL database.
        self.data_add = set() # Data to be added to genome from ticket.
        self.data_parse = set() # Data that should be parsed from the flat file.
        self.data_dict = {} # All data from import table.

        # Used to check the structure of the ticket data.
        self.evaluations = []



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

    def set_eval_mode(self, value):
        """Set the eval_mode.

        :param value: Value to be set as the eval_mode.
        :type value: str
        """
        self.eval_mode = value.lower()


    # Evaluations
    def set_eval(self, eval_id, definition, result, status):
        """Constructs and adds an Eval object to the evaluations list.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        :param definition: Description of the evaluation.
        :type definition: str
        :param result: Description of the outcome of the evaluation.
        :type result: str
        :param status: Outcome of the evaluation.
        :type status: str
        """
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_attribute(self, attribute, check_set, expect=False, eval_id=None,
                        success="correct", fail="error", eval_def=None):
        """Check that the id is valid.

        :param attribute: Name of the ImportTicket object attribute to evaluate.
        :type attribute: str
        :param check_set: Set of reference ids.
        :type check_set: set
        :param expect:
            Indicates whether the attribute value is expected to be present
            in the check set.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        :param success: Default status if the outcome is a success.
        :type success: str
        :param fail: Default status if the outcome is not a success.
        :type fail: str
        :param eval_def: Description of the evaluation.
        :type eval_def: str
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
                status = success
            else:
                result = f"The {attribute} is not valid."
                status = fail
        else:
            result = f"The {attribute} was not evaluated."
            status = "untested"
        definition = f"Check the {attribute} attribute."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_eval_flags(self, expect=True, eval_id=None,
                         success="correct", fail="error", eval_def=None):
        """Check that the eval_flags is valid.

        :param expect:
            Indicates whether the eval_flags is expected to contain data.
        :type expect: bool
        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        keys = len(self.eval_flags)
        msg = f"There are {keys} eval flags present, which is "
        if keys > 0 and expect:
            output = True
        elif keys == 0 and not expect:
            output = True
        else:
            output = False

        if output:
            result = msg + "expected."
            status = success
        else:
            result = msg + "not expected."
            status = fail
        definition = "Check if there are eval flags."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_compatible_type_and_data_retain(self, eval_id=None,
                          success="correct", fail="error", eval_def=None):
        """Check if the ticket type and data_retain are compatible.

        If the ticket type is 'add', then the data_retain set is not
        expected to have any data.

        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        msg = (f"The ticket type is '{self.type}' and "
               f"there are {len(self.data_retain)} values in the database "
               "that are set to be retained, which is ")
        if (self.type == "add" and len(self.data_retain) > 0):
            result = msg + "not expected."
            status = fail
        else:
            result = msg + "expected."
            status = success
        definition = ("Check if the ticket type and data_retain setting "
                      "are compatible.")
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def check_valid_data_source(self, ref_set_attr, check_set, eval_id=None,
                                success="correct", fail="error", eval_def=None):
        """Check that the values in the specified attribute are valid.

        :param ref_set_attr:
            Name of the data_dict in the ticket to be evaluated
            (data_add, data_retain, data_retrieve, data_parse)
        :type ref_set_attr: str
        :param check_set: Set of valid field names.
        :type check_set: set
        :param eval_id: same as for check_attribute().
        :param success: same as for check_attribute().
        :param fail: same as for check_attribute().
        :param eval_def: same as for check_attribute().
        """
        if ref_set_attr == "data_add":
            ref_set = self.data_add
        elif ref_set_attr == "data_retain":
            ref_set = self.data_retain
        elif ref_set_attr == "data_retrieve":
            ref_set = self.data_retrieve
        elif ref_set_attr == "data_parse":
            ref_set = self.data_parse
        else:
            ref_set = None
        if ref_set is not None:
            invalid_values = ref_set - check_set
            msg = f"The '{ref_set_attr}' field is ."
            if len(invalid_values) == 0:
                result = "populated correctly."
                status = success
            else:
                invalid_string = basic.join_strings(invalid_values,
                                                    delimiter=", ")
                result = (
                    "not populated correctly. The following "
                    f"values are not permitted in '{ref_set_attr}': "
                    f"{invalid_string}")
                status = fail
        else:
            result = f"'{ref_set_attr}' is not a valid attribute to be evaluated."
            status = fail

        definition = f"Check if {ref_set_attr} field is correctly populated."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)
