"""Represents a collection of data about a Source feature that is commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

from pdm_utils.functions import basic
from pdm_utils.classes import eval

class Source:

    def __init__(self):

        self.id = ""
        self.name = ""
        self.seqfeature = None # Biopython SeqFeature object.
        self.start = -1
        self.stop = -1
        self.organism = ""
        self.host = ""
        self.lab_host = ""

        # Data about the parent genome.
        self.genome_id = ""

        self._organism_name = "" # Parsed from organism.
        self._organism_host_genus = "" # Parsed from organism.
        self._host_host_genus = "" # Parsed from host.
        self._lab_host_host_genus = "" # Parsed from lab_host.
        self.evaluations = [] # List of warnings and errors about source feature.
        self.type = ""


    def parse_organism(self):
        """Retrieve the phage and host_genus names from the 'organism' field."""
        output = basic.parse_names_from_record_field(self.organism)
        self._organism_name = output[0]
        self._organism_host_genus = output[1]

    def parse_host(self):
        """Retrieve the host_genus name from the 'host' field."""
        self._host_host_genus = self.host.split(" ")[0]
        # Note: 'host' field is only expected to contain host information,
        # so a different strategy than parse_organism is used.

    def parse_lab_host(self):
        """Retrieve the host_genus name from the 'lab_host' field."""
        self._lab_host_host_genus = self.lab_host.split(" ")[0]
        # Note: 'host' field is only expected to contain host information,
        # so a different strategy than parse_organism is used.


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
        """Check that the attribute value is valid.

        :param attribute: Name of the Source feature object attribute to evaluate.
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
            value1_short = basic.truncate_value(str(value1), 30, "...")
            result = f"The {attribute} value '{value1_short}' is "
            value2 = basic.check_value_expected_in_set(
                        value1, check_set, expect)
            if value2:
                result = result + "valid."
                status = success
            else:
                result = result + "not valid."
                status = fail
        else:
            result = f"'{attribute}' is not a valid attribute to be evaluated."
            status = "untested"
        definition = f"Check the value of the '{attribute}' attribute."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)
