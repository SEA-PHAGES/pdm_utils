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
        self.left = -1
        self.right = -1
        self.organism = ""
        self.host = ""
        self.lab_host = ""

        # Data about the parent genome.
        self.genome_id = ""
        self.genome_host_genus = ""

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
    def check_attribute(self, attribute, check_set, expect=False, eval_id=None):
        """Check that the id is valid.

        :param attribute: Name of the source feature object attribute to evaluate.
        :type attribute: str
        :param check_set:
            Set of reference values.
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



    # TODO this is probably no longer needed.
    # def check_organism_name(self, eval_id=None):
    #     """Check phage name spelling in the organism field.
    #
    #     :param eval_id: Unique identifier for the evaluation.
    #     :type eval_id: str
    #     """
    #
    #     if self.genome_id != self._organism_name:
    #         result = ("The phage name in the organism field "
    #                     "does not match the genome_id.")
    #         status = "error"
    #
    #     else:
    #         result = "The phage name is spelled correctly."
    #         status = "correct"
    #
    #     definition = "Check phage name spelling in the organism field."
    #     evl = eval.Eval(eval_id, definition, result, status)
    #     self.evaluations.append(evl)


    # TODO this is probably no longer needed.
    # def check_organism_host_genus(self, eval_id=None):
    #     """Check host_genus name spelling in the organism field.
    #
    #     :param eval_id: Unique identifier for the evaluation.
    #     :type eval_id: str
    #     """
    #
    #     if self.genome_host_genus != self._organism_host_genus:
    #         result = ("The host_genus name in the organism field "
    #                     "does not match the genome_host_genus.")
    #         status = "error"
    #
    #     else:
    #         result = "The host_genus name is spelled correctly."
    #         status = "correct"
    #
    #     definition = "Check host_genus name spelling in the organism field."
    #     evl = eval.Eval(eval_id, definition, result, status)
    #     self.evaluations.append(evl)


    # TODO this is probably no longer needed.
    # def check_host_host_genus(self, eval_id=None):
    #     """Check host_genus name spelling in the host field.
    #
    #     :param eval_id: Unique identifier for the evaluation.
    #     :type eval_id: str
    #     """
    #
    #     if self.genome_host_genus != self._host_host_genus:
    #         result = ("The host_genus name in the host field "
    #                     "does not match the genome_host_genus.")
    #         status = "error"
    #
    #     else:
    #         result = "The host_genus name is spelled correctly."
    #         status = "correct"
    #
    #     definition = "Check host_genus name spelling in the host field."
    #     evl = eval.Eval(eval_id, definition, result, status)
    #     self.evaluations.append(evl)


    # TODO this is probably no longer needed.
    # def check_lab_host_host_genus(self, eval_id=None):
    #     """Check host_genus name spelling in the lab_host field.
    #
    #     :param eval_id: Unique identifier for the evaluation.
    #     :type eval_id: str
    #     """
    #
    #     if self.genome_host_genus != self._lab_host_host_genus:
    #         result = ("The host_genus name in the lab_host field "
    #                     "does not match the genome_host_genus.")
    #         status = "error"
    #
    #     else:
    #         result = "The host_genus name is spelled correctly."
    #         status = "correct"
    #
    #     definition = "Check host_genus name spelling in the lab_host field."
    #     evl = eval.Eval(eval_id, definition, result, status)
    #     self.evaluations.append(evl)
