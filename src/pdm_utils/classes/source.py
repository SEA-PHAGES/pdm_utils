"""Represents a collection of data about a Source feature that is commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

from pdm_utils.functions import basic
from pdm_utils.classes import eval
import re




class Source:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.id = ""
        self.name = ""
        self.seqfeature = None
        self.type = ""
        self.left = -1 # TODO implement this.
        self.right = -1 # TODO implement this.
        self.organism = ""
        self.host = ""
        self.lab_host = ""

        # Common to Phamerator.
        self.genome_id = ""
        self.parent_host_genus = ""


        # Computed data fields.
        self._organism_name = ""
        self._organism_host_genus = ""
        self._host_host_genus = ""
        self._lab_host_host_genus = ""

        self.evaluations = [] # List of warnings and errors about source feature




    def parse_organism(self):
        """Retrieve the phage and host_genus names from the 'organism' field."""
        self._organism_name, \
        self._organism_host_genus = \
            basic.parse_names_from_record_field(self.organism)

    def parse_host(self):
        """Retrieve the host_genus name from the 'host' field."""
        name, \
        self._host_host_genus = \
            basic.parse_names_from_record_field(self.host)
        # Note: no need to assign phage name, since this field is only
        # expected to contain host information.

    def parse_lab_host(self):
        """Retrieve the host_genus name from the 'lab_host' field."""
        name, \
        self._lab_host_host_genus = \
            basic.parse_names_from_record_field(self.lab_host)
        # Note: no need to assign phage name, since this field is only
        # expected to contain host information.







    # Evalutions


    def check_organism_name(self, eval_id=None):
        """Check phage name spelling in the organism field.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """

        if self.genome_id != self._organism_name:
            result = "The phage name in the organism field " + \
                        "does not match the genome_id."
            status = "error"

        else:
            result = "The phage name is spelled correctly."
            status = "correct"

        definition = "Check phage name spelling in the organism field."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_organism_host_genus(self, eval_id=None):
        """Check host_genus name spelling in the organism field.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """

        if self.parent_host_genus != self._organism_host_genus:
            result = "The host_genus name in the organism field " + \
                        "does not match the parent_host_genus."
            status = "error"

        else:
            result = "The host_genus name is spelled correctly."
            status = "correct"

        definition = "Check host_genus name spelling in the organism field."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_host_host_genus(self, eval_id=None):
        """Check host_genus name spelling in the host field.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """

        if self.parent_host_genus != self._host_host_genus:
            result = "The host_genus name in the host field " + \
                        "does not match the parent_host_genus."
            status = "error"

        else:
            result = "The host_genus name is spelled correctly."
            status = "correct"

        definition = "Check host_genus name spelling in the host field."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def check_lab_host_host_genus(self, eval_id=None):
        """Check host_genus name spelling in the lab_host field.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """

        if self.parent_host_genus != self._lab_host_host_genus:
            result = "The host_genus name in the lab_host field " + \
                        "does not match the parent_host_genus."
            status = "error"

        else:
            result = "The host_genus name is spelled correctly."
            status = "correct"

        definition = "Check host_genus name spelling in the lab_host field."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)






























###
