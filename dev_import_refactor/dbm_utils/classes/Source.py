"""Represents a collection of data about a Source feature that is commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

from functions import basic
from classes import Eval
import re




class SourceFeature:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.type = "source"
        self.left_boundary = "" # TODO implement this.
        self.right_boundary = "" # TODO implement this.
        self.organism = ""
        self.host = ""
        self.lab_host = ""


        # Common to Phamerator.
        self.parent_genome_id = ""
        self.parent_host_genus = ""


        # Computed data fields.
        self._organism_name = ""
        self._organism_host_genus = ""
        self._host_host_genus = ""
        self._lab_host_host_genus = ""



        self.evaluations = [] # List of warnings and errors about source feature



    # TODO this is probably no longer needed.
    # def set_evaluation(self, type, message1 = None, message2 = None):
    #     """Creates an EvalResult object and adds it to the list of all
    #     evaluations."""
    #
    #     if type == "warning":
    #         eval_object = Eval.construct_warning(message1, message2)
    #     elif type == "error":
    #         eval_object = Eval.construct_error(message1)
    #     else:
    #         eval_object = Eval.EvalResult()
    #     self.evaluations.append(eval_object)

    def parse_organism(self):
        """Retrieve the phage name and host_genus name from the 'organism' field."""
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

    def check_organism_name(self):
        """Check phage name spelling in the organism field."""

        if self.parent_genome_id != self._organism_name:
            result = "The phage name in the organism field " + \
                        "does not match the parent_genome_id."
            status = "error"

        else:
            result = "The phage name is spelled correctly."
            status = "correct"

        definition = "Check phage name spelling in the organism field."
        eval = Eval.Eval(id = "SRC0001", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)



    def check_organism_host_genus(self):
        """Check host_genus name spelling in the organism field."""

        if self.parent_host_genus != self._organism_host_genus:
            result = "The host_genus name in the organism field " + \
                        "does not match the parent_host_genus."
            status = "error"

        else:
            result = "The host_genus name is spelled correctly."
            status = "correct"

        definition = "Check host_genus name spelling in the organism field."
        eval = Eval.Eval(id = "SRC0002", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)

    def check_host_host_genus(self):
        """Check host_genus name spelling in the host field."""

        if self.parent_host_genus != self._host_host_genus:
            result = "The host_genus name in the host field " + \
                        "does not match the parent_host_genus."
            status = "error"

        else:
            result = "The host_genus name is spelled correctly."
            status = "correct"

        definition = "Check host_genus name spelling in the host field."
        eval = Eval.Eval(id = "SRC0003", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)

    def check_lab_host_host_genus(self):
        """Check host_genus name spelling in the lab_host field."""

        if self.parent_host_genus != self._lab_host_host_genus:
            result = "The host_genus name in the lab_host field " + \
                        "does not match the parent_host_genus."
            status = "error"

        else:
            result = "The host_genus name is spelled correctly."
            status = "correct"

        definition = "Check host_genus name spelling in the lab_host field."
        eval = Eval.Eval(id = "SRC0004", \
                        definition = definition, \
                        result = result, \
                        status = status)
        self.evaluations.append(eval)






























###
