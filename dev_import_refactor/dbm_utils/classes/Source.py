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
        self.type_id = "source"
        self.left_boundary = "" # TODO implement this.
        self.right_boundary = "" # TODO implement this.
        self.organism = ""
        self.host = ""
        self.lab_host = ""


        # Common to Phamerator.
        self.parent_phage_id = ""
        self.parent_host = ""


        # Computed data fields.
        self._organism_phage_name = ""
        self._organism_host_name = ""
        self._host_host_name = ""
        self._lab_host_host_name = ""



        self.evaluations = [] # List of warnings and errors about source feature




    def set_evaluation(self, type, message1 = None, message2 = None):
        """Creates an EvalResult object and adds it to the list of all
        evaluations."""

        if type == "warning":
            eval_object = Eval.construct_warning(message1, message2)
        elif type == "error":
            eval_object = Eval.construct_error(message1)
        else:
            eval_object = Eval.EvalResult()
        self.evaluations.append(eval_object)

    def parse_organism(self):
        """Retrieve the phage name and host name from the 'organism' field."""
        self._organism_phage_name, \
        self._organism_host_name = \
            basic.parse_names_from_record_field(self.organism)

    def parse_host(self):
        """Retrieve the host name from the 'host' field."""
        phage_name, \
        self._host_host_name = \
            basic.parse_names_from_record_field(self.host)
        # Note: no need to assign phage name, since this field is only
        # expected to contain host information.

    def parse_lab_host(self):
        """Retrieve the host name from the 'lab_host' field."""
        phage_name, \
        self._lab_host_host_name = \
            basic.parse_names_from_record_field(self.lab_host)
        # Note: no need to assign phage name, since this field is only
        # expected to contain host information.


    def check_organism_phage_name(self):
        """Check phage name spelling in the organism field."""

        if self.parent_phage_id != self._organism_phage_name:
            message1 = "The phage name in the organism field " + \
                        "does not match the parent_phage_id."
            self.set_evaluation("warning", message1, message1)

    def check_organism_host_name(self):
        """Check host name spelling in the organism field."""

        if self.parent_host != self._organism_host_name:
            message1 = "The host name in the organism field " + \
                        "does not match the parent_host."
            self.set_evaluation("warning", message1, message1)

    def check_host_host_name(self):
        """Check host name spelling in the host field."""

        if self.parent_host != self._host_host_name:
            message1 = "The host name in the host field " + \
                        "does not match the parent_host."
            self.set_evaluation("warning", message1, message1)

    def check_lab_host_host_name(self):
        """Check host name spelling in the lab_host field."""

        if self.parent_host != self._lab_host_host_name:
            message1 = "The host name in the lab_host field " + \
                        "does not match the parent_host."
            self.set_evaluation("warning", message1, message1)































###
