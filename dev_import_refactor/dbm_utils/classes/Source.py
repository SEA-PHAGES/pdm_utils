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
        self.phage_id = ""


        # Computed data fields.
        self._organism_phage_name = ""
        self._organism_host_name = ""
        self._host_host_name = ""
        self._lab_host_host_name = ""














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












    # TODO implement.
    # TODO unit test.
    def check_phage_name_typos(self, phage_name):
        """Check phage name spelling in various fields."""
        pass

    # TODO implement.
    # TODO unit test.
    def check_host_name_typos(self, host_name):
        """Check host name spelling in various fields."""
        pass










###
