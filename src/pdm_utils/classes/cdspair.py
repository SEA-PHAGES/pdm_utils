"""Represents a structure to directly compare data between two or more CDS
features.
"""

# TODO this class needs to be refactored, with attributes and methods
# simplified. The class is used in the compare pipeline, which has only
# been partially refactored since integrating into pdm_utils.
# It relies on extra Cds object attributes that are not present
# in the base Cds class and that are added during compare pipeline.

# Variables are prefixed to indicate genome type:
# GenBank =  "gbk", "g"
# MySQL = "mysql", "m"
# PhagesDB = "pdb", "p"


# TODO refactor and test.
class CdsPair:

    # Initialize all attributes:
    def __init__(self):

        self.type = ""

        # Initialize all non-calculated attributes:
        self._m_feature = ""
        self._g_feature = ""

        # Matched data comparison results
        self._m_g_different_translations = False
        self._m_g_different_start_sites = False
        self._m_g_different_descriptions = False

        # Total errors summary
        self._total_errors = 0


    # Define all attribute setters:

    # TODO refactor and test.
    def compare_mysql_gbk_cds_ftrs(self):

        if self._m_feature.orientation == "forward":
            if str(self._m_feature.start) != str(self._g_feature.start):
                self._m_g_different_start_sites = True
        elif self._m_feature.orientation == "reverse":
            if str(self._m_feature.stop) != str(self._g_feature.stop):
                self._m_g_different_start_sites = True
        else:
            pass


        product_description_set = set()
        product_description_set.add(self._m_feature.description)
        product_description_set.add(self._g_feature.product)


        if len(product_description_set) != 1:
            self._m_g_different_descriptions = True

        if self._m_feature.translation != self._g_feature.translation:
            self._m_g_different_translations = True

        # Compute total errors
        # First add all matched feature errors
        if self._m_g_different_translations:
            self._total_errors += 1

        if self._m_g_different_start_sites:
            self._total_errors += 1

        if self._m_g_different_descriptions:
            self._total_errors += 1

        # Now add all errors from each individual feature
        # You first compute errors for each individual feature.
        # This step is performed here instead of in the mainline code
        # because you need to wait for the feature matching step
        # after the genome matching step.
        self._m_feature.compute_total_cds_errors()
        self._g_feature.compute_total_cds_errors()
        self._total_errors += self._m_feature._total_errors
        self._total_errors += self._g_feature._total_errors
