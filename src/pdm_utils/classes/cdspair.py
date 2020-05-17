"""Represents a structure to directly compare data between two or more CDS
features.
"""

# TODO this class needs to be refactored, with attributes and methods
# simplified. The class is used in the compare pipeline, which has only
# been partially refactored since integrating into pdm_utils.
# It relies on extra Cds object attributes that are NOT present
# in the base Cds class and that are added during compare pipeline.
# So do NOT use this class for anything other than in the compare pipeline
# until it has been properly refactored.

# TODO refactor and test.
class CdsPair:

    # Initialize all attributes:
    def __init__(self):

        self.type = ""
        self.cds1 = None
        self.cds2 = None

        self.different_translation = False
        self.different_start_site = False
        self.different_description = False
        self._total_errors = 0


    # Define all attribute setters:

    # TODO refactor and test.
    def compare_cds(self):

        if self.cds1.orientation == "forward":
            if str(self.cds1.start) != str(self.cds2.start):
                self.different_start_site = True
        elif self.cds1.orientation == "reverse":
            if str(self.cds1.stop) != str(self.cds2.stop):
                self.different_start_site = True
        else:
            pass


        product_description_set = set()
        product_description_set.add(self.cds1.description)
        product_description_set.add(self.cds2.product)


        if len(product_description_set) != 1:
            self.different_description = True

        if self.cds1.translation != self.cds2.translation:
            self.different_translation = True

        # Compute total errors
        # First add all matched feature errors
        if self.different_translation:
            self._total_errors += 1

        if self.different_start_site:
            self._total_errors += 1

        if self.different_description:
            self._total_errors += 1

        # Now add all errors from each individual feature
        # You first compute errors for each individual feature.
        # This step is performed here instead of in the mainline code
        # because you need to wait for the feature matching step
        # after the genome matching step.
        self.cds1.check_for_errors()
        self.cds2.check_for_errors()
        self._total_errors += self.cds1._total_errors
        self._total_errors += self.cds2._total_errors
