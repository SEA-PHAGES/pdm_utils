"""Represents a structure to pair two Genome objects and
perform comparisons between them to identify inconsistencies."""

from pdm_utils.classes import evaluation
from pdm_utils.functions import basic

# TODO methods should be added to match perfect/imperfect features between the
# two genomes as in GenomeTriad.

class GenomePair:

    # Initialize all attributes:
    def __init__(self):

        self.genome1 = None
        self.genome2 = None
        self.evaluations = []


    # Evaluations
    def set_eval(self, eval_id, definition, result, status):
        """Constructs and adds an Evaluation object to the evaluations list.

        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        :param definition: Description of the evaluation.
        :type definition: str
        :param result: Description of the outcome of the evaluation.
        :type result: str
        :param status: Outcome of the evaluation.
        :type status: str
        """
        evl = evaluation.Evaluation(eval_id, definition, result, status)
        self.evaluations.append(evl)


    def compare_attribute(self, attribute, expect_same=False, eval_id=None,
                          success="correct", fail="error", eval_def=None):
        """Compare values of the specified attribute in each genome.

        :param attribute: Name of the GenomePair object attribute to evaluate.
        :type attribute: str
        :param expect_same:
            Indicates whether the two attribute values are expected to be
            the same.
        :type expect_same: bool
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
            value1 = getattr(self.genome1, attribute)
            value2 = getattr(self.genome2, attribute)
        except:
            test = False
            value1 = None
            value2 = None
        if test:

            if value1 == value2:
                actual_same = True
            else:
                actual_same = False

            v1_short = basic.truncate_value(str(value1), 30, "...")
            v2_short = basic.truncate_value(str(value2), 30, "...")
            result = (f"The first genome is ID: {self.genome1.id}, "
                      f"Type: {self.genome1.type}. The '{attribute}' attribute "
                      f" contains: '{v1_short}'. "
                      f"The second genome is ID: {self.genome2.id}, "
                      f"Type: {self.genome2.type}. The '{attribute}' attribute "
                      f" contains: '{v2_short}'. These two values are ")

            if actual_same:
                result = result + "identical, "
            else:
                result = result + "different, "

            if actual_same and expect_same:
                result = result + "as expected."
                status = success
            elif not actual_same and not expect_same:
                result = result + "as expected."
                status = success
            else:
                result = result + "which is not expected."
                status = fail
        else:
            result = f"'{attribute}' is not a valid field to be compared."
            status = "untested"
        definition = f"Compare values of the '{attribute}' attribute in each genome."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)


    def compare_date(self, expect, eval_id=None, success="correct",
                     fail="error", eval_def=None):
        """Compare the date of each genome.

        :param expect:
            Is the first genome expected to be "newer", "equal", or "older"
            than the second genome.
        :type expect: str
        :param eval_id: same as for compare_attribute().
        :param success: same as for compare_attribute().
        :param fail: same as for compare_attribute().
        :param eval_def: same as for compare_attribute().
        """

        if expect in set(["newer", "equal", "older"]):
            if self.genome1.date > self.genome2.date:
                actual = "newer"
                actual2 = actual + "than"
            elif self.genome1.date == self.genome2.date:
                actual = "equal"
                actual2 = actual + "to"
            else:
                actual = "older"
                actual2 = actual + "than"

            msg = (f"The query genome '{self.genome1.id}' date "
                   f"is '{self.genome1.date}'."
                   f"The reference genome '{self.genome2.id}' date "
                   f"is '{self.genome2.date}'."
                   f"The date of query genome is {actual2} the "
                   "date of the reference genome, which is ")

            if actual == expect:
                result = msg + "expected."
                status = success
            else:
                result = msg + "not expected."
                status = fail
        else:
            result = f"'{expect}' is an invalid comparison."
            status = "untested"
        definition = "Compare the date of both genomes."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)
