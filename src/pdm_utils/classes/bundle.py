"""Represents a structure to directly compare data between two or more genomes."""

from pdm_utils.classes import eval
from pdm_utils.classes import ticket
from pdm_utils.functions import basic

class Bundle:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.id = ""
        self.ticket = None
        self.genome_dict = {}
        self.genome_pair_dict = {}
        self.email = None # Email object for automating responses


        # Initialize calculated attributes
        self.evaluations = []
        self.sql_statements = [] # All SQL data needed to implement ticket.

        self._errors = 0





    # TODO should this be changed to generate the GenomePair object directly?
    def set_genome_pair(self, genome_pair, key1, key2):
        """Pair two genomes and add to the paired genome dictionary.

        :param genome_pair:
            An empty GenomePair object to stored paried genomes.
        :type genome_pair: GenomePair
        :param key1:
            A valid key in the Bundle object's 'genome_dict'
            that indicates the first genome to be paired.
        :type key1: str
        :param key2:
            A valid key in the Bundle object's 'genome_dict'
            that indicates the second genome to be paired.
        :type key2: str
        """

        try:
            genome_pair.genome1 = self.genome_dict[key1]
            genome_pair.genome2 = self.genome_dict[key2]
            paired_key = key1 + "_" + key2
            self.genome_pair_dict[paired_key] = genome_pair
        except:
            pass



    # Evaluations.
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

    def check_ticket(self, eval_id=None, success="correct", fail="error",
                     eval_def=None):
        """Check for whether a Ticket object is present.

        :param eval_id: same as for check_genome_dict().
        :param success: same as for check_genome_dict().
        :param fail: same as for check_genome_dict().
        :param eval_def: same as for check_genome_dict().
        """
        if self.ticket is not None:
            result = ("A ticket is present. "
                      f"ID: {self.ticket.id}. "
                      f"Type: {self.ticket.type}. "
                      f"PhageID: {self.ticket.phage_id}.")
            status = success
        else:
            result = "A ticket is not present."
            status = fail
        definition = "Check if a ticket is present."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_genome_dict(self, key, expect=True, eval_id=None,
                          success="correct", fail="error", eval_def=None):
        """Check if a genome is present in the genome dictionary.

        :param key:
            The value to be evaluated if it is a valid key
            in the genome dictionary.
        :type key: str
        :param expect:
            Indicates whether the key is expected
            to be a valid key in the genome dictionary.
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
        if key in self.genome_dict.keys():
            result = f"The '{key}' genome is present."
            if expect:
                status = success
            else:
                status = fail
        else:
            result = f"The '{key}' genome is not present."
            if not expect:
                status = success
            else:
                status = fail

        definition = "Check if a genome is present."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_genome_pair_dict(self, key, expect=True, eval_id=None,
                               success="correct", fail="error", eval_def=None):
        """Check if a genome_pair is present in the genome_pair dictionary.

        :param key: same as for check_genome_dict().
        :param expect: same as for check_genome_dict().
        :param eval_id: same as for check_genome_dict().
        :param success: same as for check_genome_dict().
        :param fail: same as for check_genome_dict().
        :param eval_def: same as for check_genome_dict().
        """

        if key in self.genome_pair_dict.keys():
            result = f"The '{key}' genome_pair is present."
            if expect:
                status = success
            else:
                status = fail
        else:
            result = f"The '{key}' genome_pair is not present."
            if not expect:
                status = success
            else:
                status = fail

        definition = "Check if a genome_pair is present."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, result, status)

    def check_statements(self, execute_result, execute_msg, eval_id=None,
                         success="correct", fail="error", eval_def=None):
        """Check if MySQL statements were successfully executed.

        :param execute_result:
            Indication if MySQL statements were successfully execute.
        :type execute_result: int
        :param execute_msg:
            Description of MySQL statement execution result.
        :type execute_msg: str
        :param eval_id: same as for check_genome_dict().
        :param success: same as for check_genome_dict().
        :param fail: same as for check_genome_dict().
        :param eval_def: same as for check_genome_dict().
        """
        if execute_result == 0:
            status = success
        else:
            status = fail
        definition = "Check if MySQL statements were successfully executed."
        definition = basic.join_strings([definition, eval_def])
        self.set_eval(eval_id, definition, execute_msg, status)

    def check_for_errors(self):
        """Check evaluation lists of all objects contained in the Bundle
        and determine how many errors there are."""
        self._errors = 0
        for evl in self.evaluations:
            if evl.status == "error":
                self._errors += 1

        if self.ticket is not None:
            for evl in self.ticket.evaluations:
                if evl.status == "error":
                    self._errors += 1

        for key in self.genome_dict.keys():
            gnm = self.genome_dict[key]
            for evl in gnm.evaluations:
                if evl.status == "error":
                    self._errors += 1

            for src_ftr in gnm.source_features:
                for evl in src_ftr.evaluations:
                    if evl.status == "error":
                        self._errors += 1

            for cds_ftr in gnm.cds_features:
                for evl in cds_ftr.evaluations:
                    if evl.status == "error":
                        self._errors += 1

            # TODO need to implement this once this class is implemented.
            # for trna_ftr in gnm.trna_features:
            #     for evl in trna_ftr.evaluations:
            #         if evl.status == "error":
            #             self._errors += 1

        for key in self.genome_pair_dict.keys():
            genome_pair = self.genome_pair_dict[key]
            for evl in genome_pair.evaluations:
                if evl.status == "error":
                    self._errors += 1


            # TODO need to implement this once this class is implemented.
            # for cds_pair in genome_pair.matched_cds_list:
            #     for evl in cds_pair.evaluations:
            #         if evl.status == "error":
            #             self._errors += 1


    def get_evaluations(self):
        """Get all evaluations for all objects stored in the Bundle.

        :returns: Dictionary of evaluation lists for each feature.
        :rtype: dict
        """
        eval_dict = {}
        if len(self.evaluations) > 0:
            eval_dict["bundle"] = self.evaluations
        if isinstance(self.ticket, ticket.ImportTicket):
            if len(self.ticket.evaluations) > 0:
                eval_dict["ticket"] = self.ticket.evaluations
        for key in self.genome_dict.keys():
            gnm = self.genome_dict[key]
            genome_key = "genome_" + key
            if len(gnm.evaluations) > 0:
                eval_dict[genome_key] = gnm.evaluations
            for src_ftr in gnm.source_features:
                src_key = "src_" + src_ftr.id
                if len(src_ftr.evaluations) > 0:
                    eval_dict[src_key] = src_ftr.evaluations
            for cds_ftr in gnm.cds_features:
                cds_key = "cds_" + cds_ftr.id
                if len(cds_ftr.evaluations) > 0:
                    eval_dict[cds_key] = cds_ftr.evaluations

            # TODO get evaluations for trna and tmrna features


        for key in self.genome_pair_dict.keys():
            genome_pair = self.genome_pair_dict[key]
            genome_pair_key = "genome_pair_" + key
            if len(genome_pair.evaluations) > 0:
                eval_dict[genome_pair_key] = genome_pair.evaluations
        return eval_dict
