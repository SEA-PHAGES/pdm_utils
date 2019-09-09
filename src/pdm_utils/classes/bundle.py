"""Represents a structure to directly compare data between two or more genomes."""

from pdm_utils.classes import eval
from pdm_utils.functions import phamerator
from pdm_utils.classes import ticket

class Bundle:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.id = ""
        self.ticket = ""
        self.genome_dict = {}
        self.genome_pair_dict = {}
        self.email = None # Email object for automating responses


        # Initialize calculated attributes
        self.evaluations = []
        self.sql_queries = [] # All SQL data needed to implement ticket.

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





    # TODO implement.
    # TODO unit test.
    def create_sql_statements(self):
        """Create list of MySQL statements based on the ticket type."""

        if self.ticket.type == "replace" or self.ticket.type == "remove":
            gnm = self.genome_dict["remove"]
            statement = phamerator.create_genome_delete_statement(gnm)
            self.sql_queries.append(statement)

        if self.ticket.type == "replace" or self.ticket.type == "add":
            gnm = self.genome_dict["add"]
            statement = phamerator.create_genome_insert_statement(gnm)
            self.sql_queries.append(statement)

            for cds_ftr in gnm.cds_features:
                statement = create_cds_insert_statement(cds_ftr)
                self.sql_queries.append(statement)









    # Evaluations.


    # TODO is this needed? Seems to have been replaced by
    # check_genome_dictionary.
    def check_matched_genome(self, key, eval_id=None):
        """Check for whether a certain type of genome has been added to the
        genome dictionary.

        :param key:
            The value to be evaluated if it is a valid key
            in the Bundle object's 'genome_dict'.
        :type key: str
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """
        if key in self.genome_dict.keys():
            result = "The %s genome has been matched." % key
            status = "correct"
        else:
            result = "No %s genome has been matched." % key
            status = "error"

        definition = "Check if a %s genome type has been " % key + \
                        "matched to the ticket."
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)

    def check_genome_dictionary(self, key, expect=True, eval_id=None):
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
        """

        if key in self.genome_dict.keys():
            result = "The %s genome is present." % key
            if expect:
                status = "correct"
            else:
                status = "error"
        else:
            result = "The %s genome is not present." % key
            if not expect:
                status = "correct"
            else:
                status = "error"

        definition = "Check if the %s genome is present." % key
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)

    def check_genome_pair_dictionary(self, key, expect=True, eval_id=None):
        """Check if a genome_pair is present in the genome_pair dictionary.

        :param key:
            The value to be evaluated if it is a valid key
            in the genome_pair dictionary.
        :type key: str
        :param expect:
            Indicates whether the key is expected
            to be a valid key in the genome_pair dictionary.
        :type expect: bool
        :param eval_id: Unique identifier for the evaluation.
        :type eval_id: str
        """

        if key in self.genome_pair_dict.keys():
            result = "The %s genome_pair is present." % key
            if expect:
                status = "correct"
            else:
                status = "error"
        else:
            result = "The %s genome_pair is not present." % key
            if not expect:
                status = "correct"
            else:
                status = "error"

        definition = "Check if the %s genome_pair is present." % key
        evl = eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(evl)

    def check_for_errors(self):
        """Check evaluation lists of all objects contained in the Bundle
        and determine how many errors there are."""

        for evl in self.evaluations:
            if evl.status == "error":
                self._errors += 1

        for evl in self.ticket.evaluations:
            if evl.status == "error":
                self._errors += 1

        for key in self.genome_dict.keys():
            gnm = self.genome_dict[key]
            for evl in gnm.evaluations:
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
        """Iterate through the various objects stored in the Bundle
        object and return a dictionary of evaluation lists.
        """
        eval_dict = {}
        if len(self.evaluations) > 0:
            eval_dict["bundle"] = self.evaluations
        if isinstance(self.ticket, ticket.GenomeTicket):
            if len(self.ticket.evaluations) > 0:
                eval_dict["ticket"] = self.ticket.evaluations
        for key in self.genome_dict.keys():
            gnm = self.genome_dict[key]
            genome_key = "genome_" + key
            if len(gnm.evaluations) > 0:
                eval_dict[genome_key] = gnm.evaluations
            for cds_ftr in gnm.cds_features:
                cds_key = "cds_" + cds_ftr.id
                if len(cds_ftr.evaluations) > 0:
                    eval_dict[cds_key] = cds_ftr.evaluations
        for key in self.genome_pair_dict.keys():
            genome_pair = self.genome_pair_dict[key]
            genome_pair_key = "genome_pair_" + key
            if len(genome_pair.evaluations) > 0:
                eval_dict[genome_pair_key] = genome_pair.evaluations
        return eval_dict




###
