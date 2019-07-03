"""Represents a structure to directly compare data between two or more genomes."""

from classes import Eval
from classes import GenomePair

class DataGroup:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.ticket = ""

        self.genomes_dict = {}
        self.genome_pairs_dict = {}

        self.email = None # Email object for automating responses




        # Initialize calculated attributes
        self.evaluations_dict = {} # TODO may not be needed.
        self.sql_queries = [] # All SQL data needed to implement ticket.





    #
    # def set_evaluation(self, eval_list, result_type, message1 = None, message2 = None):
    #
    #     if result_type == "warning":
    #         eval_object = Eval.construct_warning(message1, message2)
    #
    #     elif result_type == "error":
    #         eval_object = Eval.construct_error(message1)
    #
    #     else:
    #         eval_object = Eval.EvalResult()
    #
    #     self.evaluations_dict[eval_list].append(eval_object)




    # TODO not sure if I need this. It requires importing the GenomePair
    # class, which may not be necessary, since this can be performed in
    # the main script.
    # TODO unit test.
    def add_to_genome_pair_dictionary(self, key1, key2):
        """Pair two genomes and add to the paired genome dictionary."""

        genome_pair = GenomePair.GenomePair()
        genome_pair.genome1 = self.genomes_dict[key1]
        genome_pair.genome2 = self.genomes_dict[key2]
        paired_key = key1 + "_" + key2
        self.genome_pairs_dict[paired_key] = genome_pair



###
