"""Represents a structure to directly compare data between two or more genomes."""

import Eval

class DataGroup:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.ticket = ""

        # self.genome1 = None
        # self.genome2 = None
        # self.genome3 = None
        # self.genome4 = None

        self.genomes_dict = {}
        self.genome_pairs_dict = {}




        self.email = ""




        # Initialize calculated attributes
        self.evaluations_dict = {}
        self.sql_queries = []






    def set_evaluation(self, eval_list, result_type, message1 = None, message2 = None):

        if result_type == "warning":
            eval_object = Eval.construct_warning(message1, message2)

        elif result_type == "error":
            eval_object = Eval.construct_error(message1)

        else:
            eval_object = Eval.EvalResult()

        self.evaluations_dict[eval_list].append(eval_object)
