"""Represents a structure to directly compare data between two or more genomes."""

from classes import Eval

class DataGroup:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.ticket = ""
        self.genome_dict = {}
        self.genome_pair_dict = {}
        self.email = None # Email object for automating responses


        # Initialize calculated attributes
        self.evaluations = []
        self.sql_queries = [] # All SQL data needed to implement ticket.







    def set_genome_pair(self, genome_pair, key1, key2):
        """Pair two genomes and add to the paired genome dictionary."""

        try:
            genome_pair.genome1 = self.genome_dict[key1]
            genome_pair.genome2 = self.genome_dict[key2]
            paired_key = key1 + "_" + key2
            self.genome_pair_dict[paired_key] = genome_pair
        except:
            pass


    def check_matched_genome(self, key):
        """Check for whether a certain type of genome has been added to the
        genome dictionary."""
        if key in self.genome_dict.keys():
            result = "The %s genome has been matched." % key
            status = "correct"
        else:
            result = "No %s genome has been matched." % key
            status = "error"

        definition = "Check if a %s genome type has been " % key + \
                        "matched to the ticket."
        eval = Eval.Eval("DATAGROUP", definition, result, status)
        self.evaluations.append(eval)





###
