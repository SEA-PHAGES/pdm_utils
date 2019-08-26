"""Represents a structure to directly compare data between two or more genomes."""

from classes import Eval
from functions import phamerator
from classes import Ticket

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





    def set_genome_pair(self, genome_pair, key1, key2):
        """Pair two genomes and add to the paired genome dictionary."""

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
            genome = self.genome_dict["remove"]
            statement = phamerator.create_genome_delete_statement(genome)
            self.sql_queries.append(statement)

        if self.ticket.type == "replace" or self.ticket.type == "add":
            genome = self.genome_dict["add"]
            statement = phamerator.create_genome_insert_statement(genome)
            self.sql_queries.append(statement)

            for cds in genome.cds_features:
                statement = create_cds_insert_statement(cds)
                self.sql_queries.append(statement)









    # Evaluations.

    def check_matched_genome(self, key, eval_id=None):
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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)

    def check_genome_dictionary(self, key, expect=True, eval_id=None):
        """Check if a genome is present in the genome dictionary.
        The 'key' parameter indicates how the genome is expected to be
        stored in the dictionary.
        The 'expect' parameter indicates whether the genome is expected
        to be present or not."""

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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)

    def check_genome_pair_dictionary(self, key, expect=True, eval_id=None):
        """Check if a genome_pair is present in the genome_pair dictionary.
        The 'key' parameter indicates how the genome_pair is expected to be
        stored in the dictionary.
        The 'expect' parameter indicates whether the genome_pair is expected
        to be present or not."""

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
        eval = Eval.Eval(eval_id, definition, result, status)
        self.evaluations.append(eval)

    def check_for_errors(self):
        """Check evaluation lists of all objects contained in the Bundle
        and determine how many errors there are."""

        for eval in self.evaluations:
            if eval.status == "error":
                self._errors += 1

        for eval in self.ticket.evaluations:
            if eval.status == "error":
                self._errors += 1

        for key in self.genome_dict.keys():
            genome = self.genome_dict[key]
            for eval in genome.evaluations:
                if eval.status == "error":
                    self._errors += 1

            for cds in genome.cds_features:
                for eval in cds.evaluations:
                    if eval.status == "error":
                        self._errors += 1

            # TODO need to implement this once this class is implemented.
            # for trna in genome.trna_features:
            #     for eval in trna.evaluations:
            #         if eval.status == "error":
            #             self._errors += 1

        for key in self.genome_pair_dict.keys():
            genome_pair = self.genome_pair_dict[key]
            for eval in genome_pair.evaluations:
                if eval.status == "error":
                    self._errors += 1


            # TODO need to implement this once this class is implemented.
            # for cds_pair in genome_pair.matched_cds_list:
            #     for eval in cds_pair.evaluations:
            #         if eval.status == "error":
            #             self._errors += 1


    def get_evaluations(self):
        """Iterate through the various objects stored in the Bundle
        object and return a dictionary of evaluation lists.
        """
        eval_dict = {}
        if len(self.evaluations) > 0:
            eval_dict["bundle"] = self.evaluations
        if isinstance(self.ticket, Ticket.GenomeTicket):
            if len(self.ticket.evaluations) > 0:
                eval_dict["ticket"] = self.ticket.evaluations
        for key in self.genome_dict.keys():
            genome = self.genome_dict[key]
            genome_key = "genome_" + key
            if len(genome.evaluations) > 0:
                eval_dict[genome_key] = genome.evaluations
            for cds in genome.cds_features:
                cds_key = "cds_" + cds.id
                if len(cds.evaluations) > 0:
                    eval_dict[cds_key] = cds.evaluations
        for key in self.genome_pair_dict.keys():
            genome_pair = self.genome_pair_dict[key]
            genome_pair_key = "genome_pair_" + key
            if len(genome_pair.evaluations) > 0:
                eval_dict[genome_pair_key] = genome_pair.evaluations
        return eval_dict




###
