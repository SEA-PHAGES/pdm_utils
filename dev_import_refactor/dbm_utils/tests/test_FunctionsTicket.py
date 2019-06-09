""" Unit tests for misc. ticket functions"""

from classes import DataGroup
from classes import Genome
from classes import Ticket
from classes import Eval
from functions import FunctionsTicket
import unittest





class TestTicketFunctions1(unittest.TestCase):


    def setUp(self):
        self.normal_ticket_list = ["add",
                                "Trixe",
                                "Mycobacterium",
                                "A",
                                "A2",
                                "Final",
                                "Hatfull",
                                "Product",
                                "ABC123",
                                "phagesdb",
                                "none"]

        self.short_ticket_list = ["add",
                                "Trixe",
                                "Mycobacterium",
                                "A",
                                "A2",
                                "Final",
                                "Hatfull"]

        self.long_ticket_list = ["add",
                                "Trixe",
                                "Mycobacterium",
                                "A",
                                "A2",
                                "Final",
                                "Hatfull",
                                "Product",
                                "ABC123",
                                "phagesdb",
                                "none",
                                "extra"]


        self.remove_ticket_list = ["remove",
                                    "none",
                                    "none",
                                    "none",
                                    "none",
                                    "none",
                                    "none",
                                    "none",
                                    "none",
                                    "none",
                                    "Trixe"]


    def test_parse_import_ticket_1(self):
        """Verify properly structured data is parsed correctly."""
        ticket, eval = \
            FunctionsTicket.parse_import_ticket(self.normal_ticket_list)
        with self.subTest():
            self.assertIsNone(eval)
        with self.subTest():
            self.assertIsNotNone(ticket)
        with self.subTest():
            self.assertEqual(ticket.type, "add")
        with self.subTest():
            self.assertEqual(ticket.primary_phage_id, "Trixe")
        with self.subTest():
            self.assertEqual(ticket.host, "Mycobacterium")
        with self.subTest():
            self.assertEqual(ticket.cluster, "A")
        with self.subTest():
            self.assertEqual(ticket.subcluster, "A2")
        with self.subTest():
            self.assertEqual(ticket.status, "Final")
        with self.subTest():
            self.assertEqual(ticket.annotation_author, "Hatfull")
        with self.subTest():
            self.assertEqual(ticket.description_field, "Product")
        with self.subTest():
            self.assertEqual(ticket.accession, "ABC123")
        with self.subTest():
            self.assertEqual(ticket.run_mode, "phagesdb")
        with self.subTest():
            self.assertEqual(ticket.secondary_phage_id, "none")

    def test_parse_import_ticket_2(self):
        """Verify improperly structured data is not parsed."""
        ticket, eval = \
            FunctionsTicket.parse_import_ticket(self.short_ticket_list)
        with self.subTest():
            self.assertIsNotNone(eval)
        with self.subTest():
            self.assertIsNone(ticket)
        with self.subTest():
            self.assertEqual(eval.status, "error")

    def test_parse_import_ticket_3(self):
        """Verify improperly structured data is not parsed."""
        ticket, eval = \
            FunctionsTicket.parse_import_ticket(self.long_ticket_list)
        with self.subTest():
            self.assertIsNotNone(eval)
        with self.subTest():
            self.assertIsNone(ticket)
        with self.subTest():
            self.assertEqual(eval.status, "error")












    def test_parse_import_tickets_1(self):
        """Verify two lists of correct data are parsed."""
        list_of_lists = [self.normal_ticket_list,
                         self.remove_ticket_list]
        list_of_tickets, list_of_evals = \
            FunctionsTicket.parse_import_tickets(list_of_lists)
        with self.subTest():
            self.assertEqual(len(list_of_evals), 0)
        with self.subTest():
            self.assertEqual(len(list_of_tickets), 2)
        with self.subTest():
            type = list_of_tickets[0].type
            self.assertEqual(type, "add")
        with self.subTest():
            type = list_of_tickets[1].type
            self.assertEqual(type, "remove")


    def test_parse_import_tickets_2(self):
        """Verify two lists of incorrect data are not parsed."""
        list_of_lists = [self.short_ticket_list,
                         self.long_ticket_list]
        list_of_tickets, list_of_evals = \
            FunctionsTicket.parse_import_tickets(list_of_lists)
        with self.subTest():
            self.assertEqual(len(list_of_evals), 2)
        with self.subTest():
            self.assertEqual(len(list_of_tickets), 0)

    def test_parse_import_tickets_3(self):
        """Verify mixed lists of correct and incorrect data are
        parsed correctly."""
        list_of_lists = [self.short_ticket_list,
                        self.normal_ticket_list,
                         self.long_ticket_list]
        list_of_tickets, list_of_evals = \
            FunctionsTicket.parse_import_tickets(list_of_lists)
        with self.subTest():
            self.assertEqual(len(list_of_evals), 2)
        with self.subTest():
            self.assertEqual(len(list_of_tickets), 1)







    def test_identify_duplicates1_1(self):
        """Verify duplicate items generates an error."""
        eval_object = \
            FunctionsTicket.identify_duplicates1(["Trixie", "Trixie"], "temp")
        with self.subTest():
            self.assertIsNotNone(eval_object)
        with self.subTest():
            self.assertIsInstance(eval_object,Eval.EvalResult)

    def test_identify_duplicates1_2(self):
        """Verify non-duplicate items do not generate an error."""
        eval_object = \
            FunctionsTicket.identify_duplicates1(["Trixie", "L5"], "temp")
        with self.subTest():
            self.assertIsNone(eval_object)
        with self.subTest():
            self.assertNotIsInstance(eval_object,Eval.EvalResult)











    def test_identify_duplicates2_1(self):
        """Verify duplicate items in two lists generate an error."""
        list1 = ["Trixie", "L5"]
        list2 = ["Trixie", "RedRock"]
        eval_object = \
            FunctionsTicket.identify_duplicates2(list1, list2, "temp")
        with self.subTest():
            self.assertIsNotNone(eval_object)
        with self.subTest():
            self.assertIsInstance(eval_object,Eval.EvalResult)




    def test_identify_duplicates2_2(self):
        """Verify non-duplicate items in two lists do not generate an error."""
        list1 = ["Trixie", "L5"]
        list2 = ["D29", "RedRock"]
        eval_object = \
            FunctionsTicket.identify_duplicates2(list1, list2, "temp")
        with self.subTest():
            self.assertIsNone(eval_object)
        with self.subTest():
            self.assertNotIsInstance(eval_object,Eval.EvalResult)




    def test_validate_tickets_1(self):
        """Verify no duplicates do not generate an error."""

        ticket1 = Ticket.ImportTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "Trixie"
        ticket1.secondary_phage_id = "Trixie"
        ticket1.accession = "ABC123"

        ticket2 = Ticket.ImportTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "L5"
        ticket2.secondary_phage_id = "L5"
        ticket2.accession = "EFG456"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            FunctionsTicket.validate_tickets(list_of_tickets)

        self.assertEqual(len(result_list), 0)

    def test_validate_tickets_2(self):
        """Verify 'none' duplicates do not generate an error."""

        ticket1 = Ticket.ImportTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "none"
        ticket1.secondary_phage_id = "none"
        ticket1.accession = "none"

        ticket2 = Ticket.ImportTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "none"
        ticket2.secondary_phage_id = "none"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            FunctionsTicket.validate_tickets(list_of_tickets)

        self.assertEqual(len(result_list), 0)




    def test_validate_tickets_3(self):
        """Verify Primary Phage ID duplicates do generate an error."""

        ticket1 = Ticket.ImportTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "Trixie"
        ticket1.secondary_phage_id = "none"
        ticket1.accession = "none"

        ticket2 = Ticket.ImportTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "Trixie"
        ticket2.secondary_phage_id = "none"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            FunctionsTicket.validate_tickets(list_of_tickets)

        self.assertEqual(len(result_list), 1)


    def test_validate_tickets_4(self):
        """Verify Secondary Phage ID duplicates generate an error."""

        ticket1 = Ticket.ImportTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "none"
        ticket1.secondary_phage_id = "Trixie"
        ticket1.accession = "none"

        ticket2 = Ticket.ImportTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "none"
        ticket2.secondary_phage_id = "Trixie"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            FunctionsTicket.validate_tickets(list_of_tickets)

        self.assertEqual(len(result_list), 1)


    def test_validate_tickets_5(self):
        """Verify Accession duplicates generate an error."""

        ticket1 = Ticket.ImportTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "none"
        ticket1.secondary_phage_id = "none"
        ticket1.accession = "ABC123"

        ticket2 = Ticket.ImportTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "none"
        ticket2.secondary_phage_id = "none"
        ticket2.accession = "ABC123"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            FunctionsTicket.validate_tickets(list_of_tickets)

        self.assertEqual(len(result_list), 1)


    def test_validate_tickets_6(self):
        """Verify multiple duplicates generate multiple errors."""

        ticket1 = Ticket.ImportTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "Trixie"
        ticket1.secondary_phage_id = "Trixie"
        ticket1.accession = "ABC123"

        ticket2 = Ticket.ImportTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "Trixie"
        ticket2.secondary_phage_id = "Trixie"
        ticket2.accession = "ABC123"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            FunctionsTicket.validate_tickets(list_of_tickets)

        self.assertEqual(len(result_list), 3)








    def test_validate_tickets_7(self):
        """Verify duplicate Primary Phage ID from non-standard ticket type
        does not generate an error."""

        ticket1 = Ticket.ImportTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "Trixie"
        ticket1.secondary_phage_id = "none"
        ticket1.accession = "none"

        ticket2 = Ticket.ImportTicket()
        ticket2.type = "other"
        ticket2.primary_phage_id = "Trixie"
        ticket2.secondary_phage_id = "none"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            FunctionsTicket.validate_tickets(list_of_tickets)

        self.assertEqual(len(result_list), 0)

    def test_validate_tickets_8(self):
        """Verify duplicate Secondary Phage ID from non-standard ticket type
        does not generate an error."""

        ticket1 = Ticket.ImportTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "none"
        ticket1.secondary_phage_id = "Trixie"
        ticket1.accession = "none"

        ticket2 = Ticket.ImportTicket()
        ticket2.type = "other"
        ticket2.primary_phage_id = "none"
        ticket2.secondary_phage_id = "Trixie"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            FunctionsTicket.validate_tickets(list_of_tickets)

        self.assertEqual(len(result_list), 0)

    def test_validate_tickets_9(self):
        """Verify Primary Phage ID duplicates from different ticket
        types generate an error."""

        ticket1 = Ticket.ImportTicket()
        ticket1.type = "add"
        ticket1.primary_phage_id = "Trixie"
        ticket1.secondary_phage_id = "none"
        ticket1.accession = "none"

        ticket2 = Ticket.ImportTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "Trixie"
        ticket2.secondary_phage_id = "none"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            FunctionsTicket.validate_tickets(list_of_tickets)

        self.assertEqual(len(result_list), 1)

    def test_validate_tickets_10(self):
        """Verify Secondary Phage ID duplicates from different ticket
        types generate an error."""

        ticket1 = Ticket.ImportTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "none"
        ticket1.secondary_phage_id = "Trixie"
        ticket1.accession = "none"

        ticket2 = Ticket.ImportTicket()
        ticket2.type = "remove"
        ticket2.primary_phage_id = "none"
        ticket2.secondary_phage_id = "Trixie"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            FunctionsTicket.validate_tickets(list_of_tickets)

        self.assertEqual(len(result_list), 1)

    def test_validate_tickets_11(self):
        """Verify Accession duplicates from different ticket
        types generate an error."""

        ticket1 = Ticket.ImportTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "none"
        ticket1.secondary_phage_id = "none"
        ticket1.accession = "ABC123"

        ticket2 = Ticket.ImportTicket()
        ticket2.type = "remove"
        ticket2.primary_phage_id = "none"
        ticket2.secondary_phage_id = "none"
        ticket2.accession = "ABC123"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            FunctionsTicket.validate_tickets(list_of_tickets)
        self.assertEqual(len(result_list), 1)

    def test_validate_tickets_12(self):
        """Verify a conflict between an update/add Primary Phage ID and all
        Secondary Phage IDs generate an error."""

        ticket1 = Ticket.ImportTicket()
        ticket1.type = "add"
        ticket1.primary_phage_id = "Trixie"
        ticket1.secondary_phage_id = "none"
        ticket1.accession = "none"

        ticket2 = Ticket.ImportTicket()
        ticket2.type = "remove"
        ticket2.primary_phage_id = "none"
        ticket2.secondary_phage_id = "Trixie"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            FunctionsTicket.validate_tickets(list_of_tickets)

        self.assertEqual(len(result_list), 1)

    def test_validate_tickets_13(self):
        """Verify a conflict between a replace Primary Phage ID and a remove
        Secondary Phage ID generate an error."""

        ticket1 = Ticket.ImportTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "Trixie"
        ticket1.secondary_phage_id = "none"
        ticket1.accession = "none"

        ticket2 = Ticket.ImportTicket()
        ticket2.type = "remove"
        ticket2.primary_phage_id = "none"
        ticket2.secondary_phage_id = "Trixie"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            FunctionsTicket.validate_tickets(list_of_tickets)

        self.assertEqual(len(result_list), 1)







class TestTicketFunctions2(unittest.TestCase):

    def setUp(self):

        self.genome1 = Genome.Genome()
        self.genome2 = Genome.Genome()
        self.genome3 = Genome.Genome()

        self.genome1.phage_id = "Trixie"
        self.genome2.phage_id = "L5"
        self.genome3.phage_id = "D29"

        self.genome1.filename = "trixie_flat_file"
        self.genome2.filename = "l5_flat_file"
        self.genome3.filename = "d29_flat_file"

        self.genome_dict = {self.genome1.phage_id:self.genome1,
                            self.genome2.phage_id:self.genome2,
                            self.genome3.phage_id:self.genome3}

        self.genome_list = [self.genome1, self.genome2, self.genome3]


        self.ticket1 = Ticket.ImportTicket()
        self.ticket2 = Ticket.ImportTicket()
        self.ticket3 = Ticket.ImportTicket()
        self.ticket4 = Ticket.ImportTicket()
        self.ticket5 = Ticket.ImportTicket()


        self.ticket1.primary_phage_id = "Trixie"
        self.ticket2.primary_phage_id = "L5"
        self.ticket3.primary_phage_id = "D29"
        self.ticket4.primary_phage_id = "RedRock"
        self.ticket5.primary_phage_id = "EagleEye"


        self.datagroup1 = DataGroup.DataGroup()
        self.datagroup2 = DataGroup.DataGroup()
        self.datagroup3 = DataGroup.DataGroup()
        self.datagroup4 = DataGroup.DataGroup()
        self.datagroup5 = DataGroup.DataGroup()

        self.datagroup1.ticket = self.ticket1
        self.datagroup2.ticket = self.ticket2
        self.datagroup3.ticket = self.ticket3
        self.datagroup4.ticket = self.ticket4
        self.datagroup5.ticket = self.ticket5




    def test_match_genomes_to_tickets1_1(self):
        """Verify that one genome is matched to ticket correctly."""

        list1 = [self.datagroup1] # Trixie
        eval_list = \
            FunctionsTicket.match_genomes_to_tickets1(list1,
                                                    self.genome_dict,
                                                    "phamerator")

        matched_genome = list1[0].genomes_dict["phamerator"]
        id = matched_genome.phage_id
        expected_id = "Trixie"
        with self.subTest():
            self.assertEqual(id, expected_id)
        with self.subTest():
            self.assertEqual(len(eval_list), 0)

    def test_match_genomes_to_tickets1_2(self):
        """Verify that no genome is matched to empty ticket list."""

        list1 = []
        eval_list = \
            FunctionsTicket.match_genomes_to_tickets1(list1,
                                                    self.genome_dict,
                                                    "phamerator")

        self.assertEqual(len(eval_list), 0)

    def test_match_genomes_to_tickets1_3(self):
        """Verify that genome is not matched to ticket."""

        list1 = [self.datagroup4] # RedRock

        eval_list = \
            FunctionsTicket.match_genomes_to_tickets1(list1,
                                                    self.genome_dict,
                                                    "phamerator")

        num_dict_keys = len(list1[0].genomes_dict.keys())

        with self.subTest():
            self.assertEqual(num_dict_keys, 0)
        with self.subTest():
            self.assertEqual(len(eval_list), 1)

    def test_match_genomes_to_tickets1_4(self):
        """Verify that two genomes are matched to tickets correctly."""

        list1 = [self.datagroup1, self.datagroup2] # Trixie, L5

        eval_list = \
            FunctionsTicket.match_genomes_to_tickets1(list1,
                                                    self.genome_dict,
                                                    "phamerator")

        matched_genome1 = list1[0].genomes_dict["phamerator"]
        matched_genome2 = list1[1].genomes_dict["phamerator"]

        id1 = matched_genome1.phage_id
        expected_id1 = "Trixie"

        id2 = matched_genome2.phage_id
        expected_id2 = "L5"

        with self.subTest():
            self.assertEqual(id1, expected_id1)
        with self.subTest():
            self.assertEqual(id2, expected_id2)
        with self.subTest():
            self.assertEqual(len(eval_list), 0)

    def test_match_genomes_to_tickets1_5(self):
        """Verify that first genome is not matched to ticket correctly."""

        list1 = [self.datagroup4, self.datagroup2] # RedRock, L5

        eval_list = \
            FunctionsTicket.match_genomes_to_tickets1(list1,
                                                    self.genome_dict,
                                                    "phamerator")

        num_dict_keys1 = len(list1[0].genomes_dict.keys())
        matched_genome2 = list1[1].genomes_dict["phamerator"]

        id2 = matched_genome2.phage_id
        expected_id2 = "L5"

        with self.subTest():
            self.assertEqual(num_dict_keys1, 0)
        with self.subTest():
            self.assertEqual(id2, expected_id2)
        with self.subTest():
            self.assertEqual(len(eval_list), 1)

    def test_match_genomes_to_tickets1_6(self):
        """Verify that second genome is not matched to ticket correctly."""

        list1 = [self.datagroup2, self.datagroup4] # L5, RedRock

        eval_list = \
            FunctionsTicket.match_genomes_to_tickets1(list1,
                                                    self.genome_dict,
                                                    "phamerator")

        matched_genome1 = list1[0].genomes_dict["phamerator"]
        num_dict_keys2 = len(list1[1].genomes_dict.keys())

        id1 = matched_genome1.phage_id
        expected_id1 = "L5"

        with self.subTest():
            self.assertEqual(num_dict_keys2, 0)
        with self.subTest():
            self.assertEqual(id1, expected_id1)
        with self.subTest():
            self.assertEqual(len(eval_list), 1)

    def test_match_genomes_to_tickets1_7(self):
        """Verify that two genomes are not matched to tickets."""

        list1 = [self.datagroup4, self.datagroup5] # RedRock, EagleEye

        eval_list = \
            FunctionsTicket.match_genomes_to_tickets1(list1,
                                                    self.genome_dict,
                                                    "phamerator")

        num_dict_keys1 = len(list1[0].genomes_dict.keys())
        num_dict_keys2 = len(list1[1].genomes_dict.keys())

        with self.subTest():
            self.assertEqual(num_dict_keys1, 0)
        with self.subTest():
            self.assertEqual(num_dict_keys2, 0)
        with self.subTest():
            self.assertEqual(len(eval_list), 2)




    def test_assign_match_strategy_1(self):
        """Verify strategy is assigned with no error produced."""
        input_strategy = "phage_id"
        self.datagroup1.ticket.match_strategy = input_strategy
        self.datagroup2.ticket.match_strategy = input_strategy
        list1 = [self.datagroup1, self.datagroup2] # Trixie, L5
        output_strategy, eval_result = \
            FunctionsTicket.assign_match_strategy(list1)
        with self.subTest():
            self.assertEqual(output_strategy, input_strategy)
        with self.subTest():
            self.assertIsNone(eval_result)

    def test_assign_match_strategy_2(self):
        """Verify no strategy is assigned and an error is produced."""
        input_strategy1 = "phage_id"
        input_strategy2 = "filename"
        self.datagroup1.ticket.match_strategy = input_strategy1
        self.datagroup2.ticket.match_strategy = input_strategy2
        list1 = [self.datagroup1, self.datagroup2] # Trixie, L5
        output_strategy, eval_result = \
            FunctionsTicket.assign_match_strategy(list1)
        expected_strategy = ""
        with self.subTest():
            self.assertEqual(output_strategy, expected_strategy)
        with self.subTest():
            self.assertIsNotNone(eval_result)










###

class TestTicketFunctions3(unittest.TestCase):

    def setUp(self):

        self.genome1 = Genome.Genome()
        self.genome2 = Genome.Genome()
        self.genome3 = Genome.Genome()
        self.genome4 = Genome.Genome()

        self.ticket1 = Ticket.ImportTicket()
        self.ticket2 = Ticket.ImportTicket()
        self.ticket3 = Ticket.ImportTicket()
        self.ticket4 = Ticket.ImportTicket()

        self.datagroup1 = DataGroup.DataGroup()
        self.datagroup2 = DataGroup.DataGroup()
        self.datagroup3 = DataGroup.DataGroup()
        self.datagroup4 = DataGroup.DataGroup()

        self.datagroup1.ticket = self.ticket1
        self.datagroup2.ticket = self.ticket2
        self.datagroup3.ticket = self.ticket3
        self.datagroup4.ticket = self.ticket4


    def test_match_genomes_to_tickets2_1(self):
        """Verify that one genome is matched to ticket using phage_id."""

        self.datagroup1.ticket.primary_phage_id = "Trixie"
        self.datagroup1.ticket.match_strategy = "phage_id"
        self.genome1.phage_id = "Trixie"


        list1 = [self.datagroup1] # Trixie
        list2 = [self.genome1] # Trixie
        eval_list = \
            FunctionsTicket.match_genomes_to_tickets2(list1,
                                                    list2,
                                                    "phamerator")

        matched_genome = list1[0].genomes_dict["phamerator"]
        id = matched_genome.phage_id
        expected_id = "Trixie"
        with self.subTest():
            self.assertEqual(id, expected_id)
        with self.subTest():
            self.assertEqual(len(eval_list), 0)

    def test_match_genomes_to_tickets2_2(self):
        """Verify that one genome is matched to ticket using filename."""

        self.datagroup1.ticket.primary_phage_id = "Trixie"
        self.datagroup1.ticket.match_strategy = "filename"
        self.genome1.filename = "Trixie"

        list1 = [self.datagroup1] # Trixie
        list2 = [self.genome1] # Trixie
        eval_list = \
            FunctionsTicket.match_genomes_to_tickets2(list1,
                                                    list2,
                                                    "flat_file")

        matched_genome = list1[0].genomes_dict["flat_file"]
        id = matched_genome.filename
        expected_id = "Trixie"
        with self.subTest():
            self.assertEqual(id, expected_id)
        with self.subTest():
            self.assertEqual(len(eval_list), 0)

    def test_match_genomes_to_tickets2_3(self):
        """Verify that one genome is matched to ticket,
        and one genome is not matched (no matching ticket), using filename."""

        self.datagroup1.ticket.primary_phage_id = "Trixie"
        self.datagroup1.ticket.match_strategy = "filename"
        self.genome1.filename = "Trixie"
        self.genome2.filename = "L5"


        list1 = [self.datagroup1] # Trixie
        list2 = [self.genome1, self.genome2] # Trixie, L5
        eval_list = \
            FunctionsTicket.match_genomes_to_tickets2(list1,
                                                    list2,
                                                    "flat_file")

        matched_genome = list1[0].genomes_dict["flat_file"]
        id = matched_genome.filename
        expected_id = "Trixie"
        with self.subTest():
            self.assertEqual(id, expected_id)
        with self.subTest():
            self.assertEqual(len(eval_list), 1)

    def test_match_genomes_to_tickets2_4(self):
        """Verify that one ticket is matched to genome,
        and one ticket is not matched (no matching genome), using filename."""

        self.datagroup1.ticket.primary_phage_id = "Trixie"
        self.datagroup1.ticket.match_strategy = "filename"
        self.genome1.filename = "Trixie"

        self.datagroup2.ticket.primary_phage_id = "L5"
        self.datagroup2.ticket.match_strategy = "filename"

        list1 = [self.datagroup1, self.datagroup2] # Trixie, L5
        list2 = [self.genome1] # Trixie

        eval_list = \
            FunctionsTicket.match_genomes_to_tickets2(list1,
                                                    list2,
                                                    "flat_file")

        matched_genome = list1[0].genomes_dict["flat_file"]
        id = matched_genome.filename
        expected_id = "Trixie"
        with self.subTest():
            self.assertEqual(id, expected_id)
        with self.subTest():
            self.assertEqual(len(eval_list), 1)

    def test_match_genomes_to_tickets2_5(self):
        """Verify that two genomes are matched to tickets,
        using filename."""

        self.datagroup1.ticket.primary_phage_id = "Trixie"
        self.datagroup1.ticket.match_strategy = "filename"
        self.genome1.filename = "Trixie"

        self.datagroup2.ticket.primary_phage_id = "L5"
        self.datagroup2.ticket.match_strategy = "filename"
        self.genome2.filename = "L5"


        list1 = [self.datagroup1, self.datagroup2] # Trixie, L5
        list2 = [self.genome1, self.genome2] # Trixie, L5
        eval_list = \
            FunctionsTicket.match_genomes_to_tickets2(list1,
                                                    list2,
                                                    "flat_file")

        matched_genome1 = list1[0].genomes_dict["flat_file"]
        id1 = matched_genome1.filename
        expected_id1 = "Trixie"

        matched_genome2 = list1[1].genomes_dict["flat_file"]
        id2 = matched_genome2.filename
        expected_id2 = "L5"

        with self.subTest():
            self.assertEqual(id1, expected_id1)
        with self.subTest():
            self.assertEqual(id2, expected_id2)
        with self.subTest():
            self.assertEqual(len(eval_list), 0)

    def test_match_genomes_to_tickets2_6(self):
        """Verify that no genomes are matched to tickets due to
        conflicting strategies."""

        self.datagroup1.ticket.primary_phage_id = "Trixie"
        self.datagroup1.ticket.match_strategy = "phage_id"
        self.genome1.filename = "Trixie"

        self.datagroup2.ticket.primary_phage_id = "L5"
        self.datagroup2.ticket.match_strategy = "filename"
        self.genome2.filename = "L5"

        list1 = [self.datagroup1, self.datagroup2] # Trixie, L5
        list2 = [self.genome1, self.genome2] # Trixie, L5
        eval_list = \
            FunctionsTicket.match_genomes_to_tickets2(list1,
                                                    list2,
                                                    "flat_file")

        with self.subTest():
            self.assertEqual(len(list1[0].genomes_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(len(list1[0].genomes_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(len(eval_list), 3)

    def test_match_genomes_to_tickets2_7(self):
        """Verify that three genomes are matched to tickets,
        using filename."""

        self.datagroup1.ticket.primary_phage_id = "Trixie"
        self.datagroup1.ticket.match_strategy = "filename"
        self.genome1.filename = "Trixie"

        self.datagroup2.ticket.primary_phage_id = "L5"
        self.datagroup2.ticket.match_strategy = "filename"
        self.genome2.filename = "L5"

        self.datagroup3.ticket.primary_phage_id = "D29"
        self.datagroup3.ticket.match_strategy = "filename"
        self.genome3.filename = "D29"

        list1 = [self.datagroup1, self.datagroup2, self.datagroup3]
        list2 = [self.genome1, self.genome2, self.genome3]
        eval_list = \
            FunctionsTicket.match_genomes_to_tickets2(list1,
                                                    list2,
                                                    "flat_file")

        matched_genome1 = list1[0].genomes_dict["flat_file"]
        id1 = matched_genome1.filename
        expected_id1 = "Trixie"

        matched_genome2 = list1[1].genomes_dict["flat_file"]
        id2 = matched_genome2.filename
        expected_id2 = "L5"

        matched_genome3 = list1[2].genomes_dict["flat_file"]
        id3 = matched_genome3.filename
        expected_id3 = "D29"

        with self.subTest():
            self.assertEqual(id1, expected_id1)
        with self.subTest():
            self.assertEqual(id2, expected_id2)
        with self.subTest():
            self.assertEqual(id3, expected_id3)
        with self.subTest():
            self.assertEqual(len(eval_list), 0)

    def test_match_genomes_to_tickets2_8(self):
        """Verify that two tickets are matched to genomes,
        and two tickets are not matched (same identifier), using filename."""


        self.datagroup1.ticket.primary_phage_id = "Trixie"
        self.datagroup1.ticket.match_strategy = "filename"
        self.genome1.filename = "Trixie"

        self.datagroup2.ticket.primary_phage_id = "L5"
        self.datagroup2.ticket.match_strategy = "filename"
        self.genome2.filename = "L5"

        self.datagroup3.ticket.primary_phage_id = "D29"
        self.datagroup3.ticket.match_strategy = "filename"

        self.datagroup4.ticket.primary_phage_id = "D29"
        self.datagroup4.ticket.match_strategy = "filename"

        list1 = [self.datagroup1, self.datagroup2,
                self.datagroup3, self.datagroup4]
        list2 = [self.genome1, self.genome2]
        eval_list = \
            FunctionsTicket.match_genomes_to_tickets2(list1,
                                                    list2,
                                                    "flat_file")

        matched_genome1 = list1[0].genomes_dict["flat_file"]
        id1 = matched_genome1.filename
        expected_id1 = "Trixie"

        matched_genome2 = list1[1].genomes_dict["flat_file"]
        id2 = matched_genome2.filename
        expected_id2 = "L5"

        with self.subTest():
            self.assertEqual(id1, expected_id1)
        with self.subTest():
            self.assertEqual(id2, expected_id2)
        with self.subTest():
            self.assertEqual(len(list1[2].genomes_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(len(list1[3].genomes_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(len(eval_list), 1)

    def test_match_genomes_to_tickets2_9(self):
        """Verify that two genomes are matched to tickets,
        and two genomes are not matched (same identifier), using filename."""


        self.datagroup1.ticket.primary_phage_id = "Trixie"
        self.datagroup1.ticket.match_strategy = "filename"
        self.genome1.filename = "Trixie"

        self.datagroup2.ticket.primary_phage_id = "L5"
        self.datagroup2.ticket.match_strategy = "filename"
        self.genome2.filename = "L5"

        self.genome3.filename = "D29"
        self.genome4.filename = "D29"

        list1 = [self.datagroup1, self.datagroup2]
        list2 = [self.genome1, self.genome2, self.genome3, self.genome4]
        eval_list = \
            FunctionsTicket.match_genomes_to_tickets2(list1,
                                                    list2,
                                                    "flat_file")

        matched_genome1 = list1[0].genomes_dict["flat_file"]
        id1 = matched_genome1.filename
        expected_id1 = "Trixie"

        matched_genome2 = list1[1].genomes_dict["flat_file"]
        id2 = matched_genome2.filename
        expected_id2 = "L5"

        with self.subTest():
            self.assertEqual(id1, expected_id1)
        with self.subTest():
            self.assertEqual(id2, expected_id2)
        with self.subTest():
            self.assertEqual(len(eval_list), 1)
















class TestTicketFunctions4(unittest.TestCase):

    def setUp(self):

        self.ticket1 = Ticket.ImportTicket()
        self.ticket2 = Ticket.ImportTicket()
        self.ticket3 = Ticket.ImportTicket()
        self.ticket4 = Ticket.ImportTicket()

        self.datagroup1 = DataGroup.DataGroup()
        self.datagroup2 = DataGroup.DataGroup()
        self.datagroup3 = DataGroup.DataGroup()
        self.datagroup4 = DataGroup.DataGroup()

        self.datagroup1.ticket = self.ticket1
        self.datagroup2.ticket = self.ticket2
        self.datagroup3.ticket = self.ticket3
        self.datagroup4.ticket = self.ticket4


    def test_create_matched_object_dict_1(self):
        """Verify dictionary of update tickets is created."""

        self.datagroup1.ticket.type = "update"
        self.datagroup2.ticket.type = "update"
        self.datagroup3.ticket.type = "update"
        self.datagroup4.ticket.type = "update"
        list = [self.datagroup1, self.datagroup2,
                self.datagroup3, self.datagroup4]

        result_dict = FunctionsTicket.create_matched_object_dict(list)
        update_list = result_dict["update"]

        with self.subTest():
            self.assertEqual(len(result_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(len(update_list), 4)

    def test_create_matched_object_dict_2(self):
        """Verify dictionary of empty tickets is created."""

        list = [self.datagroup1, self.datagroup2,
                self.datagroup3, self.datagroup4]

        result_dict = FunctionsTicket.create_matched_object_dict(list)
        update_list = result_dict[""]

        with self.subTest():
            self.assertEqual(len(result_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(len(update_list), 4)


    def test_create_matched_object_dict_3(self):
        """Verify dictionary of multiple ticket types is created."""

        self.datagroup1.ticket.type = "update"
        self.datagroup2.ticket.type = "add"
        self.datagroup3.ticket.type = "update"
        self.datagroup4.ticket.type = "replace"

        list = [self.datagroup1, self.datagroup2,
                self.datagroup3, self.datagroup4]

        result_dict = FunctionsTicket.create_matched_object_dict(list)
        update_list = result_dict["update"]
        add_list = result_dict["add"]
        replace_list = result_dict["replace"]

        with self.subTest():
            self.assertEqual(len(result_dict.keys()), 3)
        with self.subTest():
            self.assertEqual(len(update_list), 2)
        with self.subTest():
            self.assertEqual(len(add_list), 1)
        with self.subTest():
            self.assertEqual(len(replace_list), 1)







if __name__ == '__main__':
    unittest.main()
