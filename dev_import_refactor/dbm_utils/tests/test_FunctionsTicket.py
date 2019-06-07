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

        self.genome_dict = {self.genome1.phage_id:self.genome1,
                            self.genome2.phage_id:self.genome2,
                            self.genome3.phage_id:self.genome3}

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


    def test_match_genomes_to_tickets_1(self):
        """Verify that one genome is matched to ticket correctly."""

        list1 = [self.datagroup1] # Trixie

        eval_list = \
            FunctionsTicket.match_genomes_to_tickets(list1,
                                                    self.genome_dict,
                                                    "phamerator")

        matched_genome = list1[0].genomes_dict["phamerator"]
        id = matched_genome.phage_id
        expected_id = "Trixie"
        with self.subTest():
            self.assertEqual(id, expected_id)
        with self.subTest():
            self.assertEqual(len(eval_list), 0)





    def test_match_genomes_to_tickets_2(self):
        """Verify that no genome is matched to empty ticket list."""

        list1 = []

        eval_list = \
            FunctionsTicket.match_genomes_to_tickets(list1,
                                                    self.genome_dict,
                                                    "phamerator")

        self.assertEqual(len(eval_list), 0)




    def test_match_genomes_to_tickets_3(self):
        """Verify that genome is not matched to ticket."""

        list1 = [self.datagroup4] # RedRock

        eval_list = \
            FunctionsTicket.match_genomes_to_tickets(list1,
                                                    self.genome_dict,
                                                    "phamerator")

        num_dict_keys = len(list1[0].genomes_dict.keys())

        with self.subTest():
            self.assertEqual(num_dict_keys, 0)
        with self.subTest():
            self.assertEqual(len(eval_list), 1)






    def test_match_genomes_to_tickets_4(self):
        """Verify that two genomes are matched to tickets correctly."""

        list1 = [self.datagroup1, self.datagroup2] # Trixie, L5

        eval_list = \
            FunctionsTicket.match_genomes_to_tickets(list1,
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




    def test_match_genomes_to_tickets_4(self):
        """Verify that first genome is not matched to ticket correctly."""

        list1 = [self.datagroup4, self.datagroup2] # RedRock, L5

        eval_list = \
            FunctionsTicket.match_genomes_to_tickets(list1,
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



    def test_match_genomes_to_tickets_5(self):
        """Verify that second genome is not matched to ticket correctly."""

        list1 = [self.datagroup2, self.datagroup4] # L5, RedRock

        eval_list = \
            FunctionsTicket.match_genomes_to_tickets(list1,
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





    def test_match_genomes_to_tickets_6(self):
        """Verify that two genomes are not matched to tickets."""

        list1 = [self.datagroup4, self.datagroup5] # RedRock, EagleEye

        eval_list = \
            FunctionsTicket.match_genomes_to_tickets(list1,
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









if __name__ == '__main__':
    unittest.main()
