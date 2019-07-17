""" Unit tests for misc. ticket functions"""

from classes import DataGroup
from classes import Genome
from classes import Ticket
from classes import Eval
from functions import tickets
import unittest





class TestTicketFunctions1(unittest.TestCase):


    def setUp(self):
        self.ticket = Ticket.GenomeTicket()

        self.normal_ticket_list = ["add",
                                "Trixie",
                                "Mycobacterium",
                                "A",
                                "A2",
                                "Final",
                                "Hatfull",
                                "Product",
                                "ABC123",
                                1,
                                1,
                                "phagesdb",
                                "none"]

        self.short_ticket_list = ["add",
                                "Trixie",
                                "Mycobacterium",
                                "A",
                                "A2",
                                "Final",
                                "Hatfull"]

        self.long_ticket_list = ["add",
                                "Trixie",
                                "Mycobacterium",
                                "A",
                                "A2",
                                "Final",
                                "Hatfull",
                                "Product",
                                "ABC123",
                                1,
                                1,
                                "phagesdb",
                                "none",
                                "extra"]

        self.normal_ticket_list2 = ["replace",
                                "KatherineG",
                                "Gordonia",
                                "A",
                                "A15",
                                "Final",
                                "Hatfull",
                                "Product",
                                "XYZ456",
                                1,
                                1,
                                "phagesdb",
                                "KatherineG"]

    def test_parse_import_ticket_1(self):
        """Verify properly structured data is parsed correctly."""
        tickets.parse_import_ticket(self.ticket, self.normal_ticket_list)

        with self.subTest():
            self.assertEqual(self.ticket.type, "add")
        with self.subTest():
            self.assertEqual(self.ticket.primary_phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(self.ticket.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.ticket.cluster, "A")
        with self.subTest():
            self.assertEqual(self.ticket.subcluster, "A2")
        with self.subTest():
            self.assertEqual(self.ticket.annotation_status, "final")
        with self.subTest():
            self.assertEqual(self.ticket.annotation_author, "hatfull")
        with self.subTest():
            self.assertEqual(self.ticket.annotation_qc, 1)
        with self.subTest():
            self.assertEqual(self.ticket.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(self.ticket.description_field, "product")
        with self.subTest():
            self.assertEqual(self.ticket.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.ticket.run_mode, "phagesdb")
        with self.subTest():
            self.assertEqual(self.ticket.secondary_phage_id, "none")

    def test_parse_import_ticket_2(self):
        """Verify improperly structured data is not parsed."""
        tickets.parse_import_ticket(self.ticket, self.short_ticket_list)
        self.assertEqual(self.ticket.type, "")

    def test_parse_import_ticket_3(self):
        """Verify improperly structured data is not parsed."""
        tickets.parse_import_ticket(self.ticket, self.long_ticket_list)
        self.assertEqual(self.ticket.type, "")







# TODO HERE IN PROGRESS fixing broken tickets functions.

    def test_parse_import_tickets_1(self):
        """Verify two lists of correct data are parsed."""
        list_of_lists = [self.normal_ticket_list,
                         self.normal_ticket_list2]
        list_of_tickets = tickets.parse_import_tickets(list_of_lists)
        with self.subTest():
            self.assertEqual(len(list_of_tickets), 2)
        with self.subTest():
            type = list_of_tickets[0].type
            self.assertEqual(type, "add")
        with self.subTest():
            type = list_of_tickets[1].type
            self.assertEqual(type, "replace")


    def test_parse_import_tickets_2(self):
        """Verify two lists of incorrect data are not parsed."""
        list_of_lists = [self.short_ticket_list,
                         self.long_ticket_list]
        list_of_tickets = tickets.parse_import_tickets(list_of_lists)
        with self.subTest():
            self.assertEqual(len(list_of_tickets), 2)
        with self.subTest():
            type = list_of_tickets[0].type
            self.assertEqual(type, "")
        with self.subTest():
            type = list_of_tickets[1].type
            self.assertEqual(type, "")

    def test_parse_import_tickets_3(self):
        """Verify mixed lists of correct and incorrect data are
        parsed correctly."""
        list_of_lists = [self.short_ticket_list,
                        self.normal_ticket_list,
                         self.long_ticket_list]
        list_of_tickets = tickets.parse_import_tickets(list_of_lists)
        with self.subTest():
            self.assertEqual(len(list_of_tickets), 3)
        with self.subTest():
            type = list_of_tickets[0].type
            self.assertEqual(type, "")
        with self.subTest():
            type = list_of_tickets[1].type
            self.assertEqual(type, "add")
        with self.subTest():
            type = list_of_tickets[2].type
            self.assertEqual(type, "")




    def test_compare_tickets_1(self):
        """Verify two tickets with no duplicates do not generate an error."""

        ticket1 = Ticket.GenomeTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "Trixie"
        ticket1.secondary_phage_id = "Trixie"
        ticket1.accession = "ABC123"

        ticket2 = Ticket.GenomeTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "L5"
        ticket2.secondary_phage_id = "L5"
        ticket2.accession = "EFG456"

        list_of_tickets = [ticket1, ticket2]
        tickets.compare_tickets(list_of_tickets)

        ticket1_errors = 0
        for eval in ticket1.evaluations:
            if eval.status == "error":
                ticket1_errors += 1

        ticket2_errors = 0
        for eval in ticket2.evaluations:
            if eval.status == "error":
                ticket2_errors += 1

        with self.subTest():
            self.assertEqual(len(ticket1.evaluations), 3)
        with self.subTest():
            self.assertEqual(len(ticket2.evaluations), 3)
        with self.subTest():
            self.assertEqual(ticket1_errors, 0)
        with self.subTest():
            self.assertEqual(ticket2_errors, 0)

    def test_compare_tickets_2(self):
        """Verify two tickets with 'none' duplicates do not generate an error."""

        ticket1 = Ticket.GenomeTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "none"
        ticket1.secondary_phage_id = "none"
        ticket1.accession = "none"

        ticket2 = Ticket.GenomeTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "none"
        ticket2.secondary_phage_id = "none"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        tickets.compare_tickets(list_of_tickets)

        ticket1_errors = 0
        for eval in ticket1.evaluations:
            if eval.status == "error":
                ticket1_errors += 1

        ticket2_errors = 0
        for eval in ticket2.evaluations:
            if eval.status == "error":
                ticket2_errors += 1

        with self.subTest():
            self.assertEqual(ticket1_errors, 0)
        with self.subTest():
            self.assertEqual(ticket2_errors, 0)

    def test_compare_tickets_3(self):
        """Verify two tickets with Primary Phage ID duplicates
        do generate an error."""

        ticket1 = Ticket.GenomeTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "Trixie"
        ticket1.secondary_phage_id = "none"
        ticket1.accession = "none"

        ticket2 = Ticket.GenomeTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "Trixie"
        ticket2.secondary_phage_id = "none"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        tickets.compare_tickets(list_of_tickets)

        ticket1_errors = 0
        for eval in ticket1.evaluations:
            if eval.status == "error":
                ticket1_errors += 1

        ticket2_errors = 0
        for eval in ticket2.evaluations:
            if eval.status == "error":
                ticket2_errors += 1

        with self.subTest():
            self.assertEqual(ticket1_errors, 1)
        with self.subTest():
            self.assertEqual(ticket2_errors, 1)

    def test_compare_tickets_4(self):
        """Verify two tickets with Secondary Phage ID duplicates
        do generate an error."""

        ticket1 = Ticket.GenomeTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "none"
        ticket1.secondary_phage_id = "Trixie"
        ticket1.accession = "none"

        ticket2 = Ticket.GenomeTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "none"
        ticket2.secondary_phage_id = "Trixie"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        tickets.compare_tickets(list_of_tickets)

        ticket1_errors = 0
        for eval in ticket1.evaluations:
            if eval.status == "error":
                ticket1_errors += 1

        ticket2_errors = 0
        for eval in ticket2.evaluations:
            if eval.status == "error":
                ticket2_errors += 1

        with self.subTest():
            self.assertEqual(ticket1_errors, 1)
        with self.subTest():
            self.assertEqual(ticket2_errors, 1)

    def test_compare_tickets_5(self):
        """Verify two tickets with Accession duplicates do generate an error."""

        ticket1 = Ticket.GenomeTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "none"
        ticket1.secondary_phage_id = "none"
        ticket1.accession = "ABC123"

        ticket2 = Ticket.GenomeTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "none"
        ticket2.secondary_phage_id = "none"
        ticket2.accession = "ABC123"

        list_of_tickets = [ticket1, ticket2]
        tickets.compare_tickets(list_of_tickets)

        ticket1_errors = 0
        for eval in ticket1.evaluations:
            if eval.status == "error":
                ticket1_errors += 1

        ticket2_errors = 0
        for eval in ticket2.evaluations:
            if eval.status == "error":
                ticket2_errors += 1

        with self.subTest():
            self.assertEqual(ticket1_errors, 1)
        with self.subTest():
            self.assertEqual(ticket2_errors, 1)

    def test_compare_tickets_6(self):
        """Verify two tickets with multiple duplicates
        do generate multiple errors."""

        ticket1 = Ticket.GenomeTicket()
        ticket1.type = "replace"
        ticket1.primary_phage_id = "Trixie"
        ticket1.secondary_phage_id = "none"
        ticket1.accession = "ABC123"

        ticket2 = Ticket.GenomeTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "Trixie"
        ticket2.secondary_phage_id = "none"
        ticket2.accession = "ABC123"

        list_of_tickets = [ticket1, ticket2]
        tickets.compare_tickets(list_of_tickets)

        ticket1_errors = 0
        for eval in ticket1.evaluations:
            if eval.status == "error":
                ticket1_errors += 1

        ticket2_errors = 0
        for eval in ticket2.evaluations:
            if eval.status == "error":
                ticket2_errors += 1

        with self.subTest():
            self.assertEqual(ticket1_errors, 2)
        with self.subTest():
            self.assertEqual(ticket2_errors, 2)

    def test_compare_tickets_7(self):
        """Verify a conflict between a Primary Phage ID and a
        Secondary Phage ID does generate an error."""

        ticket1 = Ticket.GenomeTicket()
        ticket1.type = "add"
        ticket1.primary_phage_id = "Trixie"
        ticket1.secondary_phage_id = "none"
        ticket1.accession = "none"

        ticket2 = Ticket.GenomeTicket()
        ticket2.type = "replace"
        ticket2.primary_phage_id = "L5"
        ticket2.secondary_phage_id = "Trixie"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        tickets.compare_tickets(list_of_tickets)

        ticket1_errors = 0
        for eval in ticket1.evaluations:
            if eval.status == "error":
                ticket1_errors += 1

        ticket2_errors = 0
        for eval in ticket2.evaluations:
            if eval.status == "error":
                ticket2_errors += 1

        with self.subTest():
            self.assertEqual(ticket1_errors, 1)
        with self.subTest():
            self.assertEqual(ticket2_errors, 1)




class TestTicketFunctions2(unittest.TestCase):

    def setUp(self):

        self.ticket1 = Ticket.GenomeTicket()
        self.ticket2 = Ticket.GenomeTicket()

        self.ticket1.primary_phage_id = "Trixie"
        self.ticket2.primary_phage_id = "L5"

        self.datagroup1 = DataGroup.DataGroup()
        self.datagroup2 = DataGroup.DataGroup()

        self.datagroup1.ticket = self.ticket1
        self.datagroup2.ticket = self.ticket2


    # TODO no longer needed probably.
    # def test_assign_match_strategy_1(self):
    #     """Verify strategy is assigned with no error produced."""
    #     input_strategy = "phage_id"
    #     self.datagroup1.ticket.match_strategy = input_strategy
    #     self.datagroup2.ticket.match_strategy = input_strategy
    #     list1 = [self.datagroup1, self.datagroup2] # Trixie, L5
    #     output_strategy, eval_result = \
    #         tickets.assign_match_strategy(list1)
    #     with self.subTest():
    #         self.assertEqual(output_strategy, input_strategy)
    #     with self.subTest():
    #         self.assertEqual(eval_result.status, "correct")
    #
    # def test_assign_match_strategy_2(self):
    #     """Verify no strategy is assigned and an error is produced."""
    #     input_strategy1 = "phage_id"
    #     input_strategy2 = "filename"
    #     self.datagroup1.ticket.match_strategy = input_strategy1
    #     self.datagroup2.ticket.match_strategy = input_strategy2
    #     list1 = [self.datagroup1, self.datagroup2] # Trixie, L5
    #     output_strategy, eval_result = \
    #         tickets.assign_match_strategy(list1)
    #     expected_strategy = ""
    #     with self.subTest():
    #         self.assertEqual(output_strategy, expected_strategy)
    #     with self.subTest():
    #         self.assertEqual(eval_result.status, "error")




class TestTicketFunctions3(unittest.TestCase):

    def setUp(self):

        self.ticket1 = Ticket.GenomeTicket()
        self.ticket2 = Ticket.GenomeTicket()
        self.ticket3 = Ticket.GenomeTicket()
        self.ticket4 = Ticket.GenomeTicket()

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

        result_dict = tickets.create_matched_object_dict(list)
        update_list = result_dict["update"]

        with self.subTest():
            self.assertEqual(len(result_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(len(update_list), 4)

    def test_create_matched_object_dict_2(self):
        """Verify dictionary of empty tickets is created."""

        list = [self.datagroup1, self.datagroup2,
                self.datagroup3, self.datagroup4]

        result_dict = tickets.create_matched_object_dict(list)
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

        result_dict = tickets.create_matched_object_dict(list)
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











class TestTicketFunctions4(unittest.TestCase):

    def setUp(self):


        # Add ticket.
        self.add_ticket = Ticket.GenomeTicket()
        self.add_ticket.type = "add"
        self.add_ticket.primary_phage_id = "Trixie_Draft"
        self.add_ticket.run_mode = "phagesdb"
        self.add_ticket.description_field = "product"
        self.add_ticket.host_genus = "Mycobacterium smegmatis"
        self.add_ticket.cluster = "A"
        self.add_ticket.subcluster = "A2"
        self.add_ticket.annotation_status = "final"
        self.add_ticket.annotation_author = "hatfull"
        self.add_ticket.annotation_qc = 1
        self.add_ticket.retrieve_record = 1
        self.add_ticket.accession = "ABC123.1"
        self.add_ticket.secondary_phage_id = "none"
        self.datagroup1 = DataGroup.DataGroup()
        self.datagroup1.ticket = self.add_ticket

        # Remove ticket.
        self.remove_ticket = Ticket.GenomeTicket()
        self.remove_ticket.type = "replace"
        self.remove_ticket.primary_phage_id = "Trixie_Draft"
        self.remove_ticket.run_mode = "phagesdb"
        self.remove_ticket.description_field = "product"
        self.remove_ticket.host_genus = "Mycobacterium smegmatis"
        self.remove_ticket.cluster = "A"
        self.remove_ticket.subcluster = "A2"
        self.remove_ticket.annotation_status = "final"
        self.remove_ticket.annotation_author = "hatfull"
        self.remove_ticket.annotation_qc = 1
        self.remove_ticket.retrieve_record = 1
        self.remove_ticket.accession = "ABC123.1"
        self.remove_ticket.secondary_phage_id = "L5"
        self.datagroup2 = DataGroup.DataGroup()
        self.datagroup2.ticket = self.remove_ticket




        # Invalid ticket.
        self.invalid_ticket = Ticket.GenomeTicket()
        self.invalid_ticket.type = "invalid"
        self.invalid_ticket.primary_phage_id = "Trixie_Draft"
        self.invalid_ticket.run_mode = "phagesdb"
        self.invalid_ticket.description_field = "product"
        self.invalid_ticket.host_genus = "Mycobacterium smegmatis"
        self.invalid_ticket.cluster = "A"
        self.invalid_ticket.subcluster = "A2"
        self.invalid_ticket.annotation_status = "final"
        self.invalid_ticket.annotation_author = "hatfull"
        self.invalid_ticket.annotation_qc = 1
        self.invalid_ticket.retrieve_record = 1
        self.invalid_ticket.accession = "ABC123.1"
        self.invalid_ticket.secondary_phage_id = "L5"
        self.datagroup3 = DataGroup.DataGroup()
        self.datagroup3.ticket = self.invalid_ticket


    def test_copy_ticket_to_genome_1(self):
        """Verify data from 'add' ticket is added to genome."""

        tickets.copy_ticket_to_genome(self.datagroup1)
        matched_genome = self.datagroup1.genome_dict["add"]
        with self.subTest():
            self.assertEqual(matched_genome.phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(matched_genome.phage_name, "Trixie_Draft")
        with self.subTest():
            self.assertEqual(matched_genome.type, "add")
        with self.subTest():
            self.assertEqual(matched_genome.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(matched_genome.cluster, "A")
        with self.subTest():
            self.assertEqual(matched_genome.subcluster, "A2")
        with self.subTest():
            self.assertEqual(matched_genome.cluster_subcluster, "A2")
        with self.subTest():
            self.assertEqual(matched_genome.annotation_status, "final")
        with self.subTest():
            self.assertEqual(matched_genome.annotation_author, 1)
        with self.subTest():
            self.assertEqual(matched_genome.annotation_qc, 1)
        with self.subTest():
            self.assertEqual(matched_genome.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(matched_genome.accession, "ABC123")

    def test_copy_ticket_to_genome_2(self):
        """Verify data from 'remove' ticket is added to genome."""

        tickets.copy_ticket_to_genome(self.datagroup2)
        matched_genome1 = self.datagroup2.genome_dict["add"]
        matched_genome2 = self.datagroup2.genome_dict["remove"]

        with self.subTest():
            self.assertEqual(matched_genome1.phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(matched_genome1.phage_name, "Trixie_Draft")
        with self.subTest():
            self.assertEqual(matched_genome2.phage_id, "L5")
        with self.subTest():
            self.assertEqual(matched_genome2.phage_name, "")

    def test_copy_ticket_to_genome_3(self):
        """Verify data from 'invalid' ticket is not added to genome."""

        tickets.copy_ticket_to_genome(self.datagroup3)
        self.assertEqual(len(self.datagroup3.genome_dict.keys()), 0)







# TODO below code is broken since I have changed how eval is structured.
#
#
#
#
# class TestTicketFunctions5(unittest.TestCase):
#
#     def setUp(self):
#
#         self.genome1 = Genome.Genome()
#         self.genome2 = Genome.Genome()
#         self.genome3 = Genome.Genome()
#         self.genome4 = Genome.Genome()
#
#         self.ticket1 = Ticket.GenomeTicket()
#         self.ticket2 = Ticket.GenomeTicket()
#         self.ticket3 = Ticket.GenomeTicket()
#         self.ticket4 = Ticket.GenomeTicket()
#
#         self.datagroup1 = DataGroup.DataGroup()
#         self.datagroup2 = DataGroup.DataGroup()
#         self.datagroup3 = DataGroup.DataGroup()
#         self.datagroup4 = DataGroup.DataGroup()
#
#         self.datagroup1.ticket = self.ticket1
#         self.datagroup2.ticket = self.ticket2
#         self.datagroup3.ticket = self.ticket3
#         self.datagroup4.ticket = self.ticket4
#
#
#     def test_match_genomes_to_tickets_1(self):
#         """Verify that one genome is matched to ticket using phage_id."""
#
#         self.datagroup1.ticket.primary_phage_id = "Trixie"
#         self.datagroup1.ticket.match_strategy = "phage_id"
#         self.genome1.phage_id = "Trixie"
#
#
#         list1 = [self.datagroup1] # Trixie
#         list2 = [self.genome1] # Trixie
#         eval_list = \
#             tickets.match_genomes_to_tickets(list1,
#                                                     list2,
#                                                     "phamerator")
#
#         matched_genome = list1[0].genome_dict["phamerator"]
#         id = matched_genome.phage_id
#         expected_id = "Trixie"
#         with self.subTest():
#             self.assertEqual(id, expected_id)
#         with self.subTest():
#             self.assertEqual(len(eval_list), 0)
#
#     def test_match_genomes_to_tickets_2(self):
#         """Verify that one genome is matched to ticket using filename."""
#
#         self.datagroup1.ticket.primary_phage_id = "Trixie"
#         self.datagroup1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#
#         list1 = [self.datagroup1] # Trixie
#         list2 = [self.genome1] # Trixie
#         eval_list = \
#             tickets.match_genomes_to_tickets(list1,
#                                                     list2,
#                                                     "flat_file")
#
#         matched_genome = list1[0].genome_dict["flat_file"]
#         id = matched_genome.filename
#         expected_id = "Trixie"
#         with self.subTest():
#             self.assertEqual(id, expected_id)
#         with self.subTest():
#             self.assertEqual(len(eval_list), 0)
#
#     def test_match_genomes_to_tickets_3(self):
#         """Verify that one genome is matched to ticket,
#         and one genome is not matched (no matching ticket), using filename."""
#
#         self.datagroup1.ticket.primary_phage_id = "Trixie"
#         self.datagroup1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#         self.genome2.filename = "L5"
#
#
#         list1 = [self.datagroup1] # Trixie
#         list2 = [self.genome1, self.genome2] # Trixie, L5
#         eval_list = \
#             tickets.match_genomes_to_tickets(list1,
#                                                     list2,
#                                                     "flat_file")
#
#         matched_genome = list1[0].genome_dict["flat_file"]
#         id = matched_genome.filename
#         expected_id = "Trixie"
#         with self.subTest():
#             self.assertEqual(id, expected_id)
#         with self.subTest():
#             self.assertEqual(len(eval_list), 1)
#
#     def test_match_genomes_to_tickets_4(self):
#         """Verify that one ticket is matched to genome,
#         and one ticket is not matched (no matching genome), using filename."""
#
#         self.datagroup1.ticket.primary_phage_id = "Trixie"
#         self.datagroup1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#
#         self.datagroup2.ticket.primary_phage_id = "L5"
#         self.datagroup2.ticket.match_strategy = "filename"
#
#         list1 = [self.datagroup1, self.datagroup2] # Trixie, L5
#         list2 = [self.genome1] # Trixie
#
#         eval_list = \
#             tickets.match_genomes_to_tickets(list1,
#                                                     list2,
#                                                     "flat_file")
#
#         matched_genome = list1[0].genome_dict["flat_file"]
#         id = matched_genome.filename
#         expected_id = "Trixie"
#         with self.subTest():
#             self.assertEqual(id, expected_id)
#         with self.subTest():
#             self.assertEqual(len(eval_list), 1)
#
#     def test_match_genomes_to_tickets_5(self):
#         """Verify that two genomes are matched to tickets,
#         using filename."""
#
#         self.datagroup1.ticket.primary_phage_id = "Trixie"
#         self.datagroup1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#
#         self.datagroup2.ticket.primary_phage_id = "L5"
#         self.datagroup2.ticket.match_strategy = "filename"
#         self.genome2.filename = "L5"
#
#
#         list1 = [self.datagroup1, self.datagroup2] # Trixie, L5
#         list2 = [self.genome1, self.genome2] # Trixie, L5
#         eval_list = \
#             tickets.match_genomes_to_tickets(list1,
#                                                     list2,
#                                                     "flat_file")
#
#         matched_genome1 = list1[0].genome_dict["flat_file"]
#         id1 = matched_genome1.filename
#         expected_id1 = "Trixie"
#
#         matched_genome2 = list1[1].genome_dict["flat_file"]
#         id2 = matched_genome2.filename
#         expected_id2 = "L5"
#
#         with self.subTest():
#             self.assertEqual(id1, expected_id1)
#         with self.subTest():
#             self.assertEqual(id2, expected_id2)
#         with self.subTest():
#             self.assertEqual(len(eval_list), 0)
#
#     def test_match_genomes_to_tickets_6(self):
#         """Verify that no genomes are matched to tickets due to
#         conflicting strategies."""
#
#         self.datagroup1.ticket.primary_phage_id = "Trixie"
#         self.datagroup1.ticket.match_strategy = "phage_id"
#         self.genome1.filename = "Trixie"
#
#         self.datagroup2.ticket.primary_phage_id = "L5"
#         self.datagroup2.ticket.match_strategy = "filename"
#         self.genome2.filename = "L5"
#
#         list1 = [self.datagroup1, self.datagroup2] # Trixie, L5
#         list2 = [self.genome1, self.genome2] # Trixie, L5
#         eval_list = \
#             tickets.match_genomes_to_tickets(list1,
#                                                     list2,
#                                                     "flat_file")
#
#         with self.subTest():
#             self.assertEqual(len(list1[0].genome_dict.keys()), 0)
#         with self.subTest():
#             self.assertEqual(len(list1[0].genome_dict.keys()), 0)
#         with self.subTest():
#             self.assertEqual(len(eval_list), 3)
#
#     def test_match_genomes_to_tickets_7(self):
#         """Verify that three genomes are matched to tickets,
#         using filename."""
#
#         self.datagroup1.ticket.primary_phage_id = "Trixie"
#         self.datagroup1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#
#         self.datagroup2.ticket.primary_phage_id = "L5"
#         self.datagroup2.ticket.match_strategy = "filename"
#         self.genome2.filename = "L5"
#
#         self.datagroup3.ticket.primary_phage_id = "D29"
#         self.datagroup3.ticket.match_strategy = "filename"
#         self.genome3.filename = "D29"
#
#         list1 = [self.datagroup1, self.datagroup2, self.datagroup3]
#         list2 = [self.genome1, self.genome2, self.genome3]
#         eval_list = \
#             tickets.match_genomes_to_tickets(list1,
#                                                     list2,
#                                                     "flat_file")
#
#         matched_genome1 = list1[0].genome_dict["flat_file"]
#         id1 = matched_genome1.filename
#         expected_id1 = "Trixie"
#
#         matched_genome2 = list1[1].genome_dict["flat_file"]
#         id2 = matched_genome2.filename
#         expected_id2 = "L5"
#
#         matched_genome3 = list1[2].genome_dict["flat_file"]
#         id3 = matched_genome3.filename
#         expected_id3 = "D29"
#
#         with self.subTest():
#             self.assertEqual(id1, expected_id1)
#         with self.subTest():
#             self.assertEqual(id2, expected_id2)
#         with self.subTest():
#             self.assertEqual(id3, expected_id3)
#         with self.subTest():
#             self.assertEqual(len(eval_list), 0)
#
#     def test_match_genomes_to_tickets_8(self):
#         """Verify that two tickets are matched to genomes,
#         and two tickets are not matched (same identifier), using filename."""
#
#
#         self.datagroup1.ticket.primary_phage_id = "Trixie"
#         self.datagroup1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#
#         self.datagroup2.ticket.primary_phage_id = "L5"
#         self.datagroup2.ticket.match_strategy = "filename"
#         self.genome2.filename = "L5"
#
#         self.datagroup3.ticket.primary_phage_id = "D29"
#         self.datagroup3.ticket.match_strategy = "filename"
#
#         self.datagroup4.ticket.primary_phage_id = "D29"
#         self.datagroup4.ticket.match_strategy = "filename"
#
#         list1 = [self.datagroup1, self.datagroup2,
#                 self.datagroup3, self.datagroup4]
#         list2 = [self.genome1, self.genome2]
#         eval_list = \
#             tickets.match_genomes_to_tickets(list1,
#                                                     list2,
#                                                     "flat_file")
#
#         matched_genome1 = list1[0].genome_dict["flat_file"]
#         id1 = matched_genome1.filename
#         expected_id1 = "Trixie"
#
#         matched_genome2 = list1[1].genome_dict["flat_file"]
#         id2 = matched_genome2.filename
#         expected_id2 = "L5"
#
#         with self.subTest():
#             self.assertEqual(id1, expected_id1)
#         with self.subTest():
#             self.assertEqual(id2, expected_id2)
#         with self.subTest():
#             self.assertEqual(len(list1[2].genome_dict.keys()), 0)
#         with self.subTest():
#             self.assertEqual(len(list1[3].genome_dict.keys()), 0)
#         with self.subTest():
#             self.assertEqual(len(eval_list), 1)
#
#     def test_match_genomes_to_tickets_9(self):
#         """Verify that two genomes are matched to tickets,
#         and two genomes are not matched (same identifier), using filename."""
#
#
#         self.datagroup1.ticket.primary_phage_id = "Trixie"
#         self.datagroup1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#
#         self.datagroup2.ticket.primary_phage_id = "L5"
#         self.datagroup2.ticket.match_strategy = "filename"
#         self.genome2.filename = "L5"
#
#         self.genome3.filename = "D29"
#         self.genome4.filename = "D29"
#
#         list1 = [self.datagroup1, self.datagroup2]
#         list2 = [self.genome1, self.genome2, self.genome3, self.genome4]
#         eval_list = \
#             tickets.match_genomes_to_tickets(list1,
#                                                     list2,
#                                                     "flat_file")
#
#         matched_genome1 = list1[0].genome_dict["flat_file"]
#         id1 = matched_genome1.filename
#         expected_id1 = "Trixie"
#
#         matched_genome2 = list1[1].genome_dict["flat_file"]
#         id2 = matched_genome2.filename
#         expected_id2 = "L5"
#
#         with self.subTest():
#             self.assertEqual(id1, expected_id1)
#         with self.subTest():
#             self.assertEqual(id2, expected_id2)
#         with self.subTest():
#             self.assertEqual(len(eval_list), 1)
#
#
#
#
#
#
#
#
#
#
#
#
#
#




if __name__ == '__main__':
    unittest.main()
