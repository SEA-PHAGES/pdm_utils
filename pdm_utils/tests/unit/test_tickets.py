""" Unit tests for misc. ticket functions"""

from classes import bundle
from classes import genome
from classes import ticket
from classes import eval
from functions import tickets
import unittest





class TestTicketFunctions1(unittest.TestCase):


    def setUp(self):
        self.tkt = ticket.GenomeTicket()

        self.normal_ticket_list = ["add",
                                "Trixie",
                                "Mycobacterium",
                                "A",
                                "A2",
                                "Final",
                                "Hatfull",
                                "Product",
                                "ABC123",
                                "1",
                                "1",
                                "phagesdb"]

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
                                "1",
                                "1",
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
                                "1",
                                "1",
                                "phagesdb"]


        self.empty_ticket_list = [None] * 12

        self.filled_ticket = ticket.GenomeTicket()
        self.filled_ticket.type = "add"
        self.filled_ticket.phage_id = "Trixie"
        self.filled_ticket.host_genus = "Mycobacterium"
        self.filled_ticket.cluster = "A"
        self.filled_ticket.subcluster = "A2"
        self.filled_ticket.annotation_status = "final"
        self.filled_ticket.annotation_author = "hatfull"
        self.filled_ticket.annotation_qc = 1
        self.filled_ticket.retrieve_record = 1
        self.filled_ticket.description_field = "product"
        self.filled_ticket.accession = "ABC123"
        self.filled_ticket.run_mode = "phagesdb"




    # def test_parse_import_ticket_1(self):
    #     """Verify properly structured data is parsed correctly."""
    #     tickets.parse_import_ticket(self.tkt, self.normal_ticket_list)
    #
    #     with self.subTest():
    #         self.assertEqual(self.tkt.type, "add")
    #     with self.subTest():
    #         self.assertEqual(self.tkt.phage_id, "Trixie")
    #     with self.subTest():
    #         self.assertEqual(self.tkt.host_genus, "Mycobacterium")
    #     with self.subTest():
    #         self.assertEqual(self.tkt.cluster, "A")
    #     with self.subTest():
    #         self.assertEqual(self.tkt.subcluster, "A2")
    #     with self.subTest():
    #         self.assertEqual(self.tkt.annotation_status, "final")
    #     with self.subTest():
    #         self.assertEqual(self.tkt.annotation_author, "hatfull")
    #     with self.subTest():
    #         self.assertEqual(self.tkt.annotation_qc, 1)
    #     with self.subTest():
    #         self.assertEqual(self.tkt.retrieve_record, 1)
    #     with self.subTest():
    #         self.assertEqual(self.tkt.description_field, "product")
    #     with self.subTest():
    #         self.assertEqual(self.tkt.accession, "ABC123")
    #     with self.subTest():
    #         self.assertEqual(self.tkt.run_mode, "phagesdb")
    #
    # def test_parse_import_ticket_2(self):
    #     """Verify id is set appropriately."""
    #     tickets.parse_import_ticket(self.tkt,
    #                                 self.normal_ticket_list,
    #                                 id = "TicketXYZ")
    #     self.assertEqual(self.tkt.id, "TicketXYZ")
    #
    # def test_parse_import_ticket_3(self):
    #     """Verify improperly structured data is not parsed."""
    #     tickets.parse_import_ticket(self.tkt, self.short_ticket_list)
    #     self.assertEqual(self.tkt.type, "")
    #
    # def test_parse_import_ticket_4(self):
    #     """Verify improperly structured data is not parsed."""
    #     tickets.parse_import_ticket(self.tkt, self.long_ticket_list)
    #     self.assertEqual(self.tkt.type, "")








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
        with self.subTest():
            self.assertEqual(list_of_tickets[0].id, 1)
        with self.subTest():
            self.assertEqual(list_of_tickets[1].id, 2)


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



















###

    def test_parse_import_ticket_data_1(self):
        """Verify properly structured data is parsed correctly."""
        tickets.parse_import_ticket_data(self.tkt, self.normal_ticket_list)

        with self.subTest():
            self.assertEqual(self.tkt.type, "add")
        with self.subTest():
            self.assertEqual(self.tkt.phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(self.tkt.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.tkt.cluster, "A")
        with self.subTest():
            self.assertEqual(self.tkt.subcluster, "A2")
        with self.subTest():
            self.assertEqual(self.tkt.annotation_status, "final")
        with self.subTest():
            self.assertEqual(self.tkt.annotation_author, "hatfull")
        with self.subTest():
            self.assertEqual(self.tkt.annotation_qc, 1)
        with self.subTest():
            self.assertEqual(self.tkt.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(self.tkt.description_field, "product")
        with self.subTest():
            self.assertEqual(self.tkt.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.tkt.run_mode, "phagesdb")

    def test_parse_import_ticket_data_2(self):
        """Verify id is set appropriately."""
        tickets.parse_import_ticket_data(self.tkt,
                                    self.normal_ticket_list,
                                    id = "TicketXYZ")
        self.assertEqual(self.tkt.id, "TicketXYZ")

    def test_parse_import_ticket_data_3(self):
        """Verify improperly structured data is not parsed."""
        tickets.parse_import_ticket_data(self.tkt, self.short_ticket_list)
        self.assertEqual(self.tkt.type, "")

    def test_parse_import_ticket_data_4(self):
        """Verify improperly structured data is not parsed."""
        tickets.parse_import_ticket_data(self.tkt, self.long_ticket_list)
        self.assertEqual(self.tkt.type, "")






    def test_parse_import_ticket_data_5(self):
        """Verify properly structured data is parsed correctly."""

        tickets.parse_import_ticket_data(self.filled_ticket,
                                         self.empty_ticket_list,
                                         direction = "ticket_to_list")

        with self.subTest():
            self.assertEqual(self.empty_ticket_list[0], "add")
        with self.subTest():
            self.assertEqual(self.empty_ticket_list[1], "Trixie")
        with self.subTest():
            self.assertEqual(self.empty_ticket_list[2], "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.empty_ticket_list[3], "A")
        with self.subTest():
            self.assertEqual(self.empty_ticket_list[4], "A2")
        with self.subTest():
            self.assertEqual(self.empty_ticket_list[5], "final")
        with self.subTest():
            self.assertEqual(self.empty_ticket_list[6], "hatfull")
        with self.subTest():
            self.assertEqual(self.empty_ticket_list[9], "1")
        with self.subTest():
            self.assertEqual(self.empty_ticket_list[10], "1")
        with self.subTest():
            self.assertEqual(self.empty_ticket_list[7], "product")
        with self.subTest():
            self.assertEqual(self.empty_ticket_list[8], "ABC123")
        with self.subTest():
            self.assertEqual(self.empty_ticket_list[11], "phagesdb")




    def test_parse_import_ticket_data_6(self):
        """Verify no changes are made if the direction is invalid."""

        tickets.parse_import_ticket_data(self.filled_ticket,
                                         self.empty_ticket_list,
                                         direction = "invalid")

        with self.subTest():
            self.assertIsNone(self.empty_ticket_list[0])
        with self.subTest():
            self.assertEqual(self.filled_ticket.type, "add")



###









    def test_compare_tickets_1(self):
        """Verify two tickets with no duplicates do not generate an error."""

        ticket1 = ticket.GenomeTicket()
        ticket1.type = "replace"
        ticket1.phage_id = "Trixie"
        ticket1.accession = "ABC123"

        ticket2 = ticket.GenomeTicket()
        ticket2.type = "replace"
        ticket2.phage_id = "L5"
        ticket2.accession = "EFG456"

        list_of_tickets = [ticket1, ticket2]
        tickets.compare_tickets(list_of_tickets)

        ticket1_errors = 0
        for evl in ticket1.evaluations:
            if evl.status == "error":
                ticket1_errors += 1

        ticket2_errors = 0
        for evl in ticket2.evaluations:
            if evl.status == "error":
                ticket2_errors += 1

        with self.subTest():
            self.assertEqual(len(ticket1.evaluations), 2)
        with self.subTest():
            self.assertEqual(len(ticket2.evaluations), 2)
        with self.subTest():
            self.assertEqual(ticket1_errors, 0)
        with self.subTest():
            self.assertEqual(ticket2_errors, 0)

    def test_compare_tickets_2(self):
        """Verify two tickets with 'none' duplicates do not generate an error."""

        ticket1 = ticket.GenomeTicket()
        ticket1.type = "replace"
        ticket1.phage_id = "none"
        ticket1.accession = "none"

        ticket2 = ticket.GenomeTicket()
        ticket2.type = "replace"
        ticket2.phage_id = "none"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        tickets.compare_tickets(list_of_tickets)

        ticket1_errors = 0
        for evl in ticket1.evaluations:
            if evl.status == "error":
                ticket1_errors += 1

        ticket2_errors = 0
        for evl in ticket2.evaluations:
            if evl.status == "error":
                ticket2_errors += 1

        with self.subTest():
            self.assertEqual(ticket1_errors, 0)
        with self.subTest():
            self.assertEqual(ticket2_errors, 0)

    def test_compare_tickets_3(self):
        """Verify two tickets with Primary Phage ID duplicates
        do generate an error."""

        ticket1 = ticket.GenomeTicket()
        ticket1.type = "replace"
        ticket1.phage_id = "Trixie"
        ticket1.accession = "none"

        ticket2 = ticket.GenomeTicket()
        ticket2.type = "replace"
        ticket2.phage_id = "Trixie"
        ticket2.accession = "none"

        list_of_tickets = [ticket1, ticket2]
        tickets.compare_tickets(list_of_tickets)

        ticket1_errors = 0
        for evl in ticket1.evaluations:
            if evl.status == "error":
                ticket1_errors += 1

        ticket2_errors = 0
        for evl in ticket2.evaluations:
            if evl.status == "error":
                ticket2_errors += 1

        with self.subTest():
            self.assertEqual(ticket1_errors, 1)
        with self.subTest():
            self.assertEqual(ticket2_errors, 1)



    def test_compare_tickets_5(self):
        """Verify two tickets with Accession duplicates do generate an error."""

        ticket1 = ticket.GenomeTicket()
        ticket1.type = "replace"
        ticket1.phage_id = "none"
        ticket1.accession = "ABC123"

        ticket2 = ticket.GenomeTicket()
        ticket2.type = "replace"
        ticket2.phage_id = "none"
        ticket2.accession = "ABC123"

        list_of_tickets = [ticket1, ticket2]
        tickets.compare_tickets(list_of_tickets)

        ticket1_errors = 0
        for evl in ticket1.evaluations:
            if evl.status == "error":
                ticket1_errors += 1

        ticket2_errors = 0
        for evl in ticket2.evaluations:
            if evl.status == "error":
                ticket2_errors += 1

        with self.subTest():
            self.assertEqual(ticket1_errors, 1)
        with self.subTest():
            self.assertEqual(ticket2_errors, 1)

    def test_compare_tickets_6(self):
        """Verify two tickets with multiple duplicates
        do generate multiple errors."""

        ticket1 = ticket.GenomeTicket()
        ticket1.type = "replace"
        ticket1.phage_id = "Trixie"
        ticket1.accession = "ABC123"

        ticket2 = ticket.GenomeTicket()
        ticket2.type = "replace"
        ticket2.phage_id = "Trixie"
        ticket2.accession = "ABC123"

        list_of_tickets = [ticket1, ticket2]
        tickets.compare_tickets(list_of_tickets)

        ticket1_errors = 0
        for evl in ticket1.evaluations:
            if evl.status == "error":
                ticket1_errors += 1

        ticket2_errors = 0
        for evl in ticket2.evaluations:
            if evl.status == "error":
                ticket2_errors += 1

        with self.subTest():
            self.assertEqual(ticket1_errors, 2)
        with self.subTest():
            self.assertEqual(ticket2_errors, 2)





class TestTicketFunctions2(unittest.TestCase):

    def setUp(self):

        self.ticket1 = ticket.GenomeTicket()
        self.ticket2 = ticket.GenomeTicket()

        self.ticket1.phage_id = "Trixie"
        self.ticket2.phage_id = "L5"

        self.bundle1 = bundle.Bundle()
        self.bundle2 = bundle.Bundle()

        self.bundle1.ticket = self.ticket1
        self.bundle2.ticket = self.ticket2


    # TODO no longer needed probably.
    # def test_assign_match_strategy_1(self):
    #     """Verify strategy is assigned with no error produced."""
    #     input_strategy = "phage_id"
    #     self.bundle1.ticket.match_strategy = input_strategy
    #     self.bundle2.ticket.match_strategy = input_strategy
    #     list1 = [self.bundle1, self.bundle2] # Trixie, L5
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
    #     self.bundle1.ticket.match_strategy = input_strategy1
    #     self.bundle2.ticket.match_strategy = input_strategy2
    #     list1 = [self.bundle1, self.bundle2] # Trixie, L5
    #     output_strategy, eval_result = \
    #         tickets.assign_match_strategy(list1)
    #     expected_strategy = ""
    #     with self.subTest():
    #         self.assertEqual(output_strategy, expected_strategy)
    #     with self.subTest():
    #         self.assertEqual(eval_result.status, "error")




class TestTicketFunctions3(unittest.TestCase):

    def setUp(self):

        self.ticket1 = ticket.GenomeTicket()
        self.ticket2 = ticket.GenomeTicket()
        self.ticket3 = ticket.GenomeTicket()
        self.ticket4 = ticket.GenomeTicket()

        self.bundle1 = bundle.Bundle()
        self.bundle2 = bundle.Bundle()
        self.bundle3 = bundle.Bundle()
        self.bundle4 = bundle.Bundle()

        self.bundle1.ticket = self.ticket1
        self.bundle2.ticket = self.ticket2
        self.bundle3.ticket = self.ticket3
        self.bundle4.ticket = self.ticket4


    def test_create_bundle_dict_1(self):
        """Verify dictionary of update tickets is created."""

        self.bundle1.ticket.type = "update"
        self.bundle2.ticket.type = "update"
        self.bundle3.ticket.type = "update"
        self.bundle4.ticket.type = "update"
        list = [self.bundle1, self.bundle2,
                self.bundle3, self.bundle4]

        result_dict = tickets.create_bundle_dict(list)
        update_list = result_dict["update"]

        with self.subTest():
            self.assertEqual(len(result_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(len(update_list), 4)

    def test_create_bundle_dict_2(self):
        """Verify dictionary of empty tickets is created."""

        list = [self.bundle1, self.bundle2,
                self.bundle3, self.bundle4]

        result_dict = tickets.create_bundle_dict(list)
        update_list = result_dict[""]

        with self.subTest():
            self.assertEqual(len(result_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(len(update_list), 4)


    def test_create_bundle_dict_3(self):
        """Verify dictionary of multiple ticket types is created."""

        self.bundle1.ticket.type = "update"
        self.bundle2.ticket.type = "add"
        self.bundle3.ticket.type = "update"
        self.bundle4.ticket.type = "replace"

        list = [self.bundle1, self.bundle2,
                self.bundle3, self.bundle4]

        result_dict = tickets.create_bundle_dict(list)
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
        self.add_ticket = ticket.GenomeTicket()
        self.add_ticket.type = "add"
        self.add_ticket.phage_id = "Trixie_Draft"
        self.add_ticket.run_mode = "phagesdb"
        self.add_ticket.description_field = "product"
        self.add_ticket.host_genus = "Mycobacterium smegmatis"
        self.add_ticket.cluster = "A"
        self.add_ticket.subcluster = "A2"
        self.add_ticket.annotation_status = "final"
        self.add_ticket.annotation_author = 1
        self.add_ticket.annotation_qc = 1
        self.add_ticket.retrieve_record = 1
        self.add_ticket.accession = "ABC123.1"
        self.bundle1 = bundle.Bundle()
        self.bundle1.ticket = self.add_ticket

        # Remove ticket.
        self.remove_ticket = ticket.GenomeTicket()
        self.remove_ticket.type = "replace"
        self.remove_ticket.phage_id = "Trixie_Draft"
        self.remove_ticket.run_mode = "phagesdb"
        self.remove_ticket.description_field = "product"
        self.remove_ticket.host_genus = "Mycobacterium smegmatis"
        self.remove_ticket.cluster = "A"
        self.remove_ticket.subcluster = "A2"
        self.remove_ticket.annotation_status = "final"
        self.remove_ticket.annotation_author = 1
        self.remove_ticket.annotation_qc = 1
        self.remove_ticket.retrieve_record = 1
        self.remove_ticket.accession = "ABC123.1"
        self.bundle2 = bundle.Bundle()
        self.bundle2.ticket = self.remove_ticket




        # Invalid ticket.
        self.invalid_ticket = ticket.GenomeTicket()
        self.invalid_ticket.type = "invalid"
        self.invalid_ticket.phage_id = "Trixie_Draft"
        self.invalid_ticket.run_mode = "phagesdb"
        self.invalid_ticket.description_field = "product"
        self.invalid_ticket.host_genus = "Mycobacterium smegmatis"
        self.invalid_ticket.cluster = "A"
        self.invalid_ticket.subcluster = "A2"
        self.invalid_ticket.annotation_status = "final"
        self.invalid_ticket.annotation_author = 1
        self.invalid_ticket.annotation_qc = 1
        self.invalid_ticket.retrieve_record = 1
        self.invalid_ticket.accession = "ABC123.1"
        self.bundle3 = bundle.Bundle()
        self.bundle3.ticket = self.invalid_ticket


    def test_copy_ticket_to_genome_1(self):
        """Verify data from 'add' ticket is added to genome."""

        tickets.copy_ticket_to_genome(self.bundle1)
        matched_genome = self.bundle1.genome_dict["add"]
        with self.subTest():
            self.assertEqual(matched_genome.id, "Trixie")
        with self.subTest():
            self.assertEqual(matched_genome.name, "Trixie_Draft")
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

    # TODO this is probably no longer needed.
    # def test_copy_ticket_to_genome_2(self):
    #     """Verify data from 'remove' ticket is added to genome."""
    #
    #     tickets.copy_ticket_to_genome(self.bundle2)
    #     matched_genome1 = self.bundle2.genome_dict["add"]
    #     # matched_genome2 = self.bundle2.genome_dict["remove"]
    #
    #     with self.subTest():
    #         self.assertEqual(matched_genome1.id, "Trixie")
    #     with self.subTest():
    #         self.assertEqual(matched_genome1.name, "Trixie_Draft")
    #     # with self.subTest():
    #     #     self.assertEqual(matched_genome2.id, "L5")
    #     # with self.subTest():
    #     #     self.assertEqual(matched_genome2.name, "")

    def test_copy_ticket_to_genome_3(self):
        """Verify data from 'invalid' ticket is not added to genome."""

        tickets.copy_ticket_to_genome(self.bundle3)
        self.assertEqual(len(self.bundle3.genome_dict.keys()), 0)







# TODO below code is broken since I have changed how eval is structured.
#
#
#
#
# class TestTicketFunctions5(unittest.TestCase):
#
#     def setUp(self):
#
#         self.genome1 = genome.Genome()
#         self.genome2 = genome.Genome()
#         self.genome3 = genome.Genome()
#         self.genome4 = genome.Genome()
#
#         self.ticket1 = ticket.GenomeTicket()
#         self.ticket2 = ticket.GenomeTicket()
#         self.ticket3 = ticket.GenomeTicket()
#         self.ticket4 = ticket.GenomeTicket()
#
#         self.bundle1 = bundle.Bundle()
#         self.bundle2 = bundle.Bundle()
#         self.bundle3 = bundle.Bundle()
#         self.bundle4 = bundle.Bundle()
#
#         self.bundle1.ticket = self.ticket1
#         self.bundle2.ticket = self.ticket2
#         self.bundle3.ticket = self.ticket3
#         self.bundle4.ticket = self.ticket4
#
#
#     def test_match_genomes_to_tickets_1(self):
#         """Verify that one genome is matched to ticket using phage_id."""
#
#         self.bundle1.ticket.phage_id = "Trixie"
#         self.bundle1.ticket.match_strategy = "phage_id"
#         self.genome1.id = "Trixie"
#
#
#         list1 = [self.bundle1] # Trixie
#         list2 = [self.genome1] # Trixie
#         eval_list = \
#             tickets.match_genomes_to_tickets(list1,
#                                                     list2,
#                                                     "phamerator")
#
#         matched_genome = list1[0].genome_dict["phamerator"]
#         id = matched_genome.id
#         expected_id = "Trixie"
#         with self.subTest():
#             self.assertEqual(id, expected_id)
#         with self.subTest():
#             self.assertEqual(len(eval_list), 0)
#
#     def test_match_genomes_to_tickets_2(self):
#         """Verify that one genome is matched to ticket using filename."""
#
#         self.bundle1.ticket.phage_id = "Trixie"
#         self.bundle1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#
#         list1 = [self.bundle1] # Trixie
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
#         self.bundle1.ticket.phage_id = "Trixie"
#         self.bundle1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#         self.genome2.filename = "L5"
#
#
#         list1 = [self.bundle1] # Trixie
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
#         self.bundle1.ticket.phage_id = "Trixie"
#         self.bundle1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#
#         self.bundle2.ticket.phage_id = "L5"
#         self.bundle2.ticket.match_strategy = "filename"
#
#         list1 = [self.bundle1, self.bundle2] # Trixie, L5
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
#         self.bundle1.ticket.phage_id = "Trixie"
#         self.bundle1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#
#         self.bundle2.ticket.phage_id = "L5"
#         self.bundle2.ticket.match_strategy = "filename"
#         self.genome2.filename = "L5"
#
#
#         list1 = [self.bundle1, self.bundle2] # Trixie, L5
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
#         self.bundle1.ticket.phage_id = "Trixie"
#         self.bundle1.ticket.match_strategy = "phage_id"
#         self.genome1.filename = "Trixie"
#
#         self.bundle2.ticket.phage_id = "L5"
#         self.bundle2.ticket.match_strategy = "filename"
#         self.genome2.filename = "L5"
#
#         list1 = [self.bundle1, self.bundle2] # Trixie, L5
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
#         self.bundle1.ticket.phage_id = "Trixie"
#         self.bundle1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#
#         self.bundle2.ticket.phage_id = "L5"
#         self.bundle2.ticket.match_strategy = "filename"
#         self.genome2.filename = "L5"
#
#         self.bundle3.ticket.phage_id = "D29"
#         self.bundle3.ticket.match_strategy = "filename"
#         self.genome3.filename = "D29"
#
#         list1 = [self.bundle1, self.bundle2, self.bundle3]
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
#         self.bundle1.ticket.phage_id = "Trixie"
#         self.bundle1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#
#         self.bundle2.ticket.phage_id = "L5"
#         self.bundle2.ticket.match_strategy = "filename"
#         self.genome2.filename = "L5"
#
#         self.bundle3.ticket.phage_id = "D29"
#         self.bundle3.ticket.match_strategy = "filename"
#
#         self.bundle4.ticket.phage_id = "D29"
#         self.bundle4.ticket.match_strategy = "filename"
#
#         list1 = [self.bundle1, self.bundle2,
#                 self.bundle3, self.bundle4]
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
#         self.bundle1.ticket.phage_id = "Trixie"
#         self.bundle1.ticket.match_strategy = "filename"
#         self.genome1.filename = "Trixie"
#
#         self.bundle2.ticket.phage_id = "L5"
#         self.bundle2.ticket.match_strategy = "filename"
#         self.genome2.filename = "L5"
#
#         self.genome3.filename = "D29"
#         self.genome4.filename = "D29"
#
#         list1 = [self.bundle1, self.bundle2]
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
