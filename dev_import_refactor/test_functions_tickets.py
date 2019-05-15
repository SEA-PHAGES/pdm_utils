""" Unit tests for misc. ticket functions"""


import Ticket
import Eval
import functions_tickets
import unittest





class TestTicketFunctionsClass(unittest.TestCase):


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

        # self.list_of_tickets = [add_ticket_list, remove_ticket_list]





    def test_parse_import_ticket_1(self):
        """Verify properly structured data is parsed correctly."""
        ticket, eval = \
            functions_tickets.parse_import_ticket(self.normal_ticket_list)
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
            functions_tickets.parse_import_ticket(self.short_ticket_list)
        with self.subTest():
            self.assertIsNotNone(eval)
        with self.subTest():
            self.assertIsNone(ticket)
        with self.subTest():
            self.assertEqual(eval.status, "error")

    def test_parse_import_ticket_3(self):
        """Verify improperly structured data is not parsed."""
        ticket, eval = \
            functions_tickets.parse_import_ticket(self.long_ticket_list)
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
            functions_tickets.parse_import_tickets(list_of_lists)
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
            functions_tickets.parse_import_tickets(list_of_lists)
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
            functions_tickets.parse_import_tickets(list_of_lists)
        with self.subTest():
            self.assertEqual(len(list_of_evals), 2)
        with self.subTest():
            self.assertEqual(len(list_of_tickets), 1)







    def test_identify_duplicates1_1(self):
        """Verify duplicate items generates an error."""
        eval_object = \
            functions_tickets.identify_duplicates1(["Trixie", "Trixie"], "temp")
        with self.subTest():
            self.assertIsNotNone(eval_object)
        with self.subTest():
            self.assertIsInstance(eval_object,Eval.EvalResult)

    def test_identify_duplicates1_2(self):
        """Verify non-duplicate items do not generate an error."""
        eval_object = \
            functions_tickets.identify_duplicates1(["Trixie", "L5"], "temp")
        with self.subTest():
            self.assertIsNone(eval_object)
        with self.subTest():
            self.assertNotIsInstance(eval_object,Eval.EvalResult)











    def test_identify_duplicates2_1(self):
        """Verify duplicate items in two lists generate an error."""
        list1 = ["Trixie", "L5"]
        list2 = ["Trixie", "RedRock"]
        eval_object = \
            functions_tickets.identify_duplicates2(list1, list2, "temp")
        with self.subTest():
            self.assertIsNotNone(eval_object)
        with self.subTest():
            self.assertIsInstance(eval_object,Eval.EvalResult)




    def test_identify_duplicates2_2(self):
        """Verify non-duplicate items in two lists do not generate an error."""
        list1 = ["Trixie", "L5"]
        list2 = ["D29", "RedRock"]
        eval_object = \
            functions_tickets.identify_duplicates2(list1, list2, "temp")
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
        ticket2.accession = "DEF456"

        list_of_tickets = [ticket1, ticket2]
        result_list = \
            functions_tickets.validate_tickets(list_of_tickets)

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
            functions_tickets.validate_tickets(list_of_tickets)

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
            functions_tickets.validate_tickets(list_of_tickets)

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
            functions_tickets.validate_tickets(list_of_tickets)

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
            functions_tickets.validate_tickets(list_of_tickets)

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
            functions_tickets.validate_tickets(list_of_tickets)

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
            functions_tickets.validate_tickets(list_of_tickets)

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
            functions_tickets.validate_tickets(list_of_tickets)

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
            functions_tickets.validate_tickets(list_of_tickets)

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
            functions_tickets.validate_tickets(list_of_tickets)

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
            functions_tickets.validate_tickets(list_of_tickets)
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
            functions_tickets.validate_tickets(list_of_tickets)

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
            functions_tickets.validate_tickets(list_of_tickets)

        self.assertEqual(len(result_list), 1)













if __name__ == '__main__':
    unittest.main()
