"""Unit tests for misc. ticket functions."""

from pdm_utils.classes import bundle
from pdm_utils.classes import genome
from pdm_utils.classes import ticket
from pdm_utils.classes import eval
from pdm_utils.functions import tickets
from pdm_utils.constants import constants
import unittest





class TestTicketFunctions1(unittest.TestCase):


    def setUp(self):
        self.required_keys = constants.IMPORT_TABLE_REQ_DICT.keys()
        self.optional_keys = constants.IMPORT_TABLE_OPT_DICT.keys()
        self.keywords = set(["retrieve", "retain", "none"])

        self.ticket_dict1 = {}
        self.ticket_dict1["type"] = "add"
        self.ticket_dict1["id"] = 1
        self.ticket_dict1["phage_id"] = "Trixie"
        self.ticket_dict1["description_field"] = "product"
        self.ticket_dict1["run_mode"] = "phagesdb"
        self.ticket_dict1["host_genus"] = "retrieve"
        self.ticket_dict1["cluster"] = "retain"

        self.ticket_dict2 = {}

        self.ticket_dict3 = {}
        self.ticket_dict3["type"] = "ADD"
        self.ticket_dict3["id"] = 1
        self.ticket_dict3["phage_id"] = "Trixie"
        self.ticket_dict3["description_field"] = "PRODUCT"
        self.ticket_dict3["run_mode"] = "PHAGESDB"
        self.ticket_dict3["host_genus"] = "RETRIEVE"
        self.ticket_dict3["subcluster"] = None


        self.ticket_dict4 = {}
        self.ticket_dict4["type"] = "ADD"
        self.ticket_dict4["id"] = 1
        self.ticket_dict4["phage_id"] = "Trixie"


    def test_modify_import_data_1(self):
        """Verify returns False if there are missing required keys."""
        result = tickets.modify_import_data(self.ticket_dict2,
                    self.required_keys, self.optional_keys, self.keywords)
        self.assertFalse(result)


    def test_modify_import_data_2(self):
        """Verify returns False if there are extra keys."""
        self.ticket_dict3["extra"] = "extra"
        result = tickets.modify_import_data(self.ticket_dict3,
                    self.required_keys, self.optional_keys, self.keywords)
        self.assertFalse(result)


    def test_modify_import_data_3(self):
        """Verify returns True with completed dictionary."""
        result = tickets.modify_import_data(self.ticket_dict3,
                    self.required_keys, self.optional_keys, self.keywords)
        with self.subTest():
            self.assertTrue(result)
        with self.subTest():
            self.assertEqual(self.ticket_dict3["subcluster"], "retrieve")
        with self.subTest():
            self.assertEqual(self.ticket_dict3["host_genus"], "retrieve")
        with self.subTest():
            self.assertEqual(self.ticket_dict3["annotation_author"], "1")
        with self.subTest():
            self.assertEqual(self.ticket_dict3["type"], "add")
        with self.subTest():
            self.assertEqual(self.ticket_dict3["description_field"], "product")
        with self.subTest():
            self.assertEqual(self.ticket_dict3["run_mode"], "phagesdb")




    def test_parse_import_ticket_data_1(self):
        """Verify ticket is generated from correct data dictionary."""
        tkt = tickets.parse_import_ticket_data(self.ticket_dict1)
        with self.subTest():
            self.assertEqual(tkt.id, 1)
        with self.subTest():
            self.assertEqual(tkt.type, "add")
        with self.subTest():
            self.assertEqual(tkt.phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(tkt.description_field, "product")
        with self.subTest():
            self.assertEqual(tkt.run_mode, "phagesdb")
        with self.subTest():
            self.assertEqual(len(tkt.data_dict.keys()), 7)
        with self.subTest():
            self.assertEqual(len(tkt.data_retrieve), 1)
        with self.subTest():
            self.assertEqual(len(tkt.data_retain), 1)
        with self.subTest():
            self.assertEqual(len(tkt.data_ticket), 5)




    def test_set_empty_1(self):
        """Verify one None value is set to ''."""
        data_dict = {"type":"add","cluster":None}
        tickets.set_empty(data_dict)
        with self.subTest():
            self.assertEqual(data_dict["type"], "add")
        with self.subTest():
            self.assertEqual(data_dict["cluster"], "")




    def test_set_keywords_1(self):
        """Verify one value is lowercased."""
        data_dict = {"type":"ADD","cluster":"RETRIEVE"}
        keywords = set(["retrieve", "retain"])
        tickets.set_keywords(data_dict, keywords)
        with self.subTest():
            self.assertEqual(data_dict["type"], "ADD")
        with self.subTest():
            self.assertEqual(data_dict["cluster"], "retrieve")




    def test_set_missing_keys_1(self):
        """Verify one missing key is added."""
        data_dict = {"type":"add", "cluster":""}
        key_set = set(["type", "host_genus"])
        tickets.set_missing_keys(data_dict, key_set)
        with self.subTest():
            self.assertEqual(len(data_dict.keys()), 3)
        with self.subTest():
            self.assertEqual(data_dict["host_genus"], "")

    def test_set_missing_keys_2(self):
        """Verify no missing key is added."""
        data_dict = {"type":"add", "cluster":""}
        key_set = set(["type", "cluster"])
        tickets.set_missing_keys(data_dict, key_set)
        self.assertEqual(len(data_dict.keys()), 2)




    def test_set_dict_value_1(self):
        """Verify empty value is replaced with first value."""
        data_dict = {"type":"add", "cluster":""}
        tickets.set_dict_value(data_dict, "cluster", "A", "B")
        self.assertEqual(data_dict["cluster"], "A")

    def test_set_dict_value_2(self):
        """Verify empty value is replaced with second value."""
        data_dict = {"type":"replace", "cluster":""}
        tickets.set_dict_value(data_dict, "cluster", "A", "B")
        self.assertEqual(data_dict["cluster"], "B")

    def test_set_dict_value_3(self):
        """Verify non-empty value is not replaced."""
        data_dict = {"type":"replace", "cluster":"C"}
        tickets.set_dict_value(data_dict, "cluster", "A", "B")
        self.assertEqual(data_dict["cluster"], "C")



    def test_construct_tickets_1(self):
        """Verify two tickets constructed correctly. The first ticket
        contains all required and optional fields. The second ticket contains
        all required fields."""
        dict_list = [self.ticket_dict1, self.ticket_dict4]
        list_of_tickets = tickets.construct_tickets(dict_list,
                "pecaan", "function", self.required_keys,
                self.optional_keys, self.keywords)
        with self.subTest():
            self.assertEqual(len(list_of_tickets), 2)
        with self.subTest():
            self.assertEqual(list_of_tickets[0].run_mode, "phagesdb")
        with self.subTest():
            self.assertEqual(list_of_tickets[0].description_field, "product")
        with self.subTest():
            self.assertTrue(list_of_tickets[0].eval_flags["check_locus_tag"])
        with self.subTest():
            self.assertEqual(list_of_tickets[1].run_mode, "pecaan")
        with self.subTest():
            self.assertEqual(list_of_tickets[1].description_field, "function")
        with self.subTest():
            self.assertFalse(list_of_tickets[1].eval_flags["check_locus_tag"])

    def test_construct_tickets_2(self):
        """Verify one ticket is constructed correctly. The second data
        dictionary is not structured correctly."""
        dict_list = [self.ticket_dict1, self.ticket_dict2]
        list_of_tickets = tickets.construct_tickets(dict_list,
                "pecaan", "function", self.required_keys,
                self.optional_keys, self.keywords)
        with self.subTest():
            self.assertEqual(len(list_of_tickets), 1)




    def test_identify_duplicates_1(self):
        """Verify no duplicates are produced."""

        ticket1 = ticket.GenomeTicket()
        ticket1.id = 1
        ticket1.type = "replace"
        ticket1.phage_id = "Trixie"

        ticket2 = ticket.GenomeTicket()
        ticket2.id = 2
        ticket2.type = "replace"
        ticket2.phage_id = "L5"

        null_set = set(["none"])
        list_of_tickets = [ticket1, ticket2]
        id_dupes, phage_id_dupes = \
            tickets.identify_duplicates(list_of_tickets, null_set=null_set)

        with self.subTest():
            self.assertEqual(len(id_dupes), 0)
        with self.subTest():
            self.assertEqual(len(phage_id_dupes), 0)


    def test_identify_duplicates_2(self):
        """Verify two tickets with 'none' duplicates
        do not generate an error."""

        ticket1 = ticket.GenomeTicket()
        ticket1.id = "none"
        ticket1.type = "replace"
        ticket1.phage_id = "none"

        ticket2 = ticket.GenomeTicket()
        ticket2.id = "none"
        ticket2.type = "replace"
        ticket2.phage_id = "none"

        null_set = set(["none"])
        list_of_tickets = [ticket1, ticket2]
        id_dupes, phage_id_dupes = \
            tickets.identify_duplicates(list_of_tickets, null_set=null_set)
        with self.subTest():
            self.assertEqual(len(id_dupes), 0)
        with self.subTest():
            self.assertEqual(len(phage_id_dupes), 0)


    def test_identify_duplicates_3(self):
        """Verify two tickets with id duplicates
        do generate an error."""

        ticket1 = ticket.GenomeTicket()
        ticket1.id = 1
        ticket1.type = "replace"
        ticket1.phage_id = "L5"

        ticket2 = ticket.GenomeTicket()
        ticket2.id = 1
        ticket2.type = "replace"
        ticket2.phage_id = "Trixie"

        null_set = set(["none"])
        list_of_tickets = [ticket1, ticket2]
        id_dupes, phage_id_dupes = \
            tickets.identify_duplicates(list_of_tickets, null_set=null_set)
        with self.subTest():
            self.assertEqual(len(id_dupes), 1)
        with self.subTest():
            self.assertEqual(len(phage_id_dupes), 0)



    def test_identify_duplicates_4(self):
        """Verify two tickets with Primary Phage ID duplicates
        do generate an error."""

        ticket1 = ticket.GenomeTicket()
        ticket1.id = 1
        ticket1.type = "replace"
        ticket1.phage_id = "Trixie"

        ticket2 = ticket.GenomeTicket()
        ticket2.id = 2
        ticket2.type = "replace"
        ticket2.phage_id = "Trixie"

        null_set = set(["none"])
        list_of_tickets = [ticket1, ticket2]
        id_dupes, phage_id_dupes = \
            tickets.identify_duplicates(list_of_tickets, null_set=null_set)
        with self.subTest():
            self.assertEqual(len(id_dupes), 0)
        with self.subTest():
            self.assertEqual(len(phage_id_dupes), 1)


    def test_identify_duplicates_6(self):
        """Verify two tickets with multiple duplicates
        do generate multiple errors."""

        ticket1 = ticket.GenomeTicket()
        ticket1.id = 1
        ticket1.type = "replace"
        ticket1.phage_id = "Trixie"

        ticket2 = ticket.GenomeTicket()
        ticket2.id = 1
        ticket2.type = "replace"
        ticket2.phage_id = "Trixie"

        null_set = set(["none"])
        list_of_tickets = [ticket1, ticket2]
        id_dupes, phage_id_dupes = \
            tickets.identify_duplicates(list_of_tickets, null_set=null_set)
        with self.subTest():
            self.assertEqual(len(id_dupes), 1)
        with self.subTest():
            self.assertEqual(len(phage_id_dupes), 1)



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
            self.assertEqual(matched_genome.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(matched_genome.accession, "ABC123")


    def test_copy_ticket_to_genome_3(self):
        """Verify data from 'invalid' ticket is not added to genome."""

        tickets.copy_ticket_to_genome(self.bundle3)
        self.assertEqual(len(self.bundle3.genome_dict.keys()), 0)





if __name__ == '__main__':
    unittest.main()
