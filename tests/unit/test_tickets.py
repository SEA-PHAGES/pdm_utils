"""Unit tests for misc. ticket functions."""

import unittest

from pdm_utils.classes import bundle
from pdm_utils.classes import genome
from pdm_utils.classes import ticket
from pdm_utils.functions import tickets
from pdm_utils.constants import constants


class TestTicketFunctions1(unittest.TestCase):


    def setUp(self):
        self.required_keys = constants.IMPORT_TABLE_STRUCTURE["required"]
        self.optional_keys = constants.IMPORT_TABLE_STRUCTURE["optional"]
        self.keywords = constants.IMPORT_TABLE_STRUCTURE["keywords"]

        self.ticket_dict1 = {}
        self.ticket_dict1["type"] = "add"
        self.ticket_dict1["phage_id"] = "Trixie"
        self.ticket_dict1["description_field"] = "product"
        self.ticket_dict1["eval_mode"] = "final"
        self.ticket_dict1["host_genus"] = "retrieve"
        self.ticket_dict1["cluster"] = "retain"
        self.ticket_dict1["subcluster"] = "A2"
        self.ticket_dict1["accession"] = "parse"


        self.ticket_dict2 = {}

        self.ticket_dict3 = {}
        self.ticket_dict3["type"] = "ADD"
        self.ticket_dict3["phage_id"] = "Trixie"
        self.ticket_dict3["description_field"] = "PRODUCT"
        self.ticket_dict3["eval_mode"] = "FINAL"
        self.ticket_dict3["host_genus"] = "RETRIEVE"
        self.ticket_dict3["subcluster"] = None
        self.ticket_dict3["accession"] = "PARSE"
        self.ticket_dict3["retrieve_record"] = "RETAIN"


        self.ticket_dict4 = {}
        self.ticket_dict4["type"] = "ADD"
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
        """Verify returns True with keywords identified and values lowercased."""
        result = tickets.modify_import_data(self.ticket_dict3,
                    self.required_keys, self.optional_keys, self.keywords)
        with self.subTest():
            self.assertTrue(result)
        with self.subTest():
            self.assertEqual(self.ticket_dict3["host_genus"], "retrieve")
        with self.subTest():
            self.assertEqual(self.ticket_dict3["retrieve_record"], "retain")
        with self.subTest():
            self.assertEqual(self.ticket_dict3["subcluster"], "retrieve")
        with self.subTest():
            self.assertEqual(self.ticket_dict3["accession"], "parse")
        with self.subTest():
            self.assertEqual(self.ticket_dict3["type"], "add")
        with self.subTest():
            self.assertEqual(self.ticket_dict3["description_field"], "product")
        with self.subTest():
            self.assertEqual(self.ticket_dict3["eval_mode"], "final")


    def test_modify_import_data_4(self):
        """Verify returns True with completed dictionary from a
        minimal add ticket."""
        self.ticket_dict4["description_field"] = "product"
        self.ticket_dict4["eval_mode"] = "final"
        result = tickets.modify_import_data(self.ticket_dict4,
                    self.required_keys, self.optional_keys, self.keywords)
        with self.subTest():
            self.assertTrue(result)
        with self.subTest():
            self.assertEqual(self.ticket_dict4["host_genus"], "retrieve")
        with self.subTest():
            self.assertEqual(self.ticket_dict4["cluster"], "retrieve")
        with self.subTest():
            self.assertEqual(self.ticket_dict4["subcluster"], "retrieve")
        with self.subTest():
            self.assertEqual(self.ticket_dict4["annotation_author"], "1")
        with self.subTest():
            self.assertEqual(self.ticket_dict4["retrieve_record"], "1")
        with self.subTest():
            self.assertEqual(self.ticket_dict4["annotation_status"], "draft")
        with self.subTest():
            self.assertEqual(self.ticket_dict4["accession"], "")


    def test_modify_import_data_5(self):
        """Verify returns True with completed dictionary from a
        minimal replace ticket."""
        self.ticket_dict4["type"] = "replace"
        self.ticket_dict4["description_field"] = "product"
        self.ticket_dict4["eval_mode"] = "final"
        result = tickets.modify_import_data(self.ticket_dict4,
                    self.required_keys, self.optional_keys, self.keywords)
        with self.subTest():
            self.assertTrue(result)
        with self.subTest():
            self.assertEqual(self.ticket_dict4["host_genus"], "retain")
        with self.subTest():
            self.assertEqual(self.ticket_dict4["cluster"], "retain")
        with self.subTest():
            self.assertEqual(self.ticket_dict4["subcluster"], "retain")
        with self.subTest():
            self.assertEqual(self.ticket_dict4["annotation_author"], "retain")
        with self.subTest():
            self.assertEqual(self.ticket_dict4["retrieve_record"], "retain")
        with self.subTest():
            self.assertEqual(self.ticket_dict4["annotation_status"], "final")
        with self.subTest():
            self.assertEqual(self.ticket_dict4["accession"], "retain")




    def test_parse_import_ticket_data_1(self):
        """Verify ticket is generated from correct data dictionary."""
        tkt = tickets.parse_import_ticket_data(self.ticket_dict1)
        with self.subTest():
            self.assertEqual(tkt.type, "add")
        with self.subTest():
            self.assertEqual(tkt.phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(tkt.description_field, "product")
        with self.subTest():
            self.assertEqual(tkt.eval_mode, "final")
        with self.subTest():
            self.assertEqual(len(tkt.data_dict.keys()), 8)
        with self.subTest():
            self.assertEqual(tkt.data_retrieve, set(["host_genus"]))
        with self.subTest():
            self.assertEqual(tkt.data_retain, set(["cluster"]))
        with self.subTest():
            self.assertEqual(tkt.data_parse, set(["accession"]))
        with self.subTest():
            self.assertEqual(tkt.data_add, set(["subcluster"]))

    def test_parse_import_ticket_data_2(self):
        """Verify ticket is generated from correct data dictionary with
        no data in 'retain', 'retrieve', or 'parse' sets."""
        self.ticket_dict1["host_genus"] = "Mycobacterium"
        self.ticket_dict1["cluster"] = "A"
        self.ticket_dict1["subcluster"] = "A2"
        self.ticket_dict1["accession"] = "ABC123"
        tkt = tickets.parse_import_ticket_data(self.ticket_dict1)
        with self.subTest():
            self.assertEqual(tkt.type, "add")
        with self.subTest():
            self.assertEqual(tkt.phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(tkt.description_field, "product")
        with self.subTest():
            self.assertEqual(tkt.eval_mode, "final")
        with self.subTest():
            self.assertEqual(len(tkt.data_dict.keys()), 8)
        with self.subTest():
            self.assertEqual(tkt.data_retrieve, set())
        with self.subTest():
            self.assertEqual(tkt.data_retain, set())
        with self.subTest():
            self.assertEqual(tkt.data_parse, set())
        with self.subTest():
            self.assertEqual(tkt.data_add, set(["subcluster", "host_genus",
                                                "cluster", "accession"]))

    def test_parse_import_ticket_data_3(self):
        """Verify ticket is generated from correct data dictionary with
        no data in 'add' sets."""
        self.ticket_dict1["host_genus"] = "retrieve"
        self.ticket_dict1["cluster"] = "retrieve"
        self.ticket_dict1["subcluster"] = "retrieve"
        self.ticket_dict1["accession"] = "retrieve"
        tkt = tickets.parse_import_ticket_data(self.ticket_dict1)
        with self.subTest():
            self.assertEqual(tkt.type, "add")
        with self.subTest():
            self.assertEqual(tkt.phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(tkt.description_field, "product")
        with self.subTest():
            self.assertEqual(tkt.eval_mode, "final")
        with self.subTest():
            self.assertEqual(len(tkt.data_dict.keys()), 8)
        with self.subTest():
            self.assertEqual(tkt.data_retrieve, set(["subcluster", "host_genus",
                                                "cluster", "accession"]))
        with self.subTest():
            self.assertEqual(tkt.data_retain, set())
        with self.subTest():
            self.assertEqual(tkt.data_parse, set())
        with self.subTest():
            self.assertEqual(tkt.data_add, set())




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
        data_dict = {"type":"ADD",
                     "cluster":"RETRIEVE",
                     "subcluster": "NONE",
                     "host_genus": "PARSE",
                     "retrieve_record": "RETAIN"}
        keywords = set(["retrieve", "retain"])
        tickets.set_keywords(data_dict, self.keywords)
        with self.subTest():
            self.assertEqual(data_dict["type"], "ADD")
        with self.subTest():
            self.assertEqual(data_dict["cluster"], "retrieve")
        with self.subTest():
            self.assertEqual(data_dict["subcluster"], "none")
        with self.subTest():
            self.assertEqual(data_dict["host_genus"], "parse")
        with self.subTest():
            self.assertEqual(data_dict["retrieve_record"], "retain")




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
        """Verify two tickets are constructed correctly.
        The first ticket contains all required and optional fields.
        The second ticket contains all required fields."""
        dict_list = [self.ticket_dict1, self.ticket_dict4]
        eval_data_dict = {"eval_mode": "custom_eval_mode",
                          "eval_flag_dict": {"check_locus_tag": False}}
        list_of_tickets = tickets.construct_tickets(dict_list,
                eval_data_dict, "function", self.required_keys,
                self.optional_keys, self.keywords)
        with self.subTest():
            self.assertEqual(len(list_of_tickets), 2)
        with self.subTest():
            self.assertEqual(list_of_tickets[0].id, 1)
        with self.subTest():
            self.assertEqual(list_of_tickets[0].eval_mode, "final")
        with self.subTest():
            self.assertEqual(list_of_tickets[0].description_field, "product")
        with self.subTest():
            self.assertTrue(list_of_tickets[0].eval_flags["check_locus_tag"])
        with self.subTest():
            self.assertEqual(list_of_tickets[1].id, 2)
        with self.subTest():
            self.assertEqual(list_of_tickets[1].eval_mode, "custom_eval_mode")
        with self.subTest():
            self.assertEqual(list_of_tickets[1].description_field, "function")
        with self.subTest():
            self.assertFalse(list_of_tickets[1].eval_flags["check_locus_tag"])

    def test_construct_tickets_2(self):
        """Verify one ticket is constructed correctly. The second data
        dictionary is not structured correctly."""
        dict_list = [self.ticket_dict1, self.ticket_dict2]
        eval_data_dict = {"eval_mode": "custom_eval_mode",
                          "eval_flag_dict": {}}
        list_of_tickets = tickets.construct_tickets(dict_list,
                eval_data_dict, "function", self.required_keys,
                self.optional_keys, self.keywords)
        with self.subTest():
            self.assertEqual(len(list_of_tickets), 1)

    def test_construct_tickets_3(self):
        """Verify four tickets constructed correctly. The first two tickets
        contain all required and optional fields. The second two tickets
        contain all required fields. Verify that each eval_flag dictionary
        is a separate object that can be modified without impacting the other
        eval_flag dictionaries."""

        tkt_dict1 = {}
        tkt_dict1["type"] = "add"
        tkt_dict1["phage_id"] = "Trixie"
        tkt_dict1["description_field"] = "product"
        tkt_dict1["eval_mode"] = "final"

        tkt_dict2 = {}
        tkt_dict2["type"] = "add"
        tkt_dict2["phage_id"] = "L5"
        tkt_dict2["description_field"] = "product"
        tkt_dict2["eval_mode"] = "final"

        tkt_dict3 = {}
        tkt_dict3["type"] = "add"
        tkt_dict3["phage_id"] = "RedRock"

        tkt_dict4 = {}
        tkt_dict4["type"] = "add"
        tkt_dict4["phage_id"] = "Bxb1"

        dict_list = [tkt_dict1, tkt_dict2, tkt_dict3, tkt_dict4]
        eval_data_dict = {"eval_mode": "custom_eval_mode",
                          "eval_flag_dict": {"check_locus_tag": False}}
        tkt_list = tickets.construct_tickets(dict_list,
                eval_data_dict, "function", self.required_keys,
                self.optional_keys, self.keywords)

        tkt_list[0].eval_flags["check_locus_tag"] = 0
        tkt_list[1].eval_flags["check_locus_tag"] = 1
        tkt_list[2].eval_flags["check_locus_tag"] = 2
        tkt_list[3].eval_flags["check_locus_tag"] = 3

        with self.subTest():
            self.assertEqual(tkt_list[0].eval_flags["check_locus_tag"], 0)
        with self.subTest():
            self.assertEqual(tkt_list[1].eval_flags["check_locus_tag"], 1)
        with self.subTest():
            self.assertEqual(tkt_list[2].eval_flags["check_locus_tag"], 2)
        with self.subTest():
            self.assertEqual(tkt_list[3].eval_flags["check_locus_tag"], 3)



    def test_identify_duplicates_1(self):
        """Verify no duplicates are produced."""

        ticket1 = ticket.ImportTicket()
        ticket1.id = 1
        ticket1.type = "replace"
        ticket1.phage_id = "Trixie"

        ticket2 = ticket.ImportTicket()
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

        ticket1 = ticket.ImportTicket()
        ticket1.id = "none"
        ticket1.type = "replace"
        ticket1.phage_id = "none"

        ticket2 = ticket.ImportTicket()
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

        ticket1 = ticket.ImportTicket()
        ticket1.id = 1
        ticket1.type = "replace"
        ticket1.phage_id = "L5"

        ticket2 = ticket.ImportTicket()
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

        ticket1 = ticket.ImportTicket()
        ticket1.id = 1
        ticket1.type = "replace"
        ticket1.phage_id = "Trixie"

        ticket2 = ticket.ImportTicket()
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

        ticket1 = ticket.ImportTicket()
        ticket1.id = 1
        ticket1.type = "replace"
        ticket1.phage_id = "Trixie"

        ticket2 = ticket.ImportTicket()
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

        self.ticket1 = ticket.ImportTicket()
        self.ticket2 = ticket.ImportTicket()

        self.ticket1.phage_id = "Trixie"
        self.ticket2.phage_id = "L5"

        self.bundle1 = bundle.Bundle()
        self.bundle2 = bundle.Bundle()

        self.bundle1.ticket = self.ticket1
        self.bundle2.ticket = self.ticket2




class TestTicketFunctions3(unittest.TestCase):

    def setUp(self):
        self.data_dict = {}
        self.data_dict["host_genus"] = "Mycobacterium smegmatis"
        self.data_dict["accession"] = "ABC123.1"
        self.data_dict["annotation_status"] = "final"
        self.data_dict["cluster"] = "A"
        self.data_dict["subcluster"] = "A2"
        self.data_dict["annotation_author"] = 1
        self.data_dict["retrieve_record"] = 1
        self.tkt1 = ticket.ImportTicket()
        self.tkt1.phage_id = "Trixie_Draft"
        self.tkt1.data_dict = self.data_dict

    def test_get_genome_1(self):
        """Verify no data from ticket is added to genome."""
        self.tkt1.data_add = set([""])
        gnm = tickets.get_genome(self.tkt1, gnm_type="add")
        with self.subTest():
            self.assertEqual(gnm.id, "Trixie")
        with self.subTest():
            self.assertEqual(gnm.name, "Trixie_Draft")
        with self.subTest():
            self.assertEqual(gnm.type, "add")
        with self.subTest():
            self.assertEqual(gnm.host_genus, "")
        with self.subTest():
            self.assertEqual(gnm.cluster, "")
        with self.subTest():
            self.assertEqual(gnm.subcluster, "")
        with self.subTest():
            self.assertEqual(gnm.annotation_status, "")
        with self.subTest():
            self.assertEqual(gnm.annotation_author, -1)
        with self.subTest():
            self.assertEqual(gnm.retrieve_record, -1)
        with self.subTest():
            self.assertEqual(gnm.accession, "")

    def test_get_genome_2(self):
        """Verify host_genus data from ticket is added to genome."""
        self.tkt1.data_add = set(["host_genus"])
        gnm = tickets.get_genome(self.tkt1, gnm_type="add")
        with self.subTest():
            self.assertEqual(gnm.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(gnm.cluster, "")

    def test_get_genome_3(self):
        """Verify cluster data from ticket is added to genome."""
        self.tkt1.data_add = set(["cluster"])
        gnm = tickets.get_genome(self.tkt1, gnm_type="add")
        with self.subTest():
            self.assertEqual(gnm.host_genus, "")
        with self.subTest():
            self.assertEqual(gnm.cluster, "A")

    def test_get_genome_4(self):
        """Verify subcluster data from ticket is added to genome."""
        self.tkt1.data_add = set(["subcluster"])
        gnm = tickets.get_genome(self.tkt1, gnm_type="add")
        with self.subTest():
            self.assertEqual(gnm.host_genus, "")
        with self.subTest():
            self.assertEqual(gnm.subcluster, "A2")

    def test_get_genome_5(self):
        """Verify annotation_status data from ticket is added to genome."""
        self.tkt1.data_add = set(["annotation_status"])
        gnm = tickets.get_genome(self.tkt1, gnm_type="add")
        with self.subTest():
            self.assertEqual(gnm.host_genus, "")
        with self.subTest():
            self.assertEqual(gnm.annotation_status, "final")

    def test_get_genome_6(self):
        """Verify annotation_author data from ticket is added to genome."""
        self.tkt1.data_add = set(["annotation_author"])
        gnm = tickets.get_genome(self.tkt1, gnm_type="add")
        with self.subTest():
            self.assertEqual(gnm.host_genus, "")
        with self.subTest():
            self.assertEqual(gnm.annotation_author, 1)

    def test_get_genome_7(self):
        """Verify retrieve_record data from ticket is added to genome."""
        self.tkt1.data_add = set(["retrieve_record"])
        gnm = tickets.get_genome(self.tkt1, gnm_type="add")
        with self.subTest():
            self.assertEqual(gnm.host_genus, "")
        with self.subTest():
            self.assertEqual(gnm.retrieve_record, 1)

    def test_get_genome_8(self):
        """Verify accession data from ticket is added to genome."""
        self.tkt1.data_add = set(["accession"])
        gnm = tickets.get_genome(self.tkt1, gnm_type="add")
        with self.subTest():
            self.assertEqual(gnm.host_genus, "")
        with self.subTest():
            self.assertEqual(gnm.accession, "ABC123")

if __name__ == '__main__':
    unittest.main()
