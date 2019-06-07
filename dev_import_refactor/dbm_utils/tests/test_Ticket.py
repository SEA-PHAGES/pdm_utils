""" Unit tests for the ticket class"""



from classes import Ticket
import unittest



class TestImportTicketClass(unittest.TestCase):


    def setUp(self):

        # Empty ticket to test simple methods
        self.ticket = Ticket.ImportTicket()

        # Standard update ticket
        self.update_ticket = Ticket.ImportTicket()
        self.update_ticket.type = "update"
        self.update_ticket.primary_phage_id = "Trixie"
        self.update_ticket.host = "Mycobacterium"
        self.update_ticket.cluster = "A"
        self.update_ticket.status = "Final"
        self.update_ticket.description_field = "none"
        self.update_ticket.secondary_phage_id = "none"
        self.update_ticket.annotation_author = "1"
        self.update_ticket.run_mode = "none"



        # Standard add ticket
        self.add_ticket = Ticket.ImportTicket()
        self.add_ticket.type = "add"
        self.add_ticket.primary_phage_id = "Trixie"
        self.add_ticket.host = "Mycobacterium"
        self.add_ticket.cluster = "A"
        self.add_ticket.status = "draft"
        self.add_ticket.description_field = "Product"
        self.add_ticket.secondary_phage_id = "none"
        self.add_ticket.annotation_author = "1"
        self.add_ticket.run_mode = "phagesdb"


        # Standard remove ticket
        self.remove_ticket = Ticket.ImportTicket()
        self.remove_ticket.type = "remove"
        self.remove_ticket.primary_phage_id = "none"
        self.remove_ticket.host = "none"
        self.remove_ticket.cluster = "none"
        self.remove_ticket.subcluster = "none"
        self.remove_ticket.status = "none"
        self.remove_ticket.description_field = "none"
        self.remove_ticket.accession = "none"
        self.remove_ticket.annotation_author = "none"
        self.remove_ticket.run_mode = "none"
        self.remove_ticket.secondary_phage_id = "Trixie"





        # Standard replace ticket
        self.replace_ticket = Ticket.ImportTicket()
        self.replace_ticket.type = "replace"
        self.replace_ticket.primary_phage_id = "Trixie"
        self.replace_ticket.host = "Mycobacterium"
        self.replace_ticket.cluster = "A"
        self.replace_ticket.status = "draft"
        self.replace_ticket.description_field = "Product"
        self.replace_ticket.secondary_phage_id = "none"
        self.replace_ticket.annotation_author = "1"
        self.replace_ticket.run_mode = "phagesdb"
        self.replace_ticket.secondary_phage_id = "Trixie"








    def test_set_case_1(self):
        self.ticket.type = "CDS"
        self.ticket.set_case()
        self.assertEqual(self.ticket.type, "cds")

    def test_set_case_2(self):
        self.ticket.primary_phage_id = "Trixie"
        self.ticket.set_case()
        self.assertEqual(self.ticket.primary_phage_id, "Trixie")

    def test_set_case_3(self):
        self.ticket.primary_phage_id = "None"
        self.ticket.set_case()
        self.assertEqual(self.ticket.primary_phage_id, "none")

    def test_set_case_4(self):
        self.ticket.host = "Mycobacterium"
        self.ticket.set_case()
        self.assertEqual(self.ticket.host, "Mycobacterium")

    def test_set_case_5(self):
        self.ticket.host = "None"
        self.ticket.set_case()
        self.assertEqual(self.ticket.host, "none")

    def test_set_case_6(self):
        self.ticket.host = "Retrieve"
        self.ticket.set_case()
        self.assertEqual(self.ticket.host, "retrieve")

    def test_set_case_7(self):
        self.ticket.cluster = "A"
        self.ticket.set_case()
        self.assertEqual(self.ticket.cluster, "A")

    def test_set_case_8(self):
        self.ticket.cluster = "None"
        self.ticket.set_case()
        self.assertEqual(self.ticket.cluster, "none")

    def test_set_case_9(self):
        self.ticket.cluster = "Retrieve"
        self.ticket.set_case()
        self.assertEqual(self.ticket.cluster, "retrieve")

    def test_set_case_10(self):
        self.ticket.subcluster = "A2"
        self.ticket.set_case()
        self.assertEqual(self.ticket.subcluster, "A2")

    def test_set_case_11(self):
        self.ticket.subcluster = "None"
        self.ticket.set_case()
        self.assertEqual(self.ticket.subcluster, "none")

    def test_set_case_12(self):
        self.ticket.subcluster = "Retrieve"
        self.ticket.set_case()
        self.assertEqual(self.ticket.subcluster, "retrieve")

    def test_set_case_13(self):
        self.ticket.status = "Final"
        self.ticket.set_case()
        self.assertEqual(self.ticket.status, "final")

    def test_set_case_14(self):
        self.ticket.description_field = "Product"
        self.ticket.set_case()
        self.assertEqual(self.ticket.description_field, "product")

    def test_set_case_15(self):
        self.ticket.accession = "ABCD1234"
        self.ticket.set_case()
        self.assertEqual(self.ticket.accession, "ABCD1234")

    def test_set_case_16(self):
        self.ticket.accession = "None"
        self.ticket.set_case()
        self.assertEqual(self.ticket.accession, "none")

    def test_set_case_17(self):
        self.ticket.accession = "Retrieve"
        self.ticket.set_case()
        self.assertEqual(self.ticket.accession, "retrieve")

    def test_set_case_18(self):
        self.ticket.annotation_author = "Hatfull"
        self.ticket.set_case()
        self.assertEqual(self.ticket.annotation_author, "hatfull")

    def test_set_case_19(self):
        self.ticket.secondary_phage_id = "Trixie"
        self.ticket.set_case()
        self.assertEqual(self.ticket.secondary_phage_id, "Trixie")

    def test_set_case_20(self):
        self.ticket.secondary_phage_id = "None"
        self.ticket.set_case()
        self.assertEqual(self.ticket.secondary_phage_id, "none")

    def test_set_case_21(self):
        self.ticket.run_mode = "Phagesdb"
        self.ticket.set_case()
        self.assertEqual(self.ticket.run_mode, "phagesdb")



    def test_set_evaluation_1(self):
        self.ticket.set_evaluation("none")
        self.assertEqual(len(self.ticket.evaluations), 1)

    def test_set_evaluation_2(self):
        self.ticket.set_evaluation("warning","message1")
        self.assertEqual(len(self.ticket.evaluations), 1)

    def test_set_evaluation_3(self):
        self.ticket.set_evaluation("error","message1","message2")
        self.assertEqual(len(self.ticket.evaluations), 1)






    def test_clear_retrieve_status_1(self):
        self.ticket.host = "retrieve"
        self.ticket.clear_retrieve_status()
        with self.subTest():
            self.assertEqual(self.ticket.host, "none")
        with self.subTest():
            self.assertEqual(len(self.ticket.evaluations), 1)
        with self.subTest():
            self.assertEqual(self.ticket.evaluations[0].status, "error")

    def test_clear_retrieve_status_2(self):
        self.ticket.cluster = "retrieve"
        self.ticket.clear_retrieve_status()
        with self.subTest():
            self.assertEqual(self.ticket.cluster, "none")
        with self.subTest():
            self.assertEqual(len(self.ticket.evaluations), 1)
        with self.subTest():
            self.assertEqual(self.ticket.evaluations[0].status, "error")

    def test_clear_retrieve_status_3(self):
        self.ticket.subcluster = "retrieve"
        self.ticket.clear_retrieve_status()
        with self.subTest():
            self.assertEqual(self.ticket.subcluster, "none")
        with self.subTest():
            self.assertEqual(len(self.ticket.evaluations), 1)
        with self.subTest():
            self.assertEqual(self.ticket.evaluations[0].status, "error")

    def test_clear_retrieve_status_4(self):
        self.ticket.accession = "retrieve"
        self.ticket.clear_retrieve_status()
        with self.subTest():
            self.assertEqual(self.ticket.accession, "none")
        with self.subTest():
            self.assertEqual(len(self.ticket.evaluations), 1)
        with self.subTest():
            self.assertEqual(self.ticket.evaluations[0].status, "error")






    def test_check_type_1(self):
        self.ticket.type = "add"
        self.ticket.check_type()
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_type_2(self):
        self.ticket.type = "remove"
        self.ticket.check_type()
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_type_3(self):
        self.ticket.type = "replace"
        self.ticket.check_type()
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_type_4(self):
        self.ticket.type = "update"
        self.ticket.check_type()
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_type_5(self):
        self.ticket.type = "abcd"
        self.ticket.check_type()
        self.assertEqual(len(self.ticket.evaluations), 1)









    def test_check_host_1(self):
        test_list = ["Mycobacterium","Gordonia"]
        self.ticket.host = "none"
        self.ticket.check_host(test_list)
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_host_2(self):
        test_list = ["Mycobacterium","Gordonia"]
        self.ticket.host = "Mycobacterium"
        self.ticket.check_host(test_list)
        with self.subTest():
            self.assertEqual(self.ticket.host, "Mycobacterium")
        with self.subTest():
            self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_host_3(self):
        test_list = ["Mycobacterium","Gordonia"]
        self.ticket.host = "Mycobacterium smegmatis"
        self.ticket.check_host(test_list)
        with self.subTest():
            self.assertEqual(self.ticket.host, "Mycobacterium")
        with self.subTest():
            self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_host_4(self):
        test_list = ["Mycobacterium","Gordonia"]
        msg1 = "The host strain Arthrobacter is not currently in the database."
        msg2 = "The host strain Arthrobacter is not correct."
        self.ticket.host = "Arthrobacter abcd"
        self.ticket.check_host(test_list)
        with self.subTest():
            self.assertEqual(self.ticket.host, "Arthrobacter")
        with self.subTest():
            self.assertEqual(len(self.ticket.evaluations), 1)
        with self.subTest():
            self.assertEqual(self.ticket.evaluations[0].messages["warning"], msg1)
        with self.subTest():
            self.assertEqual(self.ticket.evaluations[0].messages["error"], msg2)





    def test_check_subcluster_1(self):
        test_list = ["AA1", "AA2"]
        self.ticket.subcluster = "none"
        self.ticket.check_subcluster(test_list)
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_subcluster_2(self):
        test_list = ["AA1", "AA2"]
        self.ticket.subcluster = "AA1"
        self.ticket.check_subcluster(test_list)
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_subcluster_3(self):
        test_list = ["AA1", "AA2"]
        self.ticket.subcluster = "AA3"
        self.ticket.check_subcluster(test_list)
        self.assertEqual(len(self.ticket.evaluations), 1)

    def test_check_subcluster_4(self):
        test_list = ["AA1", "AA2"]
        self.ticket.subcluster = "AAAAAAA3"
        self.ticket.check_subcluster(test_list)
        self.assertEqual(len(self.ticket.evaluations), 2)

    def test_check_subcluster_5(self):
        test_list = ["AA1", "AA2", "AAAAAAA3"]
        self.ticket.subcluster = "AAAAAAA3"
        self.ticket.check_subcluster(test_list)
        self.assertEqual(len(self.ticket.evaluations), 1)








    def test_check_cluster_1(self):
        test_list = ["AA", "AA"]
        self.ticket.cluster = "none"
        self.ticket.check_cluster(test_list)
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_cluster_2(self):
        test_list = ["AA", "AA"]
        self.ticket.cluster = "AA"
        self.ticket.check_cluster(test_list)
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_cluster_3(self):
        test_list = ["AA", "AA"]
        self.ticket.cluster = "Singleton"
        self.ticket.check_cluster(test_list)
        with self.subTest():
            self.assertEqual(len(self.ticket.evaluations), 0)
        with self.subTest():
            self.assertEqual(self.ticket.cluster, "singleton")

    def test_check_cluster_4(self):
        test_list = ["AA", "AA"]
        self.ticket.cluster = "AB"
        self.ticket.check_cluster(test_list)
        self.assertEqual(len(self.ticket.evaluations), 1)

    def test_check_cluster_5(self):
        test_list = ["AA", "AA"]
        self.ticket.cluster = "AAAAAAB"
        self.ticket.check_cluster(test_list)
        self.assertEqual(len(self.ticket.evaluations), 2)

    def test_check_cluster_6(self):
        test_list = ["AA", "AA", "AAAAAAB"]
        self.ticket.cluster = "AAAAAAB"
        self.ticket.check_cluster(test_list)
        self.assertEqual(len(self.ticket.evaluations), 1)





    def test_check_cluster_subcluster_1(self):
        self.ticket.cluster = "none"
        self.ticket.subcluster = "none"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_cluster_subcluster_2(self):
        self.ticket.cluster = "none"
        self.ticket.subcluster = "AA1"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(len(self.ticket.evaluations), 1)

    def test_check_cluster_subcluster_3(self):
        self.ticket.cluster = "singleton"
        self.ticket.subcluster = "none"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_cluster_subcluster_4(self):
        self.ticket.cluster = "singleton"
        self.ticket.subcluster = "AA1"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(len(self.ticket.evaluations), 1)

    def test_check_cluster_subcluster_5(self):
        self.ticket.cluster = "UNK"
        self.ticket.subcluster = "none"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_cluster_subcluster_6(self):
        self.ticket.cluster = "UNK"
        self.ticket.subcluster = "AA1"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(len(self.ticket.evaluations), 1)

    def test_check_cluster_subcluster_7(self):
        self.ticket.cluster = "AA"
        self.ticket.subcluster = "none"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_cluster_subcluster_8(self):
        self.ticket.cluster = "AA"
        self.ticket.subcluster = "AA1"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_cluster_subcluster_9(self):
        self.ticket.cluster = "AA"
        self.ticket.subcluster = "A1"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(len(self.ticket.evaluations), 1)

    def test_check_cluster_subcluster_10(self):
        self.ticket.cluster = "AA"
        self.ticket.subcluster = "1"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(len(self.ticket.evaluations), 1)

    def test_check_cluster_subcluster_11(self):
        self.ticket.cluster = "AA"
        self.ticket.subcluster = "AAA1"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(len(self.ticket.evaluations), 1)

    def test_check_cluster_subcluster_12(self):
        self.ticket.cluster = "AA"
        self.ticket.subcluster = "A123"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(len(self.ticket.evaluations), 1)





    def test_check_status_1(self):
        test_list = ["final", "draft"]
        self.ticket.status = "none"
        self.ticket.check_status(test_list)
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_status_2(self):
        test_list = ["final", "draft"]
        self.ticket.status = "final"
        self.ticket.check_status(test_list)
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_status_3(self):
        test_list = ["final", "draft"]
        self.ticket.status = "new"
        self.ticket.check_status(test_list)
        self.assertEqual(len(self.ticket.evaluations), 1)

    def test_check_status_4(self):
        test_list = ["final", "draft"]
        self.ticket.status = "newwwwww"
        self.ticket.check_status(test_list)
        self.assertEqual(len(self.ticket.evaluations), 2)

    def test_check_status_5(self):
        test_list = ["final", "draft", "newwwwww"]
        self.ticket.status = "newwwwww"
        self.ticket.check_status(test_list)
        self.assertEqual(len(self.ticket.evaluations), 1)




    def test_check_description_field_1(self):
        test_list = ["product", "notes"]
        self.ticket.description_field = "product"
        self.ticket.check_description_field(test_list)
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_description_field_2(self):
        test_list = ["product", "notes"]
        self.ticket.description_field = "none"
        self.ticket.check_description_field(test_list)
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_description_field_3(self):
        test_list = ["product", "notes"]
        self.ticket.description_field = "new"
        self.ticket.check_description_field(test_list)
        self.assertEqual(len(self.ticket.evaluations), 1)




    #TODO unittest check_accession

    def test_check_accession_1(self):
        """Test that empty is changed to none."""
        self.ticket.accession = ""
        self.ticket.check_accession()
        self.assertEqual(self.ticket.accession, "none")

    def test_check_accession_2(self):
        """Test that whitespace is removed."""
        self.ticket.accession = "   ABC123    "
        self.ticket.check_accession()
        self.assertEqual(self.ticket.accession, "ABC123")

    def test_check_accession_3(self):
        """Test that data is split."""
        self.ticket.accession = "ABC123.456"
        self.ticket.check_accession()
        self.assertEqual(self.ticket.accession, "ABC123")

    def test_check_accession_4(self):
        """Test that data is split."""
        self.ticket.accession = "   .   "
        self.ticket.check_accession()
        self.assertEqual(self.ticket.accession, "none")









    def test_check_annotation_author_1(self):
        self.ticket.annotation_author = "hatfull"
        self.ticket.check_annotation_author()
        with self.subTest():
            self.assertEqual(self.ticket.annotation_author, "1")
        with self.subTest():
            self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_annotation_author_2(self):
        self.ticket.annotation_author = "gbk"
        self.ticket.check_annotation_author()
        with self.subTest():
            self.assertEqual(self.ticket.annotation_author, "0")
        with self.subTest():
            self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_annotation_author_3(self):
        self.ticket.annotation_author = "none"
        self.ticket.check_annotation_author()
        with self.subTest():
            self.assertEqual(self.ticket.annotation_author, "none")
        with self.subTest():
            self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check_annotation_author_4(self):
        self.ticket.annotation_author = "other"
        self.ticket.check_annotation_author()
        with self.subTest():
            self.assertEqual(self.ticket.annotation_author, "error")
        with self.subTest():
            self.assertEqual(len(self.ticket.evaluations), 1)



    def test_check__run_mode_1(self):
        run_mode_dict = {"phagesdb":"abcd","pecaan":"efgh"}
        self.ticket.run_mode = "phagesdb"
        self.ticket.check_run_mode(run_mode_dict)
        self.assertEqual(len(self.ticket.evaluations), 0)

    def test_check__run_mode_2(self):
        run_mode_dict = {"phagesdb":"abcd","pecaan":"efgh"}
        self.ticket.run_mode = "other"
        self.ticket.check_run_mode(run_mode_dict)
        self.assertEqual(len(self.ticket.evaluations), 1)

    def test_check__run_mode_3(self):
        run_mode_dict = {"phagesdb":"abcd","pecaan":"efgh"}
        self.ticket.run_mode = "custom"
        self.ticket.check_run_mode(run_mode_dict)
        self.assertEqual(len(self.ticket.evaluations), 0)





    def test_check_update_ticket_1(self):
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.update_ticket.check_update_ticket(phage_id_set)
        self.assertEqual(len(self.update_ticket.evaluations), 0)

    def test_check_update_ticket_2(self):
        """Primary Phage ID not in set."""
        phage_id_set = set(["L5","RedRock"])
        self.update_ticket.check_update_ticket(phage_id_set)
        self.assertEqual(len(self.update_ticket.evaluations), 1)

    def test_check_update_ticket_3(self):
        """Host is none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.update_ticket.host = "none"
        self.update_ticket.check_update_ticket(phage_id_set)
        self.assertEqual(len(self.update_ticket.evaluations), 1)

    def test_check_update_ticket_4(self):
        """Cluster is none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.update_ticket.cluster = "none"
        self.update_ticket.check_update_ticket(phage_id_set)
        self.assertEqual(len(self.update_ticket.evaluations), 1)

    def test_check_update_ticket_5(self):
        """Status is none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.update_ticket.status = "none"
        self.update_ticket.check_update_ticket(phage_id_set)
        self.assertEqual(len(self.update_ticket.evaluations), 1)

    def test_check_update_ticket_6(self):
        """Description Field is not none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.update_ticket.description_field = "Product"
        self.update_ticket.check_update_ticket(phage_id_set)
        self.assertEqual(len(self.update_ticket.evaluations), 1)

    def test_check_update_ticket_7(self):
        """Secondary Phage ID is not none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.update_ticket.secondary_phage_id = "Trixie"
        self.update_ticket.check_update_ticket(phage_id_set)
        self.assertEqual(len(self.update_ticket.evaluations), 1)

    def test_check_update_ticket_8(self):
        """Annotation Author is 0."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.update_ticket.annotation_author = "0"
        self.update_ticket.check_update_ticket(phage_id_set)
        self.assertEqual(len(self.update_ticket.evaluations), 0)

    def test_check_update_ticket_9(self):
        """Annotation Author is not 1 or 0."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.update_ticket.annotation_author = "none"
        self.update_ticket.check_update_ticket(phage_id_set)
        self.assertEqual(len(self.update_ticket.evaluations), 1)

    def test_check_update_ticket_10(self):
        """Run Mode is not none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.update_ticket.run_mode = "phagesdb"
        self.update_ticket.check_update_ticket(phage_id_set)
        self.assertEqual(len(self.update_ticket.evaluations), 1)








    def test_check_add_ticket_1(self):
        """Standard add ticket."""
        phage_id_set = set(["L5","RedRock"])
        self.add_ticket.check_add_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 0)

    def test_check_add_ticket_2(self):
        """Primary Phage ID already present."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.add_ticket.check_add_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 1)

    def test_check_add_ticket_3(self):
        """Primary Phage ID is none."""
        phage_id_set = set(["L5","RedRock"])
        self.add_ticket.primary_phage_id = "none"
        self.add_ticket.check_add_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 1)

    def test_check_add_ticket_4(self):
        """Host is none."""
        phage_id_set = set(["L5","RedRock"])
        self.add_ticket.host = "none"
        self.add_ticket.check_add_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 1)

    def test_check_add_ticket_5(self):
        """Cluster is none."""
        phage_id_set = set(["L5","RedRock"])
        self.add_ticket.cluster = "none"
        self.add_ticket.check_add_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 1)

    def test_check_add_ticket_6(self):
        """Status is none."""
        phage_id_set = set(["L5","RedRock"])
        self.add_ticket.status = "none"
        self.add_ticket.check_add_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 1)

    def test_check_add_ticket_7(self):
        """Status is Final."""
        phage_id_set = set(["L5","RedRock"])
        self.add_ticket.status = "final"
        self.add_ticket.check_add_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 1)

    def test_check_add_ticket_8(self):
        """Description Field is none."""
        phage_id_set = set(["L5","RedRock"])
        self.add_ticket.description_field = "none"
        self.add_ticket.check_add_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 1)

    def test_check_add_ticket_9(self):
        """Secondary Phage ID is not none."""
        phage_id_set = set(["L5","RedRock"])
        self.add_ticket.secondary_phage_id = "Trixie"
        self.add_ticket.check_add_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 1)

    def test_check_add_ticket_10(self):
        """Annotation Author is 0."""
        phage_id_set = set(["L5","RedRock"])
        self.add_ticket.annotation_author = "0"
        self.add_ticket.check_add_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 0)

    def test_check_add_ticket_11(self):
        """Annotation Author is not 1 or 0."""
        phage_id_set = set(["L5","RedRock"])
        self.add_ticket.annotation_author = "none"
        self.add_ticket.check_add_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 1)

    def test_check_add_ticket_12(self):
        """Run Mode is none."""
        phage_id_set = set(["L5","RedRock"])
        self.add_ticket.run_mode = "none"
        self.add_ticket.check_add_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 1)








    def test_check_remove_ticket_1(self):
        """Standard remove ticket."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.remove_ticket.check_remove_ticket(phage_id_set)
        self.assertEqual(len(self.remove_ticket.evaluations), 0)

    def test_check_remove_ticket_2(self):
        """Primary Phage ID is not none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.remove_ticket.primary_phage_id = "Trixie"
        self.remove_ticket.check_remove_ticket(phage_id_set)
        self.assertEqual(len(self.remove_ticket.evaluations), 1)

    def test_check_remove_ticket_3(self):
        """Host is not none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.remove_ticket.host = "Mycobacterium"
        self.remove_ticket.check_remove_ticket(phage_id_set)
        self.assertEqual(len(self.remove_ticket.evaluations), 1)

    def test_check_remove_ticket_4(self):
        """Cluster is not none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.remove_ticket.cluster = "A"
        self.remove_ticket.check_remove_ticket(phage_id_set)
        self.assertEqual(len(self.remove_ticket.evaluations), 1)

    def test_check_remove_ticket_5(self):
        """Subcluster is not none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.remove_ticket.subcluster = "A2"
        self.remove_ticket.check_remove_ticket(phage_id_set)
        self.assertEqual(len(self.remove_ticket.evaluations), 1)

    def test_check_remove_ticket_6(self):
        """Status is not none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.remove_ticket.status = "final"
        self.remove_ticket.check_remove_ticket(phage_id_set)
        self.assertEqual(len(self.remove_ticket.evaluations), 1)

    def test_check_remove_ticket_7(self):
        """Description Field is not none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.remove_ticket.description_field = "Product"
        self.remove_ticket.check_remove_ticket(phage_id_set)
        self.assertEqual(len(self.remove_ticket.evaluations), 1)

    def test_check_remove_ticket_8(self):
        """Accession is not none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.remove_ticket.accession = "ABC123"
        self.remove_ticket.check_remove_ticket(phage_id_set)
        self.assertEqual(len(self.remove_ticket.evaluations), 1)

    def test_check_remove_ticket_9(self):
        """Annotation Author is not none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.remove_ticket.annotation_author = "1"
        self.remove_ticket.check_remove_ticket(phage_id_set)
        self.assertEqual(len(self.remove_ticket.evaluations), 1)

    def test_check_remove_ticket_10(self):
        """Run Mode is not none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.remove_ticket.run_mode = "phagesdb"
        self.remove_ticket.check_remove_ticket(phage_id_set)
        self.assertEqual(len(self.remove_ticket.evaluations), 1)

    def test_check_remove_ticket_11(self):
        """Secondary Phage ID is not present."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.remove_ticket.secondary_phage_id = "D29"
        self.remove_ticket.check_remove_ticket(phage_id_set)
        self.assertEqual(len(self.remove_ticket.evaluations), 1)









    def test_check_replace_ticket_1(self):
        """Standard replace ticket."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.replace_ticket.check_replace_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 0)

    def test_check_replace_ticket_2(self):
        """Primary Phage ID is none."""
        phage_id_set = set(["Trixie","L5","RedRock","none"])
        self.replace_ticket.primary_phage_id = "none"
        self.replace_ticket.secondary_phage_id = "none"
        self.replace_ticket.check_replace_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 1)

    def test_check_replace_ticket_3(self):
        """Primary Phage ID not present and different from Secondary Phage ID."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.replace_ticket.primary_phage_id = "D29"
        self.replace_ticket.secondary_phage_id = "L5"
        self.replace_ticket.check_replace_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 1)

    def test_check_replace_ticket_4(self):
        """Primary Phage ID is present and different from Secondary Phage ID."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.replace_ticket.secondary_phage_id = "L5"
        self.replace_ticket.check_replace_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 2)

    def test_check_replace_ticket_5(self):
        """Host is none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.replace_ticket.host = "none"
        self.replace_ticket.check_replace_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 1)

    def test_check_replace_ticket_6(self):
        """Cluster is none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.replace_ticket.cluster = "none"
        self.replace_ticket.check_replace_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 1)

    def test_check_replace_ticket_7(self):
        """Status is none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.replace_ticket.status = "none"
        self.replace_ticket.check_replace_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 1)

    def test_check_replace_ticket_8(self):
        """Description Field is none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.replace_ticket.description_field = "none"
        self.replace_ticket.check_replace_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 1)

    def test_check_replace_ticket_9(self):
        """Secondary Phage ID is not present."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.replace_ticket.primary_phage_id = "D29"
        self.replace_ticket.secondary_phage_id = "D29"
        self.replace_ticket.check_replace_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 1)

    def test_check_replace_ticket_10(self):
        """Annotation Author is 0."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.replace_ticket.annotation_author = "0"
        self.replace_ticket.check_replace_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 0)

    def test_check_replace_ticket_11(self):
        """Annotation Author is not 1 or 0."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.replace_ticket.annotation_author = "none"
        self.replace_ticket.check_replace_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 1)

    def test_check_replace_ticket_12(self):
        """Run Mode is none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.replace_ticket.run_mode = "none"
        self.replace_ticket.check_replace_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 1)










    def test_check_ticket_1(self):
        """Check standard update ticket."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.update_ticket.check_ticket(phage_id_set)
        self.assertEqual(len(self.update_ticket.evaluations), 0)

    def test_check_ticket_2(self):
        """Check update ticket with Primary Phage ID not present."""
        phage_id_set = set(["L5","RedRock"])
        self.update_ticket.check_ticket(phage_id_set)
        self.assertEqual(len(self.update_ticket.evaluations), 1)

    def test_check_ticket_3(self):
        """Check standard add ticket."""
        phage_id_set = set(["L5","RedRock"])
        self.add_ticket.check_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 0)

    def test_check_ticket_4(self):
        """Check add ticket with Primary Phage ID already present."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.add_ticket.check_ticket(phage_id_set)
        self.assertEqual(len(self.add_ticket.evaluations), 1)

    def test_check_ticket_5(self):
        """Check standard remove ticket."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.remove_ticket.check_ticket(phage_id_set)
        self.assertEqual(len(self.remove_ticket.evaluations), 0)

    def test_check_ticket_6(self):
        """Check remove ticket with Primary Phage ID not none."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.remove_ticket.primary_phage_id = "Trixie"
        self.remove_ticket.check_ticket(phage_id_set)
        self.assertEqual(len(self.remove_ticket.evaluations), 1)

    def test_check_ticket_7(self):
        """Check standard replace ticket."""
        phage_id_set = set(["Trixie","L5","RedRock"])
        self.replace_ticket.check_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 0)

    def test_check_ticket_8(self):
        """Check replace ticket when Primary Phage ID is none."""
        phage_id_set = set(["Trixie","L5","RedRock","none"])
        self.replace_ticket.primary_phage_id = "none"
        self.replace_ticket.secondary_phage_id = "none"
        self.replace_ticket.check_ticket(phage_id_set)
        self.assertEqual(len(self.replace_ticket.evaluations), 1)

    def test_check_ticket_9(self):
        """Check non-standard type of ticket."""
        phage_id_set = set(["L5","RedRock"])
        self.update_ticket.type = "other"
        self.update_ticket.check_ticket(phage_id_set)
        self.assertEqual(len(self.update_ticket.evaluations), 0)


    # TODO unit test validate_ticket method


if __name__ == '__main__':
    unittest.main()
