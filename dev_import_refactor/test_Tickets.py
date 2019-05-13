""" Unit tests for the ticket class"""



import Ticket
import unittest



class TestTicketClass(unittest.TestCase):

    #
    def setUp(self):
        self.ticket = Ticket.ImportTicket()
    #     self.ticket.type = ""
    #     self.ticket.primary_phage_id = "Trixie"
    #     self.ticket.host = "Mycobacterium"
    #     self.ticket.cluster = "A"
    #     self.ticket.subcluster = "A2"
    #     self.ticket.status = "Final"
    #     self.ticket.annotation_author = "Hatfull"
    #     self.ticket.description_field = "Product"
    #     self.ticket.accession = "None"
    #     self.ticket.run_mode = "phagesdb"
    #     self.ticket.secondary_phage_id = "Trixie"
    #     self.ticket.match_strategy = "phageid"
    #     self.ticket.evaluations = []

    # def tearDown(self):
    #     self.ticket = ""

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

    def test_check_annotation_author_3(self):
        self.ticket.annotation_author = "other"
        self.ticket.check_annotation_author()
        with self.subTest():
            self.assertEqual(self.ticket.annotation_author, "error")
        with self.subTest():
            self.assertEqual(len(self.ticket.evaluations), 1)


    #TODO next = unittest for check_run_mode




if __name__ == '__main__':
    unittest.main()
