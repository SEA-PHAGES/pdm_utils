""" Unit tests for the ticket class"""



from classes import Ticket
import unittest



class TestGenomeTicketClass(unittest.TestCase):


    def setUp(self):

        # Empty ticket to test simple methods
        self.ticket = Ticket.GenomeTicket()

        # Standard update ticket
        self.update_ticket = Ticket.GenomeTicket()
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
        self.add_ticket = Ticket.GenomeTicket()
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
        self.remove_ticket = Ticket.GenomeTicket()
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
        self.replace_ticket = Ticket.GenomeTicket()
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




    def test_clear_retrieve_status_1(self):
        self.ticket.host = "retrieve"
        self.ticket.clear_retrieve_status()
        with self.subTest():
            self.assertEqual(self.ticket.host, "none")
        with self.subTest():
            self.assertEqual(self.ticket.evaluations[0].status, "error")

    def test_clear_retrieve_status_2(self):
        self.ticket.cluster = "retrieve"
        self.ticket.clear_retrieve_status()
        with self.subTest():
            self.assertEqual(self.ticket.cluster, "none")
        with self.subTest():
            self.assertEqual(self.ticket.evaluations[1].status, "error")

    def test_clear_retrieve_status_3(self):
        self.ticket.subcluster = "retrieve"
        self.ticket.clear_retrieve_status()
        with self.subTest():
            self.assertEqual(self.ticket.subcluster, "none")
        with self.subTest():
            self.assertEqual(self.ticket.evaluations[2].status, "error")

    def test_clear_retrieve_status_4(self):
        self.ticket.accession = "retrieve"
        self.ticket.clear_retrieve_status()
        with self.subTest():
            self.assertEqual(self.ticket.accession, "none")
        with self.subTest():
            self.assertEqual(self.ticket.evaluations[3].status, "error")





    def test_check_type_1(self):
        """Check that no error is produced if the
        type is present in the first set."""
        set1 = set(["add", "remove"])
        self.ticket.type = "add"
        self.ticket.check_type(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_type_2(self):
        """Check that an error is produced if the
        type is not present in the first set."""
        set1 = set(["add", "remove"])
        self.ticket.type = "none"
        self.ticket.type = "Trixie"
        self.ticket.check_type(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_primary_phage_id_1(self):
        """Check that no error is produced if the
        primary_phage_id is expected to be in the first set,
        and is present."""
        set1 = set(["Trixie", "L5"])
        set2 = set(["none"])
        self.ticket.primary_phage_id = "Trixie"
        self.ticket.check_primary_phage_id(set1, set2, "first")
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_primary_phage_id_2(self):
        """Check that an error is produced if the
        primary_phage_id is expected to be in neither set,
        and is present in one."""
        set1 = set(["Trixie", "L5"])
        set2 = set(["none"])
        self.ticket.primary_phage_id = "Trixie"
        self.ticket.check_primary_phage_id(set1, set2, "neither")
        self.assertEqual(self.ticket.evaluations[0].status, "error")



    def test_check_secondary_phage_id_1(self):
        """Check that no error is produced if the
        secondary_phage_id is expected to be in the first set,
        and is present."""
        set1 = set(["Trixie", "L5"])
        set2 = set(["none"])
        self.ticket.secondary_phage_id = "Trixie"
        self.ticket.check_secondary_phage_id(set1, set2, "first")
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_secondary_phage_id_2(self):
        """Check that an error is produced if the
        secondary_phage_id is expected to be in neither set,
        and is present in one."""
        set1 = set(["Trixie", "L5"])
        set2 = set(["none"])
        self.ticket.secondary_phage_id = "Trixie"
        self.ticket.check_secondary_phage_id(set1, set2, "neither")
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_host_1(self):
        """Check that no error is produced if the
        host is expected to be in the first set,
        and is present."""
        set1 = set(["Mycobacterium", "Gordonia"])
        set2 = set(["none"])
        self.ticket.host = "Mycobacterium"
        self.ticket.check_host(set1, set2, "first")
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_host_2(self):
        """Check that an error is produced if the
        host is expected to be in neither set,
        and is present in one."""
        set1 = set(["Mycobacterium", "Gordonia"])
        set2 = set(["none"])
        self.ticket.host = "Mycobacterium"
        self.ticket.check_host(set1, set2, "neither")
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_subcluster_1(self):
        """Check that no error is produced if the
        subcluster is expected to be in the first set,
        and is present."""
        set1 = set(["A2", "B1"])
        set2 = set(["none"])
        self.ticket.subcluster = "A2"
        self.ticket.check_subcluster(set1, set2, "first")
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_subcluster_2(self):
        """Check that an error is produced if the
        subcluster is expected to be in neither set,
        and is present in one."""
        set1 = set(["A2", "B1"])
        set2 = set(["none"])
        self.ticket.subcluster = "A2"
        self.ticket.check_subcluster(set1, set2, "neither")
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_cluster_1(self):
        """Check that no error is produced if the
        cluster is expected to be in the first set,
        and is present."""
        set1 = set(["A", "B"])
        set2 = set(["none"])
        self.ticket.cluster = "A"
        self.ticket.check_cluster(set1, set2, "first")
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_cluster_2(self):
        """Check that an error is produced if the
        cluster is expected to be in neither set,
        and is present in one."""
        set1 = set(["A", "B"])
        set2 = set(["none"])
        self.ticket.cluster = "A"
        self.ticket.check_cluster(set1, set2, "neither")
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_subcluster_structure_1(self):
        """Check that no error is produced if the
        non-empty subcluster is not structured correctly."""
        self.ticket.subcluster = "A1"
        self.ticket.check_subcluster_structure()
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_subcluster_structure_2(self):
        """Check that an error is produced if the
        non-empty subcluster is not structured correctly."""
        self.ticket.subcluster = "A"
        self.ticket.check_subcluster_structure()
        self.assertEqual(self.ticket.evaluations[0].status, "error")

    def test_check_subcluster_structure_3(self):
        """Check that no error is produced if the
        subcluster is empty."""
        self.ticket.subcluster = "none"
        self.ticket.check_subcluster_structure()
        self.assertEqual(self.ticket.evaluations[0].status, "not_evaluated")




    def test_check_cluster_structure_1(self):
        """Check that no error is produced if the
        non-empty cluster is not structured correctly."""
        self.ticket.cluster = "A"
        self.ticket.check_cluster_structure()
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_cluster_structure_2(self):
        """Check that an error is produced if the
        non-empty cluster is not structured correctly."""
        self.ticket.cluster = "A1"
        self.ticket.check_cluster_structure()
        self.assertEqual(self.ticket.evaluations[0].status, "error")

    def test_check_cluster_structure_3(self):
        """Check that no error is produced if the
        cluster is empty."""
        self.ticket.cluster = "none"
        self.ticket.check_cluster_structure()
        self.assertEqual(self.ticket.evaluations[0].status, "not_evaluated")




    def test_check_cluster_subcluster_1(self):
        """Check that compatible Cluster and subcluster
        do not produce an error."""
        self.ticket.cluster = "A"
        self.ticket.subcluster = "A1"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_cluster_subcluster_2(self):
        """Check that incompatible Cluster and subcluster
        produce an error."""
        self.ticket.cluster = "A"
        self.ticket.subcluster = "B1"
        self.ticket.check_cluster_subcluster()
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_status_1(self):
        """Check that no error is produced if the
        status is expected to be in the first set,
        and is present."""
        set1 = set(["final", "draft"])
        set2 = set(["none"])
        self.ticket.status = "final"
        self.ticket.check_status(set1, set2, "first")
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_status_2(self):
        """Check that an error is produced if the
        status is expected to be in neither set,
        and is present in one."""
        set1 = set(["final", "draft"])
        set2 = set(["none"])
        self.ticket.status = "final"
        self.ticket.check_status(set1, set2, "neither")
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_description_field_1(self):
        """Check that no error is produced if the
        description_field is expected to be in the first set,
        and is present."""
        set1 = set(["product", "function"])
        set2 = set(["none"])
        self.ticket.description_field = "product"
        self.ticket.check_description_field(set1, set2, "first")
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_description_field_2(self):
        """Check that an error is produced if the
        description_field is expected to be in neither set,
        and is present in one."""
        set1 = set(["product", "function"])
        set2 = set(["none"])
        self.ticket.description_field = "product"
        self.ticket.check_description_field(set1, set2, "neither")
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_accession_1(self):
        """Check that no error is produced if the
        accession is expected to be in the first set,
        and is present."""
        set1 = set(["product", "function"])
        set2 = set(["none"])
        self.ticket.accession = "product"
        self.ticket.check_accession(set1, set2, "first")
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_accession_2(self):
        """Check that an error is produced if the
        accession is expected to be in neither set,
        and is present in one."""
        set1 = set(["product", "function"])
        set2 = set(["none"])
        self.ticket.accession = "product"
        self.ticket.check_accession(set1, set2, "neither")
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_annotation_author_1(self):
        """Check that no error is produced if the
        annotation_author is expected to be in the first set,
        and is present."""
        set1 = set(["hatfull", "gbk"])
        set2 = set(["none"])
        self.ticket.annotation_author = "hatfull"
        self.ticket.check_annotation_author(set1, set2, "first")
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_annotation_author_2(self):
        """Check that an error is produced if the
        annotation_author is expected to be in neither set,
        and is present in one."""
        set1 = set(["hatfull", "gbk"])
        set2 = set(["none"])
        self.ticket.annotation_author = "hatfull"
        self.ticket.check_annotation_author(set1, set2, "neither")
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_run_mode_1(self):
        """Check that no error is produced if the
        run_mode is expected to be in the first set,
        and is present."""
        set1 = set(["phagesdb", "pecaan"])
        set2 = set(["none"])
        self.ticket.run_mode = "phagesdb"
        self.ticket.check_run_mode(set1, set2, "first")
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_run_mode_2(self):
        """Check that an error is produced if the
        run_mode is expected to be in neither set,
        and is present in one."""
        set1 = set(["phagesdb", "pecaan"])
        set2 = set(["none"])
        self.ticket.run_mode = "phagesdb"
        self.ticket.check_run_mode(set1, set2, "neither")
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_primary_secondary_phage_ids_1(self):
        """Check that an error is produced if the
        primary_phage_id is not the same as the secondary_phage_id."""
        self.ticket.primary_phage_id = "Trixie"
        self.ticket.secondary_phage_id = "L5"
        self.ticket.check_primary_secondary_phage_ids()
        self.assertEqual(self.ticket.evaluations[0].status, "error")

    def test_check_primary_secondary_phage_ids_2(self):
        """Check that no error is produced if the
        primary_phage_id is the same as the secondary_phage_id."""
        self.ticket.primary_phage_id = "Trixie"
        self.ticket.secondary_phage_id = "Trixie"
        self.ticket.check_primary_secondary_phage_ids()
        self.assertEqual(self.ticket.evaluations[0].status, "correct")












if __name__ == '__main__':
    unittest.main()
