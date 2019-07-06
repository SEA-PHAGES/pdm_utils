""" Unit tests for the ticket class"""



from classes import Ticket
import unittest



class TestGenomeTicketClass(unittest.TestCase):


    def setUp(self):

        # Empty ticket to test simple methods
        self.ticket = Ticket.GenomeTicket()




        ###TODO below eventually move to test_evaluate.

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
        ###TODO above eventually move to test_evaluate.




    def test_set_type_1(self):
        """Check that value is lowercased when not 'none'."""
        value = "ADD"
        self.ticket.set_type(value)
        self.assertEqual(self.ticket.type, "add")

    def test_set_type_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_type(value)
        self.assertEqual(self.ticket.type, "none")




    def test_set_primary_phage_id_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "Trixie"
        self.ticket.set_primary_phage_id(value)
        self.assertEqual(self.ticket.primary_phage_id, "Trixie")

    def test_set_primary_phage_id_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_primary_phage_id(value)
        self.assertEqual(self.ticket.primary_phage_id, "none")




    def test_set_host_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "Mycobacterium"
        self.ticket.set_host(value)
        self.assertEqual(self.ticket.host, "Mycobacterium")

    def test_set_host_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_host(value)
        self.assertEqual(self.ticket.host, "none")

    def test_set_host_3(self):
        """Check that value is lowercased when 'retrieve'."""
        value = "RETRIEVE"
        self.ticket.set_host(value)
        self.assertEqual(self.ticket.host, "retrieve")




    def test_set_cluster_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "A"
        self.ticket.set_cluster(value)
        self.assertEqual(self.ticket.cluster, "A")

    def test_set_cluster_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_cluster(value)
        self.assertEqual(self.ticket.cluster, "none")




    def test_set_subcluster_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "A1"
        self.ticket.set_subcluster(value)
        self.assertEqual(self.ticket.subcluster, "A1")

    def test_set_subcluster_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_subcluster(value)
        self.assertEqual(self.ticket.subcluster, "none")




    def test_set_status_1(self):
        """Check that value is lowercased when not 'none'."""
        value = "Final"
        self.ticket.set_status(value)
        self.assertEqual(self.ticket.status, "final")

    def test_set_status_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_status(value)
        self.assertEqual(self.ticket.status, "none")




    def test_set_description_field_1(self):
        """Check that value is lowercased when not 'none'."""
        value = "Product"
        self.ticket.set_description_field(value)
        self.assertEqual(self.ticket.description_field, "product")

    def test_set_description_field_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_description_field(value)
        self.assertEqual(self.ticket.description_field, "none")




    def test_set_accession_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "ABC123"
        self.ticket.set_accession(value)
        self.assertEqual(self.ticket.accession, "ABC123")

    def test_set_accession_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_accession(value)
        self.assertEqual(self.ticket.accession, "none")




    def test_set_annotation_author_1(self):
        """Check that value is lowercased when not 'none'."""
        value = "Hatfull"
        self.ticket.set_annotation_author(value)
        self.assertEqual(self.ticket.annotation_author, "hatfull")

    def test_set_annotation_author_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_annotation_author(value)
        self.assertEqual(self.ticket.annotation_author, "none")




    def test_set_secondary_phage_id_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "Trixie"
        self.ticket.set_secondary_phage_id(value)
        self.assertEqual(self.ticket.secondary_phage_id, "Trixie")

    def test_set_secondary_phage_id_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_secondary_phage_id(value)
        self.assertEqual(self.ticket.secondary_phage_id, "none")




    def test_set_run_mode_1(self):
        """Check that value is lowercased when not 'none'."""
        value = "PhagesDB"
        self.ticket.set_run_mode(value)
        self.assertEqual(self.ticket.run_mode, "phagesdb")

    def test_set_run_mode_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_run_mode(value)
        self.assertEqual(self.ticket.run_mode, "none")




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
        self.ticket.check_type(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_primary_phage_id_1(self):
        """Check that no error is produced if the
        primary_phage_id is not present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.primary_phage_id = "Trixie"
        self.ticket.check_primary_phage_id(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_primary_phage_id_2(self):
        """Check that no error is produced if the
        primary_phage_id is present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.primary_phage_id = "none"
        self.ticket.check_primary_phage_id(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_secondary_phage_id_1(self):
        """Check that no error is produced if the
        secondary_phage_id is not present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.secondary_phage_id = "Trixie"
        self.ticket.check_secondary_phage_id(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_secondary_phage_id_2(self):
        """Check that no error is produced if the
        secondary_phage_id is present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.secondary_phage_id = "none"
        self.ticket.check_secondary_phage_id(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_host_1(self):
        """Check that no error is produced if the
        host is not present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.host = "Mycobacterium"
        self.ticket.check_host(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_host_2(self):
        """Check that no error is produced if the
        host is present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.host = "none"
        self.ticket.check_host(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_subcluster_1(self):
        """Check that no error is produced if the
        subcluster is not present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.subcluster = "Mycobacterium"
        self.ticket.check_subcluster(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_subcluster_2(self):
        """Check that no error is produced if the
        subcluster is present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.subcluster = "none"
        self.ticket.check_subcluster(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_cluster_1(self):
        """Check that no error is produced if the
        cluster is not present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.cluster = "Mycobacterium"
        self.ticket.check_cluster(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_cluster_2(self):
        """Check that no error is produced if the
        cluster is present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.cluster = "none"
        self.ticket.check_cluster(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_status_1(self):
        """Check that no error is produced if the
        status is not present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.status = "Mycobacterium"
        self.ticket.check_status(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_status_2(self):
        """Check that no error is produced if the
        status is present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.status = "none"
        self.ticket.check_status(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_description_field_1(self):
        """Check that no error is produced if the
        description_field is not present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.description_field = "Mycobacterium"
        self.ticket.check_description_field(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_description_field_2(self):
        """Check that no error is produced if the
        description_field is present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.description_field = "none"
        self.ticket.check_description_field(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_accession_1(self):
        """Check that no error is produced if the
        accession is not present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.accession = "Mycobacterium"
        self.ticket.check_accession(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_accession_2(self):
        """Check that no error is produced if the
        accession is present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.accession = "none"
        self.ticket.check_accession(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_annotation_author_1(self):
        """Check that no error is produced if the
        annotation_author is not present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.annotation_author = "Mycobacterium"
        self.ticket.check_annotation_author(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_annotation_author_2(self):
        """Check that no error is produced if the
        annotation_author is present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.annotation_author = "none"
        self.ticket.check_annotation_author(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_run_mode_1(self):
        """Check that no error is produced if the
        run_mode is not present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.run_mode = "Mycobacterium"
        self.ticket.check_run_mode(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_run_mode_2(self):
        """Check that no error is produced if the
        run_mode is present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.run_mode = "none"
        self.ticket.check_run_mode(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "error")














if __name__ == '__main__':
    unittest.main()
