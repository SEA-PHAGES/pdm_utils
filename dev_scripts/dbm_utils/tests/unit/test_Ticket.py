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
        self.update_ticket.phage_id = "Trixie"
        self.update_ticket.host_genus = "Mycobacterium"
        self.update_ticket.cluster = "A"
        self.update_ticket.annotation_status = "Final"
        self.update_ticket.description_field = "none"
        self.update_ticket.annotation_author = "1"
        self.update_ticket.run_mode = "none"



        # Standard add ticket
        self.add_ticket = Ticket.GenomeTicket()
        self.add_ticket.type = "add"
        self.add_ticket.phage_id = "Trixie"
        self.add_ticket.host_genus = "Mycobacterium"
        self.add_ticket.cluster = "A"
        self.add_ticket.annotation_status = "draft"
        self.add_ticket.description_field = "Product"
        self.add_ticket.annotation_author = "1"
        self.add_ticket.run_mode = "phagesdb"


        # Standard remove ticket
        self.remove_ticket = Ticket.GenomeTicket()
        self.remove_ticket.type = "remove"
        self.remove_ticket.phage_id = "none"
        self.remove_ticket.host_genus = "none"
        self.remove_ticket.cluster = "none"
        self.remove_ticket.subcluster = "none"
        self.remove_ticket.annotation_status = "none"
        self.remove_ticket.description_field = "none"
        self.remove_ticket.accession = "none"
        self.remove_ticket.annotation_author = "none"
        self.remove_ticket.run_mode = "none"


        # Standard replace ticket
        self.replace_ticket = Ticket.GenomeTicket()
        self.replace_ticket.type = "replace"
        self.replace_ticket.phage_id = "Trixie"
        self.replace_ticket.host_genus = "Mycobacterium"
        self.replace_ticket.cluster = "A"
        self.replace_ticket.annotation_status = "draft"
        self.replace_ticket.description_field = "Product"
        self.replace_ticket.annotation_author = "1"
        self.replace_ticket.run_mode = "phagesdb"
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




    def test_set_phage_id_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "Trixie"
        self.ticket.set_phage_id(value)
        self.assertEqual(self.ticket.phage_id, "Trixie")

    def test_set_phage_id_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_phage_id(value)
        self.assertEqual(self.ticket.phage_id, "none")




    def test_set_host_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "Mycobacterium"
        self.ticket.set_host(value)
        self.assertEqual(self.ticket.host_genus, "Mycobacterium")

    def test_set_host_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_host(value)
        self.assertEqual(self.ticket.host_genus, "none")

    def test_set_host_3(self):
        """Check that value is lowercased when 'retrieve'."""
        value = "RETRIEVE"
        self.ticket.set_host(value)
        self.assertEqual(self.ticket.host_genus, "retrieve")




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




    def test_set_annotation_status_1(self):
        """Check that value is lowercased when not 'none'."""
        value = "Final"
        self.ticket.set_annotation_status(value)
        self.assertEqual(self.ticket.annotation_status, "final")

    def test_set_annotation_status_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_annotation_status(value)
        self.assertEqual(self.ticket.annotation_status, "none")




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




    def test_set_annotation_qc_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "ABC"
        self.ticket.set_annotation_qc(value)
        self.assertEqual(self.ticket.annotation_qc, "ABC")

    def test_set_annotation_qc_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_annotation_qc(value)
        self.assertEqual(self.ticket.annotation_qc, "none")




    def test_set_retrieve_record_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "ABC"
        self.ticket.set_retrieve_record(value)
        self.assertEqual(self.ticket.retrieve_record, "ABC")

    def test_set_retrieve_record_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.ticket.set_retrieve_record(value)
        self.assertEqual(self.ticket.retrieve_record, "none")




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




    def test_set_value_flag_1(self):
        """Verify that the 'retain' setting is set to True."""
        self.ticket._value_flag = False
        self.ticket.cluster = "retain"
        self.ticket.set_value_flag("retain")
        self.assertTrue(self.ticket._value_flag)

    def test_set_value_flag_2(self):
        """Verify that the 'retain' setting is set to False."""
        self.ticket._value_flag = True
        self.ticket.cluster = "A"
        self.ticket.set_value_flag("retain")
        self.assertFalse(self.ticket._value_flag)











    def test_check_parsed_fields_1(self):
        """Check that no error is produced if the
        correct number of fields were parsed."""
        self.ticket._parsed_fields = 11
        self.ticket.check_parsed_fields()
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_parsed_fields_2(self):
        """Check that an error is produced if the
        incorrect number of fields were parsed."""
        self.ticket._parsed_fields = 10
        self.ticket.check_parsed_fields()
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_type_1(self):
        """Check that no error is produced if the
        type is present in the set and expected to be in the set."""
        set1 = set(["add", "remove"])
        self.ticket.type = "add"
        self.ticket.check_type(set1, True)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_type_2(self):
        """Check that an error is produced if the
        type is not present in the set and expected to be in the set."""
        set1 = set(["add", "remove"])
        self.ticket.type = "none"
        self.ticket.check_type(set1)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_phage_id_1(self):
        """Check that no error is produced if the
        phage_id is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.phage_id = "Trixie"
        self.ticket.check_phage_id(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_phage_id_2(self):
        """Check that an error is produced if the
        phage_id is present in the empty/null set."""
        set1 = set(["none"])
        self.ticket.phage_id = "none"
        self.ticket.check_phage_id(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_host_genus_1(self):
        """Check that no error is produced if the
        host_genus is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.host_genus = "Mycobacterium"
        self.ticket.check_host_genus(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_host_genus_2(self):
        """Check that an error is produced if the
        host_genus is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.host_genus = "none"
        self.ticket.check_host_genus(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_subcluster_1(self):
        """Check that no error is produced if the
        subcluster is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.subcluster = "Mycobacterium"
        self.ticket.check_subcluster(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_subcluster_2(self):
        """Check that an error is produced if the
        subcluster is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.subcluster = "none"
        self.ticket.check_subcluster(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_cluster_1(self):
        """Check that no error is produced if the
        cluster is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.cluster = "Mycobacterium"
        self.ticket.check_cluster(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_cluster_2(self):
        """Check that an error is produced if the
        cluster is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.cluster = "none"
        self.ticket.check_cluster(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_annotation_status_1(self):
        """Check that no error is produced if the
        annotation_status is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.annotation_status = "Mycobacterium"
        self.ticket.check_annotation_status(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_annotation_status_2(self):
        """Check that an error is produced if the
        annotation_status is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.annotation_status = "none"
        self.ticket.check_annotation_status(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_description_field_1(self):
        """Check that no error is produced if the
        description_field is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.description_field = "Mycobacterium"
        self.ticket.check_description_field(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_description_field_2(self):
        """Check that an error is produced if the
        description_field is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.description_field = "none"
        self.ticket.check_description_field(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_accession_1(self):
        """Check that no error is produced if the
        accession is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.accession = "Mycobacterium"
        self.ticket.check_accession(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_accession_2(self):
        """Check that an error is produced if the
        accession is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.accession = "none"
        self.ticket.check_accession(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_annotation_author_1(self):
        """Check that no error is produced if the
        annotation_author is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.annotation_author = "Mycobacterium"
        self.ticket.check_annotation_author(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_annotation_author_2(self):
        """Check that an error is produced if the
        annotation_author is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.annotation_author = "none"
        self.ticket.check_annotation_author(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_run_mode_1(self):
        """Check that no error is produced if the
        run_mode is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.run_mode = "Mycobacterium"
        self.ticket.check_run_mode(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_run_mode_2(self):
        """Check that an error is produced if the
        run_mode is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.ticket.run_mode = "none"
        self.ticket.check_run_mode(set1, False)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_duplicate_phage_id_1(self):
        """Check that no error is produced if the
        phage_id is not present in the set of duplicated values."""
        dupe_set = set(["Trixie", "L5"])
        self.ticket.phage_id = "D29"
        self.ticket.check_duplicate_phage_id(dupe_set)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_duplicate_phage_id_2(self):
        """Check that an error is produced if the
        phage_id is present in the set of duplicated values."""
        dupe_set = set(["Trixie", "L5"])
        self.ticket.phage_id = "Trixie"
        self.ticket.check_duplicate_phage_id(dupe_set)
        self.assertEqual(self.ticket.evaluations[0].status, "error")







    def test_check_duplicate_accession_1(self):
        """Check that no error is produced if the
        accession is not present in the set of duplicated values."""
        dupe_set = set(["ABC123", "EFG456"])
        self.ticket.accession = "XYZ789"
        self.ticket.check_duplicate_accession(dupe_set)
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_duplicate_accession_2(self):
        """Check that an error is produced if the
        accession is present in the set of duplicated values."""
        dupe_set = set(["ABC123", "EFG456"])
        self.ticket.accession = "ABC123"
        self.ticket.check_duplicate_accession(dupe_set)
        self.assertEqual(self.ticket.evaluations[0].status, "error")




    def test_check_compatible_type_and_annotation_status_1(self):
        """Verify that no error is produced with "add" type and "draft"
        annotation_status."""
        self.ticket.type = "add"
        self.ticket.annotation_status = "draft"
        self.ticket.check_compatible_type_and_annotation_status()
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_compatible_type_and_annotation_status_2(self):
        """Verify that an error is produced with "add" type and "final"
        annotation_status."""
        self.ticket.type = "add"
        self.ticket.annotation_status = "final"
        self.ticket.check_compatible_type_and_annotation_status()
        self.assertEqual(self.ticket.evaluations[0].status, "error")

    def test_check_compatible_type_and_annotation_status_3(self):
        """Verify that no error is produced with "replace" type and "final"
        annotation_status."""
        self.ticket.type = "replace"
        self.ticket.annotation_status = "final"
        self.ticket.check_compatible_type_and_annotation_status()
        self.assertEqual(self.ticket.evaluations[0].status, "correct")

    def test_check_compatible_type_and_annotation_status_4(self):
        """Verify that no error is produced with "replace" type and "draft"
        annotation_status."""
        self.ticket.type = "replace"
        self.ticket.annotation_status = "draft"
        self.ticket.check_compatible_type_and_annotation_status()
        self.assertEqual(self.ticket.evaluations[0].status, "error")






if __name__ == '__main__':
    unittest.main()
