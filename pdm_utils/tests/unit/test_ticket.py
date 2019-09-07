""" Unit tests for the ticket class"""



from pdm_utils.classes import ticket
import unittest



class TestGenomeTicketClass(unittest.TestCase):


    def setUp(self):

        # Empty ticket to test simple methods
        self.tkt = ticket.GenomeTicket()




        ###TODO below eventually move to test_evaluate.

        # Standard update ticket
        self.update_ticket = ticket.GenomeTicket()
        self.update_ticket.type = "update"
        self.update_ticket.phage_id = "Trixie"
        self.update_ticket.host_genus = "Mycobacterium"
        self.update_ticket.cluster = "A"
        self.update_ticket.annotation_status = "Final"
        self.update_ticket.description_field = "none"
        self.update_ticket.annotation_author = "1"
        self.update_ticket.run_mode = "none"



        # Standard add ticket
        self.add_ticket = ticket.GenomeTicket()
        self.add_ticket.type = "add"
        self.add_ticket.phage_id = "Trixie"
        self.add_ticket.host_genus = "Mycobacterium"
        self.add_ticket.cluster = "A"
        self.add_ticket.annotation_status = "draft"
        self.add_ticket.description_field = "Product"
        self.add_ticket.annotation_author = "1"
        self.add_ticket.run_mode = "phagesdb"


        # Standard remove ticket
        self.remove_ticket = ticket.GenomeTicket()
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
        self.replace_ticket = ticket.GenomeTicket()
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
        self.tkt.set_type(value)
        self.assertEqual(self.tkt.type, "add")

    def test_set_type_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_type(value)
        self.assertEqual(self.tkt.type, "none")




    def test_set_phage_id_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "Trixie"
        self.tkt.set_phage_id(value)
        self.assertEqual(self.tkt.phage_id, "Trixie")

    def test_set_phage_id_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_phage_id(value)
        self.assertEqual(self.tkt.phage_id, "none")




    def test_set_host_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "Mycobacterium"
        self.tkt.set_host(value)
        self.assertEqual(self.tkt.host_genus, "Mycobacterium")

    def test_set_host_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_host(value)
        self.assertEqual(self.tkt.host_genus, "none")

    def test_set_host_3(self):
        """Check that value is lowercased when 'retrieve'."""
        value = "RETRIEVE"
        self.tkt.set_host(value)
        self.assertEqual(self.tkt.host_genus, "retrieve")




    def test_set_cluster_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "A"
        self.tkt.set_cluster(value)
        self.assertEqual(self.tkt.cluster, "A")

    def test_set_cluster_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_cluster(value)
        self.assertEqual(self.tkt.cluster, "none")




    def test_set_subcluster_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "A1"
        self.tkt.set_subcluster(value)
        self.assertEqual(self.tkt.subcluster, "A1")

    def test_set_subcluster_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_subcluster(value)
        self.assertEqual(self.tkt.subcluster, "none")




    def test_set_annotation_status_1(self):
        """Check that value is lowercased when not 'none'."""
        value = "Final"
        self.tkt.set_annotation_status(value)
        self.assertEqual(self.tkt.annotation_status, "final")

    def test_set_annotation_status_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_annotation_status(value)
        self.assertEqual(self.tkt.annotation_status, "none")




    def test_set_description_field_1(self):
        """Check that value is lowercased when not 'none'."""
        value = "Product"
        self.tkt.set_description_field(value)
        self.assertEqual(self.tkt.description_field, "product")

    def test_set_description_field_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_description_field(value)
        self.assertEqual(self.tkt.description_field, "none")




    def test_set_accession_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "ABC123"
        self.tkt.set_accession(value)
        self.assertEqual(self.tkt.accession, "ABC123")

    def test_set_accession_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_accession(value)
        self.assertEqual(self.tkt.accession, "none")




    def test_set_annotation_author_1(self):
        """Check that value is lowercased when not 'none'."""
        value = "Hatfull"
        self.tkt.set_annotation_author(value)
        self.assertEqual(self.tkt.annotation_author, "hatfull")

    def test_set_annotation_author_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_annotation_author(value)
        self.assertEqual(self.tkt.annotation_author, "none")




    def test_set_annotation_qc_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "ABC"
        self.tkt.set_annotation_qc(value)
        self.assertEqual(self.tkt.annotation_qc, "ABC")

    def test_set_annotation_qc_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_annotation_qc(value)
        self.assertEqual(self.tkt.annotation_qc, "none")




    def test_set_retrieve_record_1(self):
        """Check that value is not lowercased when not 'none'."""
        value = "ABC"
        self.tkt.set_retrieve_record(value)
        self.assertEqual(self.tkt.retrieve_record, "ABC")

    def test_set_retrieve_record_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_retrieve_record(value)
        self.assertEqual(self.tkt.retrieve_record, "none")




    def test_set_run_mode_1(self):
        """Check that value is lowercased when not 'none'."""
        value = "PhagesDB"
        self.tkt.set_run_mode(value)
        self.assertEqual(self.tkt.run_mode, "phagesdb")

    def test_set_run_mode_2(self):
        """Check that value is lowercased when 'none'."""
        value = "NONE"
        self.tkt.set_run_mode(value)
        self.assertEqual(self.tkt.run_mode, "none")




    def test_set_value_flag_1(self):
        """Verify that the 'retain' setting is set to True."""
        self.tkt._value_flag = False
        self.tkt.cluster = "retain"
        self.tkt.set_value_flag("retain")
        self.assertTrue(self.tkt._value_flag)

    def test_set_value_flag_2(self):
        """Verify that the 'retain' setting is set to False."""
        self.tkt._value_flag = True
        self.tkt.cluster = "A"
        self.tkt.set_value_flag("retain")
        self.assertFalse(self.tkt._value_flag)




    def test_check_parsed_fields_1(self):
        """Check that no error is produced if the
        correct number of fields were parsed."""
        self.tkt._parsed_fields = 11
        self.tkt.check_parsed_fields(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_parsed_fields_2(self):
        """Check that an error is produced if the
        incorrect number of fields were parsed."""
        self.tkt._parsed_fields = 10
        self.tkt.check_parsed_fields()
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_type_1(self):
        """Check that no error is produced if the
        type is present in the set and expected to be in the set."""
        set1 = set(["add", "remove"])
        self.tkt.type = "add"
        self.tkt.check_type(set1, True, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_type_2(self):
        """Check that an error is produced if the
        type is not present in the set and expected to be in the set."""
        set1 = set(["add", "remove"])
        self.tkt.type = "none"
        self.tkt.check_type(set1)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_phage_id_1(self):
        """Check that no error is produced if the
        phage_id is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.phage_id = "Trixie"
        self.tkt.check_phage_id(set1, False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")



    def test_check_phage_id_2(self):
        """Check that an error is produced if the
        phage_id is present in the empty/null set."""
        set1 = set(["none"])
        self.tkt.phage_id = "none"
        self.tkt.check_phage_id(set1, False)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_host_genus_1(self):
        """Check that no error is produced if the
        host_genus is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.host_genus = "Mycobacterium"
        self.tkt.check_host_genus(set1, False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_host_genus_2(self):
        """Check that an error is produced if the
        host_genus is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.host_genus = "none"
        self.tkt.check_host_genus(set1, False)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_subcluster_1(self):
        """Check that no error is produced if the
        subcluster is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.subcluster = "Mycobacterium"
        self.tkt.check_subcluster(set1, False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_subcluster_2(self):
        """Check that an error is produced if the
        subcluster is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.subcluster = "none"
        self.tkt.check_subcluster(set1, False)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_cluster_1(self):
        """Check that no error is produced if the
        cluster is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.cluster = "Mycobacterium"
        self.tkt.check_cluster(set1, False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_cluster_2(self):
        """Check that an error is produced if the
        cluster is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.cluster = "none"
        self.tkt.check_cluster(set1, False)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_annotation_status_1(self):
        """Check that no error is produced if the
        annotation_status is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.annotation_status = "Mycobacterium"
        self.tkt.check_annotation_status(set1, False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_annotation_status_2(self):
        """Check that an error is produced if the
        annotation_status is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.annotation_status = "none"
        self.tkt.check_annotation_status(set1, False)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_description_field_1(self):
        """Check that no error is produced if the
        description_field is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.description_field = "Mycobacterium"
        self.tkt.check_description_field(set1, False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_description_field_2(self):
        """Check that an error is produced if the
        description_field is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.description_field = "none"
        self.tkt.check_description_field(set1, False)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_accession_1(self):
        """Check that no error is produced if the
        accession is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.accession = "Mycobacterium"
        self.tkt.check_accession(set1, False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_accession_2(self):
        """Check that an error is produced if the
        accession is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.accession = "none"
        self.tkt.check_accession(set1, False)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_annotation_author_1(self):
        """Check that no error is produced if the
        annotation_author is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.annotation_author = "Mycobacterium"
        self.tkt.check_annotation_author(set1, False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_annotation_author_2(self):
        """Check that an error is produced if the
        annotation_author is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.annotation_author = "none"
        self.tkt.check_annotation_author(set1, False)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_run_mode_1(self):
        """Check that no error is produced if the
        run_mode is not present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.run_mode = "Mycobacterium"
        self.tkt.check_run_mode(set1, False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_run_mode_2(self):
        """Check that an error is produced if the
        run_mode is present in the empty/null set
        and not expected to be in the set."""
        set1 = set(["none"])
        self.tkt.run_mode = "none"
        self.tkt.check_run_mode(set1, False)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_duplicate_phage_id_1(self):
        """Check that no error is produced if the
        phage_id is not present in the set of duplicated values."""
        dupe_set = set(["Trixie", "L5"])
        self.tkt.phage_id = "D29"
        self.tkt.check_duplicate_phage_id(dupe_set, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_duplicate_phage_id_2(self):
        """Check that an error is produced if the
        phage_id is present in the set of duplicated values."""
        dupe_set = set(["Trixie", "L5"])
        self.tkt.phage_id = "Trixie"
        self.tkt.check_duplicate_phage_id(dupe_set)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_duplicate_accession_1(self):
        """Check that no error is produced if the
        accession is not present in the set of duplicated values."""
        dupe_set = set(["ABC123", "EFG456"])
        self.tkt.accession = "XYZ789"
        self.tkt.check_duplicate_accession(dupe_set, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_duplicate_accession_2(self):
        """Check that an error is produced if the
        accession is present in the set of duplicated values."""
        dupe_set = set(["ABC123", "EFG456"])
        self.tkt.accession = "ABC123"
        self.tkt.check_duplicate_accession(dupe_set)
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)




    def test_check_compatible_type_and_annotation_status_1(self):
        """Verify that no error is produced with "add" type and "draft"
        annotation_status."""
        self.tkt.type = "add"
        self.tkt.annotation_status = "draft"
        self.tkt.check_compatible_type_and_annotation_status(
            eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].id, "eval_id")

    def test_check_compatible_type_and_annotation_status_2(self):
        """Verify that an error is produced with "add" type and "final"
        annotation_status."""
        self.tkt.type = "add"
        self.tkt.annotation_status = "final"
        self.tkt.check_compatible_type_and_annotation_status()
        with self.subTest():
            self.assertEqual(self.tkt.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.tkt.evaluations[0].id)

    def test_check_compatible_type_and_annotation_status_3(self):
        """Verify that no error is produced with "replace" type and "final"
        annotation_status."""
        self.tkt.type = "replace"
        self.tkt.annotation_status = "final"
        self.tkt.check_compatible_type_and_annotation_status()
        self.assertEqual(self.tkt.evaluations[0].status, "correct")

    def test_check_compatible_type_and_annotation_status_4(self):
        """Verify that no error is produced with "replace" type and "draft"
        annotation_status."""
        self.tkt.type = "replace"
        self.tkt.annotation_status = "draft"
        self.tkt.check_compatible_type_and_annotation_status()
        self.assertEqual(self.tkt.evaluations[0].status, "error")






if __name__ == '__main__':
    unittest.main()
