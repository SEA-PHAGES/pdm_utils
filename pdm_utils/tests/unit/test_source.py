""" Unit tests for the CDS class."""

from pdm_utils.classes import source
import unittest



class TestSourceClass(unittest.TestCase):


    def setUp(self):
        self.feature = source.Source()




    def test_parse_organism_1(self):
        """Verify empty string is parsed correctly."""
        self.feature.organism = ""
        self.feature.parse_organism()
        expected_phage = ""
        expected_host_genus = ""
        with self.subTest():
            self.assertEqual(\
                self.feature._organism_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.feature._organism_host_genus, expected_host_genus)

    def test_parse_organism_2(self):
        """Verify string is parsed correctly."""
        self.feature.organism = "asdf Mycobacterium phage Trixie."
        self.feature.parse_organism()
        expected_phage = "Trixie"
        expected_host_genus = "Mycobacterium"
        with self.subTest():
            self.assertEqual(\
                self.feature._organism_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.feature._organism_host_genus, expected_host_genus)




    def test_parse_host_1(self):
        """Verify empty string is parsed correctly."""
        self.feature.host = ""
        self.feature.parse_host()
        expected_host_genus = ""
        self.assertEqual(\
            self.feature._host_host_genus, expected_host_genus)

    def test_parse_host_2(self):
        """Verify string is parsed correctly."""
        self.feature.host = "asdf Mycobacterium phage Trixie."
        self.feature.parse_host()
        expected_host_genus = "Mycobacterium"
        self.assertEqual(\
            self.feature._host_host_genus, expected_host_genus)




    def test_parse_lab_host_1(self):
        """Verify empty string is parsed correctly."""
        self.feature.lab_host = ""
        self.feature.parse_lab_host()
        expected_host_genus = ""
        self.assertEqual(\
            self.feature._lab_host_host_genus, expected_host_genus)

    def test_parse_lab_host_2(self):
        """Verify string is parsed correctly."""
        self.feature.lab_host = "asdf Mycobacterium phage Trixie."
        self.feature.parse_lab_host()
        expected_host_genus = "Mycobacterium"
        self.assertEqual(\
            self.feature._lab_host_host_genus, expected_host_genus)




    def test_check_organism_name_1(self):
        """Check that no warning is produced."""
        self.feature.genome_id = "Trixie"
        self.feature._organism_name = "Trixie"
        self.feature.check_organism_name(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_organism_name_2(self):
        """Check that a warning is produced."""
        self.feature.genome_id = "L5"
        self.feature._organism_name = "Trixie"
        self.feature.check_organism_name()
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)




    def test_check_organism_host_genus_1(self):
        """Check that no warning is produced."""
        self.feature.parent_host_genus = "Mycobacterium"
        self.feature._organism_host_genus = "Mycobacterium"
        self.feature.check_organism_host_genus(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_organism_host_genus_2(self):
        """Check that a warning is produced."""
        self.feature.parent_host_genus = "Gordonia"
        self.feature._organism_host_genus = "Mycobacterium"
        self.feature.check_organism_host_genus()
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)




    def test_check_host_host_genus_1(self):
        """Check that no warning is produced."""
        self.feature.parent_host_genus = "Mycobacterium"
        self.feature._host_host_genus = "Mycobacterium"
        self.feature.check_host_host_genus(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_host_host_genus_2(self):
        """Check that a warning is produced."""
        self.feature.parent_host_genus = "Gordonia"
        self.feature._host_host_genus = "Mycobacterium"
        self.feature.check_host_host_genus()
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)




    def test_check_lab_host_host_genus_1(self):
        """Check that no warning is produced."""
        self.feature.parent_host_genus = "Mycobacterium"
        self.feature._lab_host_host_genus = "Mycobacterium"
        self.feature.check_lab_host_host_genus(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_lab_host_host_genus_2(self):
        """Check that a warning is produced."""
        self.feature.parent_host_genus = "Gordonia"
        self.feature._lab_host_host_genus = "Mycobacterium"
        self.feature.check_lab_host_host_genus()
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)
























if __name__ == '__main__':
    unittest.main()
