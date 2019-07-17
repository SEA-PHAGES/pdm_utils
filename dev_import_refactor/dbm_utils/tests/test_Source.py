""" Unit tests for the CDS class."""

from classes import Source
import unittest



class TestSourceFeatureClass(unittest.TestCase):


    def setUp(self):
        self.feature = Source.SourceFeature()






    # TODO this is probably no longer needed.
    # def test_set_evaluation_1(self):
    #     """Set an empty evaluation object."""
    #     self.feature.set_evaluation("none")
    #     self.assertEqual(len(self.feature.evaluations), 1)
    #
    # def test_set_evaluation_2(self):
    #     """Set a warning evaluation object."""
    #     self.feature.set_evaluation("warning","message1")
    #     self.assertEqual(len(self.feature.evaluations), 1)
    #
    # def test_set_evaluation_3(self):
    #     """Set an error evaluation object."""
    #     self.feature.set_evaluation("error","message1","message2")
    #     self.assertEqual(len(self.feature.evaluations), 1)







    def test_parse_organism_1(self):
        """Verify empty string is parsed correctly."""
        self.feature.organism = ""
        self.feature.parse_organism()
        expected_phage = ""
        expected_host_genus = ""
        with self.subTest():
            self.assertEqual(\
                self.feature._organism_phage_name, expected_phage)
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
                self.feature._organism_phage_name, expected_phage)
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




    def test_check_organism_phage_name_1(self):
        """Check that no warning is produced."""
        self.feature.parent_phage_id = "Trixie"
        self.feature._organism_phage_name = "Trixie"
        self.feature.check_organism_phage_name()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_organism_phage_name_2(self):
        """Check that a warning is produced."""
        self.feature.parent_phage_id = "L5"
        self.feature._organism_phage_name = "Trixie"
        self.feature.check_organism_phage_name()
        self.assertEqual(self.feature.evaluations[0].status, "error")




    def test_check_organism_host_genus_1(self):
        """Check that no warning is produced."""
        self.feature.parent_host_genus = "Mycobacterium"
        self.feature._organism_host_genus = "Mycobacterium"
        self.feature.check_organism_host_genus()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_organism_host_genus_2(self):
        """Check that a warning is produced."""
        self.feature.parent_host_genus = "Gordonia"
        self.feature._organism_host_genus = "Mycobacterium"
        self.feature.check_organism_host_genus()
        self.assertEqual(self.feature.evaluations[0].status, "error")




    def test_check_host_host_genus_1(self):
        """Check that no warning is produced."""
        self.feature.parent_host_genus = "Mycobacterium"
        self.feature._host_host_genus = "Mycobacterium"
        self.feature.check_host_host_genus()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_host_host_genus_2(self):
        """Check that a warning is produced."""
        self.feature.parent_host_genus = "Gordonia"
        self.feature._host_host_genus = "Mycobacterium"
        self.feature.check_host_host_genus()
        self.assertEqual(self.feature.evaluations[0].status, "error")




    def test_check_lab_host_host_genus_1(self):
        """Check that no warning is produced."""
        self.feature.parent_host_genus = "Mycobacterium"
        self.feature._lab_host_host_genus = "Mycobacterium"
        self.feature.check_lab_host_host_genus()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_lab_host_host_genus_2(self):
        """Check that a warning is produced."""
        self.feature.parent_host_genus = "Gordonia"
        self.feature._lab_host_host_genus = "Mycobacterium"
        self.feature.check_lab_host_host_genus()
        self.assertEqual(self.feature.evaluations[0].status, "error")
























if __name__ == '__main__':
    unittest.main()
