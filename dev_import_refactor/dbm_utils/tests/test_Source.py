""" Unit tests for the CDS class."""

from classes import Source
import unittest



class TestSourceFeatureClass(unittest.TestCase):


    def setUp(self):
        self.feature = Source.SourceFeature()









    def test_parse_organism_1(self):
        """Verify empty string is parsed correctly."""
        self.feature.organism = ""
        self.feature.parse_organism()
        expected_phage = ""
        expected_host = ""
        with self.subTest():
            self.assertEqual(\
                self.feature._organism_phage_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.feature._organism_host_name, expected_host)

    def test_parse_organism_2(self):
        """Verify string is parsed correctly."""
        self.feature.organism = "asdf Mycobacterium phage Trixie."
        self.feature.parse_organism()
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        with self.subTest():
            self.assertEqual(\
                self.feature._organism_phage_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.feature._organism_host_name, expected_host)




    def test_parse_host_1(self):
        """Verify empty string is parsed correctly."""
        self.feature.host = ""
        self.feature.parse_host()
        expected_host = ""
        self.assertEqual(\
            self.feature._host_host_name, expected_host)

    def test_parse_host_2(self):
        """Verify string is parsed correctly."""
        self.feature.host = "asdf Mycobacterium phage Trixie."
        self.feature.parse_host()
        expected_host = "Mycobacterium"
        self.assertEqual(\
            self.feature._host_host_name, expected_host)




    def test_parse_lab_host_1(self):
        """Verify empty string is parsed correctly."""
        self.feature.lab_host = ""
        self.feature.parse_lab_host()
        expected_host = ""
        self.assertEqual(\
            self.feature._lab_host_host_name, expected_host)

    def test_parse_lab_host_2(self):
        """Verify string is parsed correctly."""
        self.feature.lab_host = "asdf Mycobacterium phage Trixie."
        self.feature.parse_lab_host()
        expected_host = "Mycobacterium"
        self.assertEqual(\
            self.feature._lab_host_host_name, expected_host)





















if __name__ == '__main__':
    unittest.main()
