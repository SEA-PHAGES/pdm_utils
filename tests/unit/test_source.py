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
        self.assertEqual(self.feature._host_host_genus, expected_host_genus)

    def test_parse_host_2(self):
        """Verify string is parsed correctly."""
        self.feature.host = "Mycobacterium smegmatis."
        self.feature.parse_host()
        expected_host_genus = "Mycobacterium"
        self.assertEqual(self.feature._host_host_genus, expected_host_genus)



    def test_parse_lab_host_1(self):
        """Verify empty string is parsed correctly."""
        self.feature.lab_host = ""
        self.feature.parse_lab_host()
        expected_host_genus = ""
        self.assertEqual(\
            self.feature._lab_host_host_genus, expected_host_genus)

    def test_parse_lab_host_2(self):
        """Verify string is parsed correctly."""
        self.feature.lab_host = "Mycobacterium Trixie."
        self.feature.parse_lab_host()
        expected_host_genus = "Mycobacterium"
        self.assertEqual(\
            self.feature._lab_host_host_genus, expected_host_genus)




    def test_check_attribute_1(self):
        """Verify no error is produced when the host_genus
        is in the check_set and is expected to be in the set."""
        check_set = set(["Mycobacterium", "Mycobacterio"])
        self.feature._organism_host_genus = "Mycobacterium"
        self.feature.check_attribute("_organism_host_genus", check_set, True, "eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_attribute_2(self):
        """Verify an error is produced when the host_genus
        is not in the check_set and is expected to be in the set."""
        check_set = set(["Mycobacterium", "Mycobacterio"])
        self.feature._organism_host_genus = "Gordonia"
        self.feature.check_attribute("_organism_host_genus", check_set, True)
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_attribute_3(self):
        """Verify no error is produced when the host_genus
        is not in the check_set and is not expected to be in the set."""
        check_set = set(["Mycobacterium", "Mycobacterio"])
        self.feature._organism_host_genus = "Gordonia"
        self.feature.check_attribute("_organism_host_genus", check_set, False)
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_attribute_4(self):
        """Verify an error is produced when the host_genus
        is in the check_set and is not expected to be in the set."""
        check_set = set(["Mycobacterium", "Mycobacterio"])
        self.feature._organism_host_genus = "Mycobacterium"
        self.feature.check_attribute("_organism_host_genus", check_set, False)
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_attribute_5(self):
        """Verify no test is performed when the attribute is invalid."""
        check_set = set([1, 0])
        self.feature.check_attribute("invalid", check_set, True)
        self.assertEqual(self.feature.evaluations[0].status, "untested")



if __name__ == '__main__':
    unittest.main()
