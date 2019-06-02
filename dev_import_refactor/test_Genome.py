""" Unit tests for the Genome class."""


import unittest
import Genome



class TestGenomeClass(unittest.TestCase):


    def setUp(self):
        self.genome = Genome.Genome()



    def test_set_evaluation_1(self):
        """Set an empty evaluation object."""
        self.genome.set_evaluation("none")
        self.assertEqual(len(self.genome.evaluations), 1)

    def test_set_evaluation_2(self):
        """Set a warning evaluation object."""
        self.genome.set_evaluation("warning","message1")
        self.assertEqual(len(self.genome.evaluations), 1)

    def test_set_evaluation_3(self):
        """Set an error evaluation object."""
        self.genome.set_evaluation("error","message1","message2")
        self.assertEqual(len(self.genome.evaluations), 1)





    def test_set_filename_1(self):
        """Confirm file path is split appropriately."""
        filepath = "/path/to/folder/Trixie.gbk"
        self.genome.set_filename(filepath)
        with self.subTest():
            self.assertEqual(self.genome.filename, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome.search_filename, "trixie")





    def test_set_host_1(self):
        """Check that host name is split appropriately."""
        host = "Mycobacterium smegmatis"
        self.genome.set_host(host)
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_2(self):
        """Check that whitespace is removed."""
        host = "  Mycobacterium smegmatis  "
        self.genome.set_host(host)
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_3(self):
        """Check that none is set appropriately."""
        host = "  none  "
        self.genome.set_host(host)
        self.assertEqual(self.genome.host, "")



    def test_set_sequence_1(self):
        """Check that sequence is set appropriately."""
        seq = "abcd"
        self.genome.set_sequence(seq)
        with self.subTest():
            self.assertEqual(self.genome.sequence, "ABCD")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)








    def test_set_accession_1(self):
        """Check that accession is set appropriately."""
        accession = "ABC123.1"
        self.genome.set_accession(accession,"filled")
        self.assertEqual(self.genome.accession, "ABC123")

    def test_set_accession_2(self):
        """Check that empty accession is set appropriately from None."""
        accession = None
        self.genome.set_accession(accession,"empty")
        self.assertEqual(self.genome.accession, "")

    def test_set_accession_3(self):
        """Check that filled accession is set appropriately from None."""
        accession = None
        self.genome.set_accession(accession,"filled")
        self.assertEqual(self.genome.accession, "none")

    def test_set_accession_4(self):
        """Check that null accession is set appropriately from None."""
        accession = None
        self.genome.set_accession(accession,"null")
        self.assertIsNone(self.genome.accession)

    def test_set_accession_5(self):
        """Check that null accession is set appropriately from empty."""
        accession = ""
        self.genome.set_accession(accession,"null")
        self.assertIsNone(self.genome.accession)





    def test_check_nucleotides_1(self):
        """All nucleotides are in the alphabet."""
        alphabet = set(["A","B","C"])
        self.genome.sequence = "AB"
        self.genome.check_nucleotides(alphabet)
        self.assertEqual(len(self.genome.evaluations), 0)

    def test_check_nucleotides_2(self):
        """Some nucleotides are not in the alphabet."""
        alphabet = set(["A","B","C"])
        self.genome.sequence = "AD"
        self.genome.check_nucleotides(alphabet)
        self.assertEqual(len(self.genome.evaluations), 1)











    def test_set_cds_features_1(self):
        """Check that feature list is set and length is computed."""
        features_list = [0,1,2,3]
        self.genome.set_cds_features(features_list)
        self.assertEqual(len(self.genome._cds_features), 4)




    def test_set_phage_id_1(self):
        """Check that name without '_Draft' suffix is not changed."""
        phage_name = "Trixie"
        self.genome.set_phage_id(phage_name)
        with self.subTest():
            self.assertEqual(self.genome.phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome.search_id, "trixie")

    def test_set_phage_id_2(self):
        """Check that '_Draft' suffix is removed."""
        phage_name = "Trixie_Draft"
        self.genome.set_phage_id(phage_name)
        with self.subTest():
            self.assertEqual(self.genome.phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome.search_id, "trixie")




if __name__ == '__main__':
    unittest.main()
