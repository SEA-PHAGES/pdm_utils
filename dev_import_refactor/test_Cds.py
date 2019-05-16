""" Unit tests for the CDS class."""

import Cds
import unittest



class TestCdsFeatureClass(unittest.TestCase):


    def setUp(self):
        self.feature = Cds.CdsFeature()





    def test_set_strand_1(self):
        self.feature.set_strand("f", "long")
        self.assertEqual(self.feature.strand, "forward")

    def test_set_strand_2(self):
        self.feature.set_strand("reverse", "short")
        self.assertEqual(self.feature.strand, "r")


    def test_set_translation_1(self):
        self.feature.set_translation("abcd")
        with self.subTest():
            self.assertEqual(self.feature.translation, "ABCD")
        with self.subTest():
            self.assertEqual(self.feature._translation_length, 4)






    def test_set_evaluation_1(self):
        self.feature.set_evaluation("none")
        self.assertEqual(len(self.feature.evaluations), 1)

    def test_set_evaluation_2(self):
        self.feature.set_evaluation("warning","message1")
        self.assertEqual(len(self.feature.evaluations), 1)

    def test_set_evaluation_3(self):
        self.feature.set_evaluation("error","message1","message2")
        self.assertEqual(len(self.feature.evaluations), 1)







    def test_check_translation_1(self):
        """All amino acids in alphabet."""
        alphabet = set(["A","B","C"])
        self.feature.translation = "AB"
        self.feature.check_translation(alphabet)
        self.assertEqual(len(self.feature.evaluations), 0)

    def test_check_translation_2(self):
        """Some amino acids not in alphabet."""
        alphabet = set(["A","B","C"])
        self.feature.translation = "AD"
        self.feature.check_translation(alphabet)
        self.assertEqual(len(self.feature.evaluations), 1)












    def test_set_start_end_1(self):
        """Forward strand feature, long format."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.strand = "forward"
        start = 5
        end = 10
        self.feature.set_start_end()
        with self.subTest():
            self.assertEqual(self.feature.start, start)
        with self.subTest():
            self.assertEqual(self.feature.end, end)

    def test_set_start_end_2(self):
        """Reverse strand feature, long format."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.strand = "reverse"
        start = 10
        end = 5
        self.feature.set_start_end()
        with self.subTest():
            self.assertEqual(self.feature.start, start)
        with self.subTest():
            self.assertEqual(self.feature.end, end)

    def test_set_start_end_3(self):
        """Reverse strand feature, short format."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.strand = "r"
        start = 10
        end = 5
        self.feature.set_start_end()
        with self.subTest():
            self.assertEqual(self.feature.start, start)
        with self.subTest():
            self.assertEqual(self.feature.end, end)

    def test_set_start_end_4(self):
        """Non-standard strand feature."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.strand = "other"
        start = ""
        end = ""
        self.feature.set_start_end()
        with self.subTest():
            self.assertEqual(self.feature.start, start)
        with self.subTest():
            self.assertEqual(self.feature.end, end)




    def test_set_location_id_1(self):
        """Forward strand feature, both values should be set."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.strand = "forward"
        self.feature.end = 10
        location_id_1 = (5, 10, "forward")
        location_id_2 = (10, "forward")
        self.feature.set_location_id()
        with self.subTest():
            self.assertEqual(self.feature._left_right_strand_id, location_id_1)
        with self.subTest():
            self.assertEqual(self.feature._end_strand_id, location_id_2)

    def test_set_location_id_2(self):
        """Reverse strand feature, both values should be set."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.strand = "reverse"
        self.feature.end = 5
        location_id_1 = (5, 10, "reverse")
        location_id_2 = (5, "reverse")
        self.feature.set_location_id()
        with self.subTest():
            self.assertEqual(self.feature._left_right_strand_id, location_id_1)
        with self.subTest():
            self.assertEqual(self.feature._end_strand_id, location_id_2)

    def test_set_location_id_3(self):
        """Test forward strand numeric format."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.strand = 1
        self.feature.end = 10
        location_id_1 = (5, 10, 1)
        location_id_2 = (10, 1)
        self.feature.set_location_id()
        with self.subTest():
            self.assertEqual(self.feature._left_right_strand_id, location_id_1)
        with self.subTest():
            self.assertEqual(self.feature._end_strand_id, location_id_2)

    def test_set_location_id_4(self):
        """Test non-standard strand format."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.strand = "abcd"
        location_id_1 = (5, 10, "abcd")
        location_id_2 = ("","abcd")
        self.feature.set_location_id()
        with self.subTest():
            self.assertEqual(self.feature._left_right_strand_id, location_id_1)
        with self.subTest():
            self.assertEqual(self.feature._end_strand_id, location_id_2)





















if __name__ == '__main__':
    unittest.main()
