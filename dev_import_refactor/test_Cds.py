""" Unit tests for the CDS class."""

import Cds
import unittest



class TestCdsFeatureClass(unittest.TestCase):


    def setUp(self):
        self.feature = Cds.CdsFeature()





    def test_set_strand_1(self):
        self.feature.set_strand("f")
        self.assertEqual(self.feature.strand, "forward")


    def test_set_translation_1(self):
        self.feature.set_translation("abcd")
        with self.subTest():
            self.assertEqual(self.feature.translation, "ABCD")
        with self.subTest():
            self.assertEqual(self.feature._translation_length, 4)



if __name__ == '__main__':
    unittest.main()
