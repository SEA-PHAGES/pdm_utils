""" Unit tests for the CDS class."""

from classes import Cds
import unittest



class TestCdsFeatureClass(unittest.TestCase):


    def setUp(self):
        self.feature = Cds.CdsFeature()





    def test_set_strand_1(self):
        self.feature.set_strand("f", "fr_long")
        self.assertEqual(self.feature.strand, "forward")

    def test_set_strand_2(self):
        self.feature.set_strand("reverse", "fr_short")
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







    def test_check_amino_acids_1(self):
        """All amino acids in alphabet."""
        alphabet = set(["A","B","C"])
        self.feature.translation = "AB"
        self.feature.check_amino_acids(alphabet)
        self.assertEqual(len(self.feature.evaluations), 0)

    def test_check_amino_acids_2(self):
        """Some amino acids not in alphabet."""
        alphabet = set(["A","B","C"])
        self.feature.translation = "AD"
        self.feature.check_amino_acids(alphabet)
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

    def test_set_start_end_5(self):
        """Operator strand feature."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.strand = "+"
        start = 5
        end = 10
        self.feature.set_start_end()
        with self.subTest():
            self.assertEqual(self.feature.start, start)
        with self.subTest():
            self.assertEqual(self.feature.end, end)

    def test_set_start_end_6(self):
        """Numeric strand feature."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.strand = -1
        start = 10
        end = 5
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
        self.feature.start = 5
        self.feature.end = 10
        location_id_1 = (5, 10, "forward")
        location_id_2 = (10, "forward")
        location_id_3 = (5, 10)
        self.feature.set_location_id()
        with self.subTest():
            self.assertEqual(self.feature._left_right_strand_id, location_id_1)
        with self.subTest():
            self.assertEqual(self.feature._end_strand_id, location_id_2)
        with self.subTest():
            self.assertEqual(self.feature._start_end_id, location_id_3)

    def test_set_location_id_2(self):
        """Reverse strand feature, both values should be set."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.strand = "reverse"
        self.feature.start = 10
        self.feature.end = 5
        location_id_1 = (5, 10, "reverse")
        location_id_2 = (5, "reverse")
        location_id_3 = (10, 5)
        self.feature.set_location_id()
        with self.subTest():
            self.assertEqual(self.feature._left_right_strand_id, location_id_1)
        with self.subTest():
            self.assertEqual(self.feature._end_strand_id, location_id_2)
        with self.subTest():
            self.assertEqual(self.feature._start_end_id, location_id_3)

    def test_set_location_id_3(self):
        """Test forward strand numeric format."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.strand = 1
        self.feature.start = 5
        self.feature.end = 10
        location_id_1 = (5, 10, 1)
        location_id_2 = (10, 1)
        location_id_3 = (5, 10)
        self.feature.set_location_id()
        with self.subTest():
            self.assertEqual(self.feature._left_right_strand_id, location_id_1)
        with self.subTest():
            self.assertEqual(self.feature._end_strand_id, location_id_2)
        with self.subTest():
            self.assertEqual(self.feature._start_end_id, location_id_3)

    def test_set_location_id_4(self):
        """Test non-standard strand format."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.strand = "abcd"
        location_id_1 = (5, 10, "abcd")
        location_id_2 = ("", "abcd")
        location_id_3 = ("", "")
        self.feature.set_location_id()
        with self.subTest():
            self.assertEqual(self.feature._left_right_strand_id, location_id_1)
        with self.subTest():
            self.assertEqual(self.feature._end_strand_id, location_id_2)
        with self.subTest():
            self.assertEqual(self.feature._start_end_id, location_id_3)














    def test_check_boundaries_1(self):
        """Test correct boundaries."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = 10
        self.feature.check_boundaries()
        self.assertEqual(len(self.feature.evaluations), 0)

    def test_check_boundaries_2(self):
        """Test incorrect left boundary."""
        self.feature.left_boundary = "a"
        self.feature.right_boundary = 10
        self.feature.check_boundaries()
        self.assertEqual(len(self.feature.evaluations), 1)

    def test_check_boundaries_3(self):
        """Test incorrect right boundary."""
        self.feature.left_boundary = 5
        self.feature.right_boundary = "a"
        self.feature.check_boundaries()
        self.assertEqual(len(self.feature.evaluations), 1)





    def test_set_phage_id_1(self):
        """Test no draft suffix."""
        phage_id = "Trixie"
        self.feature.set_phage_id(phage_id)
        with self.subTest():
            self.assertEqual(self.feature.phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(self.feature._search_id, "Trixie")

    def test_set_phage_id_2(self):
        """Test draft suffix."""
        phage_id = "Trixie_Draft"
        self.feature.set_phage_id(phage_id)
        with self.subTest():
            self.assertEqual(self.feature.phage_id, "Trixie_Draft")
        with self.subTest():
            self.assertEqual(self.feature._search_id, "Trixie")


    #TODO: remove this unit test after I decide whether or not to keep
    # the set_primary_description function.
    # de test_set_primary_description_1(self):
    #     """Test primary description."""
    #     description1 = "ABCD"
    #     description2 = "EFGH"
    #     self.feature.set_primary_description(description1, description2)
    #     with self.subTest():
    #         self.assertEqual(self.feature.primary_description, "ABCD")
    #     with self.subTest():
    #         self.assertEqual(self.feature.processed_primary_description, "EFGH")




    def test_check_locus_tag_present_1(self):
        """Check if absent locus tag is expected to be absent."""
        self.feature.locus_tag = ""
        self.feature.check_locus_tag_present("absent")
        self.assertEqual(len(self.feature.evaluations), 0)

    def test_check_locus_tag_present_2(self):
        """Check if absent locus tag is expected to be present."""
        self.feature.locus_tag = ""
        self.feature.check_locus_tag_present("present")
        self.assertEqual(len(self.feature.evaluations), 1)

    def test_check_locus_tag_present_3(self):
        """Check if present locus tag is expected to be present."""
        self.feature.locus_tag = "ABCD"
        self.feature.check_locus_tag_present("present")
        self.assertEqual(len(self.feature.evaluations), 0)

    def test_check_locus_tag_present_4(self):
        """Check if present locus tag is expected to be absent."""
        self.feature.locus_tag = "ABCD"
        self.feature.check_locus_tag_present("absent")
        self.assertEqual(len(self.feature.evaluations), 1)






    def test_check_description_1(self):
        """Product is present and function is present."""
        self.feature.processed_product_description = "ABC"
        self.feature.processed_function_description = "EFG"
        self.feature.check_description()
        self.assertEqual(len(self.feature.evaluations), 0)

    def test_check_description_2(self):
        """Product is present and function is absent."""
        self.feature.processed_product_description = "ABC"
        self.feature.check_description()
        self.assertEqual(len(self.feature.evaluations), 0)

    def test_check_description_3(self):
        """Product is absent and function is present."""
        self.feature.processed_function_description = "EFG"
        self.feature.check_description()
        self.assertEqual(len(self.feature.evaluations), 1)





    def test_check_locus_tag_typo_1(self):
        """The locus_tag does not contain a typo."""
        self.feature.locus_tag = "ABC_TRIXIE_123"
        self.feature.check_locus_tag_typo("Trixie")
        self.assertEqual(len(self.feature.evaluations), 0)

    def test_check_locus_tag_typo_2(self):
        """The locus_tag contains a typo."""
        self.feature.locus_tag = "ABC_TRIXE_123"
        self.feature.check_locus_tag_typo("Trixie")
        self.assertEqual(len(self.feature.evaluations), 1)







    def test_check_lengths_1(self):
        """The translation length is correct."""
        self.feature.compound_parts = 1
        self.feature.left_boundary = 0
        self.feature.right_boundary = 11
        self.feature.set_translation("ABC")
        self.feature.check_lengths()
        self.assertEqual(len(self.feature.evaluations), 0)

    def test_check_lengths_2(self):
        """The translation length is not correct."""
        self.feature.compound_parts = 1
        self.feature.left_boundary = 0
        self.feature.right_boundary = 12
        self.feature.set_translation("ABC")
        self.feature.check_lengths()
        self.assertEqual(len(self.feature.evaluations), 1)

    def test_check_lengths_3(self):
        """Compound feature is not computed."""
        self.feature.compound_parts = 2
        self.feature.left_boundary = 0
        self.feature.right_boundary = 12
        self.feature.set_translation("ABC")
        self.feature.check_lengths()
        self.assertEqual(len(self.feature.evaluations), 0)




    def test_set_nucleotide_length_1(self):
        """Verify the nucleotide length is correct."""
        self.feature.left_boundary = 0
        self.feature.right_boundary = 11
        self.feature.set_nucleotide_length()
        exp = 12
        self.assertEqual(self.feature._nucleotide_length, 12)



    # TODO add more set_nucleotide_length unit tests.




    def test_check_translation_length_1(self):
        """Verify a present translation does not produce an error."""
        self.feature._translation_length = 1
        self.feature.check_translation_length()
        self.assertEqual(len(self.feature.evaluations), 0)

    def test_check_translation_length_2(self):
        """Verify a present translation does not produce an error."""
        self.feature._translation_length = 100
        self.feature.check_translation_length()
        self.assertEqual(len(self.feature.evaluations), 0)

    def test_check_translation_length_3(self):
        """Verify that no translation produces an error."""
        self.feature._translation_length = 0
        self.feature.check_translation_length()
        self.assertEqual(len(self.feature.evaluations), 1)








    def test_choose_description_1(self):
        """Verify product description is assigned to primary description."""
        self.feature.product_description = "ABCD"
        self.feature.processed_primary_description = "EFGH"
        self.feature.choose_description("product")
        with self.subTest():
            self.assertEqual(self.feature.primary_description, "ABCD")
        with self.subTest():
            self.assertEqual(self.feature.processed_primary_description, "EFGH")

    def test_choose_description_2(self):
        """Verify function description is assigned to primary description."""
        self.feature.function_description = "ABCD"
        self.feature.processed_function_description = "EFGH"
        self.feature.choose_description("function")
        with self.subTest():
            self.assertEqual(self.feature.primary_description, "ABCD")
        with self.subTest():
            self.assertEqual(self.feature.processed_primary_description, "EFGH")

    def test_choose_description_3(self):
        """Verify note description is assigned to primary description."""
        self.feature.note_description = "ABCD"
        self.feature.processed_note_description = "EFGH"
        self.feature.choose_description("note")
        with self.subTest():
            self.assertEqual(self.feature.primary_description, "ABCD")
        with self.subTest():
            self.assertEqual(self.feature.processed_primary_description, "EFGH")

    def test_choose_description_4(self):
        """Verify no description is assigned to primary description."""
        self.feature.note_description = "ABCD"
        self.feature.processed_note_description = "EFGH"
        self.feature.choose_description("invalid")
        with self.subTest():
            self.assertEqual(self.feature.primary_description, "")
        with self.subTest():
            self.assertEqual(self.feature.processed_primary_description, "")










if __name__ == '__main__':
    unittest.main()
