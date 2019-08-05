""" Unit tests for the CDS class."""

from classes import Cds
from constants import constants
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import unittest



class TestCdsClass(unittest.TestCase):


    def setUp(self):
        self.feature = Cds.Cds()




    def test_set_translation_table_1(self):
        """Verify translation_table is set correctly from valid integer."""
        self.feature.set_translation_table(11)
        self.assertEqual(self.feature.translation_table, 11)

    def test_set_translation_table_2(self):
        """Verify translation_table is set correctly from valid string."""
        self.feature.set_translation_table("11")
        self.assertEqual(self.feature.translation_table, 11)

    def test_set_translation_table_2(self):
        """Verify translation_table is set correctly from invalid string."""
        self.feature.set_translation_table("x")
        self.assertEqual(self.feature.translation_table, -1)




    def test_set_strand_1(self):
        """Verify 'f' is converted correctly."""
        self.feature.set_strand("f", "fr_long")
        self.assertEqual(self.feature.strand, "forward")

    def test_set_strand_2(self):
        """Verify 'r' is converted correctly."""
        self.feature.set_strand("reverse", "fr_short")
        self.assertEqual(self.feature.strand, "r")




    def test_translate_seq_1(self):
        """Verify translation is produced from valid nucleotide sequence."""
        self.feature.seq = Seq("ATGTTTTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        translation = self.feature.translate_seq()
        self.assertEqual(translation, "MF")

    def test_translate_seq_2(self):
        """Verify translation is produced from valid nucleotide sequence
        that contains a non-standard start codon."""
        self.feature.seq = Seq("GTGTTTTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        translation = self.feature.translate_seq()
        self.assertEqual(translation, "MF")

    def test_translate_seq_3(self):
        """Verify that sequence and length is set from translating
        invalid nucleotide sequence (contains one additional codon
        after stop codon)."""
        self.feature.seq = Seq("GTGTTTTGAATG", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        translation = self.feature.translate_seq()
        self.assertEqual(translation, "")

    def test_translate_seq_4(self):
        """Verify that sequence and length is set from translating
        invalid nucleotide sequence (contains ambiguous nucleotides)."""
        self.feature.seq = Seq("GTGRRRTTTTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        translation = self.feature.translate_seq()
        self.assertEqual(translation, "")

    def test_translate_seq_5(self):
        """Verify that sequence and length is set from translating
        invalid nucleotide sequence (contains an extra nucleotide)."""
        self.feature.seq = Seq("GTGATTTTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        translation = self.feature.translate_seq()
        self.assertEqual(translation, "")

    def test_translate_seq_6(self):
        """Verify that sequence and length is set from translating
        invalid nucleotide sequence (contains an invalid start codon)."""
        self.feature.seq = Seq("GTATTTTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        translation = self.feature.translate_seq()
        self.assertEqual(translation, "")




    def test_set_translation_1(self):
        """Verify that sequence and length is set from valid Seq object."""
        self.feature.set_translation(Seq("mf"))
        with self.subTest():
            self.assertEqual(self.feature.translation, "MF")
        with self.subTest():
            self.assertEqual(self.feature._translation_length, 2)

    def test_set_translation_2(self):
        """Verify that sequence and length is set from valid input."""
        self.feature.set_translation("mf")
        with self.subTest():
            self.assertEqual(self.feature.translation, "MF")
        with self.subTest():
            self.assertEqual(self.feature._translation_length, 2)

    def test_set_translation_3(self):
        """Verify that sequence and length is set from invalid input."""
        self.feature.set_translation(1)
        with self.subTest():
            self.assertEqual(self.feature.translation, "")
        with self.subTest():
            self.assertEqual(self.feature._translation_length, 0)

    def test_set_translation_4(self):
        """Verify that sequence and length is set from translating
        valid nucleotide sequence."""
        self.feature.seq = Seq("ATGTTTTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        self.feature.set_translation(translate=True)
        with self.subTest():
            self.assertEqual(self.feature.translation, "MF")
        with self.subTest():
            self.assertEqual(self.feature._translation_length, 2)

    def test_set_translation_5(self):
        """Verify that sequence and length is set from translating
        invalid nucleotide sequence (contains one additional codon
        after stop codon)."""
        self.feature.seq = Seq("GTGTTTTGAATG", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        self.feature.set_translation(translate=True)
        with self.subTest():
            self.assertEqual(self.feature.translation, "")
        with self.subTest():
            self.assertEqual(self.feature._translation_length, 0)

    def test_set_translation_6(self):
        """Verify that sequence and length is set when no options are
        selected."""
        self.feature.set_translation()
        with self.subTest():
            self.assertEqual(self.feature.translation, "")
        with self.subTest():
            self.assertEqual(self.feature._translation_length, 0)




    def test_check_amino_acids_1(self):
        """Verify no error is produced if all amino acids
        are in the protein alphabet."""
        self.feature.translation = "ADE"
        self.feature.check_amino_acids()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_amino_acids_2(self):
        """Verify an error is produced if some amino acids
        are not in the protein alphabet."""
        self.feature.translation = "ABDE"
        self.feature.check_amino_acids()
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_amino_acids_3(self):
        """Verify no error is produced if all amino acids
        are in a custom protein alphabet."""
        alphabet = set(["A","B","C", "D", "E", "F"])
        self.feature.translation = "ABDE"
        self.feature.check_amino_acids(alphabet)
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_amino_acids_4(self):
        """Verify an error is produced if some amino acids
        are not in a custom protein alphabet."""
        alphabet = set(["A","B","C", "D", "E", "F"])
        self.feature.translation = "ABDEG"
        self.feature.check_amino_acids(alphabet)
        self.assertEqual(self.feature.evaluations[0].status, "error")




    def test_set_start_end_1(self):
        """Forward strand feature, long format."""
        self.feature.left = 5
        self.feature.right = 10
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
        self.feature.left = 5
        self.feature.right = 10
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
        self.feature.left = 5
        self.feature.right = 10
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
        self.feature.left = 5
        self.feature.right = 10
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
        self.feature.left = 5
        self.feature.right = 10
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
        self.feature.left = 5
        self.feature.right = 10
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
        self.feature.left = 5
        self.feature.right = 10
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
        self.feature.left = 5
        self.feature.right = 10
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
        self.feature.left = 5
        self.feature.right = 10
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
        self.feature.left = 5
        self.feature.right = 10
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




    def test_check_strand_1(self):
        """Verify no error is produced when the strand is
        formatted correctly using default settings."""
        self.feature.strand = "F"
        self.feature.check_strand()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_strand_2(self):
        """Verify an error is produced when the strand is
        formatted incorrectly using default settings."""
        self.feature.strand = 1
        self.feature.check_strand()
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_strand_3(self):
        """Verify no error is produced when the strand is
        formatted correctly using custom settings."""
        self.feature.strand = 1
        self.feature.check_strand(format="numeric", case=False)
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_strand_4(self):
        """Verify an error is produced when the strand is
        formatted incorrectly using custom settings."""
        self.feature.strand = "F"
        self.feature.check_strand(format="numeric", case=False)
        self.assertEqual(self.feature.evaluations[0].status, "error")




    def test_check_boundaries_1(self):
        """Test correct boundaries."""
        self.feature.left = 5
        self.feature.right = 10
        self.feature.check_boundaries()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_boundaries_2(self):
        """Test incorrect left boundary."""
        self.feature.left = "a"
        self.feature.right = 10
        self.feature.check_boundaries()
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_boundaries_3(self):
        """Test incorrect right boundary."""
        self.feature.left = 5
        self.feature.right = "a"
        self.feature.check_boundaries()
        self.assertEqual(self.feature.evaluations[0].status, "error")


    #TODO: remove this unit test after I decide whether or not to keep
    # the set_description function.
    # de test_set_description_1(self):
    #     """Test primary description."""
    #     description1 = "ABCD"
    #     description2 = "EFGH"
    #     self.feature.set_description(description1, description2)
    #     with self.subTest():
    #         self.assertEqual(self.feature.description, "ABCD")
    #     with self.subTest():
    #         self.assertEqual(self.feature.processed_description, "EFGH")




    def test_check_locus_tag_present_1(self):
        """Check if absent locus tag is expected to be absent."""
        self.feature.locus_tag = ""
        self.feature.check_locus_tag_present(False)
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_locus_tag_present_2(self):
        """Check if absent locus tag is expected to be present."""
        self.feature.locus_tag = ""
        self.feature.check_locus_tag_present(True)
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_locus_tag_present_3(self):
        """Check if present locus tag is expected to be present."""
        self.feature.locus_tag = "ABCD"
        self.feature.check_locus_tag_present(True)
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_locus_tag_present_4(self):
        """Check if present locus tag is expected to be absent."""
        self.feature.locus_tag = "ABCD"
        self.feature.check_locus_tag_present(False)
        self.assertEqual(self.feature.evaluations[0].status, "error")




    def test_check_description_1(self):
        """Verify no error is produced when a description is present in
        the processed_product as expected and the processed_function and
        the processed_note are empty."""
        self.feature.processed_product = "ABC"
        self.feature.check_description()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_description_2(self):
        """Verify no error is produced when a description is present in
        the processed_function as expected and the processed_product and
        the processed_note are empty."""
        self.feature.processed_function = "ABC"
        self.feature.check_description(description_field="function")
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_description_3(self):
        """Verify no error is produced when a description is present in
        the processed_note as expected and the processed_product and
        the processed_function are empty."""
        self.feature.processed_note = "ABC"
        self.feature.check_description(description_field="note")
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_description_4(self):
        """Verify an error is produced when a description is not present in
        the processed_product but there is a description in
        the processed_function."""
        self.feature.processed_function = "ABC"
        self.feature.check_description()
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_description_5(self):
        """Verify an error is produced when a description is not present in
        the processed_function but there is a description in
        the processed_product."""
        self.feature.processed_product = "ABC"
        self.feature.check_description(description_field="function")
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_description_6(self):
        """Verify an error is produced when a description is not present in
        the processed_function but there is a description in
        the processed_product."""
        self.feature.processed_product = "ABC"
        self.feature.check_description(description_field="note")
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_description_7(self):
        """Verify an error is produced when a description is not present due
        to an invalid input but there is a description in
        the processed_product."""
        self.feature.processed_product = "ABC"
        self.feature.check_description(description_field="invalid")
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_description_8(self):
        """Verify no error is produced when a description is present in
        the processed_product as expected as well as the processed_function."""
        self.feature.processed_product = "ABC"
        self.feature.processed_function = "FGH"
        self.feature.check_description()
        self.assertEqual(self.feature.evaluations[0].status, "correct")




    def test_check_locus_tag_typo_1(self):
        """The locus_tag does not contain a typo."""
        self.feature.locus_tag = "ABC_TRIXIE_123"
        self.feature.check_locus_tag_typo("Trixie")
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_locus_tag_typo_2(self):
        """The locus_tag contains a typo."""
        self.feature.locus_tag = "ABC_TRIXE_123"
        self.feature.check_locus_tag_typo("Trixie")
        self.assertEqual(self.feature.evaluations[0].status, "error")







    def test_check_lengths_1(self):
        """The translation length is correct."""
        self.feature.compound_parts = 1
        self.feature.left = 0
        self.feature.right = 11
        self.feature.set_translation("ABC")
        self.feature.check_lengths()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_lengths_2(self):
        """The translation length is not correct."""
        self.feature.compound_parts = 1
        self.feature.left = 0
        self.feature.right = 12
        self.feature.set_translation("ABC")
        self.feature.check_lengths()
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_lengths_3(self):
        """Compound feature is not computed."""
        self.feature.compound_parts = 2
        self.feature.left = 0
        self.feature.right = 12
        self.feature.set_translation("ABC")
        self.feature.check_lengths()
        self.assertEqual(self.feature.evaluations[0].status, "untested")




    def test_set_nucleotide_length_1(self):
        """Verify the nucleotide length is correct for a 0-based
        half-open interval."""
        self.feature.left = 0
        self.feature.right = 11
        self.feature.coordinate_format = "0_half_open"
        self.feature.set_nucleotide_length()
        self.assertEqual(self.feature._length, 11)

    def test_set_nucleotide_length_2(self):
        """Verify the nucleotide length is correct for a 1-based
        closed interval."""
        self.feature.left = 0
        self.feature.right = 11
        self.feature.coordinate_format = "1_closed"
        self.feature.set_nucleotide_length()
        self.assertEqual(self.feature._length, 12)

    def test_set_nucleotide_length_3(self):
        """Verify the nucleotide length is not set for invalid
        coordinate format."""
        self.feature.left = 0
        self.feature.right = 11
        self.feature.coordinate_format = "invalid"
        self.feature.set_nucleotide_length()
        self.assertEqual(self.feature._length, -1)







    def test_reformat_left_and_right_boundaries_1(self):
        """Verify the coordinates are converted to 1-based closed interval."""
        self.feature.left = 5
        self.feature.right = 11
        self.feature.coordinate_format = "0_half_open"
        new_format = "1_closed"
        self.feature.reformat_left_and_right_boundaries(new_format)
        with self.subTest():
            self.assertEqual(self.feature.left, 6)
        with self.subTest():
            self.assertEqual(self.feature.right, 11)
        with self.subTest():
            self.assertEqual(self.feature.coordinate_format, new_format)

    def test_reformat_left_and_right_boundaries_2(self):
        """Verify the coordinates are converted to 0-based half open interval."""
        self.feature.left = 5
        self.feature.right = 11
        self.feature.coordinate_format = "1_closed"
        new_format = "0_half_open"
        self.feature.reformat_left_and_right_boundaries(new_format)
        with self.subTest():
            self.assertEqual(self.feature.left, 4)
        with self.subTest():
            self.assertEqual(self.feature.right, 11)
        with self.subTest():
            self.assertEqual(self.feature.coordinate_format, new_format)

    def test_reformat_left_and_right_boundaries_3(self):
        """Verify the coordinates are not converted."""
        self.feature.left = 5
        self.feature.right = 11
        self.feature.coordinate_format = "1_closed"
        new_format = "invalid"
        self.feature.reformat_left_and_right_boundaries(new_format)
        with self.subTest():
            self.assertEqual(self.feature.left, 5)
        with self.subTest():
            self.assertEqual(self.feature.right, 11)
        with self.subTest():
            self.assertEqual(self.feature.coordinate_format, "1_closed")




    def test_set_nucleotide_sequence_1(self):
        """Verify that supplied Seq object is set correctly."""
        seq = Seq("aattcg")
        self.feature.set_nucleotide_sequence(value=seq)
        with self.subTest():
            self.assertEqual(self.feature.seq, "AATTCG")
        with self.subTest():
            self.assertIsInstance(self.feature.seq, Seq)

    def test_set_nucleotide_sequence_2(self):
        """Verify that supplied sequence is set correclty
        (converted to Seq object)."""
        seq = "aattcg"
        self.feature.set_nucleotide_sequence(value=seq)
        with self.subTest():
            self.assertEqual(self.feature.seq, "AATTCG")
        with self.subTest():
            self.assertIsInstance(self.feature.seq, Seq)

    def test_set_nucleotide_sequence_3(self):
        """Verify that supplied invalid sequence is set correctly
        (converted to empty Seq object)."""
        seq = 1
        self.feature.set_nucleotide_sequence(value=seq)
        with self.subTest():
            self.assertEqual(self.feature.seq, "")
        with self.subTest():
            self.assertIsInstance(self.feature.seq, Seq)

    def test_set_nucleotide_sequence_4(self):
        """Verify that expected sequence is extracted from top strand."""
        seq = Seq("AATTCG")
        self.feature.seqfeature = SeqFeature(FeatureLocation(1, 5),
                                             type="CDS",
                                             strand=1)
        self.feature.set_nucleotide_sequence(parent_genome_seq=seq)
        expected_seq = Seq("ATTC")
        self.assertEqual(self.feature.seq, expected_seq)

    def test_set_nucleotide_sequence_5(self):
        """Verify that expected sequence is extracted from bottom strand."""
        seq = Seq("AATTCG")
        self.feature.seqfeature = SeqFeature(FeatureLocation(1, 5),
                                             type="CDS",
                                             strand=-1)
        self.feature.set_nucleotide_sequence(parent_genome_seq=seq)
        expected_seq = Seq("GAAT")
        self.assertEqual(self.feature.seq, expected_seq)

    def test_set_nucleotide_sequence_6(self):
        """Verify that no sequence is extracted if the 'seqfeature'
        attribute is not a Biopython SeqFeature object."""
        seq = Seq("AATTCG")
        self.feature.seqfeature = ""
        self.feature.set_nucleotide_sequence(parent_genome_seq=seq)
        expected_seq = Seq("")
        self.assertEqual(self.feature.seq, expected_seq)

    def test_set_nucleotide_sequence_7(self):
        """Verify that expected sequence is extracted from a compound
        feature."""
        seq = Seq("AATTCGAGCT")
        self.feature.seqfeature = \
            SeqFeature(CompoundLocation([FeatureLocation(1, 5, strand=1),
                                         FeatureLocation(3, 7, strand=1)]),
                       type="CDS")
        self.feature.set_nucleotide_sequence(parent_genome_seq=seq)
        # feature #1 = ATTC
        # feature #2 = TCGA
        expected_seq = Seq("ATTCTCGA")
        self.assertEqual(self.feature.seq, expected_seq)

    def test_set_nucleotide_sequence_8(self):
        """Verify that empty Seq object set when no parameters selected."""
        self.feature.set_nucleotide_sequence()
        with self.subTest():
            self.assertEqual(self.feature.seq, "")
        with self.subTest():
            self.assertIsInstance(self.feature.seq, Seq)




    def test_check_translation_length_1(self):
        """Verify a present translation does not produce an error."""
        self.feature._translation_length = 1
        self.feature.check_translation_length()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_translation_length_2(self):
        """Verify a present translation does not produce an error."""
        self.feature._translation_length = 100
        self.feature.check_translation_length()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_translation_length_3(self):
        """Verify that no translation produces an error."""
        self.feature._translation_length = 0
        self.feature.check_translation_length()
        self.assertEqual(self.feature.evaluations[0].status, "error")








    def test_set_description_1(self):
        """Verify product description is assigned to primary description."""
        self.feature.product = "ABCD"
        self.feature.processed_product = "EFGH"
        self.feature.set_description("product")
        with self.subTest():
            self.assertEqual(self.feature.description, "ABCD")
        with self.subTest():
            self.assertEqual(self.feature.processed_description, "EFGH")

    def test_set_description_2(self):
        """Verify function description is assigned to primary description."""
        self.feature.function = "ABCD"
        self.feature.processed_function = "EFGH"
        self.feature.set_description("function")
        with self.subTest():
            self.assertEqual(self.feature.description, "ABCD")
        with self.subTest():
            self.assertEqual(self.feature.processed_description, "EFGH")

    def test_set_description_3(self):
        """Verify note description is assigned to primary description."""
        self.feature.note = "ABCD"
        self.feature.processed_note = "EFGH"
        self.feature.set_description("note")
        with self.subTest():
            self.assertEqual(self.feature.description, "ABCD")
        with self.subTest():
            self.assertEqual(self.feature.processed_description, "EFGH")

    def test_set_description_4(self):
        """Verify no description is assigned to primary description."""
        self.feature.note = "ABCD"
        self.feature.processed_note = "EFGH"
        self.feature.set_description("invalid")
        with self.subTest():
            self.assertEqual(self.feature.description, "")
        with self.subTest():
            self.assertEqual(self.feature.processed_description, "")




    def test_check_translation_table_present_1(self):
        """Verify no error is produced."""
        self.feature.translation_table = 11
        self.feature.check_translation_table_present()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_translation_table_present_2(self):
        """Verify an error is produced."""
        self.feature.translation_table = "11"
        self.feature.check_translation_table_present()
        self.assertEqual(self.feature.evaluations[0].status, "error")



    def test_check_translation_table_typo_1(self):
        """Verify no error is produced."""
        self.feature.translation_table = 11
        self.feature.check_translation_table_typo()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_translation_table_typo_2(self):
        """Verify an error is produced."""
        self.feature.translation_table = "11"
        self.feature.check_translation_table_typo()
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_translation_table_typo_3(self):
        """Verify no error is produced when a modified translation
        table is supplied."""
        self.feature.translation_table = "11"
        self.feature.check_translation_table_typo("11")
        self.assertEqual(self.feature.evaluations[0].status, "correct")




    def test_check_translation_1(self):
        """Verify no error is produced by a correct translation."""
        self.feature.translation = Seq("MF", IUPAC.protein)
        self.feature._translation_length = 2
        self.feature.seq = Seq("ATGTTTTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        self.feature.check_translation()
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_translation_2(self):
        """Verify an error is produced by a translation with an internal
        stop codon."""
        self.feature.translation = Seq("MF", IUPAC.protein)
        self.feature._translation_length = 2
        self.feature.seq = Seq("ATGTTTTGATGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        self.feature.check_translation()
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_translation_3(self):
        """Verify an error is produced by a translation shorter than
        expected."""
        self.feature.translation = Seq("MF", IUPAC.protein)
        self.feature._translation_length = 2
        self.feature.seq = Seq("ATGTTTATGTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        self.feature.check_translation()
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_translation_4(self):
        """Verify an error is produced by a translation longer than
        expected."""
        self.feature.translation = Seq("MF", IUPAC.protein)
        self.feature._translation_length = 2
        self.feature.seq = Seq("ATGTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        self.feature.check_translation()
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_translation_5(self):
        """Verify an error is produced by a translation different than
        expected (but same length)."""
        self.feature.translation = Seq("MF", IUPAC.protein)
        self.feature._translation_length = 2
        self.feature.seq = Seq("ATGATGTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        self.feature.check_translation()
        self.assertEqual(self.feature.evaluations[0].status, "error")







if __name__ == '__main__':
    unittest.main()
