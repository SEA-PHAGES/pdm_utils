""" Unit tests for the CDS class."""

from pdm_utils.classes import cds
from pdm_utils.constants import constants
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import unittest



class TestCdsClass(unittest.TestCase):


    def setUp(self):
        self.feature = cds.Cds()




    def test_set_locus_tag_1(self):
        """Verify that standard 3-part locus_tag is parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXIE_20")
        with self.subTest():
            self.assertEqual(self.feature.locus_tag, "SEA_TRIXIE_20")
        with self.subTest():
            self.assertEqual(self.feature._locus_tag_num, "20")

    def test_set_locus_tag_2(self):
        """Verify that standard 3-part locus_tag is parsed correctly
        when custom delimiter is provided."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA-TRIXIE-20", delimiter="-")
        self.assertEqual(self.feature._locus_tag_num, "20")

    def test_set_locus_tag_3(self):
        """Verify that standard 3-part locus_tag is parsed correctly
        when custom genome ID is provided."""
        self.feature.genome_id = "L5"
        self.feature.set_locus_tag(tag="SEA_TRIXIE_20", check_value="Trixie")
        self.assertEqual(self.feature._locus_tag_num, "20")

    def test_set_locus_tag_4(self):
        """Verify that non-standard 4-part locus_tag is
        parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXIE_DRAFT_20")
        self.assertEqual(self.feature._locus_tag_num, "20")

    def test_set_locus_tag_5(self):
        """Verify that non-standard 4-part locus_tag
        with no number is parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXIE_DRAFT_ABCD")
        self.assertEqual(self.feature._locus_tag_num, "")

    def test_set_locus_tag_6(self):
        """Verify that non-standard 2-part locus_tag
        with correct genome ID merged with number
        is partially parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXIE20")
        self.assertEqual(self.feature._locus_tag_num, "20")

    def test_set_locus_tag_7(self):
        """Verify that non-standard 2-part locus_tag
        with correct genome ID merged with no number
        is partially parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXIEABCD")
        self.assertEqual(self.feature._locus_tag_num, "")

    def test_set_locus_tag_8(self):
        """Verify that standard 3-part locus_tag
        with incorrect genome ID
        is partially parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXI_20")
        self.assertEqual(self.feature._locus_tag_num, "20")

    def test_set_locus_tag_9(self):
        """Verify that standard 3-part locus_tag
        with incorrect genome ID and no number
        is partially parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXI_AB20")
        self.assertEqual(self.feature._locus_tag_num, "")




    def test_set_name_1(self):
        """Verify that name is set from the parameter."""
        self.feature.gene = "1"
        self.feature._locus_tag_num = "2"
        self.feature.set_name(value="3")
        self.assertEqual(self.feature.name, "3")

    def test_set_name_2(self):
        """Verify that name is set from 'gene' attribute."""
        self.feature.gene = "1"
        self.feature._locus_tag_num = "2"
        self.feature.set_name()
        self.assertEqual(self.feature.name, "1")

    def test_set_name_3(self):
        """Verify that name is set from '_locus_tag_num' attribute."""
        self.feature.gene = ""
        self.feature._locus_tag_num = "2"
        self.feature.set_name()
        self.assertEqual(self.feature.name, "2")

    def test_set_name_4(self):
        """Verify that name is set as empty."""
        self.feature.gene = ""
        self.feature._locus_tag_num = ""
        self.feature.set_name()
        self.assertEqual(self.feature.name, "")




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
            self.assertEqual(self.feature.translation_length, 2)

    def test_set_translation_2(self):
        """Verify that sequence and length is set from valid input."""
        self.feature.set_translation("mf")
        with self.subTest():
            self.assertEqual(self.feature.translation, "MF")
        with self.subTest():
            self.assertEqual(self.feature.translation_length, 2)

    def test_set_translation_3(self):
        """Verify that sequence and length is set from invalid input."""
        self.feature.set_translation(1)
        with self.subTest():
            self.assertEqual(self.feature.translation, "")
        with self.subTest():
            self.assertEqual(self.feature.translation_length, 0)

    def test_set_translation_4(self):
        """Verify that sequence and length is set from translating
        valid nucleotide sequence."""
        self.feature.seq = Seq("ATGTTTTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        self.feature.set_translation(translate=True)
        with self.subTest():
            self.assertEqual(self.feature.translation, "MF")
        with self.subTest():
            self.assertEqual(self.feature.translation_length, 2)

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
            self.assertEqual(self.feature.translation_length, 0)

    def test_set_translation_6(self):
        """Verify that sequence and length is set when no options are
        selected."""
        self.feature.set_translation()
        with self.subTest():
            self.assertEqual(self.feature.translation, "")
        with self.subTest():
            self.assertEqual(self.feature.translation_length, 0)




    def test_check_amino_acids_1(self):
        """Verify no error is produced if all amino acids
        are in the protein alphabet."""
        self.feature.translation = "ADE"
        self.feature.check_amino_acids(
            check_set=constants.PROTEIN_ALPHABET, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_amino_acids_2(self):
        """Verify an error is produced if some amino acids
        are not in the protein alphabet."""
        self.feature.translation = "ABDE"
        self.feature.check_amino_acids(check_set=constants.PROTEIN_ALPHABET)
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_amino_acids_3(self):
        """Verify no error is produced if all amino acids
        are in a custom protein alphabet."""
        alphabet = set(["A","B","C", "D", "E", "F"])
        self.feature.translation = "ABDE"
        self.feature.check_amino_acids(check_set=alphabet)
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_amino_acids_4(self):
        """Verify an error is produced if some amino acids
        are not in a custom protein alphabet."""
        alphabet = set(["A","B","C", "D", "E", "F"])
        self.feature.translation = "ABDEG"
        self.feature.check_amino_acids(check_set=alphabet)
        self.assertEqual(self.feature.evaluations[0].status, "error")




    def test_get_start_end_1(self):
        """Forward strand feature, long format."""
        self.feature.left = 5
        self.feature.right = 10
        self.feature.strand = "forward"
        start, end = self.feature.get_start_end()
        with self.subTest():
            self.assertEqual(start, 5)
        with self.subTest():
            self.assertEqual(end, 10)

    def test_get_start_end_2(self):
        """Reverse strand feature, short format."""
        self.feature.left = 5
        self.feature.right = 10
        self.feature.strand = "r"
        start, end = self.feature.get_start_end()
        with self.subTest():
            self.assertEqual(start, 10)
        with self.subTest():
            self.assertEqual(end, 5)

    def test_get_start_end_3(self):
        """Non-standard strand feature."""
        self.feature.left = 5
        self.feature.right = 10
        self.feature.strand = "other"
        start, end = self.feature.get_start_end()
        with self.subTest():
            self.assertEqual(start, -1)
        with self.subTest():
            self.assertEqual(end, -1)





















    def test_set_location_id_1(self):
        """Forward strand feature, both values should be set."""
        self.feature.left = 5
        self.feature.right = 10
        self.feature.strand = "forward"
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
        location_id_2 = (-1, "abcd")
        location_id_3 = (-1, -1)
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
        self.feature.check_strand(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_strand_2(self):
        """Verify an error is produced when the strand is
        formatted incorrectly using default settings."""
        self.feature.strand = 1
        self.feature.check_strand()
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

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




    def test_check_coordinates_1(self):
        """Test correct boundaries."""
        self.feature.left = 5
        self.feature.right = 10
        self.feature.check_coordinates(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_coordinates_2(self):
        """Test incorrect left boundary."""
        self.feature.left = "a"
        self.feature.right = 10
        self.feature.check_coordinates()
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_coordinates_3(self):
        """Test incorrect right boundary."""
        self.feature.left = 5
        self.feature.right = "a"
        self.feature.check_coordinates()
        self.assertEqual(self.feature.evaluations[0].status, "error")




    def test_check_locus_tag_present_1(self):
        """Check if absent locus tag is expected to be absent."""
        self.feature.locus_tag = ""
        self.feature.check_locus_tag_present(False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_locus_tag_present_2(self):
        """Check if absent locus tag is expected to be present."""
        self.feature.locus_tag = ""
        self.feature.check_locus_tag_present(True)
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

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




    def test_check_gene_present_1(self):
        """Check if absent gene is expected to be absent."""
        self.feature.gene = ""
        self.feature.check_gene_present(False, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_gene_present_2(self):
        """Check if absent gene is expected to be present."""
        self.feature.gene = ""
        self.feature.check_gene_present(True)
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_gene_present_3(self):
        """Check if present gene is expected to be present."""
        self.feature.gene = "ABCD"
        self.feature.check_gene_present(True)
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_gene_present_4(self):
        """Check if present gene is expected to be absent."""
        self.feature.gene = "ABCD"
        self.feature.check_gene_present(False)
        self.assertEqual(self.feature.evaluations[0].status, "error")




    def test_check_gene_structure_1(self):
        """Verify no error is produced when gene is an integer."""
        self.feature.gene = "1"
        self.feature.check_gene_structure(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_gene_structure_2(self):
        """Verify an error is produced when gene is not integer."""
        self.feature.gene = "abcd"
        self.feature.check_gene_structure()
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)




    def test_check_compatible_gene_and_locus_tag_1(self):
        """Verify no error is produced when gene and locus_tag match."""
        self.feature.gene = "1"
        self.feature._locus_tag_num = "1"
        self.feature.check_compatible_gene_and_locus_tag(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_compatible_gene_and_locus_tag_2(self):
        """Verify an error is produced when gene and locus_tag do not match."""
        self.feature.gene = "1"
        self.feature._locus_tag_num = "10"
        self.feature.check_compatible_gene_and_locus_tag()
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)



    def test_check_description_field_1(self):
        """Verify no error is produced when a description is present in
        the processed_product as expected and the processed_function and
        the processed_note are empty."""
        self.feature.processed_product = "ABC"
        self.feature.check_description_field(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_description_field_2(self):
        """Verify no error is produced when a description is present in
        the processed_function as expected and the processed_product and
        the processed_note are empty."""
        self.feature.processed_function = "ABC"
        self.feature.check_description_field(attribute="function")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_description_field_3(self):
        """Verify no error is produced when a description is present in
        the processed_note as expected and the processed_product and
        the processed_function are empty."""
        self.feature.processed_note = "ABC"
        self.feature.check_description_field(attribute="note")
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_description_field_4(self):
        """Verify an error is produced when a description is not present in
        the processed_product but there is a description in
        the processed_function."""
        self.feature.processed_function = "ABC"
        self.feature.check_description_field()
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_description_field_5(self):
        """Verify an error is produced when a description is not present in
        the processed_function but there is a description in
        the processed_product."""
        self.feature.processed_product = "ABC"
        self.feature.check_description_field(attribute="function")
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_description_field_6(self):
        """Verify an error is produced when a description is not present in
        the processed_function but there is a description in
        the processed_product."""
        self.feature.processed_product = "ABC"
        self.feature.check_description_field(attribute="note")
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_description_field_7(self):
        """Verify an error is produced when a description is not present due
        to an invalid input but there is a description in
        the processed_product."""
        self.feature.processed_product = "ABC"
        self.feature.check_description_field(attribute="invalid")
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_description_field_8(self):
        """Verify no error is produced when a description is present in
        the processed_product as expected as well as the processed_function."""
        self.feature.processed_product = "ABC"
        self.feature.processed_function = "FGH"
        self.feature.check_description_field()
        self.assertEqual(self.feature.evaluations[0].status, "correct")




    def test_check_locus_tag_structure_1(self):
        """Verify no error is produced when the locus_tag has a
        correct structure."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "SEA_TRIXIE_123"
        self.feature.check_locus_tag_structure(
            prefix_set=constants.LOCUS_TAG_PREFIX_SET, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_locus_tag_structure_2(self):
        """Verify no error is produced when the locus_tag has a
        correct structure and the 'only_typo' parameter is chosen."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "SEATrixie123"
        self.feature.check_locus_tag_structure(only_typo=True,
            prefix_set=constants.LOCUS_TAG_PREFIX_SET)
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_locus_tag_structure_3(self):
        """Verify an error is produced when the locus_tag has an
        incorrect structure and the 'only_typo' parameter is chosen."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "SEATrixi123"
        self.feature.check_locus_tag_structure(only_typo=True,
            prefix_set=constants.LOCUS_TAG_PREFIX_SET)
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_locus_tag_structure_4(self):
        """Verify an error is produced when the locus_tag is
        not capitalized."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "sea_trixie_123"
        self.feature.check_locus_tag_structure(
            prefix_set=constants.LOCUS_TAG_PREFIX_SET)
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_locus_tag_structure_5(self):
        """Verify no error is produced when the locus_tag is
        not capitalized but the 'case' parameter is set to False."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "sea_trixie_123"
        self.feature.check_locus_tag_structure(case=False,
            prefix_set=constants.LOCUS_TAG_PREFIX_SET)
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_locus_tag_structure_6(self):
        """Verify an error is produced when the locus_tag has an
        incorrect number of parts."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "ABCTRIXIE_123"
        self.feature.check_locus_tag_structure(
            prefix_set=constants.LOCUS_TAG_PREFIX_SET)
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_locus_tag_structure_7(self):
        """Verify an error is produced when the locus_tag has an
        incorrect prefix."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "ABC_TRIXIE_123"
        self.feature.check_locus_tag_structure(
            prefix_set=constants.LOCUS_TAG_PREFIX_SET)
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_locus_tag_structure_8(self):
        """Verify no error is produced when the locus_tag has an
        incorrect prefix but 'prefix_set' parameters is set to None."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "ABC_TRIXIE_123"
        self.feature.check_locus_tag_structure(prefix_set=None)
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_locus_tag_structure_9(self):
        """Verify no error is produced when the locus_tag has an
        incorrect prefix but 'prefix_set' parameters is set to new set."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "ABC_TRIXIE_123"
        self.feature.check_locus_tag_structure(prefix_set=set(["ABC"]))
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_locus_tag_structure_10(self):
        """Verify an error is produced when the locus_tag has a
        common prefix but 'prefix_set' parameters is set to new set."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "SEA_TRIXIE_123"
        self.feature.check_locus_tag_structure(prefix_set=set(["ABC"]))
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_locus_tag_structure_11(self):
        """Verify an error is produced when the locus_tag has an
        incorrect genome due to misspelling."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "SEA_TRIXIEX_123"
        self.feature.check_locus_tag_structure(
            prefix_set=constants.LOCUS_TAG_PREFIX_SET)
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_locus_tag_structure_12(self):
        """Verify an error is produced when the locus_tag has an
        incorrect number."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "SEA_TRIXIE_123x"
        self.feature.check_locus_tag_structure(
            prefix_set=constants.LOCUS_TAG_PREFIX_SET)
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_locus_tag_structure_13(self):
        """Verify an error is produced when the locus_tag has an
        incorrect prefix, genome, and number."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "ABC_TRIXIEX_123X"
        self.feature.check_locus_tag_structure(
            prefix_set=constants.LOCUS_TAG_PREFIX_SET)
        results = self.feature.evaluations[0].result.split(".")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertTrue(results[1].strip().startswith("The prefix"))
        with self.subTest():
            self.assertTrue(results[2].strip().startswith("The genome"))
        with self.subTest():
            self.assertTrue(results[3].strip().startswith("The feature"))

    def test_check_locus_tag_structure_14(self):
        """Verify an error is produced when the locus_tag has an
        incorrect genome when the 'check_value' parameter is used."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "SEA_TRIXIE_123"
        self.feature.check_locus_tag_structure(check_value="L5",
            prefix_set=constants.LOCUS_TAG_PREFIX_SET)
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_locus_tag_structure_15(self):
        """Verify no error is produced when the locus_tag has a
        correct genome when the 'check_value' parameter is used."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "SEA_L5_123"
        self.feature.check_locus_tag_structure(check_value="L5",
            prefix_set=constants.LOCUS_TAG_PREFIX_SET)
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_locus_tag_structure_16(self):
        """Verify no error is produced when the locus_tag has a correct
        genome when 'only_typo' and 'check_value' parameters are used."""
        self.feature.genome_id = "Trixie"
        self.feature.locus_tag = "SEAL5123"
        self.feature.check_locus_tag_structure(only_typo=True,
            check_value="L5", prefix_set=constants.LOCUS_TAG_PREFIX_SET)
        self.assertEqual(self.feature.evaluations[0].status, "correct")









    def test_set_nucleotide_length_1(self):
        """Verify the nucleotide length is correct when computed directly
        from a sequence."""
        self.feature.seq = Seq("ATCG", IUPAC.unambiguous_dna)
        self.feature.set_nucleotide_length(seq=True)
        self.assertEqual(self.feature.length, 4)

    def test_set_nucleotide_length_2(self):
        """Verify the nucleotide length is correct for a 0-based
        half-open interval."""
        self.feature.left = 0
        self.feature.right = 11
        self.feature.coordinate_format = "0_half_open"
        self.feature.set_nucleotide_length()
        self.assertEqual(self.feature.length, 11)

    def test_set_nucleotide_length_3(self):
        """Verify the nucleotide length is correct for a 1-based
        closed interval."""
        self.feature.left = 0
        self.feature.right = 11
        self.feature.coordinate_format = "1_closed"
        self.feature.set_nucleotide_length()
        self.assertEqual(self.feature.length, 12)

    def test_set_nucleotide_length_4(self):
        """Verify the nucleotide length is not set for invalid
        coordinate format."""
        self.feature.left = 0
        self.feature.right = 11
        self.feature.coordinate_format = "invalid"
        self.feature.set_nucleotide_length()
        self.assertEqual(self.feature.length, -1)







    def test_reformat_left_and_right_1(self):
        """Verify the coordinates are converted to 1-based closed interval."""
        self.feature.left = 5
        self.feature.right = 11
        self.feature.coordinate_format = "0_half_open"
        new_format = "1_closed"
        self.feature.reformat_left_and_right(new_format)
        with self.subTest():
            self.assertEqual(self.feature.left, 6)
        with self.subTest():
            self.assertEqual(self.feature.right, 11)
        with self.subTest():
            self.assertEqual(self.feature.coordinate_format, new_format)

    def test_reformat_left_and_right_2(self):
        """Verify the coordinates are converted to 0-based half open interval."""
        self.feature.left = 5
        self.feature.right = 11
        self.feature.coordinate_format = "1_closed"
        new_format = "0_half_open"
        self.feature.reformat_left_and_right(new_format)
        with self.subTest():
            self.assertEqual(self.feature.left, 4)
        with self.subTest():
            self.assertEqual(self.feature.right, 11)
        with self.subTest():
            self.assertEqual(self.feature.coordinate_format, new_format)

    def test_reformat_left_and_right_3(self):
        """Verify the coordinates are not converted."""
        self.feature.left = 5
        self.feature.right = 11
        self.feature.coordinate_format = "1_closed"
        new_format = "invalid"
        self.feature.reformat_left_and_right(new_format)
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






    def test_set_seqfeature_1(self):
        """Verify seqfeature is set correctly with no changes."""
        self.feature.left = 2
        self.feature.right = 5
        self.feature.strand = 1
        self.feature.coordinate_format = "0_half_open"
        self.feature.set_seqfeature()
        with self.subTest():
            self.assertEqual(self.feature.seqfeature.strand, 1)
        with self.subTest():
            self.assertEqual(
                self.feature.seqfeature.location.start.position, 2)
        with self.subTest():
            self.assertEqual(
                self.feature.seqfeature.location.end.position, 5)

    def test_set_seqfeature_2(self):
        """Verify seqfeature is set correctly with left coordinate and
        strand reformatted."""
        self.feature.left = 2
        self.feature.right = 5
        self.feature.strand = "F"
        self.feature.coordinate_format = "1_closed"
        self.feature.set_seqfeature()
        with self.subTest():
            self.assertEqual(self.feature.seqfeature.strand, 1)
        with self.subTest():
            self.assertEqual(
                self.feature.seqfeature.location.start.position, 1)
        with self.subTest():
            self.assertEqual(
                self.feature.seqfeature.location.end.position, 5)




    def test_check_translation_length_1(self):
        """Verify a present translation does not produce an error."""
        self.feature.translation_length = 1
        self.feature.check_translation_length(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_translation_length_2(self):
        """Verify a present translation does not produce an error."""
        self.feature.translation_length = 100
        self.feature.check_translation_length()
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_translation_length_3(self):
        """Verify that no translation produces an error."""
        self.feature.translation_length = 0
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




    def test_check_translation_table_1(self):
        """Verify no error is produced."""
        self.feature.translation_table = 11
        self.feature.check_translation_table(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_translation_table_2(self):
        """Verify an error is produced."""
        self.feature.translation_table = "11"
        self.feature.check_translation_table()
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_translation_table_3(self):
        """Verify no error is produced when a modified translation
        table is supplied."""
        self.feature.translation_table = "11"
        self.feature.check_translation_table("11")
        self.assertEqual(self.feature.evaluations[0].status, "correct")




    def test_check_translation_1(self):
        """Verify no error is produced by a correct translation."""
        self.feature.translation = Seq("MF", IUPAC.protein)
        self.feature.translation_length = 2
        self.feature.seq = Seq("ATGTTTTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        self.feature.check_translation(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_translation_2(self):
        """Verify an error is produced by a translation with an internal
        stop codon."""
        self.feature.translation = Seq("MF", IUPAC.protein)
        self.feature.translation_length = 2
        self.feature.seq = Seq("ATGTTTTGATGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        self.feature.check_translation()
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_translation_3(self):
        """Verify an error is produced by a translation shorter than
        expected."""
        self.feature.translation = Seq("MF", IUPAC.protein)
        self.feature.translation_length = 2
        self.feature.seq = Seq("ATGTTTATGTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        self.feature.check_translation()
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_translation_4(self):
        """Verify an error is produced by a translation longer than
        expected."""
        self.feature.translation = Seq("MF", IUPAC.protein)
        self.feature.translation_length = 2
        self.feature.seq = Seq("ATGTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        self.feature.check_translation()
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_translation_5(self):
        """Verify an error is produced by a translation different than
        expected (but same length)."""
        self.feature.translation = Seq("MF", IUPAC.protein)
        self.feature.translation_length = 2
        self.feature.seq = Seq("ATGATGTGA", IUPAC.unambiguous_dna)
        self.feature.translation_table = 11
        self.feature.check_translation()
        self.assertEqual(self.feature.evaluations[0].status, "error")




    def test_check_generic_data_1(self):
        """Verify no error is produced if the product contains valid data."""
        self.feature.product = "terminase"
        self.feature.processed_product = "terminase"
        self.feature.check_generic_data("product", eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_generic_data_2(self):
        """Verify no error is produced if the function contains valid data."""
        self.feature.function = "terminase"
        self.feature.processed_function = "terminase"
        self.feature.check_generic_data("product")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_generic_data_3(self):
        """Verify no error is produced if the note contains valid data."""
        self.feature.note = "terminase"
        self.feature.processed_note = "terminase"
        self.feature.check_generic_data("product")
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_generic_data_4(self):
        """Verify an error is produced if the product contains invalid data."""
        self.feature.product = "gp104"
        self.feature.processed_product = ""
        self.feature.check_generic_data("product")
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_generic_data_5(self):
        """Verify no error is produced if no attribute is selected."""
        self.feature.product = "gp104"
        self.feature.processed_product = ""
        self.feature.check_generic_data()
        self.assertEqual(self.feature.evaluations[0].status, "correct")










if __name__ == '__main__':
    unittest.main()
