""" Unit tests for the CDS class."""


from pathlib import Path
import unittest
import sys

from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq

from pdm_utils.classes import cds
from pdm_utils.constants import constants

# Import helper functions to build mock database and mock flat files
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_data_utils



class TestCdsClass(unittest.TestCase):

    def setUp(self):
        self.feature = cds.Cds()
        self.seq1 = Seq("AATTCGAGCT")
        self.seqfeature1 = test_data_utils.create_1_part_seqfeature(
                                1, 5, 1, "CDS")
        self.seqfeature2 = test_data_utils.create_1_part_seqfeature(
                                1, 5, -1, "CDS")
        self.seqfeature3 = test_data_utils.create_2_part_seqfeature(
                                1, 5, 1, 3, 7, 1, "CDS")




    def test_set_locus_tag_1(self):
        """Verify that standard 3-part locus_tag is parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXIE_20", delimiter=None)
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
        self.feature.set_locus_tag(tag="SEA_TRIXIE_20", delimiter=None,
                                   check_value="Trixie")
        self.assertEqual(self.feature._locus_tag_num, "20")

    def test_set_locus_tag_4(self):
        """Verify that non-standard 4-part locus_tag is
        parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXIE_DRAFT_20", delimiter=None)
        self.assertEqual(self.feature._locus_tag_num, "20")

    def test_set_locus_tag_5(self):
        """Verify that non-standard 4-part locus_tag
        with no number is parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXIE_DRAFT_ABCD", delimiter=None)
        self.assertEqual(self.feature._locus_tag_num, "ABCD")

    def test_set_locus_tag_6(self):
        """Verify that non-standard 2-part locus_tag
        with correct genome ID merged with number
        is partially parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXIE20", delimiter=None)
        self.assertEqual(self.feature._locus_tag_num, "20")

    def test_set_locus_tag_7(self):
        """Verify that non-standard 2-part locus_tag
        with correct genome ID merged with no number
        is partially parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXIEABCD", delimiter=None)
        self.assertEqual(self.feature._locus_tag_num, "ABCD")

    def test_set_locus_tag_8(self):
        """Verify that standard 3-part locus_tag
        with incorrect genome ID
        is partially parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXI_20", delimiter=None)
        self.assertEqual(self.feature._locus_tag_num, "20")

    def test_set_locus_tag_9(self):
        """Verify that standard 3-part locus_tag
        with incorrect genome ID and no number
        is partially parsed correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="SEA_TRIXI_AB20", delimiter=None)
        self.assertEqual(self.feature._locus_tag_num, "20")

    def test_set_locus_tag_10(self):
        """Verify that empty locus_tag data from database
        is stored correctly."""
        self.feature.set_locus_tag(tag=None, delimiter=None)
        with self.subTest():
            self.assertEqual(self.feature.locus_tag, "")
        with self.subTest():
            self.assertEqual(self.feature._locus_tag_num, "")

    def test_set_locus_tag_11(self):
        """Verify that locus_tag with float is stored correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="TRIXIE_20.1", delimiter=None)
        with self.subTest():
            self.assertEqual(self.feature.locus_tag, "TRIXIE_20.1")
        with self.subTest():
            self.assertEqual(self.feature._locus_tag_num, "20.1")

    def test_set_locus_tag_12(self):
        """Verify that locus_tag with generic 'GP' prefix and
        float is stored correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="TRIXIE_GP020.1", delimiter=None)
        with self.subTest():
            self.assertEqual(self.feature.locus_tag, "TRIXIE_GP020.1")
        with self.subTest():
            self.assertEqual(self.feature._locus_tag_num, "020.1")

    def test_set_locus_tag_13(self):
        """Verify that locus_tag with generic 'ORF' prefix and
        float is stored correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="TRIXIE_ORF020.1", delimiter=None)
        with self.subTest():
            self.assertEqual(self.feature.locus_tag, "TRIXIE_ORF020.1")
        with self.subTest():
            self.assertEqual(self.feature._locus_tag_num, "020.1")

    def test_set_locus_tag_14(self):
        """Verify that locus_tag with number that includes alpha suffix
        is stored correctly."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="TRIXIE_ORF20a", delimiter=None)
        with self.subTest():
            self.assertEqual(self.feature.locus_tag, "TRIXIE_ORF20a")
        with self.subTest():
            self.assertEqual(self.feature._locus_tag_num, "20a")

    def test_set_locus_tag_15(self):
        """Verify that locus_tag with multiple potential delimiters
        is stored correctly with explicit '_' delimiter."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="TRIXIE_GP020.1", delimiter="_")
        with self.subTest():
            self.assertEqual(self.feature.locus_tag, "TRIXIE_GP020.1")
        with self.subTest():
            self.assertEqual(self.feature._locus_tag_num, "020.1")

    def test_set_locus_tag_16(self):
        """Verify that locus_tag with multiple potential delimiters
        is stored correctly with explicit '.' delimiter."""
        self.feature.genome_id = "Trixie"
        self.feature.set_locus_tag(tag="TRIXIE_GP020.1", delimiter=".")
        with self.subTest():
            self.assertEqual(self.feature.locus_tag, "TRIXIE_GP020.1")
        with self.subTest():
            self.assertEqual(self.feature._locus_tag_num, "1")




    def test_set_name_1(self):
        """Verify that name is set from the supplied parameter value."""
        self.feature._locus_tag_num = "2"
        self.feature._gene_num = "3"
        self.feature._product_num = "4"
        self.feature._note_num = "5"
        self.feature._function_num = "6"
        self.feature.gene = "7"
        self.feature.set_name(value="1")
        self.assertEqual(self.feature.name, "1")

    def test_set_name_2(self):
        """Verify that name is set from '_locus_tag_num' attribute."""
        self.feature._locus_tag_num = "2"
        self.feature._gene_num = "3"
        self.feature._product_num = "4"
        self.feature._note_num = "5"
        self.feature._function_num = "6"
        self.feature.gene = "7"
        self.feature.set_name()
        self.assertEqual(self.feature.name, "2")

    def test_set_name_3(self):
        """Verify that name is set from '_gene_num' attribute."""
        self.feature._locus_tag_num = ""
        self.feature._gene_num = "3"
        self.feature._product_num = "4"
        self.feature._note_num = "5"
        self.feature._function_num = "6"
        self.feature.gene = "7"
        self.feature.set_name()
        self.assertEqual(self.feature.name, "3")

    def test_set_name_4(self):
        """Verify that name is set from '_product_num' attribute."""
        self.feature._locus_tag_num = ""
        self.feature._gene_num = ""
        self.feature._product_num = "4"
        self.feature._note_num = "5"
        self.feature._function_num = "6"
        self.feature.gene = "7"
        self.feature.set_name()
        self.assertEqual(self.feature.name, "4")

    def test_set_name_5(self):
        """Verify that name is set from '_note_num' attribute."""
        self.feature._locus_tag_num = ""
        self.feature._gene_num = ""
        self.feature._product_num = ""
        self.feature._note_num = "5"
        self.feature._function_num = "6"
        self.feature.gene = "7"
        self.feature.set_name()
        self.assertEqual(self.feature.name, "5")

    def test_set_name_6(self):
        """Verify that name is set from '_function_num' attribute."""
        self.feature._locus_tag_num = ""
        self.feature._gene_num = ""
        self.feature._product_num = ""
        self.feature._note_num = ""
        self.feature._function_num = "6"
        self.feature.gene = "7"
        self.feature.set_name()
        self.assertEqual(self.feature.name, "6")

    def test_set_name_7(self):
        """Verify that name is set from 'gene' attribute."""
        self.feature._locus_tag_num = ""
        self.feature._gene_num = ""
        self.feature._product_num = ""
        self.feature._note_num = ""
        self.feature._function_num = ""
        self.feature.gene = "7"
        self.feature.set_name()
        self.assertEqual(self.feature.name, "7")

    def test_set_name_8(self):
        """Verify that name is set as empty."""
        self.feature._locus_tag_num = ""
        self.feature._gene_num = ""
        self.feature._product_num = ""
        self.feature._note_num = ""
        self.feature._function_num = ""
        self.feature.gene = ""
        self.feature.set_name()
        self.assertEqual(self.feature.name, "")

    def test_set_name_9(self):
        """Verify that name is set from first float attribute."""
        self.feature._locus_tag_num = "2b"
        self.feature._gene_num = "3b"
        self.feature._product_num = "4"
        self.feature._note_num = "5"
        self.feature._function_num = "6"
        self.feature.gene = "7"
        self.feature.set_name()
        self.assertEqual(self.feature.name, "4")

    def test_set_name_10(self):
        """Verify that name is set from first non-float attribute."""
        self.feature._locus_tag_num = "2b"
        self.feature._gene_num = "3b"
        self.feature._product_num = "4b"
        self.feature._note_num = "5b"
        self.feature._function_num = "6b"
        self.feature.gene = "7b"
        self.feature.set_name()
        self.assertEqual(self.feature.name, "7b")

    def test_set_name_11(self):
        """Verify that name is set from first non-float attribute."""
        self.feature._locus_tag_num = ""
        self.feature._gene_num = ""
        self.feature._product_num = ""
        self.feature._note_num = "5b"
        self.feature._function_num = "6b"
        self.feature.gene = ""
        self.feature.set_name()
        self.assertEqual(self.feature.name, "5b")




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




    def test_set_orientation_1(self):
        """Verify 'f' is converted correctly."""
        self.feature.set_orientation("f", "fr_long")
        self.assertEqual(self.feature.orientation, "forward")

    def test_set_orientation_2(self):
        """Verify 'r' is converted correctly."""
        self.feature.set_orientation("reverse", "fr_short")
        self.assertEqual(self.feature.orientation, "r")




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




    def test_get_begin_end_1(self):
        """Forward orientation feature, long format."""
        self.feature.start = 5
        self.feature.stop = 10
        self.feature.orientation = "forward"
        start, end = self.feature.get_begin_end()
        with self.subTest():
            self.assertEqual(start, 5)
        with self.subTest():
            self.assertEqual(end, 10)

    def test_get_begin_end_2(self):
        """Reverse orientation feature, short format."""
        self.feature.start = 5
        self.feature.stop = 10
        self.feature.orientation = "r"
        start, end = self.feature.get_begin_end()
        with self.subTest():
            self.assertEqual(start, 10)
        with self.subTest():
            self.assertEqual(end, 5)

    def test_get_begin_end_3(self):
        """Non-standard orientation feature."""
        self.feature.start = 5
        self.feature.stop = 10
        self.feature.orientation = "other"
        start, end = self.feature.get_begin_end()
        with self.subTest():
            self.assertEqual(start, -1)
        with self.subTest():
            self.assertEqual(end, -1)




    def test_set_location_id_1(self):
        """Forward orientation feature, both values should be set."""
        self.feature.start = 5
        self.feature.stop = 10
        self.feature.orientation = "forward"
        location_id_1 = (5, 10, "forward")
        location_id_2 = (10, "forward")
        location_id_3 = (5, 10)
        self.feature.set_location_id()
        with self.subTest():
            self.assertEqual(self.feature._start_stop_orient_id, location_id_1)
        with self.subTest():
            self.assertEqual(self.feature._end_orient_id, location_id_2)
        with self.subTest():
            self.assertEqual(self.feature._start_end_id, location_id_3)

    def test_set_location_id_2(self):
        """Reverse orientation feature, both values should be set."""
        self.feature.start = 5
        self.feature.stop = 10
        self.feature.orientation = "reverse"
        location_id_1 = (5, 10, "reverse")
        location_id_2 = (5, "reverse")
        location_id_3 = (10, 5)
        self.feature.set_location_id()
        with self.subTest():
            self.assertEqual(self.feature._start_stop_orient_id, location_id_1)
        with self.subTest():
            self.assertEqual(self.feature._end_orient_id, location_id_2)
        with self.subTest():
            self.assertEqual(self.feature._start_end_id, location_id_3)

    def test_set_location_id_3(self):
        """Test forward orientation numeric format."""
        self.feature.start = 5
        self.feature.stop = 10
        self.feature.orientation = 1
        location_id_1 = (5, 10, 1)
        location_id_2 = (10, 1)
        location_id_3 = (5, 10)
        self.feature.set_location_id()
        with self.subTest():
            self.assertEqual(self.feature._start_stop_orient_id, location_id_1)
        with self.subTest():
            self.assertEqual(self.feature._end_orient_id, location_id_2)
        with self.subTest():
            self.assertEqual(self.feature._start_end_id, location_id_3)

    def test_set_location_id_4(self):
        """Test non-standard orientation format."""
        self.feature.start = 5
        self.feature.stop = 10
        self.feature.orientation = "abcd"
        location_id_1 = (5, 10, "abcd")
        location_id_2 = (-1, "abcd")
        location_id_3 = (-1, -1)
        self.feature.set_location_id()
        with self.subTest():
            self.assertEqual(self.feature._start_stop_orient_id, location_id_1)
        with self.subTest():
            self.assertEqual(self.feature._end_orient_id, location_id_2)
        with self.subTest():
            self.assertEqual(self.feature._start_end_id, location_id_3)




    def test_check_attribute_1(self):
        """Verify no error is produced when the id
        is in the check_set and is expected to be in the set."""
        check_set = set(["Trixie_CDS_1", "L5_CDS_1"])
        self.feature.id = "Trixie_CDS_1"
        self.feature.check_attribute("id", check_set, True, "eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_attribute_2(self):
        """Verify an error is produced when the id
        is not in the check_set and is expected to be in the set."""
        check_set = set(["Trixie_CDS_1", "L5_CDS_1"])
        self.feature.id = "Trixie_CDS_2"
        self.feature.check_attribute("id", check_set, True)
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_attribute_3(self):
        """Verify no error is produced when the id
        is not in the check_set and is not expected to be in the set."""
        check_set = set(["Trixie_CDS_1", "L5_CDS_1"])
        self.feature.id = "Trixie_CDS_2"
        self.feature.check_attribute("id", check_set, False)
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_attribute_4(self):
        """Verify an error is produced when the id
        is in the check_set and is not expected to be in the set."""
        check_set = set(["Trixie_CDS_1", "L5_CDS_1"])
        self.feature.id = "Trixie_CDS_1"
        self.feature.check_attribute("id", check_set, False)
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_attribute_5(self):
        """Verify an error is produced when the translation
        is in the check_set and is expected to be in the set."""
        check_set = {Seq("MF", IUPAC.protein), Seq("MA", IUPAC.protein)}
        self.feature.translation = Seq("MF", IUPAC.protein)
        self.feature.check_attribute("translation", check_set, True)
        self.assertEqual(self.feature.evaluations[0].status, "correct")




    def test_check_magnitude_1(self):
        """Verify no error is produced when
        'start' is greater than 0, as expected."""
        self.feature.start = 1000
        self.feature.check_magnitude("start", ">", 0, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_magnitude_2(self):
        """Verify no error is produced when
        'start' is equal to 0, as expected."""
        self.feature.start = 0
        self.feature.check_magnitude("start", "=", 0)
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_magnitude_3(self):
        """Verify no error is produced when
        'start' is less than 0, as expected."""
        self.feature.start = -100
        self.feature.check_magnitude("start", "<", 0)
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_magnitude_4(self):
        """Verify an error is produced when
        'start' is greater than 0, unexpectedly."""
        self.feature.start = 100
        self.feature.check_magnitude("start", "=", 0)
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_magnitude_5(self):
        """Verify an error is produced when
        'start' is less than 0, unexpectedly."""
        self.feature.start = -100
        self.feature.check_magnitude("start", ">", 0)
        self.assertEqual(self.feature.evaluations[0].status, "error")




    def test_check_orientation_1(self):
        """Verify no error is produced when the orientation is
        formatted correctly using default settings."""
        self.feature.orientation = "F"
        self.feature.check_orientation(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_orientation_2(self):
        """Verify an error is produced when the orientation is
        formatted incorrectly using default settings."""
        self.feature.orientation = 1
        self.feature.check_orientation()
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_orientation_3(self):
        """Verify no error is produced when the orientation is
        formatted correctly using custom settings."""
        self.feature.orientation = 1
        self.feature.check_orientation(format="numeric", case=False)
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_orientation_4(self):
        """Verify an error is produced when the orientation is
        formatted incorrectly using custom settings."""
        self.feature.orientation = "F"
        self.feature.check_orientation(format="numeric", case=False)
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
        the product as expected and the function and
        the note are empty."""
        self.feature.product = "ABC"
        self.feature.check_description_field(eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_description_field_2(self):
        """Verify no error is produced when a description is present in
        the function as expected and the product and
        the note are empty."""
        self.feature.function = "ABC"
        self.feature.check_description_field(attribute="function")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_description_field_3(self):
        """Verify no error is produced when a description is present in
        the note as expected and the product and
        the function are empty."""
        self.feature.note = "ABC"
        self.feature.check_description_field(attribute="note")
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_description_field_4(self):
        """Verify an error is produced when a description is not present in
        the product but there is a description in
        the function."""
        self.feature.function = "ABC"
        self.feature.check_description_field()
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_description_field_5(self):
        """Verify an error is produced when a description is not present in
        the function but there is a description in
        the product."""
        self.feature.product = "ABC"
        self.feature.check_description_field(attribute="function")
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_description_field_6(self):
        """Verify an error is produced when a description is not present in
        the function but there is a description in
        the product."""
        self.feature.product = "ABC"
        self.feature.check_description_field(attribute="note")
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_description_field_7(self):
        """Verify an error is produced when a description is not present due
        to an invalid input but there is a description in
        the product."""
        self.feature.product = "ABC"
        self.feature.check_description_field(attribute="invalid")
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_description_field_8(self):
        """Verify no error is produced when a description is present in
        the product as expected as well as the function."""
        self.feature.product = "ABC"
        self.feature.function = "FGH"
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
            self.assertTrue(results[2].strip().startswith("The prefix"))
        with self.subTest():
            self.assertTrue(results[3].strip().startswith("The genome"))
        with self.subTest():
            self.assertTrue(results[4].strip().startswith("The feature"))

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
        """Verify the nucleotide length is correct when computed directly
        from a translation."""
        self.feature.translation = Seq("MF", IUPAC.protein)
        self.feature.set_nucleotide_length(translation=True)
        self.assertEqual(self.feature.length, 9)

    def test_set_nucleotide_length_3(self):
        """Verify the nucleotide length is correct for a 0-based
        half-open interval."""
        self.feature.start = 0
        self.feature.stop = 11
        self.feature.coordinate_format = "0_half_open"
        self.feature.set_nucleotide_length()
        self.assertEqual(self.feature.length, 11)

    def test_set_nucleotide_length_4(self):
        """Verify the nucleotide length is correct for a 1-based
        closed interval."""
        self.feature.start = 0
        self.feature.stop = 11
        self.feature.coordinate_format = "1_closed"
        self.feature.set_nucleotide_length()
        self.assertEqual(self.feature.length, 12)

    def test_set_nucleotide_length_5(self):
        """Verify the nucleotide length is not set for invalid
        coordinate format."""
        self.feature.start = 0
        self.feature.stop = 11
        self.feature.coordinate_format = "invalid"
        self.feature.set_nucleotide_length()
        self.assertEqual(self.feature.length, -1)




    def test_reformat_start_and_stop_1(self):
        """Verify the coordinates are converted to 1-based closed interval."""
        self.feature.start = 5
        self.feature.stop = 11
        self.feature.coordinate_format = "0_half_open"
        new_format = "1_closed"
        self.feature.reformat_start_and_stop(new_format)
        with self.subTest():
            self.assertEqual(self.feature.start, 6)
        with self.subTest():
            self.assertEqual(self.feature.stop, 11)
        with self.subTest():
            self.assertEqual(self.feature.coordinate_format, new_format)

    def test_reformat_start_and_stop_2(self):
        """Verify the coordinates are converted to 0-based half open interval."""
        self.feature.start = 5
        self.feature.stop = 11
        self.feature.coordinate_format = "1_closed"
        new_format = "0_half_open"
        self.feature.reformat_start_and_stop(new_format)
        with self.subTest():
            self.assertEqual(self.feature.start, 4)
        with self.subTest():
            self.assertEqual(self.feature.stop, 11)
        with self.subTest():
            self.assertEqual(self.feature.coordinate_format, new_format)

    def test_reformat_start_and_stop_3(self):
        """Verify the coordinates are not converted."""
        self.feature.start = 5
        self.feature.stop = 11
        self.feature.coordinate_format = "1_closed"
        new_format = "invalid"
        with self.assertRaises(ValueError):
            self.feature.reformat_start_and_stop(new_format)
        with self.subTest():
            self.assertEqual(self.feature.start, 5)
        with self.subTest():
            self.assertEqual(self.feature.stop, 11)
        with self.subTest():
            self.assertEqual(self.feature.coordinate_format, "1_closed")




    def test_set_nucleotide_sequence_1(self):
        """Verify that supplied Seq object is set correctly."""
        seq = Seq("aattcgagct")
        self.feature.set_nucleotide_sequence(value=seq)
        with self.subTest():
            self.assertEqual(self.feature.seq, str(self.seq1))
        with self.subTest():
            self.assertIsInstance(self.feature.seq, Seq)

    def test_set_nucleotide_sequence_2(self):
        """Verify that supplied sequence is set correclty
        (converted to Seq object)."""
        seq = "aattcgagct"
        self.feature.set_nucleotide_sequence(value=seq)
        with self.subTest():
            self.assertEqual(self.feature.seq, str(self.seq1))
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
        """Verify that expected sequence is extracted from top orientation."""
        self.feature.seqfeature = self.seqfeature1
        self.feature.set_nucleotide_sequence(parent_genome_seq=self.seq1)
        expected_seq = Seq("ATTC")
        self.assertEqual(self.feature.seq, expected_seq)

    def test_set_nucleotide_sequence_5(self):
        """Verify that expected sequence is extracted from bottom orientation."""
        self.feature.seqfeature = self.seqfeature2
        self.feature.set_nucleotide_sequence(parent_genome_seq=self.seq1)
        expected_seq = Seq("GAAT")
        self.assertEqual(self.feature.seq, expected_seq)

    def test_set_nucleotide_sequence_6(self):
        """Verify that no sequence is extracted if the 'seqfeature'
        attribute is not a Biopython SeqFeature object."""
        self.feature.seqfeature = ""
        self.feature.set_nucleotide_sequence(parent_genome_seq=self.seq1)
        expected_seq = Seq("")
        self.assertEqual(self.feature.seq, expected_seq)

    def test_set_nucleotide_sequence_7(self):
        """Verify that expected sequence is extracted from a compound
        feature."""
        self.feature.seqfeature = self.seqfeature3
        self.feature.set_nucleotide_sequence(parent_genome_seq=self.seq1)
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








    def test_set_num_1(self):
        """Verify non-generic product with no number description
        is parsed correctly."""
        attr = "_product_num"
        value = " terminase "
        self.feature.set_num(attr, value, None, None)
        self.assertEqual(self.feature._product_num, "")

    def test_set_num_2(self):
        """Verify non-generic product with number description
        is parsed correctly."""
        attr = "_product_num"
        value = " terminase10 "
        self.feature.set_num(attr, value, None, None)
        self.assertEqual(self.feature._product_num, "")

    def test_set_num_3(self):
        """Verify non-generic product description with generic 'GP' number
        is parsed correctly."""
        attr = "_product_num"
        value = " GP10; terminase "
        self.feature.set_num(attr, value, None, None)
        self.assertEqual(self.feature._product_num, "10")

    def test_set_num_4(self):
        """Verify non-generic product description with generic 'GP' number
        at end of string is parsed correctly."""
        attr = "_product_num"
        value = " terminase; GP10 "
        self.feature.set_num(attr, value, None, None)
        self.assertEqual(self.feature._product_num, "10")

    def test_set_num_5(self):
        """Verify non-generic product description with generic 'GP' number
        but with new delimiter is parsed correctly."""
        attr = "_product_num"
        value = " GP10; terminase "
        self.feature.set_num(attr, value, "-", None)
        self.assertEqual(self.feature._product_num, "")

    def test_set_num_6(self):
        """Verify non-generic product description with generic 'GP' number
        but incorrect prefix_set is parsed correctly."""
        attr = "_product_num"
        value = " GP10; terminase "
        self.feature.set_num(attr, value, None, {"xyz"})
        self.assertEqual(self.feature._product_num, "")

    def test_set_num_7(self):
        """Verify non-generic product description with generic 'XYZ' number
        and correct prefix_set is parsed correctly."""
        attr = "_product_num"
        value = " XYZ10; terminase "
        self.feature.set_num(attr, value, None, {"xyz"})
        self.assertEqual(self.feature._product_num, "10")

    def test_set_num_8(self):
        """Verify non-generic product description with generic 'ORF' number
        is parsed correctly."""
        attr = "_product_num"
        value = " ORF10; terminase "
        self.feature.set_num(attr, value, None, None)
        self.assertEqual(self.feature._product_num, "10")

    def test_set_num_9(self):
        """Verify generic product description with generic 'GP' number
        is parsed correctly."""
        attr = "_product_num"
        value = " GP10 "
        self.feature.set_num(attr, value, None, None)
        self.assertEqual(self.feature._product_num, "10")

    def test_set_num_10(self):
        """Verify generic product description with only number
        is parsed correctly."""
        attr = "_product_num"
        value = " 10 "
        self.feature.set_num(attr, value, None, None)
        self.assertEqual(self.feature._product_num, "10")

    def test_set_num_11(self):
        """Verify correct ';' delimiter is selected."""
        attr = "_product_num"
        value = " GP10; terminase, GP11 "
        self.feature.set_num(attr, value, None, None)
        self.assertEqual(self.feature._product_num, "10")

    def test_set_num_12(self):
        """Verify correct ', delimiter is selected."""
        attr = "_product_num"
        value = " GP10; terminase, lysin, GP11 "
        self.feature.set_num(attr, value, None, None)
        self.assertEqual(self.feature._product_num, "11")

    def test_set_num_13(self):
        """Verify correct number is selected when more than one is present."""
        attr = "_product_num"
        value = " GP10; GP11; terminase "
        self.feature.set_num(attr, value, None, None)
        self.assertEqual(self.feature._product_num, "10")




    def test_set_gene_1(self):
        """Verify gene and gene_num is set correctly when gene is int."""
        value = " 10 "
        self.feature.set_gene(value, None, None)
        with self.subTest():
            self.assertEqual(self.feature.gene, "10")
        with self.subTest():
            self.assertEqual(self.feature._gene_num, "10")

    def test_set_gene_2(self):
        """Verify gene and gene_num is set correctly when gene is float."""
        value = " 10.1 "
        self.feature.set_gene(value, None, None)
        with self.subTest():
            self.assertEqual(self.feature.gene, "10.1")
        with self.subTest():
            self.assertEqual(self.feature._gene_num, "10.1")

    def test_set_gene_3(self):
        """Verify gene and gene_num is set correctly when gene contains prefix."""
        value = " GP10.1 "
        self.feature.set_gene(value, None, None)
        with self.subTest():
            self.assertEqual(self.feature.gene, "GP10.1")
        with self.subTest():
            self.assertEqual(self.feature._gene_num, "10.1")

    def test_set_gene_4(self):
        """Verify gene_num is not set when gene contains invalid prefix."""
        value = " XYZ10.1 "
        self.feature.set_gene(value, None, None)
        with self.subTest():
            self.assertEqual(self.feature.gene, "XYZ10.1")
        with self.subTest():
            self.assertEqual(self.feature._gene_num, "")

    def test_set_gene_5(self):
        """Verify gene_num is set correctly when new prefix is provided."""
        value = " XYZ10.1 "
        self.feature.set_gene(value, None, {"xyz"})
        with self.subTest():
            self.assertEqual(self.feature.gene, "XYZ10.1")
        with self.subTest():
            self.assertEqual(self.feature._gene_num, "10.1")

    def test_set_gene_6(self):
        """Verify gene_num is not set when gene does not contain a number."""
        value = " terL "
        self.feature.set_gene(value, None, None)
        with self.subTest():
            self.assertEqual(self.feature.gene, "terL")
        with self.subTest():
            self.assertEqual(self.feature._gene_num, "")




    def test_set_description_field_1(self):
        """Verify non-generic product with no number description
        is parsed correctly."""
        value = " terminase "
        self.feature.set_description_field("product", value, None, None)
        with self.subTest():
            self.assertEqual(self.feature.raw_product, "terminase")
        with self.subTest():
            self.assertEqual(self.feature.product, "terminase")
        with self.subTest():
            self.assertEqual(self.feature._product_num, "")

    def test_set_description_field_2(self):
        """Verify non-generic product description with generic 'GP' number
        is parsed correctly."""
        value = " GP10; terminase "
        self.feature.set_description_field("product", value, None, None)
        with self.subTest():
            self.assertEqual(self.feature.raw_product, "GP10; terminase")
        with self.subTest():
            self.assertEqual(self.feature.product, "GP10; terminase")
        with self.subTest():
            self.assertEqual(self.feature._product_num, "10")

    def test_set_description_field_3(self):
        """Verify non-generic product description with generic 'GP' number
        but different prefix_set is parsed correctly."""
        value = " GP10; terminase "
        self.feature.set_description_field("product", value, None, {"xyz"})
        with self.subTest():
            self.assertEqual(self.feature.raw_product, "GP10; terminase")
        with self.subTest():
            self.assertEqual(self.feature.product, "GP10; terminase")
        with self.subTest():
            self.assertEqual(self.feature._product_num, "")


    def test_set_description_field_4(self):
        """Verify non-generic product description with generic 'GP' number
        but different delimiter is parsed correctly."""
        value = " GP10; terminase "
        self.feature.set_description_field("product", value, "-", None)
        with self.subTest():
            self.assertEqual(self.feature.raw_product, "GP10; terminase")
        with self.subTest():
            self.assertEqual(self.feature.product, "GP10; terminase")
        with self.subTest():
            self.assertEqual(self.feature._product_num, "")

    def test_set_description_field_5(self):
        """Verify non-generic function description with generic 'GP' number
        is parsed correctly."""
        value = " GP10; terminase "
        self.feature.set_description_field("function", value, None, None)
        with self.subTest():
            self.assertEqual(self.feature.raw_function, "GP10; terminase")
        with self.subTest():
            self.assertEqual(self.feature.function, "GP10; terminase")
        with self.subTest():
            self.assertEqual(self.feature._function_num, "10")

    def test_set_description_field_6(self):
        """Verify non-generic note description with generic 'GP' number
        is parsed correctly."""
        value = " GP10; terminase "
        self.feature.set_description_field("note", value, None, None)
        with self.subTest():
            self.assertEqual(self.feature.raw_note, "GP10; terminase")
        with self.subTest():
            self.assertEqual(self.feature.note, "GP10; terminase")
        with self.subTest():
            self.assertEqual(self.feature._note_num, "10")

    def test_set_description_field_7(self):
        """Verify generic product description with generic 'GP' number
        is parsed correctly."""
        value = " GP10 "
        self.feature.set_description_field("product", value, None, None)
        with self.subTest():
            self.assertEqual(self.feature.raw_product, "GP10")
        with self.subTest():
            self.assertEqual(self.feature.product, "")
        with self.subTest():
            self.assertEqual(self.feature._product_num, "10")

    def test_set_description_field_8(self):
        """Verify generic product description with only number
        is parsed correctly."""
        value = " 10 "
        self.feature.set_description_field("product", value, None, None)
        with self.subTest():
            self.assertEqual(self.feature.raw_product, "10")
        with self.subTest():
            self.assertEqual(self.feature.product, "")
        with self.subTest():
            self.assertEqual(self.feature._product_num, "10")




    def test_set_description_1(self):
        """Verify product description is assigned to primary description."""
        self.feature.raw_product = "ABCD"
        self.feature.product = "EFGH"
        self.feature.set_description("product")
        with self.subTest():
            self.assertEqual(self.feature.raw_description, "ABCD")
        with self.subTest():
            self.assertEqual(self.feature.description, "EFGH")

    def test_set_description_2(self):
        """Verify function description is assigned to primary description."""
        self.feature.raw_function = "ABCD"
        self.feature.function = "EFGH"
        self.feature.set_description("function")
        with self.subTest():
            self.assertEqual(self.feature.raw_description, "ABCD")
        with self.subTest():
            self.assertEqual(self.feature.description, "EFGH")

    def test_set_description_3(self):
        """Verify note description is assigned to primary description."""
        self.feature.raw_note = "ABCD"
        self.feature.note = "EFGH"
        self.feature.set_description("note")
        with self.subTest():
            self.assertEqual(self.feature.raw_description, "ABCD")
        with self.subTest():
            self.assertEqual(self.feature.description, "EFGH")

    def test_set_description_4(self):
        """Verify no description is assigned to primary description."""
        self.feature.raw_note = "ABCD"
        self.feature.note = "EFGH"
        self.feature.set_description("invalid")
        with self.subTest():
            self.assertEqual(self.feature.raw_description, "")
        with self.subTest():
            self.assertEqual(self.feature.description, "")




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
        self.feature.raw_product = "terminase"
        self.feature.product = "terminase"
        self.feature.check_generic_data("product", eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].id, "eval_id")

    def test_check_generic_data_2(self):
        """Verify no error is produced if the function contains valid data."""
        self.feature.raw_function = "terminase"
        self.feature.function = "terminase"
        self.feature.check_generic_data("product")
        with self.subTest():
            self.assertEqual(self.feature.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.feature.evaluations[0].id)

    def test_check_generic_data_3(self):
        """Verify no error is produced if the note contains valid data."""
        self.feature.raw_note = "terminase"
        self.feature.note = "terminase"
        self.feature.check_generic_data("product")
        self.assertEqual(self.feature.evaluations[0].status, "correct")

    def test_check_generic_data_4(self):
        """Verify an error is produced if the product contains invalid data."""
        self.feature.raw_product = "gp104"
        self.feature.product = ""
        self.feature.check_generic_data("product")
        self.assertEqual(self.feature.evaluations[0].status, "error")

    def test_check_generic_data_5(self):
        """Verify no error is produced if no attribute is selected."""
        self.feature.raw_product = "gp104"
        self.feature.product = ""
        self.feature.check_generic_data()
        self.assertEqual(self.feature.evaluations[0].status, "correct")




if __name__ == '__main__':
    unittest.main()
