"""Unit tests for misc. functions that interact with
GenBank-formatted flat files."""


import unittest
from datetime import datetime
from functions import basic
from functions import flat_files
from classes import Cds
from classes import Source
from classes import Genome
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqFeature import ExactPosition, Reference
from classes import Bundle

class TestFlatFileFunctions1(unittest.TestCase):


    def setUp(self):
        self.cds = Cds.CdsFeature()
        self.source = Source.SourceFeature()
        self.genome = Genome.Genome()

    def test_parse_coordinates_1(self):
        """Verify non-compound location is parsed correctly."""

        seqfeature = SeqFeature(FeatureLocation( \
            ExactPosition(2), ExactPosition(10)), \
            type = "CDS", \
            strand = 1)

        # output_left, output_right, parts, output_eval = \
        #     flat_files.parse_coordinates(seqfeature)
        output_left, output_right, parts = \
            flat_files.parse_coordinates(seqfeature)

        exp_left = 2
        exp_right = 10
        exp_parts = 1

        with self.subTest():
            self.assertEqual(output_left, exp_left)
        with self.subTest():
            self.assertEqual(output_right, exp_right)
        with self.subTest():
            self.assertEqual(parts, exp_parts)
        # with self.subTest():
        #     self.assertIsNone(output_eval)


    def test_parse_coordinates_2(self):
        """Verify 1 strand 2-part compound location is parsed correctly."""

        seqfeature = SeqFeature(CompoundLocation( \
            [FeatureLocation( \
                ExactPosition(2), ExactPosition(10)), \
            FeatureLocation( \
                ExactPosition(8), ExactPosition(20))]), \
            type = "CDS", \
            strand = 1)

        # output_left, output_right, parts, output_eval = \
        #     flat_files.parse_coordinates(seqfeature)
        output_left, output_right, parts = \
            flat_files.parse_coordinates(seqfeature)

        exp_left = 2
        exp_right = 20
        exp_parts = 2

        with self.subTest():
            self.assertEqual(output_left, exp_left)
        with self.subTest():
            self.assertEqual(output_right, exp_right)
        with self.subTest():
            self.assertEqual(parts, exp_parts)
        # with self.subTest():
        #     self.assertIsNone(output_eval)


    def test_parse_coordinates_3(self):
        """Verify -1 strand 2-part compound location is parsed correctly."""

        seqfeature = SeqFeature(CompoundLocation( \
            [FeatureLocation( \
                ExactPosition(2), ExactPosition(10)), \
            FeatureLocation( \
                ExactPosition(8), ExactPosition(20))]), \
            type = "CDS", \
            strand = -1)

        # output_left, output_right, parts, output_eval = \
        #     flat_files.parse_coordinates(seqfeature)
        output_left, output_right, parts = \
            flat_files.parse_coordinates(seqfeature)

        exp_left = 8
        exp_right = 10
        exp_parts = 2

        with self.subTest():
            self.assertEqual(output_left, exp_left)
        with self.subTest():
            self.assertEqual(output_right, exp_right)
        with self.subTest():
            self.assertEqual(parts, exp_parts)
        # with self.subTest():
        #     self.assertIsNone(output_eval)


    def test_parse_coordinates_4(self):
        """Verify undefined strand 2-part compound location is not parsed."""

        seqfeature = SeqFeature(CompoundLocation( \
            [FeatureLocation( \
                ExactPosition(2), ExactPosition(10)), \
            FeatureLocation( \
                ExactPosition(8), ExactPosition(20))]), \
            type = "CDS", \
            strand = None)

        # output_left, output_right, parts, output_eval = \
        #     flat_files.parse_coordinates(seqfeature)
        output_left, output_right, parts = \
            flat_files.parse_coordinates(seqfeature)

        exp_left = -1
        exp_right = -1
        exp_parts = 0

        with self.subTest():
            self.assertEqual(output_left, exp_left)
        with self.subTest():
            self.assertEqual(output_right, exp_right)
        with self.subTest():
            self.assertEqual(parts, exp_parts)
        # with self.subTest():
        #     self.assertIsNotNone(output_eval)


    def test_parse_coordinates_5(self):
        """Verify 1 strand 3-part compound location is not parsed."""

        seqfeature = SeqFeature(CompoundLocation( \
            [FeatureLocation( \
                ExactPosition(2), ExactPosition(10)), \
            FeatureLocation( \
                ExactPosition(8), ExactPosition(20)), \
            FeatureLocation( \
                ExactPosition(30), ExactPosition(50))]), \
            type = "CDS", \
            strand = 1)

        # output_left, output_right, parts, output_eval = \
        #     flat_files.parse_coordinates(seqfeature)
        output_left, output_right, parts = \
            flat_files.parse_coordinates(seqfeature)

        exp_left = -1
        exp_right = -1
        exp_parts = 3

        with self.subTest():
            self.assertEqual(output_left, exp_left)
        with self.subTest():
            self.assertEqual(output_right, exp_right)
        with self.subTest():
            self.assertEqual(parts, exp_parts)
        # with self.subTest():
        #     self.assertIsNotNone(output_eval)


    def test_parse_coordinates_6(self):
        """Verify location of invalid data type is not parsed."""

        seqfeature = SeqFeature(None, type = "CDS", strand = None)

        # output_left, output_right, parts, output_eval = \
        #     flat_files.parse_coordinates(seqfeature)
        output_left, output_right, parts = \
            flat_files.parse_coordinates(seqfeature)

        exp_left = -1
        exp_right = -1
        exp_parts = 0

        with self.subTest():
            self.assertEqual(output_left, exp_left)
        with self.subTest():
            self.assertEqual(output_right, exp_right)
        with self.subTest():
            self.assertEqual(parts, exp_parts)
        # with self.subTest():
        #     self.assertIsNotNone(output_eval)






    def test_parse_cds_feature_1(self):
        """Verify CDS features is parsed."""
        qualifier_dict = {"locus_tag": ["SEA_L5_1"], \
                            "translation": ["ABCDE"], \
                            "transl_table": ["11"], \
                            "product": [" unknown "], \
                            "function": [" hypothetical protein "], \
                            "note": [" gp5 "], \
                            "gene": ["1"]}

        seqfeature = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1, \
                    qualifiers = qualifier_dict)

        # eval_result = flat_files.parse_cds_feature(self.cds, seqfeature)
        flat_files.parse_cds_feature(self.cds, seqfeature)

        with self.subTest():
            self.assertEqual(self.cds.type, "")
        with self.subTest():
            self.assertEqual(self.cds.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(self.cds.strand, "F")
        with self.subTest():
            self.assertEqual(self.cds.left_boundary, 2)
        with self.subTest():
            self.assertEqual(self.cds.right_boundary, 10)
        with self.subTest():
            self.assertEqual(self.cds.compound_parts, 1)
        # with self.subTest():
        #     self.assertIsNone(eval_result)
        with self.subTest():
            self.assertEqual(self.cds.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(self.cds.translation, "ABCDE")
        with self.subTest():
            self.assertEqual(self.cds._translation_length, 5)
        with self.subTest():
            self.assertEqual(self.cds._nucleotide_length, 8)
        with self.subTest():
            self.assertEqual(self.cds.translation_table, "11")
        with self.subTest():
            self.assertEqual(self.cds.product_description, "unknown")
        with self.subTest():
            self.assertEqual(self.cds.processed_product_description, "")
        with self.subTest():
            self.assertEqual(self.cds.function_description, \
                "hypothetical protein")
        with self.subTest():
            self.assertEqual(self.cds.processed_function_description, "")
        with self.subTest():
            self.assertEqual(self.cds.note_description, "gp5")
        with self.subTest():
            self.assertEqual(self.cds.processed_note_description, "")
        with self.subTest():
            self.assertEqual(self.cds.gene_number, "1")
        with self.subTest():
            self.assertTrue(isinstance(self.cds.seqfeature, SeqFeature))




    def test_parse_cds_feature_2(self):
        """Verify CDS features is parsed with no locus tag."""
        qualifier_dict = {"locus_tag_x": ["SEA_L5_1"], \
                            "translation": ["ABCDE"], \
                            "transl_table": ["11"], \
                            "product": [" unknown "], \
                            "function": [" hypothetical protein "], \
                            "note": [" gp5 "], \
                            "gene": ["1"]}

        seqfeature = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1, \
                    qualifiers = qualifier_dict)

        # eval_result = flat_files.parse_cds_feature(self.cds, seqfeature)
        flat_files.parse_cds_feature(self.cds, seqfeature)

        with self.subTest():
            self.assertEqual(self.cds.type, "")
        with self.subTest():
            self.assertEqual(self.cds.locus_tag, "")
        with self.subTest():
            self.assertEqual(self.cds.strand, "F")
        with self.subTest():
            self.assertEqual(self.cds.left_boundary, 2)
        with self.subTest():
            self.assertEqual(self.cds.right_boundary, 10)
        with self.subTest():
            self.assertEqual(self.cds.compound_parts, 1)
        # with self.subTest():
        #     self.assertIsNone(eval_result)
        with self.subTest():
            self.assertEqual(self.cds.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(self.cds.translation, "ABCDE")
        with self.subTest():
            self.assertEqual(self.cds._translation_length, 5)
        with self.subTest():
            self.assertEqual(self.cds._nucleotide_length, 8)
        with self.subTest():
            self.assertEqual(self.cds.translation_table, "11")
        with self.subTest():
            self.assertEqual(self.cds.product_description, "unknown")
        with self.subTest():
            self.assertEqual(self.cds.processed_product_description, "")
        with self.subTest():
            self.assertEqual(self.cds.function_description, \
                "hypothetical protein")
        with self.subTest():
            self.assertEqual(self.cds.processed_function_description, "")
        with self.subTest():
            self.assertEqual(self.cds.note_description, "gp5")
        with self.subTest():
            self.assertEqual(self.cds.processed_note_description, "")
        with self.subTest():
            self.assertEqual(self.cds.gene_number, "1")
        with self.subTest():
            self.assertTrue(isinstance(self.cds.seqfeature, SeqFeature))


    def test_parse_cds_feature_3(self):
        """Verify CDS features is parsed with problematic coordinates."""
        qualifier_dict = {"locus_tag": ["SEA_L5_1"], \
                            "translation": ["ABCDE"], \
                            "transl_table": ["11"], \
                            "product": [" unknown "], \
                            "function": [" hypothetical protein "], \
                            "note": [" gp5 "], \
                            "gene": ["1"]}

        seqfeature = SeqFeature(CompoundLocation( \
            [FeatureLocation( \
                ExactPosition(2), ExactPosition(10)), \
            FeatureLocation( \
                ExactPosition(8), ExactPosition(20)), \
            FeatureLocation( \
                ExactPosition(30), ExactPosition(50))]), \
            type = "CDS", \
            strand = 1, \
            qualifiers = qualifier_dict)

        # eval_result = flat_files.parse_cds_feature(self.cds, seqfeature)
        flat_files.parse_cds_feature(self.cds, seqfeature)

        with self.subTest():
            self.assertEqual(self.cds.type, "")
        with self.subTest():
            self.assertEqual(self.cds.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(self.cds.strand, "F")
        with self.subTest():
            self.assertEqual(self.cds.left_boundary, -1)
        with self.subTest():
            self.assertEqual(self.cds.right_boundary, -1)
        with self.subTest():
            self.assertEqual(self.cds.compound_parts, 3)
        # with self.subTest():
        #     self.assertIsNotNone(eval_result)
        with self.subTest():
            self.assertEqual(self.cds.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(self.cds.translation, "ABCDE")
        with self.subTest():
            self.assertEqual(self.cds._translation_length, 5)
        with self.subTest():
            self.assertEqual(self.cds._nucleotide_length, 0)
        with self.subTest():
            self.assertEqual(self.cds.translation_table, "11")
        with self.subTest():
            self.assertEqual(self.cds.product_description, "unknown")
        with self.subTest():
            self.assertEqual(self.cds.processed_product_description, "")
        with self.subTest():
            self.assertEqual(self.cds.function_description, \
                "hypothetical protein")
        with self.subTest():
            self.assertEqual(self.cds.processed_function_description, "")
        with self.subTest():
            self.assertEqual(self.cds.note_description, "gp5")
        with self.subTest():
            self.assertEqual(self.cds.processed_note_description, "")
        with self.subTest():
            self.assertEqual(self.cds.gene_number, "1")
        with self.subTest():
            self.assertTrue(isinstance(self.cds.seqfeature, SeqFeature))


    def test_parse_cds_feature_4(self):
        """Verify CDS features is parsed with no translation."""
        qualifier_dict = {"locus_tag": ["SEA_L5_1"], \
                            "translation_x": ["ABCDE"], \
                            "transl_table": ["11"], \
                            "product": [" unknown "], \
                            "function": [" hypothetical protein "], \
                            "note": [" gp5 "], \
                            "gene": ["1"]}

        seqfeature = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1, \
                    qualifiers = qualifier_dict)

        # eval_result = flat_files.parse_cds_feature(self.cds, seqfeature)
        flat_files.parse_cds_feature(self.cds, seqfeature)

        with self.subTest():
            self.assertEqual(self.cds.type, "")
        with self.subTest():
            self.assertEqual(self.cds.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(self.cds.strand, "F")
        with self.subTest():
            self.assertEqual(self.cds.left_boundary, 2)
        with self.subTest():
            self.assertEqual(self.cds.right_boundary, 10)
        with self.subTest():
            self.assertEqual(self.cds.compound_parts, 1)
        # with self.subTest():
        #     self.assertIsNone(eval_result)
        with self.subTest():
            self.assertEqual(self.cds.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(self.cds.translation, "")
        with self.subTest():
            self.assertEqual(self.cds._translation_length, 0)
        with self.subTest():
            self.assertEqual(self.cds._nucleotide_length, 8)
        with self.subTest():
            self.assertEqual(self.cds.translation_table, "11")
        with self.subTest():
            self.assertEqual(self.cds.product_description, "unknown")
        with self.subTest():
            self.assertEqual(self.cds.processed_product_description, "")
        with self.subTest():
            self.assertEqual(self.cds.function_description, \
                "hypothetical protein")
        with self.subTest():
            self.assertEqual(self.cds.processed_function_description, "")
        with self.subTest():
            self.assertEqual(self.cds.note_description, "gp5")
        with self.subTest():
            self.assertEqual(self.cds.processed_note_description, "")
        with self.subTest():
            self.assertEqual(self.cds.gene_number, "1")
        with self.subTest():
            self.assertTrue(isinstance(self.cds.seqfeature, SeqFeature))


    def test_parse_cds_feature_5(self):
        """Verify CDS features is parsed with no translation table."""
        qualifier_dict = {"locus_tag": ["SEA_L5_1"], \
                            "translation": ["ABCDE"], \
                            "transl_table_x": ["11"], \
                            "product": [" unknown "], \
                            "function": [" hypothetical protein "], \
                            "note": [" gp5 "], \
                            "gene": ["1"]}

        seqfeature = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1, \
                    qualifiers = qualifier_dict)

        # eval_result = flat_files.parse_cds_feature(self.cds, seqfeature)
        flat_files.parse_cds_feature(self.cds, seqfeature)

        with self.subTest():
            self.assertEqual(self.cds.type, "")
        with self.subTest():
            self.assertEqual(self.cds.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(self.cds.strand, "F")
        with self.subTest():
            self.assertEqual(self.cds.left_boundary, 2)
        with self.subTest():
            self.assertEqual(self.cds.right_boundary, 10)
        with self.subTest():
            self.assertEqual(self.cds.compound_parts, 1)
        # with self.subTest():
        #     self.assertIsNone(eval_result)
        with self.subTest():
            self.assertEqual(self.cds.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(self.cds.translation, "ABCDE")
        with self.subTest():
            self.assertEqual(self.cds._translation_length, 5)
        with self.subTest():
            self.assertEqual(self.cds._nucleotide_length, 8)
        with self.subTest():
            self.assertEqual(self.cds.translation_table, "")
        with self.subTest():
            self.assertEqual(self.cds.product_description, "unknown")
        with self.subTest():
            self.assertEqual(self.cds.processed_product_description, "")
        with self.subTest():
            self.assertEqual(self.cds.function_description, \
                "hypothetical protein")
        with self.subTest():
            self.assertEqual(self.cds.processed_function_description, "")
        with self.subTest():
            self.assertEqual(self.cds.note_description, "gp5")
        with self.subTest():
            self.assertEqual(self.cds.processed_note_description, "")
        with self.subTest():
            self.assertEqual(self.cds.gene_number, "1")
        with self.subTest():
            self.assertTrue(isinstance(self.cds.seqfeature, SeqFeature))


    def test_parse_cds_feature_6(self):
        """Verify CDS features is parsed with no product."""
        qualifier_dict = {"locus_tag": ["SEA_L5_1"], \
                            "translation": ["ABCDE"], \
                            "transl_table": ["11"], \
                            "product_x": [" unknown "], \
                            "function": [" hypothetical protein "], \
                            "note": [" gp5 "], \
                            "gene": ["1"]}

        seqfeature = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1, \
                    qualifiers = qualifier_dict)

        # eval_result = flat_files.parse_cds_feature(self.cds, seqfeature)
        flat_files.parse_cds_feature(self.cds, seqfeature)

        with self.subTest():
            self.assertEqual(self.cds.type, "")
        with self.subTest():
            self.assertEqual(self.cds.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(self.cds.strand, "F")
        with self.subTest():
            self.assertEqual(self.cds.left_boundary, 2)
        with self.subTest():
            self.assertEqual(self.cds.right_boundary, 10)
        with self.subTest():
            self.assertEqual(self.cds.compound_parts, 1)
        # with self.subTest():
        #     self.assertIsNone(eval_result)
        with self.subTest():
            self.assertEqual(self.cds.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(self.cds.translation, "ABCDE")
        with self.subTest():
            self.assertEqual(self.cds._translation_length, 5)
        with self.subTest():
            self.assertEqual(self.cds._nucleotide_length, 8)
        with self.subTest():
            self.assertEqual(self.cds.translation_table, "11")
        with self.subTest():
            self.assertEqual(self.cds.product_description, "")
        with self.subTest():
            self.assertEqual(self.cds.processed_product_description, "")
        with self.subTest():
            self.assertEqual(self.cds.function_description, \
                "hypothetical protein")
        with self.subTest():
            self.assertEqual(self.cds.processed_function_description, "")
        with self.subTest():
            self.assertEqual(self.cds.note_description, "gp5")
        with self.subTest():
            self.assertEqual(self.cds.processed_note_description, "")
        with self.subTest():
            self.assertEqual(self.cds.gene_number, "1")
        with self.subTest():
            self.assertTrue(isinstance(self.cds.seqfeature, SeqFeature))


    def test_parse_cds_feature_7(self):
        """Verify CDS features is parsed with no function."""
        qualifier_dict = {"locus_tag": ["SEA_L5_1"], \
                            "translation": ["ABCDE"], \
                            "transl_table": ["11"], \
                            "product": [" unknown "], \
                            "function_x": [" hypothetical protein "], \
                            "note": [" gp5 "], \
                            "gene": ["1"]}

        seqfeature = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1, \
                    qualifiers = qualifier_dict)

        # eval_result = flat_files.parse_cds_feature(self.cds, seqfeature)
        flat_files.parse_cds_feature(self.cds, seqfeature)

        with self.subTest():
            self.assertEqual(self.cds.type, "")
        with self.subTest():
            self.assertEqual(self.cds.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(self.cds.strand, "F")
        with self.subTest():
            self.assertEqual(self.cds.left_boundary, 2)
        with self.subTest():
            self.assertEqual(self.cds.right_boundary, 10)
        with self.subTest():
            self.assertEqual(self.cds.compound_parts, 1)
        # with self.subTest():
        #     self.assertIsNone(eval_result)
        with self.subTest():
            self.assertEqual(self.cds.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(self.cds.translation, "ABCDE")
        with self.subTest():
            self.assertEqual(self.cds._translation_length, 5)
        with self.subTest():
            self.assertEqual(self.cds._nucleotide_length, 8)
        with self.subTest():
            self.assertEqual(self.cds.translation_table, "11")
        with self.subTest():
            self.assertEqual(self.cds.product_description, "unknown")
        with self.subTest():
            self.assertEqual(self.cds.processed_product_description, "")
        with self.subTest():
            self.assertEqual(self.cds.function_description, \
                "")
        with self.subTest():
            self.assertEqual(self.cds.processed_function_description, "")
        with self.subTest():
            self.assertEqual(self.cds.note_description, "gp5")
        with self.subTest():
            self.assertEqual(self.cds.processed_note_description, "")
        with self.subTest():
            self.assertEqual(self.cds.gene_number, "1")
        with self.subTest():
            self.assertTrue(isinstance(self.cds.seqfeature, SeqFeature))


    def test_parse_cds_feature_8(self):
        """Verify CDS features is parsed with no note."""
        qualifier_dict = {"locus_tag": ["SEA_L5_1"], \
                            "translation": ["ABCDE"], \
                            "transl_table": ["11"], \
                            "product": [" unknown "], \
                            "function": [" hypothetical protein "], \
                            "note_x": [" gp5 "], \
                            "gene": ["1"]}

        seqfeature = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1, \
                    qualifiers = qualifier_dict)

        # eval_result = flat_files.parse_cds_feature(self.cds, seqfeature)
        flat_files.parse_cds_feature(self.cds, seqfeature)

        with self.subTest():
            self.assertEqual(self.cds.type, "")
        with self.subTest():
            self.assertEqual(self.cds.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(self.cds.strand, "F")
        with self.subTest():
            self.assertEqual(self.cds.left_boundary, 2)
        with self.subTest():
            self.assertEqual(self.cds.right_boundary, 10)
        with self.subTest():
            self.assertEqual(self.cds.compound_parts, 1)
        # with self.subTest():
        #     self.assertIsNone(eval_result)
        with self.subTest():
            self.assertEqual(self.cds.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(self.cds.translation, "ABCDE")
        with self.subTest():
            self.assertEqual(self.cds._translation_length, 5)
        with self.subTest():
            self.assertEqual(self.cds._nucleotide_length, 8)
        with self.subTest():
            self.assertEqual(self.cds.translation_table, "11")
        with self.subTest():
            self.assertEqual(self.cds.product_description, "unknown")
        with self.subTest():
            self.assertEqual(self.cds.processed_product_description, "")
        with self.subTest():
            self.assertEqual(self.cds.function_description, \
                "hypothetical protein")
        with self.subTest():
            self.assertEqual(self.cds.processed_function_description, "")
        with self.subTest():
            self.assertEqual(self.cds.note_description, "")
        with self.subTest():
            self.assertEqual(self.cds.processed_note_description, "")
        with self.subTest():
            self.assertEqual(self.cds.gene_number, "1")
        with self.subTest():
            self.assertTrue(isinstance(self.cds.seqfeature, SeqFeature))


    def test_parse_cds_feature_9(self):
        """Verify CDS features is parsed with no gene."""
        qualifier_dict = {"locus_tag": ["SEA_L5_1"], \
                            "translation": ["ABCDE"], \
                            "transl_table": ["11"], \
                            "product": [" unknown "], \
                            "function": [" hypothetical protein "], \
                            "note": [" gp5 "], \
                            "gene_x": ["1"]}

        seqfeature = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1, \
                    qualifiers = qualifier_dict)

        # eval_result = flat_files.parse_cds_feature(self.cds, seqfeature)
        flat_files.parse_cds_feature(self.cds, seqfeature)

        with self.subTest():
            self.assertEqual(self.cds.type, "")
        with self.subTest():
            self.assertEqual(self.cds.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(self.cds.strand, "F")
        with self.subTest():
            self.assertEqual(self.cds.left_boundary, 2)
        with self.subTest():
            self.assertEqual(self.cds.right_boundary, 10)
        with self.subTest():
            self.assertEqual(self.cds.compound_parts, 1)
        # with self.subTest():
        #     self.assertIsNone(eval_result)
        with self.subTest():
            self.assertEqual(self.cds.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(self.cds.translation, "ABCDE")
        with self.subTest():
            self.assertEqual(self.cds._translation_length, 5)
        with self.subTest():
            self.assertEqual(self.cds._nucleotide_length, 8)
        with self.subTest():
            self.assertEqual(self.cds.translation_table, "11")
        with self.subTest():
            self.assertEqual(self.cds.product_description, "unknown")
        with self.subTest():
            self.assertEqual(self.cds.processed_product_description, "")
        with self.subTest():
            self.assertEqual(self.cds.function_description, \
                "hypothetical protein")
        with self.subTest():
            self.assertEqual(self.cds.processed_function_description, "")
        with self.subTest():
            self.assertEqual(self.cds.note_description, "gp5")
        with self.subTest():
            self.assertEqual(self.cds.processed_note_description, "")
        with self.subTest():
            self.assertEqual(self.cds.gene_number, "")
        with self.subTest():
            self.assertTrue(isinstance(self.cds.seqfeature, SeqFeature))



    def test_create_feature_dictionary_1(self):
        """Verify feature dictionary is constructed correctly with
        one feature."""

        feature_list = [SeqFeature(type = "CDS")]

        feature_dict = flat_files.create_feature_dictionary(feature_list)

        with self.subTest():
            self.assertEqual(len(feature_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(len(feature_dict["CDS"]), 1)


    def test_create_feature_dictionary_2(self):
        """Verify feature dictionary is constructed correctly with
        no features."""

        feature_list = []

        feature_dict = flat_files.create_feature_dictionary(feature_list)

        with self.subTest():
            self.assertEqual(len(feature_dict.keys()), 0)


    def test_create_feature_dictionary_3(self):
        """Verify feature dictionary is constructed correctly with
        several different feature types."""

        feature_list = [ \
            SeqFeature(type = "CDS"), \
            SeqFeature(type = "CDS"), \
            SeqFeature(type = "tRNA"), \
            SeqFeature(type = "tmRNA"), \
            SeqFeature(type = "other"), \
            SeqFeature(type = "gene")]

        feature_dict = flat_files.create_feature_dictionary(feature_list)

        with self.subTest():
            self.assertEqual(len(feature_dict.keys()), 5)
        with self.subTest():
            self.assertEqual(len(feature_dict["CDS"]), 2)
        with self.subTest():
            self.assertEqual(len(feature_dict["tRNA"]), 1)
        with self.subTest():
            self.assertEqual(len(feature_dict["tmRNA"]), 1)
        with self.subTest():
            self.assertEqual(len(feature_dict["other"]), 1)
        with self.subTest():
            self.assertEqual(len(feature_dict["gene"]), 1)













    def test_create_cds_objects_1(self):
        """Verify cds objects list is constructed from empty Biopython
        CDS feature list."""
        biopython_feature_list = []
        cds_object_list = flat_files.create_cds_objects(biopython_feature_list)
        with self.subTest():
            self.assertEqual(len(cds_object_list), 0)

    def test_create_cds_objects_2(self):
        """Verify cds objects list is constructed from list of one Biopython
        CDS features."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        biopython_feature_list = [seqfeature1]

        cds_object_list = flat_files.create_cds_objects(biopython_feature_list)
        with self.subTest():
            self.assertEqual(len(cds_object_list), 1)



    def test_create_cds_objects_3(self):
        """Verify cds objects list is constructed from list of three Biopython
        CDS features."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)


        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(80)), \
                    type = "CDS", \
                    strand = -1)


        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(5), ExactPosition(6)), \
                    type = "CDS", \
                    strand = 1)

        biopython_feature_list = [seqfeature1, seqfeature2, seqfeature3]

        cds_object_list = flat_files.create_cds_objects(biopython_feature_list)
        with self.subTest():
            self.assertEqual(len(cds_object_list), 3)





    #
    # def test_create_cds_objects_4(self):
    #     """Verify cds objects list is constructed from list of two Biopython
    #     CDS features when a third has an error."""
    #
    #     seqfeature1 = SeqFeature(FeatureLocation( \
    #                 ExactPosition(2), ExactPosition(10)), \
    #                 type = "CDS", \
    #                 strand = 1)
    #
    #
    #     seqfeature2 = SeqFeature(FeatureLocation( \
    #                 ExactPosition(50), ExactPosition(80)), \
    #                 type = "CDS", \
    #                 strand = -1)
    #
    #     seqfeature3 = SeqFeature(FeatureLocation( \
    #                 ExactPosition(5), ExactPosition(6)), \
    #                 type = "CDS", \
    #                 strand = None)
    #
    #
    #     biopython_feature_list = [seqfeature1, seqfeature2, seqfeature3]
    #
    #     cds_object_list = flat_files.create_cds_objects(biopython_feature_list)
    #     self.assertEqual(len(cds_object_list), 2)














    def test_parse_source_feature_1(self):
        """Verify source feature is parsed."""

        string1 = "Mycobacterium phage Trixie"
        string2 = "Mycobacterium smegmatis"
        string3 = "Gordonia terrae"

        qualifier_dict = {"organism": [string1], \
                            "host": [string2], \
                            "lab_host": [string3]}

        seqfeature = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "source", \
                    qualifiers = qualifier_dict)

        flat_files.parse_source_feature(self.source, seqfeature)

        with self.subTest():
            self.assertEqual(self.source.organism, string1)
        with self.subTest():
            self.assertEqual(self.source.host, string2)
        with self.subTest():
            self.assertEqual(self.source.lab_host, string3)


    def test_parse_source_feature_2(self):
        """Verify source feature is parsed with no organism qualifier."""

        string1 = "Mycobacterium phage Trixie"
        string2 = "Mycobacterium smegmatis"
        string3 = "Gordonia terrae"

        qualifier_dict = {"organism_x": [string1], \
                            "host": [string2], \
                            "lab_host": [string3]}

        seqfeature = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "source", \
                    qualifiers = qualifier_dict)

        flat_files.parse_source_feature(self.source, seqfeature)

        with self.subTest():
            self.assertEqual(self.source.organism, "")
        with self.subTest():
            self.assertEqual(self.source.host, string2)
        with self.subTest():
            self.assertEqual(self.source.lab_host, string3)


    def test_parse_source_feature_3(self):
        """Verify source feature is parsed with no host qualifier."""

        string1 = "Mycobacterium phage Trixie"
        string2 = "Mycobacterium smegmatis"
        string3 = "Gordonia terrae"

        qualifier_dict = {"organism": [string1], \
                            "host_x": [string2], \
                            "lab_host": [string3]}

        seqfeature = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "source", \
                    qualifiers = qualifier_dict)

        flat_files.parse_source_feature(self.source, seqfeature)

        with self.subTest():
            self.assertEqual(self.source.organism, string1)
        with self.subTest():
            self.assertEqual(self.source.host, "")
        with self.subTest():
            self.assertEqual(self.source.lab_host, string3)


    def test_parse_source_feature_4(self):
        """Verify source feature is parsed with no lab_host qualifier."""

        string1 = "Mycobacterium phage Trixie"
        string2 = "Mycobacterium smegmatis"
        string3 = "Gordonia terrae"

        qualifier_dict = {"organism": [string1], \
                            "host": [string2], \
                            "lab_host_x": [string3]}

        seqfeature = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "source", \
                    qualifiers = qualifier_dict)

        flat_files.parse_source_feature(self.source, seqfeature)

        with self.subTest():
            self.assertEqual(self.source.organism, string1)
        with self.subTest():
            self.assertEqual(self.source.host, string2)
        with self.subTest():
            self.assertEqual(self.source.lab_host, "")






    def test_create_source_objects_1(self):
        """Verify source objects list is constructed from empty Biopython
        source feature list."""
        biopython_feature_list = []
        source_object_list = \
            flat_files.create_source_objects(biopython_feature_list)
        self.assertEqual(len(source_object_list), 0)




    def test_create_source_objects_2(self):
        """Verify source objects list is constructed from list of
        one Biopython source feature."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "source", \
                    strand = 1)

        biopython_feature_list = [seqfeature1]

        source_object_list = \
            flat_files.create_source_objects(biopython_feature_list)
        self.assertEqual(len(source_object_list), 1)



    def test_create_source_objects_3(self):
        """Verify source objects list is constructed from list of
        three Biopython source features."""


        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "source", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(80)), \
                    type = "source", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(5), ExactPosition(6)), \
                    type = "source", \
                    strand = 1)


        biopython_feature_list = [seqfeature1, seqfeature2, seqfeature3]

        source_object_list = \
            flat_files.create_source_objects(biopython_feature_list)
        self.assertEqual(len(source_object_list), 3)




    def test_parse_flat_file_data_1(self):
        """Verify retrieved flat file data is parsed correctly."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature3, seqfeature4]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )


        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.id,"KatherineG")
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")



    def test_parse_flat_file_data_2(self):
        """Verify retrieved flat file data is parsed correctly with no
        record name."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature3, seqfeature4]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )



        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_3(self):
        """Verify retrieved flat file data is parsed correctly with no
        record organism."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature3, seqfeature4]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )

        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, "")
        with self.subTest():
            self.assertEqual(self.genome._organism_name, "")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, "")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_4(self):
        """Verify retrieved flat file data is parsed correctly with no
        record id."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature3, seqfeature4]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )



        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.id, "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_5(self):
        """Verify retrieved flat file data is parsed correctly with no
        accession."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature3, seqfeature4]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = { \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )

        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_6(self):
        """Verify retrieved flat file data is parsed correctly with more
        than one accession."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature3, seqfeature4]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 ", "TUV456"], \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )

        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_7(self):
        """Verify retrieved flat file data is parsed correctly with no
        record description."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature3, seqfeature4]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            features = feature_list, \
                            )

        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, "")
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_8(self):
        """Verify retrieved flat file data is parsed correctly with no
        record source."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature3, seqfeature4]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )

        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, "")
        with self.subTest():
            self.assertEqual(self.genome._source_name, "")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, "")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_9(self):
        """Verify retrieved flat file data is parsed correctly with no
        references."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature3, seqfeature4]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )

        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_10(self):
        """Verify retrieved flat file data is parsed correctly with
        references that contain no authors."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature3, seqfeature4]


        reference1 = Reference()
        reference2 = Reference()
        reference3 = Reference()

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )

        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_11(self):
        """Verify retrieved flat file data is parsed correctly with an
        empty sequence."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature3, seqfeature4]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq(""), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )

        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "")
        with self.subTest():
            self.assertEqual(self.genome._length, 0)
        with self.subTest():
            self.assertEqual(self.genome._gc, -1)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_12(self):
        """Verify retrieved flat file data is parsed correctly with no
        CDS features."""

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        feature_list = [seqfeature2, seqfeature3]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )

        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 0)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 0)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_13(self):
        """Verify retrieved flat file data is parsed correctly with no
        source features."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature4]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )

        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 0)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 0)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_14(self):
        """Verify retrieved flat file data is parsed correctly with no
        tRNA features."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature3, seqfeature4]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )

        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 0)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 0)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_15(self):
        """Verify retrieved flat file data is parsed correctly with no
        features."""


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            )

        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 0)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 0)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 0)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 0)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 0)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 0)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")



    def test_parse_flat_file_data_16(self):
        """Verify retrieved flat file data is parsed correctly with no
        filepath provided."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature3, seqfeature4]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        date = "23-JAN-2014"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            "date": date}

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )

        flat_files.parse_flat_file_data(self.genome, record)

        exp_date = datetime.strptime(date,'%d-%b-%Y')

        with self.subTest():
            self.assertEqual(self.genome.filename, "")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_17(self):
        """Verify retrieved flat file data is parsed correctly with no
        date provided."""

        seqfeature1 = SeqFeature(FeatureLocation( \
                    ExactPosition(2), ExactPosition(10)), \
                    type = "CDS", \
                    strand = 1)

        seqfeature2 = SeqFeature(FeatureLocation( \
                    ExactPosition(50), ExactPosition(55)), \
                    type = "tRNA", \
                    strand = 1)

        seqfeature3 = SeqFeature(FeatureLocation( \
                    ExactPosition(20), ExactPosition(30)), \
                    type = "source", \
                    strand = 1)

        seqfeature4 = SeqFeature(FeatureLocation( \
                    ExactPosition(100), ExactPosition(1000)), \
                    type = "CDS", \
                    strand = 1)

        feature_list = [seqfeature1, seqfeature2, seqfeature3, seqfeature4]


        reference1 = Reference()
        reference1.authors = "Jane"

        reference2 = Reference()
        reference2.authors = "Doe"

        reference3 = Reference()
        reference3.authors = "Smith"

        refs_list = [reference1, reference2, reference3]

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            "references": refs_list, \
                            }

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            features = feature_list, \
                            )


        filepath = "/path/to/file/Phage_ZZZ.gb"
        flat_files.parse_flat_file_data(self.genome, record, filepath)

        exp_date = basic.convert_empty("","empty_datetime_obj")

        with self.subTest():
            self.assertEqual(self.genome.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(self.genome.name, "XYZ123")
        with self.subTest():
            self.assertEqual(self.genome.organism, organism)
        with self.subTest():
            self.assertEqual(self.genome._organism_name, \
                                "KatherineG")
        with self.subTest():
            self.assertEqual(self.genome._organism_host_genus, \
                                "Gordonia")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.description, description)
        with self.subTest():
            self.assertEqual(self.genome._description_name, \
                                "L5")
        with self.subTest():
            self.assertEqual(self.genome._description_host_genus, \
                                "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.source, source)
        with self.subTest():
            self.assertEqual(self.genome._source_name, \
                                "phiC31")
        with self.subTest():
            self.assertEqual(self.genome._source_host_genus, \
                                "Streptomyces")
        with self.subTest():
            self.assertEqual(self.genome.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(self.genome.seq, "ATGC")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)
        with self.subTest():
            self.assertEqual(self.genome._gc, 50.00)
        with self.subTest():
            self.assertEqual(self.genome.date, exp_date)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 2)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 1)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 1)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 2)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(self.genome.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome.type, "flat_file")


    def test_parse_flat_file_data_18(self):
        """Verify retrieved flat file data is parsed correctly with
        modified translation table."""

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            )
        flat_files.parse_flat_file_data(self.genome, record, \
                                        translation_table = 1)
        self.assertEqual(self.genome.translation_table, 1)





    def test_parse_flat_file_data_19(self):
        """Verify retrieved flat file data is parsed correctly with
        id field specified as non-standard field."""

        description = "Mycobacterium phage L5 complete genome"
        organism = "Gordonia phage KatherineG"
        source = "Streptomyces phage phiC31"

        annotation_dict = {"accessions": [" ABC123.1 "], \
                            "source": source, \
                            "organism": organism, \
                            }

        record = SeqRecord(seq = Seq("atgc"), \
                            id = "OPQ123.1", \
                            name = "XYZ123", \
                            annotations = annotation_dict, \
                            description = description, \
                            )

        flat_files.parse_flat_file_data(
                            self.genome,
                            record,
                            phage_id_field = "description_name")

        self.assertEqual(self.genome.id, "L5")








    def test_check_flat_file_type_1(self):
        """Verify valid file does not produce an error."""
        filepath = "/path/to/file/l5.gb"
        result = flat_files.check_flat_file_type(filepath)
        self.assertTrue(result)

    def test_check_flat_file_type_2(self):
        """Verify invalid file produces an error."""
        filepath = "/path/to/file/l5.exe"
        result = flat_files.check_flat_file_type(filepath)
        self.assertFalse(result)






### Pasted below from phagesdb
class TestFlatFileFunctions2(unittest.TestCase):


    def setUp(self):


        self.genome1 = Genome.Genome()
        self.genome1.id = "L5"
        self.genome1.cluster = "B"
        self.genome1.type = "flat_file"
        self.genome1._value_flag = False
        self.genome1.translation_table = "empty"

        self.bundle1 = Bundle.Bundle()


        self.genome2 = Genome.Genome()
        self.genome2.id = "L5"
        self.genome2.type = "add"
        self.genome2.cluster = "A"
        self.genome2.subcluster = "A2"
        self.genome2.name = "L5_Draft"
        self.genome2.host_genus = "Mycobacterium"
        self.genome2.accession = "ABC123"
        self.genome2.cluster_subcluster = "C"
        self.genome2.annotation_status = "final"
        self.genome2.annotation_author = 1
        self.genome2.annotation_qc = 2
        self.genome2.retrieve_record = 3
        self.genome2.translation_table = 11


    def test_copy_data_to_flat_file_1(self):
        """Check that a "flat_file" genome is successfully populated."""

        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.bundle1.genome_dict[self.genome2.type] = self.genome2
        flat_files.copy_data_to_flat_file(self.bundle1, "add")
        genome1 = self.bundle1.genome_dict["flat_file"]
        with self.subTest():
            self.assertFalse(genome1._value_flag)
        with self.subTest():
            self.assertEqual(genome1.cluster, "A")
        with self.subTest():
            self.assertEqual(genome1.subcluster, "A2")
        with self.subTest():
            self.assertEqual(genome1.name, "L5_Draft")
        with self.subTest():
            self.assertEqual(genome1.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(genome1.accession, "ABC123")
        with self.subTest():
            self.assertEqual(genome1.cluster_subcluster, "C")
        with self.subTest():
            self.assertEqual(genome1.annotation_status, "final")
        with self.subTest():
            self.assertEqual(genome1.annotation_author, 1)
        with self.subTest():
            self.assertEqual(genome1.annotation_qc, 2)
        with self.subTest():
            self.assertEqual(genome1.retrieve_record, 3)
        with self.subTest():
            self.assertEqual(genome1.translation_table, "empty")
        with self.subTest():
            self.assertEqual(genome1.evaluations[0].status, "correct")


    def test_copy_data_to_flat_file_2(self):
        """Check that the function can handle a missing "flat_file" genome."""

        self.bundle1.genome_dict[self.genome2.type] = self.genome2
        flat_files.copy_data_to_flat_file(self.bundle1, "add")
        self.assertEqual(len(self.bundle1.genome_pair_dict.keys()), 0)


    def test_copy_data_to_flat_file_3(self):
        """Check that a "flat_file" genome is not successfully populated
        when a second genome is absent."""

        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        flat_files.copy_data_to_flat_file(self.bundle1, "add")
        genome1 = self.bundle1.genome_dict["flat_file"]
        with self.subTest():
            self.assertTrue(genome1._value_flag)
        with self.subTest():
            self.assertEqual(genome1.evaluations[0].status, "error")


    def test_copy_data_to_flat_file_4(self):
        """Check that a "flat_file" genome is successfully populated when
        a non-standard field flag is used."""

        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.bundle1.genome_dict[self.genome2.type] = self.genome2
        flat_files.copy_data_to_flat_file(self.bundle1, "add", "empty")
        genome1 = self.bundle1.genome_dict["flat_file"]
        with self.subTest():
            self.assertFalse(genome1._value_flag)
        with self.subTest():
            self.assertEqual(genome1.cluster, "A")
        with self.subTest():
            self.assertEqual(genome1.translation_table, 11)
        with self.subTest():
            self.assertEqual(genome1.evaluations[0].status, "correct")













if __name__ == '__main__':
    unittest.main()
