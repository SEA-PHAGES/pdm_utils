import sys
import unittest
from pathlib import Path
from unittest.mock import patch

from pdm_utils.classes.alchemyhandler import AlchemyHandler

from pdm_utils.functions import annotation

unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))

import test_db_utils
import test_data_utils

class TestAnnotationRetrieval(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

    @classmethod
    def tearDownClass(self):
        test_db_utils.remove_db()
    
    def setUp(self):
        self.alchemist = AlchemyHandler()
        self.alchemist.username = test_db_utils.USER
        self.alchemist.password = test_db_utils.PWD
        self.alchemist.database = test_db_utils.DB
        self.alchemist.connect(ask_database=True, login_attempts=0)

    def test_get_relative_gene_1(self):
        """Verify get_relative_gene() returns GeneID string as expected."""

        rel_geneid = annotation.get_relative_gene(self.alchemist, "Trixie_CDS_2", 
                                                                        -1)

        self.assertEqual(rel_geneid, "Trixie_CDS_1")

    def test_get_relative_gene_2(self):
        """Verify get_relative_gene() returns None when expected."""

        rel_geneid = annotation.get_relative_gene(self.alchemist, "Trixie_CDS_1", 
                                                                        -1)

        self.assertEqual(rel_geneid, None)

    def test_get_relative_gene_3(self):
        """Verify get_relative_gene() raises ValueError from bad GeneID."""
        with self.assertRaises(ValueError):
            annotation.get_relative_gene(self.alchemist, "NOT A GENE", 8675309)

    def test_get_adjacent_genes_1(self):
        """Verify get_adjacent_phams() returns get_relative_gene() results."""
        adjacent_genes = annotation.get_adjacent_genes(self.alchemist, 
                                                            "Trixie_CDS_2")

        self.assertEqual(adjacent_genes[0], "Trixie_CDS_1")
        self.assertEqual(adjacent_genes[1], "Trixie_CDS_3")

    def test_get_adjacent_phams_1(self):
        """Verify get_adjacent_phams() returns expected data type."""

        adjacent_phams = annotation.get_distinct_adjacent_phams(
                                                        self.alchemist, 42006)

        self.assertTrue(isinstance(adjacent_phams, tuple))
        self.assertTrue(len(adjacent_phams) == 2)

        self.assertTrue(isinstance(adjacent_phams[0], list))
        self.assertTrue(isinstance(adjacent_phams[1], list))

        for left_pham in adjacent_phams[0]:
            with self.subTest(pham=left_pham):
                self.assertTrue(isinstance(left_pham, int))

        for right_pham in adjacent_phams[1]:
            with self.subTest(pham=right_pham):
                self.assertTrue(isinstance(right_pham, int))

    def test_get_count_annotations_in_pham_1(self):
        """Verify get_count_annotations_in_pham() returns expected data type."""

        annotation_counts = annotation.get_count_annotations_in_pham(
                                                        self.alchemist, 42006)

        self.assertTrue(isinstance(annotation_counts, dict))

        for key in annotation_counts.keys():
            with self.subTest(annotation=key):
                self.assertTrue(isinstance(annotation_counts, str))

                self.assertTrue(isinstance(annotation_counts[key], int))


if __name__ == "__main__":
    unittest.main()
