"""
Unit tests for functions in phameration.py
"""

import os
import shutil
import unittest
import subprocess
import sys
from pathlib import Path

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.functions.phameration import *
from pdm_utils.functions import mysqldb_basic

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

# Standard pdm_anon user/pwd and pdm_test_db
DB = test_db_utils.DB
USER = test_db_utils.USER
PWD = test_db_utils.PWD


class TestPhamerationFunctions(unittest.TestCase):
    def setUp(self):

        # Create test database that contains data for several phages.
        test_db_utils.create_filled_test_db()

        self.alchemist = AlchemyHandler(database=DB, username=USER, password=PWD)
        self.alchemist.build_engine()
        self.engine = self.alchemist.engine
        self.temp_dir = "/tmp/pdm_utils_tests_phamerate"

    def tearDown(self):
        self.engine.dispose()
        test_db_utils.remove_db()

        run_dir = Path.cwd()
        err_file = run_dir.joinpath("error.log")
        if err_file.exists():
            print("Found leftover blastclust file... removing")
            err_file.unlink()

    def test_1_get_pham_geneids(self):
        """Verify we get back a dictionary"""
        old_phams = get_pham_geneids(self.engine)
        # old_phams should be a dict
        self.assertEqual(type(old_phams), type(dict()))

    def test_2_get_pham_colors(self):
        """Verify we get back a dictionary"""
        old_colors = get_pham_colors(self.engine)
        # old_colors should be a dict
        self.assertEqual(type(old_colors), type(dict()))

    def test_3_get_pham_geneids_and_colors(self):
        """Verify both dictionaries have the same keys"""
        old_phams = get_pham_geneids(self.engine)
        old_colors = get_pham_colors(self.engine)

        # Can't have same keys without the same number of keys...
        with self.subTest():
            self.assertEqual(len(old_phams), len(old_colors))

        # Intersection should be equal to either set of keys - check against old_phams
        with self.subTest():
            self.assertEqual(set(old_phams.keys()).intersection(set(old_colors.keys())), set(old_phams.keys()))

    def test_4_get_unphamerated_genes(self):
        """Verify we get back a set of length 0"""
        unphamerated = get_new_geneids(self.engine)
        # unphamerated should be a set
        with self.subTest():
            self.assertEqual(type(unphamerated), type(set()))
        # pdm_test_db has 0 unphamerated genes
        with self.subTest():
            self.assertEqual(len(unphamerated), 0)

    def test_5_get_geneids_and_translations(self):
        """Verify we get back a dictionary"""
        gs_to_ts = get_geneids_and_translations(self.engine)

        command = "SELECT distinct(GeneID) FROM gene"
        results = mysqldb_basic.query_dict_list(self.engine, command)

        # gs_to_ts should be a dictionary
        with self.subTest():
            self.assertEqual(type(gs_to_ts), type(dict()))
        # gs_to_ts should have the right number of geneids
        with self.subTest():
            self.assertEqual(len(gs_to_ts), len(results))

    def test_6_get_translation_groups(self):
        """Verify we get back a dictionary"""
        ts_to_gs = get_translation_groups(self.engine)

        command = "SELECT distinct(CONVERT(Translation USING utf8)) FROM gene"
        results = mysqldb_basic.query_dict_list(self.engine, command)

        # ts_to_gs should be a dictionary
        with self.subTest():
            self.assertEqual(type(ts_to_gs), type(dict()))
        # ts_to_gs should have the right number of translations
        with self.subTest():
            self.assertEqual(len(ts_to_gs), len(results))

    def test_7_refresh_tempdir_1(self):
        """Verify if no temp_dir, refresh can make one"""
        if not os.path.exists(self.temp_dir):
            refresh_tempdir(self.temp_dir)
        self.assertTrue(os.path.exists(self.temp_dir))

    def test_8_refresh_tempdir_2(self):
        """Verify if temp_dir with something, refresh makes new empty one"""
        filename = f"{self.temp_dir}/test.txt"
        if not os.path.exists(self.temp_dir):
            refresh_tempdir(self.temp_dir)
        f = open(filename, "w")
        f.write("test\n")
        f.close()

        # Our test file should now exist
        with self.subTest():
            self.assertTrue(os.path.exists(filename))

        # Refresh temp_dir
        refresh_tempdir(self.temp_dir)

        # temp_dir should now exist, but test file should not
        with self.subTest():
            self.assertTrue(os.path.exists(self.temp_dir))
        with self.subTest():
            self.assertFalse(os.path.exists(filename))

    def test_9_write_fasta(self):
        """Verify file gets written properly"""
        filename = f"{self.temp_dir}/input.fasta"

        # refresh_tempdir
        refresh_tempdir(self.temp_dir)

        # Get translations to geneid mappings
        ts_to_gs = get_translation_groups(self.engine)

        # Write fasta
        write_fasta(ts_to_gs, filename)

        # Read fasta, make sure number of lines is 2x number of unique translations
        with open(filename, "r") as fh:
            lines = fh.readlines()

        with self.subTest():
            self.assertEqual(len(lines), 2 * len(ts_to_gs))

        # all odd-index lines should map to a key in ts_to_gs
        for i in range(len(lines)):
            if i % 2 == 1:
                with self.subTest():
                    self.assertTrue(lines[i].lstrip(">").rstrip() in ts_to_gs.keys())

    def test_10_create_blastdb(self):
        """Verify blast protein database gets made"""
        filename = f"{self.temp_dir}/input.fasta"
        db_name = "blastdb"
        db_path = f"{self.temp_dir}/blastdb"

        refresh_tempdir(self.temp_dir)

        ts_to_gs = get_translation_groups(self.engine)
        write_fasta(ts_to_gs, filename)

        create_blastdb(filename, db_name, db_path)

        # Check that database files were made
        for ext in ["phr", "pin", "pog", "psd", "psi", "psq"]:
            with self.subTest():
                self.assertTrue(os.path.exists(f"{db_path}.{ext}"))

    def test_11_create_mmseqsdb(self):
        """Verify mmseqs database gets made"""
        filename = f"{self.temp_dir}/input.fasta"
        db_file = f"{self.temp_dir}/sequenceDB"

        refresh_tempdir(self.temp_dir)

        ts_to_gs = get_translation_groups(self.engine)
        write_fasta(ts_to_gs, filename)

        mmseqs_createdb(filename, db_file)

        # Check that database file was made
        self.assertTrue(os.path.exists(db_file))


def refresh_tempdir(tmpdir):
    """
    Recursively deletes tmpdir if it exists, otherwise makes it
    :param tmpdir: directory to refresh
    :return:
    """
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)


if __name__ == '__main__':
    unittest.main()
