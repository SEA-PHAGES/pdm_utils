"""Integration tests for main run script."""

import unittest
from unittest.mock import patch
import argparse
import os
from pdm_utils.classes import mysqlconnectionhandler as mch
from pdm_utils import run

user = "tester"
pwd = "tester"
db = "test_db"

class TestRunFunctions(unittest.TestCase):

    def setUp(self):
        self.test_filepath1 = \
            os.path.join(os.path.dirname(__file__), \
            "test_files/test_flat_file_1.gb")

        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("-p", "--pipeline")
        self.parser.add_argument("-db", "--database")
        self.parser.add_argument("-g", "--genome_folder")
        self.parser.add_argument("-t", "--table")
        self.parser.add_argument("-f", "--filename")
        self.parser.add_argument("-tr", "--testrun")
        self.parser.add_argument("-rm", "--run_mode")
        self.parser.add_argument("-df", "--description_field")

        self.sql_handle = mch.MySQLConnectionHandler()

    # Using the patch decorator to bypass the call to import_io.
    @patch("pdm_utils.pipelines.db_import.import_genome.import_io")
    def test_run_import_1(self, import_io_mock):
        """Verify that correct args calls import_io."""
        args_list = ["--database=Actino_Draft",
                     "--genome_folder=/path/to/folder",
                     "--table=/path/to/table.csv",
                     "--run_mode=phagesdb"]
        args = self.parser.parse_args(args=args_list)
        run.run_import(args, self.sql_handle)
        self.assertTrue(import_io_mock.called)


    # Note: the import_io() patch is added because if sys.exit() is mocked,
    # then the function continues and calls import_io(), which will
    # cause an error.
    @patch("sys.exit", return_value=True)
    @patch("pdm_utils.pipelines.db_import.import_genome.import_io")
    def test_run_import_2(self, sys_exit_mock, import_io_mock):
        """Verify that import_io is not called due to
        a missing database arg."""
        args_list = ["--genome_folder=/path/to/folder",
                     "--table=/path/to/table.csv",
                     "--run_mode=phagesdb"]
        args = self.parser.parse_args(args=args_list)
        run.run_import(args, self.sql_handle)
        self.assertTrue(sys_exit_mock.called)




    # Note: the import_io() patch is added because if sys.exit() is mocked,
    # then the function continues and calls import_io(), which will
    # cause an error.
    @patch("sys.exit", return_value=True)
    @patch("pdm_utils.pipelines.db_import.import_genome.import_io")
    def test_run_import_3(self, sys_exit_mock, import_io_mock):
        """Verify that import_io is not called due to
        a missing genome_folder arg."""
        args_list = ["--database=Actino_Draft",
                     "--table=/path/to/table.csv",
                     "--run_mode=phagesdb"]
        args = self.parser.parse_args(args=args_list)
        run.run_import(args, self.sql_handle)
        self.assertTrue(sys_exit_mock.called)


    # Note: the import_io() patch is added because if sys.exit() is mocked,
    # then the function continues and calls import_io(), which will
    # cause an error.
    @patch("sys.exit", return_value=True)
    @patch("pdm_utils.pipelines.db_import.import_genome.import_io")
    def test_run_import_4(self, sys_exit_mock, import_io_mock):
        """Verify that import_io is not called due to
        a missing table arg."""
        args_list = ["--database=Actino_Draft",
                     "--genome_folder=/path/to/folder",
                     "--run_mode=phagesdb"]
        args = self.parser.parse_args(args=args_list)
        run.run_import(args, self.sql_handle)
        self.assertTrue(sys_exit_mock.called)



    # Note: the import_io() patch is added because if sys.exit() is mocked,
    # then the function continues and calls import_io(), which will
    # cause an error.
    @patch("sys.exit", return_value=True)
    @patch("pdm_utils.pipelines.db_import.import_genome.import_io")
    def test_run_import_5(self, sys_exit_mock, import_io_mock):
        """Verify that import_io is not called due to
        a missing run_mode arg."""
        args_list = ["--database=Actino_Draft",
                     "--genome_folder=/path/to/folder",
                     "--table=/path/to/table.csv"]
        args = self.parser.parse_args(args=args_list)
        run.run_import(args, self.sql_handle)
        self.assertTrue(sys_exit_mock.called)


    def test_create_parser_1(self):
        """Verify that args are returned."""
        parser = run.create_parser()
        args = parser.parse_args(["--pipeline", "import"])
        with self.subTest():
            self.assertIsInstance(parser, argparse.ArgumentParser)
        with self.subTest():
            self.assertEqual(args.pipeline, "import")



    # In progress
    # def test_get_sql_handle1(self):
    #     """."""
    #     args = self.parser.parse_args(["--pipeline", "import"])
    #     mch_mock.return_value
    #     parse_args_mock.return_value = args
    #     run.main()



    # @patch("argparse.ArgumentParser.parse_args")
    # @patch("sys.exit", return_value=True)
    # @patch("mysqlconnectionhandler.MySQLConnectionHandler")
    # def test_main_1(self, parse_args_mock, sys_exit_mock, mch_mock):
    #     """Verify that."""
    #     args = self.parser.parse_args(["--pipeline", "import"])
    #     mch_mock.return_value
    #     parse_args_mock.return_value = args
    #     run.main()












if __name__ == '__main__':
    unittest.main()
