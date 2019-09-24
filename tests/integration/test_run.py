"""Integration tests for main run script."""

import unittest
from unittest.mock import patch
import argparse
import os
import shutil
from pdm_utils.classes import mysqlconnectionhandler as mch
from pdm_utils import run

user = "tester"
pwd = "tester"
db = "test_db"

class TestRunFunctions1(unittest.TestCase):

    def setUp(self):
        self.test_filepath1 = \
            os.path.join(os.path.dirname(__file__), \
            "test_files/test_flat_file_1.gb")

        self.test_directory1 = \
            os.path.join(os.path.dirname(__file__),
            "test_wd/input_folder")
        os.mkdir(self.test_directory1)

        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("pipeline")
        self.parser.add_argument("-db", "--database")
        self.parser.add_argument("-if", "--input_folder")
        self.parser.add_argument("-it", "--import_table")
        self.parser.add_argument("-gf", "--genome_id_field")
        self.parser.add_argument("-tr", "--test_run")
        self.parser.add_argument("-rm", "--run_mode")
        self.parser.add_argument("-df", "--description_field")
        self.parser.add_argument("--sql_handle")

        self.sql_handle = mch.MySQLConnectionHandler()

        self.args1 = self.parser.parse_args(["import"])
        self.args2 = self.parser.parse_args(
                        ["import",
                        "--database", "Actino_Draft"])
        self.args3 = self.parser.parse_args([
                        "import",
                        "--input_folder", self.test_directory1])
        self.args4 = self.parser.parse_args([
                        "import",
                        "--input_folder",
                        self.test_directory1 + "/invalid_folder"])
        self.args5 = self.parser.parse_args([
                        "import",
                        "--import_table",
                        self.test_filepath1])
        self.args6 = self.parser.parse_args([
                        "import",
                        "--import_table",
                        self.test_filepath1 + "abcd"])

        self.parser2 = argparse.ArgumentParser()
        self.parser2.add_argument("-db", "--database")




    @patch("pdm_utils.run.create_parser")
    def test_get_args_1(self, create_parser_mock):
        """Verify that 'sql_handle' arg is added correctly."""
        create_parser_mock.return_value = self.parser2
        args_list = ["--database", "Actino_Draft"]
        args = run.get_args(args_list)
        self.assertIsNone(args.sql_handle)











    # Using the patch decorator to bypass the call to import_io.
    @patch("pdm_utils.pipelines.db_import.import_genome.import_io")
    def test_run_import_1(self, import_io_mock):
        """Verify that correct args calls import_io."""
        args_list = ["import",
                     "--database=Actino_Draft",
                     "--input_folder=/path/to/folder",
                     "--import_table=/path/to/table.csv",
                     "--run_mode=phagesdb"]
        args = self.parser.parse_args(args=args_list)
        run.run_import(args)
        self.assertTrue(import_io_mock.called)

    # Note: the import_io() patch is added because if sys.exit() is mocked,
    # then the function continues and calls import_io(), which will
    # cause an error.
    @patch("pdm_utils.pipelines.db_import.import_genome.import_io")
    @patch("sys.exit")
    def test_run_import_2(self, sys_exit_mock, import_io_mock):
        """Verify that import_io is not called due to
        a missing database arg."""
        args_list = ["--input_folder=/path/to/folder",
                     "--import_table=/path/to/table.csv",
                     "--run_mode=phagesdb"]
        args = self.parser.parse_args(args=args_list)
        run.run_import(args)
        self.assertTrue(sys_exit_mock.called)

    @patch("pdm_utils.pipelines.db_import.import_genome.import_io")
    @patch("sys.exit")
    def test_run_import_3(self, sys_exit_mock, import_io_mock):
        """Verify that import_io is not called due to
        a missing input_folder arg."""
        args_list = ["--database=Actino_Draft",
                     "--import_table=/path/to/table.csv",
                     "--run_mode=phagesdb"]
        args = self.parser.parse_args(args=args_list)
        run.run_import(args)
        self.assertTrue(sys_exit_mock.called)

    @patch("pdm_utils.pipelines.db_import.import_genome.import_io")
    @patch("sys.exit")
    def test_run_import_4(self, sys_exit_mock, import_io_mock):
        """Verify that import_io is not called due to
        a missing import_table arg."""
        args_list = ["--database=Actino_Draft",
                     "--input_folder=/path/to/folder",
                     "--run_mode=phagesdb"]
        args = self.parser.parse_args(args=args_list)
        run.run_import(args)
        self.assertTrue(sys_exit_mock.called)

    @patch("pdm_utils.pipelines.db_import.import_genome.import_io")
    @patch("sys.exit")
    def test_run_import_5(self, sys_exit_mock, import_io_mock):
        """Verify that import_io is not called due to
        a missing run_mode arg."""
        args_list = ["--database=Actino_Draft",
                     "--input_folder=/path/to/folder",
                     "--import_table=/path/to/table.csv"]
        args = self.parser.parse_args(args=args_list)
        run.run_import(args)
        self.assertTrue(sys_exit_mock.called)




    def test_create_parser_1(self):
        """Verify that args are returned."""
        parser = run.create_parser()
        args = parser.parse_args(["import"])
        with self.subTest():
            self.assertIsInstance(parser, argparse.ArgumentParser)
        with self.subTest():
            self.assertEqual(args.pipeline, "import")




    @patch("pdm_utils.run.run_import")
    @patch("argparse.ArgumentParser.parse_args")
    def test_main_1(self, parse_args_mock, run_pipeline_mock):
        """Verify that main run correctly with:
        valid pipeline and valid default options."""
        parse_args_mock.return_value = self.args1
        run.main()
        with self.subTest():
            self.assertTrue(parse_args_mock.called)
        with self.subTest():
            self.assertTrue(run_pipeline_mock.called)

    @patch("pdm_utils.run.run_import")
    @patch("pdm_utils.run.get_sql_handle")
    @patch("argparse.ArgumentParser.parse_args")
    def test_main_2(self, parse_args_mock, get_sql_mock, run_pipeline_mock):
        """Verify that main run correctly with:
        valid pipeline, database provided, and sql_handle provided."""
        parse_args_mock.return_value = self.args2
        get_sql_mock.return_value = "temp"
        run.main()
        with self.subTest():
            self.assertTrue(parse_args_mock.called)
        with self.subTest():
            self.assertTrue(get_sql_mock.called)
        with self.subTest():
            self.assertTrue(run_pipeline_mock.called)

    @patch("pdm_utils.run.run_import")
    @patch("sys.exit")
    @patch("pdm_utils.run.get_sql_handle")
    @patch("argparse.ArgumentParser.parse_args")
    def test_main_3(self, parse_args_mock, get_sql_mock,
                     sys_exit_mock, run_pipeline_mock):
        """Verify that main run correctly with:
        valid pipeline, database provided, and no sql_handle provided."""
        parse_args_mock.return_value = self.args2
        get_sql_mock.return_value = None
        run.main()
        with self.subTest():
            self.assertTrue(parse_args_mock.called)
        with self.subTest():
            self.assertTrue(get_sql_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)
        with self.subTest():
            self.assertTrue(run_pipeline_mock.called)

    @patch("pdm_utils.run.run_import")
    @patch("argparse.ArgumentParser.parse_args")
    def test_main_4(self, parse_args_mock, run_pipeline_mock):
        """Verify that main runs correctly with:
        valid pipeline and valid input_folder provided."""
        parse_args_mock.return_value = self.args3
        run.main()
        with self.subTest():
            self.assertTrue(parse_args_mock.called)
        with self.subTest():
            self.assertTrue(run_pipeline_mock.called)

    @patch("pdm_utils.run.run_import")
    @patch("sys.exit")
    @patch("argparse.ArgumentParser.parse_args")
    def test_main_5(self, parse_args_mock, sys_exit_mock, run_pipeline_mock):
        """Verify that main runs correctly with:
        valid pipeline and invalid input_folder provided."""
        parse_args_mock.return_value = self.args4
        run.main()
        with self.subTest():
            self.assertTrue(parse_args_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)
        with self.subTest():
            self.assertTrue(run_pipeline_mock.called)

    @patch("pdm_utils.run.run_import")
    @patch("argparse.ArgumentParser.parse_args")
    def test_main_6(self, parse_args_mock, run_pipeline_mock):
        """Verify that main runs correctly with:
        valid pipeline and valid import_table provided."""
        parse_args_mock.return_value = self.args5
        run.main()
        with self.subTest():
            self.assertTrue(parse_args_mock.called)
        with self.subTest():
            self.assertTrue(run_pipeline_mock.called)

    @patch("pdm_utils.run.run_import")
    @patch("sys.exit")
    @patch("argparse.ArgumentParser.parse_args")
    def test_main_7(self, parse_args_mock, sys_exit_mock, run_pipeline_mock):
        """Verify that main runs correctly with:
        valid pipeline and invalid import_table provided."""
        parse_args_mock.return_value = self.args6
        run.main()
        with self.subTest():
            self.assertTrue(parse_args_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)
        with self.subTest():
            self.assertTrue(run_pipeline_mock.called)


    # In progress
    # def test_get_sql_handle_1(self):
    #     """."""
    #     args = self.parser.parse_args(["--pipeline", "import"])
    #     mch_mock.return_value
    #     parse_args_mock.return_value = args
    #     run.main()


    def tearDown(self):
        shutil.rmtree(self.test_directory1)




if __name__ == '__main__':
    unittest.main()
