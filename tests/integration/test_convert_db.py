"""Integration tests for the convert pipeline."""

from pathlib import Path
import sys
import shutil
import subprocess
import unittest
from unittest.mock import patch

from pdm_utils import run
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.pipelines import convert_db

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

pipeline = "convert"
USER = test_db_utils.USER
PWD = test_db_utils.PWD
DB = test_db_utils.DB
DB2 = "pdm_test_2"
query = test_db_utils.version_table_query

# Create the main test directory in which all files will be
# created and managed. Gets created once for all tests.
test_root_dir = Path("/tmp", "pdm_utils_tests_convert")
if test_root_dir.exists() == True:
    shutil.rmtree(test_root_dir)
test_root_dir.mkdir()

# Within main test folder, create new folder that will get created/removed
# for each test.
results_path = Path(test_root_dir, "output")

def get_unparsed_convert_args(schema_version=None, new_db=None):
    """Returns list of command line arguments to convert database."""
    unparsed_args = ["run.py", pipeline, DB]
    if schema_version is not None:
        unparsed_args.extend(["-s", schema_version])
    if new_db is not None:
        unparsed_args.extend(["-n", new_db])
    return unparsed_args


def diff_files(filepath1, filepath2, diff_filepath):
    """Compare two files using diff."""
    # diff filepath1 filepath2 > diff_filepath
    command_string = f"diff {filepath1} {filepath2}"
    command_list = command_string.split(" ")
    with diff_filepath.open("w") as handle:
        subprocess.check_call(command_list, stdout=handle)


class TestConvert(unittest.TestCase):

    def setUp(self):
        test_db_utils.create_filled_test_db()
        self.version_data = test_db_utils.get_data(query)
        self.current = self.version_data[0]["SchemaVersion"]
        self.alchemist = AlchemyHandler(database=DB, username=USER, password=PWD)
        self.alchemist.build_engine()
        results_path.mkdir()

    def tearDown(self):
        shutil.rmtree(results_path)
        # Remove 'pdm_test_db'
        test_db_utils.remove_db()

        # Remove 'pdm_test_2'
        exists = test_db_utils.check_if_exists(db=DB2)
        if exists:
            test_db_utils.remove_db(db=DB2)

    @patch("pdm_utils.pipelines.convert_db.AlchemyHandler")
    def test_main_1(self, alchemy_mock):
        """Verify step-by-step database conversion round trip is successful and
        generates identical empty schema files."""
        alchemy_mock.return_value = self.alchemist

        # Downgrade one schema version at a time and generate empty schema file.
        # If current = 4, down = [3, 2, 1, 0]
        down_versions = list(reversed(range(0, self.current)))
        for step in down_versions:
            unparsed_args = get_unparsed_convert_args(str(step))
            schema_filename = f"down_{unparsed_args[-1]}.sql"
            schema_filepath = Path(results_path, schema_filename)
            run.main(unparsed_args)
            test_db_utils.create_schema_file(schema_filepath)

        # Now upgrade one schema version at a time and generate empty schema file.
        # If current = 4, up = [1, 2, 3, 4]
        up_versions = list(range(1, self.current + 1))
        for step in up_versions:
            unparsed_args = get_unparsed_convert_args(str(step))
            schema_filename = f"up_{unparsed_args[-1]}.sql"
            schema_filepath = Path(results_path, schema_filename)
            run.main(unparsed_args)
            test_db_utils.create_schema_file(schema_filepath)

        # Compare the empty schema file representing the same schema version,
        # generated from a downgrade (e.g. 7 to 6) and an upgrade (e.g. 5 to 6).
        # If current = 4, compare = [1, 2, 3]
        compare_versions = list(range(1, self.current))
        diff_file_list = []
        for step in compare_versions:
            filepath1 = Path(results_path, f"up_{step}.sql")
            filepath2 = Path(results_path, f"down_{step}.sql")
            diff_filepath = Path(results_path, f"diff_{step}.txt")
            diff_file_list.append(diff_filepath)
            diff_files(filepath1, filepath2, diff_filepath)

        # Get the number of lines in each diff file. If there is no
        # difference, there should be 0 lines.
        results = []
        for diff_filepath in diff_file_list:
            with open(diff_filepath, "r") as handle:
                lines = handle.readlines()
            results.append(len(lines))

        # Test the results.
        for result in results:
            with self.subTest():
                self.assertEqual(result, 0)


    @patch("pdm_utils.pipelines.convert_db.AlchemyHandler")
    def test_main_2(self, alchemy_mock):
        """Verify non-step-by-step database conversion round trip is successful
        and generates identical empty schema files."""
        alchemy_mock.return_value = self.alchemist

        # Get current schema file.
        schema_filepath1 = Path(results_path, "before.sql")
        test_db_utils.create_schema_file(schema_filepath1)

        # Downgrade to schema version 0.
        step = 0
        unparsed_args = get_unparsed_convert_args(str(step))
        run.main(unparsed_args)

        # Now upgrade to the original schema version.
        unparsed_args = get_unparsed_convert_args(str(self.current))
        run.main(unparsed_args)

        # Get data and schema file after round trip.
        schema_filepath2 = Path(results_path, "after.sql")
        test_db_utils.create_schema_file(schema_filepath2)
        new_version_data = test_db_utils.get_data(query)
        new_schema_version = new_version_data[0]["SchemaVersion"]

        # Compare the two empty schema files representing the same schema version,
        # generated from a downgrade to 0 and upgrade.
        diff_filepath = Path(results_path, "diff.txt")
        diff_files(schema_filepath1, schema_filepath2, diff_filepath)

        # Get the number of lines in each diff file. If there is no
        # difference, there should be 0 lines.
        with open(diff_filepath, "r") as handle:
            lines = handle.readlines()
        result = len(lines)

        # Test the results.
        self.assertEqual(result, 0)


    @patch("pdm_utils.pipelines.convert_db.AlchemyHandler")
    def test_main_3(self, alchemy_mock):
        """Verify database with new name is created."""
        alchemy_mock.return_value = self.alchemist

        # Downgrade one step.
        step = self.current - 1
        unparsed_args = get_unparsed_convert_args(str(step), new_db=DB2)
        run.main(unparsed_args)

        # Get data from both databases.
        data1 = test_db_utils.get_data(query)
        v1 = data1[0]["SchemaVersion"]
        data2 = test_db_utils.get_data(query, db=DB2)
        v2 = data2[0]["SchemaVersion"]

        # Test
        with self.subTest():
            self.assertEqual(v1, self.current)
        with self.subTest():
            self.assertEqual(v2, step)


if __name__ == '__main__':
    unittest.main()
