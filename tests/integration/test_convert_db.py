"""Integration tests for the convert pipeline."""

from pathlib import Path
import sys
import shutil
import subprocess
import unittest
from unittest.mock import patch

import sqlalchemy

from pdm_utils import run
from pdm_utils.pipelines import convert_db

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

#sqlalchemy setup
engine_string = test_db_utils.create_engine_string()

pipeline = "convert"
user = test_db_utils.USER
pwd = test_db_utils.PWD
db = test_db_utils.DB

# Create the main test directory in which all files will be
# created and managed. Gets created once for all tests.
test_root_dir = Path("/tmp", "pdm_utils_tests_convert")
if test_root_dir.exists() == True:
    shutil.rmtree(test_root_dir)
test_root_dir.mkdir()

# Within main test folder, create new folder that will get created/removed
# for each test.
results_path = Path(test_root_dir, "output")

def get_unparsed_convert_args(schema_version=None):
    """Returns list of command line arguments to convert database."""
    unparsed_args = ["run.py", pipeline, db]
    if schema_version is not None:
        unparsed_args.extend(["-s", schema_version])
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
        self.version_data = test_db_utils.get_data(test_db_utils.version_table_query)
        self.current = self.version_data[0]["SchemaVersion"]
        self.engine = sqlalchemy.create_engine(engine_string, echo=False)
        results_path.mkdir()

    def tearDown(self):
        shutil.rmtree(results_path)
        # Remove 'pdm_test_db'
        test_db_utils.remove_db()

    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    def test_main_1(self, ctd_mock):
        """Verify step-by-step database conversion round trip is successful and
        generates identical empty schema files."""
        ctd_mock.return_value = self.engine

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


    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    def test_main_2(self, ctd_mock):
        """Verify non-step-by-step database conversion round trip is successful
        and generates identical empty schema files."""
        ctd_mock.return_value = self.engine

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
        new_version_data = test_db_utils.get_data(test_db_utils.version_table_query)
        new_schem_version = new_version_data[0]["SchemaVersion"]

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




if __name__ == '__main__':
    unittest.main()
