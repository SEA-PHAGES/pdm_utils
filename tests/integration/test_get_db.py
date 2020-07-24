"""Integration tests for the get_db pipeline."""

from pathlib import Path
import sys
import shutil
import unittest
from unittest.mock import patch

from pdm_utils import run
from pdm_utils.constants import constants
from pdm_utils.pipelines import get_db

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

PIPELINE = "get_db"
DB = test_db_utils.DB
USER = test_db_utils.USER
PWD = test_db_utils.PWD

# TODO: There are two test files staged on the server for testing purposes.
# pdm_test_db.sql and pdm_test_db.version. If get_db pipeline is revamped
# to utilize urllib3, these two test files could be removed, and the tests
# would only need to confirm the download request status attribute = 200.
DB2 = "pdm_test_db"

# Create the main test directory in which all files will be
# created and managed. Gets created once for all tests.
test_root_dir = Path("/tmp", "pdm_utils_tests_get_db")
if test_root_dir.exists() == True:
    shutil.rmtree(test_root_dir)
test_root_dir.mkdir()

# Since these tests involve creating new databases, be sure to
# remove the existing test database if it is already present.
exists = test_db_utils.check_if_exists()
if exists:
    test_db_utils.remove_db()

# Within main test folder, this new folder will be created/removed
# for each test. Within the output_folder, get_db will dynamically create
# a new folder, but only if files are downloaded.
output_path = Path(test_root_dir, "output")
results_path = Path(output_path, get_db.RESULTS_FOLDER)

def get_unparsed_args(db=None, option=None, download=False, output_folder=None,
                      version=False, url="", interactive=False):
    """Returns list of command line arguments to convert database."""

    # to make sure that file and new always have a db value - like a default
    # db is an optional argument when option == "server"


    if option == "file":
        unparsed_args = ["run.py", PIPELINE, option, DB]
        unparsed_args.extend([str(test_db_utils.TEST_DB_FILEPATH)])
    elif option == "server":
        if db is None:
        	unparsed_args = ["run.py", PIPELINE, option]
        else:
        	unparsed_args = ["run.py", PIPELINE, option, "-db", db]

        if download:
            unparsed_args.extend(["-d"])
        if version:
            unparsed_args.extend(["-v"])
        if url != "":
            unparsed_args.extend(["-u", url])
        if output_folder is not None:
            unparsed_args.extend(["-o", str(output_folder)])
        if interactive:
            unparsed_args.extend(["-i"])
    else:
        unparsed_args = ["run.py", PIPELINE, option ,DB]
        pass
    return unparsed_args


class TestGetDb(unittest.TestCase):

    def setUp(self):
        output_path.mkdir()

    def tearDown(self):
        shutil.rmtree(output_path)
        exists = test_db_utils.check_if_exists()
        if exists:
            test_db_utils.remove_db()

    # Can't mock call to AlchemyHandler since that is called twice.
    @patch("pdm_utils.classes.alchemyhandler.getpass")
    def test_main_1(self, getpass_mock):
        """Verify database is installed from file."""
        getpass_mock.side_effect = [USER, PWD]
        unparsed_args = get_unparsed_args(option="file")
        run.main(unparsed_args)
        # Query for version data. This verifies that the databases exists
        # and that it contains a pdm_utils schema with data.
        version_data = test_db_utils.get_data(test_db_utils.version_table_query)
        self.assertEqual(len(version_data), 1)

    @patch("pdm_utils.classes.alchemyhandler.getpass")
    def test_main_2(self, getpass_mock):
        """Verify new database is created."""
        # The convert pipeline gets called, which will ask for MySQL
        # username and password. So provide these credentials twice.
        getpass_mock.side_effect = [USER, PWD, USER, PWD]
        unparsed_args = get_unparsed_args(option="new")
        run.main(unparsed_args)

        # Query for version data. This verifies that the databases exists
        # and that it contains a pdm_utils schema with data.
        version_data = test_db_utils.get_data(test_db_utils.version_table_query)
        self.assertEqual(len(version_data), 1)


    @patch("pdm_utils.classes.alchemyhandler.getpass")
    def test_main_3(self, getpass_mock):
        """Verify database sql and version files are downloaded from server."""
        getpass_mock.side_effect = [USER, PWD]
        # If the entire Actinobacteriophage database is being downloaded for testing,
        # be sure to only download the SQL file and do NOT install it,
        # else it will overwrite the existing Actinobacteriophage database.
        # Since the pdm_anon user is calling this pipeline, and since
        # this user should not have MySQL privileges to do anything other
        # than select data from Actinobacteriophage, this shouldn't be a problem.
        unparsed_args = get_unparsed_args(db=DB2, option="server",
                                          download=True, version=True,
                                          output_folder=output_path)
        run.main(unparsed_args)
        file1 = Path(results_path, f"{DB2}.sql")
        file2 = Path(results_path, f"{DB2}.version")
        with self.subTest():
            self.assertTrue(file1.exists())
        with self.subTest():
            self.assertTrue(file2.exists())


    @patch("pdm_utils.classes.alchemyhandler.getpass")
    def test_main_4(self, getpass_mock):
        """Verify database sql file, but not version file, is downloaded
        from server using url option."""
        getpass_mock.side_effect = [USER, PWD]
        # If the entire Actinobacteriophage database is being downloaded for testing,
        # be sure to only download the SQL file and do NOT install it,
        # else it will overwrite the existing Actinobacteriophage database.
        # Since the pdm_anon user is calling this pipeline, and since
        # this user should not have MySQL privileges to do anything other
        # than select data from Actinobacteriophage, this shouldn't be a problem.
        unparsed_args = get_unparsed_args(db=DB2, option="server",
                                          download=True,
                                          url=constants.DB_WEBSITE,
                                          output_folder=output_path)
        run.main(unparsed_args)
        file1 = Path(results_path, f"{DB2}.sql")
        file2 = Path(results_path, f"{DB2}.version")
        with self.subTest():
            self.assertTrue(file1.exists())
        with self.subTest():
            self.assertFalse(file2.exists())


    @patch("pdm_utils.classes.alchemyhandler.getpass")
    def test_main_5(self, getpass_mock):
        """Verify database is installed from file and overwrites
        existing database."""
        getpass_mock.side_effect = [USER, PWD]
        # First install a database with data. Then delete version table.
        test_db_utils.create_filled_test_db()
        test_db_utils.execute("DROP TABLE version")
        unparsed_args = get_unparsed_args(option="file")

        run.main(unparsed_args)
        # Now query for version data. This verifies that it replaced
        # the first database.
        version_data = test_db_utils.get_data(test_db_utils.version_table_query)
        self.assertEqual(len(version_data), 1)

    @patch("pdm_utils.pipelines.get_db.input")
    @patch("pdm_utils.classes.alchemyhandler.getpass")
    def test_main_6(self, getpass_mock, input_mock):
        """Verify that interactive mode is operational"""
        getpass_mock.side_effect = [USER, PWD]
        input_mock.return_value = 1
        # Need to ensure that there is a good response from the specified request (200)
        # option has to be server for interactive mode to work
        unparsed_args = get_unparsed_args(option="server")

        run.main(unparsed_args)

        # Now, ensure that the request returns 200
        # From get_db function request_url()

        request_data = get_db.request_url().status

        self.assertEqual(request_data, 200)



if __name__ == '__main__':
    unittest.main()
