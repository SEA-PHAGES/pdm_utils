"""Integration tests for the compare pipeline."""

from pathlib import Path
import shutil
import sys
import unittest
from unittest.mock import patch

from pdm_utils import run
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.pipelines import compare_db

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

# Create the main test directory in which all files will be
# created and managed.
test_root_dir = Path("/tmp", "pdm_utils_tests_compare")
if test_root_dir.exists() == True:
    shutil.rmtree(test_root_dir)
test_root_dir.mkdir()

# New folder that will get created/removed for each test.
test_folder = Path(test_root_dir, "output")


pipeline = "compare"
USER = test_db_utils.USER
PWD = test_db_utils.PWD
DB = test_db_utils.DB

def create_update(table, field, value, phage_id=None):
    """Creates a MySQL UPDATE statement."""
    statement = f"UPDATE {table} SET {field} = '{value}'"
    if phage_id is not None:
        statement = statement + f" WHERE PhageID = '{phage_id}'"
    return statement

def get_unparsed_args():
    """Returns list of command line arguments to compare databases."""
    unparsed_args = ["run.py", pipeline, DB,
                      "-o", str(test_folder),
                      "-p", #this will retrieve all data from PhagesDB
                      "-g", "-s", "-ia", "-d", "-f", "-u"
                    ]
    return unparsed_args

def count_files(path_to_folder):
    """Count number of files in a folder."""
    count = 0
    for item in path_to_folder.iterdir():
        count += 1
    return count


TRIXIE_ACC = "JN408461"
ALICE_ACC = "JF704092"
L5_ACC = "Z18946"
D29_ACC = "AF022214"

FASTA_FILE = ">Trixie\nATCG"

def get_pdb_dict():
    """Mock PhagesDB data dictionary."""
    dict ={
        "phage_name": "Trixie",
        "isolation_host": {"genus": "Mycobacterium"},
        "genbank_accession": "ABC123",
        "pcluster": {"cluster": "A"},
        "psubcluster": {"subcluster": "A2"},
        "fasta_file": "invalid.fasta"
    }
    return dict

def get_pdb_json_data():
    """Mock PhagesDB json data."""
    dict = {"results": []}
    return dict



class TestCompare(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

    @classmethod
    def tearDownClass(self):
        # Remove 'pdm_test_db'
        test_db_utils.remove_db()

    def setUp(self):
        test_folder.mkdir()

        # Standardize values in certain fields to define the data
        stmts = []
        stmts.append(create_update("phage", "Status", "draft"))
        stmts.append(create_update("phage", "Accession", ""))
        stmts.append(create_update("phage", "AnnotationAuthor", "0"))

        stmts.append(create_update("phage", "Accession", TRIXIE_ACC, "Trixie"))
        stmts.append(create_update("phage", "Accession", ALICE_ACC, "Alice"))
        stmts.append(create_update("phage", "Accession", L5_ACC, "L5"))
        stmts.append(create_update("phage", "Accession", TRIXIE_ACC, "D29"))

        stmts.append(create_update("phage", "Status", "final", "Trixie"))
        stmts.append(create_update("phage", "Status", "final", "Alice"))
        stmts.append(create_update("phage", "Status", "final", "L5"))
        stmts.append(create_update("phage", "Status", "final", "D29"))

        stmts.append(create_update("phage", "AnnotationAuthor", "1", "Trixie"))
        stmts.append(create_update("phage", "AnnotationAuthor", "1", "Alice"))
        stmts.append(create_update("phage", "AnnotationAuthor", "1", "L5"))
        stmts.append(create_update("phage", "AnnotationAuthor", "1", "D29"))

        for stmt in stmts:
            test_db_utils.execute(stmt)

        self.unparsed_args = get_unparsed_args()

        self.alchemist = AlchemyHandler(database=DB, username=USER, password=PWD)
        self.alchemist.build_engine()

        self.pdb_data1 = get_pdb_dict()
        self.pdb_data2 = get_pdb_dict()
        self.pdb_data3 = get_pdb_dict()
        self.pdb_data1["phage_name"] = "Trixie"
        self.pdb_data2["phage_name"] = "L5"
        self.pdb_data3["phage_name"] = "unmatched"

        json_results = [self.pdb_data1, self.pdb_data2, self.pdb_data3]
        self.pdb_json_data = get_pdb_json_data()
        self.pdb_json_data["results"] = json_results
        self.pdb_json_results = json_results

    def tearDown(self):
        shutil.rmtree(test_folder)


    def test_queries_1(self):
        """Verify hard-coded SQL queries are structured correctly."""

        version_data = test_db_utils.get_data(compare_db.VERSION_QUERY)
        phage_data = test_db_utils.get_data(compare_db.PHAGE_QUERY)
        gene_data = test_db_utils.get_data(compare_db.GENE_QUERY)

        ref_phage_data = test_db_utils.get_data(test_db_utils.phage_table_query)
        ref_gene_data = test_db_utils.get_data(test_db_utils.gene_table_query)

        with self.subTest():
            self.assertEqual(len(version_data), 1)
        with self.subTest():
            self.assertEqual(len(phage_data), len(ref_phage_data))
        with self.subTest():
            self.assertEqual(len(gene_data), len(ref_gene_data))


    # Calls to PhagesDB need to be mocked, since the pipeline downloads
    # and parses all sequenced genome data.
    @patch("pdm_utils.functions.phagesdb.retrieve_url_data")
    @patch("pdm_utils.functions.phagesdb.get_phagesdb_data")
    @patch("pdm_utils.pipelines.compare_db.AlchemyHandler")
    def test_main_1(self, alchemy_mock, gpd_mock, rud_mock):
        """Verify compare runs successfully with:
        MySQL, PhagesDB, and GenBank records saved,
        a duplicate MySQL accession (for D29 and Trixie),
        an invalid accession (for L5),
        a duplicate phage.Name (for Constance and Et2Brutus),
        a PhagesDB name unmatched to MySQL (for 'unmatched')."""
        alchemy_mock.return_value = self.alchemist
        gpd_mock.return_value = self.pdb_json_results
        rud_mock.return_value = FASTA_FILE

        # Make modifications to cause errors.
        stmts = []
        stmts.append(create_update("phage", "Accession", L5_ACC + "1", "L5"))
        stmts.append(create_update("phage", "Accession", TRIXIE_ACC, "D29"))

        stmts.append(create_update("phage", "Name", "Dupe", "Constance"))
        stmts.append(create_update("phage", "Name", "Dupe", "Et2Brutus"))

        stmts.append(create_update("phage", "PhageID", "Dupe", "Constance"))
        stmts.append(create_update("phage", "PhageID", "Dupe_Draft", "Et2Brutus"))
        for stmt in stmts:
            test_db_utils.execute(stmt)


        run.main(self.unparsed_args)
        count = count_files(test_folder)
        # input("check")
        with self.subTest():
            self.assertTrue(count > 0)
        with self.subTest():
            gpd_mock.assert_called()
        with self.subTest():
            rud_mock.assert_called()

    # Calls to PhagesDB need to be mocked, since the pipeline downloads
    # and parses all sequenced genome data.
    @patch("pdm_utils.pipelines.compare_db.get_pdb_data")
    @patch("pdm_utils.pipelines.compare_db.AlchemyHandler")
    def test_main_2(self, alchemy_mock, gpd_mock):
        """Verify duplicate PhagesDB names are identified."""
        # Clear accessions so that GenBank is not queried. No need for that.
        stmt = create_update("phage", "Accession", "")
        test_db_utils.execute(stmt)

        alchemy_mock.return_value = self.alchemist
        gpd_mock.return_value = ({}, {"L5"}, {"D29"})
        run.main(self.unparsed_args)
        count = count_files(test_folder)
        # input("check")
        with self.subTest():
            self.assertTrue(count > 0)
        with self.subTest():
            gpd_mock.assert_called()

    # Calls to PhagesDB need to be mocked, since the pipeline downloads
    # and parses all sequenced genome data.
    @patch("pdm_utils.pipelines.compare_db.get_pdb_data")
    @patch("pdm_utils.pipelines.compare_db.filter_mysql_genomes")
    @patch("pdm_utils.pipelines.compare_db.AlchemyHandler")
    def test_main_3(self, alchemy_mock, fmg_mock, gpd_mock):
        """Verify duplicate MySQL names are identified."""
        # Clear accessions so that GenBank is not queried. No need for that.
        stmt = create_update("phage", "Accession", "")
        test_db_utils.execute(stmt)

        alchemy_mock.return_value = self.alchemist
        fmg_mock.return_value = ({}, set(), {"L5"}, set(), set())
        gpd_mock.return_value = ({}, set(), set())

        run.main(self.unparsed_args)
        count = count_files(test_folder)
        # input("check")
        with self.subTest():
            self.assertTrue(count > 0)
        with self.subTest():
            fmg_mock.assert_called()
        with self.subTest():
            gpd_mock.assert_called()





if __name__ == '__main__':
    unittest.main()
