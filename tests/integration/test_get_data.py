"""Integration tests for the get_data pipeline."""

from pathlib import Path
import shutil
import sqlalchemy
import sys
import unittest
from unittest.mock import patch
from pdm_utils import run
from pdm_utils.pipelines import get_data
from pdm_utils.classes import genomepair
from pdm_utils.classes import genome

# Import helper functions to build mock database
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils

#sqlalchemy setup
engine_string = test_db_utils.create_engine_string()

pipeline = "get_data"
db = test_db_utils.DB

# Create the main test directory in which all files will be
# created and managed.
test_root_dir = Path("/tmp", "pdm_utils_tests_get_data")
if test_root_dir.exists() == True:
    shutil.rmtree(test_root_dir)
test_root_dir.mkdir()

# New folder that will get created/removed for each test.
test_folder = Path(test_root_dir, "output")

# How the output folder and files are named in the get_data pipeline.
results_folder = Path(get_data.RESULTS_FOLDER)
results_path = Path(test_folder, results_folder)

pecaan_folder = Path(results_path, get_data.PECAAN_FOLDER)
genbank_folder = Path(results_path, get_data.GENBANK_FOLDER)
phagesdb_folder = Path(results_path, get_data.PHAGESDB_FOLDER)
updates_folder = Path(results_path, get_data.UPDATES_FOLDER)

phagesdb_genome_folder = Path(phagesdb_folder, get_data.GENOME_FOLDER)
pecaan_genomes_folder = Path(pecaan_folder, get_data.GENOME_FOLDER)
genbank_genomes_folder = Path(genbank_folder, get_data.GENOME_FOLDER)

phagesdb_import_table = Path(phagesdb_folder, get_data.IMPORT_TABLE)
pecaan_import_table = Path(pecaan_folder, get_data.IMPORT_TABLE)
genbank_import_table = Path(genbank_folder, get_data.IMPORT_TABLE)

genbank_results_table = Path(genbank_folder, get_data.GENBANK_RESULTS_TABLE)

update_table = Path(updates_folder, get_data.UPDATE_TABLE)

TRIXIE_ACC = "JN408461"

def create_update(table, field, value, phage_id=None):
    """Creates a MySQL UPDATE statement."""
    num_field_set = {"RetrieveRecord", "AnnotationAuthor"}
    if field not in num_field_set:
        value = f"'{value}'"
    statement = f"UPDATE {table} SET {field} = {value}"
    if phage_id is not None:
        statement = statement + f" WHERE PhageID = '{phage_id}'"
    return statement

def get_unparsed_args(draft=False, final=False, genbank=False, update=False,
                      force_download=False, genbank_results=True):
    """Returns list of command line arguments to get data."""
    unparsed_args = ["run.py", pipeline, db, str(test_folder)]
    if draft:
        unparsed_args.append("-d")
    if final:
        unparsed_args.append("-f")
    if genbank:
        unparsed_args.append("-g")
    if update:
        unparsed_args.append("-u")
    if force_download:
        unparsed_args.append("-fd")
    if genbank_results:
        unparsed_args.append("-gr")
    return unparsed_args

def count_files(path_to_folder):
    """Count number of files in a folder."""
    count = 0
    for item in path_to_folder.iterdir():
        count += 1
    return count


def count_rows(filepath, header=True):
    """Count number of rows in a file."""
    with open(filepath, "r") as handle:
        lines = handle.readlines()
    count = len(lines)
    if header:
        count -= 1
    return count


def create_matched_genomes():
    """Create list of GenomePair objects."""

    gnm1 = genome.Genome()
    gnm1.id = "Trixie"
    gnm1.annotation_status = "draft"

    gnm2 = genome.Genome()
    gnm2.id = "Trixie"

    gnm_pair1 = genomepair.GenomePair()
    gnm_pair1.genome1 = gnm1
    gnm_pair1.genome2 = gnm2

    gnm3 = genome.Genome()
    gnm3.id = "Alice"
    gnm3.annotation_status = "final"

    gnm4 = genome.Genome()
    gnm4.id = "Alice"

    gnm_pair2 = genomepair.GenomePair()
    gnm_pair2.genome1 = gnm3
    gnm_pair2.genome2 = gnm4

    matched_genomes = [gnm_pair1, gnm_pair2]
    return matched_genomes




class TestGetData(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        test_db_utils.create_filled_test_db()

    @classmethod
    def tearDownClass(self):
        test_db_utils.remove_db()

    def setUp(self):
        self.engine = sqlalchemy.create_engine(engine_string, echo=False)
        test_folder.mkdir()

        # Standardize values in certain fields to define the data
        stmt1 = create_update("phage", "Status", "unknown")
        test_db_utils.execute(stmt1)
        stmt2 = create_update("phage", "HostGenus", "unknown")
        test_db_utils.execute(stmt2)
        stmt3 = create_update("phage", "Accession", "")
        test_db_utils.execute(stmt3)
        stmt4 = create_update("phage", "DateLastModified", "1900-01-01")
        test_db_utils.execute(stmt4)
        stmt5 = create_update("phage", "RetrieveRecord", "0")
        test_db_utils.execute(stmt5)

    def tearDown(self):
        shutil.rmtree(test_folder)




    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    def test_main_1(self, ctd_mock):
        """Verify update data and final data are retrieved."""
        # Testing the update flag and final flag have been merged so that
        # PhagesDB is only queried once for all data in the genome, since
        # it is time-intensive.
        ctd_mock.return_value = self.engine
        # If final=True, any genome in database will be checked on PhagesDB
        # regardless of AnnotationAuthor
        unparsed_args = get_unparsed_args(update=True, final=True)
        run.main(unparsed_args)
        count1 = count_files(updates_folder)
        count2 = count_files(phagesdb_genome_folder)
        count3 = count_rows(update_table)
        count4 = count_rows(phagesdb_import_table)
        with self.subTest():
            self.assertEqual(count1, 1)
        with self.subTest():
            # It's not clear how stable the storage of any particular final
            # flat file is on PhagesDB. There is possibility that for any
            # genome the associated flat file will be removed. So this is
            # one area of testing that could be improved. For now, simply
            # verify that 1 or more files have been retrieved.
            self.assertTrue(count2 > 0)
        with self.subTest():
            # There should be several rows of updates.
            self.assertTrue(count3 > 0)
        with self.subTest():
            self.assertEqual(count2, count4)

    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    def test_main_2(self, ctd_mock):
        """Verify genbank data is retrieved."""
        ctd_mock.return_value = self.engine
        stmt1 = create_update("phage", "RetrieveRecord", "1", phage_id="Trixie")
        test_db_utils.execute(stmt1)
        stmt2 = create_update("phage", "Accession", TRIXIE_ACC, phage_id="Trixie")
        test_db_utils.execute(stmt2)
        unparsed_args = get_unparsed_args(genbank=True, genbank_results=True)
        run.main(unparsed_args)
        query = "SELECT COUNT(*) FROM phage"
        count1 = test_db_utils.get_data(query)[0]["COUNT(*)"]
        count2 = count_files(genbank_genomes_folder)
        count3 = count_rows(genbank_import_table)
        count4 = count_rows(genbank_results_table)
        with self.subTest():
            self.assertEqual(count2, 1)
        with self.subTest():
            self.assertEqual(count2, count3)
        with self.subTest():
            self.assertEqual(count1, count4)

    @patch("pdm_utils.pipelines.get_data.match_genomes")
    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    def test_main_3(self, ctd_mock, mg_mock):
        """Verify draft data is retrieved."""
        matched_genomes = create_matched_genomes()
        ctd_mock.return_value = self.engine
        mg_mock.return_value = (matched_genomes, {"EagleEye"})
        unparsed_args = get_unparsed_args(draft=True)
        run.main(unparsed_args)
        count1 = count_files(pecaan_genomes_folder)
        count2 = count_rows(pecaan_import_table)
        with self.subTest():
            self.assertEqual(count1, 1)
        with self.subTest():
            self.assertEqual(count1, count2)

    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    def test_main_4(self, ctd_mock):
        """Verify final data with very recent date are retrieved
        with force_download."""
        ctd_mock.return_value = self.engine
        stmt = create_update("phage", "DateLastModified", "2200-01-01")
        test_db_utils.execute(stmt)
        unparsed_args = get_unparsed_args(final=True, force_download=True)
        run.main(unparsed_args)
        count = count_files(phagesdb_genome_folder)
        with self.subTest():
            # It's not clear how stable the storage of any particular final
            # flat file is on PhagesDB. There is possibility that for any
            # genome the associated flat file will be removed. So this is
            # one area of testing that could be improved. For now, simply
            # verify that 1 or more files have been retrieved.
            self.assertTrue(count > 0)

    @patch("pdm_utils.pipelines.get_data.match_genomes")
    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    def test_main_5(self, ctd_mock, mg_mock):
        """Verify draft data already in database is retrieved
        with force_download."""
        # Create a list of 2 matched genomes, only one of which has
        # status = draft.
        matched_genomes = create_matched_genomes()
        ctd_mock.return_value = self.engine
        mg_mock.return_value = (matched_genomes, {"EagleEye"})
        unparsed_args = get_unparsed_args(draft=True, force_download=True)
        run.main(unparsed_args)
        count = count_files(pecaan_genomes_folder)
        self.assertEqual(count, 2)




if __name__ == '__main__':
    unittest.main()
