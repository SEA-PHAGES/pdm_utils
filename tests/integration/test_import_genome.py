"""Integration tests for the main import pipeline."""

import csv
import getpass
from pathlib import Path
import pymysql
import shutil
import subprocess
import time
import unittest
from unittest.mock import patch
import sqlalchemy
from pdm_utils.classes import bundle, genome, genomepair, ticket, eval, cds, source
from pdm_utils.constants import constants
from pdm_utils.functions import basic, mysqldb
from pdm_utils.pipelines import import_genome


# Create the main test directory in which all files will be
# created and managed.
test_root_dir = Path("/tmp", "pdm_utils_tests_import1")
if test_root_dir.exists() == True:
    shutil.rmtree(test_root_dir)
test_root_dir.mkdir()

def count_status(item, *args):
    count = 0
    status_set = set()
    for arg in args:
        status_set.add(arg)
    for evl in item.evaluations:
        if evl.status in status_set:
            count += 1
    return count

def count_status_from_dict(dict, *args):
    """Iterate through a dictionary produced from bundle.get_evaluations()."""
    status_set = set()
    for arg in args:
        status_set.add(arg)
    count = 0
    for key in dict:
        sub_dict = dict[key]
        for sub_key in sub_dict:
            evl_list = sub_dict[sub_key]
            for evl in evl_list:
                if evl.status in status_set:
                    count += 1
    return count

# The following integration tests user the 'pdm_anon' MySQL user.
# It is expected that this user has all privileges for 'test_db' database.
user = "pdm_anon"
pwd = "pdm_anon"
db = "test_db"
db2 = "Actinobacteriophage"
#sqlalchemy setup
engine_string1 = f"mysql+pymysql://{user}:{pwd}@localhost/{db}"
engine_string2 = f"mysql+pymysql://{user}:{pwd}@localhost/{db2}"

unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
test_file_dir = Path(test_dir, "test_files")
schema_version = constants.CODE_SCHEMA_VERSION
schema_file = f"test_schema_{schema_version}.sql"
schema_filepath = Path(test_file_dir, schema_file)
version_table_data = {"Version":1, "SchemaVersion":schema_version}


def create_new_db(schema_filepath, db, user, pwd):
    """Creates a new, empty database."""
    connection = pymysql.connect(host = "localhost",
                                 user = user,
                                 password = pwd,
                                 cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()

    # First, test if a test database already exists within mysql.
    # If there is, delete it so that a fresh test database is installed.
    sql = ("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA "
          f"WHERE SCHEMA_NAME = '{db}'")
    cur.execute(sql)
    result = cur.fetchall()
    if len(result) != 0:
        cur.execute(f"DROP DATABASE {db}")
        connection.commit()

    # Next, create the database within mysql.
    cur.execute(f"CREATE DATABASE {db}")
    connection.commit()
    connection.close()

    # Now import the empty schema from file.
    # Seems like pymysql has trouble with this step, so use subprocess.
    handle = open(schema_filepath, "r")
    command_string = f"mysql -u {user} -p{pwd} {db}"
    command_list = command_string.split(" ")
    proc = subprocess.check_call(command_list, stdin = handle)
    handle.close()


def remove_db(db, user, pwd):
    """Remove the MySQL database created for the test."""
    connection = pymysql.connect(host="localhost",
                                 user=user,
                                 password=pwd,
                                 cursorclass=pymysql.cursors.DictCursor)
    cur = connection.cursor()
    cur.execute(f"DROP DATABASE {db}")
    connection.commit()
    connection.close()

def insert_data_into_version_table(db, user, pwd, data_dict):
    """Insert data into the version table."""
    connection = pymysql.connect(host = "localhost",
                                 user = user,
                                 password = pwd,
                                 database = db,
                                 cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()
    sql = (
        "INSERT INTO version "
        "(Version, SchemaVersion) "
        "VALUES ("
        f"'{data_dict['Version']}', '{data_dict['SchemaVersion']}');"
        )
    cur.execute(sql)
    connection.commit()
    connection.close()


def count_contents(path_to_folder):
    count = 0
    for item in path_to_folder.iterdir():
        count += 1
    return count


class TestImportGenomeMain1(unittest.TestCase):


    def setUp(self):
        create_new_db(schema_filepath, db, user, pwd)
        insert_data_into_version_table(db, user, pwd, version_table_data)

        self.base_dir = Path(test_root_dir, "test_import")
        self.base_dir.mkdir()

        self.genome_folder = Path(self.base_dir,"genome_folder")
        self.genome_folder.mkdir()

        self.test_flat_file1 = Path(test_file_dir, "test_flat_file_1.gb")
        self.test_flat_file2 = Path(test_file_dir, "test_flat_file_2.gb")

        self.engine = sqlalchemy.create_engine(engine_string1, echo=False)

        self.eval_flags = {
            "check_locus_tag":True,
            "check_description_field":True,
            "check_replace":True,
            "check_trna":True,
            "import_locus_tag":True,
            "check_id_typo":True,
            "check_host_typo":True,
            "check_author":True,
            "check_description":True,
            "check_gene":True
            }

        self.data_dict = {}
        self.data_dict["host_genus"] = "Arthrobacter"
        self.data_dict["cluster"] = "B"
        self.data_dict["subcluster"] = "B2"
        self.data_dict["annotation_status"] = "draft"
        self.data_dict["annotation_author"] = 1
        self.data_dict["retrieve_record"] = 1
        self.data_dict["accession"] = "ABC123"

        self.tkt1 = ticket.ImportTicket()
        self.tkt1.id = 1
        self.tkt1.type = "add"
        self.tkt1.phage_id = "L5"
        self.tkt1.eval_mode = "final"
        self.tkt1.description_field = "product"
        self.tkt1.eval_flags = self.eval_flags
        self.tkt1.data_dict = self.data_dict

        self.tkt2 = ticket.ImportTicket()


    def tearDown(self):

        # Remove all contents in the directory created for the test.
        shutil.rmtree(self.base_dir)

        # Close all open connections.
        self.engine.dispose()

        # Remove the MySQL database created for the test.
        remove_db(db, user, pwd)

    def test_prepare_bundle_1(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'add' ticket, no PhagesDB data."""
        # Omit cluster to verify that only attributes in the data_add set
        # are copied.
        self.tkt1.data_add = set(["host_genus", "subcluster",
                                  "annotation_status", "annotation_author",
                                  "retrieve_record", "accession"])

        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    host_genus_field="_organism_host_genus",
                    file_ref="flat_file",
                    ticket_ref="ticket")
        ff_gnm = bndl.genome_dict["flat_file"]
        tkt_gnm = bndl.genome_dict["ticket"]
        bndl_tkt = bndl.ticket
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 2)
        with self.subTest():
            self.assertEqual(bndl.id, 1)
        with self.subTest():
            self.assertEqual(ff_gnm.id, "L5")
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(ff_gnm.host_genus, "Arthrobacter")
        with self.subTest():
            self.assertEqual(ff_gnm.cluster, "")
        with self.subTest():
            self.assertEqual(bndl_tkt.phage_id, "L5")
        with self.subTest():
            self.assertTrue(ff_gnm._cds_descriptions_tally > 0)
        with self.subTest():
            self.assertEqual(ff_gnm._cds_descriptions_tally,
                             ff_gnm._cds_products_tally)


    def test_prepare_bundle_2(self):
        """Verify bundle is returned from a flat file with:
        no record."""
        self.tkt1.data_add = set(["host_genus", "cluster", "subcluster",
                                  "annotation_status", "annotation_author",
                                  "retrieve_record", "accession"])
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file2,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name")
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 0)
        with self.subTest():
            self.assertIsNone(bndl.ticket)


    def test_prepare_bundle_3(self):
        """Verify bundle is returned from a flat file with:
        one record, no ticket."""
        self.tkt1.data_add = set(["host_genus", "cluster", "subcluster",
                                  "annotation_status", "annotation_author",
                                  "retrieve_record", "accession"])
        tkt_dict = {"L5x":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    file_ref="flat_file")
        ff_gnm = bndl.genome_dict["flat_file"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 1)
        with self.subTest():
            self.assertIsNone(bndl.ticket)
        with self.subTest():
            self.assertEqual(ff_gnm.id, "L5")
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, -1)


    def test_prepare_bundle_4(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'add' ticket, no data_add, no PhagesDB data."""
        self.tkt1.data_add = set()
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    file_ref="flat_file")
        ff_gnm = bndl.genome_dict["flat_file"]
        bndl_tkt = bndl.ticket
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(ff_gnm.retrieve_record, -1)


    def test_prepare_bundle_5(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'add' ticket, with PhagesDB data."""
        # Use cluster and host_genus to confirm that only attributes
        # within the data_retrieve set are copied.
        self.tkt1.data_add = set(["cluster", "subcluster",
                                  "annotation_status", "annotation_author",
                                  "retrieve_record", "accession"])
        self.tkt1.data_dict["host_genus"] = "retrieve"
        self.tkt1.data_retrieve = set(["host_genus"])
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    file_ref="flat_file",
                    retrieve_ref="phagesdb")
        ff_gnm = bndl.genome_dict["flat_file"]
        pdb_gnm = bndl.genome_dict["phagesdb"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 3)
        with self.subTest():
            self.assertEqual(ff_gnm.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(ff_gnm.cluster, "B")
        with self.subTest():
            self.assertEqual(pdb_gnm.cluster, "A")


    @patch("pdm_utils.functions.phagesdb.get_genome")
    def test_prepare_bundle_6(self, gg_mock):
        """Verify bundle is returned from a flat file with:
        one record, one 'add' ticket, with PhagesDB data to be retrieved,
        but unable to retrieve PhagesDB data."""
        # Use cluster and host_genus to confirm that only attributes
        # within the data_retrieve set are copied.
        gg_mock.return_value = None
        self.tkt1.data_add = set(["cluster", "subcluster",
                                  "annotation_status", "annotation_author",
                                  "retrieve_record", "accession"])
        self.tkt1.data_dict["host_genus"] = "retrieve"
        self.tkt1.data_retrieve = set(["host_genus"])
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    file_ref="flat_file",
                    retrieve_ref="phagesdb")
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 2)
        with self.subTest():
            self.assertFalse("phagesdb" in bndl.genome_dict.keys())


    def test_prepare_bundle_7(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'replace' ticket, with MySQL data,
        and no PhagesDB data."""
        # Use host_genus and accession to confirm that only attributes
        # in the data_retain set are copied.
        l5_data = ["L5", "EFG789", "L5_Draft", "Gordonia", "ATCG",
                   4, 1, "draft", constants.EMPTY_DATE, 1, 1]
        trixie_data = ["Trixie", "XYZ456", "Trixie", "Mycobacterium", "AATTC",
                       5, 1, "final", constants.EMPTY_DATE, 1, 1]
        d29_data = ["D29", "XYZ456", "D29", "Mycobacterium", "GGCCATT",
                    7, 1, "final", constants.EMPTY_DATE, 1, 1]
        input_phage_ids_and_seqs = [l5_data, trixie_data, d29_data]
        connection = pymysql.connect(host = "localhost",
                                     user = user,
                                     password = pwd,
                                     database = db,
                                     cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for data in input_phage_ids_and_seqs:
            sql1 = (
                "INSERT INTO phage (PhageID, Accession, Name, "
                "HostGenus, Sequence, Length, GC, Status, "
                "DateLastModified, RetrieveRecord, "
                "AnnotationAuthor) VALUES ("
                f"'{data[0]}', '{data[1]}', '{data[2]}', "
                f"'{data[3]}', '{data[4]}', "
                f" {data[5]}, {data[6]}, '{data[7]}', "
                f"'{data[8]}', {data[9]}, {data[10]});"
                )
            cur.execute(sql1)
        connection.commit()
        connection.close()

        self.tkt1.type = "replace"
        self.tkt1.data_dict["host_genus"] = "retain"
        self.tkt1.data_add = set(["cluster", "subcluster",
                                     "annotation_status", "annotation_author",
                                     "retrieve_record", "accession"])
        self.tkt1.data_retain = set(["host_genus"])
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    file_ref="flat_file",
                    retain_ref="mysql",
                    engine=self.engine)
        ff_gnm = bndl.genome_dict["flat_file"]
        pmr_gnm = bndl.genome_dict["mysql"]
        ff_pmr_pair = bndl.genome_pair_dict["flat_file_mysql"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 3)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(ff_gnm.host_genus, "Gordonia")
        with self.subTest():
            self.assertEqual(ff_gnm.accession, "ABC123")
        with self.subTest():
            self.assertEqual(pmr_gnm.accession, "EFG789")


    def test_prepare_bundle_8(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'replace' ticket, no MySQL data,
        and no PhagesDB data."""

        l5_data = ["L5x", "ABC123", "L5_Draft", "Gordonia", "ATCG",
                   4, 1, "draft", constants.EMPTY_DATE, 1, 1]
        trixie_data = ["Trixie", "XYZ456", "Trixie", "Gordonia", "AATTC",
                       5, 1, "final", constants.EMPTY_DATE, 1, 1]
        d29_data = ["D29", "XYZ456", "D29", "Mycobacterium", "GGCCATT",
                    7, 1, "final", constants.EMPTY_DATE, 1, 1]
        input_phage_ids_and_seqs = [l5_data, trixie_data, d29_data]
        connection = pymysql.connect(host = "localhost",
                                     user = user,
                                     password = pwd,
                                     database = db,
                                     cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for data in input_phage_ids_and_seqs:
            sql1 = (
                "INSERT INTO phage (PhageID, Accession, Name, "
                "HostGenus, Sequence, Length, GC, Status, "
                "DateLastModified, RetrieveRecord, "
                "AnnotationAuthor) VALUES ("
                f"'{data[0]}', '{data[1]}', '{data[2]}', "
                f"'{data[3]}', '{data[4]}', "
                f" {data[5]}, {data[6]}, '{data[7]}', "
                f"'{data[8]}', {data[9]}, {data[10]});"
                )
            cur.execute(sql1)
        connection.commit()
        connection.close()

        self.tkt1.type = "replace"
        self.tkt1.data_dict["host_genus"] = "retain"
        self.tkt1.data_add = set(["cluster", "subcluster",
                                     "annotation_status", "annotation_author",
                                     "retrieve_record", "accession"])
        self.tkt1.data_retain = set(["host_genus"])
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    file_ref="flat_file",
                    engine=self.engine)
        ff_gnm = bndl.genome_dict["flat_file"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 2)
        with self.subTest():
            self.assertEqual(ff_gnm.host_genus, "")


    def test_prepare_bundle_9(self):
        """Verify bundle is returned from a flat file with:
        one record, one 'replace' ticket, with MySQL data,
        no MySQL engine, and no PhagesDB data."""

        l5_data = ["L5", "ABC123", "L5_Draft", "Gordonia", "ATCG",
                   4, 1, "draft", constants.EMPTY_DATE, 1, 1]
        trixie_data = ["Trixie", "XYZ456", "Trixie", "Mycobacterium", "AATTC",
                       5, 1, "final", constants.EMPTY_DATE, 1, 1]
        d29_data = ["D29", "XYZ456", "D29", "Mycobacterium", "GGCCATT",
                    7, 1, "final", constants.EMPTY_DATE, 1, 1]
        input_phage_ids_and_seqs = [l5_data, trixie_data, d29_data]
        connection = pymysql.connect(host = "localhost",
                                     user = user,
                                     password = pwd,
                                     database = db,
                                     cursorclass = pymysql.cursors.DictCursor)
        cur = connection.cursor()
        for data in input_phage_ids_and_seqs:
            sql1 = (
                "INSERT INTO phage (PhageID, Accession, Name, "
                "HostGenus, Sequence, Length, GC, Status, "
                "DateLastModified, RetrieveRecord, "
                "AnnotationAuthor) VALUES ("
                f"'{data[0]}', '{data[1]}', '{data[2]}', "
                f"'{data[3]}', '{data[4]}', "
                f" {data[5]}, {data[6]}, '{data[7]}', "
                f"'{data[8]}', {data[9]}, {data[10]});"
                )
            cur.execute(sql1)
        connection.commit()
        connection.close()

        self.tkt1.type = "replace"
        self.tkt1.data_dict["host_genus"] = "retain"
        self.tkt1.data_add = set(["cluster", "subcluster",
                                     "annotation_status", "annotation_author",
                                     "retrieve_record", "accession"])
        self.tkt1.data_retain = set(["host_genus"])
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    file_ref="flat_file")
        ff_gnm = bndl.genome_dict["flat_file"]
        with self.subTest():
            self.assertEqual(len(bndl.genome_dict.keys()), 2)
        with self.subTest():
            self.assertEqual(len(bndl.genome_pair_dict.keys()), 0)


    def test_prepare_bundle_10(self):
        """Verify bundle is returned with the genome id converted using
        the id dictionary."""
        id_dict = {"L5": "new_id"}
        self.tkt1.phage_id = "new_id"
        self.tkt1.data_add = set(["host_genus", "cluster", "subcluster",
                                  "annotation_status", "annotation_author",
                                  "retrieve_record", "accession"])

        tkt_dict = {"new_id":self.tkt1, "Trixie":self.tkt2}
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    host_genus_field="_organism_host_genus",
                    file_ref="flat_file",
                    ticket_ref="ticket",
                    id_conversion_dict=id_dict)
        ff_gnm = bndl.genome_dict["flat_file"]
        tkt_gnm = bndl.genome_dict["ticket"]
        bndl_tkt = bndl.ticket
        with self.subTest():
            self.assertEqual(ff_gnm.id, "new_id")
        with self.subTest():
            self.assertEqual(ff_gnm.name, "L5")
        with self.subTest():
            self.assertEqual(ff_gnm.cluster, "B")
        with self.subTest():
            self.assertEqual(ff_gnm.subcluster, "B2")
        with self.subTest():
            self.assertEqual(bndl_tkt.phage_id, "new_id")


    @patch("pdm_utils.functions.basic.choose_from_list")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_prepare_bundle_11(self, ask_mock, choose_mock):
        """Verify bundle is returned with the CDS descriptions derived
        from 'product' instead of 'function', after it is
        interactively selected using interactive = True.
        There are 0 'function' descriptions.
        There are ~30 of 'product' descriptions."""
        self.tkt1.data_add = set(["host_genus", "subcluster",
                                  "annotation_status", "annotation_author",
                                  "retrieve_record", "accession"])
        self.tkt1.description_field = "function"
        tkt_dict = {"L5":self.tkt1, "Trixie":self.tkt2}
        ask_mock.side_effect = [False]
        choose_mock.side_effect = ["product"]
        bndl = import_genome.prepare_bundle(
                    filepath=self.test_flat_file1,
                    ticket_dict=tkt_dict, id=1,
                    genome_id_field="_organism_name",
                    host_genus_field="_organism_host_genus",
                    file_ref="flat_file",
                    ticket_ref="ticket",
                    interactive=True)
        ff_gnm = bndl.genome_dict["flat_file"]
        tkt_gnm = bndl.genome_dict["ticket"]
        bndl_tkt = bndl.ticket
        with self.subTest():
            self.assertTrue(ff_gnm._cds_descriptions_tally != \
                            ff_gnm._cds_functions_tally)
        with self.subTest():
            self.assertEqual(ff_gnm._cds_descriptions_tally,
                             ff_gnm._cds_products_tally)




class TestImportGenomeMain2(unittest.TestCase):

    def setUp(self):

        self.test_filepath1 = Path(test_file_dir, "test_flat_file_1.gb")
        self.test_directory1 = Path(test_root_dir, "test_dir")
        self.test_directory1.mkdir()

        # Minimum args list
        self.args_list = ["run.py",
                          "import",
                          "Actinobacteriophage",
                          str(self.test_directory1),
                          str(self.test_filepath1)]


    def tearDown(self):
        shutil.rmtree(self.test_directory1)




    def test_parse_args_1(self):
        """Verify args when minimum args_list is provided."""
        args = import_genome.parse_args(self.args_list)
        with self.subTest():
            self.assertEqual(args.database, "Actinobacteriophage")
        with self.subTest():
            self.assertEqual(args.input_folder, self.test_directory1)
        with self.subTest():
            self.assertEqual(args.import_table, self.test_filepath1)
        with self.subTest():
            self.assertEqual(args.genome_id_field, "_organism_name")
        with self.subTest():
            self.assertFalse(args.prod_run)
        with self.subTest():
            self.assertEqual(args.eval_mode, "final")
        with self.subTest():
            self.assertEqual(args.description_field, "product")
        with self.subTest():
            self.assertEqual(args.output_folder, Path("/tmp/"))
        with self.subTest():
            self.assertEqual(args.log_file, "import.log")
        with self.subTest():
            self.assertFalse(args.interactive)


    @patch("sys.exit")
    def test_parse_args_2(self, sys_exit_mock):
        """Verify sys exit when too few arguments are provided."""
        self.args_list.pop()
        args = import_genome.parse_args(self.args_list)
        self.assertTrue(sys_exit_mock.called)


    def test_parse_args_3(self):
        """Verify parsed args when all args are explicitly provided
        using short name."""
        output_folder = "/path/to/output/"
        log_file = "logfile.txt"
        self.args_list.extend(["-g", "FILENAME",
                               "-p",
                               "-e", "DRAFT",
                               "-d", "FUNCTION",
                               "-o", output_folder,
                               "-l", log_file,
                               "-i"
                               ])
        args = import_genome.parse_args(self.args_list)
        with self.subTest():
            self.assertEqual(args.genome_id_field, "filename")
        with self.subTest():
            self.assertTrue(args.prod_run)
        with self.subTest():
            self.assertEqual(args.eval_mode, "draft")
        with self.subTest():
            self.assertEqual(args.description_field, "function")
        with self.subTest():
            self.assertEqual(args.output_folder, Path(output_folder))
        with self.subTest():
            self.assertEqual(args.log_file, log_file)
        with self.subTest():
            self.assertTrue(args.interactive)


    def test_parse_args_4(self):
        """Verify parsed args when all args are explicitly provided
        using long name."""
        output_folder = "/path/to/output/"
        log_file = "logfile.txt"
        self.args_list.extend(["--genome_id_field", "FILENAME",
                               "--prod_run",
                               "--eval_mode", "DRAFT",
                               "--description_field", "FUNCTION",
                               "--output_folder", output_folder,
                               "--log_file", log_file,
                               "--interactive"
                               ])
        args = import_genome.parse_args(self.args_list)
        with self.subTest():
            self.assertEqual(args.genome_id_field, "filename")
        with self.subTest():
            self.assertTrue(args.prod_run)
        with self.subTest():
            self.assertEqual(args.eval_mode, "draft")
        with self.subTest():
            self.assertEqual(args.description_field, "function")
        with self.subTest():
            self.assertEqual(args.output_folder, Path(output_folder))
        with self.subTest():
            self.assertEqual(args.log_file, log_file)
        with self.subTest():
            self.assertTrue(args.interactive)




class TestImportGenomeMain3(unittest.TestCase):

    def setUp(self):
        # The test_db is only needed to test shema compatibility.
        # Otherwise, Actinobacteriophage is sufficient.
        create_new_db(schema_filepath, db, user, pwd)

        self.import_table = Path(test_file_dir, "test_import_table_1.csv")
        self.base_dir = Path(test_root_dir, "test_folder")
        self.base_dir.mkdir()

        self.input_folder = Path(self.base_dir, "input_folder")
        self.output_folder = Path(self.base_dir, "output_folder")
        self.log_file = Path(self.base_dir, "test_log.txt")

        self.engine = sqlalchemy.create_engine(engine_string1, echo=False)

        self.args_list = ["run.py",
                          "import",
                          db,
                          str(self.input_folder),
                          str(self.import_table),
                          "-g", "FILENAME",
                          "-p",
                          "-e", "DRAFT",
                          "-d", "FUNCTION",
                          "-o", str(self.output_folder),
                          "-l", str(self.log_file)
                          ]

    def tearDown(self):
        shutil.rmtree(self.base_dir)
        self.engine.dispose()
        remove_db(db, user, pwd)




    @patch("pdm_utils.pipelines.import_genome.data_io")
    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    def test_main_1(self, ctd_mock, data_io_mock):
        """Verify that correct args calls data_io."""
        insert_data_into_version_table(db, user, pwd, version_table_data)
        self.input_folder.mkdir()
        self.output_folder.mkdir()
        ctd_mock.return_value = self.engine
        import_genome.main(self.args_list)
        self.assertTrue(data_io_mock.called)


    @patch("pdm_utils.pipelines.import_genome.data_io")
    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    @patch("sys.exit")
    def test_main_2(self, sys_exit_mock, ctd_mock, data_io_mock):
        """Verify that invalid input folder calls sys exit."""
        insert_data_into_version_table(db, user, pwd, version_table_data)
        self.output_folder.mkdir()
        ctd_mock.return_value = self.engine
        import_genome.main(self.args_list)
        self.assertTrue(sys_exit_mock.called)


    @patch("pdm_utils.pipelines.import_genome.data_io")
    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    @patch("sys.exit")
    def test_main_3(self, sys_exit_mock, ctd_mock, data_io_mock):
        """Verify that invalid import file calls sys exit."""
        insert_data_into_version_table(db, user, pwd, version_table_data)
        self.input_folder.mkdir()
        self.output_folder.mkdir()
        self.args_list[4] = ""
        ctd_mock.return_value = self.engine
        import_genome.main(self.args_list)
        self.assertTrue(sys_exit_mock.called)


    @patch("pdm_utils.pipelines.import_genome.data_io")
    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    @patch("pdm_utils.functions.basic.make_new_dir")
    @patch("sys.exit")
    def test_main_4(self, sys_exit_mock, mnd_mock, ctd_mock, data_io_mock):
        """Verify that invalid output folder calls sys exit."""
        insert_data_into_version_table(db, user, pwd, version_table_data)
        self.input_folder.mkdir()
        # Need to provide filename in a valid directory to create log file
        # since output folder is invalid.
        mnd_mock.return_value = Path(self.input_folder, "temp")
        ctd_mock.return_value = self.engine
        import_genome.main(self.args_list)
        self.assertTrue(sys_exit_mock.called)


    @patch("pdm_utils.pipelines.import_genome.data_io")
    @patch("pdm_utils.functions.mysqldb.check_schema_compatibility")
    @patch("sys.exit")
    @patch("getpass.getpass")
    def test_main_5(self, getpass_mock, sys_exit_mock, csc_mock, data_io_mock):
        """Verify that invalid database calls sys exit."""
        self.input_folder.mkdir()
        self.output_folder.mkdir()
        self.args_list[2] = "Actinobacteriophage_x"
        # Assumes that mysqldb.get_engine(attempts=5)
        getpass_mock.side_effect = [user, pwd,
                                    user, pwd,
                                    user, pwd,
                                    user, pwd,
                                    user, pwd]
        import_genome.main(self.args_list)
        self.assertTrue(sys_exit_mock.called)


    @patch("pdm_utils.pipelines.import_genome.data_io")
    @patch("sys.exit")
    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    def test_main_6(self, ctd_mock, sys_exit_mock, data_io_mock):
        """Verify that invalid database schema version calls sys exit."""
        new_version_table_data = {"Version":1, "SchemaVersion":0}
        insert_data_into_version_table(db, user, pwd, new_version_table_data)
        self.input_folder.mkdir()
        self.output_folder.mkdir()
        ctd_mock.return_value = self.engine
        import_genome.main(self.args_list)
        self.assertTrue(sys_exit_mock.called)




class TestImportGenomeMain4(unittest.TestCase):

    def setUp(self):

        self.eval_data_dict = {"eval_mode": "custom_eval_mode",
                                   "eval_flag_dict": {"a":1}}

        self.table_structure_dict = constants.IMPORT_TABLE_STRUCTURE

        self.test_import_table_1 = Path(test_file_dir, "test_import_table_1.csv")

        # Valid data dictionary.
        self.data_dict1 = {}
        self.data_dict1["type"] = "replace"
        self.data_dict1["phage_id"] = "Trixie"
        self.data_dict1["description_field"] = "product"
        self.data_dict1["eval_mode"] = "final"
        self.data_dict1["host_genus"] = "retrieve"
        self.data_dict1["cluster"] = "retain"

        # Valid data dictionary.
        self.data_dict2 = {}
        self.data_dict2["type"] = "replace"
        self.data_dict2["phage_id"] = "L5"
        self.data_dict2["description_field"] = "product"
        self.data_dict2["eval_mode"] = "final"
        self.data_dict2["host_genus"] = "retrieve"
        self.data_dict2["cluster"] = "retain"


        self.dict1_dict2 = [self.data_dict1, self.data_dict2]




    def test_prepare_tickets_1(self):
        """Verify dictionary is returned from a correct import table file."""

        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table_1,
                        eval_data_dict=self.eval_data_dict,
                        description_field="product",
                        table_structure_dict=self.table_structure_dict)
        self.assertEqual(len(tkt_dict.keys()), 2)


    # Patch so that a variety of different types of files don't need
    # to be created just to test this function.
    @patch("pdm_utils.functions.basic.retrieve_data_dict")
    def test_prepare_tickets_2(self, mock_retrieve_tickets):
        """Verify dictionary is returned from two correct
        import data dictionaries."""

        mock_retrieve_tickets.return_value = self.dict1_dict2
        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table_1,
                        eval_data_dict=self.eval_data_dict,
                        description_field="product",
                        table_structure_dict=self.table_structure_dict)
        self.assertEqual(len(tkt_dict.keys()), 2)


    @patch("pdm_utils.functions.basic.retrieve_data_dict")
    def test_prepare_tickets_3(self, mock_retrieve_tickets):
        """Verify no dictionary is returned from one correct
        and one incorrect import data dictionaries."""

        self.data_dict2.pop("phage_id")
        mock_retrieve_tickets.return_value = self.dict1_dict2
        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table_1,
                        eval_data_dict=self.eval_data_dict,
                        description_field="product",
                        table_structure_dict=self.table_structure_dict)
        self.assertIsNone(tkt_dict)




    @patch("pdm_utils.functions.basic.retrieve_data_dict")
    def test_prepare_tickets_4(self, mock_retrieve_tickets):
        """Verify no dictionary is returned from one correct
        data dictionary and one correct data dictionary with
        duplicated phage_id."""

        self.data_dict2["phage_id"] = "Trixie"
        mock_retrieve_tickets.return_value = self.dict1_dict2
        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table_1,
                        eval_data_dict=self.eval_data_dict,
                        description_field="product",
                        table_structure_dict=self.table_structure_dict)
        self.assertIsNone(tkt_dict)


    @patch("pdm_utils.functions.basic.retrieve_data_dict")
    def test_prepare_tickets_5(self, mock_retrieve_tickets):
        """Verify no dictionary is returned from one correct
        data dictionary and one incorrect data dictionary with
        invalid ticket type."""

        self.data_dict2["type"] = "invalid"
        mock_retrieve_tickets.return_value = self.dict1_dict2
        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table_1,
                        eval_data_dict=self.eval_data_dict,
                        description_field="product",
                        table_structure_dict=self.table_structure_dict)
        self.assertIsNone(tkt_dict)


    @patch("pdm_utils.functions.basic.retrieve_data_dict")
    def test_prepare_tickets_6(self, mock_retrieve_tickets):
        """Verify no dictionary is returned from one correct
        data dictionary and one incorrect data dictionary with
        invalid cluster (set to 'parse')."""

        self.data_dict2["cluster"] = "parse"
        mock_retrieve_tickets.return_value = self.dict1_dict2
        tkt_dict = import_genome.prepare_tickets(
                        import_table_file=self.test_import_table_1,
                        eval_data_dict=self.eval_data_dict,
                        description_field="product",
                        table_structure_dict=self.table_structure_dict)
        self.assertIsNone(tkt_dict)




class TestImportGenomeMain5(unittest.TestCase):


    def setUp(self):

        self.flat_file_l5 = Path(test_file_dir, "test_flat_file_1.gb")
        self.flat_file_trixie = Path(test_file_dir, "test_flat_file_6.gb")

        self.engine = sqlalchemy.create_engine(engine_string2, echo=False)

        self.eval_flags = {}
        self.eval_flags["check_seq"] = True
        self.eval_flags["check_id_typo"] = True
        self.eval_flags["check_host_typo"] = True
        self.eval_flags["check_author"] = True
        self.eval_flags["check_trna"] = True
        self.eval_flags["check_gene"] = True
        self.eval_flags["check_locus_tag"] = True
        self.eval_flags["check_description"] = True
        self.eval_flags["check_description_field"] = True
        self.eval_flags["import_locus_tag"] = True

        self.data_dict1 = {}
        self.data_dict1["type"] = "replace"
        self.data_dict1["phage_id"] = "L5"
        self.data_dict1["description_field"] = "product"
        self.data_dict1["eval_mode"] = "final"
        self.data_dict1["host_genus"] = "Mycobacterium"
        self.data_dict1["cluster"] = "A"
        self.data_dict1["annotation_status"] = "draft"

        self.tkt1 = ticket.ImportTicket()
        self.tkt1.id = 1
        self.tkt1.phage_id = "L5"
        self.tkt1.eval_mode = "final"
        self.tkt1.description_field = "product"
        self.tkt1.eval_flags = self.eval_flags
        self.tkt1.data_dict = self.data_dict1

        self.data_dict2 = {}
        self.data_dict2["type"] = "replace"
        self.data_dict2["id"] = 1
        self.data_dict2["phage_id"] = "Trixie"
        self.data_dict2["description_field"] = "product"
        self.data_dict2["eval_mode"] = "final"
        self.data_dict2["host_genus"] = "Gordonia"
        self.data_dict2["cluster"] = "B"

        self.tkt2 = ticket.ImportTicket()
        self.tkt2.id = 2
        self.tkt2.phage_id = "Trixie"
        self.tkt2.eval_mode = "final"
        self.tkt2.description_field = "product"
        self.tkt2.eval_flags = self.eval_flags
        self.tkt2.data_dict = self.data_dict2

        # To test how log files are managed:
        self.base_dir = Path(test_root_dir,"test_import")
        self.base_dir.mkdir()
        self.success_path = Path(self.base_dir, "success_folder")
        self.success_path.mkdir()
        self.fail_path = Path(self.base_dir, "fail_folder")
        self.fail_path.mkdir()
        self.paths_dict = {"success": self.success_path,
                           "fail": self.fail_path}

    def tearDown(self):
        self.engine.dispose()
        shutil.rmtree(self.base_dir)




    def test_process_files_and_tickets_1(self):
        """Verify correct output using:
        no files,
        no unmatched tickets."""
        ticket_dict = {}
        files = []
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, engine=self.engine,
                            prod_run=False, genome_id_field="_organism_name")
        success_ticket_list = results_tuple[0]
        failed_ticket_list = results_tuple[1]
        success_filename_list = results_tuple[2]
        failed_filename_list = results_tuple[3]
        evaluation_dict = results_tuple[4]
        with self.subTest():
            self.assertEqual(len(success_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(success_filename_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_filename_list), 0)
        with self.subTest():
            self.assertEqual(len(evaluation_dict.keys()), 0)

    # Patching to avoid an attempt to add data to the database.
    @patch("pdm_utils.pipelines.import_genome.import_into_db")
    # Patching glp since the bundled data is incomplete,
    # bndl.check_for_errors() will always throw an error.
    @patch("pdm_utils.pipelines.import_genome.get_logfile_path")
    def test_process_files_and_tickets_2(self, glp_mock, import_into_db_mock):
        """Verify correct output using:
        two files with no tickets,
        unsuccessful import,
        no unmatched tickets.
        Also testing that file-specific log files are generated."""
        ticket_dict = {}
        files = [self.flat_file_l5, self.flat_file_trixie]
        log_l5 = Path(self.success_path, self.flat_file_l5.stem)
        log_trixie = Path(self.success_path, self.flat_file_trixie.stem)
        glp_mock.side_effect = [log_l5, log_trixie]
        import_into_db_mock.side_effect = [False, False]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, engine=self.engine,
                            prod_run=False, genome_id_field="_organism_name",
                            log_folder_paths_dict=self.paths_dict)
        success_ticket_list = results_tuple[0]
        failed_ticket_list = results_tuple[1]
        success_filename_list = results_tuple[2]
        failed_filename_list = results_tuple[3]
        evaluation_dict = results_tuple[4]
        s_count = count_contents(self.success_path)
        f_count = count_contents(self.fail_path)
        with self.subTest():
            self.assertEqual(len(success_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(success_filename_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_filename_list), 2)
        with self.subTest():
            self.assertEqual(evaluation_dict.keys(), set([1, 2]))
        with self.subTest():
            self.assertEqual(s_count, 2)
        with self.subTest():
            self.assertEqual(f_count, 0)


    # Patching to avoid an attempt to add data to the database.
    @patch("pdm_utils.pipelines.import_genome.import_into_db")
    def test_process_files_and_tickets_3(self, import_into_db_mock):
        """Verify correct output using:
        two files with matched tickets,
        successful import,
        no unmatched tickets.
        Also testing that NO file-specific log files are generated when
        no paths_dict is provided."""
        ticket_dict = {self.tkt1.phage_id: self.tkt1,
                       self.tkt2.phage_id: self.tkt2}
        files = [self.flat_file_l5, self.flat_file_trixie]
        import_into_db_mock.side_effect = [True, True]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, engine=self.engine,
                            prod_run=False, genome_id_field="_organism_name",
                            log_folder_paths_dict=None)
        success_ticket_list = results_tuple[0]
        failed_ticket_list = results_tuple[1]
        success_filename_list = results_tuple[2]
        failed_filename_list = results_tuple[3]
        evaluation_dict = results_tuple[4]
        s_count = count_contents(self.success_path)
        f_count = count_contents(self.fail_path)
        with self.subTest():
            self.assertEqual(len(success_ticket_list), 2)
        with self.subTest():
            self.assertEqual(len(failed_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(success_filename_list), 2)
        with self.subTest():
            self.assertEqual(len(failed_filename_list), 0)
        with self.subTest():
            self.assertEqual(evaluation_dict.keys(), set([1, 2]))
        with self.subTest():
            self.assertEqual(s_count, 0)
        with self.subTest():
            self.assertEqual(f_count, 0)


    # Patching to avoid an attempt to add data to the database.
    @patch("pdm_utils.pipelines.import_genome.import_into_db")
    def test_process_files_and_tickets_4(self, import_into_db_mock):
        """Verify correct output using:
        two files with matched tickets,
        unsuccessful import,
        no unmatched tickets."""
        ticket_dict = {self.tkt1.phage_id: self.tkt1,
                       self.tkt2.phage_id: self.tkt2}
        files = [self.flat_file_l5, self.flat_file_trixie]
        import_into_db_mock.side_effect = [False, False]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, engine=self.engine,
                            prod_run=False, genome_id_field="_organism_name")
        success_ticket_list = results_tuple[0]
        failed_ticket_list = results_tuple[1]
        success_filename_list = results_tuple[2]
        failed_filename_list = results_tuple[3]
        evaluation_dict = results_tuple[4]
        with self.subTest():
            self.assertEqual(len(success_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_ticket_list), 2)
        with self.subTest():
            self.assertEqual(len(success_filename_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_filename_list), 2)
        with self.subTest():
            self.assertEqual(evaluation_dict.keys(), set([1, 2]))

    def test_process_files_and_tickets_5(self):
        """Verify correct output using:
        no files,
        two unmatched tickets.
        Also testing that ticket-specific log files are generated when
        paths_dict is provided."""
        ticket_dict = {self.tkt1.phage_id: self.tkt1,
                       self.tkt2.phage_id: self.tkt2}
        files = []
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, engine=self.engine,
                            prod_run=False, genome_id_field="_organism_name",
                            log_folder_paths_dict=self.paths_dict)
        success_ticket_list = results_tuple[0]
        failed_ticket_list = results_tuple[1]
        success_filename_list = results_tuple[2]
        failed_filename_list = results_tuple[3]
        evaluation_dict = results_tuple[4]
        s_count = count_contents(self.success_path)
        f_count = count_contents(self.fail_path)
        with self.subTest():
            self.assertEqual(len(success_ticket_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_ticket_list), 2)
        with self.subTest():
            self.assertEqual(len(success_filename_list), 0)
        with self.subTest():
            self.assertEqual(len(failed_filename_list), 0)
        with self.subTest():
            self.assertEqual(evaluation_dict.keys(), set([1, 2]))
        with self.subTest():
            self.assertEqual(s_count, 0)
        with self.subTest():
            self.assertEqual(f_count, 2)


    def test_process_files_and_tickets_6(self):
        """Same test as test_process_files_and_tickets_5, except that
        NO ticket-specific log files are generated when
        paths_dict is NOT provided."""
        ticket_dict = {self.tkt1.phage_id: self.tkt1,
                       self.tkt2.phage_id: self.tkt2}
        files = []
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, engine=self.engine,
                            prod_run=False, genome_id_field="_organism_name",
                            log_folder_paths_dict=None)
        s_count = count_contents(self.success_path)
        f_count = count_contents(self.fail_path)
        with self.subTest():
            self.assertEqual(s_count, 0)
        with self.subTest():
            self.assertEqual(f_count, 0)



    # Patching to avoid an attempt to add data to the database.
    @patch("pdm_utils.pipelines.import_genome.import_into_db")
    def test_process_files_and_tickets_7(self, import_into_db_mock):
        """Verify correct output using:
        one file matched to ticket with successful import,
        one file unmatched to ticket with unsuccessful import,
        one unmatched ticket."""
        self.tkt2.phage_id = "Trixie_x"
        ticket_dict = {self.tkt1.phage_id: self.tkt1,
                       self.tkt2.phage_id: self.tkt2}
        files = [self.flat_file_l5, self.flat_file_trixie]
        import_into_db_mock.side_effect = [True, False]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, engine=self.engine,
                            prod_run=False, genome_id_field="_organism_name")
        success_ticket_list = results_tuple[0]
        failed_ticket_list = results_tuple[1]
        success_filename_list = results_tuple[2]
        failed_filename_list = results_tuple[3]
        evaluation_dict = results_tuple[4]
        with self.subTest():
            self.assertEqual(len(success_ticket_list), 1)
        with self.subTest():
            self.assertEqual(len(failed_ticket_list), 1)
        with self.subTest():
            self.assertEqual(len(success_filename_list), 1)
        with self.subTest():
            self.assertEqual(len(failed_filename_list), 1)
        with self.subTest():
            self.assertEqual(evaluation_dict.keys(), set([1, 2, 3]))


    # Patching to avoid an attempt to add data to the database.
    @patch("pdm_utils.functions.mysqldb.execute_transaction")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_process_files_and_tickets_8(self, ask_mock, execute_mock):
        """Verify correct output using:
        one file with matched ticket,
        one error evaluation in ticket (ensuring at least one error
        in all evaluations),
        interactive = False."""
        self.tkt1.evaluations = [eval.Eval(status="error")]
        ticket_dict = {self.tkt1.phage_id: self.tkt1}
        files = [self.flat_file_l5]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, engine=self.engine,
                            prod_run=True, genome_id_field="_organism_name",
                            interactive=False)
        evaluation_dict = results_tuple[4]
        error_count = count_status_from_dict(evaluation_dict, "error")
        with self.subTest():
            self.assertFalse(ask_mock.called)
        with self.subTest():
            self.assertFalse(execute_mock.called)
        with self.subTest():
            self.assertTrue(error_count > 0)


    # Patching to avoid an attempt to add data to the database.
    @patch("pdm_utils.functions.mysqldb.execute_transaction")
    # Since interactive=True and there are evals with 'error' status,
    # need to patch input to proceed.
    @patch("builtins.input")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_process_files_and_tickets_9(self, ask_mock, input_mock,
                                         execute_mock):
        """Verify correct output using:
        one file with matched ticket,
        one error evaluation in ticket (ensuring at least one error
        in all evaluations),
        interactive = True,
        with no 'warning' status corrections."""
        self.tkt1.evaluations = [eval.Eval(status="error")]
        ask_mock.return_value = True
        ticket_dict = {self.tkt1.phage_id: self.tkt1}
        files = [self.flat_file_l5]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, engine=self.engine,
                            prod_run=True, genome_id_field="_organism_name",
                            interactive=True)
        evaluation_dict = results_tuple[4]
        error_count = count_status_from_dict(evaluation_dict, "error")
        warning_count = count_status_from_dict(evaluation_dict, "warning")
        with self.subTest():
            self.assertTrue(ask_mock.called)
        with self.subTest():
            self.assertFalse(execute_mock.called)
        with self.subTest():
            self.assertTrue(error_count > 0)
        with self.subTest():
            self.assertTrue(warning_count > 0)


    # Patching to avoid an attempt to add data to the database.
    @patch("pdm_utils.functions.mysqldb.execute_transaction")
    # Since interactive=True and there are evals with 'error' status,
    # need to patch input to proceed.
    @patch("builtins.input")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_process_files_and_tickets_10(self, ask_mock, input_mock,
                                         execute_mock):
        """Verify correct output using:
        one file with matched ticket,
        one error evaluation in ticket (ensuring at least one error
        in all evaluations),
        interactive = True,
        with all 'warning' status changes to 'error'."""
        self.tkt1.evaluations = [eval.Eval(status="error")]
        ask_mock.return_value = False
        ticket_dict = {self.tkt1.phage_id: self.tkt1}
        files = [self.flat_file_l5]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, engine=self.engine,
                            prod_run=True, genome_id_field="_organism_name",
                            interactive=True)
        evaluation_dict = results_tuple[4]
        error_count = count_status_from_dict(evaluation_dict, "error")
        warning_count = count_status_from_dict(evaluation_dict, "warning")
        with self.subTest():
            self.assertTrue(ask_mock.called)
        with self.subTest():
            self.assertFalse(execute_mock.called)
        with self.subTest():
            self.assertEqual(warning_count, 0)
        with self.subTest():
            self.assertTrue(error_count > 0)


    # Patching to avoid an attempt to add data to the database.
    @patch("pdm_utils.functions.mysqldb.execute_transaction")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    @patch("pdm_utils.pipelines.import_genome.run_checks")
    @patch("pdm_utils.pipelines.import_genome.prepare_bundle")
    def test_process_files_and_tickets_11(self, prep_mock, run_checks_mock,
                                          ask_mock, execute_mock):
        """Verify correct output using:
        two files with matched tickets,
        one warning evaluation in ticket (ensuring at least one warning
        in all evaluations),
        interactive = True,
        with all 'warning' status corrections to 'error' in first file
        but no corrections in second file."""
        self.tkt1.evaluations = [eval.Eval(status="warning")]
        self.tkt2.evaluations = [eval.Eval(status="warning")]
        ticket_dict = {self.tkt1.phage_id: self.tkt1,
                       self.tkt2.phage_id: self.tkt2}
        files = [self.flat_file_l5, self.flat_file_trixie]
        gnm1 = genome.Genome()
        bndl1 = bundle.Bundle()
        bndl1.id = 1
        bndl1.ticket = self.tkt1
        bndl1.genome_dict["flat_file"] = gnm1

        gnm2 = genome.Genome()
        bndl2 = bundle.Bundle()
        bndl2.id = 2
        bndl2.ticket = self.tkt2
        bndl2.genome_dict["flat_file"] = gnm2

        prep_mock.side_effect = [bndl1, bndl2]
        ask_mock.side_effect = [False, True]
        results_tuple = import_genome.process_files_and_tickets(ticket_dict,
                            files, engine=self.engine,
                            prod_run=True, genome_id_field="_organism_name",
                            interactive=True)
        with self.subTest():
            self.assertTrue(ask_mock.called)
        with self.subTest():
            self.assertTrue(execute_mock.called)
        with self.subTest():
            self.assertEqual(bndl1._errors, 1)
        with self.subTest():
            self.assertEqual(bndl2._errors, 0)




class TestImportGenomeMain6(unittest.TestCase):

    def setUp(self):

        self.base_dir = Path(test_root_dir,"test_import")
        self.base_dir.mkdir()

        self.genome_folder = Path(self.base_dir, "genomes")

        self.flat_file1 = Path(self.genome_folder, "flat_file1.txt")
        self.flat_file2 = Path(self.genome_folder, "flat_file2.txt")

        self.valid_import_table_file = Path(test_file_dir,
                                            "test_import_table_1.csv")
        self.invalid_import_table_file = Path(test_file_dir,
                                            "test_import_table_2.csv")
        self.output_folder = Path(self.base_dir, import_genome.RESULTS_FOLDER)


        self.engine = sqlalchemy.create_engine(engine_string1, echo=False)

        self.tkt1 = ticket.ImportTicket()
        self.tkt1.phage_id = "L5"
        self.tkt1.host_genus = "Mycobacterium"

        self.tkt2 = ticket.ImportTicket()
        self.tkt2.phage_id = "Trixie"
        self.tkt2.host_genus = "Mycobacterium"

        self.tkt_dict1 = {"type": "add", "phage_id": "L5",
                          "host_genus": "Mycobacterium"}
        self.tkt_dict2 = {"type": "replace", "phage_id": "Trixie",
                          "host_genus": "Mycobacterium"}

        self.date = time.strftime("%Y%m%d")

        self.exp_success = Path(self.output_folder, "success")
        self.exp_success_tkt_table = Path(self.exp_success, "import_tickets.csv")
        self.exp_success_genomes = Path(self.exp_success, "genomes")
        self.exp_success_logs = Path(self.exp_success, "logs")

        self.exp_fail = Path(self.output_folder, "fail")
        self.exp_fail_tkt_table = Path(self.exp_fail, "import_tickets.csv")
        self.exp_fail_genomes = Path(self.exp_fail, "genomes")
        self.exp_fail_logs = Path(self.exp_fail, "logs")


    def tearDown(self):
        shutil.rmtree(self.base_dir)
        self.engine.dispose()


    @patch("sys.exit")
    @patch("pdm_utils.pipelines.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.mysqldb.get_phage_table_count")
    def test_data_io_1(self, get_count_mock, pft_mock, sys_exit_mock):
        """Verify data_io runs correctly when there are no errors."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        self.flat_file1.touch()
        self.flat_file2.touch()

        success_ticket_list = [self.tkt_dict1]
        failed_ticket_list = [self.tkt_dict1, self.tkt_dict2]
        success_filename_list = [self.flat_file1]
        failed_filename_list = [self.flat_file2]
        evaluation_dict = {}
        get_count_mock.return_value = 0
        pft_mock.return_value = (success_ticket_list,
                                 failed_ticket_list,
                                 success_filename_list,
                                 failed_filename_list,
                                 evaluation_dict)
        import_genome.data_io(engine=self.engine,
            genome_folder=self.genome_folder,
            import_table_file=self.valid_import_table_file,
            output_folder=self.output_folder)

        exp_success_tkts = []
        with open(self.exp_success_tkt_table,'r') as file:
            file_reader = csv.DictReader(file)
            for dict in file_reader:
                exp_success_tkts.append(dict)

        exp_fail_tkts = []
        with open(self.exp_fail_tkt_table,'r') as file:
            file_reader = csv.DictReader(file)
            for dict in file_reader:
                exp_fail_tkts.append(dict)

        input_genomes_count = 0
        for item in self.genome_folder.iterdir():
            input_genomes_count += 1

        success_genomes_count = 0
        for item in self.exp_success_genomes.iterdir():
            success_genomes_count += 1
        fail_genomes_count = 0
        for item in self.exp_fail_genomes.iterdir():
            fail_genomes_count += 1

        # Unable to check how many log files are present -
        # the log folders don't exist because process_files_and_tickets()
        # was patched, so no log files created, so data_io() removes the folders.
        # success_log_files_count = 0
        # for item in self.exp_success_logs.iterdir():
        #     success_log_files_count += 1
        # fail_log_files_count = 0
        # for item in self.exp_fail_logs.iterdir():
        #     fail_log_files_count += 1

        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertFalse(sys_exit_mock.called)
        with self.subTest():
            self.assertEqual(len(exp_success_tkts), 1)
        with self.subTest():
            self.assertEqual(len(exp_success_tkts[0].keys()), 11)
        with self.subTest():
            self.assertEqual(len(exp_fail_tkts), 2)
        with self.subTest():
            self.assertEqual(len(exp_fail_tkts[0].keys()), 11)
        with self.subTest():
            self.assertEqual(input_genomes_count, 2)
        with self.subTest():
            self.assertEqual(success_genomes_count, 1)
        with self.subTest():
            self.assertEqual(fail_genomes_count, 1)
        # Since pft was patched, file-specific log files were not created,
        # so the log folders should be removed:
        with self.subTest():
            self.assertFalse(self.exp_fail_logs.is_dir())
        with self.subTest():
            self.assertFalse(self.exp_success_logs.is_dir())


    @patch("sys.exit")
    @patch("pdm_utils.pipelines.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.mysqldb.get_phage_table_count")
    def test_data_io_2(self, get_count_mock, pft_mock, sys_exit_mock):
        """Verify data_io is successful with
        success tickets but no success files, and
        fail tickets but no fail files."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        self.flat_file1.touch()
        self.flat_file2.touch()

        success_ticket_list = [self.tkt_dict1]
        failed_ticket_list = [self.tkt_dict1, self.tkt_dict2]
        success_filename_list = []
        failed_filename_list = []
        evaluation_dict = {}
        get_count_mock.return_value = 0
        pft_mock.return_value = (success_ticket_list,
                                 failed_ticket_list,
                                 success_filename_list,
                                 failed_filename_list,
                                 evaluation_dict)
        import_genome.data_io(engine=self.engine,
            genome_folder=self.genome_folder,
            import_table_file=self.valid_import_table_file,
            output_folder=self.output_folder)

        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertFalse(sys_exit_mock.called)
        with self.subTest():
            self.assertTrue(self.exp_success_tkt_table.exists())
        with self.subTest():
            self.assertTrue(self.exp_fail_tkt_table.exists())
        with self.subTest():
            self.assertFalse(self.exp_success_genomes.exists())
        with self.subTest():
            self.assertFalse(self.exp_fail_genomes.exists())


    @patch("sys.exit")
    @patch("pdm_utils.pipelines.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.mysqldb.get_phage_table_count")
    def test_data_io_3(self, get_count_mock, pft_mock, sys_exit_mock):
        """Verify data_io is successful with
        success files but no success tickets, and
        fail files but no fail tickets."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        self.flat_file1.touch()
        self.flat_file2.touch()
        success_ticket_list = []
        failed_ticket_list = []
        success_filename_list = [self.flat_file1]
        failed_filename_list = [self.flat_file2]
        evaluation_dict = {}
        get_count_mock.return_value = 0
        pft_mock.return_value = (success_ticket_list,
                                 failed_ticket_list,
                                 success_filename_list,
                                 failed_filename_list,
                                 evaluation_dict)
        import_genome.data_io(engine=self.engine,
            genome_folder=self.genome_folder,
            import_table_file=self.valid_import_table_file,
            output_folder=self.output_folder)
        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertFalse(sys_exit_mock.called)
        with self.subTest():
            self.assertFalse(self.exp_success_tkt_table.exists())
        with self.subTest():
            self.assertFalse(self.exp_fail_tkt_table.exists())
        with self.subTest():
            self.assertTrue(self.exp_success_genomes.exists())
        with self.subTest():
            self.assertTrue(self.exp_fail_genomes.exists())


    @patch("pdm_utils.pipelines.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.mysqldb.get_phage_table_count")
    @patch("sys.exit")
    def test_data_io_4(self, sys_exit_mock, get_count_mock, pft_mock):
        """Verify data_io is not successful when there
        are no files to process."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        get_count_mock.return_value = 0
        pft_mock.return_value = ([], [], [], [], {})
        import_genome.data_io(engine=self.engine,
            genome_folder=self.genome_folder,
            import_table_file=self.valid_import_table_file,
            output_folder=self.output_folder)
        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)


    @patch("sys.exit")
    @patch("pdm_utils.pipelines.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.mysqldb.get_phage_table_count")
    def test_data_io_5(self, get_count_mock, pft_mock, sys_exit_mock):
        """Verify data_io is not successful when there are no tickets to process."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        self.flat_file1.touch()
        self.flat_file2.touch()
        get_count_mock.return_value = 0
        pft_mock.return_value = ([], [], [], [], {})
        import_genome.data_io(engine=self.engine,
            genome_folder=self.genome_folder,
            import_table_file=self.invalid_import_table_file,
            output_folder=self.output_folder)
        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertTrue(sys_exit_mock.called)


    @patch("sys.exit")
    @patch("pdm_utils.pipelines.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.mysqldb.get_phage_table_count")
    def test_data_io_6(self, get_count_mock, pft_mock, sys_exit_mock):
        """Verify data_io runs correctly when there are
        successful tickets and genomes, and no failed tickets or genomes."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        self.flat_file1.touch()
        success_ticket_list = [self.tkt_dict1]
        failed_ticket_list = []
        success_filename_list = [self.flat_file1]
        failed_filename_list = []
        evaluation_dict = {}
        get_count_mock.return_value = 0
        pft_mock.return_value = (success_ticket_list,
                                 failed_ticket_list,
                                 success_filename_list,
                                 failed_filename_list,
                                 evaluation_dict)
        import_genome.data_io(engine=self.engine,
            genome_folder=self.genome_folder,
            import_table_file=self.valid_import_table_file,
            output_folder=self.output_folder)

        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertFalse(sys_exit_mock.called)
        with self.subTest():
            self.assertTrue(self.exp_success_tkt_table.exists())
        with self.subTest():
            self.assertTrue(self.exp_success_genomes.exists())
        with self.subTest():
            self.assertFalse(self.exp_fail.exists())


    @patch("sys.exit")
    @patch("pdm_utils.pipelines.import_genome.process_files_and_tickets")
    @patch("pdm_utils.functions.mysqldb.get_phage_table_count")
    def test_data_io_7(self, get_count_mock, pft_mock, sys_exit_mock):
        """Verify data_io runs correctly when there are
        failed tickets and genomes, and no successful tickets or genomes."""
        self.genome_folder.mkdir()
        self.output_folder.mkdir()
        self.flat_file1.touch()
        success_ticket_list = []
        failed_ticket_list = [self.tkt_dict1]
        success_filename_list = []
        failed_filename_list = [self.flat_file1]
        evaluation_dict = {}
        get_count_mock.return_value = 0
        pft_mock.return_value = (success_ticket_list,
                                 failed_ticket_list,
                                 success_filename_list,
                                 failed_filename_list,
                                 evaluation_dict)
        import_genome.data_io(engine=self.engine,
            genome_folder=self.genome_folder,
            import_table_file=self.valid_import_table_file,
            output_folder=self.output_folder)

        with self.subTest():
            self.assertTrue(pft_mock.called)
        with self.subTest():
            self.assertFalse(sys_exit_mock.called)
        with self.subTest():
            self.assertTrue(self.exp_fail_tkt_table.exists())
        with self.subTest():
            self.assertTrue(self.exp_fail_genomes.exists())
        with self.subTest():
            self.assertFalse(self.exp_success.exists())



class TestImportGenomeMain7(unittest.TestCase):


    def setUp(self):

        self.eval_warning1 = eval.Eval(status="warning", result="temp")
        self.eval_warning2 = eval.Eval(status="warning")
        self.eval_warning3 = eval.Eval(status="warning")
        self.eval_warning4 = eval.Eval(status="warning")
        self.eval_warning5 = eval.Eval(status="warning")
        self.eval_warning6 = eval.Eval(status="warning")
        self.eval_warning7 = eval.Eval(status="warning")
        self.eval_correct1 = eval.Eval(status="correct")
        self.eval_correct2 = eval.Eval(status="correct")
        self.eval_correct3 = eval.Eval(status="correct")
        self.eval_error1 = eval.Eval(status="error")


        self.tkt = ticket.ImportTicket()

        self.cds1 = cds.Cds()
        self.cds1.id = "CDS1"
        self.cds1.start = 1
        self.cds1.stop = 2
        self.cds1.orientation = "F"

        self.cds2 = cds.Cds()
        self.cds2.id = "CDS2"
        self.cds2.start = 3
        self.cds2.stop = 4
        self.cds2.orientation = "F"

        self.cds3 = cds.Cds()
        self.cds3.id = "CDS3"
        self.cds3.start = 5
        self.cds3.stop = 6
        self.cds3.orientation = "R"

        self.cds4 = cds.Cds()
        self.cds4.id = "CDS4"
        self.cds4.start = 7
        self.cds4.stop = 8
        self.cds4.orientation = "R"



        self.src1 = source.Source()
        self.src2 = source.Source()
        self.src3 = source.Source()

        self.gnm1 = genome.Genome()
        self.gnm2 = genome.Genome()
        self.gnm3 = genome.Genome()

        self.genome_pair1 = genomepair.GenomePair()
        self.genome_pair2 = genomepair.GenomePair()
        self.genome_pair3 = genomepair.GenomePair()

        self.bndl = bundle.Bundle()




    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluation_1(self, ask_mock):
        """Verify values and no need for user input with:
        interactive=False and 'warning' eval."""
        exit, correct = import_genome.review_evaluation(self.eval_warning1,
                                                        interactive=False)
        with self.subTest():
            self.assertEqual(self.eval_warning1.status, "error")
        with self.subTest():
            self.assertFalse(exit)
        with self.subTest():
            self.assertFalse(correct)
        with self.subTest():
            self.assertFalse(ask_mock.called)


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluation_2(self, ask_mock):
        """Verify values and need for user input with:
        interactive=True, 'warning' eval, response=True."""
        ask_mock.side_effect = [True]
        exit, correct = import_genome.review_evaluation(self.eval_warning1,
                                                        interactive=True)
        with self.subTest():
            self.assertEqual(self.eval_warning1.status, "warning")
        with self.subTest():
            self.assertFalse(exit)
        with self.subTest():
            self.assertFalse(correct)
        with self.subTest():
            self.assertTrue(ask_mock.called)
        with self.subTest():
            self.assertTrue(self.eval_warning1.result.endswith("temp"))


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluation_3(self, ask_mock):
        """Verify values and need for user input with:
        interactive=True, 'warning' eval, response=False."""
        ask_mock.side_effect = [False]
        exit, correct = import_genome.review_evaluation(self.eval_warning1,
                                                        interactive=True)
        with self.subTest():
            self.assertEqual(self.eval_warning1.status, "error")
        with self.subTest():
            self.assertFalse(exit)
        with self.subTest():
            self.assertFalse(correct)
        with self.subTest():
            self.assertTrue(ask_mock.called)
        with self.subTest():
            self.assertFalse(self.eval_warning1.result.startswith("Status"))
        with self.subTest():
            self.assertTrue(self.eval_warning1.result.endswith("manually."))


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluation_4(self, ask_mock):
        """Verify values and need for user input with:
        interactive=True, 'warning' eval, response=None."""
        ask_mock.side_effect = [None]
        exit, correct = import_genome.review_evaluation(self.eval_warning1,
                                                        interactive=True)
        with self.subTest():
            self.assertEqual(self.eval_warning1.status, "error")
        with self.subTest():
            self.assertTrue(exit)
        with self.subTest():
            self.assertFalse(correct)
        with self.subTest():
            self.assertTrue(ask_mock.called)
        with self.subTest():
            self.assertFalse(self.eval_warning1.result.startswith("Status"))
        with self.subTest():
            self.assertTrue(self.eval_warning1.result.endswith("exit."))


    @patch("builtins.input")
    def test_review_evaluation_5(self, input_mock):
        """Verify values and need for user input with:
        interactive=True, 'warning' eval, 3 invalid responses."""
        input_mock.side_effect = ["invalid", "invalid", "invalid"]
        exit, correct = import_genome.review_evaluation(self.eval_warning1,
                                                        interactive=True)
        with self.subTest():
            self.assertEqual(self.eval_warning1.status, "error")
        with self.subTest():
            self.assertTrue(exit)
        with self.subTest():
            self.assertFalse(correct)
        with self.subTest():
            self.assertTrue(input_mock.called)
        with self.subTest():
            self.assertFalse(self.eval_warning1.result.startswith("Status"))
        with self.subTest():
            self.assertTrue(self.eval_warning1.result.endswith("exit."))


    @patch("builtins.input")
    def test_review_evaluation_6(self, input_mock):
        """Verify values and need for user input with:
        interactive=True and 'error' eval."""
        exit, correct = import_genome.review_evaluation(self.eval_error1,
                                                        interactive=True)
        with self.subTest():
            self.assertEqual(self.eval_error1.status, "error")
        with self.subTest():
            self.assertFalse(exit)
        with self.subTest():
            self.assertFalse(correct)
        with self.subTest():
            self.assertTrue(input_mock.called)


    @patch("builtins.input")
    def test_review_evaluation_7(self, input_mock):
        """Verify values and no need for user input with:
        interactive=False and 'error' eval."""
        exit, correct = import_genome.review_evaluation(self.eval_error1,
                                                        interactive=False)
        with self.subTest():
            self.assertEqual(self.eval_error1.status, "error")
        with self.subTest():
            self.assertFalse(exit)
        with self.subTest():
            self.assertFalse(correct)
        with self.subTest():
            self.assertFalse(input_mock.called)


    @patch("builtins.input")
    def test_review_evaluation_8(self, input_mock):
        """Verify values and no need for user input with
        interactive=True and 'correct' eval."""
        exit, correct = import_genome.review_evaluation(self.eval_correct1,
                                                        interactive=True)
        with self.subTest():
            self.assertEqual(self.eval_correct1.status, "correct")
        with self.subTest():
            self.assertFalse(exit)
        with self.subTest():
            self.assertTrue(correct)
        with self.subTest():
            self.assertFalse(input_mock.called)




    @patch("pdm_utils.pipelines.import_genome.log_and_print")
    def test_review_evaluation_list_1(self, lp_mock):
        """Verify no need for user input if there are "
        "no 'warning' or 'error' evals."""
        eval_list = [self.eval_correct1, self.eval_correct2, self.eval_correct3]
        exit = import_genome.review_evaluation_list(eval_list, interactive=True)
        with self.subTest():
            self.assertEqual(self.eval_correct1.status, "correct")
        with self.subTest():
            self.assertEqual(self.eval_correct2.status, "correct")
        with self.subTest():
            self.assertEqual(self.eval_correct3.status, "correct")
        with self.subTest():
            self.assertFalse(exit)
        with self.subTest():
            self.assertTrue(lp_mock.called)


    @patch("pdm_utils.pipelines.import_genome.log_and_print")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluation_list_2(self, ask_mock, lp_mock):
        """Verify values when interactive=True, two 'warning' evals,
        and response=True."""
        ask_mock.side_effect = [True, True]
        eval_list = [self.eval_warning1, self.eval_correct2, self.eval_warning2]
        exit = import_genome.review_evaluation_list(eval_list, interactive=True)
        with self.subTest():
            self.assertEqual(self.eval_warning1.status, "warning")
        with self.subTest():
            self.assertEqual(self.eval_correct2.status, "correct")
        with self.subTest():
            self.assertEqual(self.eval_warning2.status, "warning")
        with self.subTest():
            self.assertFalse(exit)
        with self.subTest():
            self.assertFalse(lp_mock.called)


    @patch("pdm_utils.pipelines.import_genome.log_and_print")
    def test_review_evaluation_list_3(self, lp_mock):
        """Verify values when interactive=False, two 'warning' evals."""
        eval_list = [self.eval_warning1, self.eval_correct2, self.eval_warning2]
        exit = import_genome.review_evaluation_list(eval_list, interactive=False)
        with self.subTest():
            self.assertEqual(self.eval_warning1.status, "error")
        with self.subTest():
            self.assertEqual(self.eval_correct2.status, "correct")
        with self.subTest():
            self.assertEqual(self.eval_warning2.status, "error")
        with self.subTest():
            self.assertFalse(exit)
        with self.subTest():
            self.assertFalse(lp_mock.called)


    @patch("pdm_utils.pipelines.import_genome.log_and_print")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluation_list_4(self, ask_mock, lp_mock):
        """Verify values when interactive=True, two 'warning' evals,
        and response=False and True."""
        ask_mock.side_effect = [False, True]
        eval_list = [self.eval_warning1, self.eval_correct2, self.eval_warning2]
        exit = import_genome.review_evaluation_list(eval_list, interactive=True)
        with self.subTest():
            self.assertEqual(self.eval_warning1.status, "error")
        with self.subTest():
            self.assertEqual(self.eval_correct2.status, "correct")
        with self.subTest():
            self.assertEqual(self.eval_warning2.status, "warning")
        with self.subTest():
            self.assertFalse(exit)
        with self.subTest():
            self.assertFalse(lp_mock.called)


    @patch("pdm_utils.pipelines.import_genome.log_and_print")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_evaluation_list_5(self, ask_mock, lp_mock):
        """Verify values when interactive=True, three 'warning' evals,
        and response=True, exit."""
        ask_mock.side_effect = [True, None]
        eval_list = [self.eval_warning1, self.eval_warning2, self.eval_warning3]
        exit = import_genome.review_evaluation_list(eval_list, interactive=True)
        with self.subTest():
            self.assertEqual(self.eval_warning1.status, "warning")
        with self.subTest():
            self.assertEqual(self.eval_warning2.status, "error")
        with self.subTest():
            self.assertEqual(self.eval_warning3.status, "error")
        with self.subTest():
            self.assertTrue(exit)
        with self.subTest():
            self.assertFalse(lp_mock.called)




    @patch("pdm_utils.pipelines.import_genome.review_evaluation_list")
    def test_review_object_list_1(self, rel_mock):
        """Verify results for list of 0 features."""
        import_genome.review_object_list([], "CDS",
            ["id", "start", "stop", "orientation"], interactive=False)
        self.assertFalse(rel_mock.called)


    @patch("pdm_utils.pipelines.import_genome.review_evaluation_list")
    def test_review_object_list_2(self, rel_mock):
        """Verify results for list of 1 None feature."""
        import_genome.review_object_list([None], "CDS",
            ["id", "start", "stop", "orientation"], interactive=False)
        self.assertFalse(rel_mock.called)


    @patch("pdm_utils.pipelines.import_genome.review_evaluation_list")
    def test_review_object_list_3(self, rel_mock):
        """Verify results for list of three CDS features,
        each with 0 evaluations, interactive=False."""
        cds_features = [self.cds1, self.cds2, self.cds3]
        import_genome.review_object_list(cds_features, "CDS",
            ["id", "start", "stop", "orientation"], interactive=False)
        self.assertFalse(rel_mock.called)


    def test_review_object_list_4(self):
        """Verify results for list of three CDS features,
        one with three warning evaluations,
        and two with one warning evaluation, interactive=False."""
        self.cds1.evaluations = [self.eval_warning1,
                                 self.eval_warning2,
                                 self.eval_warning3]
        self.cds2.evaluations = [self.eval_warning4]
        self.cds3.evaluations = [self.eval_warning5]
        cds_features = [self.cds1, self.cds2, self.cds3]
        import_genome.review_object_list(cds_features, "CDS",
            ["id", "start", "stop", "orientation"], interactive=False)
        w1 = cds_features[0].evaluations[0]
        w2 = cds_features[0].evaluations[1]
        w3 = cds_features[0].evaluations[2]
        w4 = cds_features[1].evaluations[0]
        w5 = cds_features[2].evaluations[0]
        with self.subTest():
            self.assertEqual(w1.status, "error")
        with self.subTest():
            self.assertEqual(w2.status, "error")
        with self.subTest():
            self.assertEqual(w3.status, "error")
        with self.subTest():
            self.assertEqual(w4.status, "error")
        with self.subTest():
            self.assertEqual(w5.status, "error")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_object_list_5(self, ask_mock):
        """Verify results for list of three CDS features,
        one with three warning evaluations,
        and two with one warning evaluation, interactive=True,
        response=True, False, True, False, True."""
        ask_mock.side_effect = [True, False, True, False, True]
        self.cds1.evaluations = [self.eval_warning1,
                                 self.eval_warning2,
                                 self.eval_warning3]
        self.cds2.evaluations = [self.eval_warning4]
        self.cds3.evaluations = [self.eval_warning5]
        cds_features = [self.cds1, self.cds2, self.cds3]
        import_genome.review_object_list(cds_features, "CDS",
            ["id", "start", "stop", "orientation"], interactive=True)
        w1 = cds_features[0].evaluations[0]
        w2 = cds_features[0].evaluations[1]
        w3 = cds_features[0].evaluations[2]
        w4 = cds_features[1].evaluations[0]
        w5 = cds_features[2].evaluations[0]
        with self.subTest():
            self.assertEqual(w1.status, "warning")
        with self.subTest():
            self.assertEqual(w2.status, "error")
        with self.subTest():
            self.assertEqual(w3.status, "warning")
        with self.subTest():
            self.assertEqual(w4.status, "error")
        with self.subTest():
            self.assertEqual(w5.status, "warning")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_object_list_6(self, ask_mock):
        """Verify results for list of three CDS features,
        one with three warning evaluations,
        and two with one warning evaluation, interactive=True,
        response=True, exit."""
        ask_mock.side_effect = [True, None]
        self.cds1.evaluations = [self.eval_warning1,
                                 self.eval_warning2,
                                 self.eval_warning3]
        self.cds2.evaluations = [self.eval_warning4]
        self.cds3.evaluations = [self.eval_warning5]
        cds_features = [self.cds1, self.cds2, self.cds3]
        import_genome.review_object_list(cds_features, "CDS",
            ["id", "start", "stop", "orientation"], interactive=True)
        w1 = cds_features[0].evaluations[0]
        w2 = cds_features[0].evaluations[1]
        w3 = cds_features[0].evaluations[2]
        w4 = cds_features[1].evaluations[0]
        w5 = cds_features[2].evaluations[0]
        with self.subTest():
            self.assertEqual(w1.status, "warning")
        with self.subTest():
            self.assertEqual(w2.status, "error")
        with self.subTest():
            self.assertEqual(w3.status, "error")
        with self.subTest():
            self.assertEqual(w4.status, "error")
        with self.subTest():
            self.assertEqual(w5.status, "error")




    @patch("pdm_utils.pipelines.import_genome.review_evaluation_list")
    def test_review_bundled_objects_1(self, rel_mock):
        """Verify results when bundle has:
        no bundle evaluations,
        no ticket,
        no genomes in genome_dict,
        no genome_pairs in genome_pair_dict."""
        import_genome.review_bundled_objects(self.bndl, interactive=True)
        self.assertFalse(rel_mock.called)


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_bundled_objects_2(self, ask_mock):
        """Verify results when bundle has:
        one warning evaluation,
        interactive=True."""
        ask_mock.side_effect = [False]
        self.bndl.evaluations = [self.eval_warning1]
        import_genome.review_bundled_objects(self.bndl, interactive=True)
        self.assertEqual(self.bndl.evaluations[0].status, "error")


    def test_review_bundled_objects_3(self):
        """Verify results when bundle has:
        one warning evaluation,
        interactive=False."""
        self.bndl.evaluations = [self.eval_warning1]
        import_genome.review_bundled_objects(self.bndl, interactive=False)
        self.assertEqual(self.bndl.evaluations[0].status, "error")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_bundled_objects_4(self, ask_mock):
        """Verify results when bundle has:
        a ticket with one warning evaluation,
        interactive=False."""
        ask_mock.side_effect = [False]
        self.tkt.evaluations = [self.eval_warning1]
        self.bndl.ticket = self.tkt
        import_genome.review_bundled_objects(self.bndl, interactive=True)
        self.assertEqual(self.bndl.ticket.evaluations[0].status, "error")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_bundled_objects_5(self, ask_mock):
        """Verify results when bundle has:
        three genomes, two with one warning evaluation."""
        ask_mock.side_effect = [False, False]
        self.gnm1.evaluations = [self.eval_warning1]
        self.gnm2.evaluations = [self.eval_correct1]
        self.gnm3.evaluations = [self.eval_warning2]
        self.bndl.genome_dict["gnm1"] = self.gnm1
        self.bndl.genome_dict["gnm2"] = self.gnm2
        self.bndl.genome_dict["gnm3"] = self.gnm3
        import_genome.review_bundled_objects(self.bndl, interactive=True)
        evl_list1 = self.bndl.genome_dict["gnm1"].evaluations
        evl_list2 = self.bndl.genome_dict["gnm2"].evaluations
        evl_list3 = self.bndl.genome_dict["gnm3"].evaluations
        with self.subTest():
            self.assertEqual(evl_list1[0].status, "error")
        with self.subTest():
            self.assertEqual(evl_list2[0].status, "correct")
        with self.subTest():
            self.assertEqual(evl_list3[0].status, "error")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_bundled_objects_6(self, ask_mock):
        """Verify results when bundle has:
        one genome with three CDS features, two with one warning evaluation,
        one with two warning evaluations."""
        ask_mock.side_effect = [True, None]
        self.cds1.evaluations = [self.eval_warning1]
        self.cds2.evaluations = [self.eval_warning2, self.eval_warning3]
        self.cds3.evaluations = [self.eval_warning3]
        self.gnm1.cds_features = [self.cds1, self.cds2, self.cds3]
        self.bndl.genome_dict["gnm1"] = self.gnm1
        import_genome.review_bundled_objects(self.bndl, interactive=True)
        evl_list1 = self.bndl.genome_dict["gnm1"].cds_features[0].evaluations
        evl_list2 = self.bndl.genome_dict["gnm1"].cds_features[1].evaluations
        evl_list3 = self.bndl.genome_dict["gnm1"].cds_features[2].evaluations
        with self.subTest():
            self.assertEqual(evl_list1[0].status, "warning")
        with self.subTest():
            self.assertEqual(evl_list2[0].status, "error")
        with self.subTest():
            self.assertEqual(evl_list2[1].status, "error")
        with self.subTest():
            self.assertEqual(evl_list3[0].status, "error")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_bundled_objects_7(self, ask_mock):
        """Verify results when bundle has:
        one genome with three source features, two with one warning evaluation."""
        ask_mock.side_effect = [True, False]
        self.src1.evaluations = [self.eval_warning1]
        self.src2.evaluations = [self.eval_correct1]
        self.src3.evaluations = [self.eval_warning2]
        self.gnm1.source_features = [self.src1, self.src2, self.src3]
        self.bndl.genome_dict["gnm1"] = self.gnm1
        import_genome.review_bundled_objects(self.bndl, interactive=True)
        evl_list1 = self.bndl.genome_dict["gnm1"].source_features[0].evaluations
        evl_list2 = self.bndl.genome_dict["gnm1"].source_features[1].evaluations
        evl_list3 = self.bndl.genome_dict["gnm1"].source_features[2].evaluations
        with self.subTest():
            self.assertEqual(evl_list1[0].status, "warning")
        with self.subTest():
            self.assertEqual(evl_list2[0].status, "correct")
        with self.subTest():
            self.assertEqual(evl_list3[0].status, "error")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_bundled_objects_8(self, ask_mock):
        """Verify results when bundle has:
        three genome_pairs, two with one warning evaluation."""
        ask_mock.side_effect = [False, True]
        self.genome_pair1.evaluations = [self.eval_warning1]
        self.genome_pair2.evaluations = [self.eval_correct1]
        self.genome_pair3.evaluations = [self.eval_warning2]
        self.bndl.genome_pair_dict["gnm_pr1"] = self.genome_pair1
        self.bndl.genome_pair_dict["gnm_pr2"] = self.genome_pair2
        self.bndl.genome_pair_dict["gnm_pr3"] = self.genome_pair3
        import_genome.review_bundled_objects(self.bndl, interactive=True)
        evl_list1 = self.bndl.genome_pair_dict["gnm_pr1"].evaluations
        evl_list2 = self.bndl.genome_pair_dict["gnm_pr2"].evaluations
        evl_list3 = self.bndl.genome_pair_dict["gnm_pr3"].evaluations
        with self.subTest():
            self.assertEqual(evl_list1[0].status, "error")
        with self.subTest():
            self.assertEqual(evl_list2[0].status, "correct")
        with self.subTest():
            self.assertEqual(evl_list3[0].status, "warning")


    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_bundled_objects_9(self, ask_mock):
        """Verify results when bundle has:
        one ticket with four evaluations, three that are errors,
        three CDS features, each with one error evaluation, and
        one source feature with one error,
        but with user exiting review after first ticket error review
        and exiting after first CDS feature error review."""
        self.tkt.evaluations = [self.eval_warning1, self.eval_correct1,
                                self.eval_warning2, self.eval_warning3]
        self.cds1.evaluations = [self.eval_warning4]
        self.cds2.evaluations = [self.eval_warning5]
        self.cds3.evaluations = [self.eval_warning6]
        self.src1.evaluations = [self.eval_warning7]
        self.gnm1.cds_features = [self.cds1, self.cds2, self.cds3]
        self.gnm1.source_features = [self.src1]
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["gnm1"] = self.gnm1
        ask_mock.side_effect = [False, # tkt evl1
                                None, # tkt evl3
                                True, # cds1 evl1
                                None, # cds2 evl1
                                True] # src evl1
        import_genome.review_bundled_objects(self.bndl, interactive=True)
        tkt_list = self.bndl.ticket.evaluations
        cds1_list = self.bndl.genome_dict["gnm1"].cds_features[0].evaluations
        cds2_list = self.bndl.genome_dict["gnm1"].cds_features[1].evaluations
        cds3_list = self.bndl.genome_dict["gnm1"].cds_features[2].evaluations
        src1_list = self.bndl.genome_dict["gnm1"].source_features[0].evaluations
        with self.subTest():
            self.assertEqual(tkt_list[0].status, "error")
        with self.subTest():
            self.assertEqual(tkt_list[1].status, "correct")
        with self.subTest():
            self.assertEqual(tkt_list[2].status, "error")
        with self.subTest():
            self.assertEqual(tkt_list[3].status, "error")
        with self.subTest():
            self.assertEqual(cds1_list[0].status, "warning")
        with self.subTest():
            self.assertEqual(cds2_list[0].status, "error")
        with self.subTest():
            self.assertEqual(cds3_list[0].status, "error")
        with self.subTest():
            self.assertEqual(src1_list[0].status, "warning")




class TestImportGenomeMain8(unittest.TestCase):

    def setUp(self):

        self.cds1 = cds.Cds()
        self.cds1.id = "L5_001"
        self.cds1.start = 10
        self.cds1.stop = 20
        self.cds1.orientation = "F"
        self.cds1.product = "int"
        self.cds1.function = ""
        self.cds1.note = ""

        self.cds2 = cds.Cds()
        self.cds2.id = "L5_002"
        self.cds2.start = 12345
        self.cds2.stop = 12445
        self.cds2.orientation = "F"
        self.cds2.product = "capsid"
        self.cds2.function = ""
        self.cds2.note = "random"

        self.cds3 = cds.Cds()
        self.cds3.id = "L5_001"
        self.cds3.start = 150000
        self.cds3.stop = 151000
        self.cds3.orientation = "R"
        self.cds3.product = "rep"
        self.cds3.function = "lysB"
        self.cds3.note = ""

        self.field = "product"

        self.gnm = genome.Genome()
        self.tkt = ticket.ImportTicket()




    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_cds_descriptions_1(self, ask_mock):
        """Verify no change if there are no CDS features."""
        ask_mock.return_value = True
        features = []
        new_field = import_genome.review_cds_descriptions(features, self.field)
        with self.subTest():
            self.assertFalse(ask_mock.called)
        with self.subTest():
            self.assertEqual(new_field, "product")

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_cds_descriptions_2(self, ask_mock):
        """Verify no change if no descriptions in fields other than
        description_field."""
        ask_mock.return_value = True
        features = [self.cds1]
        new_field = import_genome.review_cds_descriptions(features, self.field)
        with self.subTest():
            self.assertFalse(ask_mock.called)
        with self.subTest():
            self.assertEqual(new_field, "product")

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_cds_descriptions_3(self, ask_mock):
        """Verify no change if there is a description in 'product' field and
        description_field = 'function', and
        description_field is not changed."""
        ask_mock.return_value = True
        features = [self.cds1]
        new_field = import_genome.review_cds_descriptions(features, "function")
        with self.subTest():
            self.assertTrue(ask_mock.called)
        with self.subTest():
            self.assertEqual(new_field, "function")

    @patch("pdm_utils.functions.basic.choose_from_list")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_cds_descriptions_4(self, ask_mock, cfl_mock):
        """Verify no change if there is a description in 'product' field and
        description_field = 'function', and
        description_field is changed."""
        ask_mock.return_value = False
        cfl_mock.return_value = "product"
        features = [self.cds1]
        new_field = import_genome.review_cds_descriptions(features, "function")
        self.assertEqual(new_field, "product")

    @patch("pdm_utils.functions.basic.choose_from_list")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_cds_descriptions_5(self, ask_mock, cfl_mock):
        """Verify no change if there is a description in 'product' field and
        description_field = 'function', and
        description_field was not correctly chosen from list."""
        ask_mock.return_value = False
        cfl_mock.return_value = None
        features = [self.cds1]
        new_field = import_genome.review_cds_descriptions(features, "function")
        self.assertEqual(new_field, "function")

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_cds_descriptions_6(self, ask_mock):
        """Verify no change if there are descriptions in fields in
        multiple CDS features other than description_field."""
        ask_mock.return_value = True
        features = [self.cds1, self.cds2, self.cds3]
        new_field = import_genome.review_cds_descriptions(features, self.field)
        with self.subTest():
            self.assertTrue(ask_mock.called)
        with self.subTest():
            self.assertEqual(new_field, "product")

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_cds_descriptions_7(self, ask_mock):
        """Verify no change if description_field is not correct,
        and no new field is selected."""
        ask_mock.side_effect = [False, False, False]
        features = [self.cds1, self.cds2, self.cds3]
        new_field = import_genome.review_cds_descriptions(features, self.field)
        self.assertEqual(new_field, "product")

    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_review_cds_descriptions_8(self, ask_mock):
        """Verify no change if description_field is not correct,
        and new field is selected."""
        ask_mock.side_effect = [False, False, True]
        features = [self.cds1, self.cds2, self.cds3]
        new_field = import_genome.review_cds_descriptions(features, self.field)
        self.assertTrue(new_field in {"function", "note"})




    def test_set_cds_descriptions_1(self):
        """Verify descriptions are set from product field,
        with no interactivity."""
        self.gnm.cds_features = [self.cds1, self.cds2, self.cds3]
        self.tkt.description_field = "product"
        import_genome.set_cds_descriptions(self.gnm, self.tkt, interactive=False)
        with self.subTest():
            self.assertEqual(self.gnm._cds_descriptions_tally, 3)
        with self.subTest():
            self.assertEqual(self.cds1.description, "int")
        with self.subTest():
            self.assertEqual(self.cds2.description, "capsid")
        with self.subTest():
            self.assertEqual(self.cds3.description, "rep")

    @patch("pdm_utils.functions.basic.choose_from_list")
    @patch("pdm_utils.functions.basic.ask_yes_no")
    def test_set_cds_descriptions_2(self, ask_mock, choose_mock):
        """Verify descriptions are set from product field,
        with interactivity = True and user changes description_field."""
        self.gnm.cds_features = [self.cds1, self.cds2, self.cds3]
        self.tkt.description_field = "product"
        ask_mock.side_effect = [False]
        choose_mock.side_effect = ["function"]
        import_genome.set_cds_descriptions(self.gnm, self.tkt, interactive=True)
        with self.subTest():
            self.assertEqual(self.gnm._cds_descriptions_tally, 1)
        with self.subTest():
            self.assertEqual(self.cds1.description, "")
        with self.subTest():
            self.assertEqual(self.cds2.description, "")
        with self.subTest():
            self.assertEqual(self.cds3.description, "lysB")
        with self.subTest():
            self.assertEqual(self.tkt.description_field, "function")




class TestImportGenomeClass9(unittest.TestCase):
    def setUp(self):
        self.evl1 = eval.Eval()
        self.evl1.id = "GNM0001"
        self.evl1.definition = "temp"
        self.evl1.status = "error"
        self.evl1.result = "Failed evaluation."

        self.evl2 = eval.Eval()
        self.evl2.id = "GNM0002"
        self.evl2.definition = "temp"
        self.evl2.status = "error"
        self.evl2.result = "Failed evaluation."

        self.evl3 = eval.Eval()
        self.evl3.id = "GNM0003"
        self.evl3.definition = "temp"
        self.evl3.status = "correct"
        self.evl3.result = "Failed evaluation."

        self.base_dir = Path(test_root_dir, "test_folder")
        self.base_dir.mkdir()
        self.log_file1 = Path(self.base_dir, "test_log1.txt")
        self.log_file2 = Path(self.base_dir, "test_log2.txt")

    def tearDown(self):
        shutil.rmtree(self.base_dir)


    def test_log_evaluations_1(self):
        """Verify function executes when logfile_path is None."""
        # Nothing really to test.
        evaluation_dict = {1:{"bundle": [self.evl1],
                              "ticket": [self.evl2]},
                           2:{"genome": [self.evl3]}}
        import_genome.log_evaluations(evaluation_dict, logfile_path=None)
        self.assertFalse(self.log_file1.exists())


    def test_log_evaluations_2(self):
        """Verify file-specific log file is created when
        logfile_path is not None."""
        # Nothing really to test.
        evaluation_dict = {1:{"bundle": [self.evl1],
                              "ticket": [self.evl2]},
                           2:{"genome": [self.evl3]}}
        import_genome.log_evaluations(evaluation_dict,
                                      logfile_path=self.log_file1)

        # Confirm how many evaluations were logged.
        with open(self.log_file1,'r') as file:
            lines = file.readlines()

        with self.subTest():
            self.assertTrue(self.log_file1.exists())
        with self.subTest():
            self.assertEqual(len(lines), 2)


    def test_log_evaluations_3(self):
        """Verify file-specific log file is created new each time the
        function is called."""
        # Nothing really to test.
        evaluation_dict1 = {1:{"bundle": [self.evl1],
                              "ticket": [self.evl2]},
                           2:{"genome": [self.evl3]}}
        import_genome.log_evaluations(evaluation_dict1,
                                      logfile_path=self.log_file1)

        evaluation_dict2 = {1:{"bundle": [self.evl1]}}
        import_genome.log_evaluations(evaluation_dict2,
                                      logfile_path=self.log_file2)

        # Confirm how many evaluations were logged.
        with open(self.log_file1,'r') as file:
            lines1 = file.readlines()
        with open(self.log_file2,'r') as file:
            lines2 = file.readlines()

        with self.subTest():
            self.assertTrue(self.log_file1.exists())
        with self.subTest():
            self.assertTrue(self.log_file2.exists())
        with self.subTest():
            self.assertEqual(len(lines1), 2)
        with self.subTest():
            self.assertEqual(len(lines2), 1)










if __name__ == '__main__':
    unittest.main()
