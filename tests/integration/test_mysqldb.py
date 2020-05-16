"""Integration tests for misc. functions that interact with a MySQL database."""

from pathlib import Path
import unittest
import subprocess
import sys
from unittest.mock import patch

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import sqlalchemy

from pdm_utils.classes import cds, genome
from pdm_utils.constants import constants
from pdm_utils.functions import mysqldb

# Import helper functions to build mock database and mock flat files
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils
import test_data_utils

# The following integration tests user the 'pdm_anon' MySQL user.
# It is expected that this user has all privileges for 'pdm_test_db' database.
user = test_db_utils.USER
pwd = test_db_utils.PWD
db = test_db_utils.DB
db2 = "Actinobacteriophage"
engine_string1 = test_db_utils.create_engine_string()


class TestMysqldbFunctions1(unittest.TestCase):

    def setUp(self):
        self.engine = sqlalchemy.create_engine(engine_string1, echo=False)
        test_db_utils.create_empty_test_db()
        phage_data1 = test_data_utils.get_trixie_phage_data()
        phage_data2 = test_data_utils.get_trixie_phage_data()
        phage_data3 = test_data_utils.get_trixie_phage_data()

        phage_data1["PhageID"] = "L5"
        phage_data2["PhageID"] = "Trixie"
        phage_data3["PhageID"] = "D29"

        phage_data1["HostGenus"] = "Mycobacterium"
        phage_data2["HostGenus"] = "Mycobacterium"
        phage_data3["HostGenus"] = "Gordonia"

        phage_data1["Accession"] = "ABC123"
        phage_data2["Accession"] = "XYZ456"
        phage_data3["Accession"] = ""

        phage_data1["Cluster"] = "A"
        phage_data2["Cluster"] = "B"
        phage_data3["Cluster"] = "NULL"

        phage_data1["Subcluster"] = "A1"
        phage_data2["Subcluster"] = "NULL"
        phage_data3["Subcluster"] = "NULL"

        phage_data1["Sequence"] = "atcg"
        phage_data2["Sequence"] = "AATT"
        phage_data3["Sequence"] = "GGCC"

        phage_data1["Length"] = 6
        phage_data2["Length"] = 4
        phage_data3["Length"] = 5

        phage_data1["DateLastModified"] = constants.EMPTY_DATE
        phage_data2["DateLastModified"] = constants.EMPTY_DATE
        phage_data3["DateLastModified"] = constants.EMPTY_DATE

        phage_data_list = [phage_data1, phage_data2, phage_data3]
        for phage_data in phage_data_list:
            test_db_utils.insert_phage_data(phage_data)

        gene_data1 = test_data_utils.get_trixie_gene_data()
        gene_data2 = test_data_utils.get_trixie_gene_data()
        gene_data3 = test_data_utils.get_trixie_gene_data()
        gene_data4 = test_data_utils.get_trixie_gene_data()

        gene_data1["PhageID"] = "Trixie"
        gene_data2["PhageID"] = "Trixie"
        gene_data3["PhageID"] = "Trixie"
        gene_data4["PhageID"] = "D29"

        gene_data1["GeneID"] = "Trixie_1"
        gene_data2["GeneID"] = "Trixie_2"
        gene_data3["GeneID"] = "Trixie_3"
        gene_data4["GeneID"] = "D29_1"

        gene_data_list = [gene_data1, gene_data2, gene_data3, gene_data4]
        for gene_data in gene_data_list:
            test_db_utils.insert_gene_data(gene_data)

        self.phage_query = "SELECT * FROM phage"
        self.gene_query = "SELECT * FROM gene"


    def tearDown(self):
        self.engine.dispose()
        test_db_utils.remove_db()




    # TODO not sure if this is needed.
    def test_verify_db_setup(self):
        """Confirm that the database was setup correctly for the tests."""
        phage_data = test_db_utils.get_data(test_db_utils.phage_table_query)
        gene_data = test_db_utils.get_data(test_db_utils.gene_table_query)
        with self.subTest():
            self.assertEqual(len(phage_data), 3)
        with self.subTest():
            self.assertEqual(len(gene_data), 4)


    def test_create_seq_set_1(self):
        """Retrieve a set of all data from Sequence column."""
        result = mysqldb.create_seq_set(self.engine)
        with self.subTest():
            self.assertEqual(len(result), 3)
        with self.subTest():
            self.assertTrue(Seq("ATCG", IUPAC.ambiguous_dna) in result)


    def test_parse_genome_data_1(self):
        """Verify that a Genome object is constructed correctly for a
        valid PhageID."""
        genome_list = mysqldb.parse_genome_data(self.engine,
                        phage_id_list=["L5"], phage_query=self.phage_query,
                        gnm_type="mysql")
        with self.subTest():
            self.assertEqual(len(genome_list), 1)
        with self.subTest():
            self.assertEqual(genome_list[0].id, "L5")
        with self.subTest():
            self.assertEqual(genome_list[0].seq, "ATCG")
        with self.subTest():
            self.assertEqual(genome_list[0].type, "mysql")
        with self.subTest():
            self.assertEqual(genome_list[0].date, constants.EMPTY_DATE)
        with self.subTest():
            self.assertEqual(len(genome_list[0].cds_features), 0)

    def test_parse_genome_data_2(self):
        """Verify that an empty Genome object list is constructed for an
        invalid PhageID."""
        genome_list = mysqldb.parse_genome_data(
                          self.engine, phage_id_list=["EagleEye"],
                          phage_query=self.phage_query)
        self.assertEqual(len(genome_list), 0)

    def test_parse_genome_data_3(self):
        """Verify that a Genome object with CDS features
        is constructed correctly for a valid PhageID."""
        genome_list = mysqldb.parse_genome_data(
                        self.engine, phage_id_list=["Trixie"],
                        phage_query=self.phage_query,
                        gene_query=self.gene_query)
        with self.subTest():
            self.assertEqual(len(genome_list), 1)
        with self.subTest():
            self.assertEqual(genome_list[0].id, "Trixie")
        with self.subTest():
            self.assertEqual(genome_list[0].seq, "AATT")
        with self.subTest():
            self.assertEqual(genome_list[0].type, "")
        with self.subTest():
            self.assertEqual(genome_list[0].date, constants.EMPTY_DATE)
        with self.subTest():
            self.assertEqual(len(genome_list[0].cds_features), 3)
        with self.subTest():
            self.assertEqual(genome_list[0].cds_features[0].genome_length, 4)

    def test_parse_genome_data_4(self):
        """Verify that multiple Genome objects with CDS features
        are constructed correctly for multiple valid PhageIDs."""
        genome_list = mysqldb.parse_genome_data(
                        self.engine, phage_query=self.phage_query,
                        gene_query=self.gene_query)

        genome_dict = {}
        for gnm in genome_list:
            genome_dict[gnm.id] = gnm

        with self.subTest():
            self.assertEqual(len(genome_list), 3)
        with self.subTest():
            self.assertEqual(genome_dict["Trixie"].seq, "AATT")
        with self.subTest():
            self.assertEqual(len(genome_dict["Trixie"].cds_features), 3)
        with self.subTest():
            self.assertEqual(
                genome_dict["Trixie"].cds_features[0].genome_length, 4)
        with self.subTest():
            self.assertEqual(
                genome_dict["Trixie"].cds_features[1].genome_length, 4)
        with self.subTest():
            self.assertEqual(len(genome_dict["D29"].cds_features), 1)
        with self.subTest():
            self.assertEqual(
                genome_dict["D29"].cds_features[0].genome_length, 5)
        with self.subTest():
            self.assertEqual(len(genome_dict["L5"].cds_features), 0)




    def test_parse_cds_data_1(self):
        """Verify that a Cds object is constructed correctly for a
        valid PhageID."""
        cds_list = mysqldb.parse_cds_data(self.engine, column="PhageID",
                                             phage_id_list=["Trixie"],
                                             query=self.gene_query)
        with self.subTest():
            self.assertEqual(len(cds_list), 3)
        with self.subTest():
            self.assertEqual(cds_list[0].id, "Trixie_1")
        with self.subTest():
            self.assertEqual(cds_list[0].genome_id, "Trixie")

    def test_parse_cds_data_2(self):
        """Verify that an empty Cds object list is constructed for an
        invalid PhageID."""
        cds_list = mysqldb.parse_cds_data(self.engine, column="PhageID",
                                          phage_id_list=["L5"],
                                          query=self.gene_query)
        self.assertEqual(len(cds_list), 0)

    def test_parse_cds_data_3(self):
        """Verify that Cds objects are constructed correctly for multiple
        valid PhageID when the phage_id parameter is not specified."""
        cds_list = mysqldb.parse_cds_data(self.engine, query=self.gene_query)
        with self.subTest():
            self.assertEqual(len(cds_list), 4)




class TestMysqldbFunctions2(unittest.TestCase):

    def setUp(self):
        self.engine = sqlalchemy.create_engine(engine_string1, echo=False)
        test_db_utils.create_empty_test_db()
        test_db_utils.execute("TRUNCATE version")

    def tearDown(self):
        self.engine.dispose()
        test_db_utils.remove_db()




    def test_create_phage_table_insert_1(self):
        """Verify phage table INSERT statement is created correctly."""
        # Note: even though this function returns a string and doesn't
        # actually utilize a MySQL database, this test ensures
        # that the returned statement will function properly in MySQL.
        gnm = genome.Genome()
        gnm.id = "L5"
        gnm.name = "L5_Draft"
        gnm.host_genus = "Mycobacterium"
        gnm.annotation_status = "final"
        gnm.accession = "ABC123"
        gnm.seq = Seq("ATCG", IUPAC.ambiguous_dna)
        gnm.length = 4
        gnm.gc = 0.5001
        gnm.date = constants.EMPTY_DATE
        gnm.retrieve_record = 1
        gnm.annotation_author = 1
        gnm.cluster = "Singleton"
        gnm.subcluster = "A2"
        statement = mysqldb.create_phage_table_insert(gnm)
        test_db_utils.execute(statement)
        phage_data = test_db_utils.get_data(test_db_utils.phage_table_query)
        results = phage_data[0]
        exp = ("INSERT INTO phage "
               "(PhageID, Accession, Name, HostGenus, Sequence, "
               "Length, GC, Status, DateLastModified, RetrieveRecord, "
               "AnnotationAuthor, Cluster, Subcluster) "
               "VALUES "
               "('L5', 'ABC123', 'L5_Draft', 'Mycobacterium', 'ATCG', "
               f"4, 0.5001, 'final', '{constants.EMPTY_DATE}', 1, "
               "1, NULL, 'A2');")
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(results["PhageID"], "L5")
        with self.subTest():
            self.assertEqual(results["Accession"], "ABC123")
        with self.subTest():
            self.assertEqual(results["Name"], "L5_Draft")
        with self.subTest():
            self.assertEqual(results["HostGenus"], "Mycobacterium")
        with self.subTest():
            self.assertEqual(results["Sequence"].decode("utf-8"), "ATCG")
        with self.subTest():
            self.assertEqual(results["Length"], 4)
        with self.subTest():
            self.assertEqual(results["GC"], 0.5001)
        with self.subTest():
            self.assertEqual(results["Status"], "final")
        with self.subTest():
            self.assertEqual(results["DateLastModified"], constants.EMPTY_DATE)
        with self.subTest():
            self.assertEqual(results["RetrieveRecord"], 1)
        with self.subTest():
            self.assertEqual(results["AnnotationAuthor"], 1)
        with self.subTest():
            self.assertIsNone(results["Cluster"])
        with self.subTest():
            self.assertEqual(results["Subcluster"], "A2")




    def test_change_version_1(self):
        """Verify the version is incremented by 1."""
        data = {"Version": 10, "SchemaVersion": 1}
        test_db_utils.insert_version_data(data)
        mysqldb.change_version(self.engine)
        result = test_db_utils.get_data(test_db_utils.version_table_query)
        output_value = result[0]["Version"]
        self.assertEqual(output_value, 11)

    def test_change_version_2(self):
        """Verify the version is incremented by 5."""
        data = {"Version": 10, "SchemaVersion": 1}
        test_db_utils.insert_version_data(data)
        mysqldb.change_version(self.engine, amount=5)
        result = test_db_utils.get_data(test_db_utils.version_table_query)
        output_value = result[0]["Version"]
        self.assertEqual(output_value, 15)

    def test_change_version_3(self):
        """Verify the version is decremented by 5."""
        data = {"Version": 10, "SchemaVersion": 1}
        test_db_utils.insert_version_data(data)
        mysqldb.change_version(self.engine, amount=-5)
        result = test_db_utils.get_data(test_db_utils.version_table_query)
        output_value = result[0]["Version"]
        self.assertEqual(output_value, 5)




class TestMysqldbFunctions3(unittest.TestCase):

    def setUp(self):
        self.engine = sqlalchemy.create_engine(engine_string1, echo=False)
        test_db_utils.create_empty_test_db()
        phage_data = test_data_utils.get_trixie_phage_data()
        phage_data["PhageID"] = "Trixie"
        phage_data["HostGenus"] = "Mycobacterium"
        phage_data["Accession"] = "ABC123"
        phage_data["Cluster"] = "A"
        phage_data["Subcluster"] = "A1"
        phage_data["Sequence"] = "atcg"
        phage_data["Length"] = 6
        phage_data["DateLastModified"] = constants.EMPTY_DATE
        test_db_utils.insert_phage_data(phage_data)
        self.phage_query = "SELECT * FROM phage WHERE PhageID = 'Trixie'"
        self.gene_query = "SELECT * FROM gene WHERE PhageID = 'Trixie'"


    def tearDown(self):
        self.engine.dispose()
        test_db_utils.remove_db()




    def test_create_gene_table_insert_1(self):
        """Verify gene table INSERT statement is created correctly when
        locus_tag is not empty and description contains a "'"."""
        # Note: even though this function returns a string and doesn't
        # actually utilize a MySQL database, this test ensures
        # that the returned statement will function properly in MySQL.
        cds1 = cds.Cds()
        cds1.id = "SEA_TRIXIE_123"
        cds1.genome_id = "Trixie"
        cds1.start = 5
        cds1.stop = 10
        cds1.parts = 1
        cds1.translation_length = 20
        cds1.length = 200
        cds1.name = "Int"
        cds1.type = "CDS"
        cds1.translation = Seq("ACKLG", IUPAC.protein)
        cds1.orientation = "F"
        cds1.description = "5' nucleotide phosphatase"
        cds1.locus_tag = "TAG1"
        statement = mysqldb.create_gene_table_insert(cds1)
        test_db_utils.execute(statement)
        result = test_db_utils.get_data(self.gene_query)
        results = result[0]
        exp = ("""INSERT INTO gene """
               """(GeneID, PhageID, Start, Stop, Length, Name, """
               """Translation, Orientation, Notes, LocusTag, Parts) """
               """VALUES """
               """("SEA_TRIXIE_123", "Trixie", 5, 10, 200, "Int", """
               """"ACKLG", "F", "5' nucleotide phosphatase", "TAG1", 1);""")
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(results["GeneID"], "SEA_TRIXIE_123")
        with self.subTest():
            self.assertEqual(results["PhageID"], "Trixie")
        with self.subTest():
            self.assertEqual(results["Start"], 5)
        with self.subTest():
            self.assertEqual(results["Stop"], 10)
        with self.subTest():
            self.assertEqual(results["Parts"], 1)
        with self.subTest():
            self.assertEqual(results["Length"], 200)
        with self.subTest():
            self.assertEqual(results["Name"], "Int")
        with self.subTest():
            self.assertEqual(results["Translation"].decode("utf-8"), "ACKLG")
        with self.subTest():
            self.assertEqual(results["Orientation"], "F")
        with self.subTest():
            self.assertEqual(results["Notes"].decode("utf-8"),
                             "5' nucleotide phosphatase")
        with self.subTest():
            self.assertEqual(results["LocusTag"], "TAG1")

    def test_create_gene_table_insert_2(self):
        """Verify gene table INSERT statement is created correctly when
        locus_tag is empty."""
        # Note: even though this function returns a string and doesn't
        # actually utilize a MySQL database, this test ensures
        # that the returned statement will function properly in MySQL.
        cds1 = cds.Cds()
        cds1.id = "SEA_TRIXIE_123"
        cds1.genome_id = "Trixie"
        cds1.start = 5
        cds1.stop = 10
        cds1.parts = 1
        cds1.translation_length = 20
        cds1.length = 200
        cds1.name = "Int"
        cds1.type = "CDS"
        cds1.translation = Seq("ACKLG", IUPAC.protein)
        cds1.orientation = "F"
        cds1.description = "integrase"
        cds1.locus_tag = ""
        statement = mysqldb.create_gene_table_insert(cds1)
        test_db_utils.execute(statement)
        result = test_db_utils.get_data(self.gene_query)
        results = result[0]
        exp = ("""INSERT INTO gene """
               """(GeneID, PhageID, Start, Stop, Length, Name, """
               """Translation, Orientation, Notes, LocusTag, Parts) """
               """VALUES """
               """("SEA_TRIXIE_123", "Trixie", 5, 10, 200, "Int", """
               """"ACKLG", "F", "integrase", NULL, 1);""")
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(results["GeneID"], "SEA_TRIXIE_123")
        with self.subTest():
            self.assertEqual(results["LocusTag"], None)




    def test_create_update_1(self):
        """Verify correct Cluster statement is created for a non-singleton."""
        statement = mysqldb.create_update(
            "phage", "Cluster", "B", "PhageID", "Trixie")
        test_db_utils.execute(statement)
        result = test_db_utils.get_data(self.phage_query)
        results = result[0]
        exp = "UPDATE phage SET Cluster = 'B' WHERE PhageID = 'Trixie';"
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(results["Cluster"], "B")

    def test_create_update_2(self):
        """Verify correct Cluster statement is created for a singleton."""
        statement = mysqldb.create_update(
            "phage", "Cluster", "Singleton", "PhageID", "Trixie")
        test_db_utils.execute(statement)
        result = test_db_utils.get_data(self.phage_query)
        results = result[0]
        exp = "UPDATE phage SET Cluster = NULL WHERE PhageID = 'Trixie';"
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertIsNone(results["Cluster"])

    def test_create_update_3(self):
        """Verify correct Subcluster statement is created for a
        non-empty value."""
        statement = mysqldb.create_update(
            "phage", "Subcluster", "A2", "PhageID", "Trixie")
        test_db_utils.execute(statement)
        result = test_db_utils.get_data(self.phage_query)
        results = result[0]
        exp = "UPDATE phage SET Subcluster = 'A2' WHERE PhageID = 'Trixie';"
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(results["Subcluster"], "A2")

    def test_create_update_4(self):
        """Verify Gene table statement is created correctly."""
        # First add gene Trixie_1 to the database.
        gene_data = test_data_utils.get_trixie_gene_data()
        gene_data["PhageID"] = "Trixie"
        gene_data["Notes"] = "none"
        gene_data["GeneID"] = "Trixie_1"
        test_db_utils.insert_gene_data(gene_data)

        # Second run the update statement.
        statement = mysqldb.create_update(
            "gene", "Notes", "Repressor", "GeneID", "Trixie_1")
        test_db_utils.execute(statement)
        result = test_db_utils.get_data(self.gene_query)
        results = result[0]
        exp = "UPDATE gene SET Notes = 'Repressor' WHERE GeneID = 'Trixie_1';"
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(results["Notes"].decode("utf-8"), "Repressor")




class TestMysqldbFunctions4(unittest.TestCase):

    def setUp(self):
        self.engine = sqlalchemy.create_engine(engine_string1, echo=False)
        test_db_utils.create_empty_test_db()

        phage_data1 = test_data_utils.get_trixie_phage_data()
        phage_data1["PhageID"] = "Trixie"
        phage_data1["HostGenus"] = "Mycobacterium"
        phage_data1["Accession"] = "ABC123"
        phage_data1["Cluster"] = "A"
        phage_data1["Subcluster"] = "A1"
        phage_data1["Sequence"] = "atcg"
        phage_data1["Length"] = 6
        phage_data1["DateLastModified"] = constants.EMPTY_DATE
        test_db_utils.insert_phage_data(phage_data1)

        gene_data1 = test_data_utils.get_trixie_gene_data()
        gene_data2 = test_data_utils.get_trixie_gene_data()

        gene_data1["PhageID"] = "Trixie"
        gene_data2["PhageID"] = "Trixie"

        gene_data1["GeneID"] = "Trixie_1"
        gene_data2["GeneID"] = "Trixie_2"

        test_db_utils.insert_gene_data(gene_data1)
        test_db_utils.insert_gene_data(gene_data2)

        self.phage_query = "SELECT * FROM phage where PhageID = 'Trixie'"
        self.gene_query = "SELECT * FROM gene where PhageID = 'Trixie'"

    def tearDown(self):
        self.engine.dispose()
        test_db_utils.remove_db()




    def test_create_delete_1(self):
        """Verify correct DELETE statement is created
        for a PhageID in the phage table."""
        # First, retrieve the current state of the phage table.
        results1 = test_db_utils.get_data(self.phage_query)
        results1_phageids = set()
        for dict in results1:
            results1_phageids.add(dict["PhageID"])

        # Second, execute the DELETE statement.
        statement = mysqldb.create_delete("phage", "PhageID", "Trixie")
        test_db_utils.execute(statement)

        # Third, retrieve the current state of the phage table.
        results2 = test_db_utils.get_data(self.phage_query)
        results2_phageids = set()
        for dict in results2:
            results2_phageids.add(dict["PhageID"])

        exp = "DELETE FROM phage WHERE PhageID = 'Trixie';"
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(len(results1_phageids), 1)
        with self.subTest():
            self.assertEqual(len(results2_phageids), 0)

    def test_create_delete_2(self):
        """Verify correct DELETE statement is created
        for a single GeneID in the gene table."""

        # First, retrieve the current state of the gene table.
        results1 = test_db_utils.get_data(self.gene_query)
        results1_phageids = set()
        results1_geneids = set()
        for dict in results1:
            results1_phageids.add(dict["PhageID"])
            results1_geneids.add(dict["GeneID"])

        # Second, execute the DELETE statement.
        statement = mysqldb.create_delete("gene", "GeneID", "Trixie_1")
        test_db_utils.execute(statement)

        # Third, retrieve the current state of the gene table.
        results2 = test_db_utils.get_data(self.gene_query)
        results2_phageids = set()
        results2_geneids = set()
        for dict in results2:
            results2_phageids.add(dict["PhageID"])
            results2_geneids.add(dict["GeneID"])

        exp = "DELETE FROM gene WHERE GeneID = 'Trixie_1';"
        with self.subTest():
            self.assertEqual(statement, exp)
        with self.subTest():
            self.assertEqual(len(results1_phageids), 1)
        with self.subTest():
            self.assertEqual(len(results1_geneids), 2)
        with self.subTest():
            self.assertEqual(len(results2_phageids), 1)
        with self.subTest():
            self.assertEqual(len(results2_geneids), 1)




    def test_execute_transaction_1(self):
        """Valid everything should result in creation of cursor and execution
        of all statements in the transaction - return code 0."""
        valid1 = ("INSERT INTO phage "
                   "(PhageID, Accession, Name, HostGenus, Sequence, "
                   "Length, GC, Status, DateLastModified, "
                   "RetrieveRecord, AnnotationAuthor,"
                   "Cluster, Subcluster) "
                   "VALUES ('D29', 'ABC123', 'L5_Draft', 'Mycobacterium', "
                   "'ATCG', 4, 0.5001, 'final', "
                   f"'{constants.EMPTY_DATE}', 1, 1, 'A', 'A2');")
        valid2 = ("INSERT INTO phage "
                   "(PhageID, Accession, Name, HostGenus, Sequence, "
                   "Length, GC, Status, DateLastModified, "
                   "RetrieveRecord, AnnotationAuthor,"
                   "Cluster, Subcluster) "
                   "VALUES ('L5', 'ABC123', 'L5_Draft', 'Mycobacterium', "
                   "'ATCG', 4, 0.5001, 'final', "
                   f"'{constants.EMPTY_DATE}', 1, 1, 'A', 'A2');")
        valid_stmts = [valid1, valid2]
        return_code, msg = mysqldb.execute_transaction(self.engine, valid_stmts)
        query = "SELECT COUNT(PhageID) FROM phage"
        result_list = self.engine.execute(query).fetchall()
        count = result_list[0][0]
        with self.subTest():
            self.assertEqual(count, 3)
        with self.subTest():
            self.assertEqual(return_code, 0)

    def test_execute_transaction_2(self):
        """Valid connection but invalid transaction should return code 1."""
        valid1 = ("INSERT INTO phage "
                   "(PhageID, Accession, Name, HostGenus, Sequence, "
                   "Length, GC, Status, DateLastModified, "
                   "RetrieveRecord, AnnotationAuthor,"
                   "Cluster, Subcluster) "
                   "VALUES ('D29', 'ABC123', 'L5_Draft', 'Mycobacterium', "
                   "'ATCG', 4, 0.5001, 'final', "
                   f"'{constants.EMPTY_DATE}', 1, 1, 'A', 'A2');")
        invalid1 = ("INSERT INTO phage "
                   "(PhageID, Accession, Name, HostGenus, Sequence, "
                   "Length, GC, Status, DateLastModified, "
                   "RetrieveRecord, AnnotationAuthor,"
                   "Cluster, Subcluster) "
                   "VALUES ('Trixie', 'ABC123', 'L5_Draft', 'Mycobacterium', "
                   "'ATCG', 4, 0.5001, 'final', "
                   f"'{constants.EMPTY_DATE}', 1, 1, 'A', 'A2');")
        invalid_stmts = [valid1, invalid1]
        return_code, msg = mysqldb.execute_transaction(self.engine, invalid_stmts)
        query = "SELECT COUNT(PhageID) FROM phage"
        result_list = self.engine.execute(query).fetchall()
        count = result_list[0][0]
        with self.subTest():
            self.assertEqual(count, 1)
        with self.subTest():
            self.assertEqual(return_code, 1)

    def test_execute_transaction_3(self):
        """Everything ok but no transaction should return 0."""
        return_code, msg = mysqldb.execute_transaction(self.engine)
        self.assertEqual(return_code, 0)




if __name__ == '__main__':
    unittest.main()
