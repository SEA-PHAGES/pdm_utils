""" Unit tests for misc. functions that interact with a MySQL database."""


import unittest
from pdm_utils.functions import mysqldb
from pdm_utils.classes import genome
from pdm_utils.classes import cds
from datetime import datetime
from pdm_utils.classes import bundle
from pdm_utils.constants import constants
from Bio.Seq import Seq



class TestMysqldbFunctions1(unittest.TestCase):


    def setUp(self):
        self.genome1 = genome.Genome()
        self.genome2 = genome.Genome()
        self.genome3 = genome.Genome()




    def test_parse_phage_table_data_1(self):
        """Verify standard MySQL genome data is parsed correctly
        from a data dictionary returned from a SQL query."""

        data_dict = {"PhageID":"L5",
                     "Accession":"ABC123",
                     "Name":"L5_Draft",
                     "HostGenus":"Mycobacterium",
                     "Sequence":"ATCG".encode("utf-8"),
                     "Length":10,
                     "DateLastModified":constants.EMPTY_DATE,
                     "Notes":"abc".encode("utf-8"),
                     "GC":12.12,
                     "Cluster":"A",
                     "Subcluster":"A2",
                     "Status":"final",
                     "RetrieveRecord":1,
                     "AnnotationAuthor":1}

        self.genome1 = \
            mysqldb.parse_phage_table_data(data_dict, gnm_type="mysql")

        with self.subTest():
            self.assertEqual(self.genome1.id, "L5")
        with self.subTest():
            self.assertEqual(self.genome1.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome1.name, "L5_Draft")
        with self.subTest():
            self.assertEqual(self.genome1.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome1.seq, "ATCG")
        with self.subTest():
            self.assertIsInstance(self.genome1.seq, Seq)
        with self.subTest():
            self.assertEqual(self.genome1.length, 10)
        with self.subTest():
            self.assertEqual(self.genome1.date, constants.EMPTY_DATE)
        with self.subTest():
            self.assertEqual(self.genome1.description, "abc")
        with self.subTest():
            self.assertEqual(self.genome1.gc, 12.12)
        with self.subTest():
            self.assertEqual(self.genome1.cluster, "A")
        with self.subTest():
            self.assertEqual(self.genome1.subcluster, "A2")
        with self.subTest():
            self.assertEqual(self.genome1.annotation_status, "final")
        with self.subTest():
            self.assertEqual(self.genome1.retrieve_record, 1)
        with self.subTest():
            self.assertEqual(self.genome1.annotation_author, 1)
        with self.subTest():
            self.assertEqual(self.genome1.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome1.type, "mysql")


    def test_parse_phage_table_data_2(self):
        """Verify truncated MySQL genome data is parsed correctly
        from a data dictionary returned from a SQL query."""

        data_dict = {"PhageID":"L5"}
        self.genome1 = mysqldb.parse_phage_table_data(data_dict)
        with self.subTest():
            self.assertEqual(self.genome1.id, "L5")
        with self.subTest():
            self.assertEqual(self.genome1.name, "")
        with self.subTest():
            self.assertEqual(self.genome1.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome1.type, "")




    def test_parse_gene_table_data_1(self):
        """Verify standard MySQL CDS data is parsed correctly
        from a data dictionary returned from a SQL query."""

        data_dict = {"GeneID":"L5_001",
                     "PhageID":"L5",
                     "Start":10,
                     "Stop":100,
                     "Parts":1,
                     "Length":1000,
                     "Name":"1",
                     "Translation":"AGGPT",
                     "Orientation":"F",
                     "Notes":"description".encode("utf-8"),
                     "DomainStatus":1,
                     "PhamID":1,
                     "LocusTag":"SEA_L5_001"}

        cds1 = mysqldb.parse_gene_table_data(data_dict)

        with self.subTest():
            self.assertEqual(cds1.id, "L5_001")
        with self.subTest():
            self.assertEqual(cds1.genome_id, "L5")
        with self.subTest():
            self.assertEqual(cds1.start, 10)
        with self.subTest():
            self.assertEqual(cds1.stop, 100)
        with self.subTest():
            self.assertEqual(cds1.parts, 1)
        with self.subTest():
            self.assertEqual(cds1.length, 1000)
        with self.subTest():
            self.assertEqual(cds1.name, "1")
        with self.subTest():
            self.assertEqual(cds1.type, "CDS")
        with self.subTest():
            self.assertEqual(cds1.translation, "AGGPT")
        with self.subTest():
            self.assertEqual(cds1.orientation, "F")
        with self.subTest():
            self.assertEqual(cds1.description, "description")
        with self.subTest():
            self.assertEqual(cds1.locus_tag, "SEA_L5_001")
        with self.subTest():
            self.assertEqual(cds1.translation_table, 11)
        with self.subTest():
            self.assertEqual(cds1.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(cds1.domain_status, 1)
        with self.subTest():
            self.assertEqual(cds1.pham_id, 1)


    def test_parse_gene_table_data_2(self):
        """Verify truncated MySQL CDS data is parsed correctly
        from a data dictionary returned from a SQL query."""

        data_dict = {"GeneID":"L5_001",
                     "PhageID":"L5",
                     "Start":10,
                     "LocusTag":None}
        cds1 = mysqldb.parse_gene_table_data(data_dict)
        with self.subTest():
            self.assertEqual(cds1.id, "L5_001")
        with self.subTest():
            self.assertEqual(cds1.genome_id, "L5")
        with self.subTest():
            self.assertEqual(cds1.start, 10)
        with self.subTest():
            self.assertEqual(cds1.locus_tag, "")
        with self.subTest():
            self.assertEqual(cds1.translation_table, 11)




    def test_create_genome_statements_1(self):
        """Verify list of INSERT statements is created correctly for:
        'add' ticket, and no CDS features."""

        self.genome1.id = "L5"
        self.genome1.name = "L5_Draft"
        self.genome1.host_genus = "Mycobacterium"
        self.genome1.annotation_status = "final"
        self.genome1.accession = "ABC123"
        self.genome1.seq = "ATCG"
        self.genome1.length = 4
        self.genome1.gc = 0.5001
        self.genome1.date = '1/1/2000'
        self.genome1.retrieve_record = "1"
        self.genome1.annotation_author = "1"
        self.genome1.cluster = "A"
        self.genome1.subcluster = "A2"

        statements = mysqldb.create_genome_statements(
                        self.genome1, tkt_type="add")
        self.assertEqual(len(statements), 1)


    def test_create_genome_statements_2(self):
        """Verify list of INSERT statements is created correctly for:
        'replace' ticket, and no CDS features."""

        self.genome1.id = "L5"
        self.genome1.name = "L5_Draft"
        self.genome1.host_genus = "Mycobacterium"
        self.genome1.annotation_status = "final"
        self.genome1.accession = "ABC123"
        self.genome1.seq = "ATCG"
        self.genome1.length = 4
        self.genome1.gc = 0.5001
        self.genome1.date = '1/1/2000'
        self.genome1.retrieve_record = "1"
        self.genome1.annotation_author = "1"
        self.genome1.cluster = "A"
        self.genome1.subcluster = "A2"

        statements = mysqldb.create_genome_statements(
                        self.genome1, tkt_type="replace")
        self.assertEqual(len(statements), 2)


    def test_create_genome_statements_3(self):
        """Verify list of INSERT statements is created correctly for:
        'add' ticket, and two CDS features."""

        cds1 = cds.Cds()
        cds1.genome_id = "L5"
        cds1.start = 10
        cds1.stop = 100
        cds1.parts = 1
        cds1.length = 1000
        cds1.name = "1"
        cds1.type = "CDS"
        cds1.translation = "AGGPT"
        cds1.orientation = "F"
        cds1.description = "description"
        cds1.locus_tag = "SEA_L5_001"

        cds2 = cds.Cds()
        cds2.genome_id = "L5"
        cds2.start = 100
        cds2.stop = 1000
        cds2.parts = 1
        cds2.length = 10000
        cds2.name = "2"
        cds2.type = "CDS"
        cds2.translation = "AKKQE"
        cds2.orientation = "R"
        cds2.description = "description"
        cds2.locus_tag = "SEA_L5_002"

        self.genome1.id = "L5"
        self.genome1.name = "L5_Draft"
        self.genome1.host_genus = "Mycobacterium"
        self.genome1.annotation_status = "final"
        self.genome1.accession = "ABC123"
        self.genome1.seq = "ATCG"
        self.genome1.length = 4
        self.genome1.gc = 0.5001
        self.genome1.date = '1/1/2000'
        self.genome1.retrieve_record = "1"
        self.genome1.annotation_author = "1"
        self.genome1.cluster = "A"
        self.genome1.subcluster = "A2"
        self.genome1.cds_features = [cds1, cds2]

        statements = mysqldb.create_genome_statements(
                        self.genome1, tkt_type="add")
        self.assertEqual(len(statements), 3)


if __name__ == '__main__':
    unittest.main()
