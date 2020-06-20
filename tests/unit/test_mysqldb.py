""" Unit tests for misc. functions that interact with a MySQL database."""


from datetime import datetime
from pathlib import Path
import sys
import unittest
from unittest.mock import Mock, patch, PropertyMock

from Bio.Seq import Seq

from pdm_utils.classes import bundle
from pdm_utils.classes import cds, trna, tmrna
from pdm_utils.classes import genome
from pdm_utils.constants import constants
from pdm_utils.functions import mysqldb


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
        data_dict = {"GeneID": "L5_001",
                     "PhageID": "L5",
                     "Start": 10,
                     "Stop": 100,
                     "Parts": 1,
                     "Length": 1000,
                     "Name": "1",
                     "Translation": "AGGPT".encode("utf-8"),
                     "Orientation": "F",
                     "Notes": "description".encode("utf-8"),
                     "DomainStatus": 1,
                     "PhamID": 1,
                     "LocusTag": "SEA_L5_001"}

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




    def test_parse_trna_table_data_1(self):
        """Verify standard MySQL tRNA data is parsed correctly
        from a data dictionary returned from a SQL query."""
        data_dict = {
            "GeneID": "TRIXIE_0001",
            "PhageID": "Trixie",
            "Start": 100,
            "Stop": 1100,
            "Length": 1000,
            "Name": "1",
            "Orientation": "F",
            "Note": "misc".encode("utf-8"),
            "LocusTag": "SEA_TRIXIE_0001",
            "AminoAcid": "Ala",
            "Anticodon": "AAA",
            "Structure": "AAAAAAAA".encode("utf-8"),
            "Source": "aragorn"
            }

        trna1 = mysqldb.parse_trna_table_data(data_dict)

        with self.subTest():
            self.assertEqual(trna1.id, "TRIXIE_0001")
        with self.subTest():
            self.assertEqual(trna1.genome_id, "Trixie")
        with self.subTest():
            self.assertEqual(trna1.start, 100)
        with self.subTest():
            self.assertEqual(trna1.stop, 1100)
        with self.subTest():
            self.assertEqual(trna1.parts, 1)
        with self.subTest():
            self.assertEqual(trna1.length, 1000)
        with self.subTest():
            self.assertEqual(trna1.name, "1")
        with self.subTest():
            self.assertEqual(trna1.type, "tRNA")
        with self.subTest():
            self.assertEqual(trna1.orientation, "F")
        with self.subTest():
            self.assertEqual(trna1.locus_tag, "SEA_TRIXIE_0001")
        with self.subTest():
            self.assertEqual(trna1.note, "misc")
        with self.subTest():
            self.assertEqual(trna1.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(trna1.amino_acid, "Ala")
        with self.subTest():
            self.assertEqual(trna1.anticodon, "AAA")
        with self.subTest():
            self.assertEqual(trna1.structure, "AAAAAAAA")
        with self.subTest():
            self.assertEqual(trna1.use, "aragorn")

    def test_parse_trna_table_data_2(self):
        """Verify truncated MySQL tRNA data is parsed correctly
        from a data dictionary returned from a SQL query."""
        data_dict = {
            "GeneID": "TRIXIE_0001",
            "PhageID": "Trixie",
            "Start": 100,
            "LocusTag": "SEA_TRIXIE_0001"
            }

        trna1 = mysqldb.parse_trna_table_data(data_dict)

        with self.subTest():
            self.assertEqual(trna1.id, "TRIXIE_0001")
        with self.subTest():
            self.assertEqual(trna1.genome_id, "Trixie")
        with self.subTest():
            self.assertEqual(trna1.start, 100)
        with self.subTest():
            self.assertEqual(trna1.locus_tag, "SEA_TRIXIE_0001")




    def test_parse_tmrna_table_data_1(self):
        """Verify standard MySQL tmRNA data is parsed correctly
        from a data dictionary returned from a SQL query."""
        data_dict = {
            "GeneID": "TRIXIE_0001",
            "PhageID": "Trixie",
            "Start": 100,
            "Stop": 1100,
            "Length": 1000,
            "Name": "1",
            "Orientation": "F",
            "Note": "misc".encode("utf-8"),
            "LocusTag": "SEA_TRIXIE_0001",
            "PeptideTag": "random"
            }

        tmrna1 = mysqldb.parse_tmrna_table_data(data_dict)

        with self.subTest():
            self.assertEqual(tmrna1.id, "TRIXIE_0001")
        with self.subTest():
            self.assertEqual(tmrna1.genome_id, "Trixie")
        with self.subTest():
            self.assertEqual(tmrna1.start, 100)
        with self.subTest():
            self.assertEqual(tmrna1.stop, 1100)
        with self.subTest():
            self.assertEqual(tmrna1.parts, 1)
        with self.subTest():
            self.assertEqual(tmrna1.length, 1000)
        with self.subTest():
            self.assertEqual(tmrna1.name, "1")
        with self.subTest():
            self.assertEqual(tmrna1.type, "tmRNA")
        with self.subTest():
            self.assertEqual(tmrna1.orientation, "F")
        with self.subTest():
            self.assertEqual(tmrna1.locus_tag, "SEA_TRIXIE_0001")
        with self.subTest():
            self.assertEqual(tmrna1.note, "misc")
        with self.subTest():
            self.assertEqual(tmrna1.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(tmrna1.peptide_tag, "random")

    def test_parse_tmrna_table_data_2(self):
        """Verify truncated MySQL tmRNA data is parsed correctly
        from a data dictionary returned from a SQL query."""
        data_dict = {
            "GeneID": "TRIXIE_0001",
            "PhageID": "Trixie",
            "Start": 100,
            "LocusTag": "SEA_TRIXIE_0001",
            }

        tmrna1 = mysqldb.parse_tmrna_table_data(data_dict)

        with self.subTest():
            self.assertEqual(tmrna1.id, "TRIXIE_0001")
        with self.subTest():
            self.assertEqual(tmrna1.genome_id, "Trixie")
        with self.subTest():
            self.assertEqual(tmrna1.start, 100)
        with self.subTest():
            self.assertEqual(tmrna1.locus_tag, "SEA_TRIXIE_0001")




class TestMysqldbFunctions2(unittest.TestCase):


    def setUp(self):
        self.genome1 = genome.Genome()
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

        self.cds1 = cds.Cds()
        self.cds1.genome_id = "L5"
        self.cds1.start = 10
        self.cds1.stop = 100
        self.cds1.parts = 1
        self.cds1.length = 1000
        self.cds1.name = "1"
        self.cds1.type = "CDS"
        self.cds1.translation = "AGGPT"
        self.cds1.orientation = "F"
        self.cds1.description = "description"
        self.cds1.locus_tag = "SEA_L5_001"

        self.cds2 = cds.Cds()
        self.cds2.genome_id = "L5"
        self.cds2.start = 100
        self.cds2.stop = 1000
        self.cds2.parts = 1
        self.cds2.length = 10000
        self.cds2.name = "2"
        self.cds2.type = "CDS"
        self.cds2.translation = "AKKQE"
        self.cds2.orientation = "R"
        self.cds2.description = "description"
        self.cds2.locus_tag = "SEA_L5_002"

        self.cds_features = [self.cds1, self.cds2]

        self.trna1 = trna.Trna()
        self.trna1.id = "Trixie_1"
        self.trna1.genome_id = "Trixie"
        self.trna1.name = "1"
        self.trna1.locus_tag = "TAG1"
        self.trna1.start = 5
        self.trna1.stop = 10
        self.trna1.length = 200
        self.trna1.orientation = "F"
        self.trna1.note = "misc"
        self.trna1.amino_acid = "Ala"
        self.trna1.anticodon = "AAA"
        self.trna1.structure = "random"
        self.trna1.use = "aragorn"

        self.trna2 = trna.Trna()
        self.trna2.id = "Trixie_1"
        self.trna2.genome_id = "Trixie"
        self.trna2.name = "1"
        self.trna2.locus_tag = "TAG1"
        self.trna2.start = 5
        self.trna2.stop = 10
        self.trna2.length = 200
        self.trna2.orientation = "F"
        self.trna2.note = "misc"
        self.trna2.amino_acid = "Ala"
        self.trna2.anticodon = "AAA"
        self.trna2.structure = "random"
        self.trna2.use = "aragorn"

        self.trna_features = [self.trna1, self.trna2]

        self.tmrna1 = tmrna.Tmrna()
        self.tmrna1.id = "Trixie_1"
        self.tmrna1.genome_id = "Trixie"
        self.tmrna1.name = "1"
        self.tmrna1.locus_tag = "TAG1"
        self.tmrna1.start = 5
        self.tmrna1.stop = 10
        self.tmrna1.length = 200
        self.tmrna1.orientation = "F"
        self.tmrna1.note = "misc"
        self.tmrna1.peptide_tag = "random"

        self.tmrna2 = tmrna.Tmrna()
        self.tmrna2.id = "Trixie_1"
        self.tmrna2.genome_id = "Trixie"
        self.tmrna2.name = "1"
        self.tmrna2.locus_tag = "TAG1"
        self.tmrna2.start = 5
        self.tmrna2.stop = 10
        self.tmrna2.length = 200
        self.tmrna2.orientation = "F"
        self.tmrna2.note = "misc"
        self.tmrna2.peptide_tag = "random"

        self.tmrna_features = [self.tmrna1, self.tmrna2]

    def test_create_genome_statements_1(self):
        """Verify list of INSERT statements is created correctly for:
        'add' ticket, and no CDS features."""
        statements = mysqldb.create_genome_statements(
                        self.genome1, tkt_type="add")
        self.assertEqual(len(statements), 1)

    def test_create_genome_statements_2(self):
        """Verify list of INSERT statements is created correctly for:
        'replace' ticket, and no CDS features."""
        statements = mysqldb.create_genome_statements(
                        self.genome1, tkt_type="replace")
        self.assertEqual(len(statements), 2)

    def test_create_genome_statements_3(self):
        """Verify list of INSERT statements is created correctly for:
        'add' ticket, two CDS features, two tRNA features, and
        two tmRNA features."""
        self.genome1.cds_features = self.cds_features
        self.genome1.trna_features = self.trna_features
        self.genome1.tmrna_features = self.tmrna_features
        statements = mysqldb.create_genome_statements(
                        self.genome1, tkt_type="add")
        self.assertEqual(len(statements), 7)

if __name__ == '__main__':
    unittest.main()
