""" Unit tests for misc. functions that interact with a MySQL database."""


import unittest
from pdm_utils.functions import phamerator
from pdm_utils.classes import genome
from pdm_utils.classes import cds
from datetime import datetime
from pdm_utils.classes import bundle
from pdm_utils.constants import constants
from Bio.Seq import Seq



class TestPhameratorFunctions1(unittest.TestCase):


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
                     "HostStrain":"Mycobacterium",
                     "Sequence":"ATCG".encode("utf-8"),
                     "SequenceLength":10,
                     "DateLastModified":constants.EMPTY_DATE,
                     "Notes":"abc".encode("utf-8"),
                     "GC":12.12,
                     "Cluster":"B",
                     "Cluster2":"A",
                     "Subcluster2":"A2",
                     "Status":"final",
                     "RetrieveRecord":1,
                     "AnnotationAuthor":1}

        self.genome1 = \
            phamerator.parse_phage_table_data(data_dict, gnm_type="phamerator")

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
            self.assertEqual(self.genome1.cluster_subcluster, "B")
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
            self.assertEqual(self.genome1.type, "phamerator")


    def test_parse_phage_table_data_2(self):
        """Verify truncated MySQL genome data is parsed correctly
        from a data dictionary returned from a SQL query."""

        data_dict = {"PhageID":"L5"}
        self.genome1 = phamerator.parse_phage_table_data(data_dict)
        with self.subTest():
            self.assertEqual(self.genome1.id, "L5")
        with self.subTest():
            self.assertEqual(self.genome1.name, "")
        with self.subTest():
            self.assertEqual(self.genome1.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome1.type, "")




    def test_convert_for_sql_1(self):
        """Verify non-empy value returned contains ''."""
        value = phamerator.convert_for_sql("A")
        self.assertEqual(value, "'A'")

    def test_convert_for_sql_2(self):
        """Verify empty value returned is NULL."""
        value = phamerator.convert_for_sql("")
        self.assertEqual(value, "NULL")

    def test_convert_for_sql_3(self):
        """Verify 'singleton' value returned is NULL."""
        value = phamerator.convert_for_sql("SINGLETON")
        self.assertEqual(value, "NULL")




    def test_parse_gene_table_data_1(self):
        """Verify standard MySQL CDS data is parsed correctly
        from a data dictionary returned from a SQL query."""

        data_dict = {"GeneID":"L5_001",
                     "PhageID":"L5",
                     "Start":10,
                     "Stop":100,
                     "Length":1000,
                     "Name":"1",
                     "Translation":"AGGPT",
                     "Orientation":"F",
                     "Notes":"description".encode("utf-8"),
                     "LocusTag":"SEA_L5_001"}

        cds1 = phamerator.parse_gene_table_data(data_dict)

        with self.subTest():
            self.assertEqual(cds1.id, "L5_001")
        with self.subTest():
            self.assertEqual(cds1.genome_id, "L5")
        with self.subTest():
            self.assertEqual(cds1.left, 10)
        with self.subTest():
            self.assertEqual(cds1.right, 100)
        with self.subTest():
            self.assertEqual(cds1.length, 1000)
        with self.subTest():
            self.assertEqual(cds1.name, "1")
        with self.subTest():
            self.assertEqual(cds1.type, "CDS")
        with self.subTest():
            self.assertEqual(cds1.translation, "AGGPT")
        with self.subTest():
            self.assertEqual(cds1.strand, "F")
        with self.subTest():
            self.assertEqual(cds1.description, "description")
        with self.subTest():
            self.assertEqual(cds1.locus_tag, "SEA_L5_001")
        with self.subTest():
            self.assertEqual(cds1.translation_table, 11)
        with self.subTest():
            self.assertEqual(cds1.coordinate_format, "0_half_open")


    def test_parse_gene_table_data_2(self):
        """Verify truncated MySQL CDS data is parsed correctly
        from a data dictionary returned from a SQL query."""

        data_dict = {"GeneID":"L5_001",
                     "PhageID":"L5",
                     "Start":10}
        cds1 = phamerator.parse_gene_table_data(data_dict)
        with self.subTest():
            self.assertEqual(cds1.id, "L5_001")
        with self.subTest():
            self.assertEqual(cds1.genome_id, "L5")
        with self.subTest():
            self.assertEqual(cds1.left, 10)
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
        self.genome1.cluster_subcluster = "A123"
        self.genome1.cluster = "A"
        self.genome1.subcluster = "A2"

        statements = phamerator.create_genome_statements(
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
        self.genome1.cluster_subcluster = "A123"
        self.genome1.cluster = "A"
        self.genome1.subcluster = "A2"

        statements = phamerator.create_genome_statements(
                        self.genome1, tkt_type="replace")
        self.assertEqual(len(statements), 2)


    def test_create_genome_statements_3(self):
        """Verify list of INSERT statements is created correctly for:
        'add' ticket, and two CDS features."""

        cds1 = cds.Cds()
        cds1.genome_id = "L5"
        cds1.left = 10
        cds1.right = 100
        cds1.length = 1000
        cds1.name = "1"
        cds1.type = "CDS"
        cds1.translation = "AGGPT"
        cds1.strand = "F"
        cds1.description = "description"
        cds1.locus_tag = "SEA_L5_001"

        cds2 = cds.Cds()
        cds2.genome_id = "L5"
        cds2.left = 100
        cds2.right = 1000
        cds2.length = 10000
        cds2.name = "2"
        cds2.type = "CDS"
        cds2.translation = "AKKQE"
        cds2.strand = "R"
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
        self.genome1.cluster_subcluster = "A123"
        self.genome1.cluster = "A"
        self.genome1.subcluster = "A2"
        self.genome1.cds_features = [cds1, cds2]

        statements = phamerator.create_genome_statements(
                        self.genome1, tkt_type="add")
        self.assertEqual(len(statements), 3)




# class TestPhameratorFunctions2(unittest.TestCase):
#
#     def setUp(self):
#
#
#         self.genome1 = genome.Genome()
#         self.genome1.id = "L5"
#         self.genome1.type = "add"
#         self.genome1.host_genus = "Gordonia"
#         self.genome1.cluster = "B"
#         self.genome1._value_flag = True
#
#         self.bundle1 = bundle.Bundle()
#
#
#         self.genome2 = genome.Genome()
#         self.genome2.id = "L5"
#         self.genome2.type = "phamerator"
#         self.genome2.host_genus = "Mycobacterium"
#         self.genome2.cluster = "A"
#
#     def test_copy_data_1(self):
#         """Check that an "add" genome with no fields set to 'retain' is
#         not impacted."""
#
#         self.bundle1.genome_dict[self.genome1.type] = self.genome1
#         phamerator.copy_data(self.bundle1, "phamerator", "add")
#         genome1 = self.bundle1.genome_dict["add"]
#         with self.subTest():
#             self.assertFalse(genome1._value_flag)
#         with self.subTest():
#             self.assertEqual(genome1.host_genus, "Gordonia")
#         with self.subTest():
#             self.assertEqual(genome1.cluster, "B")
#
#     def test_copy_data_2(self):
#         """Check that an "add" genome with host_genus field set to 'retain' is
#         populated correctly."""
#
#         self.bundle1.genome_dict[self.genome1.type] = self.genome1
#         self.genome1.host_genus = "retain"
#         self.bundle1.genome_dict[self.genome2.type] = self.genome2
#         phamerator.copy_data(self.bundle1, "phamerator", "add")
#         genome1 = self.bundle1.genome_dict["add"]
#         with self.subTest():
#             self.assertFalse(genome1._value_flag)
#         with self.subTest():
#             self.assertEqual(genome1.host_genus, "Mycobacterium")
#         with self.subTest():
#             self.assertEqual(genome1.cluster, "B")
#
#     def test_copy_data_3(self):
#         """Check that an "invalid" genome with host_genus field set to 'retain' is
#         not populated correctly."""
#
#         self.genome1.type = "invalid"
#         self.bundle1.genome_dict[self.genome1.type] = self.genome1
#         self.genome1.host_genus = "retain"
#         phamerator.copy_data(self.bundle1, "phamerator", "add")
#         with self.subTest():
#             self.assertEqual(
#                 len(self.bundle1.genome_pair_dict.keys()), 0)
#         with self.subTest():
#             self.assertEqual(self.genome1.host_genus, "retain")
#
#     def test_copy_data_4(self):
#         """Check that an "add" genome with host_genus field set to 'retain' is
#         not populated correctly when "invalid" type is requested."""
#
#         self.bundle1.genome_dict[self.genome1.type] = self.genome1
#         self.genome1.host_genus = "retain"
#         phamerator.copy_data(self.bundle1, "phamerator", "invalid")
#         with self.subTest():
#             self.assertEqual(
#                 len(self.bundle1.genome_pair_dict.keys()), 0)
#         with self.subTest():
#             self.assertEqual(self.genome1.host_genus, "retain")
#
#     def test_copy_data_5(self):
#         """Check that an "add" genome with host_genus field set to 'retain' is
#         not populated correctly when there is no matching "phamerator"
#         genomet type."""
#
#         self.bundle1.genome_dict[self.genome1.type] = self.genome1
#         self.genome1.host_genus = "retain"
#         self.genome1._value_flag = False
#         phamerator.copy_data(self.bundle1, "phamerator", "add")
#         with self.subTest():
#             self.assertTrue(self.genome1._value_flag)
#         with self.subTest():
#             self.assertEqual(
#                 len(self.bundle1.genome_pair_dict.keys()), 0)
#         with self.subTest():
#             self.assertEqual(self.genome1.host_genus, "retain")













if __name__ == '__main__':
    unittest.main()
