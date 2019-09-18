""" Unit tests for misc. functions that interact with PhameratorDB."""


import unittest
from pdm_utils.functions import phamerator
from pdm_utils.classes import genome
from pdm_utils.classes import cds
from datetime import datetime
from pdm_utils.classes import bundle
from pdm_utils.constants import constants
from Bio.Seq import Seq



class TestPhameratorFunctions(unittest.TestCase):


    def setUp(self):
        self.genome1 = genome.Genome()
        self.genome2 = genome.Genome()
        self.genome3 = genome.Genome()

        self.data_tuple1 = ("L5",
                            "L5",
                            "Mycobacterium",
                            "ATCG",
                            "final",
                            "A",
                            "1/1/1900",
                            "ABC123",
                            "A2",
                            "1",
                            "1",
                            "1")

        self.data_tuple2 = ("Trixie",
                            "Trixie",
                            "Mycobacterium",
                            "TTTT",
                            "final",
                            "A",
                            "1/1/1900",
                            "EFG123",
                            "A2",
                            "1",
                            "1",
                            "1")

        self.data_tuple3 = ("KatherineG",
                            "KatherineG",
                            "Gordonia",
                            "CCCC",
                            "final",
                            "A",
                            "1/1/1900",
                            "XYZ123",
                            "A15",
                            "1",
                            "1",
                            "1")

        self.data_tuple4 = ("XYZ",
                            "XYZ_Draft",
                            "Arthrobacter",
                            "GGGG",
                            "gbk",
                            "X",
                            "1/1/1900",
                            "none",
                            "",
                            "0",
                            "0",
                            "0")

        self.cds1 = cds.Cds()
        self.cds2 = cds.Cds()
        self.cds3 = cds.Cds()




    # # TODO this may no longer be needed now that
    # # parse_phage_table_data() is available.
    # def test_parse_phamerator_data_1(self):
    #     """Verify standard Phamerator genome data is parsed correctly."""
    #
    #     input_phage_id = "Trixie_Draft"
    #     input_phage_name = "Trixie_Draft"
    #     input_host = "   Mycobacterium smegmatis  "
    #     input_sequence = "atcg"
    #     input_status = "final"
    #     input_cluster = "A"
    #     input_subcluster = "A2"
    #     input_date = datetime.strptime('1/1/2000', '%m/%d/%Y')
    #     input_accession = "  ABC123.1  "
    #     input_annotation_author = "1"
    #     input_annotation_qc = "1"
    #     input_retrieve_record = "1"
    #
    #     data_tuple = (input_phage_id,
    #                     input_phage_name,
    #                     input_host,
    #                     input_sequence,
    #                     input_status,
    #                     input_cluster,
    #                     input_date,
    #                     input_accession,
    #                     input_subcluster,
    #                     input_annotation_author,
    #                     input_annotation_qc,
    #                     input_retrieve_record)
    #
    #     phamerator.parse_phamerator_data(self.genome1, data_tuple)
    #
    #     output_phage_id = "Trixie"
    #     output_phage_name = "Trixie_Draft"
    #     output_host = "Mycobacterium"
    #     output_sequence = "ATCG"
    #     output_status = "final"
    #     output_cluster = "A"
    #     output_subcluster = "A2"
    #     output_date = datetime.strptime('1/1/2000', '%m/%d/%Y')
    #     output_accession = "ABC123"
    #     output_annotation_author = "1"
    #     output_annotation_qc = "1"
    #     output_retrieve_record = "1"
    #     output_type = "phamerator"
    #     output_seq_length = 4
    #
    #     with self.subTest():
    #         self.assertEqual(self.genome1.id, output_phage_id)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.name, output_phage_name)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.host_genus, output_host)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.seq, output_sequence)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.annotation_status, output_status)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.cluster, output_cluster)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.subcluster, output_subcluster)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.date, \
    #             output_date)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.accession, output_accession)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.annotation_author, \
    #             output_annotation_author)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.annotation_qc, output_annotation_qc)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.retrieve_record, \
    #             output_retrieve_record)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.length, output_seq_length)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.type, output_type)
    #
    #
    # def test_parse_phamerator_data_2(self):
    #     """Verify empty Phamerator genome data is parsed correctly."""
    #
    #     input_phage_id = "Trixie_Draft"
    #     input_phage_name = "Trixie_Draft"
    #     input_host = ""
    #     input_sequence = "atcg"
    #     input_status = "final"
    #     input_cluster = "SINGLETON"
    #     input_subcluster = "NONE"
    #     input_date = None
    #     input_accession = None
    #     input_annotation_author = "1"
    #     input_annotation_qc = "1"
    #     input_retrieve_record = "1"
    #
    #     data_tuple = (input_phage_id,
    #                     input_phage_name,
    #                     input_host,
    #                     input_sequence,
    #                     input_status,
    #                     input_cluster,
    #                     input_date,
    #                     input_accession,
    #                     input_subcluster,
    #                     input_annotation_author,
    #                     input_annotation_qc,
    #                     input_retrieve_record)
    #
    #     phamerator.parse_phamerator_data(self.genome1, data_tuple)
    #
    #     output_phage_id = "Trixie"
    #     output_phage_name = "Trixie_Draft"
    #     output_host = ""
    #     output_sequence = "ATCG"
    #     output_status = "final"
    #     output_cluster = "singleton"
    #     output_subcluster = ""
    #     output_date = datetime.strptime('1/1/0001', '%m/%d/%Y')
    #     output_accession = ""
    #     output_annotation_author = "1"
    #     output_annotation_qc = "1"
    #     output_retrieve_record = "1"
    #     output_type = "phamerator"
    #     output_seq_length = 4
    #
    #     with self.subTest():
    #         self.assertEqual(self.genome1.id, output_phage_id)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.name, output_phage_name)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.host_genus, output_host)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.seq, output_sequence)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.annotation_status, output_status)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.cluster, output_cluster)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.subcluster, output_subcluster)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.date, \
    #             output_date)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.accession, output_accession)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.annotation_author, \
    #             output_annotation_author)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.annotation_qc, output_annotation_qc)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.retrieve_record, \
    #             output_retrieve_record)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.length, output_seq_length)
    #     with self.subTest():
    #         self.assertEqual(self.genome1.type, output_type)



















    def test_parse_phage_table_data_1(self):
        """Verify standard Phamerator genome data is parsed correctly
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
                     "status":"final",
                     "RetrieveRecord":1,
                     "AnnotationQC":1,
                     "AnnotationAuthor":1}

        self.genome1 = phamerator.parse_phage_table_data(data_dict)

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
            self.assertEqual(self.genome1.annotation_qc, 1)
        with self.subTest():
            self.assertEqual(self.genome1.annotation_author, 1)
        with self.subTest():
            self.assertEqual(self.genome1.translation_table, 11)
        with self.subTest():
            self.assertEqual(self.genome1.type, "phamerator")


    def test_parse_phage_table_data_2(self):
        """Verify truncated Phamerator genome data is parsed correctly
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
            self.assertEqual(self.genome1.type, "phamerator")














    def test_parse_gene_table_data_1(self):
        """Verify standard Phamerator CDS data is parsed correctly
        from a data dictionary returned from a SQL query."""

        data_dict = {"GeneID":"L5_001",
                     "PhageID":"L5",
                     "Start":10,
                     "Stop":100,
                     "Length":1000,
                     "Name":"1",
                     "TypeID":"CDS",
                     "translation":"AGGPT",
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


    def test_parse_gene_table_data_2(self):
        """Verify truncated Phamerator CDS data is parsed correctly
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









    ## TODO this may no longer be needed.
    #
    # def test_create_phamerator_dict_1(self):
    #     """Verify Phamerator MySQL query output is parsed correctly."""
    #
    #     data_tuples = (self.data_tuple1, self.data_tuple2, self.data_tuple3)
    #     genome_dict = phamerator.create_phamerator_dict(data_tuples)
    #     genome_l5 = genome_dict["L5"]
    #     genome_trixie = genome_dict["Trixie"]
    #     genome_katherineg = genome_dict["KatherineG"]
    #
    #     with self.subTest():
    #         self.assertEqual(len(genome_dict.keys()), 3)
    #     with self.subTest():
    #         self.assertEqual(genome_l5.accession, "ABC123")
    #     with self.subTest():
    #         self.assertEqual(genome_trixie.accession, "EFG123")
    #     with self.subTest():
    #         self.assertEqual(genome_katherineg.accession, "XYZ123")
    #
    #
    #
    #
    #
    #
    # def test_create_data_sets_1(self):
    #     """Verify multiple sets of unique Phamerator data are produced.
    #     Verify that empty accession is not added.
    #     Verify that empty subcluster is not added."""
    #
    #     data_tuples = (self.data_tuple1,
    #                     self.data_tuple2,
    #                     self.data_tuple3,
    #                     self.data_tuple4)
    #     genome_dict = phamerator.create_phamerator_dict(data_tuples)
    #     returned_dict = phamerator.create_data_sets(genome_dict)
    #
    #     exp_ids = set(["L5","Trixie","KatherineG","XYZ"])
    #     exp_host = set(["Mycobacterium","Gordonia","Arthrobacter"])
    #     exp_status = set(["final","gbk"])
    #     exp_cluster = set(["A","X"])
    #     exp_subcluster = set(["A2","A15"])
    #     exp_accession = set(["ABC123","EFG123","XYZ123"])
    #
    #     with self.subTest():
    #         self.assertEqual(len(returned_dict.keys()), 6)
    #     with self.subTest():
    #         self.assertEqual(returned_dict["phage_id"], exp_ids)
    #     with self.subTest():
    #         self.assertEqual(returned_dict["host_genus"], exp_host)
    #     with self.subTest():
    #         self.assertEqual(returned_dict["annotation_status"], exp_status)
    #     with self.subTest():
    #         self.assertEqual(returned_dict["cluster"], exp_cluster)
    #     with self.subTest():
    #         self.assertEqual(returned_dict["subcluster"], exp_subcluster)
    #     with self.subTest():
    #         self.assertEqual(returned_dict["accession"], exp_accession)




    def test_create_update_statement_1(self):
        """Verify correct Cluster statement is created for a non-singleton."""
        statement = phamerator.create_update_statement( \
            "phage", "PhageID", "L5", "Cluster", "A")
        exp_statement = \
            "UPDATE phage SET Cluster = 'A' WHERE PhageID = 'L5';"
        with self.subTest():
            self.assertEqual(statement, exp_statement)

    def test_create_update_statement_2(self):
        """Verify correct Cluster statement is created for a singleton."""
        statement = phamerator.create_update_statement( \
                "phage", "PhageID", "L5", "Cluster", "SINGLETON")
        exp_statement = \
            "UPDATE phage SET Cluster = NULL WHERE PhageID = 'L5';"
        with self.subTest():
            self.assertEqual(statement, exp_statement)

    def test_create_update_statement_3(self):
        """Verify correct Cluster2 statement is created for a singleton."""
        statement = phamerator.create_update_statement( \
                "phage", "PhageID", "L5", "Cluster2", "SINGLETON")
        exp_statement = \
            "UPDATE phage SET Cluster2 = NULL WHERE PhageID = 'L5';"
        with self.subTest():
            self.assertEqual(statement, exp_statement)

    def test_create_update_statement_4(self):
        """Verify correct Subcluster2 statement is created for a
        non-empty value."""
        statement = phamerator.create_update_statement( \
                "phage", "PhageID", "L5", "Subcluster2", "A2")
        exp_statement = \
            "UPDATE phage SET Subcluster2 = 'A2' WHERE PhageID = 'L5';"
        with self.subTest():
            self.assertEqual(statement, exp_statement)

    def test_create_update_statement_5(self):
        """Verify correct Subcluster2 statement is created for an
        empty value."""
        statement = phamerator.create_update_statement( \
                "phage", "PhageID", "L5", "Subcluster2", "none")
        exp_statement = \
            "UPDATE phage SET Subcluster2 = NULL WHERE PhageID = 'L5';"
        with self.subTest():
            self.assertEqual(statement, exp_statement)

    def test_create_cluster_statement_6(self):
        """Verify Gene statement is created correctly."""
        statement = phamerator.create_update_statement( \
            "gene", "GeneID", "SEA_L5_123", "Notes", "Integrase")
        exp_statement = \
            "UPDATE gene SET Notes = 'Integrase' WHERE GeneID = 'SEA_L5_123';"
        with self.subTest():
            self.assertEqual(statement, exp_statement)




    def test_create_genome_update_statements_1(self):
        """Verify list of statements are created correctly."""

        self.genome1.id = "L5"
        self.genome1.host_genus = "Mycobacterium"
        self.genome1.annotation_status = "final"
        self.genome1.accession = "ABC123"
        self.genome1.author = "1"
        self.genome1.cluster_subcluster = "A123"
        self.genome1.cluster = "A"
        self.genome1.subcluster = "A2"

        statements = phamerator.create_genome_update_statements(self.genome1)

        end = " WHERE PhageID = 'L5';"
        exp_statements = 7
        exp_host = "UPDATE phage SET HostStrain = 'Mycobacterium'" + end
        exp_status = "UPDATE phage SET status = 'final'" + end
        exp_accession = "UPDATE phage SET Accession = 'ABC123'" + end
        exp_author = "UPDATE phage SET AnnotationAuthor = '1'" + end
        exp_clust_sub = "UPDATE phage SET Cluster = 'A123'" + end
        exp_cluster = "UPDATE phage SET Cluster2 = 'A'" + end
        exp_subcluster = "UPDATE phage SET Subcluster2 = 'A2'" + end


        with self.subTest():
            self.assertEqual(len(statements), 7)
        with self.subTest():
            self.assertEqual(statements[0], exp_host)
        with self.subTest():
            self.assertEqual(statements[1], exp_status)
        with self.subTest():
            self.assertEqual(statements[2], exp_accession)
        with self.subTest():
            self.assertEqual(statements[3], exp_author)
        with self.subTest():
            self.assertEqual(statements[4], exp_clust_sub)
        with self.subTest():
            self.assertEqual(statements[5], exp_cluster)
        with self.subTest():
            self.assertEqual(statements[6], exp_subcluster)




    def test_create_delete_statement_1(self):
        """Verify correct DELETE statement is created."""
        statement = \
            phamerator.create_delete_statement("phage", "PhageID", "L5")
        exp_statement = "DELETE FROM phage WHERE PhageID = 'L5';"
        with self.subTest():
            self.assertEqual(statement, exp_statement)




    def test_create_genome_delete_statement_1(self):
        """Verify list of DELETE statements are created correctly."""

        self.genome1.id = "L5"
        statements = phamerator.create_genome_delete_statement(self.genome1)
        exp_statements = 1
        exp = "DELETE FROM phage WHERE PhageID = 'L5';"

        with self.subTest():
            self.assertEqual(len(statements), 1)
        with self.subTest():
            self.assertEqual(statements[0], exp)




    def test_create_cds_insert_statement_1(self):
        """Verify CDS INSERT statement is created correctly."""

        self.cds1.id = "SEA_L5_123"
        self.cds1.genome_id = "L5"
        self.cds1.left = 5
        self.cds1.right = 10
        self.cds1._translation_length = 20
        self.cds1.name = "Int"
        self.cds1.type = "CDS"
        self.cds1.translation = "ACKLG"
        self.cds1.strand = "forward"
        self.cds1.processed_description = "integrase"
        self.cds1.locus_tag = "TAG1"

        statement = phamerator.create_cds_insert_statement(self.cds1)

        exp = "INSERT INTO gene " + \
        "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
        "translation, Orientation, Notes, LocusTag) " + \
        "VALUES " + \
        "('SEA_L5_123', 'L5', 5, 10, 20, 'Int', 'CDS', " + \
        "'ACKLG', 'forward', 'integrase', 'TAG1');"

        self.assertEqual(statement, exp)




    def test_create_cds_insert_statements_1(self):
        """Verify list of CDS INSERT statements is created correctly."""

        self.cds1.id = "SEA_L5_123"
        self.cds1.genome_id = "L5"
        self.cds1.left = 5
        self.cds1.right = 10
        self.cds1._translation_length = 20
        self.cds1.name = "Int"
        self.cds1.type = "CDS"
        self.cds1.translation = "ACKLG"
        self.cds1.strand = "forward"
        self.cds1.processed_description = "integrase"
        self.cds1.locus_tag = "TAG1"

        self.cds2.id = "SEA_TRIXIE_123"
        self.cds2.genome_id = "Trixie"
        self.cds2.left = 1
        self.cds2.right = 100
        self.cds2._translation_length = 15
        self.cds2.name = "parA"
        self.cds2.type = "CDS"
        self.cds2.translation = "PPGLA"
        self.cds2.strand = "reverse"
        self.cds2.processed_description = "parA"
        self.cds2.locus_tag = "TAG2"

        self.cds3.id = "SEA_D29_123"
        self.cds3.genome_id = "D29"
        self.cds3.left = 42
        self.cds3.right = 62
        self.cds3._translation_length = 3
        self.cds3.name = "D29_123"
        self.cds3.type = "CDS"
        self.cds3.translation = "AVY"
        self.cds3.strand = "reverse"
        self.cds3.processed_description = ""
        self.cds3.locus_tag = ""

        cds_features = [self.cds1, self.cds2, self.cds3]

        statements = phamerator.create_cds_insert_statements(cds_features)

        exp1 = "INSERT INTO gene " + \
        "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
        "translation, Orientation, Notes, LocusTag) " + \
        "VALUES " + \
        "('SEA_L5_123', 'L5', 5, 10, 20, 'Int', 'CDS', " + \
        "'ACKLG', 'forward', 'integrase', 'TAG1');"

        exp2 = "INSERT INTO gene " + \
        "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
        "translation, Orientation, Notes, LocusTag) " + \
        "VALUES " + \
        "('SEA_TRIXIE_123', 'Trixie', 1, 100, 15, 'parA', 'CDS', " + \
        "'PPGLA', 'reverse', 'parA', 'TAG2');"

        exp3 = "INSERT INTO gene " + \
        "(GeneID, PhageID, Start, Stop, Length, Name, TypeID, " + \
        "translation, Orientation, Notes, LocusTag) " + \
        "VALUES " + \
        "('SEA_D29_123', 'D29', 42, 62, 3, 'D29_123', 'CDS', " + \
        "'AVY', 'reverse', '', '');"

        with self.subTest():
            self.assertEqual(len(statements), 3)
        with self.subTest():
            self.assertEqual(statements[0], exp1)
        with self.subTest():
            self.assertEqual(statements[1], exp2)
        with self.subTest():
            self.assertEqual(statements[2], exp3)




    def test_create_genome_insert_statement_1(self):
        """Verify genome INSERT statement is created correctly."""

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
        self.genome1.annotation_qc = "1"
        self.genome1.annotation_author = "1"

        statement = phamerator.create_genome_insert_statement(self.genome1)

        exp = "INSERT INTO phage " + \
        "(PhageID, Accession, Name, HostStrain, Sequence, " + \
        "SequenceLength, GC, status, DateLastModified, RetrieveRecord, " + \
        "AnnotationQC, AnnotationAuthor) " + \
        "VALUES ('L5', 'ABC123', 'L5_Draft', 'Mycobacterium', 'ATCG', 4, " + \
        "0.5001, 'final', '1/1/2000', '1', '1', '1');"

        self.assertEqual(statement, exp)




    def test_create_genome_insert_statements_1(self):
        """Verify list of genome INSERT statements is created correctly."""

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
        self.genome1.annotation_qc = "1"
        self.genome1.annotation_author = "1"

        self.genome1.cluster_subcluster = "A123"
        self.genome1.cluster = "A"
        self.genome1.subcluster = "A2"

        statements = phamerator.create_genome_insert_statements(self.genome1)
        self.assertEqual(len(statements), 4)








class TestPhameratorFunctions2(unittest.TestCase):

    def setUp(self):


        self.genome1 = genome.Genome()
        self.genome1.id = "L5"
        self.genome1.type = "add"
        self.genome1.host_genus = "Gordonia"
        self.genome1.cluster = "B"
        self.genome1._value_flag = True

        self.bundle1 = bundle.Bundle()


        self.genome2 = genome.Genome()
        self.genome2.id = "L5"
        self.genome2.type = "phamerator"
        self.genome2.host_genus = "Mycobacterium"
        self.genome2.cluster = "A"

    def test_copy_data_from_1(self):
        """Check that an "add" genome with no fields set to 'retain' is
        not impacted."""

        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        phamerator.copy_data_from(self.bundle1, "add")
        genome1 = self.bundle1.genome_dict["add"]
        with self.subTest():
            self.assertFalse(genome1._value_flag)
        with self.subTest():
            self.assertEqual(genome1.host_genus, "Gordonia")
        with self.subTest():
            self.assertEqual(genome1.cluster, "B")
        with self.subTest():
            self.assertEqual(genome1.evaluations[0].status, "correct")

    def test_copy_data_from_2(self):
        """Check that an "add" genome with host_genus field set to 'retain' is
        populated correctly."""

        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.genome1.host_genus = "retain"
        self.bundle1.genome_dict[self.genome2.type] = self.genome2
        phamerator.copy_data_from(self.bundle1, "add")
        genome1 = self.bundle1.genome_dict["add"]
        with self.subTest():
            self.assertFalse(genome1._value_flag)
        with self.subTest():
            self.assertEqual(genome1.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(genome1.cluster, "B")
        with self.subTest():
            self.assertEqual(genome1.evaluations[0].status, "correct")

    def test_copy_data_from_3(self):
        """Check that an "invalid" genome with host_genus field set to 'retain' is
        not populated correctly."""

        self.genome1.type = "invalid"
        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.genome1.host_genus = "retain"
        phamerator.copy_data_from(self.bundle1, "add")
        with self.subTest():
            self.assertEqual(
                len(self.bundle1.genome_pair_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(self.genome1.host_genus, "retain")
        with self.subTest():
            self.assertEqual(len(self.genome1.evaluations), 0)

    def test_copy_data_from_4(self):
        """Check that an "add" genome with host_genus field set to 'retain' is
        not populated correctly when "invalid" type is requested."""

        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.genome1.host_genus = "retain"
        phamerator.copy_data_from(self.bundle1, "invalid")
        with self.subTest():
            self.assertEqual(
                len(self.bundle1.genome_pair_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(self.genome1.host_genus, "retain")
        with self.subTest():
            self.assertEqual(len(self.genome1.evaluations), 0)

    def test_copy_data_from_5(self):
        """Check that an "add" genome with host_genus field set to 'retain' is
        not populated correctly when there is no matching "phamerator"
        genomet type."""

        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.genome1.host_genus = "retain"
        self.genome1._value_flag = False
        phamerator.copy_data_from(self.bundle1, "add")
        with self.subTest():
            self.assertTrue(self.genome1._value_flag)
        with self.subTest():
            self.assertEqual(
                len(self.bundle1.genome_pair_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(self.genome1.host_genus, "retain")
        with self.subTest():
            self.assertEqual(self.genome1.evaluations[0].status, "error")













if __name__ == '__main__':
    unittest.main()
