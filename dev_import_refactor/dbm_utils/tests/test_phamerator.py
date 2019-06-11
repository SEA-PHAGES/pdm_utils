""" Unit tests for misc. functions that interact with PhameratorDB."""


import unittest
from functions import phamerator
from classes import Genome
from datetime import datetime


class TestTicketFunctions1(unittest.TestCase):


    def setUp(self):
        self.genome1 = Genome.Genome()
        self.genome2 = Genome.Genome()
        self.genome3 = Genome.Genome()

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


    def test_parse_phamerator_data_1(self):
        """Verify standard Phamerator genome data is parsed correctly."""

        input_phage_id = "Trixie_Draft"
        input_phage_name = "Trixie_Draft"
        input_host = "   Mycobacterium smegmatis  "
        input_sequence = "atcg"
        input_status = "final"
        input_cluster = "A"
        input_subcluster = "A2"
        input_date_last_modified = datetime.strptime('1/1/2000', '%m/%d/%Y')
        input_accession = "  ABC123.1  "
        input_annotation_author = "1"
        input_annotation_qc = "1"
        input_retrieve_record = "1"

        data_tuple = (input_phage_id,
                        input_phage_name,
                        input_host,
                        input_sequence,
                        input_status,
                        input_cluster,
                        input_date_last_modified,
                        input_accession,
                        input_subcluster,
                        input_annotation_author,
                        input_annotation_qc,
                        input_retrieve_record)

        phamerator.parse_phamerator_data(self.genome1, data_tuple)

        output_phage_id = "Trixie"
        output_phage_name = "Trixie_Draft"
        output_host = "Mycobacterium"
        output_sequence = "ATCG"
        output_status = "final"
        output_cluster = "A"
        output_subcluster = "A2"
        output_date_last_modified = datetime.strptime('1/1/2000', '%m/%d/%Y')
        output_accession = "ABC123"
        output_annotation_author = "1"
        output_annotation_qc = "1"
        output_retrieve_record = "1"

        output_search_id = "trixie"
        output_seq_length = 4

        with self.subTest():
            self.assertEqual(self.genome1.phage_id, output_phage_id)
        with self.subTest():
            self.assertEqual(self.genome1.phage_name, output_phage_name)
        with self.subTest():
            self.assertEqual(self.genome1.host, output_host)
        with self.subTest():
            self.assertEqual(self.genome1.sequence, output_sequence)
        with self.subTest():
            self.assertEqual(self.genome1.status, output_status)
        with self.subTest():
            self.assertEqual(self.genome1.cluster, output_cluster)
        with self.subTest():
            self.assertEqual(self.genome1.subcluster, output_subcluster)
        with self.subTest():
            self.assertEqual(self.genome1.date_last_modified, \
                output_date_last_modified)
        with self.subTest():
            self.assertEqual(self.genome1.accession, output_accession)
        with self.subTest():
            self.assertEqual(self.genome1.annotation_author, \
                output_annotation_author)
        with self.subTest():
            self.assertEqual(self.genome1.annotation_qc, output_annotation_qc)
        with self.subTest():
            self.assertEqual(self.genome1.retrieve_record, \
                output_retrieve_record)
        with self.subTest():
            self.assertEqual(self.genome1.search_id, output_search_id)
        with self.subTest():
            self.assertEqual(self.genome1._length, output_seq_length)


    def test_parse_phamerator_data_2(self):
        """Verify empty Phamerator genome data is parsed correctly."""

        input_phage_id = "Trixie_Draft"
        input_phage_name = "Trixie_Draft"
        input_host = ""
        input_sequence = "atcg"
        input_status = "final"
        input_cluster = "SINGLETON"
        input_subcluster = "NONE"
        input_date_last_modified = None
        input_accession = None
        input_annotation_author = "1"
        input_annotation_qc = "1"
        input_retrieve_record = "1"

        data_tuple = (input_phage_id,
                        input_phage_name,
                        input_host,
                        input_sequence,
                        input_status,
                        input_cluster,
                        input_date_last_modified,
                        input_accession,
                        input_subcluster,
                        input_annotation_author,
                        input_annotation_qc,
                        input_retrieve_record)

        phamerator.parse_phamerator_data(self.genome1, data_tuple)

        output_phage_id = "Trixie"
        output_phage_name = "Trixie_Draft"
        output_host = ""
        output_sequence = "ATCG"
        output_status = "final"
        output_cluster = "singleton"
        output_subcluster = ""
        output_date_last_modified = datetime.strptime('1/1/0001', '%m/%d/%Y')
        output_accession = ""
        output_annotation_author = "1"
        output_annotation_qc = "1"
        output_retrieve_record = "1"

        output_search_id = "trixie"
        output_seq_length = 4

        with self.subTest():
            self.assertEqual(self.genome1.phage_id, output_phage_id)
        with self.subTest():
            self.assertEqual(self.genome1.phage_name, output_phage_name)
        with self.subTest():
            self.assertEqual(self.genome1.host, output_host)
        with self.subTest():
            self.assertEqual(self.genome1.sequence, output_sequence)
        with self.subTest():
            self.assertEqual(self.genome1.status, output_status)
        with self.subTest():
            self.assertEqual(self.genome1.cluster, output_cluster)
        with self.subTest():
            self.assertEqual(self.genome1.subcluster, output_subcluster)
        with self.subTest():
            self.assertEqual(self.genome1.date_last_modified, \
                output_date_last_modified)
        with self.subTest():
            self.assertEqual(self.genome1.accession, output_accession)
        with self.subTest():
            self.assertEqual(self.genome1.annotation_author, \
                output_annotation_author)
        with self.subTest():
            self.assertEqual(self.genome1.annotation_qc, output_annotation_qc)
        with self.subTest():
            self.assertEqual(self.genome1.retrieve_record, \
                output_retrieve_record)
        with self.subTest():
            self.assertEqual(self.genome1.search_id, output_search_id)
        with self.subTest():
            self.assertEqual(self.genome1._length, output_seq_length)









    def test_create_phamerator_dict_1(self):
        """Verify Phamerator MySQL query output is parsed correctly."""

        data_tuples = (self.data_tuple1, self.data_tuple2, self.data_tuple3)
        genome_dict = phamerator.create_phamerator_dict(data_tuples)
        genome_l5 = genome_dict["L5"]
        genome_trixie = genome_dict["Trixie"]
        genome_katherineg = genome_dict["KatherineG"]

        with self.subTest():
            self.assertEqual(len(genome_dict.keys()), 3)
        with self.subTest():
            self.assertEqual(genome_l5.accession, "ABC123")
        with self.subTest():
            self.assertEqual(genome_trixie.accession, "EFG123")
        with self.subTest():
            self.assertEqual(genome_katherineg.accession, "XYZ123")






    def test_create_data_sets_1(self):
        """Verify multiple sets of unique Phamerator data are produced.
        Verify that empty accession is not added.
        Verify that empty subcluster is not added."""

        data_tuples = (self.data_tuple1,
                        self.data_tuple2,
                        self.data_tuple3,
                        self.data_tuple4)
        genome_dict = phamerator.create_phamerator_dict(data_tuples)
        returned_dict = phamerator.create_data_sets(genome_dict)

        exp_ids = set(["L5","Trixie","KatherineG","XYZ"])
        exp_host = set(["Mycobacterium","Gordonia","Arthrobacter"])
        exp_status = set(["final","gbk"])
        exp_cluster = set(["A","X"])
        exp_subcluster = set(["A2","A15"])
        exp_accession = set(["ABC123","EFG123","XYZ123"])

        with self.subTest():
            self.assertEqual(len(returned_dict.keys()), 6)
        with self.subTest():
            self.assertEqual(returned_dict["phage_id"], exp_ids)
        with self.subTest():
            self.assertEqual(returned_dict["host"], exp_host)
        with self.subTest():
            self.assertEqual(returned_dict["status"], exp_status)
        with self.subTest():
            self.assertEqual(returned_dict["cluster"], exp_cluster)
        with self.subTest():
            self.assertEqual(returned_dict["subcluster"], exp_subcluster)
        with self.subTest():
            self.assertEqual(returned_dict["accession"], exp_accession)




    def test_create_cluster_statement_1(self):
        """Verify correct cluster statement is created for a non-singleton."""
        statement = phamerator.create_cluster_statement("L5", "A")
        exp_statement = \
            """UPDATE phage SET Cluster = 'A' WHERE PhageID = 'L5';"""
        with self.subTest():
            self.assertEqual(statement, exp_statement)

    def test_create_cluster_statement_2(self):
        """Verify correct cluster statement is created for a singleton."""
        statement = phamerator.create_cluster_statement("L5", "SINGLETON")
        exp_statement = \
            """UPDATE phage SET Cluster = NULL WHERE PhageID = 'L5';"""
        with self.subTest():
            self.assertEqual(statement, exp_statement)









if __name__ == '__main__':
    unittest.main()
