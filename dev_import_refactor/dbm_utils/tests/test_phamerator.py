""" Unit tests for misc. functions that interact with PhameratorDB."""


import unittest
from functions import phamerator
from classes import Genome
from datetime import datetime


class TestTicketFunctions1(unittest.TestCase):


    def setUp(self):
        self.genome = Genome.Genome()







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

        phamerator.parse_phamerator_data(self.genome, data_tuple)

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
            self.assertEqual(self.genome.phage_id, output_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.phage_name, output_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.host, output_host)
        with self.subTest():
            self.assertEqual(self.genome.sequence, output_sequence)
        with self.subTest():
            self.assertEqual(self.genome.status, output_status)
        with self.subTest():
            self.assertEqual(self.genome.cluster, output_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, output_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.date_last_modified, \
                output_date_last_modified)
        with self.subTest():
            self.assertEqual(self.genome.accession, output_accession)
        with self.subTest():
            self.assertEqual(self.genome.annotation_author, \
                output_annotation_author)
        with self.subTest():
            self.assertEqual(self.genome.annotation_qc, output_annotation_qc)
        with self.subTest():
            self.assertEqual(self.genome.retrieve_record, \
                output_retrieve_record)
        with self.subTest():
            self.assertEqual(self.genome.search_id, output_search_id)
        with self.subTest():
            self.assertEqual(self.genome._length, output_seq_length)


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

        phamerator.parse_phamerator_data(self.genome, data_tuple)

        output_phage_id = "Trixie"
        output_phage_name = "Trixie_Draft"
        output_host = "none"
        output_sequence = "ATCG"
        output_status = "final"
        output_cluster = "singleton"
        output_subcluster = "none"
        output_date_last_modified = datetime.strptime('1/1/1900', '%m/%d/%Y')
        output_accession = "none"
        output_annotation_author = "1"
        output_annotation_qc = "1"
        output_retrieve_record = "1"

        output_search_id = "trixie"
        output_seq_length = 4

        with self.subTest():
            self.assertEqual(self.genome.phage_id, output_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.phage_name, output_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.host, output_host)
        with self.subTest():
            self.assertEqual(self.genome.sequence, output_sequence)
        with self.subTest():
            self.assertEqual(self.genome.status, output_status)
        with self.subTest():
            self.assertEqual(self.genome.cluster, output_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, output_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.date_last_modified, \
                output_date_last_modified)
        with self.subTest():
            self.assertEqual(self.genome.accession, output_accession)
        with self.subTest():
            self.assertEqual(self.genome.annotation_author, \
                output_annotation_author)
        with self.subTest():
            self.assertEqual(self.genome.annotation_qc, output_annotation_qc)
        with self.subTest():
            self.assertEqual(self.genome.retrieve_record, \
                output_retrieve_record)
        with self.subTest():
            self.assertEqual(self.genome.search_id, output_search_id)
        with self.subTest():
            self.assertEqual(self.genome._length, output_seq_length)














if __name__ == '__main__':
    unittest.main()
