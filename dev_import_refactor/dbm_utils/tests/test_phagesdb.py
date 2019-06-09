""" Unit tests for misc. functions that interact with PhagesDB."""

# from classes import DataGroup
# from classes import Genome
# from classes import Ticket
# from classes import Eval
from functions import phagesdb
import unittest
import urllib.request
import json



class TestPhagesDBFunctions(unittest.TestCase):


    def setUp(self):
        self.API_PREFIX = "https://phagesdb.org/api/phages/"
        self.API_SUFFIX = "/?format=json"




    def test_parse_phagesdb_phage_name_1(self):
        """Verify name is retrieved and no error is produced."""

        data_dict = {"phage_name":"Trixie"}

        phage, eval_result = phagesdb.parse_phagesdb_phage_name(data_dict)
        expected_phage = "Trixie"
        with self.subTest():
            self.assertEqual(phage, expected_phage)
        with self.subTest():
            self.assertIsNone(eval_result)

    def test_parse_phagesdb_phage_name_2(self):
        """Verify name is not retrieved and error is produced."""

        data_dict = {"phage_id":"Trixie"}
        phage, eval_result = phagesdb.parse_phagesdb_phage_name(data_dict)
        expected_phage = ""
        with self.subTest():
            self.assertEqual(phage, expected_phage)
        with self.subTest():
            self.assertIsNotNone(eval_result)




    def test_parse_phagesdb_cluster_1(self):
        """Verify missing cluster is retrieved and no error is produced."""

        data_dict = {"pcluster": None}
        cluster, eval_result = phagesdb.parse_phagesdb_cluster(data_dict)
        expected_cluster = "UNK"
        with self.subTest():
            self.assertEqual(cluster, expected_cluster)
        with self.subTest():
            self.assertIsNone(eval_result)

    def test_parse_phagesdb_cluster_2(self):
        """Verify no cluster is retrieved and error is produced."""

        data_dict = {"pcluster_x": None}
        cluster, eval_result = phagesdb.parse_phagesdb_cluster(data_dict)
        expected_cluster = ""
        with self.subTest():
            self.assertEqual(cluster, expected_cluster)
        with self.subTest():
            self.assertIsNotNone(eval_result)

    def test_parse_phagesdb_cluster_3(self):
        """Verify standard cluster is retrieved and no error is produced."""

        data_dict = {"pcluster": {"cluster": "A"}}
        cluster, eval_result = phagesdb.parse_phagesdb_cluster(data_dict)
        expected_cluster = "A"
        with self.subTest():
            self.assertEqual(cluster, expected_cluster)
        with self.subTest():
            self.assertIsNone(eval_result)




    def test_parse_phagesdb_subcluster_1(self):
        """Verify missing subcluster is retrieved and no error is produced."""

        data_dict = {"psubcluster": None}
        subcluster, eval_result = phagesdb.parse_phagesdb_subcluster(data_dict)
        expected_subcluster = "none"
        with self.subTest():
            self.assertEqual(subcluster, expected_subcluster)
        with self.subTest():
            self.assertIsNone(eval_result)

    def test_parse_phagesdb_subcluster_2(self):
        """Verify no subcluster is retrieved and an error is produced."""

        data_dict = {"psubcluster_x": None}
        subcluster, eval_result = phagesdb.parse_phagesdb_subcluster(data_dict)
        expected_subcluster = ""
        with self.subTest():
            self.assertEqual(subcluster, expected_subcluster)
        with self.subTest():
            self.assertIsNotNone(eval_result)

    def test_parse_phagesdb_subcluster_3(self):
        """Verify standard subcluster is retrieved and no error is produced."""

        data_dict = {"psubcluster": {"subcluster": "A2"}}
        subcluster, eval_result = phagesdb.parse_phagesdb_subcluster(data_dict)
        expected_subcluster = "A2"
        with self.subTest():
            self.assertEqual(subcluster, expected_subcluster)
        with self.subTest():
            self.assertIsNone(eval_result)




    def test_parse_phagesdb_host_1(self):
        """Verify host genus is retrieved and no error is produced."""

        data_dict = {"isolation_host": {"genus": "Mycobacterium"}}
        host, eval_result = phagesdb.parse_phagesdb_host(data_dict)
        expected_host = "Mycobacterium"
        with self.subTest():
            self.assertEqual(host, expected_host)
        with self.subTest():
            self.assertIsNone(eval_result)

    def test_parse_phagesdb_host_2(self):
        """Verify host genus is not retrieved and an error is produced."""

        data_dict = {"isolation_host_x": {"genus": "Mycobacterium"}}
        host, eval_result = phagesdb.parse_phagesdb_host(data_dict)
        expected_host = ""
        with self.subTest():
            self.assertEqual(host, expected_host)
        with self.subTest():
            self.assertIsNotNone(eval_result)




    def test_parse_phagesdb_accession_1(self):
        """Verify accession is retrieved and no error is produced."""

        data_dict = {"genbank_accession": "ABC123"}
        accession, eval_result = phagesdb.parse_phagesdb_accession(data_dict)
        expected_accession = "ABC123"
        with self.subTest():
            self.assertEqual(accession, expected_accession)
        with self.subTest():
            self.assertIsNone(eval_result)

    def test_parse_phagesdb_accession_2(self):
        """Verify accession is not retrieved and an error is produced."""

        data_dict = {"genbank_accession_x": "ABC123"}
        accession, eval_result = phagesdb.parse_phagesdb_accession(data_dict)
        expected_accession = ""
        with self.subTest():
            self.assertEqual(accession, expected_accession)
        with self.subTest():
            self.assertIsNotNone(eval_result)




    def test_retrieve_phagesdb_fasta_1(self):
        """Verify fasta data is retrieved and no error is produced."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        fasta_data, eval_result = phagesdb.retrieve_phagesdb_fasta(url)
        expected_fasta_data_header = ">Mycobacterium phage L5"
        with self.subTest():
            self.assertEqual(fasta_data[:23], expected_fasta_data_header)
        with self.subTest():
            self.assertIsNone(eval_result)

    def test_retrieve_phagesdb_fasta_2(self):
        """Verify fasta data is not retrieved and an error is produced."""

        url = "https://phagesdb.org/media/fastas/L5_x.fasta"
        fasta_data, eval_result = phagesdb.retrieve_phagesdb_fasta(url)
        expected_fasta_data_header = ""
        with self.subTest():
            self.assertEqual(fasta_data, expected_fasta_data_header)
        with self.subTest():
            self.assertIsNotNone(eval_result)




if __name__ == '__main__':
    unittest.main()
