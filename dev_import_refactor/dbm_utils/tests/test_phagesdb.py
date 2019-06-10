""" Unit tests for misc. functions that interact with PhagesDB."""

# from classes import DataGroup
# from classes import Genome
# from classes import Ticket
# from classes import Eval
from functions import phagesdb
from classes import Genome
from constants import constants
import unittest
import urllib.request
import json



class TestPhagesDBFunctions(unittest.TestCase):


    def setUp(self):

        self.API_PREFIX = constants.API_PREFIX
        self.API_SUFFIX = constants.API_SUFFIX


        # self.API_PREFIX = "https://phagesdb.org/api/phages/"
        # self.API_SUFFIX = "/?format=json"

        self.genome = Genome.Genome()


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




    def test_parse_phagesdb_filename_1(self):
        """Verify fasta filename is retrieved and no error is produced."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"fasta_file": url}
        filename, eval_result = phagesdb.parse_phagesdb_filename(data_dict)
        expected_filename = url
        with self.subTest():
            self.assertEqual(filename, expected_filename)
        with self.subTest():
            self.assertIsNone(eval_result)

    def test_parse_phagesdb_filename_2(self):
        """Verify fasta filename is not retrieved and an error is produced."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"fasta_file_x": url}
        filename, eval_result = phagesdb.parse_phagesdb_filename(data_dict)
        expected_filename = ""
        with self.subTest():
            self.assertEqual(filename, expected_filename)
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




    def test_parse_fasta_file_1(self):
        """Verify it fasta file data is parsed correctly."""
        fasta_data = ">Trixie  \nAAAAAAAAAA   \nTTTTTTT \nCCC\nGGGGGGGGGGG\n\n"
        expected_header = "Trixie"
        expected_sequence = "AAAAAAAAAATTTTTTTCCCGGGGGGGGGGG"
        result_list, eval_result = phagesdb.parse_fasta_file(fasta_data)
        with self.subTest():
            self.assertEqual(result_list[0], expected_header)
        with self.subTest():
            self.assertEqual(result_list[1], expected_sequence)
        with self.subTest():
            self.assertIsNone(eval_result)

    def test_parse_fasta_file_2(self):
        """Verify it incorrect fasta file format (no ">") produces an error."""
        fasta_data = "Trixie  \nAAAAAAAAAA   \nTTTTTTT \nCCC\nGGGGGGGGGGG\n\n"
        expected_header = "Trixie"
        expected_sequence = "AAAAAAAAAATTTTTTTCCCGGGGGGGGGGG"
        result_list, eval_result = phagesdb.parse_fasta_file(fasta_data)
        with self.subTest():
            self.assertEqual(result_list[0], expected_header)
        with self.subTest():
            self.assertEqual(result_list[1], expected_sequence)
        with self.subTest():
            self.assertIsNotNone(eval_result)

    def test_parse_fasta_file_3(self):
        """Verify it incorrect fasta file format (no new lines) produces
        an error."""
        fasta_data = "Trixie  AAAAAAAAAA   TTTTTTT CCCGGGGGGGGGGG"
        expected_header = ""
        expected_sequence = ""
        result_list, eval_result = phagesdb.parse_fasta_file(fasta_data)
        with self.subTest():
            self.assertEqual(result_list[0], expected_header)
        with self.subTest():
            self.assertEqual(result_list[1], expected_sequence)
        with self.subTest():
            self.assertIsNotNone(eval_result)




    def test_parse_phagesdb_data_1(self):
        """Verify genome object is parsed from PhagesDB."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        eval_results = phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_search_id = "trixie"
        expected_cluster = "A"
        expected_subcluster = "A2"
        expected_host = "Mycobacterium"
        expected_accession = "ABC123"
        expected_filename = url
        expected_seq_start = "GGTCGGTT"
        expected_seq_end =   "GTCGGTTA"

        with self.subTest():
            self.assertEqual(self.genome.phage_name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.phage_id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.search_id, expected_search_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.filename, expected_filename)
        with self.subTest():
            self.assertEqual(len(self.genome.parsed_record), 2)
        with self.subTest():
            self.assertEqual(self.genome.sequence[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.genome.sequence[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(len(eval_results), 0)


    def test_parse_phagesdb_data_2(self):
        """Verify error is produced from no phage_name key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name_x":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        eval_results = phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = ""
        expected_phage_id = ""
        expected_search_id = ""
        expected_cluster = "A"
        expected_subcluster = "A2"
        expected_host = "Mycobacterium"
        expected_accession = "ABC123"
        expected_filename = url
        expected_seq_start = "GGTCGGTT"
        expected_seq_end =   "GTCGGTTA"

        with self.subTest():
            self.assertEqual(self.genome.phage_name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.phage_id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.search_id, expected_search_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.filename, expected_filename)
        with self.subTest():
            self.assertEqual(len(self.genome.parsed_record), 2)
        with self.subTest():
            self.assertEqual(self.genome.sequence[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.genome.sequence[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(len(eval_results), 1)


    def test_parse_phagesdb_data_3(self):
        """Verify error is produced from no pcluster key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster_x": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        eval_results = phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_search_id = "trixie"
        expected_cluster = ""
        expected_subcluster = "A2"
        expected_host = "Mycobacterium"
        expected_accession = "ABC123"
        expected_filename = url
        expected_seq_start = "GGTCGGTT"
        expected_seq_end =   "GTCGGTTA"

        with self.subTest():
            self.assertEqual(self.genome.phage_name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.phage_id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.search_id, expected_search_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.filename, expected_filename)
        with self.subTest():
            self.assertEqual(len(self.genome.parsed_record), 2)
        with self.subTest():
            self.assertEqual(self.genome.sequence[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.genome.sequence[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(len(eval_results), 1)


    def test_parse_phagesdb_data_4(self):
        """Verify error is produced from no psubcluster key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster_x": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        eval_results = phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_search_id = "trixie"
        expected_cluster = "A"
        expected_subcluster = "none"
        expected_host = "Mycobacterium"
        expected_accession = "ABC123"
        expected_filename = url
        expected_seq_start = "GGTCGGTT"
        expected_seq_end =   "GTCGGTTA"

        with self.subTest():
            self.assertEqual(self.genome.phage_name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.phage_id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.search_id, expected_search_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.filename, expected_filename)
        with self.subTest():
            self.assertEqual(len(self.genome.parsed_record), 2)
        with self.subTest():
            self.assertEqual(self.genome.sequence[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.genome.sequence[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(len(eval_results), 1)


    def test_parse_phagesdb_data_5(self):
        """Verify error is produced from no isolation_host key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host_x": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        eval_results = phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_search_id = "trixie"
        expected_cluster = "A"
        expected_subcluster = "A2"
        expected_host = ""
        expected_accession = "ABC123"
        expected_filename = url
        expected_seq_start = "GGTCGGTT"
        expected_seq_end =   "GTCGGTTA"

        with self.subTest():
            self.assertEqual(self.genome.phage_name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.phage_id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.search_id, expected_search_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.filename, expected_filename)
        with self.subTest():
            self.assertEqual(len(self.genome.parsed_record), 2)
        with self.subTest():
            self.assertEqual(self.genome.sequence[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.genome.sequence[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(len(eval_results), 1)


    def test_parse_phagesdb_data_6(self):
        """Verify error is produced from no genbank_accession key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession_x": "ABC123",
                    "fasta_file": url}
        eval_results = phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_search_id = "trixie"
        expected_cluster = "A"
        expected_subcluster = "A2"
        expected_host = "Mycobacterium"
        expected_accession = "none"
        expected_filename = url
        expected_seq_start = "GGTCGGTT"
        expected_seq_end =   "GTCGGTTA"

        with self.subTest():
            self.assertEqual(self.genome.phage_name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.phage_id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.search_id, expected_search_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.filename, expected_filename)
        with self.subTest():
            self.assertEqual(len(self.genome.parsed_record), 2)
        with self.subTest():
            self.assertEqual(self.genome.sequence[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.genome.sequence[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(len(eval_results), 1)


    def test_parse_phagesdb_data_7(self):
        """Verify error is produced from no fasta_file key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file_x": url}
        eval_results = phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_search_id = "trixie"
        expected_cluster = "A"
        expected_subcluster = "A2"
        expected_host = "Mycobacterium"
        expected_accession = "ABC123"
        expected_filename = ""
        expected_seq = ""

        with self.subTest():
            self.assertEqual(self.genome.phage_name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.phage_id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.search_id, expected_search_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.filename, expected_filename)
        with self.subTest():
            self.assertEqual(len(self.genome.parsed_record), 0)
        with self.subTest():
            self.assertEqual(self.genome.sequence, expected_seq)
        with self.subTest():
            self.assertEqual(len(eval_results), 1)


    def test_parse_phagesdb_data_8(self):
        """Verify error is produced from incorrect fasta_file URL."""

        url = "https://phagesdb.org/media/fastas/L5_x.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        eval_results = phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_search_id = "trixie"
        expected_cluster = "A"
        expected_subcluster = "A2"
        expected_host = "Mycobacterium"
        expected_accession = "ABC123"
        expected_filename = url
        expected_seq = ""

        with self.subTest():
            self.assertEqual(self.genome.phage_name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.phage_id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.search_id, expected_search_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.filename, expected_filename)
        with self.subTest():
            self.assertEqual(len(self.genome.parsed_record), 0)
        with self.subTest():
            self.assertEqual(self.genome.sequence, expected_seq)
        with self.subTest():
            self.assertEqual(len(eval_results), 1)


    def test_parse_phagesdb_data_9(self):
        """Verify that errors accumulate."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name_x":"Trixie",
                    "pcluster_x": {"cluster": "A"},
                    "psubcluster_x": {"subcluster": "A2"},
                    "isolation_host_x": {"genus": "Mycobacterium"},
                    "genbank_accession_x": "ABC123",
                    "fasta_file_x": url}
        eval_results = phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = ""
        expected_phage_id = ""
        expected_search_id = ""
        expected_cluster = ""
        expected_subcluster = "none"
        expected_host = ""
        expected_accession = "none"
        expected_filename = ""
        expected_seq = ""

        with self.subTest():
            self.assertEqual(self.genome.phage_name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.phage_id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.search_id, expected_search_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.filename, expected_filename)
        with self.subTest():
            self.assertEqual(len(self.genome.parsed_record), 0)
        with self.subTest():
            self.assertEqual(self.genome.sequence, expected_seq)
        with self.subTest():
            self.assertEqual(len(eval_results), 6)




    def test_retrieve_phagesdb_data_1(self):
        """Verify data is retrieved from PhagesDB with no error produced."""

        url = self.API_PREFIX + "L5" + self.API_SUFFIX
        data_dict, eval_result = phagesdb.retrieve_phagesdb_data(url)
        expected_phage_name = "L5"
        with self.subTest():
            self.assertEqual(data_dict["phage_name"], expected_phage_name)
        with self.subTest():
            self.assertIsNone(eval_result)

    def test_retrieve_phagesdb_data_2(self):
        """Verify data is not retrieved from PhagesDB and an error produced."""

        url = self.API_PREFIX + "L5_x" + self.API_SUFFIX
        data_dict, eval_result = phagesdb.retrieve_phagesdb_data(url)
        with self.subTest():
            self.assertEqual(len(data_dict.keys()), 0)
        with self.subTest():
            self.assertIsNotNone(eval_result)




    def test_construct_phage_url_1(self):
        """Verify URL is constructed correctly."""

        phage_name = "L5"
        expected_url = self.API_PREFIX + "L5" + self.API_SUFFIX
        url = phagesdb.construct_phage_url(phage_name)
        self.assertEqual(url, expected_url)






if __name__ == '__main__':
    unittest.main()
