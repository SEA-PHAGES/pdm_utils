""" Unit tests for misc. functions that interact with PhagesDB."""

from classes import Bundle
from functions import phagesdb
from classes import Genome
from constants import constants
import unittest
import urllib.request
import json
from Bio.SeqRecord import SeqRecord



class TestPhagesDBFunctions(unittest.TestCase):


    def setUp(self):

        self.API_PREFIX = constants.API_PREFIX
        self.API_SUFFIX = constants.API_SUFFIX

        self.genome = Genome.Genome()


    def test_parse_phagesdb_phage_name_1(self):
        """Verify name is retrieved and no error is produced."""

        data_dict = {"phage_name":"Trixie"}

        phage = phagesdb.parse_phagesdb_phage_name(data_dict)
        expected_phage = "Trixie"
        self.assertEqual(phage, expected_phage)

    def test_parse_phagesdb_phage_name_2(self):
        """Verify name is not retrieved and error is produced."""

        data_dict = {"phage_id":"Trixie"}
        phage = phagesdb.parse_phagesdb_phage_name(data_dict)
        expected_phage = ""
        self.assertEqual(phage, expected_phage)




    def test_parse_phagesdb_cluster_1(self):
        """Verify missing cluster is retrieved and no error is produced."""

        data_dict = {"pcluster": None}
        cluster = phagesdb.parse_phagesdb_cluster(data_dict)
        expected_cluster = "UNK"
        self.assertEqual(cluster, expected_cluster)

    def test_parse_phagesdb_cluster_2(self):
        """Verify no cluster is retrieved and error is produced."""

        data_dict = {"pcluster_x": None}
        cluster = phagesdb.parse_phagesdb_cluster(data_dict)
        expected_cluster = ""
        self.assertEqual(cluster, expected_cluster)

    def test_parse_phagesdb_cluster_3(self):
        """Verify standard cluster is retrieved and no error is produced."""

        data_dict = {"pcluster": {"cluster": "A"}}
        cluster = phagesdb.parse_phagesdb_cluster(data_dict)
        expected_cluster = "A"
        self.assertEqual(cluster, expected_cluster)




    def test_parse_phagesdb_subcluster_1(self):
        """Verify missing subcluster is retrieved and no error is produced."""

        data_dict = {"psubcluster": None}
        subcluster = phagesdb.parse_phagesdb_subcluster(data_dict)
        expected_subcluster = "none"
        self.assertEqual(subcluster, expected_subcluster)

    def test_parse_phagesdb_subcluster_2(self):
        """Verify no subcluster is retrieved and an error is produced."""

        data_dict = {"psubcluster_x": None}
        subcluster = phagesdb.parse_phagesdb_subcluster(data_dict)
        expected_subcluster = ""
        self.assertEqual(subcluster, expected_subcluster)

    def test_parse_phagesdb_subcluster_3(self):
        """Verify standard subcluster is retrieved and no error is produced."""

        data_dict = {"psubcluster": {"subcluster": "A2"}}
        subcluster = phagesdb.parse_phagesdb_subcluster(data_dict)
        expected_subcluster = "A2"
        self.assertEqual(subcluster, expected_subcluster)




    def test_parse_phagesdb_host_genus_1(self):
        """Verify host genus is retrieved and no error is produced."""

        data_dict = {"isolation_host": {"genus": "Mycobacterium"}}
        host = phagesdb.parse_phagesdb_host_genus(data_dict)
        expected_host = "Mycobacterium"
        self.assertEqual(host, expected_host)

    def test_parse_phagesdb_host_genus_2(self):
        """Verify host genus is not retrieved and an error is produced."""

        data_dict = {"isolation_host_x": {"genus": "Mycobacterium"}}
        host = phagesdb.parse_phagesdb_host_genus(data_dict)
        expected_host = ""
        self.assertEqual(host, expected_host)




    def test_parse_phagesdb_accession_1(self):
        """Verify accession is retrieved and no error is produced."""

        data_dict = {"genbank_accession": "ABC123"}
        accession = phagesdb.parse_phagesdb_accession(data_dict)
        expected_accession = "ABC123"
        self.assertEqual(accession, expected_accession)

    def test_parse_phagesdb_accession_2(self):
        """Verify accession is not retrieved and an error is produced."""

        data_dict = {"genbank_accession_x": "ABC123"}
        accession = phagesdb.parse_phagesdb_accession(data_dict)
        expected_accession = ""
        self.assertEqual(accession, expected_accession)




    def test_parse_phagesdb_filename_1(self):
        """Verify fasta filename is retrieved and no error is produced."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"fasta_file": url}
        filename = phagesdb.parse_phagesdb_filename(data_dict)
        expected_filename = url
        self.assertEqual(filename, expected_filename)

    def test_parse_phagesdb_filename_2(self):
        """Verify fasta filename is not retrieved and an error is produced."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"fasta_file_x": url}
        filename = phagesdb.parse_phagesdb_filename(data_dict)
        expected_filename = ""
        self.assertEqual(filename, expected_filename)




    def test_retrieve_phagesdb_fasta_1(self):
        """Verify fasta data is retrieved and no error is produced."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        fasta_data = phagesdb.retrieve_phagesdb_fasta(url)
        expected_fasta_data_header = ">Mycobacterium phage L5"
        self.assertEqual(fasta_data[:23], expected_fasta_data_header)

    def test_retrieve_phagesdb_fasta_2(self):
        """Verify fasta data is not retrieved and an error is produced."""

        url = "https://phagesdb.org/media/fastas/L5_x.fasta"
        fasta_data = phagesdb.retrieve_phagesdb_fasta(url)
        expected_fasta_data_header = ""
        self.assertEqual(fasta_data, expected_fasta_data_header)




    def test_parse_fasta_file_1(self):
        """Verify it fasta file data is parsed correctly."""
        fasta_data = ">Trixie  \nAAAAAAAAAA   \nTTTTTTT \nCCC\nGGGGGGGGGGG\n\n"
        expected_header = "Trixie"
        expected_sequence = "AAAAAAAAAATTTTTTTCCCGGGGGGGGGGG"
        result_tuple = phagesdb.parse_fasta_file(fasta_data)
        with self.subTest():
            self.assertEqual(result_tuple[0], expected_header)
        with self.subTest():
            self.assertEqual(result_tuple[1], expected_sequence)

    def test_parse_fasta_file_2(self):
        """Verify it incorrect fasta file format (no ">") produces an error."""
        fasta_data = "Trixie  \nAAAAAAAAAA   \nTTTTTTT \nCCC\nGGGGGGGGGGG\n\n"
        expected_header = "Trixie"
        expected_sequence = "AAAAAAAAAATTTTTTTCCCGGGGGGGGGGG"
        result_tuple = phagesdb.parse_fasta_file(fasta_data)
        with self.subTest():
            self.assertEqual(result_tuple[0], expected_header)
        with self.subTest():
            self.assertEqual(result_tuple[1], expected_sequence)

    def test_parse_fasta_file_3(self):
        """Verify it incorrect fasta file format (no new lines) produces
        an error."""
        fasta_data = "Trixie  AAAAAAAAAA   TTTTTTT CCCGGGGGGGGGGG"
        expected_header = ""
        expected_sequence = ""
        result_tuple = phagesdb.parse_fasta_file(fasta_data)
        with self.subTest():
            self.assertEqual(result_tuple[0], expected_header)
        with self.subTest():
            self.assertEqual(result_tuple[1], expected_sequence)




    def test_parse_phagesdb_data_1(self):
        """Verify genome object is parsed from PhagesDB."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        phagesdb.parse_phagesdb_data(self.genome, data_dict)


        description = "Mycobacterium phage L5"
        errors = 0
        for eval in self.genome.evaluations:
            if eval.status == "error":
                errors += 1

        with self.subTest():
            self.assertEqual(self.genome.name, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome.id, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome.cluster, "A")
        with self.subTest():
            self.assertEqual(self.genome.subcluster, "A2")
        with self.subTest():
            self.assertEqual(self.genome.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.genome.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.genome.record_filename, url)
        with self.subTest():
            self.assertTrue(isinstance(self.genome.record, SeqRecord))
        with self.subTest():
            self.assertEqual(self.genome.sequence[:8], "GGTCGGTT")
        with self.subTest():
            self.assertEqual(self.genome.sequence[-8:], "GTCGGTTA")
        with self.subTest():
            self.assertEqual(self.genome.record_description, description)
        with self.subTest():
            self.assertEqual(self.genome.type, "phagesdb")
        with self.subTest():
            self.assertEqual(errors, 0)


    def test_parse_phagesdb_data_2(self):
        """Verify error is produced from no phage_name key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name_x":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = ""
        expected_phage_id = ""
        expected_cluster = "A"
        expected_subcluster = "A2"
        expected_host = "Mycobacterium"
        expected_accession = "ABC123"
        expected_filename = url
        expected_seq_start = "GGTCGGTT"
        expected_seq_end =   "GTCGGTTA"
        expected_type = "phagesdb"

        errors = 0
        for eval in self.genome.evaluations:
            if eval.status == "error":
                errors += 1

        with self.subTest():
            self.assertEqual(self.genome.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.record_filename, expected_filename)
        with self.subTest():
            self.assertTrue(isinstance(self.genome.record, SeqRecord))
        with self.subTest():
            self.assertEqual(self.genome.sequence[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.genome.sequence[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(self.genome.type, expected_type)
        with self.subTest():
            self.assertEqual(errors, 2)


    def test_parse_phagesdb_data_3(self):
        """Verify error is produced from no pcluster key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster_x": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_cluster = ""
        expected_subcluster = "A2"
        expected_host = "Mycobacterium"
        expected_accession = "ABC123"
        expected_filename = url
        expected_seq_start = "GGTCGGTT"
        expected_seq_end =   "GTCGGTTA"
        expected_type = "phagesdb"

        errors = 0
        for eval in self.genome.evaluations:
            if eval.status == "error":
                errors += 1

        with self.subTest():
            self.assertEqual(self.genome.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.record_filename, expected_filename)
        with self.subTest():
            self.assertTrue(isinstance(self.genome.record, SeqRecord))
        with self.subTest():
            self.assertEqual(self.genome.sequence[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.genome.sequence[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(self.genome.type, expected_type)
        with self.subTest():
            self.assertEqual(errors, 1)


    def test_parse_phagesdb_data_4(self):
        """Verify error is produced from no psubcluster key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster_x": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_cluster = "A"
        expected_subcluster = ""
        expected_host = "Mycobacterium"
        expected_accession = "ABC123"
        expected_filename = url
        expected_seq_start = "GGTCGGTT"
        expected_seq_end =   "GTCGGTTA"
        expected_type = "phagesdb"

        errors = 0
        for eval in self.genome.evaluations:
            if eval.status == "error":
                errors += 1

        with self.subTest():
            self.assertEqual(self.genome.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.record_filename, expected_filename)
        with self.subTest():
            self.assertTrue(isinstance(self.genome.record, SeqRecord))
        with self.subTest():
            self.assertEqual(self.genome.sequence[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.genome.sequence[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(self.genome.type, expected_type)
        with self.subTest():
            self.assertEqual(errors, 1)


    def test_parse_phagesdb_data_5(self):
        """Verify error is produced from no isolation_host key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host_x": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_cluster = "A"
        expected_subcluster = "A2"
        expected_host = ""
        expected_accession = "ABC123"
        expected_filename = url
        expected_seq_start = "GGTCGGTT"
        expected_seq_end =   "GTCGGTTA"
        expected_type = "phagesdb"

        errors = 0
        for eval in self.genome.evaluations:
            if eval.status == "error":
                errors += 1

        with self.subTest():
            self.assertEqual(self.genome.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.record_filename, expected_filename)
        with self.subTest():
            self.assertTrue(isinstance(self.genome.record, SeqRecord))
        with self.subTest():
            self.assertEqual(self.genome.sequence[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.genome.sequence[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(self.genome.type, expected_type)
        with self.subTest():
            self.assertEqual(errors, 1)


    def test_parse_phagesdb_data_6(self):
        """Verify error is produced from no genbank_accession key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession_x": "ABC123",
                    "fasta_file": url}
        phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_cluster = "A"
        expected_subcluster = "A2"
        expected_host = "Mycobacterium"
        expected_accession = ""
        expected_filename = url
        expected_seq_start = "GGTCGGTT"
        expected_seq_end =   "GTCGGTTA"
        expected_type = "phagesdb"

        errors = 0
        for eval in self.genome.evaluations:
            if eval.status == "error":
                errors += 1


        with self.subTest():
            self.assertEqual(self.genome.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.record_filename, expected_filename)
        with self.subTest():
            self.assertTrue(isinstance(self.genome.record, SeqRecord))
        with self.subTest():
            self.assertEqual(self.genome.sequence[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.genome.sequence[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(self.genome.type, expected_type)
        with self.subTest():
            self.assertEqual(errors, 1)


    def test_parse_phagesdb_data_7(self):
        """Verify error is produced from no fasta_file key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file_x": url}
        phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_cluster = "A"
        expected_subcluster = "A2"
        expected_host = "Mycobacterium"
        expected_accession = "ABC123"
        expected_filename = ""
        expected_seq = ""
        expected_type = "phagesdb"

        errors = 0
        for eval in self.genome.evaluations:
            if eval.status == "error":
                errors += 1

        with self.subTest():
            self.assertEqual(self.genome.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.record_filename, expected_filename)
        with self.subTest():
            self.assertEqual(self.genome.record, "")
        with self.subTest():
            self.assertEqual(self.genome.sequence, expected_seq)
        with self.subTest():
            self.assertEqual(self.genome.type, expected_type)
        with self.subTest():
            self.assertEqual(errors, 2)


    def test_parse_phagesdb_data_8(self):
        """Verify error is produced from incorrect fasta_file URL."""

        url = "https://phagesdb.org/media/fastas/L5_x.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_cluster = "A"
        expected_subcluster = "A2"
        expected_host = "Mycobacterium"
        expected_accession = "ABC123"
        expected_filename = url
        expected_seq = ""
        expected_type = "phagesdb"

        errors = 0
        for eval in self.genome.evaluations:
            if eval.status == "error":
                errors += 1

        with self.subTest():
            self.assertEqual(self.genome.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.record_filename, expected_filename)
        with self.subTest():
            self.assertEqual(self.genome.record, "")
        with self.subTest():
            self.assertEqual(self.genome.sequence, expected_seq)
        with self.subTest():
            self.assertEqual(self.genome.type, expected_type)
        with self.subTest():
            self.assertEqual(errors, 1)


    def test_parse_phagesdb_data_9(self):
        """Verify that errors accumulate."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name_x":"Trixie",
                    "pcluster_x": {"cluster": "A"},
                    "psubcluster_x": {"subcluster": "A2"},
                    "isolation_host_x": {"genus": "Mycobacterium"},
                    "genbank_accession_x": "ABC123",
                    "fasta_file_x": url}
        phagesdb.parse_phagesdb_data(self.genome, data_dict)

        expected_phage_name = ""
        expected_phage_id = ""
        expected_cluster = ""
        expected_subcluster = ""
        expected_host = ""
        expected_accession = ""
        expected_filename = ""
        expected_seq = ""
        expected_type = "phagesdb"

        errors = 0
        for eval in self.genome.evaluations:
            if eval.status == "error":
                errors += 1

        with self.subTest():
            self.assertEqual(self.genome.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.genome.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.genome.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.genome.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.genome.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.genome.record_filename, expected_filename)
        with self.subTest():
            self.assertEqual(self.genome.record, "")
        with self.subTest():
            self.assertEqual(self.genome.sequence, expected_seq)
        with self.subTest():
            self.assertEqual(self.genome.type, expected_type)
        with self.subTest():
            self.assertEqual(errors, 8)




    def test_retrieve_phagesdb_data_1(self):
        """Verify data is retrieved from PhagesDB with no error produced."""

        url = self.API_PREFIX + "L5" + self.API_SUFFIX
        data_dict = phagesdb.retrieve_phagesdb_data(url)
        expected_phage_name = "L5"
        self.assertEqual(data_dict["phage_name"], expected_phage_name)

    def test_retrieve_phagesdb_data_2(self):
        """Verify data is not retrieved from PhagesDB and an error produced."""

        url = self.API_PREFIX + "L5_x" + self.API_SUFFIX
        data_dict = phagesdb.retrieve_phagesdb_data(url)
        self.assertEqual(len(data_dict.keys()), 0)




    def test_construct_phage_url_1(self):
        """Verify URL is constructed correctly."""

        phage_name = "L5"
        expected_url = self.API_PREFIX + "L5" + self.API_SUFFIX
        url = phagesdb.construct_phage_url(phage_name)
        self.assertEqual(url, expected_url)





    def test_retrieve_phagesdb_data_list_1(self):
        """Confirm that data is successfully retrieved."""
        url = constants.API_CLUSTERS
        data_list = phagesdb.retrieve_phagesdb_data_list(url)
        self.assertTrue(len(data_list) > 0)

    def test_retrieve_phagesdb_data_list_2(self):
        """Confirm that data is not successfully retrieved."""
        url = 'invalid url'
        data_list = phagesdb.retrieve_phagesdb_data_list(url)
        self.assertTrue(len(data_list) == 0)




    def test_create_host_genus_set_1(self):
        """Confirm that host genera data is successfully retrieved."""
        host_genera_set = phagesdb.create_host_genus_set()
        with self.subTest():
            self.assertTrue(len(host_genera_set) > 0)
        with self.subTest():
            self.assertTrue("Mycobacterium" in host_genera_set)


    def test_create_host_genus_set_2(self):
        """Confirm that host genera data are not retrieved due to
        invalid url."""
        host_genera_set = phagesdb.create_host_genus_set("invalid_url")
        self.assertTrue(len(host_genera_set) == 0)




    def test_create_cluster_subcluster_sets_1(self):
        """Confirm that cluster and subcluster data are successfully
        retrieved."""
        cluster_set, \
        subcluster_set = phagesdb.create_cluster_subcluster_sets()
        with self.subTest():
            self.assertTrue("A" in cluster_set)
        with self.subTest():
            self.assertTrue("A1" in subcluster_set)

    def test_create_cluster_subcluster_sets_2(self):
        """Confirm that cluster and subcluster data are not retrieved due
        to invalid url."""
        cluster_set, \
        subcluster_set = phagesdb.create_cluster_subcluster_sets("invalid_url")
        with self.subTest():
            self.assertTrue(len(cluster_set) == 0)
        with self.subTest():
            self.assertTrue(len(subcluster_set) == 0)









class TestPhagesDBFunctions2(unittest.TestCase):

    def setUp(self):


        self.API_PREFIX = constants.API_PREFIX
        self.API_SUFFIX = constants.API_SUFFIX


        self.genome1 = Genome.Genome()
        self.genome1.id = "L5"
        self.genome1.type = "add"
        self.genome1.host_genus = "Gordonia"
        self.genome1.cluster = "B"
        self.genome1._empty_fields = True

        self.bundle1 = Bundle.Bundle()




    def test_copy_data_from_phagesdb_1(self):
        """Check that an "add" genome with no fields set to 'retrieve' is
        not impacted."""

        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        phagesdb.copy_data_from_phagesdb(self.bundle1, "add")
        genome1 = self.bundle1.genome_dict["add"]
        with self.subTest():
            self.assertFalse(genome1._empty_fields)
        with self.subTest():
            self.assertEqual(genome1.host_genus, "Gordonia")
        with self.subTest():
            self.assertEqual(genome1.cluster, "B")
        with self.subTest():
            self.assertEqual(genome1.evaluations[0].status, "correct")

    def test_copy_data_from_phagesdb_2(self):
        """Check that an "add" genome with host_genus field set to 'retrieve' is
        populated correctly."""

        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.genome1.host_genus = "retrieve"
        phagesdb.copy_data_from_phagesdb(self.bundle1, "add")
        genome1 = self.bundle1.genome_dict["add"]
        with self.subTest():
            self.assertFalse(genome1._empty_fields)
        with self.subTest():
            self.assertEqual(genome1.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(genome1.cluster, "B")
        with self.subTest():
            self.assertEqual(genome1.evaluations[0].status, "correct")

    def test_copy_data_from_phagesdb_3(self):
        """Check that an "invalid" genome with host_genus field set to 'retrieve' is
        not populated correctly."""

        self.genome1.type = "invalid"
        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.genome1.host_genus = "retrieve"
        phagesdb.copy_data_from_phagesdb(self.bundle1, "add")
        with self.subTest():
            self.assertEqual(
                len(self.bundle1.genome_pair_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(self.genome1.host_genus, "retrieve")
        with self.subTest():
            self.assertEqual(len(self.genome1.evaluations), 0)

    def test_copy_data_from_phagesdb_4(self):
        """Check that an "add" genome with host_genus field set to 'retrieve' is
        not populated correctly when "invalid" type is requrested."""

        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.genome1.host_genus = "retrieve"
        phagesdb.copy_data_from_phagesdb(self.bundle1, "invalid")
        with self.subTest():
            self.assertEqual(
                len(self.bundle1.genome_pair_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(self.genome1.host_genus, "retrieve")
        with self.subTest():
            self.assertEqual(len(self.genome1.evaluations), 0)

    def test_copy_data_from_phagesdb_5(self):
        """Check that an "add" genome with host_genus field set to 'retrieve' is
        not populated correctly when id is not valid."""

        self.genome1.id = "invalid"
        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.genome1.host_genus = "retrieve"
        self.genome1._empty_fields = False
        phagesdb.copy_data_from_phagesdb(self.bundle1, "add")
        with self.subTest():
            self.assertTrue(self.genome1._empty_fields)
        with self.subTest():
            self.assertEqual(
                len(self.bundle1.genome_pair_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(self.genome1.host_genus, "retrieve")
        with self.subTest():
            self.assertEqual(self.genome1.evaluations[0].status, "error")

















if __name__ == '__main__':
    unittest.main()
