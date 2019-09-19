"""Integration tests for misc. functions that interact with PhagesDB."""

from pdm_utils.classes import bundle
from pdm_utils.functions import phagesdb
from pdm_utils.classes import genome
from pdm_utils.constants import constants
import unittest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq




class TestPhagesDBFunctions(unittest.TestCase):


    def setUp(self):

        self.API_PREFIX = constants.API_PREFIX
        self.API_SUFFIX = constants.API_SUFFIX

        self.gnm = genome.Genome()




    def test_retrieve_fasta_data_1(self):
        """Verify fasta data is retrieved and no error is produced."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        fasta_data = phagesdb.retrieve_fasta_data(url)
        expected_fasta_data_header = ">Mycobacterium phage L5"
        self.assertEqual(fasta_data[:23], expected_fasta_data_header)

    def test_retrieve_fasta_data_2(self):
        """Verify fasta data is not retrieved and an error is produced."""

        url = "https://phagesdb.org/media/fastas/L5_x.fasta"
        fasta_data = phagesdb.retrieve_fasta_data(url)
        expected_fasta_data_header = ""
        self.assertEqual(fasta_data, expected_fasta_data_header)




    def test_parse_genome_data_1(self):
        """Verify genome object is parsed from PhagesDB."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict)


        description = "Mycobacterium phage L5"
        errors = 0
        for evl in self.gnm.evaluations:
            if evl.status == "error":
                errors += 1

        with self.subTest():
            self.assertEqual(self.gnm.name, "Trixie")
        with self.subTest():
            self.assertEqual(self.gnm.id, "Trixie")
        with self.subTest():
            self.assertEqual(self.gnm.cluster, "A")
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, "A2")
        with self.subTest():
            self.assertEqual(self.gnm.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.gnm.accession, "ABC123")
        with self.subTest():
            self.assertEqual(self.gnm.filename, url)
        with self.subTest():
            self.assertEqual(self.gnm.seq[:8], "GGTCGGTT")
        with self.subTest():
            self.assertEqual(self.gnm.seq[-8:], "GTCGGTTA")
        with self.subTest():
            self.assertIsInstance(self.gnm.seq, Seq)
        with self.subTest():
            self.assertEqual(self.gnm.description, description)
        with self.subTest():
            self.assertEqual(self.gnm.type, "phagesdb")
        with self.subTest():
            self.assertEqual(errors, 0)


    def test_parse_genome_data_2(self):
        """Verify output when there is no phage_name key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name_x":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict)

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

        with self.subTest():
            self.assertEqual(self.gnm.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.gnm.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.gnm.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.gnm.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.gnm.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.gnm.filename, expected_filename)
        with self.subTest():
            self.assertEqual(self.gnm.seq[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.gnm.seq[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(self.gnm.type, expected_type)


    def test_parse_genome_data_3(self):
        """Verify output when there is no pcluster key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster_x": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict)

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

        with self.subTest():
            self.assertEqual(self.gnm.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.gnm.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.gnm.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.gnm.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.gnm.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.gnm.filename, expected_filename)
        with self.subTest():
            self.assertEqual(self.gnm.seq[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.gnm.seq[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(self.gnm.type, expected_type)


    def test_parse_genome_data_4(self):
        """Verify output when there is no psubcluster key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster_x": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict)

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

        with self.subTest():
            self.assertEqual(self.gnm.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.gnm.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.gnm.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.gnm.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.gnm.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.gnm.filename, expected_filename)
        with self.subTest():
            self.assertEqual(self.gnm.seq[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.gnm.seq[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(self.gnm.type, expected_type)


    def test_parse_genome_data_5(self):
        """Verify output when there is no isolation_host key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host_x": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict)

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

        with self.subTest():
            self.assertEqual(self.gnm.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.gnm.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.gnm.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.gnm.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.gnm.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.gnm.filename, expected_filename)
        with self.subTest():
            self.assertEqual(self.gnm.seq[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.gnm.seq[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(self.gnm.type, expected_type)


    def test_parse_genome_data_6(self):
        """Verify output when there is no genbank_accession key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession_x": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict)

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

        with self.subTest():
            self.assertEqual(self.gnm.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.gnm.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.gnm.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.gnm.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.gnm.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.gnm.filename, expected_filename)
        with self.subTest():
            self.assertEqual(self.gnm.seq[:8], expected_seq_start)
        with self.subTest():
            self.assertEqual(self.gnm.seq[-8:], expected_seq_end)
        with self.subTest():
            self.assertEqual(self.gnm.type, expected_type)


    def test_parse_genome_data_7(self):
        """Verify output when there is no fasta_file key."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file_x": url}
        self.gnm = phagesdb.parse_genome_data(data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_cluster = "A"
        expected_subcluster = "A2"
        expected_host = "Mycobacterium"
        expected_accession = "ABC123"
        expected_filename = ""
        expected_seq = ""
        expected_type = "phagesdb"

        with self.subTest():
            self.assertEqual(self.gnm.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.gnm.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.gnm.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.gnm.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.gnm.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.gnm.filename, expected_filename)
        with self.subTest():
            self.assertEqual(self.gnm.seq, expected_seq)
        with self.subTest():
            self.assertEqual(self.gnm.type, expected_type)


    def test_parse_genome_data_8(self):
        """Verify output when there is incorrect fasta_file URL."""

        url = "https://phagesdb.org/media/fastas/L5_x.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict)

        expected_phage_name = "Trixie"
        expected_phage_id = "Trixie"
        expected_cluster = "A"
        expected_subcluster = "A2"
        expected_host = "Mycobacterium"
        expected_accession = "ABC123"
        expected_filename = url
        expected_seq = ""
        expected_type = "phagesdb"

        with self.subTest():
            self.assertEqual(self.gnm.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.gnm.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.gnm.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.gnm.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.gnm.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.gnm.filename, expected_filename)
        with self.subTest():
            self.assertEqual(self.gnm.seq, expected_seq)
        with self.subTest():
            self.assertEqual(self.gnm.type, expected_type)


    def test_parse_genome_data_9(self):
        """Verify output when there are multiple errors."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name_x":"Trixie",
                    "pcluster_x": {"cluster": "A"},
                    "psubcluster_x": {"subcluster": "A2"},
                    "isolation_host_x": {"genus": "Mycobacterium"},
                    "genbank_accession_x": "ABC123",
                    "fasta_file_x": url}
        self.gnm = phagesdb.parse_genome_data(data_dict)

        expected_phage_name = ""
        expected_phage_id = ""
        expected_cluster = ""
        expected_subcluster = ""
        expected_host = ""
        expected_accession = ""
        expected_filename = ""
        expected_seq = ""
        expected_type = "phagesdb"

        with self.subTest():
            self.assertEqual(self.gnm.name, expected_phage_name)
        with self.subTest():
            self.assertEqual(self.gnm.id, expected_phage_id)
        with self.subTest():
            self.assertEqual(self.gnm.cluster, expected_cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, expected_subcluster)
        with self.subTest():
            self.assertEqual(self.gnm.host_genus, expected_host)
        with self.subTest():
            self.assertEqual(self.gnm.accession, expected_accession)
        with self.subTest():
            self.assertEqual(self.gnm.filename, expected_filename)
        with self.subTest():
            self.assertEqual(self.gnm.seq, expected_seq)
        with self.subTest():
            self.assertEqual(self.gnm.type, expected_type)




    def test_retrieve_genome_data_1(self):
        """Verify data is retrieved from PhagesDB with no error produced."""

        url = self.API_PREFIX + "L5" + self.API_SUFFIX
        data_dict = phagesdb.retrieve_genome_data(url)
        expected_phage_name = "L5"
        self.assertEqual(data_dict["phage_name"], expected_phage_name)

    def test_retrieve_genome_data_2(self):
        """Verify data is not retrieved from PhagesDB and an error produced."""

        url = self.API_PREFIX + "L5_x" + self.API_SUFFIX
        data_dict = phagesdb.retrieve_genome_data(url)
        self.assertEqual(len(data_dict.keys()), 0)




    def test_retrieve_data_list_1(self):
        """Confirm that data is successfully retrieved."""
        url = constants.API_CLUSTERS
        data_list = phagesdb.retrieve_data_list(url)
        self.assertTrue(len(data_list) > 0)

    def test_retrieve_data_list_2(self):
        """Confirm that data is not successfully retrieved."""
        url = 'invalid url'
        data_list = phagesdb.retrieve_data_list(url)
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


        self.genome1 = genome.Genome()
        self.genome1.id = "L5"
        self.genome1.type = "add"
        self.genome1.host_genus = "Gordonia"
        self.genome1.cluster = "B"
        self.genome1._value_flag = True

        self.bundle1 = bundle.Bundle()




    def test_copy_data_from_1(self):
        """Check that an "add" genome with no fields set to 'retrieve' is
        not impacted."""

        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        phagesdb.copy_data_from(self.bundle1, "add")
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
        """Check that an "add" genome with host_genus field set to 'retrieve' is
        populated correctly."""

        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.genome1.host_genus = "retrieve"
        phagesdb.copy_data_from(self.bundle1, "add")
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
        """Check that an "invalid" genome with host_genus field set to 'retrieve' is
        not populated correctly."""

        self.genome1.type = "invalid"
        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.genome1.host_genus = "retrieve"
        phagesdb.copy_data_from(self.bundle1, "add")
        with self.subTest():
            self.assertEqual(
                len(self.bundle1.genome_pair_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(self.genome1.host_genus, "retrieve")
        with self.subTest():
            self.assertEqual(len(self.genome1.evaluations), 0)

    def test_copy_data_from_4(self):
        """Check that an "add" genome with host_genus field set to 'retrieve' is
        not populated correctly when "invalid" type is requrested."""

        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.genome1.host_genus = "retrieve"
        phagesdb.copy_data_from(self.bundle1, "invalid")
        with self.subTest():
            self.assertEqual(
                len(self.bundle1.genome_pair_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(self.genome1.host_genus, "retrieve")
        with self.subTest():
            self.assertEqual(len(self.genome1.evaluations), 0)

    def test_copy_data_from_5(self):
        """Check that an "add" genome with host_genus field set to 'retrieve' is
        not populated correctly when id is not valid."""

        self.genome1.id = "invalid"
        self.bundle1.genome_dict[self.genome1.type] = self.genome1
        self.genome1.host_genus = "retrieve"
        self.genome1._value_flag = False
        phagesdb.copy_data_from(self.bundle1, "add")
        with self.subTest():
            self.assertTrue(self.genome1._value_flag)
        with self.subTest():
            self.assertEqual(
                len(self.bundle1.genome_pair_dict.keys()), 0)
        with self.subTest():
            self.assertEqual(self.genome1.host_genus, "retrieve")
        with self.subTest():
            self.assertEqual(self.genome1.evaluations[0].status, "error")





if __name__ == '__main__':
    unittest.main()
