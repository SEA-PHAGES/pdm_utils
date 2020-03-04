"""Integration tests for misc. functions that interact with PhagesDB."""

from pdm_utils.classes import bundle
from pdm_utils.functions import phagesdb
from pdm_utils.classes import genome
from pdm_utils.constants import constants
import unittest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path



class TestPhagesDBFunctions(unittest.TestCase):


    def setUp(self):

        self.API_PREFIX = constants.API_PREFIX
        self.API_SUFFIX = constants.API_SUFFIX

        self.gnm = genome.Genome()




    def test_retrieve_url_data_1(self):
        """Verify fasta data is retrieved and no error is produced."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        fasta_data = phagesdb.retrieve_url_data(url)
        expected_fasta_data_header = ">Mycobacterium phage L5"
        self.assertEqual(fasta_data[:23], expected_fasta_data_header)

    def test_retrieve_url_data_2(self):
        """Verify fasta data is not retrieved and an error is produced."""

        url = "https://phagesdb.org/media/fastas/L5_x.fasta"
        fasta_data = phagesdb.retrieve_url_data(url)
        expected_fasta_data_header = ""
        self.assertEqual(fasta_data, expected_fasta_data_header)




    def test_parse_genome_data_1(self):
        """Verify genome object is parsed from PhagesDB."""
        url = "https://phagesdb.org/media/fastas/L5.fasta"
        filename = Path(url).stem
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict,
                                              gnm_type="phagesdb",
                                              seq=True)
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
            self.assertEqual(self.gnm.filename, filename)
        with self.subTest():
            self.assertEqual(self.gnm.seq[:8], "GGTCGGTT")
        with self.subTest():
            self.assertEqual(self.gnm.seq[-8:], "GTCGGTTA")
        with self.subTest():
            self.assertIsInstance(self.gnm.seq, Seq)
        with self.subTest():
            self.assertEqual(self.gnm.description, "Mycobacterium phage L5")
        with self.subTest():
            self.assertEqual(self.gnm._description_name, "L5")
        with self.subTest():
            self.assertEqual(self.gnm._description_host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(self.gnm.type, "phagesdb")
        with self.subTest():
            self.assertIsInstance(self.gnm.misc, dict)


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
        with self.subTest():
            self.assertEqual(self.gnm.name, "")
        with self.subTest():
            self.assertEqual(self.gnm.id, "")



    def test_parse_genome_data_3(self):
        """Verify output when there is no pcluster key."""
        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster_x": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict, gnm_type="phagesdb")
        with self.subTest():
            self.assertEqual(self.gnm.name, "Trixie")
        with self.subTest():
            self.assertEqual(self.gnm.cluster, "")


    def test_parse_genome_data_4(self):
        """Verify output when there is no psubcluster key."""
        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster_x": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict, gnm_type="phagesdb")
        self.assertEqual(self.gnm.subcluster, "")


    def test_parse_genome_data_5(self):
        """Verify output when there is no isolation_host key."""
        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host_x": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict, gnm_type="phagesdb")
        self.assertEqual(self.gnm.host_genus, "")


    def test_parse_genome_data_6(self):
        """Verify output when there is no genbank_accession key."""
        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession_x": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict, gnm_type="phagesdb")
        self.assertEqual(self.gnm.accession, "")


    def test_parse_genome_data_7(self):
        """Verify output when there is no fasta_file key."""
        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file_x": url}
        self.gnm = phagesdb.parse_genome_data(data_dict,
                                              gnm_type="phagesdb",
                                              seq=True)
        with self.subTest():
            self.assertEqual(self.gnm.filename, "")
        with self.subTest():
            self.assertEqual(self.gnm.seq, "")
        with self.subTest():
            self.assertEqual(self.gnm.description, "")
        with self.subTest():
            self.assertEqual(self.gnm._description_name, "")
        with self.subTest():
            self.assertEqual(self.gnm._description_host_genus, "")


    def test_parse_genome_data_8(self):
        """Verify output when there is incorrect fasta_file URL."""
        url = "https://phagesdb.org/media/fastas/L5_x.fasta"
        filename = Path(url).stem
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict,
                                              gnm_type="phagesdb",
                                              seq=True)
        with self.subTest():
            self.assertEqual(self.gnm.filename, filename)
        with self.subTest():
            self.assertEqual(self.gnm.seq, "")
        with self.subTest():
            self.assertEqual(self.gnm.description, "")
        with self.subTest():
            self.assertEqual(self.gnm._description_name, "")
        with self.subTest():
            self.assertEqual(self.gnm._description_host_genus, "")


    def test_parse_genome_data_9(self):
        """Verify output when there are multiple errors."""
        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"phage_name":"Trixie",
                    "pcluster_x": {"cluster": "A"},
                    "psubcluster_x": {"subcluster": "A2"},
                    "isolation_host_x": {"genus": "Mycobacterium"},
                    "genbank_accession_x": "ABC123",
                    "fasta_file_x": url}
        self.gnm = phagesdb.parse_genome_data(data_dict, seq=True)
        with self.subTest():
            self.assertEqual(self.gnm.name, "Trixie")
        with self.subTest():
            self.assertEqual(self.gnm.id, "Trixie")
        with self.subTest():
            self.assertEqual(self.gnm.cluster, "")
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, "")
        with self.subTest():
            self.assertEqual(self.gnm.host_genus, "")
        with self.subTest():
            self.assertEqual(self.gnm.accession, "")
        with self.subTest():
            self.assertEqual(self.gnm.filename, "")
        with self.subTest():
            self.assertEqual(self.gnm.seq, "")
        with self.subTest():
            self.assertEqual(self.gnm.description, "")
        with self.subTest():
            self.assertEqual(self.gnm._description_name, "")
        with self.subTest():
            self.assertEqual(self.gnm._description_host_genus, "")
        with self.subTest():
            self.assertEqual(self.gnm.type, "")


    def test_parse_genome_data_10(self):
        """Verify output when seq is False."""
        url = "https://phagesdb.org/media/fastas/L5.fasta"
        filename = Path(url).stem
        data_dict = {"phage_name":"Trixie",
                    "pcluster": {"cluster": "A"},
                    "psubcluster": {"subcluster": "A2"},
                    "isolation_host": {"genus": "Mycobacterium"},
                    "genbank_accession": "ABC123",
                    "fasta_file": url}
        self.gnm = phagesdb.parse_genome_data(data_dict,
                                              gnm_type="phagesdb",
                                              seq=False)
        with self.subTest():
            self.assertEqual(self.gnm.filename, filename)
        with self.subTest():
            self.assertEqual(self.gnm.seq, "")
        with self.subTest():
            self.assertEqual(self.gnm.description, "")
        with self.subTest():
            self.assertEqual(self.gnm._description_name, "")
        with self.subTest():
            self.assertEqual(self.gnm._description_host_genus, "")




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

        self.bundle1 = bundle.Bundle()




    def test_get_genome_1(self):
        """Check that a genome is successfully retrieved."""
        gnm = phagesdb.get_genome("Trixie", gnm_type="phagesdb")
        with self.subTest():
            self.assertEqual(gnm.id, "Trixie")
        with self.subTest():
            self.assertEqual(gnm.host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(gnm.cluster, "A")
        with self.subTest():
            self.assertEqual(gnm.type, "phagesdb")

    def test_get_genome_2(self):
        """Check that a genome is not retrieved."""
        gnm = phagesdb.get_genome("invalid", gnm_type="phagesdb")
        self.assertIsNone(gnm)


if __name__ == '__main__':
    unittest.main()
