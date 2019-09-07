"""Unit tests for misc. functions that interact with PhagesDB."""

from pdm_utils.functions import phagesdb
from pdm_utils.constants import constants
import unittest



class TestPhagesDBFunctions(unittest.TestCase):


    def setUp(self):

        self.API_PREFIX = constants.API_PREFIX
        self.API_SUFFIX = constants.API_SUFFIX


    def test_parse_phage_name_1(self):
        """Verify name is retrieved and no error is produced."""

        data_dict = {"phage_name":"Trixie"}

        phage = phagesdb.parse_phage_name(data_dict)
        expected_phage = "Trixie"
        self.assertEqual(phage, expected_phage)

    def test_parse_phage_name_2(self):
        """Verify name is not retrieved and error is produced."""

        data_dict = {"phage_id":"Trixie"}
        phage = phagesdb.parse_phage_name(data_dict)
        expected_phage = ""
        self.assertEqual(phage, expected_phage)




    def test_parse_cluster_1(self):
        """Verify missing cluster is retrieved and no error is produced."""

        data_dict = {"pcluster": None}
        cluster = phagesdb.parse_cluster(data_dict)
        expected_cluster = "UNK"
        self.assertEqual(cluster, expected_cluster)

    def test_parse_cluster_2(self):
        """Verify no cluster is retrieved and error is produced."""

        data_dict = {"pcluster_x": None}
        cluster = phagesdb.parse_cluster(data_dict)
        expected_cluster = ""
        self.assertEqual(cluster, expected_cluster)

    def test_parse_cluster_3(self):
        """Verify standard cluster is retrieved and no error is produced."""

        data_dict = {"pcluster": {"cluster": "A"}}
        cluster = phagesdb.parse_cluster(data_dict)
        expected_cluster = "A"
        self.assertEqual(cluster, expected_cluster)




    def test_parse_subcluster_1(self):
        """Verify missing subcluster is retrieved and no error is produced."""

        data_dict = {"psubcluster": None}
        subcluster = phagesdb.parse_subcluster(data_dict)
        expected_subcluster = "none"
        self.assertEqual(subcluster, expected_subcluster)

    def test_parse_subcluster_2(self):
        """Verify no subcluster is retrieved and an error is produced."""

        data_dict = {"psubcluster_x": None}
        subcluster = phagesdb.parse_subcluster(data_dict)
        expected_subcluster = ""
        self.assertEqual(subcluster, expected_subcluster)

    def test_parse_subcluster_3(self):
        """Verify standard subcluster is retrieved and no error is produced."""

        data_dict = {"psubcluster": {"subcluster": "A2"}}
        subcluster = phagesdb.parse_subcluster(data_dict)
        expected_subcluster = "A2"
        self.assertEqual(subcluster, expected_subcluster)




    def test_parse_host_genus_1(self):
        """Verify host genus is retrieved and no error is produced."""

        data_dict = {"isolation_host": {"genus": "Mycobacterium"}}
        host = phagesdb.parse_host_genus(data_dict)
        expected_host = "Mycobacterium"
        self.assertEqual(host, expected_host)

    def test_parse_host_genus_2(self):
        """Verify host genus is not retrieved and an error is produced."""

        data_dict = {"isolation_host_x": {"genus": "Mycobacterium"}}
        host = phagesdb.parse_host_genus(data_dict)
        expected_host = ""
        self.assertEqual(host, expected_host)




    def test_parse_accession_1(self):
        """Verify accession is retrieved and no error is produced."""

        data_dict = {"genbank_accession": "ABC123"}
        accession = phagesdb.parse_accession(data_dict)
        expected_accession = "ABC123"
        self.assertEqual(accession, expected_accession)

    def test_parse_accession_2(self):
        """Verify accession is not retrieved and an error is produced."""

        data_dict = {"genbank_accession_x": "ABC123"}
        accession = phagesdb.parse_accession(data_dict)
        expected_accession = ""
        self.assertEqual(accession, expected_accession)




    def test_parse_fasta_filename_1(self):
        """Verify fasta filename is retrieved and no error is produced."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"fasta_file": url}
        filename = phagesdb.parse_fasta_filename(data_dict)
        expected_filename = url
        self.assertEqual(filename, expected_filename)

    def test_parse_fasta_filename_2(self):
        """Verify fasta filename is not retrieved and an error is produced."""

        url = "https://phagesdb.org/media/fastas/L5.fasta"
        data_dict = {"fasta_file_x": url}
        filename = phagesdb.parse_fasta_filename(data_dict)
        expected_filename = ""
        self.assertEqual(filename, expected_filename)




    def test_parse_fasta_data_1(self):
        """Verify it fasta file data is parsed correctly."""
        fasta_data = ">Trixie  \nAAAAAAAAAA   \nTTTTTTT \nCCC\nGGGGGGGGGGG\n\n"
        expected_header = "Trixie"
        expected_sequence = "AAAAAAAAAATTTTTTTCCCGGGGGGGGGGG"
        result_tuple = phagesdb.parse_fasta_data(fasta_data)
        with self.subTest():
            self.assertEqual(result_tuple[0], expected_header)
        with self.subTest():
            self.assertEqual(result_tuple[1], expected_sequence)

    def test_parse_fasta_data_2(self):
        """Verify it incorrect fasta file format (no ">") produces an error."""
        fasta_data = "Trixie  \nAAAAAAAAAA   \nTTTTTTT \nCCC\nGGGGGGGGGGG\n\n"
        expected_header = "Trixie"
        expected_sequence = "AAAAAAAAAATTTTTTTCCCGGGGGGGGGGG"
        result_tuple = phagesdb.parse_fasta_data(fasta_data)
        with self.subTest():
            self.assertEqual(result_tuple[0], expected_header)
        with self.subTest():
            self.assertEqual(result_tuple[1], expected_sequence)

    def test_parse_fasta_data_3(self):
        """Verify it incorrect fasta file format (no new lines) produces
        an error."""
        fasta_data = "Trixie  AAAAAAAAAA   TTTTTTT CCCGGGGGGGGGGG"
        expected_header = ""
        expected_sequence = ""
        result_tuple = phagesdb.parse_fasta_data(fasta_data)
        with self.subTest():
            self.assertEqual(result_tuple[0], expected_header)
        with self.subTest():
            self.assertEqual(result_tuple[1], expected_sequence)




    def test_construct_phage_url_1(self):
        """Verify URL is constructed correctly."""

        phage_name = "L5"
        expected_url = self.API_PREFIX + "L5" + self.API_SUFFIX
        url = phagesdb.construct_phage_url(phage_name)
        self.assertEqual(url, expected_url)




if __name__ == '__main__':
    unittest.main()
