""" Unit tests for the GenomePair Class."""


from classes import GenomePair
from classes import Genome
from classes import Ticket
from classes import Eval
import unittest


class TestGenomePairClass(unittest.TestCase):


    def setUp(self):
        self.genome1 = Genome.Genome()
        self.genome2 = Genome.Genome()
        self.ticket = Ticket.GenomeTicket()
        self.genome_pair = GenomePair.GenomePair()
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2











    def test_copy_data_1(self):
        """Check that phage_name is copied from genome1 to genome2,
        type is the identifying attribute, using 'copy' keyword."""
        self.genome1.phage_id = "L5"
        self.genome1.phage_name = "Trixie"
        self.genome1.type = "phamerator"

        self.genome2.phage_id = "D29"
        self.genome2.phage_name = "copy"
        self.genome2.type = "import"

        self.genome_pair.copy_data("type", "phamerator", "import", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.phage_name, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.phage_id, "D29")

    def test_copy_data_2(self):
        """Check that phage_name is copied from genome2 to genome1,
        type is the identifying attribute, using 'copy' keyword."""

        self.genome1.phage_id = "L5"
        self.genome1.phage_name = "copy"
        self.genome1.type = "phamerator"

        self.genome2.phage_id = "D29"
        self.genome2.phage_name = "Trixie"
        self.genome2.type = "import"

        self.genome_pair.copy_data("type", "import", "phamerator", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.phage_name, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.phage_id, "L5")

    def test_copy_data_3(self):
        """Check that no phage_name is copied from genome1 to genome2,
        type is the identifying attribute, using 'copy' keyword,
        since genome2 does not contain keyword."""

        self.genome1.phage_id = "L5"
        self.genome1.phage_name = "copy"
        self.genome1.type = "phamerator"

        self.genome2.phage_id = "D29"
        self.genome2.phage_name = "Trixie"
        self.genome2.type = "import"

        self.genome_pair.copy_data("type", "phamerator", "import", "copy")
        self.assertEqual(self.genome_pair.genome2.phage_name, "Trixie")

    def test_copy_data_4(self):
        """Check that phage_name is copied from genome1 to genome2,
        phage_id is the identifying attribute, using 'copy' keyword."""
        self.genome1.phage_id = "L5"
        self.genome1.phage_name = "Trixie"
        self.genome1.type = "phamerator"

        self.genome2.phage_id = "D29"
        self.genome2.phage_name = "copy"
        self.genome2.type = "import"

        self.genome_pair.copy_data("phage_id", "L5", "D29", "copy")
        self.assertEqual(self.genome_pair.genome2.phage_name, "Trixie")

    def test_copy_data_5(self):
        """Check that phage_name is copied from genome2 to genome1,
        phage_id is the identifying attribute, using 'copy' keyword."""

        self.genome1.phage_id = "L5"
        self.genome1.phage_name = "copy"
        self.genome1.type = "phamerator"

        self.genome2.phage_id = "D29"
        self.genome2.phage_name = "Trixie"
        self.genome2.type = "import"

        self.genome_pair.copy_data("phage_id", "D29", "L5", "copy")
        self.assertEqual(self.genome_pair.genome1.phage_name, "Trixie")

    def test_copy_data_6(self):
        """Check that phage_name is copied from genome2 to genome1,
        phage_id is the identifying attribute, using 'copy' keyword."""

        self.genome1.phage_id = "L5"
        self.genome1.phage_name = "copy"
        self.genome1.type = "phamerator"

        self.genome2.phage_id = "D29"
        self.genome2.phage_name = "Trixie"
        self.genome2.type = "import"

        self.genome_pair.copy_data("phage_id", "D29", "L5", "copy")
        self.assertEqual(self.genome_pair.genome1.phage_name, "Trixie")

    def test_copy_data_7(self):
        """Check that no data is copied if direction of copy is
        not determined since the 'first' and 'second' parameter values
        are not found in the 'attr' attribute."""
        self.genome1.phage_id = "L5"
        self.genome1.phage_name = "Trixie"
        self.genome1.type = "phamerator"

        self.genome2.phage_id = "D29"
        self.genome2.phage_name = "copy"
        self.genome2.type = "phamerator"

        self.genome_pair.copy_data("type", "L5", "D29", "copy")
        self.assertEqual(self.genome_pair.genome2.phage_name, "copy")


    def test_copy_data_8(self):
        """Check that no data is copied if direction of copy is
        not determined since the 'first' and 'second' parameter values
        are identical."""
        self.genome1.phage_id = "L5"
        self.genome1.phage_name = "Trixie"
        self.genome1.host = "copy"
        self.genome1.type = "phamerator"

        self.genome2.phage_id = "L5"
        self.genome2.phage_name = "copy"
        self.genome2.host = "Mycobacterium"
        self.genome2.type = "phamerator"

        self.genome_pair.copy_data("phage_id", "L5", "L5", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.host, "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.phage_name, "copy")

    def test_copy_data_9(self):
        """Check that no data is copied if direction of copy is
        not determined since the 'first' parameter value is not found
        in either genome, using 'phage_id' attribute."""
        self.genome1.phage_id = "L5"
        self.genome1.phage_name = "Trixie"
        self.genome1.host = "copy"
        self.genome1.type = "phamerator"

        self.genome2.phage_id = "D29"
        self.genome2.phage_name = "copy"
        self.genome2.host = "Mycobacterium"
        self.genome2.type = "phamerator"

        self.genome_pair.copy_data("phage_id", "Trixie", "L5", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.host, "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.phage_name, "copy")

    def test_copy_data_10(self):
        """Check that no data is copied if direction of copy is
        not determined since the 'second' parameter value is not found
        in either genome, using 'phage_id' attribute."""
        self.genome1.phage_id = "L5"
        self.genome1.phage_name = "Trixie"
        self.genome1.host = "copy"
        self.genome1.type = "phamerator"

        self.genome2.phage_id = "D29"
        self.genome2.phage_name = "copy"
        self.genome2.host = "Mycobacterium"
        self.genome2.type = "phamerator"

        self.genome_pair.copy_data("phage_id", "L5", "Trixie", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.host, "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.phage_name, "copy")

    def test_copy_data_11(self):
        """Check that no data is copied if direction of copy is
        not determined since the 'first' parameter value is not found
        in either genome, using 'type' attribute."""
        self.genome1.phage_id = "L5"
        self.genome1.phage_name = "Trixie"
        self.genome1.host = "copy"
        self.genome1.type = "phamerator"

        self.genome2.phage_id = "D29"
        self.genome2.phage_name = "copy"
        self.genome2.host = "Mycobacterium"
        self.genome2.type = "import"

        self.genome_pair.copy_data("type", "L5", "import", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.host, "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.phage_name, "copy")

    def test_copy_data_12(self):
        """Check that no data is copied if direction of copy is
        not determined since the 'second' parameter value is not found
        in either genome, using 'type' attribute."""
        self.genome1.phage_id = "L5"
        self.genome1.phage_name = "Trixie"
        self.genome1.host = "copy"
        self.genome1.type = "phamerator"

        self.genome2.phage_id = "D29"
        self.genome2.phage_name = "copy"
        self.genome2.host = "Mycobacterium"
        self.genome2.type = "import"

        self.genome_pair.copy_data("type", "phamerator", "Trixie", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.host, "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.phage_name, "copy")

    def test_copy_data_13(self):
        """Check that all data (except for 'type' attribute)
        is copied if no keyword is provided."""

        self.genome1.type = "phamerator"
        self.genome1.phage_id = 1
        self.genome1.phage_name = 2
        self.genome1.host = 3
        self.genome1.sequence = 4
        self.genome1.accession = 5
        self.genome1.author = 6
        self.genome1.status = 7
        self.genome1.cluster = 8
        self.genome1.subcluster = 9
        self.genome1.cluster_subcluster = 10
        self.genome1.date_last_modified = 12
        self.genome1.annotation_author = 13
        self.genome1.annotation_qc = 14
        self.genome1.retrieve_record = 15
        self.genome1.record_name = 16
        self.genome1.record_id = 17
        self.genome1.record_accession = 18
        self.genome1.record_description = 19
        self.genome1.record_source = 20
        self.genome1.record_organism = 21
        self.genome1.record_authors = 22
        self.genome1.record_date = 23
        self.genome1.record_filename = 24
        self.genome1.translation_table = 25
        self.genome1.record = 26
        self.genome1.search_id = 27
        self.genome1.search_name = 28
        self.genome1._length = 29
        self.genome1._gc = 30
        self.genome1.evaluations = 31
        self.genome1.cds_features = 32
        self.genome1._cds_features_tally = 33
        self.genome1._cds_start_end_ids = 34
        self.genome1._cds_end_strand_ids = 35
        self.genome1._cds_processed_primary_descriptions_tally = 36
        self.genome1.trna_features = 37
        self.genome1._trna_features_tally = 38
        self.genome1.source_features = 39
        self.genome1._source_features_tally = 40
        self.genome1.search_record_filename = 41
        self.genome1._record_description_phage_name = 42
        self.genome1._record_source_phage_name = 43
        self.genome1._record_organism_phage_name = 44
        self.genome1._record_description_host_name = 45
        self.genome1._record_source_host_name = 46
        self.genome1._record_organism_host_name = 47
        self.genome1._cds_processed_product_descriptions_tally = 48
        self.genome1._cds_processed_function_descriptions_tally = 49
        self.genome1._cds_processed_note_descriptions_tally = 50
        self.genome1._cds_unique_start_end_ids = 51
        self.genome1._cds_duplicate_start_end_ids = 52
        self.genome1._cds_unique_end_strand_ids = 53
        self.genome1._cds_duplicate_end_strand_ids = 54

        self.genome2.type = "import"

        self.genome_pair.copy_data("type", "phamerator", "import")

        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.type, "import")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.phage_id, 1)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.phage_name, 2)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.host, 3)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.sequence, 4)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.accession, 5)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.author, 6)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.status, 7)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.cluster, 8)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.subcluster, 9)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.cluster_subcluster, 10)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.date_last_modified, 12)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.annotation_author, 13)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.annotation_qc, 14)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.retrieve_record, 15)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.record_name, 16)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.record_id, 17)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.record_accession, 18)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.record_description, 19)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.record_source, 20)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.record_organism, 21)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.record_authors, 22)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.record_date, 23)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.record_filename, 24)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.translation_table, 25)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.record, 26)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.search_id, 27)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.search_name, 28)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2._length, 29)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2._gc, 30)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.evaluations, 31)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.cds_features, 32)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2._cds_features_tally, 33)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2._cds_start_end_ids, 34)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2._cds_end_strand_ids, 35)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._cds_processed_primary_descriptions_tally,
                36)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.trna_features, 37)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2._trna_features_tally, 38)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.source_features, 39)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2._source_features_tally, 40)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.search_record_filename, 41)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._record_description_phage_name, 42)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._record_source_phage_name, 43)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._record_organism_phage_name, 44)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._record_description_host_name, 45)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._record_source_host_name, 46)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._record_organism_host_name, 47)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._cds_processed_product_descriptions_tally,
                48)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._cds_processed_function_descriptions_tally,
                49)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._cds_processed_note_descriptions_tally,
                50)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._cds_unique_start_end_ids, 51)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._cds_duplicate_start_end_ids, 52)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._cds_unique_end_strand_ids, 53)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._cds_duplicate_end_strand_ids, 54)

    def test_copy_data_14(self):
        """Check that all data (except for 'phage_id' attribute)
        is copied if no keyword is provided."""

        self.genome1.phage_id = "L5"
        self.genome1.type = 1
        self.genome1.phage_name = 2

        self.genome2.phage_id = "Trixie"

        self.genome_pair.copy_data("phage_id", "L5", "Trixie")

        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.type, 1)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.phage_name, 2)




    def test_compare_phage_id_1(self):
        """Check that an error is produced if the
        phage_id is not the same."""
        self.genome1.phage_id = "Trixie"
        self.genome2.phage_id = "L5"
        self.genome_pair.compare_phage_id()
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")

    def test_compare_phage_id_2(self):
        """Check that no error is produced if the
        phage_id is the same."""
        self.genome1.phage_id = "Trixie"
        self.genome2.phage_id = "Trixie"
        self.genome_pair.compare_phage_id()
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")




    def test_compare_genome_sequence_1(self):
        """Check that identical sequences produce no warning."""
        self.genome1.sequence = "ABCD"
        self.genome2.sequence = "ABCD"
        self.genome_pair.compare_genome_sequence()
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_genome_sequence_2(self):
        """Check that different sequences produce a warning."""
        self.genome1.sequence = "ABCD"
        self.genome2.sequence = "ABCDE"
        self.genome_pair.compare_genome_sequence()
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")




    def test_compare_genome_length_1(self):
        """Check that identical sequence lengths produce no warning."""
        self.genome1._length = 5
        self.genome2._length = 5
        self.genome_pair.compare_genome_length()
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_genome_length_2(self):
        """Check that different sequence lengths produce a warning."""
        self.genome1._length = 5
        self.genome2._length = 6
        self.genome_pair.compare_genome_length()
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")




    def test_compare_cluster_1(self):
        """Check that identical clusters produce no warning."""
        self.genome1.cluster = "A"
        self.genome2.cluster = "A"
        self.genome_pair.compare_cluster()
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_cluster_2(self):
        """Check that different clusters produce a warning."""
        self.genome1.cluster = "A"
        self.genome2.cluster = "B"
        self.genome_pair.compare_cluster()
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")




    def test_compare_subcluster_1(self):
        """Check that identical subclusters produce no warning."""
        self.genome1.subcluster = "A1"
        self.genome2.subcluster = "A1"
        self.genome_pair.compare_subcluster()
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_subcluster_2(self):
        """Check that different subclusters produce a warning."""
        self.genome1.subcluster = "A1"
        self.genome2.subcluster = "A2"
        self.genome_pair.compare_subcluster()
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")




    def test_compare_accession_1(self):
        """Check that identical accessions produce no warning."""
        self.genome1.accession = "ABC123"
        self.genome2.accession = "ABC123"
        self.genome_pair.compare_accession()
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_accession_2(self):
        """Check that different accessions produce a warning."""
        self.genome1.accession = "ABC1234"
        self.genome2.accession = "ABC123"
        self.genome_pair.compare_accession()
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")




    def test_compare_host_1(self):
        """Check that identical hosts produce no warning."""
        self.genome1.host = "Mycobacterium"
        self.genome2.host = "Mycobacterium"
        self.genome_pair.compare_host()
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_host_2(self):
        """Check that different hosts produce a warning."""
        self.genome1.host = "Mycobacterium"
        self.genome2.host = "Mycobacteriums"
        self.genome_pair.compare_host()
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")




if __name__ == '__main__':
    unittest.main()
