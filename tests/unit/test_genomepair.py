""" Unit tests for the GenomePair Class."""

from datetime import datetime
from pdm_utils.classes import genomepair
from pdm_utils.classes import genome
from pdm_utils.classes import ticket
from pdm_utils.classes import eval
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class TestGenomePairClass(unittest.TestCase):


    def setUp(self):
        self.genome1 = genome.Genome()
        self.genome2 = genome.Genome()
        self.tkt = ticket.GenomeTicket()
        self.genome_pair = genomepair.GenomePair()
        self.genome_pair.genome1 = self.genome1
        self.genome_pair.genome2 = self.genome2

        self.date_jan1 = datetime.strptime('1/1/2000', '%m/%d/%Y')
        self.date_feb1 = datetime.strptime('2/1/2000', '%m/%d/%Y')
        self.date_feb1_b = datetime.strptime('2/1/2000', '%m/%d/%Y')












    def test_copy_data_1(self):
        """Check that name is copied from genome1 to genome2,
        type is the identifying attribute, using 'copy' keyword."""
        self.genome1.id = "L5"
        self.genome1.name = "Trixie"
        self.genome1.type = "mysql"

        self.genome2.id = "D29"
        self.genome2.name = "copy"
        self.genome2.type = "import"

        self.genome_pair.copy_data("type", "mysql", "import", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.name, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.id, "D29")

    def test_copy_data_2(self):
        """Check that name is copied from genome2 to genome1,
        type is the identifying attribute, using 'copy' keyword."""

        self.genome1.id = "L5"
        self.genome1.name = "copy"
        self.genome1.type = "mysql"

        self.genome2.id = "D29"
        self.genome2.name = "Trixie"
        self.genome2.type = "import"

        self.genome_pair.copy_data("type", "import", "mysql", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.name, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.id, "L5")

    def test_copy_data_3(self):
        """Check that no name is copied from genome1 to genome2,
        type is the identifying attribute, using 'copy' keyword,
        since genome2 does not contain keyword."""

        self.genome1.id = "L5"
        self.genome1.name = "copy"
        self.genome1.type = "mysql"

        self.genome2.id = "D29"
        self.genome2.name = "Trixie"
        self.genome2.type = "import"

        self.genome_pair.copy_data("type", "mysql", "import", "copy")
        self.assertEqual(self.genome_pair.genome2.name, "Trixie")

    def test_copy_data_4(self):
        """Check that name is copied from genome1 to genome2,
        id is the identifying attribute, using 'copy' keyword."""
        self.genome1.id = "L5"
        self.genome1.name = "Trixie"
        self.genome1.type = "mysql"

        self.genome2.id = "D29"
        self.genome2.name = "copy"
        self.genome2.type = "import"

        self.genome_pair.copy_data("id", "L5", "D29", "copy")
        self.assertEqual(self.genome_pair.genome2.name, "Trixie")

    def test_copy_data_5(self):
        """Check that name is copied from genome2 to genome1,
        id is the identifying attribute, using 'copy' keyword."""

        self.genome1.id = "L5"
        self.genome1.name = "copy"
        self.genome1.type = "mysql"

        self.genome2.id = "D29"
        self.genome2.name = "Trixie"
        self.genome2.type = "import"

        self.genome_pair.copy_data("id", "D29", "L5", "copy")
        self.assertEqual(self.genome_pair.genome1.name, "Trixie")

    def test_copy_data_6(self):
        """Check that name is copied from genome2 to genome1,
        id is the identifying attribute, using 'copy' keyword."""

        self.genome1.id = "L5"
        self.genome1.name = "copy"
        self.genome1.type = "mysql"

        self.genome2.id = "D29"
        self.genome2.name = "Trixie"
        self.genome2.type = "import"

        self.genome_pair.copy_data("id", "D29", "L5", "copy")
        self.assertEqual(self.genome_pair.genome1.name, "Trixie")

    def test_copy_data_7(self):
        """Check that no data is copied if direction of copy is
        not determined since the 'first' and 'second' parameter values
        are not found in the 'attr' attribute."""
        self.genome1.id = "L5"
        self.genome1.name = "Trixie"
        self.genome1.type = "mysql"

        self.genome2.id = "D29"
        self.genome2.name = "copy"
        self.genome2.type = "mysql"

        self.genome_pair.copy_data("type", "L5", "D29", "copy")
        self.assertEqual(self.genome_pair.genome2.name, "copy")


    def test_copy_data_8(self):
        """Check that no data is copied if direction of copy is
        not determined since the 'first' and 'second' parameter values
        are identical."""
        self.genome1.id = "L5"
        self.genome1.name = "Trixie"
        self.genome1.host_genus = "copy"
        self.genome1.type = "mysql"

        self.genome2.id = "L5"
        self.genome2.name = "copy"
        self.genome2.host_genus = "Mycobacterium"
        self.genome2.type = "mysql"

        self.genome_pair.copy_data("id", "L5", "L5", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.host_genus, "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.name, "copy")

    def test_copy_data_9(self):
        """Check that no data is copied if direction of copy is
        not determined since the 'first' parameter value is not found
        in either genome, using 'id' attribute."""
        self.genome1.id = "L5"
        self.genome1.name = "Trixie"
        self.genome1.host_genus = "copy"
        self.genome1.type = "mysql"

        self.genome2.id = "D29"
        self.genome2.name = "copy"
        self.genome2.host_genus = "Mycobacterium"
        self.genome2.type = "mysql"

        self.genome_pair.copy_data("id", "Trixie", "L5", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.host_genus, "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.name, "copy")

    def test_copy_data_10(self):
        """Check that no data is copied if direction of copy is
        not determined since the 'second' parameter value is not found
        in either genome, using 'id' attribute."""
        self.genome1.id = "L5"
        self.genome1.name = "Trixie"
        self.genome1.host_genus = "copy"
        self.genome1.type = "mysql"

        self.genome2.id = "D29"
        self.genome2.name = "copy"
        self.genome2.host_genus = "Mycobacterium"
        self.genome2.type = "mysql"

        self.genome_pair.copy_data("id", "L5", "Trixie", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.host_genus, "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.name, "copy")

    def test_copy_data_11(self):
        """Check that no data is copied if direction of copy is
        not determined since the 'first' parameter value is not found
        in either genome, using 'type' attribute."""
        self.genome1.id = "L5"
        self.genome1.name = "Trixie"
        self.genome1.host_genus = "copy"
        self.genome1.type = "mysql"

        self.genome2.id = "D29"
        self.genome2.name = "copy"
        self.genome2.host_genus = "Mycobacterium"
        self.genome2.type = "import"

        self.genome_pair.copy_data("type", "L5", "import", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.host_genus, "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.name, "copy")

    def test_copy_data_12(self):
        """Check that no data is copied if direction of copy is
        not determined since the 'second' parameter value is not found
        in either genome, using 'type' attribute."""
        self.genome1.id = "L5"
        self.genome1.name = "Trixie"
        self.genome1.host_genus = "copy"
        self.genome1.type = "mysql"

        self.genome2.id = "D29"
        self.genome2.name = "copy"
        self.genome2.host_genus = "Mycobacterium"
        self.genome2.type = "import"

        self.genome_pair.copy_data("type", "mysql", "Trixie", "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome1.host_genus, "copy")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.name, "copy")

    def test_copy_data_13(self):
        """Check that all data (except for 'type' attribute)
        is copied if no keyword is provided."""

        self.genome1.type = "mysql"
        self.genome1.id = 1
        self.genome1.name = 2
        self.genome1.host_genus = 3
        self.genome1.seq = 4
        self.genome1.accession = 5
        self.genome1.annotation_status = 7
        self.genome1.cluster = 8
        self.genome1.subcluster = 9
        self.genome1.cluster_subcluster = 10
        self.genome1.date = 12
        self.genome1.annotation_author = 13
        self.genome1.retrieve_record = 15
        self.genome1.description = 19
        self.genome1.source = 20
        self.genome1.organism = 21
        self.genome1.authors = 22
        self.genome1.filename = 24
        self.genome1.translation_table = 25
        self.genome1.length = 29
        self.genome1.gc = 30
        self.genome1.evaluations = 31
        self.genome1.cds_features = 32
        self.genome1._cds_features_tally = 33
        self.genome1._cds_start_end_ids = 34
        self.genome1._cds_end_strand_ids = 35
        self.genome1._cds_processed_descriptions_tally = 36
        self.genome1.trna_features = 37
        self.genome1._trna_features_tally = 38
        self.genome1.source_features = 39
        self.genome1._source_features_tally = 40
        self.genome1._description_name = 42
        self.genome1._source_name = 43
        self.genome1._organism_name = 44
        self.genome1._description_host_genus = 45
        self.genome1._source_host_genus = 46
        self.genome1._organism_host_genus = 47
        self.genome1._cds_processed_products_tally = 48
        self.genome1._cds_processed_functions_tally = 49
        self.genome1._cds_processed_notes_tally = 50
        self.genome1._cds_unique_start_end_ids = 51
        self.genome1._cds_duplicate_start_end_ids = 52
        self.genome1._cds_unique_end_strand_ids = 53
        self.genome1._cds_duplicate_end_strand_ids = 54

        self.genome2.type = "import"

        self.genome_pair.copy_data("type", "mysql", "import")

        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.type, "import")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.id, 1)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.name, 2)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.host_genus, 3)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.seq, 4)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.accession, 5)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.annotation_status, 7)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.cluster, 8)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.subcluster, 9)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.cluster_subcluster, 10)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.date, 12)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.annotation_author, 13)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.retrieve_record, 15)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.description, 19)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.source, 20)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.organism, 21)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.authors, 22)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.filename, 24)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.translation_table, 25)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.length, 29)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.gc, 30)
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
                self.genome_pair.genome2._cds_processed_descriptions_tally,
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
            self.assertEqual(
                self.genome_pair.genome2._description_name, 42)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._source_name, 43)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._organism_name, 44)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._description_host_genus, 45)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._source_host_genus, 46)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._organism_host_genus, 47)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._cds_processed_products_tally,
                48)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._cds_processed_functions_tally,
                49)
        with self.subTest():
            self.assertEqual(
                self.genome_pair.genome2._cds_processed_notes_tally,
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
        """Check that all data (except for 'id' attribute)
        is copied if no keyword is provided."""

        self.genome1.id = "L5"
        self.genome1.type = 1
        self.genome1.name = 2

        self.genome2.id = "Trixie"

        self.genome_pair.copy_data("id", "L5", "Trixie")

        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.id, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.type, 1)
        with self.subTest():
            self.assertEqual(self.genome_pair.genome2.name, 2)




    def test_compare_attribute_1(self):
        """Verify no error is produced when both genomes have
        identical ids as expected."""
        self.genome1.id = "Trixie"
        self.genome2.id = "Trixie"
        self.genome_pair.compare_attribute("id", expect_same=True,
                                           eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")


    def test_compare_attribute_2(self):
        """Verify no error is produced when both genomes have
        different ids as expected."""
        self.genome1.id = "Trixie"
        self.genome2.id = "L5"
        self.genome_pair.compare_attribute("id", expect_same=False)
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.genome_pair.evaluations[0].id)


    def test_compare_attribute_3(self):
        """Verify an error is produced when both genomes have
        different ids unexpectedly."""
        self.genome1.id = "Trixie"
        self.genome2.id = "L5"
        self.genome_pair.compare_attribute("id", expect_same=True)
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")

    def test_compare_attribute_4(self):
        """Verify an error is produced when both genomes have
        identical ids unexpectedly."""
        self.genome1.id = "Trixie"
        self.genome2.id = "Trixie"
        self.genome_pair.compare_attribute("id", expect_same=False)
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")

    def test_compare_attribute_5(self):
        """Verify no error is produced when both genomes have
        different retrieve_records as expected."""
        self.genome1.retrieve_record = 1
        self.genome2.retrieve_record = 0
        self.genome_pair.compare_attribute("retrieve_record", expect_same=False)
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")


    def test_compare_attribute_6(self):
        """Verify no error is produced when both genomes have
        identical seqs as expected."""
        self.genome1.seq = Seq("AAAA", IUPAC.ambiguous_dna)
        self.genome2.seq = Seq("AAAA", IUPAC.ambiguous_dna)
        self.genome_pair.compare_attribute("seq", expect_same=True)
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_attribute_7(self):
        """Verify no error is produced when both genomes have
        different seqs as expected."""
        self.genome1.seq = Seq("AAAA", IUPAC.ambiguous_dna)
        self.genome2.seq = Seq("AAAT", IUPAC.ambiguous_dna)
        self.genome_pair.compare_attribute("seq", expect_same=False)
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_attribute_8(self):
        """Verify no test is performed when the attribute is invalid."""
        self.genome_pair.compare_attribute("invalid", expect_same=False)
        self.assertEqual(self.genome_pair.evaluations[0].status, "untested")







    # TODO the following tests are probably no longer needed since
    # compare_attribute() is functioning.
    # def test_compare_id_1(self):
    #     """Check that an error is produced if the
    #     id is not the same."""
    #     self.genome1.id = "Trixie"
    #     self.genome2.id = "L5"
    #     self.genome_pair.compare_id(eval_id="eval_id")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "error")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")
    #
    # def test_compare_id_2(self):
    #     """Check that no error is produced if the
    #     id is the same."""
    #     self.genome1.id = "Trixie"
    #     self.genome2.id = "Trixie"
    #     self.genome_pair.compare_id()
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
    #     with self.subTest():
    #         self.assertIsNone(self.genome_pair.evaluations[0].id)
    #
    #
    #
    #
    # def test_compare_genome_sequence_1(self):
    #     """Check that identical sequences produce no warning."""
    #     self.genome1.seq = Seq("AAAA", IUPAC.ambiguous_dna)
    #     self.genome2.seq = Seq("AAAA", IUPAC.ambiguous_dna)
    #     self.genome_pair.compare_genome_sequence(eval_id="eval_id")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")
    #
    # def test_compare_genome_sequence_2(self):
    #     """Check that different sequences produce a warning."""
    #     self.genome1.seq = Seq("AAAA", IUPAC.ambiguous_dna)
    #     self.genome2.seq = Seq("TTTT", IUPAC.ambiguous_dna)
    #     self.genome_pair.compare_genome_sequence()
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "error")
    #     with self.subTest():
    #         self.assertIsNone(self.genome_pair.evaluations[0].id)
    #
    #
    #
    #
    # def test_compare_genome_length_1(self):
    #     """Check that identical sequence lengths produce no warning."""
    #     self.genome1.length = 5
    #     self.genome2.length = 5
    #     self.genome_pair.compare_genome_length(eval_id="eval_id")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")
    #
    # def test_compare_genome_length_2(self):
    #     """Check that different sequence lengths produce a warning."""
    #     self.genome1.length = 5
    #     self.genome2.length = 6
    #     self.genome_pair.compare_genome_length()
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "error")
    #     with self.subTest():
    #         self.assertIsNone(self.genome_pair.evaluations[0].id)
    #
    #
    #
    #
    # def test_compare_cluster_1(self):
    #     """Check that identical clusters produce no warning."""
    #     self.genome1.cluster = "A"
    #     self.genome2.cluster = "A"
    #     self.genome_pair.compare_cluster(eval_id="eval_id")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")
    #
    # def test_compare_cluster_2(self):
    #     """Check that different clusters produce a warning."""
    #     self.genome1.cluster = "A"
    #     self.genome2.cluster = "B"
    #     self.genome_pair.compare_cluster()
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "error")
    #     with self.subTest():
    #         self.assertIsNone(self.genome_pair.evaluations[0].id)
    #
    #
    #
    #
    # def test_compare_subcluster_1(self):
    #     """Check that identical subclusters produce no warning."""
    #     self.genome1.subcluster = "A1"
    #     self.genome2.subcluster = "A1"
    #     self.genome_pair.compare_subcluster(eval_id="eval_id")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")
    #
    # def test_compare_subcluster_2(self):
    #     """Check that different subclusters produce a warning."""
    #     self.genome1.subcluster = "A1"
    #     self.genome2.subcluster = "A2"
    #     self.genome_pair.compare_subcluster()
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "error")
    #     with self.subTest():
    #         self.assertIsNone(self.genome_pair.evaluations[0].id)
    #
    #
    #
    #
    # def test_compare_accession_1(self):
    #     """Check that identical accessions produce no warning."""
    #     self.genome1.accession = "ABC123"
    #     self.genome2.accession = "ABC123"
    #     self.genome_pair.compare_accession(eval_id="eval_id")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")
    #
    # def test_compare_accession_2(self):
    #     """Check that different accessions produce a warning."""
    #     self.genome1.accession = "ABC1234"
    #     self.genome2.accession = "ABC123"
    #     self.genome_pair.compare_accession()
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "error")
    #     with self.subTest():
    #         self.assertIsNone(self.genome_pair.evaluations[0].id)
    #
    #
    #
    #
    # def test_compare_host_genus_1(self):
    #     """Check that identical hosts produce no warning."""
    #     self.genome1.host_genus = "Mycobacterium"
    #     self.genome2.host_genus = "Mycobacterium"
    #     self.genome_pair.compare_host_genus(eval_id="eval_id")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")
    #
    # def test_compare_host_genus_2(self):
    #     """Check that different hosts produce a warning."""
    #     self.genome1.host_genus = "Mycobacterium"
    #     self.genome2.host_genus = "Mycobacteriums"
    #     self.genome_pair.compare_host_genus()
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "error")
    #     with self.subTest():
    #         self.assertIsNone(self.genome_pair.evaluations[0].id)
    #
    #
    #
    #
    # def test_compare_annotation_author_1(self):
    #     """Verify no error is produced with
    #     identical annotation_author values."""
    #     self.genome1.annotation_author = 1
    #     self.genome2.annotation_author = 1
    #     self.genome_pair.compare_annotation_author(eval_id="eval_id")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")
    #
    # def test_compare_annotation_author_2(self):
    #     """Verify an error is produced with
    #     different annotation_author values."""
    #     self.genome1.annotation_author = 1
    #     self.genome2.annotation_author = 0
    #     self.genome_pair.compare_annotation_author()
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "error")
    #     with self.subTest():
    #         self.assertIsNone(self.genome_pair.evaluations[0].id)
    #



    def test_compare_annotation_status_1(self):
        """Verify no error is produced when both genomes have an
        expected annotation_status."""
        self.genome1.type = "mysql"
        self.genome1.annotation_status = "draft"
        self.genome2.type = "flat_file"
        self.genome2.annotation_status = "final"
        self.genome_pair.compare_annotation_status(
            "type", "mysql", "flat_file", "draft", "final", "eval_id")
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")

    def test_compare_annotation_status_2(self):
        """Verify no test is performed when the attribute is invalid."""
        self.genome1.type = "mysql"
        self.genome1.annotation_status = "draft"
        self.genome2.type = "flat_file"
        self.genome2.annotation_status = "final"
        self.genome_pair.compare_annotation_status(
            "type_x", "mysql", "flat_file", "draft", "final")
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].status, "untested")
        with self.subTest():
            self.assertIsNone(self.genome_pair.evaluations[0].id)

    def test_compare_annotation_status_3(self):
        """Verify no test is performed when the ref_name is invalid."""
        self.genome1.type = "mysql"
        self.genome1.annotation_status = "draft"
        self.genome2.type = "flat_file"
        self.genome2.annotation_status = "final"
        self.genome_pair.compare_annotation_status(
            "type", "mysql_x", "flat_file", "draft", "final")
        self.assertEqual(self.genome_pair.evaluations[0].status, "untested")

    def test_compare_annotation_status_4(self):
        """Verify no test is performed when the query_name is invalid."""
        self.genome1.type = "mysql"
        self.genome1.annotation_status = "draft"
        self.genome2.type = "flat_file"
        self.genome2.annotation_status = "final"
        self.genome_pair.compare_annotation_status(
            "type", "mysql", "flat_file_x", "draft", "final")
        self.assertEqual(self.genome_pair.evaluations[0].status, "untested")

    def test_compare_annotation_status_5(self):
        """Verify no test is performed when the ref_name and query_name
        are the same."""
        self.genome1.type = "mysql"
        self.genome1.annotation_status = "draft"
        self.genome2.type = "flat_file"
        self.genome2.annotation_status = "final"
        self.genome_pair.compare_annotation_status(
            "type", "mysql", "mysql", "draft", "final")
        self.assertEqual(self.genome_pair.evaluations[0].status, "untested")

    def test_compare_annotation_status_6(self):
        """Verify no error is produced when both genomes have an
        expected annotation_status with order reverse."""
        self.genome1.type = "mysql"
        self.genome1.annotation_status = "final"
        self.genome2.type = "flat_file"
        self.genome2.annotation_status = "draft"
        self.genome_pair.compare_annotation_status(
            "type", "flat_file", "mysql", "draft", "final")
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_annotation_status_7(self):
        """Verify an error is produced when genomes do not have the
        expected annotation_status."""
        self.genome1.type = "mysql"
        self.genome1.annotation_status = "final"
        self.genome2.type = "flat_file"
        self.genome2.annotation_status = "draft"
        self.genome_pair.compare_annotation_status(
            "type", "mysql", "flat_file", "draft", "final")
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")




    def test_compare_date_1(self):
        """Verify no error is produced when
        genome1 is older than genome2 as expected."""
        self.genome1.date = self.date_jan1
        self.genome2.date = self.date_feb1
        self.genome_pair.compare_date("older", "eval_id")
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")

    def test_compare_date_2(self):
        """Verify no error is produced when
        genome1 is newer than genome2 as expected."""
        self.genome1.date = self.date_feb1
        self.genome2.date = self.date_jan1
        self.genome_pair.compare_date("newer")
        with self.subTest():
            self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.genome_pair.evaluations[0].id)

    def test_compare_date_3(self):
        """Verify no error is produced when
        genome1 is equal to genome2 as expected."""
        self.genome1.date = self.date_feb1
        self.genome2.date = self.date_feb1_b
        self.genome_pair.compare_date("equal")
        self.assertEqual(self.genome_pair.evaluations[0].status, "correct")

    def test_compare_date_4(self):
        """Verify no error is produced when
        genome1 is older than genome2 unexpectedly."""
        self.genome1.date = self.date_jan1
        self.genome2.date = self.date_feb1
        self.genome_pair.compare_date("newer")
        self.assertEqual(self.genome_pair.evaluations[0].status, "error")

    def test_compare_date_5(self):
        """Verify no test is performed when an
        invalid comparison is selected."""
        self.genome1.date = self.date_feb1
        self.genome2.date = self.date_jan1
        self.genome_pair.compare_date("invalid")
        self.assertEqual(self.genome_pair.evaluations[0].status, "untested")






    # TODO these may no longer be needed. Older version of the method.
    # def test_compare_date_1(self):
    #     """Verify no error is produced when
    #     query is newer than ref as expected."""
    #     self.genome1.type = "mysql"
    #     self.genome1.date = self.date_jan1
    #     self.genome2.type = "flat_file"
    #     self.genome2.date = self.date_feb1
    #     self.genome_pair.compare_date(
    #         "type", "mysql", "flat_file", "newer", "eval_id")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")
    #
    # def test_compare_date_2(self):
    #     """Verify no error is produced when
    #     query is older than ref as expected."""
    #     self.genome1.type = "mysql"
    #     self.genome1.date = self.date_feb1
    #     self.genome2.type = "flat_file"
    #     self.genome2.date = self.date_jan1
    #     self.genome_pair.compare_date(
    #         "type", "mysql", "flat_file", "older", "eval_id")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")
    #
    # def test_compare_date_3(self):
    #     """Verify no error is produced when
    #     query is equal to ref as expected."""
    #     self.genome1.type = "mysql"
    #     self.genome1.date = self.date_feb1
    #     self.genome2.type = "flat_file"
    #     self.genome2.date = self.date_feb1_b
    #     self.genome_pair.compare_date(
    #         "type", "mysql", "flat_file", "equal", "eval_id")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")
    #
    # def test_compare_date_4(self):
    #     """Verify no error is produced when
    #     query and ref are reversed."""
    #     self.genome1.type = "mysql"
    #     self.genome1.date = self.date_feb1
    #     self.genome2.type = "flat_file"
    #     self.genome2.date = self.date_jan1
    #     self.genome_pair.compare_date(
    #         "type", "flat_file", "mysql", "newer", "eval_id")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "correct")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")
    #
    # def test_compare_date_5(self):
    #     """Verify an error is produced when
    #     query is older to ref unexpectedly."""
    #     self.genome1.type = "mysql"
    #     self.genome1.date = self.date_feb1
    #     self.genome2.type = "flat_file"
    #     self.genome2.date = self.date_jan1
    #     self.genome_pair.compare_date(
    #         "type", "flat_file", "mysql", "older", "eval_id")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "error")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].id, "eval_id")
    #
    # def test_compare_date_6(self):
    #     """Verify no test is performed when the attribute is invalid."""
    #     self.genome1.type = "mysql"
    #     self.genome1.date = self.date_feb1
    #     self.genome2.type = "flat_file"
    #     self.genome2.date = self.date_jan1
    #     self.genome_pair.compare_date(
    #         "invalid", "flat_file", "mysql", "older")
    #     with self.subTest():
    #         self.assertEqual(self.genome_pair.evaluations[0].status, "untested")
    #     with self.subTest():
    #         self.assertIsNone(self.genome_pair.evaluations[0].id)
    #
    # def test_compare_date_7(self):
    #     """Verify no test is performed when the ref_name is invalid."""
    #     self.genome1.type = "mysql"
    #     self.genome1.date = self.date_feb1
    #     self.genome2.type = "flat_file"
    #     self.genome2.date = self.date_jan1
    #     self.genome_pair.compare_date(
    #         "type", "invalid", "mysql", "older")
    #     self.assertEqual(self.genome_pair.evaluations[0].status, "untested")
    #
    # def test_compare_date_8(self):
    #     """Verify no test is performed when the query_name is invalid."""
    #     self.genome1.type = "mysql"
    #     self.genome1.date = self.date_feb1
    #     self.genome2.type = "flat_file"
    #     self.genome2.date = self.date_jan1
    #     self.genome_pair.compare_date(
    #         "type", "mysql", "invalid", "older")
    #     self.assertEqual(self.genome_pair.evaluations[0].status, "untested")
    #
    # def test_compare_date_9(self):
    #     """Verify no test is performed when the ref_name and query_name
    #     are the same."""
    #     self.genome1.type = "mysql"
    #     self.genome1.date = self.date_feb1
    #     self.genome2.type = "flat_file"
    #     self.genome2.date = self.date_jan1
    #     self.genome_pair.compare_date(
    #         "type", "mysql", "mysql", "older")
    #     self.assertEqual(self.genome_pair.evaluations[0].status, "untested")
    #
    # def test_compare_date_10(self):
    #     """Verify no test is performed when an
    #     invalid comparison is selected."""
    #     self.genome1.type = "mysql"
    #     self.genome1.date = self.date_feb1
    #     self.genome2.type = "flat_file"
    #     self.genome2.date = self.date_jan1
    #     self.genome_pair.compare_date(
    #         "type", "mysql", "flat_file", "invalid")
    #     self.assertEqual(self.genome_pair.evaluations[0].status, "untested")






if __name__ == '__main__':
    unittest.main()
