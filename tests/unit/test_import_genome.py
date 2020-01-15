""" Unit tests for import functions."""

from datetime import datetime
from pdm_utils.classes import bundle
from pdm_utils.classes import genome
from pdm_utils.classes import source
from pdm_utils.classes import cds
from pdm_utils.classes import genomepair
from pdm_utils.constants import constants
from pdm_utils.functions import run_modes
from pdm_utils.pipelines import import_genome
from pdm_utils.classes import ticket, eval
import unittest
from Bio.Seq import Seq
from unittest.mock import patch
from pdm_utils.classes import mysqlconnectionhandler as mch
from Bio.Alphabet import IUPAC

def get_errors(item):
    errors = 0
    for evl in item.evaluations:
        if evl.status == "error":
            errors += 1
    return errors


class TestImportGenomeClass1(unittest.TestCase):


    def setUp(self):

        # self.null_set = constants.EMPTY_SET
        self.type_set = constants.IMPORT_TICKET_TYPE_SET
        self.run_mode_set = run_modes.RUN_MODES.keys()
        self.description_field_set = constants.DESCRIPTION_FIELD_SET

        self.retain_set = constants.IMPORT_TABLE_STRUCTURE["valid_retain"]
        self.retrieve_set = constants.IMPORT_TABLE_STRUCTURE["valid_retrieve"]
        self.add_set = constants.IMPORT_TABLE_STRUCTURE["valid_add"]


        self.data_dict = {
            "host_genus": "Mycobacterium smegmatis",
            "cluster": "A",
            "subcluster": "A2",
            "annotation_status": "draft",
            "annotation_author": 1,
            "retrieve_record": 1,
            "accession": "ABC123.1"
            }

        self.tkt = ticket.GenomeTicket()
        self.tkt.id = 1
        self.tkt.type = "replace"
        self.tkt.phage_id = "Trixie"
        self.tkt.run_mode = "phagesdb"
        self.tkt.description_field = "product"
        self.tkt.data_dict = self.data_dict
        self.tkt.eval_flags = {"a":1, "b":2}
        self.tkt.data_retain = set(["host_genus"])
        self.tkt.data_retrieve = set(["cluster"])
        self.tkt.data_add = set(["retrieve_record"])


    def test_check_ticket_1(self):
        """Verify no error is produced with a correctly structured
        'replace' ticket."""
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            run_mode_set=self.run_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set)
        errors = get_errors(self.tkt)
        with self.subTest():
            self.assertEqual(len(self.tkt.evaluations), 11)
        with self.subTest():
            self.assertEqual(errors, 0)

    def test_check_ticket_2(self):
        """Verify correct number of errors is produced with
        a duplicate id."""
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            run_mode_set=self.run_mode_set,
            id_dupe_set=set([1]), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set)
        errors = get_errors(self.tkt)
        self.assertEqual(errors, 1)

    def test_check_ticket_3(self):
        """Verify correct number of errors is produced with
        a duplicate phage_id."""
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            run_mode_set=self.run_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(["Trixie"]),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set)
        errors = get_errors(self.tkt)
        self.assertEqual(errors, 1)

    def test_check_ticket_4(self):
        """Verify correct number of errors is produced with
        an invalid type."""
        self.tkt.type = "invalid"
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            run_mode_set=self.run_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set)
        errors = get_errors(self.tkt)
        self.assertEqual(errors, 1)

    def test_check_ticket_5(self):
        """Verify correct number of errors is produced with
        an invalid description_field."""
        self.tkt.description_field = "invalid"
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            run_mode_set=self.run_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set)
        errors = get_errors(self.tkt)
        self.assertEqual(errors, 1)

    def test_check_ticket_6(self):
        """Verify correct number of errors is produced with
        an invalid run_mode."""
        self.tkt.run_mode = "invalid"
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            run_mode_set=self.run_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set)
        errors = get_errors(self.tkt)
        self.assertEqual(errors, 1)

    def test_check_ticket_7(self):
        """Verify correct number of errors is produced with
        an invalid eval_flag dictionary."""
        self.tkt.eval_flags = {}
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            run_mode_set=self.run_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set)
        errors = get_errors(self.tkt)
        self.assertEqual(errors, 1)

    def test_check_ticket_8(self):
        """Verify correct number of errors is produced with
        an invalid phage_id."""
        self.tkt.phage_id = ""
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            run_mode_set=self.run_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set)
        errors = get_errors(self.tkt)
        self.assertEqual(errors, 1)

    def test_check_ticket_9(self):
        """Verify correct number of errors is produced with
        an incompatible "add" type and data_retain."""
        self.tkt.type = "add"
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            run_mode_set=self.run_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set)
        errors = get_errors(self.tkt)
        self.assertEqual(errors, 1)

    def test_check_ticket_10(self):
        """Verify correct number of errors is produced with
        an data in data_retain."""
        self.tkt.data_retain = set(["invalid"])
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            run_mode_set=self.run_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set)
        errors = get_errors(self.tkt)
        self.assertEqual(errors, 1)

    def test_check_ticket_11(self):
        """Verify correct number of errors is produced with
        an data in data_retrieve."""
        self.tkt.data_retrieve = set(["invalid"])
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            run_mode_set=self.run_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set)
        errors = get_errors(self.tkt)
        self.assertEqual(errors, 1)

    def test_check_ticket_12(self):
        """Verify correct number of errors is produced with
        an data in data_add."""
        self.tkt.data_add = set(["invalid"])
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            run_mode_set=self.run_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set)
        errors = get_errors(self.tkt)
        self.assertEqual(errors, 1)






class TestImportGenomeClass2(unittest.TestCase):


    def setUp(self):

        self.tkt = ticket.GenomeTicket()
        self.tkt.type = "add"
        self.tkt.eval_flags["check_seq"] = True
        self.tkt.eval_flags["check_id_typo"] = True
        self.tkt.eval_flags["check_host_typo"] = True
        self.tkt.eval_flags["check_author"] = True
        self.tkt.eval_flags["check_trna"] = True
        self.tkt.eval_flags["check_gene"] = True
        self.tkt.eval_flags["check_locus_tag"] = True
        self.tkt.eval_flags["check_description"] = True
        self.tkt.eval_flags["check_description_field"] = True

        self.cds1 = cds.Cds()
        self.cds2 = cds.Cds()
        self.source1 = source.Source()

        self.gnm = genome.Genome()
        self.gnm.id = "Trixie"
        self.gnm.name = "Trixie_Draft"
        self.gnm.host_genus = "Mycobacterium"
        self.gnm.cluster = "A"
        self.gnm.subcluster = "A2"
        self.gnm.accession = "ABC123"
        self.gnm.filename = "Trixie.gb"
        self.gnm.seq = "ATCG"
        self.gnm.annotation_status = "final"
        self.gnm.annotation_author = 1
        self.gnm.retrieve_record = 1
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.source_features = [self.source1]

        self.null_set = set([""])
        self.id_set = set(["Trixie"])
        self.seq_set = set(["AATTAA"])
        self.host_set = set(["Mycobacterium"])
        self.cluster_set = set(["A", "B"])
        self.subcluster_set = set(["A1, A2"])




class TestImportGenomeClass3(unittest.TestCase):

    def setUp(self):

        self.date_jan1 = datetime.strptime('1/1/2000', '%m/%d/%Y')
        self.date_feb1 = datetime.strptime('2/1/2000', '%m/%d/%Y')
        self.date_feb1_b = datetime.strptime('2/1/2000', '%m/%d/%Y')

        self.ff_gnm = genome.Genome()
        self.ff_gnm.type = "flat_file"
        self.ff_gnm.id = "Trixie"
        self.ff_gnm.name = "Trixie"
        self.ff_gnm.date = self.date_feb1
        self.ff_gnm.annotation_status = "final"
        self.ff_gnm.seq = Seq("AAAA", IUPAC.ambiguous_dna)
        self.ff_gnm.length = 4
        self.ff_gnm.cluster = "A"
        self.ff_gnm.subcluster = "A2"
        self.ff_gnm.accession = "ABC123"
        self.ff_gnm.host_genus = "Mycobacterium"
        self.ff_gnm.annotation_author = 1
        self.ff_gnm.retrieve_record = 1
        self.ff_gnm.translation_table = 11

        self.pmr_gnm = genome.Genome()
        self.pmr_gnm.type = "mysql"
        self.pmr_gnm.id = "Trixie"
        self.pmr_gnm.name = "Trixie_Draft"
        self.pmr_gnm.date = self.date_jan1
        self.pmr_gnm.annotation_status = "draft"
        self.pmr_gnm.seq = Seq("AAAA", IUPAC.ambiguous_dna)
        self.pmr_gnm.length = 4
        self.pmr_gnm.cluster = "A"
        self.pmr_gnm.subcluster = "A2"
        self.pmr_gnm.accession = "ABC123"
        self.pmr_gnm.host_genus = "Mycobacterium"
        self.pmr_gnm.annotation_author = 1
        self.pmr_gnm.retrieve_record = 1
        self.pmr_gnm.translation_table = 11

        self.genome_pair = genomepair.GenomePair()
        self.genome_pair.genome1 = self.ff_gnm
        self.genome_pair.genome2 = self.pmr_gnm

        self.eval_flags = {"check_replace": True}


    def test_compare_genomes_1(self):
        """Verify correct number of evaluations are produced and
        the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft'."""
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        with self.subTest():
            self.assertEqual(len(self.genome_pair.evaluations), 12)
        with self.subTest():
            self.assertEqual(errors, 0)

    def test_compare_genomes_2(self):
        """Verify correct number of evaluations are produced and
        the correct number of errors when:
        'check_replace' is True, annotation_status = 'final'."""
        self.pmr_gnm.annotation_status = "final"
        self.pmr_gnm.name = "Trixie"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        with self.subTest():
            self.assertEqual(len(self.genome_pair.evaluations), 13)
        with self.subTest():
            self.assertEqual(errors, 0)

    def test_compare_genomes_3(self):
        """Verify correct number of evaluations are produced and
        the correct number of errors when:
        'check_replace' is False."""
        self.eval_flags = {"check_replace": False}
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        with self.subTest():
            self.assertEqual(len(self.genome_pair.evaluations), 9)
        with self.subTest():
            self.assertEqual(errors, 0)

    def test_compare_genomes_4(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'id' values are different."""
        self.pmr_gnm.id = "L5"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)

    def test_compare_genomes_5(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'seq' values are different."""
        self.pmr_gnm.seq = Seq("AAAT", IUPAC.ambiguous_dna)
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)

    def test_compare_genomes_6(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'length' values are different."""
        self.pmr_gnm.length = 5
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)

    def test_compare_genomes_7(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'cluster' values are different."""
        self.pmr_gnm.cluster = "B"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)

    def test_compare_genomes_8(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'subcluster' values are different."""
        self.pmr_gnm.subcluster = "B2"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)

    def test_compare_genomes_9(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'accession' values are different."""
        self.pmr_gnm.accession = "XYZ456"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 0)

    def test_compare_genomes_10(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'host_genus' values are different."""
        self.pmr_gnm.host_genus = "Gordonia"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)

    def test_compare_genomes_11(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'annotation_author' values are different."""
        self.pmr_gnm.annotation_author = 0
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)

    def test_compare_genomes_12(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'translation_table' values are different."""
        self.pmr_gnm.translation_table = 1
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)

    def test_compare_genomes_13(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'retrieve_record' values are different."""
        self.pmr_gnm.retrieve_record = 0
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)

    def test_compare_genomes_14(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'date' values are not expected."""
        self.pmr_gnm.date = self.date_feb1
        self.ff_gnm.date = self.date_jan1
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)

    def test_compare_genomes_15(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'annotation_status' values are the same."""
        self.ff_gnm.annotation_status = "draft"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)

    def test_compare_genomes_16(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'name' values are the same."""
        self.ff_gnm.name = "Trixie_Draft"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)

    def test_compare_genomes_17(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'final',
        and 'annotation_status' values are not the same."""
        self.pmr_gnm.annotation_status = "final"
        self.pmr_gnm.name = "Trixie"
        self.ff_gnm.annotation_status = "unknown"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)

    def test_compare_genomes_18(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'final',
        and 'name' values are not the same."""
        self.pmr_gnm.annotation_status = "final"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        errors = get_errors(self.genome_pair)
        self.assertEqual(errors, 1)




class TestImportGenomeClass4(unittest.TestCase):

    def setUp(self):
        self.date_jan1 = datetime.strptime('1/1/2000', '%m/%d/%Y')

        self.eval_dict = {"check_replace": True,
                         "import_locus_tag": True,
                         "check_locus_tag": True,
                         "check_description_field": True,
                         "check_description": True,
                         "check_trna": True,
                         "check_id_typo": True,
                         "check_host_typo": True,
                         "check_author": True,
                         "check_gene": True,
                         "check_seq": True}

        self.tkt = ticket.GenomeTicket()
        self.tkt.type = "add"
        self.tkt.phage_id = "Trixie"
        self.tkt.eval_flags = self.eval_dict
        self.tkt.description_field = "product"

        self.cds1 = cds.Cds()
        self.cds1.id = "L5_1"
        self.cds1.start = 10
        self.cds1.stop = 20
        self.cds1.coordinate_format = "0_half_open"
        self.cds1.orientation = "F"

        self.cds2 = cds.Cds()
        self.cds2.id = "L5_2"
        self.cds2.start = 100
        self.cds2.stop = 200
        self.cds2.coordinate_format = "0_half_open"
        self.cds2.orientation = "F"

        self.gnm = genome.Genome()
        self.gnm.type = "flat_file"
        self.gnm.id = "Trixie"
        self.gnm.name = "Trixie_Draft"
        self.gnm.date = self.date_jan1
        self.gnm.annotation_status = "draft"
        self.gnm.seq = Seq("AAAA", IUPAC.ambiguous_dna)
        self.gnm.length = 4
        self.gnm.gc = 0
        self.gnm.cluster = "A"
        self.gnm.subcluster = "A2"
        self.gnm.accession = ""
        self.gnm.host_genus = "Mycobacterium"
        self.gnm.annotation_author = 1
        self.gnm.retrieve_record = 1
        self.gnm.translation_table = 11
        self.gnm._cds_processed_descriptions_tally = 0
        self.gnm._cds_features_tally = 2
        self.gnm.description = "Mycobacterium phage Trixie"
        self.gnm._description_name = "Trixie"
        self.gnm._description_host_genus = "Mycobacterium"
        self.gnm.source = "Mycobacterium phage Trixie"
        self.gnm._source_name = "Trixie"
        self.gnm._source_host_genus = "Mycobacterium"
        self.gnm.organism = "Mycobacterium phage Trixie"
        self.gnm._organism_name = "Trixie"
        self.gnm._organism_host_genus = "Mycobacterium"
        self.gnm.authors = "Doe,J., Hatfull,G.F., John;R., Smith; Jones,B."
        self.gnm.cds_features = [self.cds1, self.cds2]

        self.id_set = set(["L5", "RedRock"])
        self.seq_set = set([Seq("ATGC", IUPAC.ambiguous_dna),
                           Seq("TTTT", IUPAC.ambiguous_dna)])
        self.host_set = set(["Mycobacterium", "Gordonia"])
        self.cluster_set = set(["A", "B", "C"])
        self.subcluster_set = set(["A1", "A2", "A3"])
        self.accession_set = set(["AAAAA", "XYZ123", "JKL123"])

        self.check_sum = 34

    def test_check_genome_1(self):
        """Verify correct number of evaluations are produced using:
        ticket type = 'add',
        'annotation_status' = 'draft',
        all eval_flags = True,
        'accession' != '', and
        'annotation_author' = 1."""
        self.gnm.accession = "ABC123"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        self.assertEqual(len(self.gnm.evaluations), self.check_sum)

    def test_check_genome_2(self):
        """Verify correct number of evaluations are produced using:
        ticket type = 'add',
        'annotation_status' = 'draft',
        all eval_flags = True,
        'accession' == '', and
        'annotation_author' = 1."""
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        self.assertEqual(len(self.gnm.evaluations), self.check_sum - 1)

    def test_check_genome_3(self):
        """Verify correct number of evaluations are produced using:
        ticket type = 'replace',
        'annotation_status' = 'draft',
        all eval_flags = True,
        'accession' != '', and
        'annotation_author' = 1."""
        self.tkt.type = "replace"
        self.gnm.accession = "ABC123"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        self.assertEqual(len(self.gnm.evaluations), self.check_sum - 2)

    def test_check_genome_4(self):
        """Verify correct number of evaluations are produced using:
        ticket type = 'add',
        'annotation_status' = 'final',
        all eval_flags = True,
        'accession' != '', and
        'annotation_author' = 1."""
        self.gnm.annotation_status = "final"
        self.gnm.accession = "ABC123"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        self.assertEqual(len(self.gnm.evaluations), self.check_sum - 1)

    def test_check_genome_5(self):
        """Verify correct number of evaluations are produced using:
        ticket type = 'add',
        'annotation_status' = 'unknown',
        all eval_flags = True,
        'accession' != '', and
        'annotation_author' = 1."""
        self.gnm.annotation_status = "unknown"
        self.gnm.accession = "ABC123"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        self.assertEqual(len(self.gnm.evaluations), self.check_sum - 3)

    def test_check_genome_6(self):
        """Verify correct number of evaluations are produced using:
        ticket type = 'add',
        'annotation_status' = 'draft',
        'check_seq' = False,
        'accession' != '', and
        'annotation_author' = 1."""
        self.eval_dict["check_seq"] = False
        self.gnm.accession = "ABC123"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        self.assertEqual(len(self.gnm.evaluations), self.check_sum - 1)

    def test_check_genome_7(self):
        """Verify correct number of evaluations are produced using:
        ticket type = 'add',
        'annotation_status' = 'draft',
        'check_id_typo' = False,
        'accession' != '', and
        'annotation_author' = 1."""
        self.eval_dict["check_id_typo"] = False
        self.gnm.accession = "ABC123"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        self.assertEqual(len(self.gnm.evaluations), self.check_sum - 3)

    def test_check_genome_8(self):
        """Verify correct number of evaluations are produced using:
        ticket type = 'add',
        'annotation_status' = 'draft',
        'check_host_typo' = False,
        'accession' != '', and
        'annotation_author' = 1."""
        self.eval_dict["check_host_typo"] = False
        self.gnm.accession = "ABC123"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        self.assertEqual(len(self.gnm.evaluations), self.check_sum - 3)

    def test_check_genome_9(self):
        """Verify correct number of evaluations are produced using:
        ticket type = 'add',
        'annotation_status' = 'draft',
        'check_author' = False,
        'accession' != '', and
        'annotation_author' = 1."""
        self.eval_dict["check_author"] = False
        self.gnm.accession = "ABC123"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        self.assertEqual(len(self.gnm.evaluations), self.check_sum - 2)

    def test_check_genome_10(self):
        """Verify correct number of evaluations are produced using:
        ticket type = 'add',
        'annotation_status' = 'draft',
        all eval_flags = True,
        'accession' != '', and
        'annotation_author' = 0."""
        self.gnm.annotation_author = 0
        self.gnm.accession = "ABC123"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        self.assertEqual(len(self.gnm.evaluations), self.check_sum - 1)

    def test_check_genome_11(self):
        """Verify correct number of errors are produced using:
        ticket type = 'add',
        'annotation_status' = 'draft',
        all eval_flags = True,
        'accession' == '', and
        'annotation_author' = 1."""
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 0)

    def test_check_genome_12(self):
        """Verify correct number of errors are produced using:
        'add' ticket type and 'id' in id_set."""
        self.id_set.add("Trixie")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_13(self):
        """Verify correct number of errors are produced using:
        'add' ticket type and 'id' == ''."""
        self.gnm.id = ""
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 4)

    def test_check_genome_14(self):
        """Verify correct number of errors are produced using:
        'add' ticket type and 'name' in id_set."""
        self.id_set.add("Trixie_Draft")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_15(self):
        """Verify correct number of errors are produced using:
        'add' ticket type and 'name' == ''."""
        self.gnm.name = ""
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 2)

    def test_check_genome_16(self):
        """Verify correct number of errors are produced using:
        'add' ticket type and 'seq' in seq_set."""
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_17(self):
        """Verify correct number of errors are produced using:
        'add' ticket type and 'seq' == ''."""
        self.gnm.seq = Seq("", IUPAC.ambiguous_dna)
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_18(self):
        """Verify correct number of errors are produced using:
        'add' ticket type, 'accession' != '' and not in accession_set."""
        self.gnm.annotation_status = "unknown"
        self.gnm.accession = "BBBBB"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 0)

    def test_check_genome_19(self):
        """Verify correct number of errors are produced using:
        'add' ticket type, 'accession' != '' and in accession_set."""
        self.gnm.annotation_status = "unknown"
        self.gnm.accession = "AAAAA"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_20(self):
        """Verify correct number of errors are produced using:
        'replace' ticket type, and 'id' in id_set."""
        self.tkt.type = "replace"
        self.gnm.accession = "ABC123"
        self.gnm.annotation_status = "final"
        self.gnm._cds_processed_descriptions_tally = 1
        self.gnm.name = "Trixie"
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        self.id_set.add("Trixie")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 0)

    def test_check_genome_21(self):
        """Verify correct number of errors are produced using:
        'replace' ticket type, and 'id' not in id_set."""
        self.tkt.type = "replace"
        self.gnm.accession = "ABC123"
        self.gnm.annotation_status = "final"
        self.gnm._cds_processed_descriptions_tally = 1
        self.gnm.name = "L5"
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_22(self):
        """Verify correct number of errors are produced using:
        'replace' ticket type, and 'seq' not in seq_set."""
        self.tkt.type = "replace"
        self.gnm.name = "Trixie"
        self.gnm.accession = "ABC123"
        self.gnm.annotation_status = "final"
        self.gnm._cds_processed_descriptions_tally = 1
        self.gnm.seq = Seq("GGGG", IUPAC.ambiguous_dna)
        self.id_set.add("Trixie")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_23(self):
        """Verify correct number of errors are produced using:
        'draft' annotation_status, and
        'name' does not contain '_Draft' suffix."""
        self.gnm.name = "Trixie"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_24(self):
        """Verify correct number of errors are produced using:
        'draft' annotation_status, and
        '_cds_processed_descriptions_tally' > 0."""
        self.gnm._cds_processed_descriptions_tally = 1
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_25(self):
        """Verify correct number of errors are produced using:
        'draft' annotation_status, and 'accession' != ''."""
        self.gnm.accession = "ZZZZ"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_26(self):
        """Verify correct number of errors are produced using:
        'final' annotation_status, and
        'name' contains '_Draft' suffix."""
        self.tkt.type = "replace"
        self.gnm.accession = "ABC123"
        self.gnm.annotation_status = "final"
        self.gnm._cds_processed_descriptions_tally = 1
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        self.id_set.add("Trixie")
        self.id_set.add("Trixie_Draft")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_27(self):
        """Verify correct number of errors are produced using:
        'final' annotation_status, and
        '_cds_processed_descriptions_tally' == 0."""
        self.tkt.type = "replace"
        self.gnm.accession = "ABC123"
        self.gnm.annotation_status = "final"
        self.gnm._cds_processed_descriptions_tally = 0
        self.gnm.name = "Trixie"
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        self.id_set.add("Trixie")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_28(self):
        """Verify correct number of errors are produced using:
        ticket type = 'add', id contains '_Draft' suffix.
        'annotation_status' = 'draft',
        all eval_flags = True,
        'accession' == '', and
        'annotation_author' = 1."""
        self.gnm.id = "Trixie_Draft"
        self.gnm._description_name = "Trixie_Draft"
        self.gnm._source_name = "Trixie_Draft"
        self.gnm._organism_name = "Trixie_Draft"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_29(self):
        """Verify correct number of errors are produced using:
        invalid 'annotation_status'."""
        self.gnm.annotation_status = "invalid"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_30(self):
        """Verify correct number of errors are produced using:
        invalid 'annotation_author'."""
        self.gnm.annotation_author = -1
        self.gnm.authors = "none"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_31(self):
        """Verify correct number of errors are produced using:
        invalid 'retrieve_record'."""
        self.gnm.retrieve_record = -1
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_32(self):
        """Verify correct number of errors are produced using:
        'cluster' not in cluster_set."""
        self.gnm.cluster = "Z"
        self.gnm.subcluster = "Z1"
        self.subcluster_set.add("Z1")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_33(self):
        """Verify correct number of errors are produced using:
        'subcluster' not in subcluster_set."""
        self.gnm.cluster = "Z"
        self.gnm.subcluster = "Z1"
        self.cluster_set.add("Z")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_34(self):
        """Verify correct number of errors are produced using:
        'subcluster' = 'none'."""
        self.gnm.subcluster = "none"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 0)


    def test_check_genome_36(self):
        """Verify correct number of errors are produced using:
        invalid 'translation_table'."""
        self.gnm.translation_table = 1
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_37(self):
        """Verify correct number of errors are produced using:
        'host_genus' not in host_set."""
        self.host_set = {"Gordonia"}
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_38(self):
        """Verify correct number of errors are produced using:
        'cluster' structured incorrectly."""
        self.gnm.cluster = "Z1"
        self.gnm.subcluster = "Z1"
        self.cluster_set.add("Z1")
        self.subcluster_set.add("Z1")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 2)

    def test_check_genome_39(self):
        """Verify correct number of errors are produced using:
        'subcluster' structured incorrectly."""
        self.gnm.cluster = "Z"
        self.gnm.subcluster = "Z"
        self.cluster_set.add("Z")
        self.subcluster_set.add("Z")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 2)

    def test_check_genome_40(self):
        """Verify correct number of errors are produced using:
        incompatible 'cluster' and 'subcluster'."""
        self.gnm.cluster = "Z"
        self.gnm.subcluster = "X1"
        self.cluster_set.add("Z")
        self.subcluster_set.add("X1")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_41(self):
        """Verify correct number of errors are produced using:
        invalid 'date'."""
        self.gnm.date = constants.EMPTY_DATE
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_42(self):
        """Verify correct number of errors are produced using:
        'gc' < 0."""
        self.gnm.gc = -1
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_43(self):
        """Verify correct number of errors are produced using:
        'gc' > 100."""
        self.gnm.gc = 101
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_44(self):
        """Verify correct number of errors are produced using:
        'length' = 0."""
        self.gnm.length = 0
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_45(self):
        """Verify correct number of errors are produced using:
        '_cds_features_tally' = 0."""
        self.gnm._cds_features_tally = 0
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_46(self):
        """Verify correct number of errors are produced using:
        duplicated feature coordinates."""
        self.cds2.start = 10
        self.cds2.stop = 20
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_47(self):
        """Verify correct number of errors are produced using:
        invalid nucleotides."""
        self.gnm.seq = Seq("CCCCC-C", IUPAC.ambiguous_dna)
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_48(self):
        """Verify correct number of errors are produced using:
        incompatible 'id' and '_description_name'."""
        self.gnm._description_name = "L5"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_49(self):
        """Verify correct number of errors are produced using:
        incompatible 'id' and '_source_name'."""
        self.gnm._source_name = "L5"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_50(self):
        """Verify correct number of errors are produced using:
        incompatible 'id' and '_organism_name'."""
        self.gnm._organism_name = "L5"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_51(self):
        """Verify correct number of errors are produced using:
        incompatible 'host_genus' and '_description_host_genus'."""
        self.gnm._description_host_genus = "Gordonia"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_52(self):
        """Verify correct number of errors are produced using:
        incompatible 'host_genus' and '_source_host_genus'."""
        self.gnm._source_host_genus = "Gordonia"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_53(self):
        """Verify correct number of errors are produced using:
        incompatible 'host_genus' and '_organism_host_genus'."""
        self.gnm._organism_host_genus = "Gordonia"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_54(self):
        """Verify correct number of errors are produced using:
        'annotation_author' = 1 and missing author."""
        self.gnm.authors = "Doe,J., Hatful,G.F., John;R., Smith;."
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_55(self):
        """Verify correct number of errors are produced using:
        'annotation_author' = 1 and invalid 'LASTNAME 'author."""
        self.gnm.authors = "Doe,J., Hatfull,G.F., LASTNAME, John;R., Smith;."
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)

    def test_check_genome_56(self):
        """Verify correct number of errors are produced using:
        'annotation_author' = 0 and invalid author."""
        self.gnm.annotation_author = 0
        self.gnm.authors = "Doe,J., Hatfull,G.F., John;R., Smith;."
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        errors = get_errors(self.gnm)
        self.assertEqual(errors, 1)




class TestImportGenomeClass5(unittest.TestCase):

    def setUp(self):
        self.tkt = ticket.GenomeTicket()

        self.tkt.type = "replace"
        self.tkt.data_retrieve = set(["cluster"])
        self.tkt.data_retain = set(["subcluster"])
        self.tkt.data_add = set(["host_genus"])

        self.gnm1 = genome.Genome()
        self.gnm2 = genome.Genome()
        self.gnm2.annotation_status = "final"
        self.gnm3 = genome.Genome()
        self.gnm4 = genome.Genome()

        self.genome_pair = genomepair.GenomePair()

        self.bndl = bundle.Bundle()
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["ticket"] = self.gnm1
        self.bndl.genome_dict["flat_file"] = self.gnm2
        self.bndl.genome_dict["phagesdb"] = self.gnm3
        self.bndl.genome_dict["mysql"] = self.gnm4
        self.bndl.set_genome_pair(self.genome_pair, "flat_file", "mysql")




    def test_check_bundle_1(self):
        """Verify all check methods are called with no errors."""
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql")
        errors = get_errors(self.bndl)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 7)
        with self.subTest():
            self.assertEqual(errors, 0)

    def test_check_bundle_2(self):
        """Verify correct number of errors and check methods called
        with no ticket."""
        self.bndl.ticket = None
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql")
        errors = get_errors(self.bndl)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 1)
        with self.subTest():
            self.assertEqual(errors, 1)

    def test_check_bundle_3(self):
        """Verify correct number of errors and check methods called
        with incorrect flat_file key."""
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file_x",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql")
        errors = get_errors(self.bndl)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 6)
        with self.subTest():
            self.assertEqual(errors, 2)

    def test_check_bundle_4(self):
        """Verify correct number of errors and check methods called
        with incompatible "replace" ticket type and "draft" annotation status."""
        self.gnm2.annotation_status = "draft"
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql")
        errors = get_errors(self.bndl)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 7)
        with self.subTest():
            self.assertEqual(errors, 1)

    def test_check_bundle_5(self):
        """Verify correct number of errors and check methods called
        with no data in data_add."""
        self.tkt.data_add = set()
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql")
        errors = get_errors(self.bndl)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 6)
        with self.subTest():
            self.assertEqual(errors, 0)

    def test_check_bundle_6(self):
        """Verify correct number of errors and check methods called
        with incorrect ticket key."""
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket_x",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql")
        errors = get_errors(self.bndl)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 7)
        with self.subTest():
            self.assertEqual(errors, 1)

    def test_check_bundle_7(self):
        """Verify correct number of errors and check methods called
        with no data in data_retrieve."""
        self.tkt.data_retrieve = set()
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql")
        errors = get_errors(self.bndl)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 6)
        with self.subTest():
            self.assertEqual(errors, 0)

    def test_check_bundle_8(self):
        """Verify correct number of errors and check methods called
        with incorrect retrieve key."""
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb_x",
                                   retain_ref="mysql")
        errors = get_errors(self.bndl)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 7)
        with self.subTest():
            self.assertEqual(errors, 1)

    def test_check_bundle_9(self):
        """Verify correct number of errors and check methods called
        with "add" ticket type and "draft" annotation status."""
        self.tkt.type = "add"
        self.gnm2.annotation_status = "draft"
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql")
        errors = get_errors(self.bndl)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 5)
        with self.subTest():
            self.assertEqual(errors, 0)

    def test_check_bundle_10(self):
        """Verify correct number of errors and check methods called
        with incorrect retain key."""
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql_x")
        errors = get_errors(self.bndl)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 7)
        with self.subTest():
            self.assertEqual(errors, 2)




class TestImportGenomeClass6(unittest.TestCase):

    def setUp(self):
        self.eval_dict = {"check_replace": True,
                         "import_locus_tag": True,
                         "check_locus_tag": True,
                         "check_description_field": True,
                         "check_description": True,
                         "check_trna": True,
                         "check_id_typo": True,
                         "check_host_typo": True,
                         "check_author": True,
                         "check_gene": True,
                         "check_seq": True}

        self.tkt = ticket.GenomeTicket()
        self.tkt.phage_id = "Trixie"
        self.tkt.eval_flags = self.eval_dict
        self.tkt.description_field = "product"

        self.cds1 = cds.Cds()
        self.cds2 = cds.Cds()

        # TODO using a Cds object since tRNA object is not available yet.
        self.trna1 = cds.Cds()
        self.trna2 = cds.Cds()

        self.src1 = source.Source()
        self.src2 = source.Source()

        self.gnm1 = genome.Genome()
        self.gnm1.id = "Trixie"

        self.gnm2 = genome.Genome()
        self.gnm2.type = "mysql"
        self.gnm2.id = "Trixie"

        self.genome_pair = genomepair.GenomePair()
        self.genome_pair.genome1 = self.gnm1
        self.genome_pair.genome2 = self.gnm2

        self.bndl = bundle.Bundle()

        self.null_set = set(["", "none", None])
        self.accession_set = set(["ABC123", "XYZ456"])
        self.phage_id_set = set(["L5", "Trixie"])
        self.seq_set = set(["AATTGG", "ATGC"])
        self.host_genus_set = set(["Mycobacterium", "Gordonia"])
        self.cluster_set = set(["A", "B"])
        self.subcluster_set = set(["A2", "B2"])

        self.sql_handle = mch.MySQLConnectionHandler()



    def test_run_checks_1(self):
        """Verify run_checks works using a bundle with:
        no ticket, no "flat_file" genome."""
        import_genome.run_checks(
                self.bndl,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                file_ref="flat_file", ticket_ref="ticket",
                retrieve_ref="phagesdb", retain_ref="mysql")
        with self.subTest():
            self.assertTrue(len(self.bndl.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.gnm1.evaluations) == 0)


    def test_run_checks_2(self):
        """Verify run_checks works using a bundle with:
        'add' ticket, no "flat_file" genome."""
        self.tkt.type = "add"
        self.bndl.ticket = self.tkt
        import_genome.run_checks(
                self.bndl,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                file_ref="flat_file", ticket_ref="ticket",
                retrieve_ref="phagesdb", retain_ref="mysql")
        with self.subTest():
            self.assertTrue(len(self.bndl.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.gnm1.evaluations) == 0)


    def test_run_checks_3(self):
        """Verify run_checks works using a bundle with:
        'add' ticket, 'flat_file' genome with no features."""
        self.tkt.type = "add"
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["flat_file"] = self.gnm1
        import_genome.run_checks(
                self.bndl,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                file_ref="flat_file", ticket_ref="ticket",
                retrieve_ref="phagesdb", retain_ref="mysql")
        with self.subTest():
            self.assertTrue(len(self.bndl.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.gnm1.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.cds1.evaluations) == 0)
        with self.subTest():
            self.assertTrue(len(self.cds2.evaluations) == 0)
        with self.subTest():
            self.assertTrue(len(self.src1.evaluations) == 0)
        with self.subTest():
            self.assertTrue(len(self.src2.evaluations) == 0)


    def test_run_checks_4(self):
        """Verify run_checks works using a bundle with:
        'add' ticket, 'flat_file' genome with two CDS features,
        two Source features, and two tRNA features."""
        self.tkt.type = "add"
        self.bndl.ticket = self.tkt
        self.gnm1.cds_features = [self.cds1, self.cds2]
        self.gnm1.source_features = [self.src1, self.src2]
        self.gnm1.trna_features = [self.trna1, self.trna2]
        self.bndl.genome_dict["flat_file"] = self.gnm1
        import_genome.run_checks(
                self.bndl,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                file_ref="flat_file", ticket_ref="ticket",
                retrieve_ref="phagesdb", retain_ref="mysql")
        with self.subTest():
            self.assertTrue(len(self.bndl.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.gnm1.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.cds1.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.cds2.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.src1.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.src2.evaluations) > 0)


    def test_run_checks_5(self):
        """Verify run_checks works using a bundle with:
        'add' ticket, no genome, 'flat_file_mysql' genome_pair."""
        self.tkt.type = "add"
        self.bndl.ticket = self.tkt
        self.bndl.genome_pair_dict["flat_file_mysql"] = self.genome_pair
        import_genome.run_checks(
                self.bndl,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                file_ref="flat_file", ticket_ref="ticket",
                retrieve_ref="phagesdb", retain_ref="mysql")
        with self.subTest():
            self.assertTrue(len(self.bndl.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.genome_pair.evaluations) == 0)
        with self.subTest():
            self.assertTrue(len(self.gnm1.evaluations) == 0)


    def test_run_checks_6(self):
        """Verify run_checks works using a bundle with:
        'replace' ticket, no genome, 'flat_file_mysql' genome_pair."""
        self.tkt.type = "replace"
        self.bndl.ticket = self.tkt
        self.bndl.genome_pair_dict["flat_file_mysql"] = self.genome_pair
        import_genome.run_checks(
                self.bndl,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                file_ref="flat_file", ticket_ref="ticket",
                retrieve_ref="phagesdb", retain_ref="mysql")
        with self.subTest():
            self.assertTrue(len(self.bndl.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.genome_pair.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.gnm1.evaluations) == 0)


    def test_run_checks_7(self):
        """Verify run_checks works using a bundle with:
        'replace' ticket, no matching genome, no matching genome_pair."""
        self.tkt.type = "replace"
        self.bndl.ticket = self.tkt
        self.gnm1.cds_features = [self.cds1, self.cds2]
        self.gnm1.source_features = [self.src1, self.src2]
        self.gnm1.trna_features = [self.trna1, self.trna2]
        self.bndl.genome_dict["flat_file_x"] = self.gnm1
        self.bndl.genome_pair_dict["flat_file_mysql_x"] = self.genome_pair
        import_genome.run_checks(
                self.bndl,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                file_ref="flat_file", ticket_ref="ticket",
                retrieve_ref="phagesdb", retain_ref="mysql")
        with self.subTest():
            self.assertTrue(len(self.bndl.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.genome_pair.evaluations) == 0)
        with self.subTest():
            self.assertTrue(len(self.gnm1.evaluations) == 0)
        with self.subTest():
            self.assertTrue(len(self.cds1.evaluations) == 0)
        with self.subTest():
            self.assertTrue(len(self.cds2.evaluations) == 0)
        with self.subTest():
            self.assertTrue(len(self.src1.evaluations) == 0)
        with self.subTest():
            self.assertTrue(len(self.src2.evaluations) == 0)



    @patch("pdm_utils.classes.genome.Genome.clear_locus_tags")
    def test_import_into_db_1(self, clear_mock):
        """Verify import_into_db works using a bundle with:
        1 error, prod_run = True, import_locus_tag = True."""
        self.bndl._errors = 1
        result = import_genome.import_into_db(self.bndl, self.sql_handle,
                    gnm_key="", prod_run=True)
        with self.subTest():
            self.assertFalse(result)
        with self.subTest():
            self.assertFalse(clear_mock.called)


    def test_import_into_db_2(self):
        """Verify import_into_db works using a bundle with:
        0 errors, genome present, prod_run = False."""
        self.bndl._errors = 0
        self.tkt.type = "replace"
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["flat_file"] = self.gnm1
        result = import_genome.import_into_db(self.bndl, self.sql_handle,
                    gnm_key="flat_file", prod_run=False)
        self.assertTrue(result)


    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler.execute_transaction")
    def test_import_into_db_3(self, execute_transaction_mock):
        """Verify import_into_db works using a bundle with:
        0 errors, genome present, prod_run = True, execution = failed."""
        execute_transaction_mock.return_value = 1
        self.bndl._errors = 0
        self.tkt.type = "replace"
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["flat_file"] = self.gnm1
        result = import_genome.import_into_db(self.bndl, self.sql_handle,
                    gnm_key="flat_file", prod_run=True)
        self.assertFalse(result)


    @patch("pdm_utils.classes.mysqlconnectionhandler.MySQLConnectionHandler.execute_transaction")
    def test_import_into_db_4(self, execute_transaction_mock):
        """Verify import_into_db works using a bundle with:
        0 errors, genome present, prod_run = True, execution = successful."""
        execute_transaction_mock.return_value = 0
        self.bndl._errors = 0
        self.tkt.type = "replace"
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["flat_file"] = self.gnm1
        result = import_genome.import_into_db(self.bndl, self.sql_handle,
                    gnm_key="flat_file", prod_run=True)
        self.assertTrue(result)


    @patch("pdm_utils.classes.genome.Genome.clear_locus_tags")
    def test_import_into_db_5(self, clear_mock):
        """Verify import_into_db works using a bundle with:
        0 errors, genome present, prod_run = False, import_locus_tag = False."""
        self.bndl._errors = 0
        self.tkt.type = "replace"
        self.tkt.eval_flags["import_locus_tag"] = False
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["flat_file"] = self.gnm1
        result = import_genome.import_into_db(self.bndl, self.sql_handle,
                    gnm_key="flat_file", prod_run=False)
        self.assertTrue(clear_mock.called)



class TestImportGenomeClass7(unittest.TestCase):
    def setUp(self):
        self.evl1 = eval.Eval()
        self.evl1.id = "GNM0001"
        self.evl1.definition = "temp"
        self.evl1.status = "error"
        self.evl1.result = "Failed evaluation."

        self.evl2 = eval.Eval()
        self.evl2.id = "GNM0002"
        self.evl2.definition = "temp"
        self.evl2.status = "error"
        self.evl2.result = "Failed evaluation."

        self.evl3 = eval.Eval()
        self.evl3.id = "GNM0003"
        self.evl3.definition = "temp"
        self.evl3.status = "correct"
        self.evl3.result = "Failed evaluation."


    def test_log_evaluations(self):
        """Verify function executes."""
        evaluation_dict = {1:{"bundle": [self.evl1],
                              "ticket": [self.evl2]},
                           2:{"genome": [self.evl3]}}
        import_genome.log_evaluations(evaluation_dict)











class TestImportGenomeClass8(unittest.TestCase):


    def setUp(self):

        self.cds1 = cds.Cds()
        self.cds1.id = "L5_1"
        self.cds1.name = "1"
        self.cds1.translation = Seq("MF", IUPAC.protein)
        self.cds1.translation_length = 2
        self.cds1.seq = Seq("ATGTTTTGA", IUPAC.unambiguous_dna)
        self.cds1.translation_table = 11
        self.cds1.start = 10
        self.cds1.stop = 20
        self.cds1.coordinate_format = "0_half_open"
        self.cds1.orientation = "F"
        self.cds1.parts = 1
        self.cds1.length = 9
        self.cds1.genome_id = "L5"
        self.cds1.genome_length = 50000
        self.cds1.pham = 100
        self.cds1.description = "repressor protein"
        self.cds1.processed_description = "repressor"
        self.cds1.locus_tag = "SEA_L5_1"
        self.cds1._locus_tag_num = "1"
        self.cds1.gene = "1"
        self.cds1.product = "repressor protein"
        self.cds1.function = "hypothetical protein"
        self.cds1.note = "protein"
        self.cds1.processed_product = "repressor"
        self.cds1.processed_function = ""
        self.cds1.processed_note = ""
        self.cds1.type = "CDS"

        self.eval_flags = {"check_locus_tag": True,
                           "check_gene": True,
                           "check_description": True,
                           "check_description_field": True}

    def test_check_cds_1(self):
        """Verify correct number of evaluations are produced when
        none are False."""
        import_genome.check_cds(self.cds1, self.eval_flags)
        self.assertEqual(len(self.cds1.evaluations), 12)

    def test_check_cds_2(self):
        """Verify correct number of evaluations are produced when
        check_locus_tag = False."""
        self.eval_flags["check_locus_tag"] = False
        import_genome.check_cds(self.cds1, self.eval_flags)
        self.assertEqual(len(self.cds1.evaluations), 9)

    def test_check_cds_3(self):
        """Verify correct number of evaluations are produced when
        check_gene = False."""
        self.eval_flags["check_gene"] = False
        import_genome.check_cds(self.cds1, self.eval_flags)
        self.assertEqual(len(self.cds1.evaluations), 9)

    def test_check_cds_4(self):
        """Verify correct number of evaluations are produced when
        check_description = False."""
        self.eval_flags["check_description"] = False
        import_genome.check_cds(self.cds1, self.eval_flags)
        self.assertEqual(len(self.cds1.evaluations), 12)

    def test_check_cds_5(self):
        """Verify correct number of evaluations are produced when
        all False."""
        self.eval_flags["check_locus_tag"] = False
        self.eval_flags["check_gene"] = False
        self.eval_flags["check_description"] = False
        self.eval_flags["check_description_field"] = False
        import_genome.check_cds(self.cds1, self.eval_flags)
        self.assertEqual(len(self.cds1.evaluations), 6)

    def test_check_cds_6(self):
        """Verify correct number of errors with correct CDS feature."""
        import_genome.check_cds(self.cds1, self.eval_flags)
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 0)

    def test_check_cds_7(self):
        """Verify correct number of errors with incorrect amino acids."""
        self.cds1.translation = Seq("MB", IUPAC.protein)
        import_genome.check_cds(self.cds1, self.eval_flags)
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 2)

    def test_check_cds_8(self):
        """Verify correct number of errors with incorrect translation."""
        self.cds1.translation = Seq("MM", IUPAC.protein)
        import_genome.check_cds(self.cds1, self.eval_flags)
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 1)

    def test_check_cds_9(self):
        """Verify correct number of errors with missing translation."""
        self.cds1.translation_length = 0
        import_genome.check_cds(self.cds1, self.eval_flags)
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 2)

    def test_check_cds_10(self):
        """Verify correct number of errors with incorrect translation table."""
        self.cds1.translation_table = 1
        import_genome.check_cds(self.cds1, self.eval_flags)
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 1)

    def test_check_cds_11(self):
        """Verify correct number of errors with incorrect start coordinate."""
        self.cds1.start = -1
        import_genome.check_cds(self.cds1, self.eval_flags)
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 1)

    def test_check_cds_12(self):
        """Verify correct number of errors with incorrect orientation."""
        self.cds1.orientation = "f"
        import_genome.check_cds(self.cds1, self.eval_flags)
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 1)

    def test_check_cds_13(self):
        """Verify correct number of errors with missing locus_tag."""
        self.cds1.locus_tag = ""
        import_genome.check_cds(self.cds1, self.eval_flags)
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 2)

    def test_check_cds_14(self):
        """Verify correct number of errors with incorrect locus_tag."""
        self.cds1.locus_tag = "ABCXYZ"
        import_genome.check_cds(self.cds1, self.eval_flags)
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 1)

    def test_check_cds_15(self):
        """Verify correct number of errors with missing gene qualifier."""
        self.cds1.gene = ""
        import_genome.check_cds(self.cds1, self.eval_flags)
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 3)

    def test_check_cds_16(self):
        """Verify correct number of errors with non-integer gene qualifier."""
        self.cds1.gene = "A"
        import_genome.check_cds(self.cds1, self.eval_flags)
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 2)

    def test_check_cds_17(self):
        """Verify correct number of errors with non-integer gene qualifier."""
        self.cds1.gene = "11"
        import_genome.check_cds(self.cds1, self.eval_flags)
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 1)

    def test_check_cds_18(self):
        """Verify correct number of errors with non-matching integer in
        gene qualifier and locus_tag."""
        self.cds1.gene = "11"
        import_genome.check_cds(self.cds1, self.eval_flags)
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 1)

    def test_check_cds_19(self):
        """Verify correct number of errors with non-matching integer in
        gene qualifier and locus_tag."""
        import_genome.check_cds(self.cds1, self.eval_flags,
                                description_field = "function")
        errors = get_errors(self.cds1)
        self.assertEqual(errors, 1)



class TestImportGenomeClass9(unittest.TestCase):


    def setUp(self):

        self.src1 = source.Source()
        self.src1.id = 1
        self.src1.name = ""
        self.src1.type = "source"
        self.src1.start = 1
        self.src1.stop = 50000
        self.src1.organism = "Mycobacterium phage L5"
        self.src1.host = "Mycobacterium smegmatis mc1255"
        self.src1.lab_host = "Mycobacterium smegmatis mc1255"
        self.src1.genome_id = "L5"
        self.src1.genome_host_genus = "Mycobacterium"
        self.src1._organism_name = "L5"
        self.src1._organism_host_genus = "Mycobacterium"
        self.src1._host_host_genus = "Mycobacterium"
        self.src1._lab_host_host_genus = "Mycobacterium"

        self.eval_flags = {"check_id_typo": True,
                           "check_host_typo": True}

    def test_check_source_1(self):
        """Verify correct number of evaluations are produced when
        none are False."""
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        self.assertEqual(len(self.src1.evaluations), 4)

    def test_check_source_2(self):
        """Verify correct number of evaluations are produced when
        check_id_typo = False."""
        self.eval_flags["check_id_typo"] = False
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        self.assertEqual(len(self.src1.evaluations), 3)

    def test_check_source_3(self):
        """Verify correct number of evaluations are produced when
        check_host_typo = False."""
        self.eval_flags["check_host_typo"] = False
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        self.assertEqual(len(self.src1.evaluations), 1)

    def test_check_source_4(self):
        """Verify correct number of evaluations are produced when
        organism is empty."""
        self.src1.organism = ""
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        self.assertEqual(len(self.src1.evaluations), 3)

    def test_check_source_5(self):
        """Verify correct number of evaluations are produced when
        host is empty."""
        self.src1.host = ""
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        self.assertEqual(len(self.src1.evaluations), 3)

    def test_check_source_6(self):
        """Verify correct number of evaluations are produced when
        lab_host is empty."""
        self.src1.lab_host = ""
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        self.assertEqual(len(self.src1.evaluations), 3)

    def test_check_source_7(self):
        """Verify correct number of errors with incorrect organism name."""
        self.src1._organism_name = "Trixie"
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        errors = get_errors(self.src1)
        self.assertEqual(errors, 1)

    def test_check_source_8(self):
        """Verify correct number of errors with incorrect organism host genus."""
        self.src1._organism_host_genus = "Gordonia"
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        errors = get_errors(self.src1)
        self.assertEqual(errors, 1)

    def test_check_source_9(self):
        """Verify correct number of errors with incorrect host host genus."""
        self.src1._host_host_genus = "Gordonia"
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        errors = get_errors(self.src1)
        self.assertEqual(errors, 1)

    def test_check_source_10(self):
        """Verify correct number of errors with incorrect lab host host genus."""
        self.src1._lab_host_host_genus = "Gordonia"
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        errors = get_errors(self.src1)
        self.assertEqual(errors, 1)

    def test_check_source_11(self):
        """Verify correct number of errors with permissible synonym
        in organism host genus."""
        self.src1._organism_host_genus = "Mycobacterio"
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        errors = get_errors(self.src1)
        self.assertEqual(errors, 0)

    def test_check_source_12(self):
        """Verify correct number of errors with permissible synonym
        in organism host genus but not in lab host host genus."""
        self.src1._organism_host_genus = "Mycobacterio"
        self.src1._lab_host_host_genus = "Mycobacterio"
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        errors = get_errors(self.src1)
        self.assertEqual(errors, 1)











if __name__ == '__main__':
    unittest.main()
