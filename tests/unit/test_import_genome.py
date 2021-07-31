""" Unit tests for import functions."""

import pathlib
import unittest
from unittest.mock import patch

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from datetime import datetime
import sqlalchemy

from pdm_utils.classes import bundle
from pdm_utils.classes import genome
from pdm_utils.classes import source
from pdm_utils.classes import cds, trna, tmrna
from pdm_utils.classes import genomepair
from pdm_utils.classes import ticket
from pdm_utils.constants import constants
from pdm_utils.functions import eval_modes
from pdm_utils.pipelines import import_genome

def count_status(item, *args):
    count = 0
    status_set = set()
    for arg in args:
        status_set.add(arg)
    for evl in item.evaluations:
        if evl.status in status_set:
            count += 1
    return count


class TestImportGenome1(unittest.TestCase):


    def setUp(self):

        self.type_set = constants.IMPORT_TICKET_TYPE_SET
        self.eval_mode_set = eval_modes.EVAL_MODES.keys()
        self.description_field_set = constants.DESCRIPTION_FIELD_SET

        self.retain_set = constants.IMPORT_TABLE_STRUCTURE["valid_retain"]
        self.retrieve_set = constants.IMPORT_TABLE_STRUCTURE["valid_retrieve"]
        self.add_set = constants.IMPORT_TABLE_STRUCTURE["valid_add"]
        self.parse_set = constants.IMPORT_TABLE_STRUCTURE["valid_parse"]


        self.data_dict = {
            "host_genus": "Mycobacterium smegmatis",
            "cluster": "A",
            "subcluster": "A2",
            "annotation_status": "draft",
            "annotation_author": 1,
            "retrieve_record": 1,
            "accession": "ABC123.1"
            }

        self.tkt = ticket.ImportTicket()
        self.tkt.id = 1
        self.tkt.type = "replace"
        self.tkt.phage_id = "Trixie"
        self.tkt.eval_mode = "final"
        self.tkt.description_field = "product"
        self.tkt.data_dict = self.data_dict
        self.tkt.eval_flags = {"a":1, "b":2}
        self.tkt.data_retain = set(["host_genus"])
        self.tkt.data_retrieve = set(["cluster"])
        self.tkt.data_add = set(["retrieve_record"])
        self.tkt.data_parse = set(["accession"])


    @patch("builtins.print")
    def test_log_and_print_1(self, p_mock):
        """Verify print is called."""
        import_genome.log_and_print("empty", terminal=True)
        self.assertTrue(p_mock.called)

    @patch("builtins.print")
    def test_log_and_print_2(self, p_mock):
        """Verify print is not called."""
        import_genome.log_and_print("empty", terminal=False)
        self.assertFalse(p_mock.called)




    def test_get_result_string_1(self):
        """Verify string is constructed correctly."""
        attr_list = ["type", "phage_id", "eval_mode"]
        string = import_genome.get_result_string(self.tkt, attr_list)
        exp = "type: replace, phage_id: Trixie, eval_mode: final"
        self.assertEqual(string, exp)




    def test_check_ticket_1(self):
        """Verify no error is produced with a correctly structured
        'replace' ticket."""
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            eval_mode_set=self.eval_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set, parse_set=self.parse_set)
        count = count_status(self.tkt, "error")
        with self.subTest():
            self.assertEqual(len(self.tkt.evaluations), 12)
        with self.subTest():
            self.assertEqual(count, 0)

    def test_check_ticket_2(self):
        """Verify correct number of errors is produced with
        a duplicate id."""
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            eval_mode_set=self.eval_mode_set,
            id_dupe_set=set([1]), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set, parse_set=self.parse_set)
        count = count_status(self.tkt, "error")
        self.assertEqual(count, 1)

    def test_check_ticket_3(self):
        """Verify correct number of errors is produced with
        a duplicate phage_id."""
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            eval_mode_set=self.eval_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(["Trixie"]),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set, parse_set=self.parse_set)
        count = count_status(self.tkt, "error")
        self.assertEqual(count, 1)

    def test_check_ticket_4(self):
        """Verify correct number of errors is produced with
        an invalid type."""
        self.tkt.type = "invalid"
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            eval_mode_set=self.eval_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set, parse_set=self.parse_set)
        count = count_status(self.tkt, "error")
        self.assertEqual(count, 1)

    def test_check_ticket_5(self):
        """Verify correct number of errors is produced with
        an invalid description_field."""
        self.tkt.description_field = "invalid"
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            eval_mode_set=self.eval_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set, parse_set=self.parse_set)
        count = count_status(self.tkt, "error")
        self.assertEqual(count, 1)

    def test_check_ticket_6(self):
        """Verify correct number of errors is produced with
        an invalid eval_mode."""
        self.tkt.eval_mode = "invalid"
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            eval_mode_set=self.eval_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set, parse_set=self.parse_set)
        count = count_status(self.tkt, "error")
        self.assertEqual(count, 1)

    def test_check_ticket_7(self):
        """Verify correct number of errors is produced with
        an invalid eval_flag dictionary."""
        self.tkt.eval_flags = {}
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            eval_mode_set=self.eval_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set, parse_set=self.parse_set)
        count = count_status(self.tkt, "error")
        self.assertEqual(count, 1)

    def test_check_ticket_8(self):
        """Verify correct number of errors is produced with
        an invalid phage_id."""
        self.tkt.phage_id = ""
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            eval_mode_set=self.eval_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set, parse_set=self.parse_set)
        count = count_status(self.tkt, "error")
        self.assertEqual(count, 1)

    def test_check_ticket_9(self):
        """Verify correct number of errors is produced with
        an incompatible "add" type and data_retain."""
        self.tkt.type = "add"
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            eval_mode_set=self.eval_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set, parse_set=self.parse_set)
        count = count_status(self.tkt, "error")
        self.assertEqual(count, 1)

    def test_check_ticket_10(self):
        """Verify correct number of errors is produced with
        invalid data in data_retain."""
        self.tkt.data_retain = {"invalid"}
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            eval_mode_set=self.eval_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set, parse_set=self.parse_set)
        count = count_status(self.tkt, "error")
        self.assertEqual(count, 1)

    def test_check_ticket_11(self):
        """Verify correct number of errors is produced with
        invalid data in data_retrieve."""
        self.tkt.data_retrieve = {"invalid"}
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            eval_mode_set=self.eval_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set, parse_set=self.parse_set)
        count = count_status(self.tkt, "error")
        self.assertEqual(count, 1)

    def test_check_ticket_12(self):
        """Verify correct number of errors is produced with
        invalid data in data_add."""
        self.tkt.data_add = {"invalid"}
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            eval_mode_set=self.eval_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set, parse_set=self.parse_set)
        count = count_status(self.tkt, "error")
        self.assertEqual(count, 1)

    def test_check_ticket_13(self):
        """Verify correct number of errors is produced with
        invalid data in data_parse."""
        self.tkt.data_parse = {"invalid"}
        import_genome.check_ticket(
            self.tkt, type_set=self.type_set,
            description_field_set=self.description_field_set,
            eval_mode_set=self.eval_mode_set,
            id_dupe_set=set(), phage_id_dupe_set=set(),
            retain_set=self.retain_set, retrieve_set=self.retrieve_set,
            add_set=self.add_set, parse_set=self.parse_set)
        count = count_status(self.tkt, "error")
        self.assertEqual(count, 1)



class TestImportGenome2(unittest.TestCase):


    def setUp(self):

        self.tkt = ticket.ImportTicket()
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




class TestImportGenome3(unittest.TestCase):

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
        count = count_status(self.genome_pair, "error", "warning")
        with self.subTest():
            self.assertEqual(len(self.genome_pair.evaluations), 12)
        with self.subTest():
            self.assertEqual(count, 0)

    def test_compare_genomes_2(self):
        """Verify correct number of evaluations are produced and
        the correct number of errors when:
        'check_replace' is True, annotation_status = 'final'."""
        self.pmr_gnm.annotation_status = "final"
        self.pmr_gnm.name = "Trixie"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        with self.subTest():
            self.assertEqual(len(self.genome_pair.evaluations), 13)
        with self.subTest():
            self.assertEqual(count, 0)

    def test_compare_genomes_3(self):
        """Verify correct number of evaluations are produced and
        the correct number of errors when:
        'check_replace' is False."""
        self.eval_flags = {"check_replace": False}
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        with self.subTest():
            self.assertEqual(len(self.genome_pair.evaluations), 9)
        with self.subTest():
            self.assertEqual(count, 0)

    def test_compare_genomes_4(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'id' values are different."""
        self.pmr_gnm.id = "L5"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)

    def test_compare_genomes_5(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'seq' values are different."""
        self.pmr_gnm.seq = Seq("AAAT", IUPAC.ambiguous_dna)
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)

    def test_compare_genomes_6(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'length' values are different."""
        self.pmr_gnm.length = 5
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)

    def test_compare_genomes_7(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'cluster' values are different."""
        self.pmr_gnm.cluster = "B"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)

    def test_compare_genomes_8(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'subcluster' values are different."""
        self.pmr_gnm.subcluster = "B2"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)

    def test_compare_genomes_9(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'accession' values are different."""
        self.pmr_gnm.accession = "XYZ456"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 0)

    def test_compare_genomes_10(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'host_genus' values are different."""
        self.pmr_gnm.host_genus = "Gordonia"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)

    def test_compare_genomes_11(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'annotation_author' values are different."""
        self.pmr_gnm.annotation_author = 0
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)

    def test_compare_genomes_12(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'translation_table' values are different."""
        self.pmr_gnm.translation_table = 1
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)

    def test_compare_genomes_13(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'retrieve_record' values are different."""
        self.pmr_gnm.retrieve_record = 0
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)

    def test_compare_genomes_14(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'date' values are not expected."""
        self.pmr_gnm.date = self.date_feb1
        self.ff_gnm.date = self.date_jan1
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)

    def test_compare_genomes_15(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'annotation_status' values are the same."""
        self.ff_gnm.annotation_status = "draft"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)

    def test_compare_genomes_16(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'draft',
        and 'name' values are the same."""
        self.ff_gnm.name = "Trixie_Draft"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)

    def test_compare_genomes_17(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'final',
        and 'annotation_status' values are not the same."""
        self.pmr_gnm.annotation_status = "final"
        self.pmr_gnm.name = "Trixie"
        self.ff_gnm.annotation_status = "unknown"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)

    def test_compare_genomes_18(self):
        """Verify the correct number of errors when:
        'check_replace' is True, annotation_status = 'final',
        and 'name' values are not the same."""
        self.pmr_gnm.annotation_status = "final"
        import_genome.compare_genomes(self.genome_pair, self.eval_flags)
        count = count_status(self.genome_pair, "error", "warning")
        self.assertEqual(count, 1)




class TestImportGenome4(unittest.TestCase):

    def setUp(self):
        self.date_jan1 = datetime.strptime('1/1/2000', '%m/%d/%Y')

        # Evaluation dict with all flags = True.
        self.eval_dict = eval_modes.get_eval_flag_dict("base")

        self.tkt = ticket.ImportTicket()
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

        self.trna1 = trna.Trna()
        self.trna1.id = "L5_3"
        self.trna1.start = 300
        self.trna1.stop = 400
        self.trna1.coordinate_format = "0_half_open"
        self.trna1.orientation = "F"

        self.trna2 = trna.Trna()
        self.trna2.id = "L5_4"
        self.trna2.start = 500
        self.trna2.stop = 600
        self.trna2.coordinate_format = "0_half_open"
        self.trna2.orientation = "F"

        self.tmrna1 = tmrna.Tmrna()
        self.tmrna1.id = "L5_4"
        self.tmrna1.start = 700
        self.tmrna1.stop = 800
        self.tmrna1.coordinate_format = "0_half_open"
        self.tmrna1.orientation = "F"

        self.tmrna2 = tmrna.Tmrna()
        self.tmrna2.id = "L5_5"
        self.tmrna2.start = 900
        self.tmrna2.stop = 1000
        self.tmrna2.coordinate_format = "0_half_open"
        self.tmrna2.orientation = "F"

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
        self.gnm._cds_descriptions_tally = 0
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
        self.gnm.trna_features = [self.trna1, self.trna2]
        self.gnm.tmrna_features = [self.tmrna1, self.tmrna2]

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
        'annotation_status' = 'final',
        'check_description_tally' = False,
        'accession' != '', and
        'annotation_author' = 1."""
        self.eval_dict["check_description_tally"] = False
        self.gnm.annotation_status = "final"
        self.gnm.accession = "ABC123"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        self.assertEqual(len(self.gnm.evaluations), self.check_sum - 2)

    def test_check_genome_6(self):
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

    def test_check_genome_7(self):
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

    def test_check_genome_8(self):
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

    def test_check_genome_9(self):
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

    def test_check_genome_10(self):
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

    def test_check_genome_11(self):
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

    def test_check_genome_12(self):
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
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 0)

    def test_check_genome_13(self):
        """Verify correct number of errors are produced using:
        'add' ticket type and 'id' in id_set."""
        self.id_set.add("Trixie")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_14(self):
        """Verify correct number of errors are produced using:
        'add' ticket type and 'id' == ''."""
        self.gnm.id = ""
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 4)

    def test_check_genome_15(self):
        """Verify correct number of errors are produced using:
        'add' ticket type and 'name' in id_set."""
        self.id_set.add("Trixie_Draft")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_16(self):
        """Verify correct number of errors are produced using:
        'add' ticket type and 'name' == ''."""
        self.gnm.name = ""
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 2)

    def test_check_genome_17(self):
        """Verify correct number of errors are produced using:
        'add' ticket type and 'seq' in seq_set."""
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_18(self):
        """Verify correct number of errors are produced using:
        'add' ticket type and 'seq' == ''."""
        self.gnm.seq = Seq("", IUPAC.ambiguous_dna)
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_19(self):
        """Verify correct number of errors are produced using:
        'add' ticket type and 'annotation_status' not valid."""
        self.gnm.annotation_status = "final"
        # Since status is now 'final', the cds tally and name need to be
        # changed as well, so that only the status check will throw an error.
        self.gnm._cds_descriptions_tally = 1
        self.gnm.name = "Trixie"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_20(self):
        """Verify correct number of errors are produced using:
        'add' ticket type, 'accession' != '' and not in accession_set."""
        self.gnm.annotation_status = "unknown"
        self.gnm.accession = "BBBBB"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 0)

    def test_check_genome_21(self):
        """Verify correct number of errors are produced using:
        'add' ticket type, 'accession' != '' and in accession_set."""
        self.gnm.annotation_status = "unknown"
        self.gnm.accession = "AAAAA"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_22(self):
        """Verify correct number of errors are produced using:
        'replace' ticket type, and 'id' in id_set."""
        self.tkt.type = "replace"
        self.gnm.accession = "ABC123"
        self.gnm.annotation_status = "final"
        self.gnm._cds_descriptions_tally = 1
        self.gnm.name = "Trixie"
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        self.id_set.add("Trixie")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 0)

    def test_check_genome_23(self):
        """Verify correct number of errors are produced using:
        'replace' ticket type, and 'id' not in id_set."""
        self.tkt.type = "replace"
        self.gnm.accession = "ABC123"
        self.gnm.annotation_status = "final"
        self.gnm._cds_descriptions_tally = 1
        self.gnm.name = "L5"
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_24(self):
        """Verify correct number of errors are produced using:
        'replace' ticket type, and 'seq' not in seq_set."""
        self.tkt.type = "replace"
        self.gnm.name = "Trixie"
        self.gnm.accession = "ABC123"
        self.gnm.annotation_status = "final"
        self.gnm._cds_descriptions_tally = 1
        self.gnm.seq = Seq("GGGG", IUPAC.ambiguous_dna)
        self.id_set.add("Trixie")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_25(self):
        """Verify correct number of errors are produced using:
        'replace' ticket type, and 'annotation_status' not valid."""
        self.gnm.annotation_status = "draft"
        # Since status is now 'draft', several other attributes need to be
        # changed as well, so that only the status check will throw an error.
        self.tkt.type = "replace"
        self.gnm.accession = ""
        self.gnm._cds_descriptions_tally = 0
        self.gnm.name = "Trixie_Draft"
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        self.id_set.add("Trixie")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_26(self):
        """Verify correct number of errors are produced using:
        'draft' annotation_status, and
        'name' does not contain '_Draft' suffix."""
        self.gnm.name = "Trixie"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_27(self):
        """Verify correct number of errors are produced using:
        'draft' annotation_status, and
        '_cds_descriptions_tally' > 0."""
        self.gnm._cds_descriptions_tally = 1
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_28(self):
        """Verify correct number of errors are produced using:
        'draft' annotation_status, and 'accession' != ''."""
        self.gnm.accession = "ZZZZ"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_29(self):
        """Verify correct number of errors are produced using:
        'final' annotation_status, and
        duplicated CDS coordinates."""
        self.tkt.type = "replace"
        self.gnm.accession = "ABC123"
        self.gnm.annotation_status = "final"
        self.gnm._cds_descriptions_tally = 1
        self.gnm.name = "Trixie"
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        self.id_set.add("Trixie")
        self.cds1.start = 10
        self.cds1.stop = 20
        self.cds1.orientation = "F"
        self.cds2.start = 10
        self.cds2.stop = 20
        self.cds2.orientation = "R"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_30(self):
        """Verify correct number of errors are produced using:
        'final' annotation_status, and
        duplicated CDS and tRNA coordinates."""
        self.tkt.type = "replace"
        self.gnm.accession = "ABC123"
        self.gnm.annotation_status = "final"
        self.gnm._cds_descriptions_tally = 1
        self.gnm.name = "Trixie"
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        self.id_set.add("Trixie")
        self.cds1.start = 10
        self.cds1.stop = 20
        self.cds1.orientation = "F"
        self.trna1.start = 10
        self.trna1.stop = 20
        self.trna1.orientation = "R"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_31(self):
        """Verify correct number of errors are produced using:
        'final' annotation_status, and
        duplicated tRNA and tmRNA coordinates."""
        self.tkt.type = "replace"
        self.gnm.accession = "ABC123"
        self.gnm.annotation_status = "final"
        self.gnm._cds_descriptions_tally = 1
        self.gnm.name = "Trixie"
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        self.id_set.add("Trixie")
        self.tmrna1.start = 10
        self.tmrna1.stop = 20
        self.tmrna1.orientation = "F"
        self.trna1.start = 10
        self.trna1.stop = 20
        self.trna1.orientation = "R"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_32(self):
        """Verify correct number of errors are produced using:
        'final' annotation_status, and
        'name' contains '_Draft' suffix."""
        self.tkt.type = "replace"
        self.gnm.accession = "ABC123"
        self.gnm.annotation_status = "final"
        self.gnm._cds_descriptions_tally = 1
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        self.id_set.add("Trixie")
        self.id_set.add("Trixie_Draft")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_33(self):
        """Verify correct number of errors are produced using:
        'final' annotation_status, and
        '_cds_descriptions_tally' == 0."""
        self.tkt.type = "replace"
        self.gnm.accession = "ABC123"
        self.gnm.annotation_status = "final"
        self.gnm._cds_descriptions_tally = 0
        self.gnm.name = "Trixie"
        self.gnm.seq = Seq("ATGC", IUPAC.ambiguous_dna)
        self.id_set.add("Trixie")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_34(self):
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
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_35(self):
        """Verify correct number of errors are produced using:
        invalid 'annotation_status'."""
        self.gnm.annotation_status = "invalid"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_36(self):
        """Verify correct number of errors are produced using:
        invalid 'annotation_author'."""
        self.gnm.annotation_author = -1
        self.gnm.authors = "none"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_37(self):
        """Verify correct number of errors are produced using:
        invalid 'retrieve_record'."""
        self.gnm.retrieve_record = -1
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_38(self):
        """Verify correct number of errors are produced using:
        'cluster' not in cluster_set."""
        self.gnm.cluster = "Z"
        self.gnm.subcluster = "Z1"
        self.subcluster_set.add("Z1")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_39(self):
        """Verify correct number of errors are produced using:
        'subcluster' not in subcluster_set."""
        self.gnm.cluster = "Z"
        self.gnm.subcluster = "Z1"
        self.cluster_set.add("Z")
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_40(self):
        """Verify correct number of errors are produced using:
        'subcluster' = 'none'."""
        self.gnm.subcluster = "none"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 0)

    def test_check_genome_41(self):
        """Verify correct number of errors are produced using:
        invalid 'translation_table'."""
        self.gnm.translation_table = 1
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_42(self):
        """Verify correct number of errors are produced using:
        'host_genus' not in host_set."""
        self.host_set = {"Gordonia"}
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_43(self):
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
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 2)

    def test_check_genome_44(self):
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
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 2)

    def test_check_genome_45(self):
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
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_46(self):
        """Verify correct number of errors are produced using:
        invalid 'date'."""
        self.gnm.date = constants.EMPTY_DATE
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_47(self):
        """Verify correct number of errors are produced using:
        'gc' < 0."""
        self.gnm.gc = -1
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_48(self):
        """Verify correct number of errors are produced using:
        'gc' > 100."""
        self.gnm.gc = 101
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_49(self):
        """Verify correct number of errors are produced using:
        'length' = 0."""
        self.gnm.length = 0
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_50(self):
        """Verify correct number of errors are produced using:
        '_cds_features_tally' = 0."""
        self.gnm._cds_features_tally = 0
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_51(self):
        """Verify correct number of errors are produced using:
        invalid nucleotides."""
        self.gnm.seq = Seq("CCCCC-C", IUPAC.ambiguous_dna)
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_52(self):
        """Verify correct number of errors are produced using:
        incompatible 'id' and '_description_name'."""
        self.gnm._description_name = "L5"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_53(self):
        """Verify correct number of errors are produced using:
        incompatible 'id' and '_source_name'."""
        self.gnm._source_name = "L5"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_54(self):
        """Verify correct number of errors are produced using:
        incompatible 'id' and '_organism_name'."""
        self.gnm._organism_name = "L5"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_55(self):
        """Verify correct number of errors are produced using:
        incompatible 'host_genus' and '_description_host_genus'."""
        self.gnm._description_host_genus = "Gordonia"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_56(self):
        """Verify correct number of errors are produced using:
        incompatible 'host_genus' and '_source_host_genus'."""
        self.gnm._source_host_genus = "Gordonia"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_57(self):
        """Verify correct number of errors are produced using:
        incompatible 'host_genus' and '_organism_host_genus'."""
        self.gnm._organism_host_genus = "Gordonia"
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_58(self):
        """Verify correct number of errors are produced using:
        'annotation_author' = 1 and missing author."""
        self.gnm.authors = "Doe,J., Hatful,G.F., John;R., Smith;."
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_59(self):
        """Verify correct number of errors are produced using:
        'annotation_author' = 1 and invalid 'LASTNAME 'author."""
        self.gnm.authors = "Doe,J., Hatfull,G.F., LASTNAME, John;R., Smith;."
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_genome_60(self):
        """Verify correct number of errors are produced using:
        'annotation_author' = 0 and invalid author."""
        self.gnm.annotation_author = 0
        self.gnm.authors = "Doe,J., Hatfull,G.F., John;R., Smith;."
        import_genome.check_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set, self.accession_set)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)




    def test_check_retain_genome_1(self):
        """Verify correct number of evaluations are produced using:
        all eval_flags = True."""
        import_genome.check_retain_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags)
        self.assertEqual(len(self.gnm.evaluations), 1)

    def test_check_retain_genome_2(self):
        """Verify correct number of evaluations are produced using:
        check_replace = False."""
        self.eval_dict["check_replace"] = False
        import_genome.check_retain_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags)
        self.assertEqual(len(self.gnm.evaluations), 0)

    def test_check_retain_genome_3(self):
        """Verify correct number of errors are produced using:
        all eval_flags = True,
        'annotation_status' = 'draft'."""
        import_genome.check_retain_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 0)

    def test_check_retain_genome_4(self):
        """Verify correct number of errors are produced using:
        all eval_flags = True,
        'annotation_status' = 'final'."""
        self.gnm.annotation_status = "final"
        import_genome.check_retain_genome(
            self.gnm, self.tkt.type, self.tkt.eval_flags)
        count = count_status(self.gnm, "error", "warning")
        self.assertEqual(count, 1)



class TestImportGenome5(unittest.TestCase):

    def setUp(self):
        self.tkt = ticket.ImportTicket()

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
        count = count_status(self.bndl, "error", "warning")
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 6)
        with self.subTest():
            self.assertEqual(count, 0)

    def test_check_bundle_2(self):
        """Verify correct number of errors and check methods called
        with no ticket."""
        self.bndl.ticket = None
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql")
        count = count_status(self.bndl, "error", "warning")
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 1)
        with self.subTest():
            self.assertEqual(count, 1)

    def test_check_bundle_3(self):
        """Verify correct number of errors and check methods called
        with incorrect flat_file key."""
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file_x",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql")
        count = count_status(self.bndl, "error", "warning")
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 6)
        with self.subTest():
            self.assertEqual(count, 2)

    def test_check_bundle_4(self):
        """Verify correct number of errors and check methods called
        with no data in data_add."""
        self.tkt.data_add = set()
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql")
        count = count_status(self.bndl, "error", "warning")
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 5)
        with self.subTest():
            self.assertEqual(count, 0)

    def test_check_bundle_5(self):
        """Verify correct number of errors and check methods called
        with incorrect ticket key."""
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket_x",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql")
        count = count_status(self.bndl, "error", "warning")
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 6)
        with self.subTest():
            self.assertEqual(count, 1)

    def test_check_bundle_6(self):
        """Verify correct number of errors and check methods called
        with no data in data_retrieve."""
        self.tkt.data_retrieve = set()
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql")
        count = count_status(self.bndl, "error", "warning")
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 5)
        with self.subTest():
            self.assertEqual(count, 0)

    def test_check_bundle_7(self):
        """Verify correct number of errors and check methods called
        with incorrect retrieve key."""
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb_x",
                                   retain_ref="mysql")
        count = count_status(self.bndl, "error", "warning")
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 6)
        with self.subTest():
            self.assertEqual(count, 1)

    def test_check_bundle_8(self):
        """Verify correct number of errors and check methods called
        with incorrect retain key."""
        import_genome.check_bundle(self.bndl,
                                   ticket_ref="ticket",
                                   file_ref="flat_file",
                                   retrieve_ref="phagesdb",
                                   retain_ref="mysql_x")
        count = count_status(self.bndl, "error", "warning")
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 6)
        with self.subTest():
            self.assertEqual(count, 2)




class TestImportGenome6(unittest.TestCase):

    def setUp(self):
        # Evaluation dict with all flags = True.
        self.eval_dict = eval_modes.get_eval_flag_dict("base")

        self.tkt = ticket.ImportTicket()
        self.tkt.phage_id = "Trixie"
        self.tkt.eval_flags = self.eval_dict
        self.tkt.description_field = "product"

        self.cds1 = cds.Cds()
        self.cds2 = cds.Cds()

        # Replaced trna1 and trna2 with TrnaFeatures instead of Cds
        self.trna1 = trna.Trna()
        self.trna2 = trna.Trna()

        self.tmrna1 = tmrna.Tmrna()
        self.tmrna2 = tmrna.Tmrna()

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

        engine_string = ("mysql+pymysql://"
                         "pdm_anon:pdm_anon@localhost/"
                         "actino_draft")
        self.engine = sqlalchemy.create_engine(engine_string, echo=True)


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
        two Source features, two tRNA features, and two tmRNA features."""
        self.tkt.type = "add"
        self.bndl.ticket = self.tkt
        self.gnm1.cds_features = [self.cds1, self.cds2]
        self.gnm1.source_features = [self.src1, self.src2]
        self.gnm1.trna_features = [self.trna1, self.trna2]
        self.gnm1.tmrna_features = [self.tmrna1, self.tmrna2]
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
        with self.subTest():
            self.assertTrue(len(self.trna1.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.trna1.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.tmrna1.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.tmrna2.evaluations) > 0)


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


    def test_run_checks_8(self):
        """Verify run_checks works using a bundle with:
        'replace' ticket, no file_ref genome, no matching genome_pair,
        mysql genome."""
        self.bndl.genome_dict["mysql"] = self.gnm1
        self.tkt.type = "replace"
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
            self.assertTrue(len(self.gnm1.evaluations) > 0)




    @patch("pdm_utils.functions.mysqldb.execute_transaction")
    @patch("pdm_utils.classes.genome.Genome.clear_locus_tags")
    def test_import_into_db_1(self, clear_mock, execute_mock):
        """Verify import_into_db works using a bundle with:
        1 error, prod_run = True, import_locus_tag = True."""
        self.bndl._errors = 1
        result = import_genome.import_into_db(self.bndl, self.engine,
                    gnm_key="", prod_run=True)
        with self.subTest():
            self.assertFalse(result)
        with self.subTest():
            self.assertFalse(clear_mock.called)
        with self.subTest():
            self.assertFalse(execute_mock.called)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 0)


    @patch("pdm_utils.functions.mysqldb.execute_transaction")
    def test_import_into_db_2(self, execute_mock):
        """Verify import_into_db works using a bundle with:
        0 errors, genome present, prod_run = False."""
        self.bndl._errors = 0
        self.tkt.type = "replace"
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["flat_file"] = self.gnm1
        result = import_genome.import_into_db(self.bndl, self.engine,
                    gnm_key="flat_file", prod_run=False)
        with self.subTest():
            self.assertTrue(result)
        with self.subTest():
            self.assertFalse(execute_mock.called)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 0)


    @patch("pdm_utils.functions.mysqldb.execute_transaction")
    def test_import_into_db_3(self, execute_mock):
        """Verify import_into_db works using a bundle with:
        0 errors, genome present, prod_run = True, execution = failed."""
        execute_mock.return_value = (1, "Fail")
        self.bndl._errors = 0
        self.tkt.type = "replace"
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["flat_file"] = self.gnm1
        result = import_genome.import_into_db(self.bndl, self.engine,
                    gnm_key="flat_file", prod_run=True)
        with self.subTest():
            self.assertFalse(result)
        with self.subTest():
            self.assertTrue(execute_mock.called)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 1)


    @patch("pdm_utils.functions.mysqldb.execute_transaction")
    def test_import_into_db_4(self, execute_mock):
        """Verify import_into_db works using a bundle with:
        0 errors, genome present, prod_run = True, execution = successful."""
        execute_mock.return_value = (0, "Success")
        self.bndl._errors = 0
        self.tkt.type = "replace"
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["flat_file"] = self.gnm1
        result = import_genome.import_into_db(self.bndl, self.engine,
                    gnm_key="flat_file", prod_run=True)
        with self.subTest():
            self.assertTrue(result)
        with self.subTest():
            self.assertTrue(execute_mock.called)
        with self.subTest():
            self.assertEqual(len(self.bndl.evaluations), 1)


    @patch("pdm_utils.classes.genome.Genome.clear_locus_tags")
    def test_import_into_db_5(self, clear_mock):
        """Verify import_into_db works using a bundle with:
        0 errors, genome present, prod_run = False, import_locus_tag = False."""
        self.bndl._errors = 0
        self.tkt.type = "replace"
        self.tkt.eval_flags["import_locus_tag"] = False
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["flat_file"] = self.gnm1
        result = import_genome.import_into_db(self.bndl, self.engine,
                    gnm_key="flat_file", prod_run=False)
        self.assertTrue(clear_mock.called)




class TestImportGenome7(unittest.TestCase):


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
        self.cds1.pham_id = 100
        self.cds1.raw_description = "repressor protein"
        self.cds1.description = "repressor"
        self.cds1.locus_tag = "SEA_L5_1"
        self.cds1._locus_tag_num = "1"
        self.cds1.gene = "1"
        self.cds1.raw_product = "repressor protein"
        self.cds1.raw_function = "hypothetical protein"
        self.cds1.raw_note = "protein"
        self.cds1.product = "repressor"
        self.cds1.function = ""
        self.cds1.note = ""
        self.cds1.type = "CDS"

        self.eval_flags = {"check_locus_tag": True,
                           "check_gene": True,
                           "check_description": True,
                           "check_description_field": True}

    def test_check_cds_1(self):
        """Verify correct number of evaluations are produced when
        none are False."""
        import_genome.check_cds(self.cds1, self.eval_flags)
        self.assertEqual(len(self.cds1.evaluations), 14)

    def test_check_cds_2(self):
        """Verify correct number of evaluations are produced when
        check_locus_tag = False."""
        self.eval_flags["check_locus_tag"] = False
        import_genome.check_cds(self.cds1, self.eval_flags)
        self.assertEqual(len(self.cds1.evaluations), 11)

    def test_check_cds_3(self):
        """Verify correct number of evaluations are produced when
        check_gene = False."""
        self.eval_flags["check_gene"] = False
        import_genome.check_cds(self.cds1, self.eval_flags)
        self.assertEqual(len(self.cds1.evaluations), 11)

    def test_check_cds_4(self):
        """Verify correct number of evaluations are produced when
        check_description = False."""
        self.eval_flags["check_description"] = False
        import_genome.check_cds(self.cds1, self.eval_flags)
        self.assertEqual(len(self.cds1.evaluations), 14)

    def test_check_cds_5(self):
        """Verify correct number of evaluations are produced when
        all False."""
        self.eval_flags["check_locus_tag"] = False
        self.eval_flags["check_gene"] = False
        self.eval_flags["check_description"] = False
        self.eval_flags["check_description_field"] = False
        import_genome.check_cds(self.cds1, self.eval_flags)
        self.assertEqual(len(self.cds1.evaluations), 8)

    def test_check_cds_6(self):
        """Verify correct number of errors with correct CDS feature."""
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error")
        self.assertEqual(count, 0)

    def test_check_cds_7(self):
        """Verify correct number of errors with incorrect amino acids."""
        self.cds1.translation = Seq("MB", IUPAC.protein)
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error", "warning")
        self.assertEqual(count, 2)

    def test_check_cds_8(self):
        """Verify correct number of errors with incorrect translation."""
        self.cds1.translation = Seq("MM", IUPAC.protein)
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error")
        self.assertEqual(count, 1)

    def test_check_cds_9(self):
        """Verify correct number of errors with missing translation."""
        self.cds1.translation = constants.EMPTY_PROTEIN_SEQ
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error")
        self.assertEqual(count, 2)

    def test_check_cds_10(self):
        """Verify correct number of errors with incorrect translation table."""
        self.cds1.translation_table = 1
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_cds_11(self):
        """Verify correct number of errors with incorrect start coordinate."""
        self.cds1.start = -1
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error")
        self.assertEqual(count, 1)

    def test_check_cds_12(self):
        """Verify correct number of errors with incorrect stop coordinate."""
        self.cds1.stop = -1
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error")
        self.assertEqual(count, 1)

    def test_check_cds_13(self):
        """Verify correct number of errors with incorrect parts."""
        self.cds1.parts = 0
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error")
        self.assertEqual(count, 1)

    def test_check_cds_14(self):
        """Verify correct number of errors with incorrect orientation."""
        self.cds1.orientation = "f"
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error")
        self.assertEqual(count, 1)

    def test_check_cds_15(self):
        """Verify correct number of errors with missing locus_tag."""
        self.cds1.locus_tag = ""
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error", "warning")
        self.assertEqual(count, 2)

    def test_check_cds_16(self):
        """Verify correct number of errors with incorrect locus_tag."""
        self.cds1.locus_tag = "ABCXYZ"
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_cds_17(self):
        """Verify correct number of errors with missing gene qualifier."""
        self.cds1.gene = ""
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error", "warning")
        self.assertEqual(count, 3)

    def test_check_cds_18(self):
        """Verify correct number of errors with non-integer gene qualifier."""
        self.cds1.gene = "A"
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error", "warning")
        self.assertEqual(count, 2)

    def test_check_cds_19(self):
        """Verify correct number of errors with non-integer gene qualifier."""
        self.cds1.gene = "11"
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_cds_20(self):
        """Verify correct number of errors with non-matching integer in
        gene qualifier and locus_tag."""
        self.cds1.gene = "11"
        import_genome.check_cds(self.cds1, self.eval_flags)
        count = count_status(self.cds1, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_cds_21(self):
        """Verify correct number of errors with non-matching integer in
        gene qualifier and locus_tag."""
        import_genome.check_cds(self.cds1, self.eval_flags,
                                description_field = "function")
        count = count_status(self.cds1, "error", "warning")
        self.assertEqual(count, 1)



class TestImportGenome8(unittest.TestCase):


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
        count = count_status(self.src1, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_source_8(self):
        """Verify correct number of errors with incorrect organism host genus."""
        self.src1._organism_host_genus = "Gordonia"
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        count = count_status(self.src1, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_source_9(self):
        """Verify correct number of errors with incorrect host host genus."""
        self.src1._host_host_genus = "Gordonia"
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        count = count_status(self.src1, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_source_10(self):
        """Verify correct number of errors with incorrect lab host host genus."""
        self.src1._lab_host_host_genus = "Gordonia"
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        count = count_status(self.src1, "error", "warning")
        self.assertEqual(count, 1)

    def test_check_source_11(self):
        """Verify correct number of errors with permissible synonym
        in organism host genus."""
        self.src1._organism_host_genus = "Mycobacterio"
        import_genome.check_source(self.src1, self.eval_flags,
                                   host_genus="Mycobacterium")
        count = count_status(self.src1, "error")
        self.assertEqual(count, 0)




class TestImportGenome9(unittest.TestCase):

    def setUp(self):
        self.tkt = ticket.ImportTicket()
        self.tkt.type = "replace"
        self.tkt.phage_id = "Trixie_tkt"
        self.gnm = genome.Genome()
        self.gnm.id = "Trixie_gnm"
        self.gnm.filename = "Trixie_file"
        self.bndl = bundle.Bundle()
        self.bndl._errors = 0
        self.success_path = pathlib.Path("success_folder")
        self.fail_path = pathlib.Path("fail_folder")
        self.paths_dict = {"success": self.success_path,
                           "fail": self.fail_path}
        self.flatfile_path = pathlib.Path("/folder/to/Trixie_flatfile.txt")




    def test_get_logfile_path_1(self):
        """Verify returned path when:
        bundle has 0 errors, paths_dict is provided, genome is provided,
        filepath is provided, file_ref is provided,
        and ticket is not provided."""
        self.bndl.genome_dict["flat_file"] = self.gnm
        logfile_path = import_genome.get_logfile_path(self.bndl,
                            paths_dict=self.paths_dict,
                            filepath=self.flatfile_path,
                            file_ref="flat_file")
        exp_name = self.gnm.id + "__" + self.flatfile_path.stem + ".log"
        exp_path = pathlib.Path(self.success_path, exp_name)
        self.assertEqual(logfile_path, exp_path)

    def test_get_logfile_path_2(self):
        """Verify returned path when:
        bundle has 1 errors, paths_dict is provided, genome is provided,
        filepath is provided, file_ref is provided,
        and ticket is not provided."""
        self.bndl._errors = 1
        self.bndl.genome_dict["flat_file"] = self.gnm
        logfile_path = import_genome.get_logfile_path(self.bndl,
                            paths_dict=self.paths_dict,
                            filepath=self.flatfile_path,
                            file_ref="flat_file")
        exp_name = self.gnm.id + "__" + self.flatfile_path.stem + ".log"
        exp_path = pathlib.Path(self.fail_path, exp_name)
        self.assertEqual(logfile_path, exp_path)

    def test_get_logfile_path_3(self):
        """Verify returned path when:
        bundle has 1 errors, paths_dict is provided, ticket is provided,
        filepath is not provided, file_ref is not provided,
        and genome is not provided."""
        self.bndl._errors = 1
        self.bndl.ticket = self.tkt
        logfile_path = import_genome.get_logfile_path(self.bndl,
                            paths_dict=self.paths_dict, filepath=None,
                            file_ref=None)
        exp_name = (self.tkt.phage_id + "__" + "no_file" + ".log")
        exp_path = pathlib.Path(self.fail_path, exp_name)
        self.assertEqual(logfile_path, exp_path)

    def test_get_logfile_path_4(self):
        """Verify returned path when:
        bundle has 1 errors, paths_dict is provided, genome is not provided,
        filepath is provided, file_ref is provided,
        and ticket is not provided."""
        self.bndl._errors = 1
        logfile_path = import_genome.get_logfile_path(self.bndl,
                            paths_dict=self.paths_dict,
                            filepath=self.flatfile_path,
                            file_ref="flat_file")
        exp_name = "no_id" + "__" + self.flatfile_path.stem + ".log"
        exp_path = pathlib.Path(self.fail_path, exp_name)
        self.assertEqual(logfile_path, exp_path)

    def test_get_logfile_path_5(self):
        """Verify returned path when:
        bundle has 1 errors, paths_dict is provided, genome is provided,
        filepath is provided, file_ref is provided, and ticket is provided."""
        self.bndl._errors = 1
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["flat_file"] = self.gnm
        logfile_path = import_genome.get_logfile_path(self.bndl,
                            paths_dict=self.paths_dict,
                            filepath=self.flatfile_path,
                            file_ref="flat_file")
        exp_name = self.gnm.id + "__" + self.flatfile_path.stem + ".log"
        exp_path = pathlib.Path(self.fail_path, exp_name)
        self.assertEqual(logfile_path, exp_path)

    def test_get_logfile_path_6(self):
        """Verify returned path when:
        bundle has 1 errors, genome is provided, paths_dict is provided,
        filepath is not provided, file_ref is provided,
        and ticket is provided."""
        self.bndl._errors = 1
        self.bndl.ticket = self.tkt
        self.bndl.genome_dict["flat_file"] = self.gnm
        logfile_path = import_genome.get_logfile_path(self.bndl,
                            paths_dict=self.paths_dict, filepath=None,
                            file_ref="flat_file")
        exp_name = (self.tkt.phage_id + "__" + "no_file" + ".log")
        exp_path = pathlib.Path(self.fail_path, exp_name)
        self.assertEqual(logfile_path, exp_path)

    def test_get_logfile_path_7(self):
        """Verify returned path when:
        bundle has 1 errors, paths_dict is not provided, genome is not provided,
        filepath is not provided, file_ref is not provided,
        and ticket is not provided."""
        self.bndl._errors = 1
        logfile_path = import_genome.get_logfile_path(self.bndl,
                            paths_dict=None, filepath=None,
                            file_ref=None)
        self.assertIsNone(logfile_path)




if __name__ == '__main__':
    unittest.main()
