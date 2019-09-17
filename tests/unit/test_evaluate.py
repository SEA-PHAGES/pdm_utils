""" Unit tests for evaluate functions."""


from pdm_utils.classes import genome
from pdm_utils.classes import source
from pdm_utils.classes import cds
from pdm_utils.classes import genomepair
from pdm_utils.constants import constants
from pdm_utils.pipelines.db_import import evaluate
from pdm_utils.classes import ticket
import unittest
from Bio.Seq import Seq


class TestEvaluateClass(unittest.TestCase):


    def setUp(self):

        self.null_set = constants.EMPTY_SET
        self.type_set = constants.IMPORT_TICKET_TYPE_SET
        self.run_mode_set = constants.RUN_MODES.keys()
        self.description_field_set = constants.DESCRIPTION_FIELD_SET

        self.add_ticket1 = ticket.GenomeTicket()
        self.add_ticket1.id = 1
        self.add_ticket1.type = "add"
        self.add_ticket1.phage_id = "Trixie_Draft"
        self.add_ticket1.run_mode = "phagesdb"
        self.add_ticket1.description_field = "product"
        self.add_ticket1.host_genus = "Mycobacterium smegmatis"
        self.add_ticket1.cluster = "A"
        self.add_ticket1.subcluster = "A2"
        self.add_ticket1.annotation_status = "final"
        self.add_ticket1.annotation_author = "hatfull"
        self.add_ticket1.annotation_qc = 1
        self.add_ticket1.retrieve_record = 1
        self.add_ticket1.accession = "ABC123.1"

        self.id_dupe_set = set([1])



    def test_check_ticket_structure_1(self):
        """Verify no error is produced with a correctly structured
        'add' ticket."""
        evaluate.check_ticket_structure(
            self.add_ticket1, type_set=self.type_set,
            description_field_set=self.description_field_set,
            null_set=self.null_set, run_mode_set=self.run_mode_set)
        errors = 0
        for evl in self.add_ticket1.evaluations:
            if evl.status == "error":
                errors += 1
        with self.subTest():
            self.assertEqual(len(self.add_ticket1.evaluations), 12)
        with self.subTest():
            self.assertEqual(errors, 0)


    def test_check_ticket_structure_2(self):
        """Verify an error is produced with an incorrectly structured
        'invalid' ticket 'type' field and duplicate id."""

        tkt = self.add_ticket1
        tkt.type = "invalid"
        evaluate.check_ticket_structure(
            self.add_ticket1, type_set=self.type_set,
            description_field_set=self.description_field_set,
            null_set=self.null_set, run_mode_set=self.run_mode_set,
            id_dupe_set=self.id_dupe_set)
        errors = 0
        for evl in tkt.evaluations:
            if evl.status == "error":
                errors += 1
        with self.subTest():
            self.assertEqual(len(tkt.evaluations), 12)
        with self.subTest():
            self.assertEqual(errors, 2)













    ###Below = pasted from test.Ticket.py, since some Ticket methods
    ###were moved to evaluate.py. These will now probably need to be implemented
    # in the 'compare_add_replace_ticket()'' evaluate functions.
    #
    # def test_check_update_ticket_1(self):
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[0].status, "correct")
    #
    # def test_check_update_ticket_2(self):
    #     """Primary Phage ID not in set."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[0].status, "error")
    #
    # def test_check_update_ticket_3(self):
    #     """host_genus is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.host_genus = "none"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[1].status, "error")
    #
    # def test_check_update_ticket_4(self):
    #     """Cluster is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.cluster = "none"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[2].status, "error")
    #
    # def test_check_update_ticket_5(self):
    #     """Status is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.status = "none"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[3].status, "error")
    #
    # def test_check_update_ticket_6(self):
    #     """Description Field is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.description_field = "Product"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[4].status, "error")
    #
    # def test_check_update_ticket_7(self):
    #     """Secondary Phage ID is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.secondary_phage_id = "Trixie"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[5].status, "error")
    #
    # def test_check_update_ticket_8(self):
    #     """Annotation Author is 0."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.annotation_author = "0"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[6].status, "correct")
    #
    # def test_check_update_ticket_9(self):
    #     """Annotation Author is not 1 or 0."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.annotation_author = "none"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[6].status, "error")
    #
    # def test_check_update_ticket_10(self):
    #     """Run Mode is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.run_mode = "phagesdb"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[7].status, "error")
    #






    #
    # def test_check_add_ticket_1(self):
    #     """Standard add ticket."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 0)
    #
    # def test_check_add_ticket_2(self):
    #     """Primary Phage ID already present."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_3(self):
    #     """Primary Phage ID is none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.phage_id = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_4(self):
    #     """host_genus is none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.host_genus = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_5(self):
    #     """Cluster is none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.cluster = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_6(self):
    #     """Status is none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.status = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_7(self):
    #     """Status is Final."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.status = "final"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_8(self):
    #     """Description Field is none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.description_field = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_9(self):
    #     """Secondary Phage ID is not none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.secondary_phage_id = "Trixie"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_10(self):
    #     """Annotation Author is 0."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.annotation_author = "0"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 0)
    #
    # def test_check_add_ticket_11(self):
    #     """Annotation Author is not 1 or 0."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.annotation_author = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_12(self):
    #     """Run Mode is none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.run_mode = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    #
    #
    #
    #
    #
    #
    #
    # def test_check_remove_ticket_1(self):
    #     """Standard remove ticket."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 0)
    #
    # def test_check_remove_ticket_2(self):
    #     """Primary Phage ID is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.phage_id = "Trixie"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_3(self):
    #     """host_genus is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.host_genus = "Mycobacterium"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_4(self):
    #     """Cluster is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.cluster = "A"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_5(self):
    #     """Subcluster is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.subcluster = "A2"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_6(self):
    #     """Status is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.status = "final"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_7(self):
    #     """Description Field is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.description_field = "Product"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_8(self):
    #     """Accession is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.accession = "ABC123"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_9(self):
    #     """Annotation Author is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.annotation_author = "1"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_10(self):
    #     """Run Mode is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.run_mode = "phagesdb"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_11(self):
    #     """Secondary Phage ID is not present."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.secondary_phage_id = "D29"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    #
    #
    #
    #
    #
    #
    #
    #
    # def test_check_replace_ticket_1(self):
    #     """Standard replace ticket."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 0)
    #
    # def test_check_replace_ticket_2(self):
    #     """Primary Phage ID is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock","none"])
    #     self.replace_ticket.phage_id = "none"
    #     self.replace_ticket.secondary_phage_id = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_3(self):
    #     """Primary Phage ID not present and different from Secondary Phage ID."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.phage_id = "D29"
    #     self.replace_ticket.secondary_phage_id = "L5"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_4(self):
    #     """Primary Phage ID is present and different from Secondary Phage ID."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.secondary_phage_id = "L5"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 2)
    #
    # def test_check_replace_ticket_5(self):
    #     """host_genus is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.host_genus = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_6(self):
    #     """Cluster is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.cluster = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_7(self):
    #     """Status is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.status = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_8(self):
    #     """Description Field is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.description_field = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_9(self):
    #     """Secondary Phage ID is not present."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.phage_id = "D29"
    #     self.replace_ticket.secondary_phage_id = "D29"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_10(self):
    #     """Annotation Author is 0."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.annotation_author = "0"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 0)
    #
    # def test_check_replace_ticket_11(self):
    #     """Annotation Author is not 1 or 0."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.annotation_author = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_12(self):
    #     """Run Mode is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.run_mode = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    # def test_check_ticket_1(self):
    #     """Check standard update ticket."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.update_ticket.evaluations), 0)
    #
    # def test_check_ticket_2(self):
    #     """Check update ticket with Primary Phage ID not present."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.update_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.update_ticket.evaluations), 1)
    #
    # def test_check_ticket_3(self):
    #     """Check standard add ticket."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 0)
    #
    # def test_check_ticket_4(self):
    #     """Check add ticket with Primary Phage ID already present."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.add_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_ticket_5(self):
    #     """Check standard remove ticket."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 0)
    #
    # def test_check_ticket_6(self):
    #     """Check remove ticket with Primary Phage ID not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.phage_id = "Trixie"
    #     self.remove_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_ticket_7(self):
    #     """Check standard replace ticket."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 0)
    #
    # def test_check_ticket_8(self):
    #     """Check replace ticket when Primary Phage ID is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock","none"])
    #     self.replace_ticket.phage_id = "none"
    #     self.replace_ticket.secondary_phage_id = "none"
    #     self.replace_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_ticket_9(self):
    #     """Check non-standard type of ticket."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.update_ticket.type = "other"
    #     self.update_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.update_ticket.evaluations), 0)
    #
    ###Above = pasted from test.Ticket.py, since some Ticket methods
    ###were moved to evaluate.py











class TestEvaluateClass2(unittest.TestCase):


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


    def test_check_phagesdb_genome_1(self):
        """Verify no error is produced with a correctly structured
        PhagesDB genome."""

        evaluate.check_phagesdb_genome(self.gnm, self.null_set)
        errors = 0
        for evl in self.gnm.evaluations:
            if evl.status == "error":
                errors += 1
        with self.subTest():
            self.assertEqual(len(self.gnm.evaluations), 8)
        with self.subTest():
            self.assertEqual(errors, 0)

    def test_check_phagesdb_genome_2(self):
        """Verify an error is produced with a PhagesDB genome with
        no id."""

        self.gnm.id = ""
        evaluate.check_phagesdb_genome(self.gnm, self.null_set)
        errors = 0
        for evl in self.gnm.evaluations:
            if evl.status == "error":
                errors += 1
        self.assertEqual(errors, 1)


    def test_check_phagesdb_genome_3(self):
        """Verify an error is produced with a PhagesDB genome with
        no name."""

        self.gnm.name = ""
        evaluate.check_phagesdb_genome(self.gnm, self.null_set)
        errors = 0
        for evl in self.gnm.evaluations:
            if evl.status == "error":
                errors += 1
        self.assertEqual(errors, 1)


    def test_check_phagesdb_genome_4(self):
        """Verify an error is produced with a PhagesDB genome with
        no host_genus."""

        self.gnm.host_genus = ""
        evaluate.check_phagesdb_genome(self.gnm, self.null_set)
        errors = 0
        for evl in self.gnm.evaluations:
            if evl.status == "error":
                errors += 1
        self.assertEqual(errors, 1)


    def test_check_phagesdb_genome_5(self):
        """Verify an error is produced with a PhagesDB genome with
        no cluster."""

        self.gnm.cluster = ""
        evaluate.check_phagesdb_genome(self.gnm, self.null_set)
        errors = 0
        for evl in self.gnm.evaluations:
            if evl.status == "error":
                errors += 1
        self.assertEqual(errors, 1)


    def test_check_phagesdb_genome_6(self):
        """Verify an error is produced with a PhagesDB genome with
        no subcluster."""

        self.gnm.subcluster = ""
        evaluate.check_phagesdb_genome(self.gnm, self.null_set)
        errors = 0
        for evl in self.gnm.evaluations:
            if evl.status == "error":
                errors += 1
        self.assertEqual(errors, 1)


    def test_check_phagesdb_genome_7(self):
        """Verify an error is produced with a PhagesDB genome with
        no accession."""

        self.gnm.accession = ""
        evaluate.check_phagesdb_genome(self.gnm, self.null_set)
        errors = 0
        for evl in self.gnm.evaluations:
            if evl.status == "error":
                errors += 1
        self.assertEqual(errors, 1)


    def test_check_phagesdb_genome_8(self):
        """Verify an error is produced with a PhagesDB genome with
        no filename."""

        self.gnm.filename = ""
        evaluate.check_phagesdb_genome(self.gnm, self.null_set)
        errors = 0
        for evl in self.gnm.evaluations:
            if evl.status == "error":
                errors += 1
        self.assertEqual(errors, 1)


    def test_check_phagesdb_genome_9(self):
        """Verify an error is produced with a PhagesDB genome with
        no sequence."""

        self.gnm.seq = ""
        evaluate.check_phagesdb_genome(self.gnm, self.null_set)
        errors = 0
        for evl in self.gnm.evaluations:
            if evl.status == "error":
                errors += 1
        self.assertEqual(errors, 1)





    def test_check_source_for_import_1(self):
        """Verify correct number of evaluations are produced when
        check_id_typo = True and check_host_typo = True."""
        src_ftr = source.Source()
        evaluate.check_source_for_import(src_ftr)
        self.assertEqual(len(src_ftr.evaluations), 4)

    def test_check_source_for_import_2(self):
        """Verify correct number of evaluations are produced when
        check_id_typo = False and check_host_typo = True."""
        src_ftr = source.Source()
        evaluate.check_source_for_import(src_ftr, check_id_typo=False)
        self.assertEqual(len(src_ftr.evaluations), 3)

    def test_check_source_for_import_3(self):
        """Verify correct number of evaluations are produced when
        check_id_typo = True and check_host_typo = False."""
        src_ftr = source.Source()
        evaluate.check_source_for_import(src_ftr, check_host_typo=False)
        self.assertEqual(len(src_ftr.evaluations), 1)









    def test_check_cds_for_import_1(self):
        """Verify correct number of evaluations are produced when
        check_locus_tag = True, check_gene = True, and
        check_description = True."""
        cds_ftr = cds.Cds()
        evaluate.check_cds_for_import(cds_ftr)
        self.assertEqual(len(cds_ftr.evaluations), 13)

    def test_check_cds_for_import_2(self):
        """Verify correct number of evaluations are produced when
        check_locus_tag = False, check_gene = True, and
        check_description = True."""
        cds_ftr = cds.Cds()
        evaluate.check_cds_for_import(cds_ftr, check_locus_tag=False)
        self.assertEqual(len(cds_ftr.evaluations), 10)

    def test_check_cds_for_import_3(self):
        """Verify correct number of evaluations are produced when
        check_locus_tag = True, check_gene = False, and
        check_description = True."""
        cds_ftr = cds.Cds()
        evaluate.check_cds_for_import(cds_ftr, check_gene=False)
        self.assertEqual(len(cds_ftr.evaluations), 10)

    def test_check_cds_for_import_4(self):
        """Verify correct number of evaluations are produced when
        check_locus_tag = True, check_gene = True, and
        check_description = False."""
        cds_ftr = cds.Cds()
        evaluate.check_cds_for_import(cds_ftr, check_description=False)
        self.assertEqual(len(cds_ftr.evaluations), 12)


    # TODO test_check_cds_for_import_5 to test check_description_field parameter.




    def test_compare_genomes_1(self):
        """Verify correct number of evaluations are produced when."""
        genome1 = genome.Genome()
        genome2 = genome.Genome()
        genome_pair = genomepair.GenomePair()
        genome_pair.genome1 = genome1
        genome_pair.genome2 = genome2
        evaluate.compare_genomes(genome_pair)
        self.assertEqual(len(genome_pair.evaluations), 8)





    def test_check_genome_for_import_1(self):
        """Verify correct number of evaluations are produced using
        'add' ticket and all eval_flags 'True'."""
        evaluate.check_genome_for_import(self.gnm, self.tkt, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 28)

    def test_check_genome_for_import_2(self):
        """Verify correct number of evaluations are produced using
        'replace' ticket."""
        self.tkt.type = "replace"
        evaluate.check_genome_for_import(self.gnm, self.tkt, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 27)


    def test_check_genome_for_import_3(self):
        """Verify correct number of evaluations are produced using
        'check_seq' as False."""
        self.tkt.eval_flags["check_seq"] = False
        evaluate.check_genome_for_import(self.gnm, self.tkt, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 27)


    def test_check_genome_for_import_4(self):
        """Verify correct number of evaluations are produced using
        'check_id_typo' as False."""
        self.tkt.eval_flags["check_id_typo"] = False
        evaluate.check_genome_for_import(self.gnm, self.tkt, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 25)


    def test_check_genome_for_import_5(self):
        """Verify correct number of evaluations are produced using
        'check_host_typo' as False."""
        self.tkt.eval_flags["check_host_typo"] = False
        evaluate.check_genome_for_import(self.gnm, self.tkt, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 24)


    def test_check_genome_for_import_6(self):
        """Verify correct number of evaluations are produced using
        'check_author' as False."""
        self.tkt.eval_flags["check_author"] = False
        evaluate.check_genome_for_import(self.gnm, self.tkt, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 26)



        # self.tkt.eval_flags["check_trna"] = True





        # self.tkt.eval_flags["check_seq"] = True
        # self.tkt.eval_flags["check_id_typo"] = True
        # self.tkt.eval_flags["check_host_typo"] = True
        # self.tkt.eval_flags["check_author"] = True
        # self.tkt.eval_flags["check_trna"] = True



if __name__ == '__main__':
    unittest.main()
