""" Unit tests for import functions."""


from pdm_utils.classes import bundle
from pdm_utils.classes import genome
from pdm_utils.classes import source
from pdm_utils.classes import cds
from pdm_utils.classes import genomepair
from pdm_utils.constants import constants
from pdm_utils.pipelines.db_import import import_main
from pdm_utils.classes import ticket
import unittest
from Bio.Seq import Seq


class TestImportClass1(unittest.TestCase):


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
        self.add_ticket1.annotation_status = "draft"
        self.add_ticket1.annotation_author = "hatfull"
        self.add_ticket1.annotation_qc = 1
        self.add_ticket1.retrieve_record = 1
        self.add_ticket1.accession = "ABC123.1"

        self.id_dupe_set = set([1])



    def test_check_ticket_structure_1(self):
        """Verify no error is produced with a correctly structured
        'add' ticket."""
        import_main.check_ticket_structure(
            self.add_ticket1, type_set=self.type_set,
            description_field_set=self.description_field_set,
            null_set=self.null_set, run_mode_set=self.run_mode_set)
        errors = 0
        for evl in self.add_ticket1.evaluations:
            if evl.status == "error":
                errors += 1
        with self.subTest():
            self.assertEqual(len(self.add_ticket1.evaluations), 13)
        with self.subTest():
            self.assertEqual(errors, 0)


    def test_check_ticket_structure_2(self):
        """Verify an error is produced with an incorrectly structured
        'invalid' ticket 'type' field and duplicate id."""

        tkt = self.add_ticket1
        tkt.type = "invalid"
        import_main.check_ticket_structure(
            self.add_ticket1, type_set=self.type_set,
            description_field_set=self.description_field_set,
            null_set=self.null_set, run_mode_set=self.run_mode_set,
            id_dupe_set=self.id_dupe_set)
        errors = 0
        for evl in tkt.evaluations:
            if evl.status == "error":
                errors += 1
        with self.subTest():
            self.assertEqual(len(tkt.evaluations), 13)
        with self.subTest():
            self.assertEqual(errors, 2)





class TestImportClass2(unittest.TestCase):


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




    def test_check_source_1(self):
        """Verify correct number of evaluations are produced when
        check_id_typo = True and check_host_typo = True."""
        src_ftr = source.Source()
        import_main.check_source(src_ftr)
        self.assertEqual(len(src_ftr.evaluations), 4)

    def test_check_source_2(self):
        """Verify correct number of evaluations are produced when
        check_id_typo = False and check_host_typo = True."""
        src_ftr = source.Source()
        import_main.check_source(src_ftr, check_id_typo=False)
        self.assertEqual(len(src_ftr.evaluations), 3)

    def test_check_source_3(self):
        """Verify correct number of evaluations are produced when
        check_id_typo = True and check_host_typo = False."""
        src_ftr = source.Source()
        import_main.check_source(src_ftr, check_host_typo=False)
        self.assertEqual(len(src_ftr.evaluations), 1)









    def test_check_cds_1(self):
        """Verify correct number of evaluations are produced when
        check_locus_tag = True, check_gene = True, and
        check_description = True."""
        cds_ftr = cds.Cds()
        import_main.check_cds(cds_ftr)
        self.assertEqual(len(cds_ftr.evaluations), 13)

    def test_check_cds_2(self):
        """Verify correct number of evaluations are produced when
        check_locus_tag = False, check_gene = True, and
        check_description = True."""
        cds_ftr = cds.Cds()
        import_main.check_cds(cds_ftr, check_locus_tag=False)
        self.assertEqual(len(cds_ftr.evaluations), 10)

    def test_check_cds_3(self):
        """Verify correct number of evaluations are produced when
        check_locus_tag = True, check_gene = False, and
        check_description = True."""
        cds_ftr = cds.Cds()
        import_main.check_cds(cds_ftr, check_gene=False)
        self.assertEqual(len(cds_ftr.evaluations), 10)

    def test_check_cds_4(self):
        """Verify correct number of evaluations are produced when
        check_locus_tag = True, check_gene = True, and
        check_description = False."""
        cds_ftr = cds.Cds()
        import_main.check_cds(cds_ftr, check_description=False)
        self.assertEqual(len(cds_ftr.evaluations), 12)


    # TODO test_check_cds_5 to test check_description_field parameter.




    def test_compare_genomes_1(self):
        """Verify correct number of evaluations are produced when."""
        genome1 = genome.Genome()
        genome2 = genome.Genome()
        genome_pair = genomepair.GenomePair()
        genome_pair.genome1 = genome1
        genome_pair.genome2 = genome2
        import_main.compare_genomes(genome_pair)
        self.assertEqual(len(genome_pair.evaluations), 8)





    def test_check_genome_1(self):
        """Verify correct number of evaluations are produced using
        'add' ticket and all eval_flags 'True'."""
        import_main.check_genome(self.gnm, self.tkt, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 29)

    def test_check_genome_2(self):
        """Verify correct number of evaluations are produced using
        'replace' ticket."""
        self.tkt.type = "replace"
        import_main.check_genome(self.gnm, self.tkt, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 28)


    def test_check_genome_3(self):
        """Verify correct number of evaluations are produced using
        'check_seq' as False."""
        self.tkt.eval_flags["check_seq"] = False
        import_main.check_genome(self.gnm, self.tkt, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 28)


    def test_check_genome_4(self):
        """Verify correct number of evaluations are produced using
        'check_id_typo' as False."""
        self.tkt.eval_flags["check_id_typo"] = False
        import_main.check_genome(self.gnm, self.tkt, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 26)


    def test_check_genome_5(self):
        """Verify correct number of evaluations are produced using
        'check_host_typo' as False."""
        self.tkt.eval_flags["check_host_typo"] = False
        import_main.check_genome(self.gnm, self.tkt, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 25)


    def test_check_genome_6(self):
        """Verify correct number of evaluations are produced using
        'check_author' as False."""
        self.tkt.eval_flags["check_author"] = False
        import_main.check_genome(self.gnm, self.tkt, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 27)






class TestImportClass3(unittest.TestCase):

    def setUp(self):
        self.bndl = bundle.Bundle()
        self.tkt = ticket.GenomeTicket()
        self.gnm1 = genome.Genome()
        self.gnm2 = genome.Genome()




    def test_check_bundle_1(self):
        """Verify all check methods are called."""
        self.tkt.type = "replace"
        self.tkt.cluster = "retrieve"
        self.tkt.subcluster = "retain"
        self.bndl.ticket = self.tkt
        import_main.check_bundle(self.bndl)
        self.assertEqual(len(self.bndl.evaluations), 9)

    def test_check_bundle_2(self):
        """Verify some check methods are called when there is no Ticket."""
        import_main.check_bundle(self.bndl)
        self.assertEqual(len(self.bndl.evaluations), 1)

    def test_check_bundle_3(self):
        """Verify some check methods are called when there are no 'retrieve'
        Ticket attributes."""
        self.tkt.type = "replace"
        self.tkt.subcluster = "retain"
        self.bndl.ticket = self.tkt
        import_main.check_bundle(self.bndl)
        self.assertEqual(len(self.bndl.evaluations), 7)

    def test_check_bundle_4(self):
        """Verify some check methods are called when there are no 'retain'
        Ticket attributes."""
        self.tkt.type = "replace"
        self.tkt.cluster = "retrieve"
        self.bndl.ticket = self.tkt
        import_main.check_bundle(self.bndl)
        self.assertEqual(len(self.bndl.evaluations), 8)

    def test_check_bundle_5(self):
        """Verify some check methods are called when there is no 'replace'
        Ticket."""
        self.tkt.cluster = "retrieve"
        self.bndl.ticket = self.tkt
        import_main.check_bundle(self.bndl)
        self.assertEqual(len(self.bndl.evaluations), 6)









if __name__ == '__main__':
    unittest.main()
