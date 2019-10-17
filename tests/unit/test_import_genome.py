""" Unit tests for import functions."""


from pdm_utils.classes import bundle
from pdm_utils.classes import genome
from pdm_utils.classes import source
from pdm_utils.classes import cds
from pdm_utils.classes import genomepair
from pdm_utils.constants import constants
from pdm_utils.functions import run_modes
from pdm_utils.pipelines.db_import import import_genome
from pdm_utils.classes import ticket
import unittest
from Bio.Seq import Seq
from unittest.mock import patch
from pdm_utils.classes import mysqlconnectionhandler as mch


class TestImportGenomeClass1(unittest.TestCase):


    def setUp(self):

        self.null_set = constants.EMPTY_SET
        self.type_set = constants.IMPORT_TICKET_TYPE_SET
        self.run_mode_set = run_modes.RUN_MODES.keys()
        self.description_field_set = constants.DESCRIPTION_FIELD_SET

        self.data_dict = {
            "host_genus": "Mycobacterium smegmatis",
            "cluster": "A",
            "subcluster": "A2",
            "annotation_status": "draft",
            "annotation_author": 1,
            "retrieve_record": 1,
            "accession": "ABC123.1"
            }

        self.add_ticket1 = ticket.GenomeTicket()
        self.add_ticket1.id = 1
        self.add_ticket1.type = "add"
        self.add_ticket1.phage_id = "Trixie_Draft"
        self.add_ticket1.run_mode = "phagesdb"
        self.add_ticket1.description_field = "product"
        self.add_ticket1.data_dict = self.data_dict
        self.add_ticket1.eval_flags = {"a":1, "b":2}

        self.id_dupe_set = set([1])


    def test_check_ticket_1(self):
        """Verify no error is produced with a correctly structured
        'add' ticket."""
        import_genome.check_ticket(
            self.add_ticket1, type_set=self.type_set,
            description_field_set=self.description_field_set,
            null_set=self.null_set, run_mode_set=self.run_mode_set)
        errors = 0
        for evl in self.add_ticket1.evaluations:
            if evl.status == "error":
                errors += 1
        with self.subTest():
            self.assertEqual(len(self.add_ticket1.evaluations), 9)
        with self.subTest():
            self.assertEqual(errors, 0)


    def test_check_ticket_2(self):
        """Verify an error is produced with an incorrectly structured
        'invalid' ticket 'type' field and duplicate id."""

        tkt = self.add_ticket1
        tkt.type = "invalid"
        import_genome.check_ticket(
            self.add_ticket1, type_set=self.type_set,
            description_field_set=self.description_field_set,
            null_set=self.null_set, run_mode_set=self.run_mode_set,
            id_dupe_set=self.id_dupe_set)
        errors = 0
        for evl in tkt.evaluations:
            if evl.status == "error":
                errors += 1
        with self.subTest():
            self.assertEqual(len(tkt.evaluations), 9)
        with self.subTest():
            self.assertEqual(errors, 2)





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




    def test_check_source_1(self):
        """Verify correct number of evaluations are produced when
        none are False."""
        eval_flags = {"check_id_typo": True,
                      "check_host_typo": True}
        src_ftr = source.Source()
        import_genome.check_source(src_ftr, eval_flags)
        self.assertEqual(len(src_ftr.evaluations), 4)

    def test_check_source_2(self):
        """Verify correct number of evaluations are produced when
        check_id_typo = False."""
        eval_flags = {"check_id_typo": False,
                      "check_host_typo": True}
        src_ftr = source.Source()
        import_genome.check_source(src_ftr, eval_flags)
        self.assertEqual(len(src_ftr.evaluations), 3)

    def test_check_source_3(self):
        """Verify correct number of evaluations are produced when
        check_host_typo = False."""
        eval_flags = {"check_id_typo": True,
                      "check_host_typo": False}
        src_ftr = source.Source()
        import_genome.check_source(src_ftr, eval_flags)
        self.assertEqual(len(src_ftr.evaluations), 1)




    def test_check_cds_1(self):
        """Verify correct number of evaluations are produced when
        none are False."""
        eval_flags = {"check_locus_tag": True,
                      "check_gene": True,
                      "check_description": True,
                      "check_description_field": True}
        cds_ftr = cds.Cds()
        import_genome.check_cds(cds_ftr, eval_flags)
        self.assertEqual(len(cds_ftr.evaluations), 13)

    def test_check_cds_2(self):
        """Verify correct number of evaluations are produced when
        check_locus_tag = False."""
        eval_flags = {"check_locus_tag": False,
                      "check_gene": True,
                      "check_description": True,
                      "check_description_field": True}
        cds_ftr = cds.Cds()
        import_genome.check_cds(cds_ftr, eval_flags)
        self.assertEqual(len(cds_ftr.evaluations), 10)

    def test_check_cds_3(self):
        """Verify correct number of evaluations are produced when
        check_gene = False."""
        eval_flags = {"check_locus_tag": True,
                      "check_gene": False,
                      "check_description": True,
                      "check_description_field": True}
        cds_ftr = cds.Cds()
        import_genome.check_cds(cds_ftr, eval_flags)
        self.assertEqual(len(cds_ftr.evaluations), 10)

    def test_check_cds_4(self):
        """Verify correct number of evaluations are produced when
        check_description = False."""
        eval_flags = {"check_locus_tag": True,
                      "check_gene": True,
                      "check_description": False,
                      "check_description_field": True}
        cds_ftr = cds.Cds()
        import_genome.check_cds(cds_ftr, eval_flags)
        self.assertEqual(len(cds_ftr.evaluations), 12)

    def test_check_cds_5(self):
        """Verify correct number of evaluations are produced when
        all False."""
        eval_flags = {"check_locus_tag": False,
                      "check_gene": False,
                      "check_description": False,
                      "check_description_field": False}
        cds_ftr = cds.Cds()
        import_genome.check_cds(cds_ftr, eval_flags)
        self.assertEqual(len(cds_ftr.evaluations), 6)


    # TODO test_check_cds_5 to test check_description_field parameter
    # once the cds.check_description_field() has been implemented.




    def test_compare_genomes_1(self):
        """Verify correct number of evaluations are produced when
        'check_replace' is True."""
        eval_flags = {"check_replace": True}
        genome1 = genome.Genome()
        genome2 = genome.Genome()
        genome_pair = genomepair.GenomePair()
        genome_pair.genome1 = genome1
        genome_pair.genome2 = genome2
        import_genome.compare_genomes(genome_pair, eval_flags)
        self.assertEqual(len(genome_pair.evaluations), 8)

    def test_compare_genomes_2(self):
        """Verify correct number of evaluations are produced when
        'check_replace' is False."""
        eval_flags = {"check_replace": False}
        genome1 = genome.Genome()
        genome2 = genome.Genome()
        genome_pair = genomepair.GenomePair()
        genome_pair.genome1 = genome1
        genome_pair.genome2 = genome2
        import_genome.compare_genomes(genome_pair, eval_flags)
        self.assertEqual(len(genome_pair.evaluations), 7)




    def test_check_genome_1(self):
        """Verify correct number of evaluations are produced using
        'add' ticket and all eval_flags 'True'."""
        import_genome.check_genome(self.gnm, self.tkt.type,
            self.tkt.eval_flags, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 29)

    def test_check_genome_2(self):
        """Verify correct number of evaluations are produced using
        'replace' ticket."""
        self.tkt.type = "replace"
        import_genome.check_genome(self.gnm, self.tkt.type,
            self.tkt.eval_flags, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 28)

    def test_check_genome_3(self):
        """Verify correct number of evaluations are produced using
        'check_seq' as False."""
        self.tkt.eval_flags["check_seq"] = False
        import_genome.check_genome(self.gnm, self.tkt.type,
            self.tkt.eval_flags, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 28)

    def test_check_genome_4(self):
        """Verify correct number of evaluations are produced using
        'check_id_typo' as False."""
        self.tkt.eval_flags["check_id_typo"] = False
        import_genome.check_genome(self.gnm, self.tkt.type,
            self.tkt.eval_flags, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 26)

    def test_check_genome_5(self):
        """Verify correct number of evaluations are produced using
        'check_host_typo' as False."""
        self.tkt.eval_flags["check_host_typo"] = False
        import_genome.check_genome(self.gnm, self.tkt.type,
            self.tkt.eval_flags, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 25)

    def test_check_genome_6(self):
        """Verify correct number of evaluations are produced using
        'check_author' as False."""
        self.tkt.eval_flags["check_author"] = False
        import_genome.check_genome(self.gnm, self.tkt.type,
            self.tkt.eval_flags, self.null_set,
            self.id_set, self.seq_set, self.host_set,
            self.cluster_set, self.subcluster_set)
        self.assertEqual(len(self.gnm.evaluations), 27)




class TestImportGenomeClass3(unittest.TestCase):

    def setUp(self):
        self.bndl = bundle.Bundle()
        self.tkt = ticket.GenomeTicket()
        self.gnm1 = genome.Genome()
        self.gnm2 = genome.Genome()




    def test_check_bundle_1(self):
        """Verify all check methods are called."""
        self.tkt.type = "replace"
        self.tkt.data_retrieve = set(["cluster"])
        self.tkt.data_retain = set(["subcluster"])
        self.tkt.data_ticket = set(["host_genus"])
        self.bndl.ticket = self.tkt
        import_genome.check_bundle(self.bndl)
        self.assertEqual(len(self.bndl.evaluations), 6)

    def test_check_bundle_2(self):
        """Verify some check methods are called when there is no Ticket."""
        import_genome.check_bundle(self.bndl)
        self.assertEqual(len(self.bndl.evaluations), 1)

    def test_check_bundle_3(self):
        """Verify some check methods are called when there are no 'retrieve'
        Ticket attributes."""
        self.tkt.type = "replace"
        self.tkt.data_retain = set(["subcluster"])
        self.tkt.data_ticket = set(["host_genus"])
        self.bndl.ticket = self.tkt
        import_genome.check_bundle(self.bndl)
        self.assertEqual(len(self.bndl.evaluations), 5)


    def test_check_bundle_5(self):
        """Verify some check methods are called when there is no 'replace'
        Ticket."""
        self.tkt.data_retrieve = set(["cluster"])
        self.tkt.data_retain = set(["subcluster"])
        self.tkt.data_ticket = set(["host_genus"])
        self.bndl.ticket = self.tkt
        import_genome.check_bundle(self.bndl)
        self.assertEqual(len(self.bndl.evaluations), 4)







class TestImportGenomeClass4(unittest.TestCase):

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
        self.gnm2.type = "phamerator"
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
                null_set=self.null_set,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                gnm_key="flat_file",
                gnm_pair_key="flat_file_phamerator")
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
                null_set=self.null_set,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                gnm_key="flat_file",
                gnm_pair_key="flat_file_phamerator")
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
                null_set=self.null_set,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                gnm_key="flat_file",
                gnm_pair_key="flat_file_phamerator")
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
                null_set=self.null_set,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                gnm_key="flat_file",
                gnm_pair_key="flat_file_phamerator")
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
        'add' ticket, no genome, 'flat_file_phamerator' genome_pair."""
        self.tkt.type = "add"
        self.bndl.ticket = self.tkt
        self.bndl.genome_pair_dict["flat_file_phamerator"] = self.genome_pair
        import_genome.run_checks(
                self.bndl,
                null_set=self.null_set,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                gnm_key="flat_file",
                gnm_pair_key="flat_file_phamerator")
        with self.subTest():
            self.assertTrue(len(self.bndl.evaluations) > 0)
        with self.subTest():
            self.assertTrue(len(self.genome_pair.evaluations) == 0)
        with self.subTest():
            self.assertTrue(len(self.gnm1.evaluations) == 0)


    def test_run_checks_6(self):
        """Verify run_checks works using a bundle with:
        'replace' ticket, no genome, 'flat_file_phamerator' genome_pair."""
        self.tkt.type = "replace"
        self.bndl.ticket = self.tkt
        self.bndl.genome_pair_dict["flat_file_phamerator"] = self.genome_pair
        import_genome.run_checks(
                self.bndl,
                null_set=self.null_set,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                gnm_key="flat_file",
                gnm_pair_key="flat_file_phamerator")
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
        self.bndl.genome_pair_dict["flat_file_phamerator_x"] = self.genome_pair
        import_genome.run_checks(
                self.bndl,
                null_set=self.null_set,
                accession_set=self.accession_set,
                phage_id_set=self.phage_id_set,
                seq_set=self.seq_set, host_genus_set=self.host_genus_set,
                cluster_set=self.cluster_set,
                subcluster_set=self.subcluster_set,
                gnm_key="flat_file",
                gnm_pair_key="flat_file_phamerator")
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



    def test_import_into_db_1(self):
        """Verify import_into_db works using a bundle with:
        1 error, prod_run = True."""
        self.bndl._errors = 1
        result = import_genome.import_into_db(self.bndl, self.sql_handle,
                    gnm_key="", prod_run=True)
        self.assertFalse(result)


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










if __name__ == '__main__':
    unittest.main()
