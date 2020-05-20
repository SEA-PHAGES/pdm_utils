""" Unit tests for the Genome class."""

from datetime import datetime
import pathlib
import unittest

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from pdm_utils.classes import cds, trna, tmrna, source
from pdm_utils.classes import genome
from pdm_utils.constants import constants

class TestGenomeClass1(unittest.TestCase):

    def setUp(self):
        self.gnm = genome.Genome()

        self.cds1 = cds.Cds()
        self.cds1.description = ""
        self.cds1.product = ""
        self.cds1.function = ""
        self.cds1.note = ""

        self.cds2 = cds.Cds()
        self.cds2.description = ""
        self.cds2.product = ""
        self.cds2.function = ""
        self.cds2.note = ""


        self.cds3 = cds.Cds()
        self.cds4 = cds.Cds()

        self.trna1 = trna.TrnaFeature()
        self.trna2 = trna.TrnaFeature()
        self.trna3 = trna.TrnaFeature()
        self.trna4 = trna.TrnaFeature()

        self.tmrna1 = tmrna.TmrnaFeature()
        self.tmrna2 = tmrna.TmrnaFeature()
        self.tmrna3 = tmrna.TmrnaFeature()

        self.source1 = source.Source()

    def test_set_filename_1(self):
        """Confirm file path is split appropriately."""
        filepath = pathlib.Path("/path/to/folder/Trixie.gbk")
        self.gnm.set_filename(filepath)
        self.assertEqual(self.gnm.filename, "Trixie")




    def test_set_sequence_1(self):
        """Check that sequence is set appropriately."""
        self.gnm.set_sequence("aaggcga")
        with self.subTest():
            self.assertEqual(self.gnm.seq, "AAGGCGA")
        with self.subTest():
            self.assertIsInstance(self.gnm.seq, Seq)
        with self.subTest():
            self.assertEqual(self.gnm.length, 7)
        with self.subTest():
            self.assertEqual(self.gnm.gc, 57.1429)

    def test_set_sequence_2(self):
        """Check that sequence is set appropriately if it is an empty string."""
        seq = ""
        self.gnm.set_sequence(seq)
        with self.subTest():
            self.assertEqual(self.gnm.seq, "")
        with self.subTest():
            self.assertIsInstance(self.gnm.seq, Seq)
        with self.subTest():
            self.assertEqual(self.gnm.length, 0)
        with self.subTest():
            self.assertEqual(self.gnm.gc, -1)




    def test_set_accession_1(self):
        """Check that accession is set appropriately."""
        accession = "ABC123.1"
        self.gnm.set_accession(accession, "empty_string")
        self.assertEqual(self.gnm.accession, "ABC123")

    def test_set_accession_2(self):
        """Check that an empty string accession is set appropriately
        from None."""
        accession = None
        self.gnm.set_accession(accession, "empty_string")
        self.assertEqual(self.gnm.accession, "")

    def test_set_accession_3(self):
        """Check that 'none' string accession is set appropriately
        from None."""
        accession = None
        self.gnm.set_accession(accession, "none_string")
        self.assertEqual(self.gnm.accession, "none")

    def test_set_accession_4(self):
        """Check that a None accession is set appropriately from None."""
        accession = None
        self.gnm.set_accession(accession, "none_object")
        self.assertIsNone(self.gnm.accession)

    def test_set_accession_5(self):
        """Check that a None accession is set appropriately from an
        empty string."""
        accession = ""
        self.gnm.set_accession(accession, "none_object")
        self.assertIsNone(self.gnm.accession)




    def test_set_cds_features_1(self):
        """Check that CDS feature list is set and length is computed."""
        features_list = [0, 1, 2, 3]
        self.gnm.set_cds_features(features_list)
        with self.subTest():
            self.assertEqual(len(self.gnm.cds_features), 4)
        with self.subTest():
            self.assertEqual(self.gnm._cds_features_tally, 4)




    def test_set_cds_id_list_1(self):
        """Check that CDS feature identifier lists are set."""
        cds1 = cds.Cds()
        cds1._start_end_id = (1, 5)
        cds1._end_orient_id = (5, "forward")

        cds2 = cds.Cds()
        cds2._start_end_id = (21, 2)
        cds2._end_orient_id = (2, "reverse")

        features_list = [cds1, cds2]
        self.gnm.set_cds_features(features_list)
        self.gnm.set_cds_id_list()

        start_end_id_list = [(1,5), (21,2)]
        end_strand_id_list = [(5, "forward"), (2, "reverse")]

        with self.subTest():
            self.assertEqual(
                self.gnm._cds_start_end_ids, start_end_id_list)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_end_orient_ids, end_strand_id_list)




    def test_set_trna_features_1(self):
        """Check that tRNA feature list is set and length is computed."""
        features_list = [0,1,2,3]
        self.gnm.set_trna_features(features_list)
        with self.subTest():
            self.assertEqual(len(self.gnm.trna_features), 4)
        with self.subTest():
            self.assertEqual(self.gnm._trna_features_tally, 4)




    def test_set_tmrna_features_1(self):
        """Check that tmRNA feature list is set and length is computed."""
        features_list = [0,1,2,3]
        self.gnm.set_tmrna_features(features_list)
        with self.subTest():
            self.assertEqual(len(self.gnm.tmrna_features), 4)
        with self.subTest():
            self.assertEqual(self.gnm._tmrna_features_tally, 4)



    def test_set_source_features_1(self):
        """Check that source feature list is set and length is computed."""
        features_list = [0,1,2,3]
        self.gnm.set_source_features(features_list)
        with self.subTest():
            self.assertEqual(len(self.gnm.source_features), 4)
        with self.subTest():
            self.assertEqual(self.gnm._source_features_tally, 4)




    def test_set_id_1(self):
        """Check that name without '_Draft' suffix is not changed."""
        self.gnm.set_id(value="Trixie")
        self.assertEqual(self.gnm.id, "Trixie")

    def test_set_id_2(self):
        """Check that '_Draft' suffix is removed."""
        self.gnm.set_id(value="Trixie_Draft")
        self.assertEqual(self.gnm.id, "Trixie")

    def test_set_id_3(self):
        """Check that the id is set from the name attribute."""
        self.gnm.name = "Trixie_Draft"
        self.gnm.set_id(attribute="name")
        self.assertEqual(self.gnm.id, "Trixie")

    def test_set_id_4(self):
        """Check that the id is empty from an invalid attribute."""
        self.gnm.id = "not empty"
        self.gnm.set_id(attribute="invalid")
        self.assertEqual(self.gnm.id, "")

    def test_set_id_5(self):
        """Check that the id is empty from no value or attribute."""
        self.gnm.id = "not empty"
        self.gnm.set_id()
        self.assertEqual(self.gnm.id, "")




    def test_set_host_genus_1(self):
        """Check that host_genus name is split appropriately."""
        host = "Mycobacterium smegmatis"
        self.gnm.set_host_genus(value=host, format="none_string")
        self.assertEqual(self.gnm.host_genus, "Mycobacterium")

    def test_set_host_genus_2(self):
        """Check that whitespace is removed."""
        host = "  Mycobacterium smegmatis  "
        self.gnm.set_host_genus(value=host, format="none_string")
        self.assertEqual(self.gnm.host_genus, "Mycobacterium")

    def test_set_host_genus_3(self):
        """Check that none is set appropriately."""
        host = ""
        self.gnm.set_host_genus(value=host, format="none_string")
        self.assertEqual(self.gnm.host_genus, "none")

    def test_set_host_genus_4(self):
        """Check that None object is set appropriately."""
        host = ""
        self.gnm.set_host_genus(value=host, format="none_object")
        self.assertIsNone(self.gnm.host_genus)

    def test_set_host_genus_5(self):
        """Check that the host_genus is set from the name attribute."""
        self.gnm.name = "Mycobacterium smegmatis"
        self.gnm.set_host_genus(attribute="name")
        self.assertEqual(self.gnm.host_genus, "Mycobacterium")

    def test_set_host_genus_6(self):
        """Check that the host_genus is set from the
        description_host_genus attribute."""
        self.gnm._description_host_genus = "Mycobacterium smegmatis"
        self.gnm.set_host_genus(attribute="_description_host_genus")
        self.assertEqual(self.gnm.host_genus, "Mycobacterium")

    def test_set_host_genus_7(self):
        """Check that the host_genus is empty from an invalid attribute."""
        self.gnm.set_host_genus(attribute="invalid")
        self.assertEqual(self.gnm.host_genus, "")

    def test_set_host_genus_8(self):
        """Check that the host_genus is empty from an empty value and
        empty attribute."""
        self.gnm.host_genus = "not empty"
        self.gnm.set_host_genus()
        self.assertEqual(self.gnm.host_genus, "")




    def test_set_cluster_1(self):
        """Check that standard Cluster is set appropriately."""
        cluster = "A"
        self.gnm.set_cluster(cluster)
        self.assertEqual(self.gnm.cluster, "A")

    def test_set_cluster_2(self):
        """Check that 'singleton' string is set appropriately."""
        cluster = "singleton"
        self.gnm.set_cluster(cluster)
        self.assertEqual(self.gnm.cluster, "Singleton")

    def test_set_cluster_3(self):
        """Check that None is set appropriately."""
        cluster = None
        self.gnm.set_cluster(cluster)
        self.assertEqual(self.gnm.cluster, "Singleton")

    def test_set_cluster_4(self):
        """Check that whitespace is removed."""
        cluster = " A   "
        self.gnm.set_cluster(cluster)
        self.assertEqual(self.gnm.cluster, "A")




    def test_set_subcluster_1(self):
        """Check that standard Subcluster is set appropriately."""
        subcluster = "A2"
        self.gnm.set_subcluster(subcluster)
        self.assertEqual(self.gnm.subcluster, "A2")

    def test_set_subcluster_2(self):
        """Check that whitespace is removed."""
        subcluster = "    A2    "
        self.gnm.set_subcluster(subcluster)
        self.assertEqual(self.gnm.subcluster, "A2")

    def test_set_subcluster_3(self):
        """Check that 'none' subcluster is set appropriately from None."""
        subcluster = None
        self.gnm.set_subcluster(subcluster)
        self.assertEqual(self.gnm.subcluster, "none")




    def test_set_date_1(self):
        """Check that date is set appropriately."""
        date = datetime.strptime('1/1/1900', '%m/%d/%Y')
        self.gnm.set_date(date, "empty_string")
        self.assertEqual(self.gnm.date, date)

    def test_set_date_2(self):
        """Check that a None date is set appropriately
        from None."""
        date = None
        self.gnm.set_date(date, "empty_string")
        self.assertEqual(self.gnm.date, "")

    def test_set_date_3(self):
        """Check that filled date is set appropriately
        from None."""
        date1 = None
        date2 = datetime.strptime('1/1/0001', '%m/%d/%Y')
        self.gnm.set_date(date1, "empty_datetime_obj")
        self.assertEqual(self.gnm.date, date2)

    def test_set_date_4(self):
        """Check that None date is set appropriately
        from None."""
        date = None
        self.gnm.set_date(date, "none_object")
        self.assertIsNone(self.gnm.date)

    def test_set_date_5(self):
        """Check that None date is set appropriately
        from empty string."""
        date = ""
        self.gnm.set_date(date, "none_object")
        self.assertIsNone(self.gnm.date)

    def test_set_date_6(self):
        """Check that empty string date is set appropriately from
        incorrect strategy."""
        date = None
        self.gnm.set_date(date, "invalid")
        self.assertEqual(self.gnm.date, None)




    def test_set_annotation_author_1(self):
        """Check that annotation_author to set to 1."""
        value = 1
        self.gnm.set_annotation_author(value)
        self.assertEqual(self.gnm.annotation_author, 1)

    def test_set_annotation_author_2(self):
        """Check that annotation_author remains a string."""
        value = "retain"
        self.gnm.set_annotation_author(value)
        self.assertEqual(self.gnm.annotation_author, "retain")




    def test_set_retrieve_record_1(self):
        """Check that annotation_author to set to 1."""
        value = 1
        self.gnm.set_retrieve_record(value)
        self.assertEqual(self.gnm.retrieve_record, 1)

    def test_set_retrieve_record_2(self):
        """Check that annotation_author remains a string."""
        value = "retain"
        self.gnm.set_retrieve_record(value)
        self.assertEqual(self.gnm.retrieve_record, "retain")




    def test_set_cds_descriptions_1(self):
        """Check that descriptions are set from 'product' descriptions."""
        self.cds1.product = "lysB"
        self.cds1.function = "lysA"
        self.cds2.product = "rep"
        self.cds2.function = "repressor"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.set_cds_descriptions("product")
        with self.subTest():
            self.assertEqual(self.cds1.description, "lysB")
        with self.subTest():
            self.assertEqual(self.cds2.description, "rep")

    def test_set_cds_descriptions_2(self):
        """Check that descriptions are set from 'function' descriptions."""
        self.cds1.product = "lysB"
        self.cds1.function = "lysA"
        self.cds2.product = "rep"
        self.cds2.function = "repressor"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.set_cds_descriptions("function")
        with self.subTest():
            self.assertEqual(self.cds1.description, "lysA")
        with self.subTest():
            self.assertEqual(self.cds2.description, "repressor")




    def test_tally_cds_descriptions_1(self):
        """Check that all description tallies are reset to 0 and that
        no description tally is incremented."""
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm._cds_descriptions_tally = 10
        self.gnm._cds_products_tally = 10
        self.gnm._cds_functions_tally = 10
        self.gnm._cds_notes_tally = 10
        self.gnm.tally_cds_descriptions()
        with self.subTest():
            self.assertEqual(self.gnm._cds_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual(self.gnm._cds_products_tally, 0)
        with self.subTest():
            self.assertEqual(self.gnm._cds_functions_tally, 0)
        with self.subTest():
            self.assertEqual(self.gnm._cds_notes_tally, 0)

    def test_tally_cds_descriptions_2(self):
        """Check that processed primary description tally
        is incremented."""
        self.cds1.description = "abcd"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm._cds_descriptions_tally = 10
        self.gnm.tally_cds_descriptions()
        with self.subTest():
            self.assertEqual(self.gnm._cds_descriptions_tally, 1)
        with self.subTest():
            self.assertEqual(self.gnm._cds_products_tally, 0)
        with self.subTest():
            self.assertEqual(self.gnm._cds_functions_tally, 0)
        with self.subTest():
            self.assertEqual(self.gnm._cds_notes_tally, 0)

    def test_tally_cds_descriptions_3(self):
        """Check that processed product description tally
        is incremented."""
        self.cds1.product = "abcd"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.tally_cds_descriptions()

        with self.subTest():
            self.assertEqual(self.gnm._cds_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual(self.gnm._cds_products_tally, 1)
        with self.subTest():
            self.assertEqual(self.gnm._cds_functions_tally, 0)
        with self.subTest():
            self.assertEqual(self.gnm._cds_notes_tally, 0)

    def test_tally_cds_descriptions_4(self):
        """Check that processed function description tally
        is incremented."""
        self.cds1.function = "abcd"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.tally_cds_descriptions()

        with self.subTest():
            self.assertEqual(self.gnm._cds_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual(self.gnm._cds_products_tally, 0)
        with self.subTest():
            self.assertEqual(self.gnm._cds_functions_tally, 1)
        with self.subTest():
            self.assertEqual(self.gnm._cds_notes_tally, 0)

    def test_tally_cds_descriptions_5(self):
        """Check that processed note description tally
        is incremented."""
        self.cds1.note = "abcd"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.tally_cds_descriptions()

        with self.subTest():
            self.assertEqual(self.gnm._cds_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual(self.gnm._cds_products_tally, 0)
        with self.subTest():
            self.assertEqual(self.gnm._cds_functions_tally, 0)
        with self.subTest():
            self.assertEqual(self.gnm._cds_notes_tally, 1)

    def test_tally_cds_descriptions_6(self):
        """Check that all description tallies are incremented."""
        self.cds1.description = "abcd"
        self.cds1.product = "efgh"
        self.cds1.function = "ijkl"
        self.cds1.note = "mnop"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.tally_cds_descriptions()

        with self.subTest():
            self.assertEqual(self.gnm._cds_descriptions_tally, 1)
        with self.subTest():
            self.assertEqual(self.gnm._cds_products_tally, 1)
        with self.subTest():
            self.assertEqual(self.gnm._cds_functions_tally, 1)
        with self.subTest():
            self.assertEqual(self.gnm._cds_notes_tally, 1)

    def test_tally_cds_descriptions_7(self):
        """Check that all description tallies are incremented for
        all CDS features."""
        self.cds1.description = "abcd"
        self.cds1.product = "efgh"
        self.cds1.function = "ijkl"
        self.cds1.note = "mnop"

        self.cds2.description = "ab"
        self.cds2.product = "cd"
        self.cds2.function = "ef"
        self.cds2.note = "gh"

        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.tally_cds_descriptions()

        with self.subTest():
            self.assertEqual(self.gnm._cds_descriptions_tally, 2)
        with self.subTest():
            self.assertEqual(self.gnm._cds_products_tally, 2)
        with self.subTest():
            self.assertEqual(self.gnm._cds_functions_tally, 2)
        with self.subTest():
            self.assertEqual(self.gnm._cds_notes_tally, 2)




    def test_check_subcluster_structure_1(self):
        """Check that no error is produced if the
        non-empty subcluster is structured correctly."""
        self.gnm.subcluster = "A1"
        self.gnm.check_subcluster_structure("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_subcluster_structure_2(self):
        """Check that an error is produced if the
        non-empty subcluster is not structured correctly."""
        self.gnm.subcluster = "A"
        self.gnm.check_subcluster_structure()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_subcluster_structure_3(self):
        """Check that no error is produced if the
        subcluster is empty."""
        self.gnm.subcluster = "none"
        self.gnm.check_subcluster_structure()
        self.assertEqual(self.gnm.evaluations[0].status, "untested")




    def test_check_cluster_structure_1(self):
        """Check that no error is produced if the
        non-empty cluster is structured correctly."""
        self.gnm.cluster = "A"
        self.gnm.check_cluster_structure("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_cluster_structure_2(self):
        """Check that an error is produced if the
        non-empty cluster is not structured correctly."""
        self.gnm.cluster = "A1"
        self.gnm.check_cluster_structure()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_cluster_structure_3(self):
        """Check that no error is produced if the
        cluster is empty."""
        self.gnm.cluster = "none"
        self.gnm.check_cluster_structure()
        self.assertEqual(self.gnm.evaluations[0].status, "untested")




    def test_check_compatible_cluster_and_subcluster_1(self):
        """Check that compatible Cluster and subcluster
        do not produce an error."""
        self.gnm.cluster = "A"
        self.gnm.subcluster = "A1"
        self.gnm.check_compatible_cluster_and_subcluster("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_compatible_cluster_and_subcluster_2(self):
        """Check that incompatible Cluster and subcluster
        produce an error."""
        self.gnm.cluster = "A"
        self.gnm.subcluster = "B1"
        self.gnm.check_compatible_cluster_and_subcluster()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)




    def test_check_nucleotides_1(self):
        """All nucleotides are in the alphabet."""
        alphabet = set(["A","B","C"])
        self.gnm.seq = "AB"
        self.gnm.check_nucleotides(alphabet, "eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_nucleotides_2(self):
        """Some nucleotides are not in the alphabet."""
        alphabet = set(["A","B","C"])
        self.gnm.seq = "AD"
        self.gnm.check_nucleotides(alphabet)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)




    def test_check_magnitude_1(self):
        """Verify no error is produced when
        'length' is greater than 0, as expected."""
        self.gnm.length = 1000
        self.gnm.check_magnitude("length", ">", 0, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_magnitude_2(self):
        """Verify no error is produced when
        'length' is equal to 0, as expected."""
        self.gnm.length = 0
        self.gnm.check_magnitude("length", "=", 0)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_magnitude_3(self):
        """Verify no error is produced when
        'length' is less than 0, as expected."""
        self.gnm.length = -100
        self.gnm.check_magnitude("length", "<", 0)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_magnitude_4(self):
        """Verify an error is produced when
        'length' is greater than 0, unexpectedly."""
        self.gnm.length = 100
        self.gnm.check_magnitude("length", "=", 0)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_magnitude_5(self):
        """Verify an error is produced when
        'length' is less than 0, unexpectedly."""
        self.gnm.length = -100
        self.gnm.check_magnitude("length", ">", 0)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_magnitude_6(self):
        """Verify an error is produced when
        'name' is less than ref_name, unexpectedly."""
        self.gnm.name = "a"
        self.gnm.check_magnitude("name", "<", "c")
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_magnitude_7(self):
        """Verify an error is produced when
        'date' is greater than ref_date, unexpectedly."""
        ref_date =  datetime.strptime('1/1/2000', '%m/%d/%Y')
        self.gnm.date = datetime.strptime('2/1/2000', '%m/%d/%Y')
        self.gnm.check_magnitude("date", "<", ref_date)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_magnitude_8(self):
        """Verify an error is produced when
        'date' is equal to ref_date, unexpectedly."""
        ref_date =  datetime.strptime('2/1/2000', '%m/%d/%Y')
        self.gnm.date = datetime.strptime('2/1/2000', '%m/%d/%Y')
        self.gnm.check_magnitude("date", "=", ref_date)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")




    def test_check_authors_1(self):
        """Check that no warning is produced when author is expected
        and present."""
        check_set = set(["hatfull"])
        self.gnm.authors = "abcd; efgh,s,a,s; HATFULL,G.F; xyz"
        self.gnm.check_authors(check_set=check_set, eval_id="eval_id")

        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_authors_2(self):
        """Check that no warning is produced when author is expected
        and present not separated only by whitespace."""
        check_set = set(["hatfull"])
        self.gnm.authors = "abcd; efgh,s,a,s HATFULL djkf,G.F; xyz"
        self.gnm.check_authors(check_set=check_set, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_authors_3(self):
        """Check that no warning is produced when author is not expected
        and not present."""
        check_set = set(["hatfull"])
        self.gnm.authors = "abcd; efgh; xyz"
        self.gnm.check_authors(check_set=check_set, expect=False)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_authors_4(self):
        """Check that a warning is produced when author is expected
        and not present."""
        check_set = set(["hatfull"])
        self.gnm.authors = "abcd; efgh; xyz"
        self.gnm.check_authors(check_set=check_set)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_authors_5(self):
        """Check that a warning is produced when author is not expected
        and present."""
        check_set = set(["hatfull"])
        self.gnm.authors = "abcd; efgh; HATFULL; xyz"
        self.gnm.annotation_author = 0
        self.gnm.check_authors(check_set=check_set, expect=False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_authors_6(self):
        """Check that no warning is produced when author is expected
        and present in a supplied set of multiple authors."""
        check_set = set(["abcd", "efgh", "hatfull"])
        self.gnm.authors = "hatfull,g,f,s; 1234; xyz"
        self.gnm.check_authors(check_set=check_set)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")




    def test_set_unique_cds_start_end_ids_1(self):
        """Verify that both sets are computed."""
        self.gnm._cds_start_end_ids = \
            [(1, 5), (2, 10), (10, 2), (2, 10)]
        expected_unique_set = set([(1, 5), (10, 2)])
        expected_duplicate_set = set([(2, 10)])
        self.gnm.set_unique_cds_start_end_ids()
        with self.subTest():
            self.assertEqual(self.gnm._cds_unique_start_end_ids, \
                expected_unique_set)
        with self.subTest():
            self.assertEqual(self.gnm._cds_duplicate_start_end_ids, \
                expected_duplicate_set)




    def test_set_unique_cds_end_orient_ids_1(self):
        """Verify that both sets are computed."""
        self.gnm._cds_end_orient_ids = \
            [(1, "forward"), (2, "reverse"), (2, "forward"), (2, "reverse")]
        expected_unique_set = set([(1, "forward"), (2, "forward")])
        expected_duplicate_set = set([(2, "reverse")])
        self.gnm.set_unique_cds_end_orient_ids()
        with self.subTest():
            self.assertEqual(self.gnm._cds_unique_end_orient_ids,
                expected_unique_set)
        with self.subTest():
            self.assertEqual(self.gnm._cds_duplicate_end_orient_ids,
                expected_duplicate_set)




    def test_check_cds_start_end_ids_1(self):
        """Verify that no warning is produced."""
        self.gnm.check_cds_start_end_ids("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_cds_start_end_ids_2(self):
        """Verify that a warning is produced."""
        self.gnm._cds_duplicate_start_end_ids = set([(2, 10)])
        self.gnm.check_cds_start_end_ids()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)




    def test_check_cds_end_orient_ids_1(self):
        """Verify that no warning is produced."""
        self.gnm.check_cds_end_orient_ids("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_cds_end_orient_ids_2(self):
        """Verify that a warning is produced."""
        self.gnm._cds_duplicate_end_orient_ids = set([(2, "forward")])
        self.gnm.check_cds_end_orient_ids()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)




    def test_clear_locus_tags_1(self):
        """Verify that each locus tag is cleared."""
        self.cds1.locus_tag = "L5_1"
        self.cds2.locus_tag = "L5_2"
        self.cds3.locus_tag = "L5_3"
        self.gnm.cds_features = [self.cds1, self.cds2, self.cds3]

        self.trna1.locus_tag = "L5_4"
        self.gnm.trna_features = [self.trna1]

        self.tmrna1.locus_tag = "L5_5"
        self.tmrna2.locus_tag = "L5_6"
        self.gnm.tmrna_features = [self.tmrna1, self.tmrna2]

        not_empty = 0
        for ftr in self.gnm.cds_features:
            if ftr.locus_tag != "":
                not_empty += 1
        for ftr in self.gnm.trna_features:
            if ftr.locus_tag != "":
                not_empty += 1
        for ftr in self.gnm.tmrna_features:
            if ftr.locus_tag != "":
                not_empty += 1

        self.gnm.clear_locus_tags()
        empty = 0
        for ftr in self.gnm.cds_features:
            if ftr.locus_tag == "":
                empty += 1
        for ftr in self.gnm.trna_features:
            if ftr.locus_tag == "":
                empty += 1
        for ftr in self.gnm.tmrna_features:
            if ftr.locus_tag == "":
                empty += 1

        self.assertEqual(not_empty, empty)




    def test_set_feature_genome_ids_1(self):
        """Verify that the genome_id of CDS, tRNA, tmRNA, and Source features
        are all set."""
        id1 = "L5"
        id2 = "Trixie"
        self.cds1.genome_id = id1
        self.cds2.genome_id = id1
        self.gnm.cds_features = [self.cds1, self.cds2]

        self.trna1.genome_id = id1
        self.trna2.genome_id = id1
        self.gnm.trna_features = [self.trna1, self.trna2]

        self.tmrna1.genome_id = id1
        self.tmrna2.genome_id = id1
        self.gnm.tmrna_features = [self.tmrna1, self.tmrna2]

        self.source1.genome_id = id1
        self.gnm.source_features = [self.source1]

        self.gnm.set_feature_genome_ids(use_cds=True, use_trna=True,
                                        use_tmrna=True, use_source=True,
                                        value=id2)

        for ftr in self.gnm.cds_features:
            with self.subTest():
                self.assertEqual(ftr.genome_id, id2)
        for ftr in self.gnm.trna_features:
            with self.subTest():
                self.assertEqual(ftr.genome_id, id2)
        for ftr in self.gnm.tmrna_features:
            with self.subTest():
                self.assertEqual(ftr.genome_id, id2)
        for ftr in self.gnm.source_features:
            with self.subTest():
                self.assertEqual(ftr.genome_id, id2)

    def test_set_feature_genome_ids_2(self):
        """Verify that the genome_id of only CDS features are set."""
        id1 = "L5"
        id2 = "Trixie"
        self.cds1.genome_id = id1
        self.cds2.genome_id = id1
        self.gnm.cds_features = [self.cds1, self.cds2]

        self.trna1.genome_id = id1
        self.trna2.genome_id = id1
        self.gnm.trna_features = [self.trna1, self.trna2]

        self.tmrna1.genome_id = id1
        self.tmrna2.genome_id = id1
        self.gnm.tmrna_features = [self.tmrna1, self.tmrna2]

        self.source1.genome_id = id1
        self.gnm.source_features = [self.source1]

        self.gnm.set_feature_genome_ids(use_cds=True, use_trna=False,
                                        use_tmrna=False, use_source=False,
                                        value=id2)

        for ftr in self.gnm.cds_features:
            with self.subTest():
                self.assertEqual(ftr.genome_id, id2)
        for ftr in self.gnm.trna_features:
            with self.subTest():
                self.assertEqual(ftr.genome_id, id1)
        for ftr in self.gnm.tmrna_features:
            with self.subTest():
                self.assertEqual(ftr.genome_id, id1)
        for ftr in self.gnm.source_features:
            with self.subTest():
                self.assertEqual(ftr.genome_id, id1)

    def test_set_feature_genome_ids_3(self):
        """Verify that the genome_id of only CDS features are set without
        providing a value."""
        id1 = "L5"
        id2 = "Trixie"
        # Set genome id to 'Trixie' and set all feature.genome_ids to 'L5'.
        self.gnm.id = id2
        self.cds1.genome_id = id1
        self.cds2.genome_id = id1
        self.gnm.cds_features = [self.cds1, self.cds2]

        self.trna1.genome_id = id1
        self.trna2.genome_id = id1
        self.gnm.trna_features = [self.trna1, self.trna2]

        self.tmrna1.genome_id = id1
        self.tmrna2.genome_id = id1
        self.gnm.tmrna_features = [self.tmrna1, self.tmrna2]

        self.source1.genome_id = id1
        self.gnm.source_features = [self.source1]

        self.gnm.set_feature_genome_ids(use_cds=True, use_trna=False,
                                        use_tmrna=False, use_source=False)

        for ftr in self.gnm.cds_features:
            with self.subTest():
                self.assertEqual(ftr.genome_id, id2)
        for ftr in self.gnm.trna_features:
            with self.subTest():
                self.assertEqual(ftr.genome_id, id1)
        for ftr in self.gnm.tmrna_features:
            with self.subTest():
                self.assertEqual(ftr.genome_id, id1)
        for ftr in self.gnm.source_features:
            with self.subTest():
                self.assertEqual(ftr.genome_id, id1)




    def test_check_attribute_1(self):
        """Verify no error is produced when the id
        is in the check_set and is expected to be in the set."""
        check_set = set(["Trixie", "L5"])
        self.gnm.id = "Trixie"
        self.gnm.check_attribute("id", check_set, True, "eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_attribute_2(self):
        """Verify an error is produced when the id
        is not in the check_set and is expected to be in the set."""
        check_set = set(["Trixie", "L5"])
        self.gnm.id = "D29"
        self.gnm.check_attribute("id", check_set, True)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_attribute_3(self):
        """Verify no error is produced when the id
        is not in the check_set and is not expected to be in the set."""
        check_set = set(["Trixie", "L5"])
        self.gnm.id = "D29"
        self.gnm.check_attribute("id", check_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_attribute_4(self):
        """Verify an error is produced when the id
        is in the check_set and is not expected to be in the set."""
        check_set = set(["Trixie", "L5"])
        self.gnm.id = "Trixie"
        self.gnm.check_attribute("id", check_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_attribute_5(self):
        """Verify an error is produced when the host_genus
        is in the check_set and is not expected to be in the set."""
        check_set = set(["Mycobacterium", "Gordonia"])
        self.gnm.host_genus = "Mycobacterium"
        self.gnm.check_attribute("host_genus", check_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_attribute_6(self):
        """Verify an error is produced when the seq
        is in the check_set and is expected to be in the set."""
        check_set = set([Seq("AAAA", IUPAC.ambiguous_dna),
                         Seq("ATAT", IUPAC.ambiguous_dna)])
        self.gnm.seq = Seq("AAAA", IUPAC.ambiguous_dna)
        self.gnm.check_attribute("seq", check_set, True)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_attribute_7(self):
        """Verify an error is produced when the retrieve_record
        is not in the check_set and is expected to be in the set."""
        check_set = set([1, 0])
        self.gnm.retrieve_record = -1
        self.gnm.check_attribute("retrieve_record", check_set, True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_attribute_8(self):
        """Verify no test is performed when the attribute is invalid."""
        check_set = set([1, 0])
        self.gnm.check_attribute("invalid", check_set, True)
        self.assertEqual(self.gnm.evaluations[0].status, "untested")




    def test_compare_two_attributes_1(self):
        """Verify no error is produced when both attributes have
        identical values as expected."""
        self.gnm.id = "Trixie"
        self.gnm.name = "Trixie"
        self.gnm.compare_two_attributes("id", "name", expect_same=True,
                                           eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_compare_two_attributes_2(self):
        """Verify no error is produced when both attributes have
        different values as expected."""
        self.gnm.id = "Trixie"
        self.gnm.name = "Trixie_Draft"
        self.gnm.compare_two_attributes("id", "name", expect_same=False)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_compare_two_attributes_3(self):
        """Verify an error is produced when both attributes have
        different values unexpectedly."""
        self.gnm.id = "Trixie"
        self.gnm.name = "Trixie_Draft"
        self.gnm.compare_two_attributes("id", "name", expect_same=True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_compare_two_attributes_4(self):
        """Verify an error is produced when both attributes have
        identical values unexpectedly."""
        self.gnm.id = "Trixie"
        self.gnm.name = "Trixie"
        self.gnm.compare_two_attributes("id", "name", expect_same=False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_compare_two_attributes_5(self):
        """Verify no error is produced when both int attributes have
        different values as expected."""
        self.gnm.retrieve_record = 1
        self.gnm.annotation_author = 0
        self.gnm.compare_two_attributes("retrieve_record",
                                        "annotation_author", expect_same=False)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_compare_two_attributes_6(self):
        """Verify no error is produced when two attributes of
        different data types have different values as expected."""
        self.gnm.id = "Trixie"
        self.gnm.seq = Seq("AAAA", IUPAC.ambiguous_dna)
        self.gnm.compare_two_attributes("id", "seq", expect_same=False)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_compare_two_attributes_7(self):
        """Verify no test is performed when the first attribute is invalid."""
        self.gnm.id = "Trixie"
        self.gnm.name = "Trixie_Draft"
        self.gnm.compare_two_attributes("invalid", "id", expect_same=False)
        self.assertEqual(self.gnm.evaluations[0].status, "untested")

    def test_compare_two_attributes_8(self):
        """Verify no test is performed when the second attribute is invalid."""
        self.gnm.id = "Trixie"
        self.gnm.name = "Trixie_Draft"
        self.gnm.compare_two_attributes("id", "invalid", expect_same=False)
        self.assertEqual(self.gnm.evaluations[0].status, "untested")




    def test_parse_description_1(self):
        """Verify empty string is parsed correctly."""
        self.gnm.description = ""
        self.gnm.parse_description()
        expected_phage = ""
        expected_host = ""
        with self.subTest():
            self.assertEqual(
                self.gnm._description_name, expected_phage)
        with self.subTest():
            self.assertEqual(
                self.gnm._description_host_genus, expected_host)

    def test_parse_description_2(self):
        """Verify string is parsed correctly."""
        self.gnm.description = "asdf Mycobacterium phage Trixie."
        self.gnm.parse_description()
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        with self.subTest():
            self.assertEqual(
                self.gnm._description_name, expected_phage)
        with self.subTest():
            self.assertEqual(
                self.gnm._description_host_genus, expected_host)




    def test_parse_source_1(self):
        """Verify empty string is parsed correctly."""
        self.gnm.source = ""
        self.gnm.parse_source()
        expected_phage = ""
        expected_host = ""
        with self.subTest():
            self.assertEqual(
                self.gnm._source_name, expected_phage)
        with self.subTest():
            self.assertEqual(
                self.gnm._source_host_genus, expected_host)

    def test_parse_source_2(self):
        """Verify string is parsed correctly."""
        self.gnm.source = "asdf Mycobacterium phage Trixie."
        self.gnm.parse_source()
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        with self.subTest():
            self.assertEqual(
                self.gnm._source_name, expected_phage)
        with self.subTest():
            self.assertEqual(
                self.gnm._source_host_genus, expected_host)




    def test_parse_organism_1(self):
        """Verify empty string is parsed correctly."""
        self.gnm.organism = ""
        self.gnm.parse_organism()
        expected_phage = ""
        expected_host = ""
        with self.subTest():
            self.assertEqual(
                self.gnm._organism_name, expected_phage)
        with self.subTest():
            self.assertEqual(
                self.gnm._organism_host_genus, expected_host)

    def test_parse_organism_2(self):
        """Verify string is parsed correctly."""
        self.gnm.organism = "asdf Mycobacterium phage Trixie."
        self.gnm.parse_organism()
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        with self.subTest():
            self.assertEqual(
                self.gnm._organism_name, expected_phage)
        with self.subTest():
            self.assertEqual(
                self.gnm._organism_host_genus, expected_host)




    def test_check_feature_coordinates_1(self):
        """Verify no error is produced by two CDS features with same orientation."""
        self.cds1.orientation = "F"
        self.cds1.start = 5
        self.cds1.stop = 50
        self.cds2.orientation = "F"
        self.cds2.start = 20
        self.cds2.stop = 70
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_coordinates(use_cds=True,eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_feature_coordinates_2(self):
        """Verify an error is produced by two CDS features with same orientation
        with identical start and stop coordinates."""
        self.cds1.orientation = "F"
        self.cds1.start = 5
        self.cds1.stop = 50
        self.cds2.orientation = "F"
        self.cds2.start = 5
        self.cds2.stop = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_coordinates(use_cds=True)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_feature_coordinates_3(self):
        """Verify an error is produced by two CDS features with different orientation
        with identical start and stop coordinates."""
        self.cds1.orientation = "R"
        self.cds1.start = 5
        self.cds1.stop = 50
        self.cds2.orientation = "F"
        self.cds2.start = 5
        self.cds2.stop = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_coordinates(use_cds=True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_feature_coordinates_4(self):
        """Verify no error is produced by two CDS features with different orientation
        with identical start and stop coordinates when strand is True."""
        self.cds1.orientation = "R"
        self.cds1.start = 5
        self.cds1.stop = 50
        self.cds2.orientation = "F"
        self.cds2.start = 5
        self.cds2.stop = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_coordinates(use_cds=True, strand=True)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_feature_coordinates_5(self):
        """Verify an error is produced by two CDS features with same orientation
        when they have nested coordinates."""
        self.cds1.orientation = "F"
        self.cds1.start = 10
        self.cds1.stop = 20
        self.cds2.orientation = "F"
        self.cds2.start = 5
        self.cds2.stop = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_coordinates(use_cds=True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")


    def test_check_feature_coordinates_6(self):
        """Verify no error is produced by two CDS features with "F" orientation
        when they have the same start coordinates."""
        self.cds1.orientation = "F"
        self.cds1.start = 10
        self.cds1.stop = 20
        self.cds2.orientation = "F"
        self.cds2.start = 10
        self.cds2.stop = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_coordinates(use_cds=True)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_feature_coordinates_7(self):
        """Verify an error is produced by two CDS features with "F" orientation
        when they have the same stop coordinates."""
        self.cds1.orientation = "F"
        self.cds1.start = 10
        self.cds1.stop = 50
        self.cds2.orientation = "F"
        self.cds2.start = 5
        self.cds2.stop = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_coordinates(use_cds=True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_feature_coordinates_8(self):
        """Verify no error is produced by two CDS features with "R" orientation
        when they have the same stop coordinates."""
        self.cds1.orientation = "R"
        self.cds1.start = 10
        self.cds1.stop = 50
        self.cds2.orientation = "R"
        self.cds2.start = 5
        self.cds2.stop = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_coordinates(use_cds=True)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_feature_coordinates_9(self):
        """Verify no error is produced by two CDS features with "R" orientation
        when they have the same start coordinates."""
        self.cds1.orientation = "R"
        self.cds1.start = 5
        self.cds1.stop = 20
        self.cds2.orientation = "R"
        self.cds2.start = 5
        self.cds2.stop = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_coordinates(use_cds=True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_feature_coordinates_10(self):
        """Verify an error is produced by a CDS and tRNA feature
        when they have the same coordinates."""
        self.cds1.orientation = "F"
        self.cds1.start = 5
        self.cds1.stop = 20
        self.cds2.orientation = "F"
        self.cds2.start = 100
        self.cds2.stop = 200
        self.trna1.start = 5
        self.trna1.stop = 20
        self.trna1.orientation = "F"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.trna_features = [self.trna1]
        self.gnm.check_feature_coordinates(use_cds=True, use_trna=True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_feature_coordinates_11(self):
        """Verify an error is produced by a tRNA and tmRNA feature
        when they have the same coordinates."""
        self.tmrna1.orientation = "F"
        self.tmrna1.start = 5
        self.tmrna1.stop = 20
        self.trna1.start = 5
        self.trna1.stop = 20
        self.trna1.orientation = "F"
        self.gnm.trna_features = [self.trna1]
        self.gnm.tmrna_features = [self.tmrna1]
        self.gnm.check_feature_coordinates(use_tmrna=True, use_trna=True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_feature_coordinates_12(self):
        """Verify no error is produced by a CDS and tRNA feature
        when they have the same coordinates when CDS is false."""
        self.cds1.orientation = "F"
        self.cds1.start = 5
        self.cds1.stop = 20
        self.cds2.orientation = "F"
        self.cds2.start = 100
        self.cds2.stop = 200
        self.trna1.start = 5
        self.trna1.stop = 20
        self.trna1.orientation = "F"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.trna_features = [self.trna1]
        self.gnm.check_feature_coordinates(use_trna=True)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_feature_coordinates_13(self):
        """Verify an error is produced by a CDS feature stored in the
        CDS features list and a CDS feature stored in a supplementary list
        when they have the same coordinates."""
        self.cds1.orientation = "F"
        self.cds1.start = 5
        self.cds1.stop = 20
        self.cds2.orientation = "F"
        self.cds2.start = 100
        self.cds2.stop = 200
        self.cds3.orientation = "F"
        self.cds3.start = 5
        self.cds3.stop = 20
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_coordinates(use_cds=True, other=[self.cds3])
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_feature_coordinates_14(self):
        """Verify no error is produced when empty lists are passed through."""
        self.gnm.check_feature_coordinates(use_cds=True, use_trna=True,
                                           use_tmrna=True, other=[])
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_feature_coordinates_15(self):
        """Verify no error is produced by two CDS features with same orientation
        when one wraps around the end of the genome."""
        # For wrap-around genes, the 'start' coordinate is larger than
        # the 'stop' coordinate.
        self.cds1.orientation = "F"
        self.cds1.start = 50000
        self.cds1.stop = 20
        self.cds2.orientation = "F"
        self.cds2.start = 5
        self.cds2.stop = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_coordinates(use_cds=True)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")




class TestGenomeClass2(unittest.TestCase):

    def setUp(self):
        self.gnm = genome.Genome()
        self.gnm.id = "L5"

        self.cds1 = cds.Cds()
        self.cds2 = cds.Cds()
        self.cds3 = cds.Cds()
        self.cds4 = cds.Cds()

        self.trna1 = trna.TrnaFeature()
        self.trna2 = trna.TrnaFeature()
        self.trna3 = trna.TrnaFeature()
        self.trna4 = trna.TrnaFeature()

        self.src1 = source.Source()
        self.src2 = source.Source()

        self.tmrna1 = tmrna.TmrnaFeature()
        self.tmrna2 = tmrna.TmrnaFeature()

        # 1
        self.trna3.start = 10
        self.trna3.stop = 15

        # 2
        self.trna2.start = 10
        self.trna2.stop = 17

        # 3
        self.cds4.start = 10
        self.cds4.stop = 20

        # 4
        self.cds3.start = 15
        self.cds3.stop = 30

        # 5
        self.cds2.start = 15
        self.cds2.stop = 40

        # 6
        self.tmrna1.start = 17
        self.tmrna1.stop = 18

        # 7
        self.src2.start = 19
        self.src2.stop = 35

        # 8
        self.trna1.start = 100
        self.trna1.stop = 200

        # 9
        self.src1.start = 101
        self.src1.stop = 300

        # 10
        self.tmrna2.start = 200
        self.tmrna2.stop = 205

        # 11 Wrap-around feature.
        self.cds1.start = 400
        self.cds1.stop = 3

        self.gnm.cds_features = [self.cds1, self.cds2, self.cds3, self.cds4]
        self.gnm.trna_features = [self.trna1, self.trna2, self.trna3]
        self.gnm.tmrna_features = [self.tmrna1, self.tmrna2]
        self.gnm.source_features = [self.src1, self.src2]


    def test_set_feature_ids_1(self):
        """Verify that CDS features are sorted in correct order."""
        self.gnm.set_feature_ids(use_cds=True)
        with self.subTest():
            self.assertEqual(self.cds1.id, "L5_4")
        with self.subTest():
            self.assertEqual(self.cds2.id, "L5_3")
        with self.subTest():
            self.assertEqual(self.cds3.id, "L5_2")
        with self.subTest():
            self.assertEqual(self.cds4.id, "L5_1")

    def test_set_feature_ids_2(self):
        """Verify that tRNA features are sorted in correct order."""
        self.gnm.set_feature_ids(use_trna=True)
        with self.subTest():
            self.assertEqual(self.trna1.id, "L5_3")
        with self.subTest():
            self.assertEqual(self.trna2.id, "L5_2")
        with self.subTest():
            self.assertEqual(self.trna3.id, "L5_1")

    def test_set_feature_ids_3(self):
        """Verify that source features are sorted in correct order."""
        self.gnm.set_feature_ids(use_source=True)
        with self.subTest():
            self.assertEqual(self.src2.id, "L5_1")
        with self.subTest():
            self.assertEqual(self.src1.id, "L5_2")

    def test_set_feature_ids_4(self):
        """Verify that tmRNA features are sorted in correct order."""
        self.gnm.set_feature_ids(use_tmrna=True)
        with self.subTest():
            self.assertEqual(self.tmrna1.id, "L5_1")
        with self.subTest():
            self.assertEqual(self.tmrna2.id, "L5_2")

    def test_set_feature_ids_5(self):
        """Verify that CDS, tRNA, tmrna, and source features are sorted
        in correct order."""
        self.gnm.set_feature_ids(use_cds=True, use_trna=True, use_tmrna=True,
                                 use_source=True)
        with self.subTest():
            self.assertEqual(self.trna1.id, "L5_8")
        with self.subTest():
            self.assertEqual(self.trna2.id, "L5_2")
        with self.subTest():
            self.assertEqual(self.trna3.id, "L5_1")
        with self.subTest():
            self.assertEqual(self.cds1.id, "L5_11")
        with self.subTest():
            self.assertEqual(self.cds2.id, "L5_5")
        with self.subTest():
            self.assertEqual(self.cds3.id, "L5_4")
        with self.subTest():
            self.assertEqual(self.cds4.id, "L5_3")
        with self.subTest():
            self.assertEqual(self.src1.id, "L5_9")
        with self.subTest():
            self.assertEqual(self.src2.id, "L5_7")
        with self.subTest():
            self.assertEqual(self.tmrna1.id, "L5_6")
        with self.subTest():
            self.assertEqual(self.tmrna2.id, "L5_10")

    def test_set_feature_ids_6(self):
        """Verify that CDS, tRNA, tmRNA, and source features are sorted
        in correct order with type delimiter added."""
        self.gnm.cds_features = [self.cds4]
        self.gnm.trna_features = [self.trna3]
        self.gnm.tmrna_features = [self.tmrna1]
        self.gnm.source_features = [self.src2]
        self.gnm.set_feature_ids(use_type=True, use_cds=True, use_trna=True,
                                 use_tmrna=True, use_source=True)
        with self.subTest():
            self.assertEqual(self.trna3.id, "L5_TRNA_1")
        with self.subTest():
            self.assertEqual(self.cds4.id, "L5_CDS_2")
        with self.subTest():
            self.assertEqual(self.tmrna1.id, "L5_TMRNA_3")
        with self.subTest():
            self.assertEqual(self.src2.id, "L5_SRC_4")


if __name__ == '__main__':
    unittest.main()
