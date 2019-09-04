""" Unit tests for the Genome class."""


import unittest
from constants import constants
from classes import genome
from classes import cds
from datetime import datetime
from Bio.Seq import Seq
from classes import trna


class TestGenomeClass(unittest.TestCase):


    def setUp(self):
        self.gnm = genome.Genome()



        self.cds1 = cds.Cds()
        self.cds1.processed_description = ""
        self.cds1.processed_product = ""
        self.cds1.processed_function = ""
        self.cds1.processed_note = ""

        self.cds2 = cds.Cds()
        self.cds2.processed_description = ""
        self.cds2.processed_product = ""
        self.cds2.processed_function = ""
        self.cds2.processed_note = ""


        self.cds3 = cds.Cds()
        self.cds4 = cds.Cds()

        self.trna1 = trna.TrnaFeature()
        self.trna2 = trna.TrnaFeature()
        self.trna3 = trna.TrnaFeature()
        self.trna4 = trna.TrnaFeature()

    def test_set_filename_1(self):
        """Confirm file path is split appropriately."""
        filepath = "/path/to/folder/Trixie.gbk"
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
            self.assertEqual(self.gnm._length, 7)
        with self.subTest():
            self.assertEqual(self.gnm._gc, 57.1429)

    def test_set_sequence_2(self):
        """Check that sequence is set appropriately if it is an empty string."""
        seq = ""
        self.gnm.set_sequence(seq)
        with self.subTest():
            self.assertEqual(self.gnm.seq, "")
        with self.subTest():
            self.assertIsInstance(self.gnm.seq, Seq)
        with self.subTest():
            self.assertEqual(self.gnm._length, 0)
        with self.subTest():
            self.assertEqual(self.gnm._gc, -1)






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
        cds1._end_strand_id = (5, "forward")

        cds2 = cds.Cds()
        cds2._start_end_id = (21, 2)
        cds2._end_strand_id = (2, "reverse")

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
                self.gnm._cds_end_strand_ids, end_strand_id_list)




    def test_set_trna_features_1(self):
        """Check that tRNA feature list is set and length is computed."""
        features_list = [0,1,2,3]
        self.gnm.set_trna_features(features_list)
        with self.subTest():
            self.assertEqual(len(self.gnm.trna_features), 4)
        with self.subTest():
            self.assertEqual(self.gnm._trna_features_tally, 4)




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
        """Check that the id is set from the accession attribute."""
        self.gnm.accession = "Trixie_Draft"
        self.gnm.set_id(attribute="accession")
        self.assertEqual(self.gnm.id, "Trixie")

    def test_set_id_5(self):
        """Check that the id is set from the description attribute."""
        self.gnm.description = "Trixie_Draft"
        self.gnm.set_id(attribute="description")
        self.assertEqual(self.gnm.id, "Trixie")

    def test_set_id_6(self):
        """Check that the id is set from the source attribute."""
        self.gnm.source = "Trixie_Draft"
        self.gnm.set_id(attribute="source")
        self.assertEqual(self.gnm.id, "Trixie")

    def test_set_id_7(self):
        """Check that the id is set from the organism attribute."""
        self.gnm.organism = "Trixie_Draft"
        self.gnm.set_id(attribute="organism")
        self.assertEqual(self.gnm.id, "Trixie")

    def test_set_id_8(self):
        """Check that the id is set from the filename attribute."""
        self.gnm.filename = "Trixie_Draft"
        self.gnm.set_id(attribute="filename")
        self.assertEqual(self.gnm.id, "Trixie")

    def test_set_id_9(self):
        """Check that the id is set from the
        description_name attribute."""
        self.gnm._description_name = "Trixie_Draft"
        self.gnm.set_id(attribute="description_name")
        self.assertEqual(self.gnm.id, "Trixie")

    def test_set_id_10(self):
        """Check that the id is set from the
        source_name attribute."""
        self.gnm._source_name = "Trixie_Draft"
        self.gnm.set_id(attribute="source_name")
        self.assertEqual(self.gnm.id, "Trixie")

    def test_set_id_11(self):
        """Check that the id is set from the
        organism_name attribute."""
        self.gnm._organism_name = "Trixie_Draft"
        self.gnm.set_id(attribute="organism_name")
        self.assertEqual(self.gnm.id, "Trixie")

    def test_set_id_12(self):
        """Check that the id is empty from an invalid attribute."""
        self.gnm.id = "not empty"
        self.gnm.set_id(attribute="invalid")
        self.assertEqual(self.gnm.id, "")

    def test_set_id_13(self):
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
        """Check that the host_genus is set from the description attribute."""
        self.gnm.description = "Mycobacterium smegmatis"
        self.gnm.set_host_genus(attribute="description")
        self.assertEqual(self.gnm.host_genus, "Mycobacterium")

    def test_set_host_genus_7(self):
        """Check that the host_genus is set from the source attribute."""
        self.gnm.source = "Mycobacterium smegmatis"
        self.gnm.set_host_genus(attribute="source")
        self.assertEqual(self.gnm.host_genus, "Mycobacterium")

    def test_set_host_genus_8(self):
        """Check that the host_genus is set from the organism attribute."""
        self.gnm.organism = "Mycobacterium smegmatis"
        self.gnm.set_host_genus(attribute="organism")
        self.assertEqual(self.gnm.host_genus, "Mycobacterium")

    def test_set_host_genus_9(self):
        """Check that the host_genus is set from the filename attribute."""
        self.gnm.filename = "Mycobacterium smegmatis"
        self.gnm.set_host_genus(attribute="filename")
        self.assertEqual(self.gnm.host_genus, "Mycobacterium")

    def test_set_host_genus_10(self):
        """Check that the host_genus is set from the
        description_host_genus attribute."""
        self.gnm._description_host_genus = "Mycobacterium smegmatis"
        self.gnm.set_host_genus(attribute="description_host_genus")
        self.assertEqual(self.gnm.host_genus, "Mycobacterium")

    def test_set_host_genus_11(self):
        """Check that the host_genus is set from the
        source_host_genus attribute."""
        self.gnm._source_host_genus = "Mycobacterium smegmatis"
        self.gnm.set_host_genus(attribute="source_host_genus")
        self.assertEqual(self.gnm.host_genus, "Mycobacterium")

    def test_set_host_genus_12(self):
        """Check that the host_genus is set from the
        organism_host_genus attribute."""
        self.gnm._organism_host_genus = "Mycobacterium smegmatis"
        self.gnm.set_host_genus(attribute="organism_host_genus")
        self.assertEqual(self.gnm.host_genus, "Mycobacterium")

    def test_set_host_genus_13(self):
        """Check that the host_genus is empty from an invalid attribute."""
        self.gnm.set_host_genus(attribute="invalid")
        self.assertEqual(self.gnm.host_genus, "")

    def test_set_host_genus_14(self):
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
        cluster = "Singleton"
        self.gnm.set_cluster(cluster)
        self.assertEqual(self.gnm.cluster, "singleton")

    def test_set_cluster_3(self):
        """Check that None is set appropriately."""
        cluster = None
        self.gnm.set_cluster(cluster)
        self.assertEqual(self.gnm.cluster, "singleton")

    def test_set_cluster_4(self):
        """Check that whitespace is removed."""
        cluster = " A   "
        self.gnm.set_cluster(cluster)
        self.assertEqual(self.gnm.cluster, "A")




    def test_set_subcluster_1(self):
        """Check that standard Subcluster is set appropriately."""
        subcluster = "A2"
        self.gnm.set_subcluster(subcluster, "none_string")
        self.assertEqual(self.gnm.subcluster, "A2")

    def test_set_subcluster_2(self):
        """Check that whitespace is removed."""
        subcluster = "    A2    "
        self.gnm.set_subcluster(subcluster, "none_string")
        self.assertEqual(self.gnm.subcluster, "A2")

    def test_set_subcluster_3(self):
        """Check that None subcluster is set appropriately from None."""
        subcluster = None
        self.gnm.set_subcluster(subcluster, "none_object")
        self.assertIsNone(self.gnm.subcluster)

    def test_set_subcluster_4(self):
        """Check that none_string subcluster is set appropriately from None."""
        subcluster = None
        self.gnm.set_subcluster(subcluster, "none_string")
        self.assertEqual(self.gnm.subcluster, "none")

    def test_set_subcluster_5(self):
        """Check that empty_string subcluster is set appropriately from
        none_string."""
        subcluster = "none"
        self.gnm.set_subcluster(subcluster, "empty_string")
        self.assertEqual(self.gnm.subcluster, "")

    def test_set_subcluster_6(self):
        """Check that case is accounted for."""
        subcluster = "NONE"
        self.gnm.set_subcluster(subcluster, "empty_string")
        self.assertEqual(self.gnm.subcluster, "")




    def test_set_cluster_subcluster_1(self):
        """Check that None Cluster is set as singleton cluster_subcluster."""
        self.gnm.subcluster = ""
        self.gnm.cluster = None
        self.gnm.set_cluster_subcluster()
        self.assertEqual(self.gnm.cluster_subcluster, "singleton")

    def test_set_cluster_subcluster_2(self):
        """Check that singleton Cluster is set as
        singleton cluster_subcluster."""
        self.gnm.subcluster = ""
        self.gnm.cluster = "singleton"
        self.gnm.set_cluster_subcluster()
        self.assertEqual(self.gnm.cluster_subcluster, "singleton")

    def test_set_cluster_subcluster_3(self):
        """Check that Cluster is set as cluster_subcluster."""
        self.gnm.subcluster = ""
        self.gnm.cluster = "A"
        self.gnm.set_cluster_subcluster()
        self.assertEqual(self.gnm.cluster_subcluster, "A")

    def test_set_cluster_subcluster_4(self):
        """Check that Subcluster is set as cluster_subcluster."""
        self.gnm.subcluster = "A1"
        self.gnm.cluster = "A"
        self.gnm.set_cluster_subcluster()
        self.assertEqual(self.gnm.cluster_subcluster, "A1")

    def test_set_cluster_subcluster_5(self):
        """Check that Cluster is set when subcluster is None."""
        self.gnm.subcluster = None
        self.gnm.cluster = "A"
        self.gnm.set_cluster_subcluster()
        self.assertEqual(self.gnm.cluster_subcluster, "A")

    def test_set_cluster_subcluster_6(self):
        """Check that Cluster is set when subcluster is empty string."""
        self.gnm.subcluster = ""
        self.gnm.cluster = "A"
        self.gnm.set_cluster_subcluster()
        self.assertEqual(self.gnm.cluster_subcluster, "A")

    def test_set_cluster_subcluster_7(self):
        """Check that cluster_subcluster is set when provided value is None."""
        self.gnm.set_cluster_subcluster(None)
        self.assertEqual(self.gnm.cluster_subcluster, "singleton")

    def test_set_cluster_subcluster_8(self):
        """Check that cluster_subcluster is set when provided value is
        'Singleton'."""
        self.gnm.set_cluster_subcluster("Singleton")
        self.assertEqual(self.gnm.cluster_subcluster, "singleton")

    def test_set_cluster_subcluster_9(self):
        """Check that cluster_subcluster is set when provided value is
        'A2'."""
        self.gnm.set_cluster_subcluster("A2")
        self.assertEqual(self.gnm.cluster_subcluster, "A2")






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
        self.gnm.set_date(date1, \
                                            "empty_datetime_obj")
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



    # TODO this is probably no longer needed.
    # def test_set_annotation_author_1(self):
    #     """Check that annotation_author to set to 1 when author
    #     is in author_set."""
    #     input_author = "Hatfull"
    #     self.gnm.set_annotation_author(input_author)
    #     self.assertEqual(self.gnm.annotation_author, 1)
    #
    # def test_set_annotation_author_2(self):
    #     """Check that annotation_author to set to 0 when author
    #     is not in author_set."""
    #     input_author = "Unknown author"
    #     self.gnm.set_annotation_author(input_author)
    #     self.assertEqual(self.gnm.annotation_author, 0)




    def test_tally_descriptions_1(self):
        """Check that no description tally is incremented."""
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.tally_descriptions()

        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_products_tally, 0)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_functions_tally, 0)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_notes_tally, 0)

    def test_tally_descriptions_2(self):
        """Check that processed primary description tally
        is incremented."""
        self.cds1.processed_description = "abcd"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.tally_descriptions()

        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_descriptions_tally, 1)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_products_tally, 0)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_functions_tally, 0)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_notes_tally, 0)

    def test_tally_descriptions_3(self):
        """Check that processed product description tally
        is incremented."""
        self.cds1.processed_product = "abcd"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.tally_descriptions()

        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_products_tally, 1)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_functions_tally, 0)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_notes_tally, 0)

    def test_tally_descriptions_4(self):
        """Check that processed function description tally
        is incremented."""
        self.cds1.processed_function = "abcd"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.tally_descriptions()

        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_products_tally, 0)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_functions_tally, 1)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_notes_tally, 0)

    def test_tally_descriptions_5(self):
        """Check that processed note description tally
        is incremented."""
        self.cds1.processed_note = "abcd"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.tally_descriptions()

        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_products_tally, 0)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_functions_tally, 0)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_notes_tally, 1)

    def test_tally_descriptions_6(self):
        """Check that all description tallies are incremented."""
        self.cds1.processed_description = "abcd"
        self.cds1.processed_product = "efgh"
        self.cds1.processed_function = "ijkl"
        self.cds1.processed_note = "mnop"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.tally_descriptions()

        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_descriptions_tally, 1)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_products_tally, 1)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_functions_tally, 1)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_notes_tally, 1)

    def test_tally_descriptions_7(self):
        """Check that all description tallies are incremented for
        all CDS features."""
        self.cds1.processed_description = "abcd"
        self.cds1.processed_product = "efgh"
        self.cds1.processed_function = "ijkl"
        self.cds1.processed_note = "mnop"

        self.cds2.processed_description = "ab"
        self.cds2.processed_product = "cd"
        self.cds2.processed_function = "ef"
        self.cds2.processed_note = "gh"

        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.tally_descriptions()

        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_descriptions_tally, 2)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_products_tally, 2)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_functions_tally, 2)
        with self.subTest():
            self.assertEqual(
                self.gnm._cds_processed_notes_tally, 2)




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
        self.assertEqual(self.gnm.evaluations[0].status, "not_evaluated")




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
        self.assertEqual(self.gnm.evaluations[0].status, "not_evaluated")




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




    def test_check_compatible_status_and_accession_1(self):
        """Check final annotation_status with accession."""
        self.gnm.annotation_status = "final"
        self.gnm.accession = "ABC123"
        self.gnm.check_compatible_status_and_accession("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_compatible_status_and_accession_2(self):
        """Check final annotation_status with no accession."""
        self.gnm.annotation_status = "final"
        self.gnm.accession = ""
        self.gnm.check_compatible_status_and_accession()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_compatible_status_and_accession_3(self):
        """Check draft annotation_status with no accession."""
        self.gnm.annotation_status = "draft"
        self.gnm.accession = ""
        self.gnm.check_compatible_status_and_accession()
        self.assertEqual(self.gnm.evaluations[0].status, "correct")






    def test_check_compatible_status_and_descriptions_1(self):
        """Check that draft genome with no descriptions does not produce
        an error."""
        self.gnm.annotation_status = "draft"
        self.gnm.check_compatible_status_and_descriptions("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_compatible_status_and_descriptions_2(self):
        """Check that draft genome with a description produces an error."""
        self.gnm.annotation_status = "draft"
        self.gnm._cds_processed_descriptions_tally = 1
        self.gnm.check_compatible_status_and_descriptions()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_compatible_status_and_descriptions_3(self):
        """Check that final genome with a description does not produce
        an error."""
        self.gnm.annotation_status = "final"
        self.gnm._cds_processed_descriptions_tally = 1
        self.gnm.check_compatible_status_and_descriptions()
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_compatible_status_and_descriptions_4(self):
        """Check that final genome with no descriptions produces an error."""
        self.gnm.annotation_status = "final"
        self.gnm.check_compatible_status_and_descriptions()
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_compatible_status_and_descriptions_5(self):
        """Check that gbk genome with no descriptions does not produce
        an error."""
        self.gnm.annotation_status = "gbk"
        self.gnm.check_compatible_status_and_descriptions()
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_compatible_status_and_descriptions_6(self):
        """Check that gbk genome with descriptions does not produce
        an error."""
        self.gnm.annotation_status = "gbk"
        self.gnm._cds_processed_descriptions_tally = 1
        self.gnm.check_compatible_status_and_descriptions()
        self.assertEqual(self.gnm.evaluations[0].status, "correct")












    def test_check_description_name_1(self):
        """Check that no warning is produced."""
        self.gnm.id = "Trixie"
        self.gnm._description_name = "Trixie"
        self.gnm.check_description_name("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_description_name_2(self):
        """Check that a warning is produced."""
        self.gnm.id = "L5"
        self.gnm._description_name = "Trixie"
        self.gnm.check_description_name()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)




    def test_check_source_name_1(self):
        """Check that no warning is produced."""
        self.gnm.id = "Trixie"
        self.gnm._source_name = "Trixie"
        self.gnm.check_source_name("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_source_name_2(self):
        """Check that a warning is produced."""
        self.gnm.id = "L5"
        self.gnm._source_name = "Trixie"
        self.gnm.check_source_name()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)




    def test_check_organism_name_1(self):
        """Check that no warning is produced."""
        self.gnm.id = "Trixie"
        self.gnm._organism_name = "Trixie"
        self.gnm.check_organism_name("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_organism_name_2(self):
        """Check that a warning is produced."""
        self.gnm.id = "L5"
        self.gnm._organism_name = "Trixie"
        self.gnm.check_organism_name()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)




    def test_check_description_host_genus_1(self):
        """Check that no warning is produced."""
        self.gnm.host_genus = "Mycobacterium"
        self.gnm._description_host_genus = "Mycobacterium"
        self.gnm.check_description_host_genus("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_description_host_genus_2(self):
        """Check that a warning is produced."""
        self.gnm.host_genus = "Gordonia"
        self.gnm._description_host_genus = "Mycobacterium"
        self.gnm.check_description_host_genus()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)




    def test_check_source_host_genus_1(self):
        """Check that no warning is produced."""
        self.gnm.host_genus = "Mycobacterium"
        self.gnm._source_host_genus = "Mycobacterium"
        self.gnm.check_source_host_genus("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_source_host_genus_2(self):
        """Check that a warning is produced."""
        self.gnm.host_genus = "Gordonia"
        self.gnm._source_host_genus = "Mycobacterium"
        self.gnm.check_source_host_genus()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)




    def test_check_organism_host_genus_1(self):
        """Check that no warning is produced."""
        self.gnm.host_genus = "Mycobacterium"
        self.gnm._organism_host_genus = "Mycobacterium"
        self.gnm.check_organism_host_genus("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_organism_host_genus_2(self):
        """Check that a warning is produced."""
        self.gnm.host_genus = "Gordonia"
        self.gnm._organism_host_genus = "Mycobacterium"
        self.gnm.check_organism_host_genus()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)




    def test_check_authors_1(self):
        """Check that no warning is produced when author is expected
        and present."""
        check_set = set(["hatfull"])
        self.gnm.authors = "abcd; efgh; HATFULL; xyz"
        self.gnm.check_authors(check_set=check_set, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_authors_2(self):
        """Check that no warning is produced when author is not expected
        and not present."""
        check_set = set(["hatfull"])
        self.gnm.authors = "abcd; efgh; xyz"
        self.gnm.check_authors(check_set=check_set, expect=False)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_authors_3(self):
        """Check that a warning is produced when author is expected
        and not present."""
        check_set = set(["hatfull"])
        self.gnm.authors = "abcd; efgh; xyz"
        self.gnm.check_authors(check_set=check_set)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_authors_4(self):
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
        self.gnm.authors = "hatfull; 1234; xyz"
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




    def test_set_unique_cds_end_strand_ids_1(self):
        """Verify that both sets are computed."""
        self.gnm._cds_end_strand_ids = \
            [(1, "forward"), (2, "reverse"), (2, "forward"), (2, "reverse")]
        expected_unique_set = set([(1, "forward"), (2, "forward")])
        expected_duplicate_set = set([(2, "reverse")])
        self.gnm.set_unique_cds_end_strand_ids()
        with self.subTest():
            self.assertEqual(self.gnm._cds_unique_end_strand_ids, \
                expected_unique_set)
        with self.subTest():
            self.assertEqual(self.gnm._cds_duplicate_end_strand_ids, \
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




    def test_check_cds_end_strand_ids_1(self):
        """Verify that no warning is produced."""
        self.gnm.check_cds_end_strand_ids("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_cds_end_strand_ids_2(self):
        """Verify that a warning is produced."""
        self.gnm._cds_duplicate_end_strand_ids = set([(2, "forward")])
        self.gnm.check_cds_end_strand_ids()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)




    def test_set_value_flag_1(self):
        """Verify that the 'retain' setting is set to True."""
        self.gnm._value_flag = False
        self.gnm.cluster = "retain"
        self.gnm.set_value_flag("retain")
        self.assertTrue(self.gnm._value_flag)

    def test_set_value_flag_2(self):
        """Verify that the 'retain' setting is set to False."""
        self.gnm._value_flag = True
        self.gnm.cluster = "A"
        self.gnm.set_value_flag("retain")
        self.assertFalse(self.gnm._value_flag)




    def test_set_feature_ids_1(self):
        """Verify that CDS features are sorted in correct order."""

        #1
        self.trna3.left = 10
        self.trna3.right = 15

        #2
        self.trna2.left = 10
        self.trna2.right = 17

        #3
        self.cds4.left = 10
        self.cds4.right = 20

        #4
        self.cds3.left = 18
        self.cds3.right = 30

        #5
        self.cds2.left = 18
        self.cds2.right = 40

        #6
        self.trna1.left = 100
        self.trna1.right = 200

        #7 Wrap-around feature.
        self.cds1.left = 400
        self.cds1.right = 3

        self.gnm.id = "L5"
        self.gnm.cds_features = [self.cds1, self.cds2, self.cds3, self.cds4]
        self.gnm.trna_features = [self.trna1, self.trna2, self.trna3]
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

        #1
        self.trna3.left = 10
        self.trna3.right = 15

        #2
        self.trna2.left = 10
        self.trna2.right = 17

        #3
        self.cds4.left = 10
        self.cds4.right = 20

        #4
        self.cds3.left = 18
        self.cds3.right = 30

        #5
        self.cds2.left = 18
        self.cds2.right = 40

        #6
        self.trna1.left = 100
        self.trna1.right = 200

        #7 Wrap-around feature.
        self.cds1.left = 400
        self.cds1.right = 3

        self.gnm.id = "L5"
        self.gnm.cds_features = [self.cds1, self.cds2, self.cds3, self.cds4]
        self.gnm.trna_features = [self.trna1, self.trna2, self.trna3]
        self.gnm.set_feature_ids(use_trna=True)
        with self.subTest():
            self.assertEqual(self.trna1.id, "L5_3")
        with self.subTest():
            self.assertEqual(self.trna2.id, "L5_2")
        with self.subTest():
            self.assertEqual(self.trna3.id, "L5_1")


    def test_set_feature_ids_3(self):
        """Verify that CDS and tRNA features are sorted in correct order."""

        #1
        self.trna3.left = 10
        self.trna3.right = 15

        #2
        self.trna2.left = 10
        self.trna2.right = 17

        #3
        self.cds4.left = 10
        self.cds4.right = 20

        #4
        self.cds3.left = 18
        self.cds3.right = 30

        #5
        self.cds2.left = 18
        self.cds2.right = 40

        #6
        self.trna1.left = 100
        self.trna1.right = 200

        #7 Wrap-around feature.
        self.cds1.left = 400
        self.cds1.right = 3

        self.gnm.id = "L5"
        self.gnm.cds_features = [self.cds1, self.cds2, self.cds3, self.cds4]
        self.gnm.trna_features = [self.trna1, self.trna2, self.trna3]
        self.gnm.set_feature_ids(use_cds=True, use_trna=True)
        with self.subTest():
            self.assertEqual(self.trna1.id, "L5_6")
        with self.subTest():
            self.assertEqual(self.trna2.id, "L5_2")
        with self.subTest():
            self.assertEqual(self.trna3.id, "L5_1")
        with self.subTest():
            self.assertEqual(self.cds1.id, "L5_7")
        with self.subTest():
            self.assertEqual(self.cds2.id, "L5_5")
        with self.subTest():
            self.assertEqual(self.cds3.id, "L5_4")
        with self.subTest():
            self.assertEqual(self.cds4.id, "L5_3")


    def test_set_feature_ids_4(self):
        """Verify that CDS and tRNA features are sorted in correct order
        with type delimiter added."""

        #1
        self.trna3.left = 10
        self.trna3.right = 15

        #2
        self.cds4.left = 10
        self.cds4.right = 20

        self.gnm.id = "L5"
        self.gnm.cds_features = [self.cds4]
        self.gnm.trna_features = [self.trna3]
        self.gnm.set_feature_ids(use_type=True, use_cds=True, use_trna=True)
        with self.subTest():
            self.assertEqual(self.trna3.id, "L5_TRNA_1")
        with self.subTest():
            self.assertEqual(self.cds4.id, "L5_CDS_2")







    def test_check_id_1(self):
        """Verify that no error is produced when the id
        is in the id_set and is expected to be in the set."""
        value_set = set(["Trixie", "L5"])
        self.gnm.id = "Trixie"
        self.gnm.check_id(value_set, True, "eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")


    def test_check_id_2(self):
        """Verify that an error is produced when the id
        is not in the id_set and is expected to be in the set."""
        value_set = set(["Trixie", "L5"])
        self.gnm.id = "D29"
        self.gnm.check_id(value_set, True)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)



    def test_check_id_3(self):
        """Verify that no error is produced when the id
        is not in the id_set and is not expected to be in the set."""
        value_set = set(["Trixie", "L5"])
        self.gnm.id = "D29"
        self.gnm.check_id(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_id_4(self):
        """Verify that an error is produced when the id
        is in the id_set and is not expected to be in the set."""
        value_set = set(["Trixie", "L5"])
        self.gnm.id = "Trixie"
        self.gnm.check_id(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")





    def test_check_name_1(self):
        """Verify that no error is produced when the name
        is in the name_set and is expected to be in the set."""
        value_set = set(["Trixie", "L5"])
        self.gnm.name = "Trixie"
        self.gnm.check_name(value_set, True, "eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_name_2(self):
        """Verify that an error is produced when the name
        is not in the name_set and is expected to be in the set."""
        value_set = set(["Trixie", "L5"])
        self.gnm.name = "D29"
        self.gnm.check_name(value_set, True)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_name_3(self):
        """Verify that no error is produced when the name
        is not in the name_set and is not expected to be in the set."""
        value_set = set(["Trixie", "L5"])
        self.gnm.name = "D29"
        self.gnm.check_name(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_name_4(self):
        """Verify that an error is produced when the name
        is in the name_set and is not expected to be in the set."""
        value_set = set(["Trixie", "L5"])
        self.gnm.name = "Trixie"
        self.gnm.check_name(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")




    def test_check_annotation_status_1(self):
        """Verify that no error is produced when the annotation_status
        is in the status_set and is expected to be in the set."""
        self.gnm.annotation_status = "draft"
        self.gnm.check_annotation_status(
            check_set=constants.ANNOTATION_STATUS_SET,
            expect=True, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_annotation_status_2(self):
        """Verify that an error is produced when the annotation_status
        is not in the status_set and is expected to be in the set."""
        self.gnm.annotation_status = "invalid"
        self.gnm.check_annotation_status(
            check_set=constants.ANNOTATION_STATUS_SET, expect=True)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_annotation_status_3(self):
        """Verify that no error is produced when the annotation_status
        is in a non-standard status_set and is expected to be in the set."""
        check_set = set(["new_status", "final"])
        self.gnm.annotation_status = "new_status"
        self.gnm.check_annotation_status(
            check_set=check_set, expect = True)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_annotation_status_4(self):
        """Verify that an error is produced when the annotation_status
        is not a non-standard status_set and is expected to be in the set."""
        check_set = set(["new_status", "final"])
        self.gnm.annotation_status = "draft"
        self.gnm.check_annotation_status(
            check_set=check_set, expect = True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_annotation_status_5(self):
        """Verify that an error is produced when the annotation_status
        is in the set, and is not expected to be in the set."""
        check_set = set(["new_status", "final"])
        self.gnm.annotation_status = "new_status"
        self.gnm.check_annotation_status(
            check_set=check_set, expect = False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")




    def test_check_host_genus_1(self):
        """Verify that no error is produced when the host_genus
        is in the host_set and is expected to be in the set."""
        value_set = set(["Mycobacterium", "Gordonia"])
        self.gnm.host_genus = "Mycobacterium"
        self.gnm.check_host_genus(value_set, True, "eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_host_genus_2(self):
        """Verify that an error is produced when the host_genus
        is not in the host_set and is expected to be in the set."""
        value_set = set(["Mycobacterium", "Gordonia"])
        self.gnm.host_genus = "invalid"
        self.gnm.check_host_genus(value_set, True)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_host_genus_3(self):
        """Verify that no error is produced when the host_genus
        is not in the host_set and is not expected to be in the set."""
        value_set = set(["Mycobacterium", "Gordonia"])
        self.gnm.host_genus = "Arthrobacter"
        self.gnm.check_host_genus(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_host_genus_4(self):
        """Verify that an error is produced when the host_genus
        is in the host_set and is not expected to be in the set."""
        value_set = set(["Mycobacterium", "Gordonia"])
        self.gnm.host_genus = "Mycobacterium"
        self.gnm.check_host_genus(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")




    def test_check_cluster_1(self):
        """Verify that no error is produced when the cluster
        is in the cluster_set and is expected to be in the set."""
        value_set = set(["A", "B"])
        self.gnm.cluster = "A"
        self.gnm.check_cluster(value_set, True, "eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_cluster_2(self):
        """Verify that an error is produced when the cluster
        is not in the cluster_set and is expected to be in the set."""
        value_set = set(["A", "B"])
        self.gnm.cluster = "C"
        self.gnm.check_cluster(value_set, True)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_cluster_3(self):
        """Verify that no error is produced when the cluster
        is not in the cluster_set and is not expected to be in the set."""
        value_set = set(["A", "B"])
        self.gnm.cluster = "C"
        self.gnm.check_cluster(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_cluster_4(self):
        """Verify that an error is produced when the cluster
        is in the cluster_set and is not expected to be in the set."""
        value_set = set(["A", "B"])
        self.gnm.cluster = "A"
        self.gnm.check_cluster(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")




    def test_check_subcluster_1(self):
        """Verify that no error is produced when the subcluster
        is in the subcluster_set and is expected to be in the set."""
        value_set = set(["A1", "B1"])
        self.gnm.subcluster = "A1"
        self.gnm.check_subcluster(value_set, True, "eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_subcluster_2(self):
        """Verify that an error is produced when the subcluster
        is not in the subcluster_set and is expected to be in the set."""
        value_set = set(["A1", "B1"])
        self.gnm.subcluster = "C1"
        self.gnm.check_subcluster(value_set, True)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_subcluster_3(self):
        """Verify that no error is produced when the subcluster
        is not in the subcluster_set and is not expected to be in the set."""
        value_set = set(["A1", "B1"])
        self.gnm.subcluster = "A2"
        self.gnm.check_subcluster(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_subcluster_4(self):
        """Verify that an error is produced when the subcluster
        is in the subcluster_set and is not expected to be in the set."""
        value_set = set(["A1", "B1"])
        self.gnm.subcluster = "A1"
        self.gnm.check_subcluster(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")




    def test_check_sequence_1(self):
        """Verify that no error is produced when the sequence
        is in the seq_set and is expected to be in the set."""
        value_set = set([Seq("ATCG"), Seq("AACCGGTT")])
        self.gnm.seq = Seq("ATCG")
        self.gnm.check_sequence(value_set, True, "eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_sequence_2(self):
        """Verify that an error is produced when the sequence
        is not in the seq_set and is expected to be in the set."""
        value_set = set([Seq("ATCG"), Seq("AACCGGTT")])
        self.gnm.seq = Seq("TTTTT")
        self.gnm.check_sequence(value_set, True)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_sequence_3(self):
        """Verify that no error is produced when the sequence
        is not in the seq_set and is not expected to be in the set."""
        value_set = set([Seq("ATCG"), Seq("AACCGGTT")])
        self.gnm.seq = Seq("TTTTT")
        self.gnm.check_sequence(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_sequence_4(self):
        """Verify that an error is produced when the sequence
        is in the seq_set and is not expected to be in the set."""
        value_set = set([Seq("ATCG"), Seq("AACCGGTT")])
        self.gnm.seq = Seq("ATCG")
        self.gnm.check_sequence(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")




    def test_check_accession_1(self):
        """Verify that no error is produced when the accession
        is in the accession_set and is expected to be in the set."""
        value_set = set(["ABC123", "XYZ456"])
        self.gnm.accession = "ABC123"
        self.gnm.check_accession(value_set, True, "eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_accession_2(self):
        """Verify that an error is produced when the accession
        is not in the accession_set and is expected to be in the set."""
        value_set = set(["ABC123", "XYZ456"])
        self.gnm.accession = "EFG789"
        self.gnm.check_accession(value_set, True)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_accession_3(self):
        """Verify that no error is produced when the accession
        is not in the accession_set and is not expected to be in the set."""
        value_set = set(["ABC123", "XYZ456"])
        self.gnm.accession = "EFG789"
        self.gnm.check_accession(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_accession_4(self):
        """Verify that an error is produced when the accession
        is in the accession_set and is not expected to be in the set."""
        value_set = set(["ABC123", "XYZ456"])
        self.gnm.accession = "ABC123"
        self.gnm.check_accession(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")




    def test_check_annotation_author_1(self):
        """Verify that no error is produced when the annotation_author
        is valid."""
        self.gnm.annotation_author = 0
        self.gnm.check_annotation_author(
            check_set=constants.ANNOTATION_AUTHOR_SET, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_annotation_author_2(self):
        """Verify that no error is produced when the annotation_author
        is valid."""
        self.gnm.annotation_author = 1
        self.gnm.check_annotation_author(
            check_set=constants.ANNOTATION_AUTHOR_SET)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_annotation_author_3(self):
        """Verify that an error is produced when the annotation_author
        is not valid."""
        self.gnm.annotation_author = "1"
        self.gnm.check_annotation_author(
            check_set=constants.ANNOTATION_AUTHOR_SET)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_annotation_author_4(self):
        """Verify that no error is produced when the annotation_author
        is valid based on a supplied set."""
        check_set = set(["1", "2"])
        self.gnm.annotation_author = "1"
        self.gnm.check_annotation_author(check_set=check_set)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_annotation_author_5(self):
        """Verify that an error is produced when the annotation_author
        is not valid based on a supplied set."""
        check_set = set(["1", "2"])
        self.gnm.annotation_author = 1
        self.gnm.check_annotation_author(check_set=check_set)
        self.assertEqual(self.gnm.evaluations[0].status, "error")




    def test_check_annotation_qc_1(self):
        """Verify that no error is produced when the annotation_qc
        is valid."""
        self.gnm.annotation_qc = 0
        self.gnm.check_annotation_qc(
            check_set=constants.ANNOTATION_QC_SET, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_annotation_qc_2(self):
        """Verify that no error is produced when the annotation_qc
        is valid."""
        self.gnm.annotation_qc = 1
        self.gnm.check_annotation_qc(check_set=constants.ANNOTATION_QC_SET)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_annotation_qc_3(self):
        """Verify that an error is produced when the annotation_qc
        is not valid."""
        self.gnm.annotation_qc = "1"
        self.gnm.check_annotation_qc(check_set=constants.ANNOTATION_QC_SET)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_annotation_qc_4(self):
        """Verify that no error is produced when the annotation_qc
        is valid based on a supplied set."""
        check_set = set(["1", "2"])
        self.gnm.annotation_qc = "1"
        self.gnm.check_annotation_qc(check_set=check_set)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_annotation_qc_5(self):
        """Verify that an error is produced when the annotation_qc
        is invalid based on a supplied set."""
        check_set = set(["1", "2"])
        self.gnm.annotation_qc = 1
        self.gnm.check_annotation_qc(check_set=check_set)
        self.assertEqual(self.gnm.evaluations[0].status, "error")




    def test_check_retrieve_record_1(self):
        """Verify that no error is produced when the retrieve_record
        is valid."""
        self.gnm.retrieve_record = 0
        self.gnm.check_retrieve_record(
            check_set=constants.RETRIEVE_RECORD_SET, eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_retrieve_record_2(self):
        """Verify that no error is produced when the retrieve_record
        is valid."""
        self.gnm.retrieve_record = 1
        self.gnm.check_retrieve_record(
            check_set=constants.RETRIEVE_RECORD_SET)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_retrieve_record_3(self):
        """Verify that an error is produced when the retrieve_record
        is not valid."""
        self.gnm.retrieve_record = "1"
        self.gnm.check_retrieve_record(
            check_set=constants.RETRIEVE_RECORD_SET)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_retrieve_record_4(self):
        """Verify that no error is produced when the retrieve_record
        is valid based on a supplied set."""
        check_set = set(["1", "2"])
        self.gnm.retrieve_record = "1"
        self.gnm.check_retrieve_record(check_set=check_set)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_retrieve_record_5(self):
        """Verify that no error is produced when the retrieve_record
        is valid based on a supplied set."""
        check_set = set(["1", "2"])
        self.gnm.retrieve_record = 1
        self.gnm.check_retrieve_record(check_set=check_set)
        self.assertEqual(self.gnm.evaluations[0].status, "error")




    def test_check_filename_1(self):
        """Verify that no error is produced when the filename
        is in the filename_set and is expected to be in the set."""
        value_set = set(["Trixie"])
        self.gnm.filename = "Trixie"
        self.gnm.check_filename(value_set, True, "eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_filename_2(self):
        """Verify that an error is produced when the filename
        is not in the filename_set and is expected to be in the set."""
        value_set = set(["Trixie"])
        self.gnm.filename = "L5"
        self.gnm.check_filename(value_set, True)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_filename_3(self):
        """Verify that no error is produced when the filename
        is not in the filename_set and is not expected to be in the set."""
        value_set = set(["Trixie"])
        self.gnm.filename = "L5"
        self.gnm.check_filename(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_filename_4(self):
        """Verify that an error is produced when the filename
        is in the filename_set and is not expected to be in the set."""
        value_set = set(["Trixie"])
        self.gnm.filename = "Trixie"
        self.gnm.check_filename(value_set, False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")




    def test_parse_description_1(self):
        """Verify empty string is parsed correctly."""
        self.gnm.description = ""
        self.gnm.parse_description()
        expected_phage = ""
        expected_host = ""
        with self.subTest():
            self.assertEqual(\
                self.gnm._description_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.gnm._description_host_genus, expected_host)

    def test_parse_description_2(self):
        """Verify string is parsed correctly."""
        self.gnm.description = "asdf Mycobacterium phage Trixie."
        self.gnm.parse_description()
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        with self.subTest():
            self.assertEqual(\
                self.gnm._description_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.gnm._description_host_genus, expected_host)




    def test_parse_source_1(self):
        """Verify empty string is parsed correctly."""
        self.gnm.source = ""
        self.gnm.parse_source()
        expected_phage = ""
        expected_host = ""
        with self.subTest():
            self.assertEqual(\
                self.gnm._source_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.gnm._source_host_genus, expected_host)

    def test_parse_source_2(self):
        """Verify string is parsed correctly."""
        self.gnm.source = "asdf Mycobacterium phage Trixie."
        self.gnm.parse_source()
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        with self.subTest():
            self.assertEqual(\
                self.gnm._source_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.gnm._source_host_genus, expected_host)




    def test_parse_organism_1(self):
        """Verify empty string is parsed correctly."""
        self.gnm.organism = ""
        self.gnm.parse_organism()
        expected_phage = ""
        expected_host = ""
        with self.subTest():
            self.assertEqual(\
                self.gnm._organism_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.gnm._organism_host_genus, expected_host)

    def test_parse_organism_2(self):
        """Verify string is parsed correctly."""
        self.gnm.organism = "asdf Mycobacterium phage Trixie."
        self.gnm.parse_organism()
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        with self.subTest():
            self.assertEqual(\
                self.gnm._organism_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.gnm._organism_host_genus, expected_host)




    def test_split_cluster_subcluster_1(self):
        """Verify split to only cluster, and subcluster is 'none'."""
        self.gnm.cluster_subcluster = "A"
        self.gnm.split_cluster_subcluster()
        cluster = "A"
        subcluster = "none"
        with self.subTest():
            self.assertEqual(self.gnm.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, subcluster)

    def test_split_cluster_subcluster_2(self):
        """Verify split to only cluster, and subcluster is None."""
        self.gnm.cluster_subcluster = "A"
        self.gnm.split_cluster_subcluster("none_object")
        cluster = "A"
        subcluster = None
        with self.subTest():
            self.assertEqual(self.gnm.cluster, cluster)
        with self.subTest():
            self.assertIsNone(self.gnm.subcluster)

    def test_split_cluster_subcluster_3(self):
        """Verify split to only cluster, and subcluster is ''."""
        self.gnm.cluster_subcluster = "A"
        self.gnm.split_cluster_subcluster("empty_string")
        cluster = "A"
        subcluster = ""
        with self.subTest():
            self.assertEqual(self.gnm.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, subcluster)

    def test_split_cluster_subcluster_4(self):
        """Verify singleton split to only cluster."""
        self.gnm.cluster_subcluster = "singleton"
        self.gnm.split_cluster_subcluster("empty_string")
        cluster = "singleton"
        subcluster = ""
        with self.subTest():
            self.assertEqual(self.gnm.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, subcluster)

    def test_split_cluster_subcluster_5(self):
        """Verify split to both cluster and subcluster."""
        self.gnm.cluster_subcluster = "A15"
        self.gnm.split_cluster_subcluster("empty_string")
        cluster = "A"
        subcluster = "A15"
        with self.subTest():
            self.assertEqual(self.gnm.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, subcluster)

    def test_split_cluster_subcluster_6(self):
        """Verify no split, and output is 'none'."""
        self.gnm.cluster_subcluster = "A1B2"
        self.gnm.split_cluster_subcluster()
        cluster = "none"
        subcluster = "none"
        with self.subTest():
            self.assertEqual(self.gnm.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, subcluster)

    def test_split_cluster_subcluster_7(self):
        """Verify no split, and output is ''."""
        self.gnm.cluster_subcluster = "A1B2"
        self.gnm.cluster = ""
        self.gnm.subcluster = ""
        self.gnm.split_cluster_subcluster("empty_string")
        cluster = ""
        subcluster = ""
        with self.subTest():
            self.assertEqual(self.gnm.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, subcluster)

    def test_split_cluster_subcluster_8(self):
        """Verify no split, and output is None."""
        self.gnm.cluster_subcluster = "A1B2"
        self.gnm.cluster = ""
        self.gnm.subcluster = ""
        self.gnm.split_cluster_subcluster("none_object")
        cluster = None
        subcluster = None
        with self.subTest():
            self.assertEqual(self.gnm.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, subcluster)

    def test_split_cluster_subcluster_9(self):
        """Verify no change to cluster and subcluster."""
        self.gnm.cluster_subcluster = ""
        self.gnm.cluster = "B"
        self.gnm.subcluster = "B10"
        self.gnm.split_cluster_subcluster("none_object")
        cluster = "B"
        subcluster = "B10"
        with self.subTest():
            self.assertEqual(self.gnm.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, subcluster)

    def test_split_cluster_subcluster_10(self):
        """Verify no change to cluster and subcluster."""
        self.gnm.cluster_subcluster = None
        self.gnm.cluster = "B"
        self.gnm.subcluster = "B10"
        self.gnm.split_cluster_subcluster("none_object")
        cluster = "B"
        subcluster = "B10"
        with self.subTest():
            self.assertEqual(self.gnm.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.gnm.subcluster, subcluster)




    def test_check_cds_feature_tally_1(self):
        """Verify no error is encountered when there is one CDS feature."""
        self.gnm._cds_features_tally = 1
        self.gnm.check_cds_feature_tally("eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_cds_feature_tally_2(self):
        """Verify error is encountered when there is no CDS feature."""
        self.gnm._cds_features_tally = 0
        self.gnm.check_cds_feature_tally()
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)




    def test_check_value_flag_1(self):
        """Verify that no error is produced when the _value_flag
        field is True and is expected to be True."""
        self.gnm._value_flag = True
        self.gnm.check_value_flag(True, "eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_value_flag_2(self):
        """Verify that an error is produced when the _value_flag
        field is False and is expected to be True."""
        self.gnm._value_flag = False
        self.gnm.check_value_flag(True)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_value_flag_3(self):
        """Verify that no error is produced when the _value_flag
        field is False and is expected to be False."""
        self.gnm._value_flag = False
        self.gnm.check_value_flag(False)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_value_flag_4(self):
        """Verify that an error is produced when the _value_flag
        field is True and is expected to be False."""
        self.gnm._value_flag = True
        self.gnm.check_value_flag(False)
        self.assertEqual(self.gnm.evaluations[0].status, "error")




    def test_check_feature_ids_1(self):
        """Verify no error is produced by two CDS features on same strand."""
        self.cds1.strand = "F"
        self.cds1.left = 5
        self.cds1.right = 50
        self.cds2.strand = "F"
        self.cds2.left = 20
        self.cds2.right = 70
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_ids(cds_ftr=True,eval_id="eval_id")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "correct")
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].id, "eval_id")

    def test_check_feature_ids_2(self):
        """Verify an error is produced by two CDS features on same strand
        with identical left and right coordinates."""
        self.cds1.strand = "F"
        self.cds1.left = 5
        self.cds1.right = 50
        self.cds2.strand = "F"
        self.cds2.left = 5
        self.cds2.right = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_ids(cds_ftr=True)
        with self.subTest():
            self.assertEqual(self.gnm.evaluations[0].status, "error")
        with self.subTest():
            self.assertIsNone(self.gnm.evaluations[0].id)

    def test_check_feature_ids_3(self):
        """Verify an error is produced by two CDS features on different strand
        with identical left and right coordinates."""
        self.cds1.strand = "R"
        self.cds1.left = 5
        self.cds1.right = 50
        self.cds2.strand = "F"
        self.cds2.left = 5
        self.cds2.right = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_ids(cds_ftr=True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_feature_ids_4(self):
        """Verify no error is produced by two CDS features on different strand
        with identical left and right coordinates when strand is True."""
        self.cds1.strand = "R"
        self.cds1.left = 5
        self.cds1.right = 50
        self.cds2.strand = "F"
        self.cds2.left = 5
        self.cds2.right = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_ids(cds_ftr=True, strand=True)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_feature_ids_5(self):
        """Verify an error is produced by two CDS features on same strand
        when they have nested coordinates."""
        self.cds1.strand = "F"
        self.cds1.left = 10
        self.cds1.right = 20
        self.cds2.strand = "F"
        self.cds2.left = 5
        self.cds2.right = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_ids(cds_ftr=True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_feature_ids_6(self):
        """Verify no error is produced by two CDS features on "F" strand
        when they have the same left (start) coordinates."""
        self.cds1.strand = "F"
        self.cds1.left = 10
        self.cds1.right = 20
        self.cds2.strand = "F"
        self.cds2.left = 10
        self.cds2.right = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_ids(cds_ftr=True)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_feature_ids_7(self):
        """Verify an error is produced by two CDS features on "F" strand
        when they have the same right (stop) coordinates."""
        self.cds1.strand = "F"
        self.cds1.left = 10
        self.cds1.right = 50
        self.cds2.strand = "F"
        self.cds2.left = 5
        self.cds2.right = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_ids(cds_ftr=True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_feature_ids_8(self):
        """Verify no error is produced by two CDS features on "R" strand
        when they have the same right (start) coordinates."""
        self.cds1.strand = "R"
        self.cds1.left = 10
        self.cds1.right = 50
        self.cds2.strand = "R"
        self.cds2.left = 5
        self.cds2.right = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_ids(cds_ftr=True)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_feature_ids_9(self):
        """Verify no error is produced by two CDS features on "R" strand
        when they have the same left (stop) coordinates."""
        self.cds1.strand = "R"
        self.cds1.left = 5
        self.cds1.right = 20
        self.cds2.strand = "R"
        self.cds2.left = 5
        self.cds2.right = 50
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_ids(cds_ftr=True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_feature_ids_10(self):
        """Verify an error is produced by a CDS and tRNA feature
        when they have the same coordinates."""
        self.cds1.strand = "F"
        self.cds1.left = 5
        self.cds1.right = 20
        self.cds2.strand = "F"
        self.cds2.left = 100
        self.cds2.right = 200
        self.trna1.left = 5
        self.trna1.right = 20
        self.trna1.strand = "F"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.trna_features = [self.trna1]
        self.gnm.check_feature_ids(cds_ftr=True, trna_ftr=True)
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_feature_ids_11(self):
        """Verify no error is produced by a CDS and tRNA feature
        when they have the same coordinates when CDS is false."""
        self.cds1.strand = "F"
        self.cds1.left = 5
        self.cds1.right = 20
        self.cds2.strand = "F"
        self.cds2.left = 100
        self.cds2.right = 200
        self.trna1.left = 5
        self.trna1.right = 20
        self.trna1.strand = "F"
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.trna_features = [self.trna1]
        self.gnm.check_feature_ids(trna_ftr=True)
        self.assertEqual(self.gnm.evaluations[0].status, "correct")

    def test_check_feature_ids_12(self):
        """Verify an error is produced by a CDS feature stored in the
        CDS features list and a CDS feature stored in a supplementary list
        when they have the same coordinates."""
        self.cds1.strand = "F"
        self.cds1.left = 5
        self.cds1.right = 20
        self.cds2.strand = "F"
        self.cds2.left = 100
        self.cds2.right = 200
        self.cds3.strand = "F"
        self.cds3.left = 5
        self.cds3.right = 20
        self.gnm.cds_features = [self.cds1, self.cds2]
        self.gnm.check_feature_ids(cds_ftr=True, other=[self.cds3])
        self.assertEqual(self.gnm.evaluations[0].status, "error")

    def test_check_feature_ids_13(self):
        """Verify no error is produced when empty lists are passed through."""
        self.gnm.check_feature_ids(cds_ftr=True, trna_ftr=True, tmrna=True, other=[])
        self.assertEqual(self.gnm.evaluations[0].status, "correct")








if __name__ == '__main__':
    unittest.main()
