""" Unit tests for the Genome class."""


import unittest
from classes import Genome
from classes import Cds
from datetime import datetime



class TestGenomeClass(unittest.TestCase):


    def setUp(self):
        self.genome = Genome.Genome()



        self.cds1 = Cds.CdsFeature()
        self.cds1.processed_primary_description = ""
        self.cds1.processed_product_description = ""
        self.cds1.processed_function_description = ""
        self.cds1.processed_note_description = ""

        self.cds2 = Cds.CdsFeature()
        self.cds2.processed_primary_description = ""
        self.cds2.processed_product_description = ""
        self.cds2.processed_function_description = ""
        self.cds2.processed_note_description = ""





    def test_set_filename_1(self):
        """Confirm file path is split appropriately."""
        filepath = "/path/to/folder/Trixie.gbk"
        self.genome.set_filename(filepath)
        with self.subTest():
            self.assertEqual(self.genome.filename, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome.search_filename, "trixie")




    def test_set_host_1(self):
        """Check that host name is split appropriately."""
        host = "Mycobacterium smegmatis"
        self.genome.set_host(host, "none_string")
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_2(self):
        """Check that whitespace is removed."""
        host = "  Mycobacterium smegmatis  "
        self.genome.set_host(host, "none_string")
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_3(self):
        """Check that none is set appropriately."""
        host = ""
        self.genome.set_host(host, "none_string")
        self.assertEqual(self.genome.host, "none")

    def test_set_host_4(self):
        """Check that None object is set appropriately."""
        host = ""
        self.genome.set_host(host, "none_object")
        self.assertIsNone(self.genome.host)





    def test_set_sequence_1(self):
        """Check that sequence is set appropriately."""
        seq = "aaggcga"
        self.genome.set_sequence(seq)
        with self.subTest():
            self.assertEqual(self.genome.sequence, "AAGGCGA")
        with self.subTest():
            self.assertEqual(self.genome._length, 7)
        with self.subTest():
            self.assertEqual(self.genome._gc, 57.1429)

    def test_set_sequence_2(self):
        """Check that sequence is set appropriately if it is an empty string."""
        seq = ""
        self.genome.set_sequence(seq)
        with self.subTest():
            self.assertEqual(self.genome.sequence, "")
        with self.subTest():
            self.assertEqual(self.genome._length, 0)
        with self.subTest():
            self.assertEqual(self.genome._gc, -1)






    def test_set_accession_1(self):
        """Check that accession is set appropriately."""
        accession = "ABC123.1"
        self.genome.set_accession(accession, "empty_string")
        self.assertEqual(self.genome.accession, "ABC123")

    def test_set_accession_2(self):
        """Check that an empty string accession is set appropriately
        from None."""
        accession = None
        self.genome.set_accession(accession, "empty_string")
        self.assertEqual(self.genome.accession, "")

    def test_set_accession_3(self):
        """Check that 'none' string accession is set appropriately
        from None."""
        accession = None
        self.genome.set_accession(accession, "none_string")
        self.assertEqual(self.genome.accession, "none")

    def test_set_accession_4(self):
        """Check that a None accession is set appropriately from None."""
        accession = None
        self.genome.set_accession(accession, "none_object")
        self.assertIsNone(self.genome.accession)

    def test_set_accession_5(self):
        """Check that a None accession is set appropriately from an
        empty string."""
        accession = ""
        self.genome.set_accession(accession, "none_object")
        self.assertIsNone(self.genome.accession)




    def test_set_cds_features_1(self):
        """Check that CDS feature list is set and length is computed."""
        features_list = [0, 1, 2, 3]
        self.genome.set_cds_features(features_list)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 4)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 4)




    def test_set_cds_ids_1(self):
        """Check that CDS feature identifier lists are set."""
        cds1 = Cds.CdsFeature()
        cds1._start_end_id = (1, 5)
        cds1._end_strand_id = (5, "forward")

        cds2 = Cds.CdsFeature()
        cds2._start_end_id = (21, 2)
        cds2._end_strand_id = (2, "reverse")

        features_list = [cds1, cds2]
        self.genome.set_cds_features(features_list)
        self.genome.set_cds_ids()

        start_end_id_list = [(1,5), (21,2)]
        end_strand_id_list = [(5, "forward"), (2, "reverse")]

        with self.subTest():
            self.assertEqual(self.genome._cds_start_end_ids, \
            start_end_id_list)
        with self.subTest():
            self.assertEqual(self.genome._cds_end_strand_ids, \
            end_strand_id_list)




    def test_set_trna_features_1(self):
        """Check that tRNA feature list is set and length is computed."""
        features_list = [0,1,2,3]
        self.genome.set_trna_features(features_list)
        with self.subTest():
            self.assertEqual(len(self.genome.trna_features), 4)
        with self.subTest():
            self.assertEqual(self.genome._trna_features_tally, 4)




    def test_set_source_features_1(self):
        """Check that source feature list is set and length is computed."""
        features_list = [0,1,2,3]
        self.genome.set_source_features(features_list)
        with self.subTest():
            self.assertEqual(len(self.genome.source_features), 4)
        with self.subTest():
            self.assertEqual(self.genome._source_features_tally, 4)




    def test_set_phage_id_1(self):
        """Check that name without '_Draft' suffix is not changed."""
        phage_name = "Trixie"
        self.genome.set_phage_id(phage_name)
        with self.subTest():
            self.assertEqual(self.genome.phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome.search_id, "trixie")

    def test_set_phage_id_2(self):
        """Check that '_Draft' suffix is removed."""
        phage_name = "Trixie_Draft"
        self.genome.set_phage_id(phage_name)
        with self.subTest():
            self.assertEqual(self.genome.phage_id, "Trixie")
        with self.subTest():
            self.assertEqual(self.genome.search_id, "trixie")




    def test_set_phage_id_from_field_1(self):
        """Check that the phage_id is set from the phage_name field."""
        self.genome.phage_name = "Trixie_Draft"
        self.genome.set_phage_id_from_field("phage_name")
        self.assertEqual(self.genome.phage_id, "Trixie")

    def test_set_phage_id_from_field_2(self):
        """Check that the phage_id is set from the accession field."""
        self.genome.accession = "Trixie_Draft"
        self.genome.set_phage_id_from_field("accession")
        self.assertEqual(self.genome.phage_id, "Trixie")

    def test_set_phage_id_from_field_3(self):
        """Check that the phage_id is set from the record_name field."""
        self.genome.record_name = "Trixie_Draft"
        self.genome.set_phage_id_from_field("record_name")
        self.assertEqual(self.genome.phage_id, "Trixie")

    def test_set_phage_id_from_field_4(self):
        """Check that the phage_id is set from the record_id field."""
        self.genome.record_id = "Trixie_Draft"
        self.genome.set_phage_id_from_field("record_id")
        self.assertEqual(self.genome.phage_id, "Trixie")

    def test_set_phage_id_from_field_5(self):
        """Check that the phage_id is set from the record_accession field."""
        self.genome.record_accession = "Trixie_Draft"
        self.genome.set_phage_id_from_field("record_accession")
        self.assertEqual(self.genome.phage_id, "Trixie")

    def test_set_phage_id_from_field_6(self):
        """Check that the phage_id is set from the record_description field."""
        self.genome.record_description = "Trixie_Draft"
        self.genome.set_phage_id_from_field("record_description")
        self.assertEqual(self.genome.phage_id, "Trixie")

    def test_set_phage_id_from_field_7(self):
        """Check that the phage_id is set from the record_source field."""
        self.genome.record_source = "Trixie_Draft"
        self.genome.set_phage_id_from_field("record_source")
        self.assertEqual(self.genome.phage_id, "Trixie")

    def test_set_phage_id_from_field_8(self):
        """Check that the phage_id is set from the record_organism field."""
        self.genome.record_organism = "Trixie_Draft"
        self.genome.set_phage_id_from_field("record_organism")
        self.assertEqual(self.genome.phage_id, "Trixie")

    def test_set_phage_id_from_field_9(self):
        """Check that the phage_id is set from the filename field."""
        self.genome.filename = "Trixie_Draft"
        self.genome.set_phage_id_from_field("filename")
        self.assertEqual(self.genome.phage_id, "Trixie")

    def test_set_phage_id_from_field_10(self):
        """Check that the phage_id is set from the
        record_description_phage_name field."""
        self.genome._record_description_phage_name = "Trixie_Draft"
        self.genome.set_phage_id_from_field("record_description_phage_name")
        self.assertEqual(self.genome.phage_id, "Trixie")

    def test_set_phage_id_from_field_11(self):
        """Check that the phage_id is set from the
        record_source_phage_name field."""
        self.genome._record_source_phage_name = "Trixie_Draft"
        self.genome.set_phage_id_from_field("record_source_phage_name")
        self.assertEqual(self.genome.phage_id, "Trixie")

    def test_set_phage_id_from_field_12(self):
        """Check that the phage_id is set from the
        record_organism_phage_name field."""
        self.genome._record_organism_phage_name = "Trixie_Draft"
        self.genome.set_phage_id_from_field("record_organism_phage_name")
        self.assertEqual(self.genome.phage_id, "Trixie")

    def test_set_phage_id_from_field_13(self):
        """Check that the phage_id is empty from an invalid field."""
        self.genome.set_phage_id_from_field("invalid")
        self.assertEqual(self.genome.phage_id, "")




    def test_set_host_from_field_1(self):
        """Check that the host is set from the record_name field."""
        self.genome.record_name = "Mycobacterium smegmatis"
        self.genome.set_host_from_field("record_name")
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_from_field_2(self):
        """Check that the host is set from the record_id field."""
        self.genome.record_id = "Mycobacterium smegmatis"
        self.genome.set_host_from_field("record_id")
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_from_field_3(self):
        """Check that the host is set from the record_accession field."""
        self.genome.record_accession = "Mycobacterium smegmatis"
        self.genome.set_host_from_field("record_accession")
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_from_field_4(self):
        """Check that the host is set from the record_description field."""
        self.genome.record_description = "Mycobacterium smegmatis"
        self.genome.set_host_from_field("record_description")
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_from_field_5(self):
        """Check that the host is set from the record_source field."""
        self.genome.record_source = "Mycobacterium smegmatis"
        self.genome.set_host_from_field("record_source")
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_from_field_6(self):
        """Check that the host is set from the record_organism field."""
        self.genome.record_organism = "Mycobacterium smegmatis"
        self.genome.set_host_from_field("record_organism")
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_from_field_7(self):
        """Check that the host is set from the filename field."""
        self.genome.filename = "Mycobacterium smegmatis"
        self.genome.set_host_from_field("filename")
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_from_field_8(self):
        """Check that the host is set from the
        record_description_host_name field."""
        self.genome._record_description_host_name = "Mycobacterium smegmatis"
        self.genome.set_host_from_field("record_description_host_name")
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_from_field_9(self):
        """Check that the host is set from the
        record_source_host_name field."""
        self.genome._record_source_host_name = "Mycobacterium smegmatis"
        self.genome.set_host_from_field("record_source_host_name")
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_from_field_10(self):
        """Check that the host is set from the
        record_organism_host_name field."""
        self.genome._record_organism_host_name = "Mycobacterium smegmatis"
        self.genome.set_host_from_field("record_organism_host_name")
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_from_field_11(self):
        """Check that the host is empty from an invalid field."""
        self.genome.set_host_from_field("invalid")
        self.assertEqual(self.genome.host, "")




    def test_set_cluster_1(self):
        """Check that standard Cluster is set appropriately."""
        cluster = "A"
        self.genome.set_cluster(cluster)
        self.assertEqual(self.genome.cluster, "A")

    def test_set_cluster_2(self):
        """Check that singleton Cluster is set appropriately."""
        cluster = "Singleton"
        self.genome.set_cluster(cluster)
        self.assertEqual(self.genome.cluster, "singleton")

    def test_set_cluster_3(self):
        """Check that whitespace is removed."""
        cluster = " A   "
        self.genome.set_cluster(cluster)
        self.assertEqual(self.genome.cluster, "A")




    def test_set_subcluster_1(self):
        """Check that standard Subcluster is set appropriately."""
        subcluster = "A2"
        self.genome.set_subcluster(subcluster, "none_string")
        self.assertEqual(self.genome.subcluster, "A2")

    def test_set_subcluster_2(self):
        """Check that whitespace is removed."""
        subcluster = "    A2    "
        self.genome.set_subcluster(subcluster, "none_string")
        self.assertEqual(self.genome.subcluster, "A2")

    def test_set_subcluster_3(self):
        """Check that None subcluster is set appropriately from None."""
        subcluster = None
        self.genome.set_subcluster(subcluster, "none_object")
        self.assertIsNone(self.genome.subcluster)

    def test_set_subcluster_4(self):
        """Check that none_string subcluster is set appropriately from None."""
        subcluster = None
        self.genome.set_subcluster(subcluster, "none_string")
        self.assertEqual(self.genome.subcluster, "none")

    def test_set_subcluster_5(self):
        """Check that empty_string subcluster is set appropriately from
        none_string."""
        subcluster = "none"
        self.genome.set_subcluster(subcluster, "empty_string")
        self.assertEqual(self.genome.subcluster, "")

    def test_set_subcluster_6(self):
        """Check that case is accounted for."""
        subcluster = "NONE"
        self.genome.set_subcluster(subcluster, "empty_string")
        self.assertEqual(self.genome.subcluster, "")




    def test_set_cluster_subcluster_1(self):
        """Check that None Cluster is set as singleton cluster_subcluster."""
        self.genome.subcluster = ""
        self.genome.cluster = None
        self.genome.set_cluster_subcluster()
        self.assertEqual(self.genome.cluster_subcluster, "singleton")

    def test_set_cluster_subcluster_2(self):
        """Check that singleton Cluster is set as
        singleton cluster_subcluster."""
        self.genome.subcluster = ""
        self.genome.cluster = "singleton"
        self.genome.set_cluster_subcluster()
        self.assertEqual(self.genome.cluster_subcluster, "singleton")

    def test_set_cluster_subcluster_3(self):
        """Check that Cluster is set as cluster_subcluster."""
        self.genome.subcluster = ""
        self.genome.cluster = "A"
        self.genome.set_cluster_subcluster()
        self.assertEqual(self.genome.cluster_subcluster, "A")

    def test_set_cluster_subcluster_4(self):
        """Check that Subcluster is set as cluster_subcluster."""
        self.genome.subcluster = "A1"
        self.genome.cluster = "A"
        self.genome.set_cluster_subcluster()
        self.assertEqual(self.genome.cluster_subcluster, "A1")

    def test_set_cluster_subcluster_5(self):
        """Check that Cluster is set when subcluster is None."""
        self.genome.subcluster = None
        self.genome.cluster = "A"
        self.genome.set_cluster_subcluster()
        self.assertEqual(self.genome.cluster_subcluster, "A")

    def test_set_cluster_subcluster_6(self):
        """Check that Cluster is set when subcluster is empty string."""
        self.genome.subcluster = ""
        self.genome.cluster = "A"
        self.genome.set_cluster_subcluster()
        self.assertEqual(self.genome.cluster_subcluster, "A")




    def test_set_date_last_modified_1(self):
        """Check that date_last_modified is set appropriately."""
        date_last_modified = datetime.strptime('1/1/1900', '%m/%d/%Y')
        self.genome.set_date_last_modified(date_last_modified, "empty_string")
        self.assertEqual(self.genome.date_last_modified, date_last_modified)

    def test_set_date_last_modified_2(self):
        """Check that a None date_last_modified is set appropriately
        from None."""
        date_last_modified = None
        self.genome.set_date_last_modified(date_last_modified, "empty_string")
        self.assertEqual(self.genome.date_last_modified, "")

    def test_set_date_last_modified_3(self):
        """Check that filled date_last_modified is set appropriately
        from None."""
        date_last_modified1 = None
        date_last_modified2 = datetime.strptime('1/1/0001', '%m/%d/%Y')
        self.genome.set_date_last_modified(date_last_modified1, \
                                            "empty_datetime_obj")
        self.assertEqual(self.genome.date_last_modified, date_last_modified2)

    def test_set_date_last_modified_4(self):
        """Check that None date_last_modified is set appropriately
        from None."""
        date_last_modified = None
        self.genome.set_date_last_modified(date_last_modified, "none_object")
        self.assertIsNone(self.genome.date_last_modified)

    def test_set_date_last_modified_5(self):
        """Check that None date_last_modified is set appropriately
        from empty string."""
        date_last_modified = ""
        self.genome.set_date_last_modified(date_last_modified, "none_object")
        self.assertIsNone(self.genome.date_last_modified)

    def test_set_date_last_modified_6(self):
        """Check that empty string date_last_modified is set appropriately from
        incorrect strategy."""
        date_last_modified = None
        self.genome.set_date_last_modified(date_last_modified, "invalid")
        self.assertEqual(self.genome.date_last_modified, None)




    def test_set_annotation_author_1(self):
        """Check that annotation_author to set to 1 when author
        is in author_set."""
        input_author = "Hatfull"
        self.genome.set_annotation_author(input_author)
        self.assertEqual(self.genome.annotation_author, 1)

    def test_set_annotation_author_2(self):
        """Check that annotation_author to set to 0 when author
        is not in author_set."""
        input_author = "Unknown author"
        self.genome.set_annotation_author(input_author)
        self.assertEqual(self.genome.annotation_author, 0)




    def test_tally_descriptions_1(self):
        """Check that no description tally is incremented."""
        self.genome.cds_features = [self.cds1, self.cds2]
        self.genome.tally_descriptions()

        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_primary_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_product_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_function_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_note_descriptions_tally, 0)

    def test_tally_descriptions_2(self):
        """Check that processed primary description tally
        is incremented."""
        self.cds1.processed_primary_description = "abcd"
        self.genome.cds_features = [self.cds1, self.cds2]
        self.genome.tally_descriptions()

        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_primary_descriptions_tally, 1)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_product_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_function_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_note_descriptions_tally, 0)

    def test_tally_descriptions_3(self):
        """Check that processed product description tally
        is incremented."""
        self.cds1.processed_product_description = "abcd"
        self.genome.cds_features = [self.cds1, self.cds2]
        self.genome.tally_descriptions()

        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_primary_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_product_descriptions_tally, 1)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_function_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_note_descriptions_tally, 0)

    def test_tally_descriptions_4(self):
        """Check that processed function description tally
        is incremented."""
        self.cds1.processed_function_description = "abcd"
        self.genome.cds_features = [self.cds1, self.cds2]
        self.genome.tally_descriptions()

        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_primary_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_product_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_function_descriptions_tally, 1)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_note_descriptions_tally, 0)

    def test_tally_descriptions_5(self):
        """Check that processed note description tally
        is incremented."""
        self.cds1.processed_note_description = "abcd"
        self.genome.cds_features = [self.cds1, self.cds2]
        self.genome.tally_descriptions()

        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_primary_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_product_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_function_descriptions_tally, 0)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_note_descriptions_tally, 1)

    def test_tally_descriptions_6(self):
        """Check that all description tallies are incremented."""
        self.cds1.processed_primary_description = "abcd"
        self.cds1.processed_product_description = "efgh"
        self.cds1.processed_function_description = "ijkl"
        self.cds1.processed_note_description = "mnop"
        self.genome.cds_features = [self.cds1, self.cds2]
        self.genome.tally_descriptions()

        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_primary_descriptions_tally, 1)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_product_descriptions_tally, 1)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_function_descriptions_tally, 1)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_note_descriptions_tally, 1)

    def test_tally_descriptions_7(self):
        """Check that all description tallies are incremented for
        all CDS features."""
        self.cds1.processed_primary_description = "abcd"
        self.cds1.processed_product_description = "efgh"
        self.cds1.processed_function_description = "ijkl"
        self.cds1.processed_note_description = "mnop"

        self.cds2.processed_primary_description = "ab"
        self.cds2.processed_product_description = "cd"
        self.cds2.processed_function_description = "ef"
        self.cds2.processed_note_description = "gh"

        self.genome.cds_features = [self.cds1, self.cds2]
        self.genome.tally_descriptions()

        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_primary_descriptions_tally, 2)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_product_descriptions_tally, 2)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_function_descriptions_tally, 2)
        with self.subTest():
            self.assertEqual( \
                self.genome._cds_processed_note_descriptions_tally, 2)




    def test_check_subcluster_structure_1(self):
        """Check that no error is produced if the
        non-empty subcluster is structured correctly."""
        self.genome.subcluster = "A1"
        self.genome.check_subcluster_structure()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_subcluster_structure_2(self):
        """Check that an error is produced if the
        non-empty subcluster is not structured correctly."""
        self.genome.subcluster = "A"
        self.genome.check_subcluster_structure()
        self.assertEqual(self.genome.evaluations[0].status, "error")

    def test_check_subcluster_structure_3(self):
        """Check that no error is produced if the
        subcluster is empty."""
        self.genome.subcluster = "none"
        self.genome.check_subcluster_structure()
        self.assertEqual(self.genome.evaluations[0].status, "not_evaluated")




    def test_check_cluster_structure_1(self):
        """Check that no error is produced if the
        non-empty cluster is structured correctly."""
        self.genome.cluster = "A"
        self.genome.check_cluster_structure()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_cluster_structure_2(self):
        """Check that an error is produced if the
        non-empty cluster is not structured correctly."""
        self.genome.cluster = "A1"
        self.genome.check_cluster_structure()
        self.assertEqual(self.genome.evaluations[0].status, "error")

    def test_check_cluster_structure_3(self):
        """Check that no error is produced if the
        cluster is empty."""
        self.genome.cluster = "none"
        self.genome.check_cluster_structure()
        self.assertEqual(self.genome.evaluations[0].status, "not_evaluated")




    def test_compare_cluster_subcluster_structure_1(self):
        """Check that compatible Cluster and subcluster
        do not produce an error."""
        self.genome.cluster = "A"
        self.genome.subcluster = "A1"
        self.genome.compare_cluster_subcluster_structure()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_compare_cluster_subcluster_structure_2(self):
        """Check that incompatible Cluster and subcluster
        produce an error."""
        self.genome.cluster = "A"
        self.genome.subcluster = "B1"
        self.genome.compare_cluster_subcluster_structure()
        self.assertEqual(self.genome.evaluations[0].status, "error")




    def test_check_nucleotides_1(self):
        """All nucleotides are in the alphabet."""
        alphabet = set(["A","B","C"])
        self.genome.sequence = "AB"
        self.genome.check_nucleotides(alphabet)
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_nucleotides_2(self):
        """Some nucleotides are not in the alphabet."""
        alphabet = set(["A","B","C"])
        self.genome.sequence = "AD"
        self.genome.check_nucleotides(alphabet)
        self.assertEqual(self.genome.evaluations[0].status, "error")




    def test_check_status_accession_1(self):
        """Check final status with accession."""
        self.genome.status = "final"
        self.genome.accession = "ABC123"
        self.genome.check_status_accession()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_status_accession_2(self):
        """Check final status with no accession."""
        self.genome.status = "final"
        self.genome.accession = ""
        self.genome.check_status_accession()
        self.assertEqual(self.genome.evaluations[0].status, "error")

    def test_check_status_accession_3(self):
        """Check draft status with no accession."""
        self.genome.status = "draft"
        self.genome.accession = ""
        self.genome.check_status_accession()
        self.assertEqual(self.genome.evaluations[0].status, "correct")






    def test_check_status_descriptions_1(self):
        """Check that draft genome with no descriptions does not produce
        an error."""
        self.genome.status = "draft"
        self.genome.check_status_descriptions()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_status_descriptions_2(self):
        """Check that draft genome with a description produces an error."""
        self.genome.status = "draft"
        self.genome._cds_processed_primary_descriptions_tally = 1
        self.genome.check_status_descriptions()
        self.assertEqual(self.genome.evaluations[0].status, "error")

    def test_check_status_descriptions_3(self):
        """Check that final genome with a description does not produce
        an error."""
        self.genome.status = "final"
        self.genome._cds_processed_primary_descriptions_tally = 1
        self.genome.check_status_descriptions()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_status_descriptions_4(self):
        """Check that final genome with no descriptions produces an error."""
        self.genome.status = "final"
        self.genome.check_status_descriptions()
        self.assertEqual(self.genome.evaluations[0].status, "error")

    def test_check_status_descriptions_5(self):
        """Check that gbk genome with no descriptions does not produce
        an error."""
        self.genome.status = "gbk"
        self.genome.check_status_descriptions()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_status_descriptions_6(self):
        """Check that gbk genome with descriptions does not produce
        an error."""
        self.genome.status = "gbk"
        self.genome._cds_processed_primary_descriptions_tally = 1
        self.genome.check_status_descriptions()
        self.assertEqual(self.genome.evaluations[0].status, "correct")












    def test_check_record_description_phage_name_1(self):
        """Check that no warning is produced."""
        self.genome.phage_id = "Trixie"
        self.genome._record_description_phage_name = "Trixie"
        self.genome.check_record_description_phage_name()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_record_description_phage_name_2(self):
        """Check that a warning is produced."""
        self.genome.phage_id = "L5"
        self.genome._record_description_phage_name = "Trixie"
        self.genome.check_record_description_phage_name()
        self.assertEqual(self.genome.evaluations[0].status, "error")




    def test_check_record_source_phage_name_1(self):
        """Check that no warning is produced."""
        self.genome.phage_id = "Trixie"
        self.genome._record_source_phage_name = "Trixie"
        self.genome.check_record_source_phage_name()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_record_source_phage_name_2(self):
        """Check that a warning is produced."""
        self.genome.phage_id = "L5"
        self.genome._record_source_phage_name = "Trixie"
        self.genome.check_record_source_phage_name()
        self.assertEqual(self.genome.evaluations[0].status, "error")




    def test_check_record_organism_phage_name_1(self):
        """Check that no warning is produced."""
        self.genome.phage_id = "Trixie"
        self.genome._record_organism_phage_name = "Trixie"
        self.genome.check_record_organism_phage_name()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_record_organism_phage_name_2(self):
        """Check that a warning is produced."""
        self.genome.phage_id = "L5"
        self.genome._record_organism_phage_name = "Trixie"
        self.genome.check_record_organism_phage_name()
        self.assertEqual(self.genome.evaluations[0].status, "error")




    def test_check_record_description_host_name_1(self):
        """Check that no warning is produced."""
        self.genome.host = "Mycobacterium"
        self.genome._record_description_host_name = "Mycobacterium"
        self.genome.check_record_description_host_name()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_record_description_host_name_2(self):
        """Check that a warning is produced."""
        self.genome.host = "Gordonia"
        self.genome._record_description_host_name = "Mycobacterium"
        self.genome.check_record_description_host_name()
        self.assertEqual(self.genome.evaluations[0].status, "error")




    def test_check_record_source_host_name_1(self):
        """Check that no warning is produced."""
        self.genome.host = "Mycobacterium"
        self.genome._record_source_host_name = "Mycobacterium"
        self.genome.check_record_source_host_name()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_record_source_host_name_2(self):
        """Check that a warning is produced."""
        self.genome.host = "Gordonia"
        self.genome._record_source_host_name = "Mycobacterium"
        self.genome.check_record_source_host_name()
        self.assertEqual(self.genome.evaluations[0].status, "error")




    def test_check_record_organism_host_name_1(self):
        """Check that no warning is produced."""
        self.genome.host = "Mycobacterium"
        self.genome._record_organism_host_name = "Mycobacterium"
        self.genome.check_record_organism_host_name()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_record_organism_host_name_2(self):
        """Check that a warning is produced."""
        self.genome.host = "Gordonia"
        self.genome._record_organism_host_name = "Mycobacterium"
        self.genome.check_record_organism_host_name()
        self.assertEqual(self.genome.evaluations[0].status, "error")




    def test_check_author_1(self):
        """Check that no warning is produced when author is expected
        and present."""
        self.genome.author = "Hatfull"
        self.genome.record_authors = "abcd; efgh; hatfull; xyz"
        self.genome.annotation_author = 1
        self.genome.check_author()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_author_2(self):
        """Check that no warning is produced when author is not expected
        and not present."""
        self.genome.author = "Hatfull"
        self.genome.record_authors = "abcd; efgh; xyz"
        self.genome.annotation_author = 0
        self.genome.check_author()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_author_3(self):
        """Check that warning is produced when author is expected
        and not present."""
        self.genome.author = "Hatfull"
        self.genome.record_authors = "abcd; efgh; xyz"
        self.genome.annotation_author = 1
        self.genome.check_author()
        self.assertEqual(self.genome.evaluations[0].status, "error")

    def test_check_author_4(self):
        """Check that warning is produced when author is not expected
        and present."""
        self.genome.author = "Hatfull"
        self.genome.record_authors = "abcd; efgh; hatfull; xyz"
        self.genome.annotation_author = 0
        self.genome.check_author()
        self.assertEqual(self.genome.evaluations[0].status, "error")










    def test_identify_unique_cds_start_end_ids_1(self):
        """Verify that both sets are computed."""
        self.genome._cds_start_end_ids = \
            [(1, 5), (2, 10), (10, 2), (2, 10)]
        expected_unique_set = set([(1, 5), (10, 2)])
        expected_duplicate_set = set([(2, 10)])
        self.genome.identify_unique_cds_start_end_ids()
        with self.subTest():
            self.assertEqual(self.genome._cds_unique_start_end_ids, \
                expected_unique_set)
        with self.subTest():
            self.assertEqual(self.genome._cds_duplicate_start_end_ids, \
                expected_duplicate_set)




    def test_identify_unique_cds_end_strand_ids_1(self):
        """Verify that both sets are computed."""
        self.genome._cds_end_strand_ids = \
            [(1, "forward"), (2, "reverse"), (2, "forward"), (2, "reverse")]
        expected_unique_set = set([(1, "forward"), (2, "forward")])
        expected_duplicate_set = set([(2, "reverse")])
        self.genome.identify_unique_cds_end_strand_ids()
        with self.subTest():
            self.assertEqual(self.genome._cds_unique_end_strand_ids, \
                expected_unique_set)
        with self.subTest():
            self.assertEqual(self.genome._cds_duplicate_end_strand_ids, \
                expected_duplicate_set)




    def test_check_cds_start_end_ids_1(self):
        """Verify that no warning is produced."""
        self.genome.check_cds_start_end_ids()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_cds_start_end_ids_2(self):
        """Verify that a warning is produced."""
        self.genome._cds_duplicate_start_end_ids = set([(2, 10)])
        self.genome.check_cds_start_end_ids()
        self.assertEqual(self.genome.evaluations[0].status, "error")




    def test_check_cds_end_strand_ids_1(self):
        """Verify that no warning is produced."""
        self.genome.check_cds_end_strand_ids()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_cds_end_strand_ids_2(self):
        """Verify that a warning is produced."""
        self.genome._cds_duplicate_end_strand_ids = set([(2, "forward")])
        self.genome.check_cds_end_strand_ids()
        self.assertEqual(self.genome.evaluations[0].status, "error")




    def test_parse_record_description_1(self):
        """Verify empty string is parsed correctly."""
        self.genome.record_description = ""
        self.genome.parse_record_description()
        expected_phage = ""
        expected_host = ""
        with self.subTest():
            self.assertEqual(\
                self.genome._record_description_phage_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.genome._record_description_host_name, expected_host)

    def test_parse_record_description_2(self):
        """Verify string is parsed correctly."""
        self.genome.record_description = "asdf Mycobacterium phage Trixie."
        self.genome.parse_record_description()
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        with self.subTest():
            self.assertEqual(\
                self.genome._record_description_phage_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.genome._record_description_host_name, expected_host)




    def test_parse_record_source_1(self):
        """Verify empty string is parsed correctly."""
        self.genome.record_source = ""
        self.genome.parse_record_source()
        expected_phage = ""
        expected_host = ""
        with self.subTest():
            self.assertEqual(\
                self.genome._record_source_phage_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.genome._record_source_host_name, expected_host)

    def test_parse_record_source_2(self):
        """Verify string is parsed correctly."""
        self.genome.record_source = "asdf Mycobacterium phage Trixie."
        self.genome.parse_record_source()
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        with self.subTest():
            self.assertEqual(\
                self.genome._record_source_phage_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.genome._record_source_host_name, expected_host)




    def test_parse_record_organism_1(self):
        """Verify empty string is parsed correctly."""
        self.genome.record_organism = ""
        self.genome.parse_record_organism()
        expected_phage = ""
        expected_host = ""
        with self.subTest():
            self.assertEqual(\
                self.genome._record_organism_phage_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.genome._record_organism_host_name, expected_host)

    def test_parse_record_organism_2(self):
        """Verify string is parsed correctly."""
        self.genome.record_organism = "asdf Mycobacterium phage Trixie."
        self.genome.parse_record_organism()
        expected_phage = "Trixie"
        expected_host = "Mycobacterium"
        with self.subTest():
            self.assertEqual(\
                self.genome._record_organism_phage_name, expected_phage)
        with self.subTest():
            self.assertEqual(\
                self.genome._record_organism_host_name, expected_host)




    def test_split_cluster_subcluster_1(self):
        """Verify split to only cluster, and subcluster is 'none'."""
        self.genome.cluster_subcluster = "A"
        self.genome.split_cluster_subcluster()
        cluster = "A"
        subcluster = "none"
        with self.subTest():
            self.assertEqual(self.genome.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, subcluster)

    def test_split_cluster_subcluster_2(self):
        """Verify split to only cluster, and subcluster is None."""
        self.genome.cluster_subcluster = "A"
        self.genome.split_cluster_subcluster("none_object")
        cluster = "A"
        subcluster = None
        with self.subTest():
            self.assertEqual(self.genome.cluster, cluster)
        with self.subTest():
            self.assertIsNone(self.genome.subcluster)

    def test_split_cluster_subcluster_3(self):
        """Verify split to only cluster, and subcluster is ''."""
        self.genome.cluster_subcluster = "A"
        self.genome.split_cluster_subcluster("empty_string")
        cluster = "A"
        subcluster = ""
        with self.subTest():
            self.assertEqual(self.genome.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, subcluster)

    def test_split_cluster_subcluster_4(self):
        """Verify singleton split to only cluster."""
        self.genome.cluster_subcluster = "singleton"
        self.genome.split_cluster_subcluster("empty_string")
        cluster = "singleton"
        subcluster = ""
        with self.subTest():
            self.assertEqual(self.genome.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, subcluster)

    def test_split_cluster_subcluster_5(self):
        """Verify split to both cluster and subcluster."""
        self.genome.cluster_subcluster = "A15"
        self.genome.split_cluster_subcluster("empty_string")
        cluster = "A"
        subcluster = "A15"
        with self.subTest():
            self.assertEqual(self.genome.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, subcluster)

    def test_split_cluster_subcluster_6(self):
        """Verify no split, and output is 'none'."""
        self.genome.cluster_subcluster = "A1B2"
        self.genome.split_cluster_subcluster()
        cluster = "none"
        subcluster = "none"
        with self.subTest():
            self.assertEqual(self.genome.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, subcluster)

    def test_split_cluster_subcluster_7(self):
        """Verify no split, and output is ''."""
        self.genome.cluster_subcluster = "A1B2"
        self.genome.cluster = ""
        self.genome.subcluster = ""
        self.genome.split_cluster_subcluster("empty_string")
        cluster = ""
        subcluster = ""
        with self.subTest():
            self.assertEqual(self.genome.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, subcluster)

    def test_split_cluster_subcluster_8(self):
        """Verify no split, and output is None."""
        self.genome.cluster_subcluster = "A1B2"
        self.genome.cluster = ""
        self.genome.subcluster = ""
        self.genome.split_cluster_subcluster("none_object")
        cluster = None
        subcluster = None
        with self.subTest():
            self.assertEqual(self.genome.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, subcluster)

    def test_split_cluster_subcluster_9(self):
        """Verify no change to cluster and subcluster."""
        self.genome.cluster_subcluster = ""
        self.genome.cluster = "B"
        self.genome.subcluster = "B10"
        self.genome.split_cluster_subcluster("none_object")
        cluster = "B"
        subcluster = "B10"
        with self.subTest():
            self.assertEqual(self.genome.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, subcluster)

    def test_split_cluster_subcluster_10(self):
        """Verify no change to cluster and subcluster."""
        self.genome.cluster_subcluster = None
        self.genome.cluster = "B"
        self.genome.subcluster = "B10"
        self.genome.split_cluster_subcluster("none_object")
        cluster = "B"
        subcluster = "B10"
        with self.subTest():
            self.assertEqual(self.genome.cluster, cluster)
        with self.subTest():
            self.assertEqual(self.genome.subcluster, subcluster)




    def test_check_cds_feature_tally_1(self):
        """Verify no error is encountered when there is one CDS feature."""
        self.genome._cds_features_tally = 1
        self.genome.check_cds_feature_tally()
        self.assertEqual(self.genome.evaluations[0].status, "correct")

    def test_check_cds_feature_tally_2(self):
        """Verify error is encountered when there is no CDS feature."""
        self.genome._cds_features_tally = 0
        self.genome.check_cds_feature_tally()
        self.assertEqual(self.genome.evaluations[0].status, "error")









if __name__ == '__main__':
    unittest.main()
