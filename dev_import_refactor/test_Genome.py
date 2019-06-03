""" Unit tests for the Genome class."""


import unittest
import Genome
import Cds
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



    def test_set_evaluation_1(self):
        """Set an empty evaluation object."""
        self.genome.set_evaluation("none")
        self.assertEqual(len(self.genome.evaluations), 1)

    def test_set_evaluation_2(self):
        """Set a warning evaluation object."""
        self.genome.set_evaluation("warning","message1")
        self.assertEqual(len(self.genome.evaluations), 1)

    def test_set_evaluation_3(self):
        """Set an error evaluation object."""
        self.genome.set_evaluation("error","message1","message2")
        self.assertEqual(len(self.genome.evaluations), 1)





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
        self.genome.set_host(host)
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_2(self):
        """Check that whitespace is removed."""
        host = "  Mycobacterium smegmatis  "
        self.genome.set_host(host)
        self.assertEqual(self.genome.host, "Mycobacterium")

    def test_set_host_3(self):
        """Check that none is set appropriately."""
        host = "  none  "
        self.genome.set_host(host)
        self.assertEqual(self.genome.host, "")



    def test_set_sequence_1(self):
        """Check that sequence is set appropriately."""
        seq = "abcd"
        self.genome.set_sequence(seq)
        with self.subTest():
            self.assertEqual(self.genome.sequence, "ABCD")
        with self.subTest():
            self.assertEqual(self.genome._length, 4)




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
        features_list = [0,1,2,3]
        self.genome.set_cds_features(features_list)
        with self.subTest():
            self.assertEqual(len(self.genome.cds_features), 4)
        with self.subTest():
            self.assertEqual(self.genome._cds_features_tally, 4)


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







    def test_set_cluster_subcluster_1(self):
        """Check that None Cluster is set as singleton cluster_subcluster."""
        subcluster = ""
        cluster = None
        self.genome.set_cluster_subcluster(cluster, subcluster)
        self.assertEqual(self.genome.cluster_subcluster, "singleton")

    def test_set_cluster_subcluster_2(self):
        """Check that singleton Cluster is set as
        singleton cluster_subcluster."""
        subcluster = ""
        cluster = "singleton"
        self.genome.set_cluster_subcluster(cluster, subcluster)
        self.assertEqual(self.genome.cluster_subcluster, "singleton")

    def test_set_cluster_subcluster_3(self):
        """Check that Cluster is set as cluster_subcluster."""
        subcluster = ""
        cluster = "A"
        self.genome.set_cluster_subcluster(cluster, subcluster)
        self.assertEqual(self.genome.cluster_subcluster, "A")

    def test_set_cluster_subcluster_4(self):
        """Check that Subcluster is set as cluster_subcluster."""
        subcluster = "A1"
        cluster = "A"
        self.genome.set_cluster_subcluster(cluster, subcluster)
        self.assertEqual(self.genome.cluster_subcluster, "A1")





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
        date_last_modified2 = datetime.strptime('1/1/1900', '%m/%d/%Y')
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
        self.assertEqual(self.genome.date_last_modified, "")





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











    def test_check_nucleotides_1(self):
        """All nucleotides are in the alphabet."""
        alphabet = set(["A","B","C"])
        self.genome.sequence = "AB"
        self.genome.check_nucleotides(alphabet)
        self.assertEqual(len(self.genome.evaluations), 0)

    def test_check_nucleotides_2(self):
        """Some nucleotides are not in the alphabet."""
        alphabet = set(["A","B","C"])
        self.genome.sequence = "AD"
        self.genome.check_nucleotides(alphabet)
        self.assertEqual(len(self.genome.evaluations), 1)






    def test_check_status_accession_1(self):
        """Check final status with accession."""
        self.genome.status = "final"
        self.genome.accession = "ABC123"
        self.genome.check_status_accession()
        self.assertEqual(len(self.genome.evaluations), 0)

    def test_check_status_accession_2(self):
        """Check final status with no accession."""
        self.genome.status = "final"
        self.genome.accession = ""
        self.genome.check_status_accession()
        self.assertEqual(len(self.genome.evaluations), 1)

    def test_check_status_accession_3(self):
        """Check draft status with no accession."""
        self.genome.status = "draft"
        self.genome.accession = ""
        self.genome.check_status_accession()
        self.assertEqual(len(self.genome.evaluations), 0)






    def test_check_status_descriptions_1(self):
        """Check that draft genome with no descriptions does not produce
        an error."""
        self.genome.status = "draft"
        self.genome.check_status_descriptions()
        self.assertEqual(len(self.genome.evaluations), 0)

    def test_check_status_descriptions_2(self):
        """Check that draft genome with a description produces an error."""
        self.genome.status = "draft"
        self.genome._cds_processed_primary_descriptions_tally = 1
        self.genome.check_status_descriptions()
        self.assertEqual(len(self.genome.evaluations), 1)

    def test_check_status_descriptions_3(self):
        """Check that final genome with a description does not produce
        an error."""
        self.genome.status = "final"
        self.genome._cds_processed_primary_descriptions_tally = 1
        self.genome.check_status_descriptions()
        self.assertEqual(len(self.genome.evaluations), 0)

    def test_check_status_descriptions_4(self):
        """Check that final genome with no descriptions produces an error."""
        self.genome.status = "final"
        self.genome.check_status_descriptions()
        self.assertEqual(len(self.genome.evaluations), 1)

    def test_check_status_descriptions_5(self):
        """Check that gbk genome with no descriptions does not produce
        an error."""
        self.genome.status = "gbk"
        self.genome.check_status_descriptions()
        self.assertEqual(len(self.genome.evaluations), 0)

    def test_check_status_descriptions_6(self):
        """Check that gbk genome with descriptions does not produce
        an error."""
        self.genome.status = "gbk"
        self.genome._cds_processed_primary_descriptions_tally = 1
        self.genome.check_status_descriptions()
        self.assertEqual(len(self.genome.evaluations), 0)
















if __name__ == '__main__':
    unittest.main()
