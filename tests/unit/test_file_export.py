"""Tests the functionality of the unique functions in the file_export pipeline"""

import unittest
from pdm_utils.classes import genome, cds
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from pdm_utils.pipelines.db_export import file_export
from pdm_utils.functions import flat_files

class TestDatabaseToFile(unittest.TestCase):

    def setUp(self):
        self.test_genome = genome.Genome()
        self.test_genome.seq = Seq("ATA")
        self.test_seqrecord = flat_files.\
                genome_to_seqrecord(self.test_genome)
        self.cds_1 = cds.Cds()
        self.cds_1.strand = "forward"
        self.cds_1.coordinate_format = "0_half_open" 
        self.cds_1.left = 0
        self.cds_1.right = 1
        self.cds_2 = cds.Cds()
        self.cds_2.strand = "forward"
        self.cds_2.coordinate_format = "0_half_open"
        self.cds_2.left = 1
        self.cds_2.right = 2
        self.cds_3 = cds.Cds()
        self.cds_3.strand = "forward"
        self.cds_3.coordinate_format = "0_half_open"
        self.cds_3.left = 2
        self.cds_3.right = 3
        self.test_version = {"version" : 1, "schema_version" : 1}
        
    def test_set_cds_seqfeatures_1(self):
        test_cds_list = [self.cds_1, self.cds_2, self.cds_3]
        self.test_genome.cds_features = test_cds_list
        file_export.set_cds_seqfeatures\
                (self.test_genome)

        for test_cds in self.test_genome.cds_features:
            self.assertTrue(test_cds.seqfeature != None)
       
        self.assertEqual\
                (self.test_genome.cds_features[0].left, 0)
        self.assertEqual\
                (self.test_genome.cds_features[1].left, 1)
        self.assertEqual\
                (self.test_genome.cds_features[2].left, 2)

    def test_set_cds_seqfeatures_2(self):
        test_cds_list = [self.cds_3, self.cds_2, self.cds_1]
        self.test_genome.cds_features = test_cds_list
        file_export.set_cds_seqfeatures\
                (self.test_genome)

        for test_cds in self.test_genome.cds_features:
            self.assertTrue(test_cds.seqfeature != None)
        
        self.assertEqual\
                (self.test_genome.cds_features[0].left, 0)
        self.assertEqual\
                (self.test_genome.cds_features[1].left, 1)
        self.assertEqual\
                (self.test_genome.cds_features[2].left, 2)

    def test_set_cds_seqfeatures_3(self):
        with self.assertRaises(TypeError):
            file_export.set_cds_seqfeatures(None)

    def test_append_database_version_1(self):
        file_export.append_database_version\
                (self.test_seqrecord, self.test_version)

        comment_tuple = self.test_seqrecord.annotations["comment"]
        self.assertEqual(comment_tuple[4],\
                "Database Version: 1; Schema Version: 1")

    def test_append_database_version_2(self):
        with self.assertRaises(ValueError):
            file_export.append_database_version\
                    (self.test_seqrecord, {})

        with self.assertRaises(TypeError):
            file_export.append_database_version\
                    (None, self.test_version)
