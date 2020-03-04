"""Unit tests for misc. functions that interact with
GenBank-formatted flat files."""


from datetime import datetime
import pathlib
import unittest

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqFeature import ExactPosition, BeforePosition, Reference

from pdm_utils.constants import constants
from pdm_utils.classes import cds, source, genome
from pdm_utils.functions import basic, flat_files

class TestFlatFileFunctions1(unittest.TestCase):


    def setUp(self):
        self.cds_ftr = cds.Cds()
        self.src_ftr = source.Source()
        self.gnm = genome.Genome()

    def test_parse_coordinates_1(self):
        """Verify non-compound location is parsed correctly."""

        seqfeature = SeqFeature(FeatureLocation(
            ExactPosition(2), ExactPosition(10)),
            type = "CDS",
            strand = 1)
        output_start, output_stop, parts = \
            flat_files.parse_coordinates(seqfeature)
        with self.subTest():
            self.assertEqual(output_start, 2)
        with self.subTest():
            self.assertEqual(output_stop, 10)
        with self.subTest():
            self.assertEqual(parts, 1)

    def test_parse_coordinates_2(self):
        """Verify 1 strand 2-part compound location is parsed correctly."""

        seqfeature = SeqFeature(CompoundLocation([
                        FeatureLocation(
                            ExactPosition(2),
                            ExactPosition(10),
                            strand=1),
                        FeatureLocation(
                            ExactPosition(8),
                            ExactPosition(20),
                            strand=1)],
                        'join'),
                        type='CDS',
                        location_operator='join')
        output_start, output_stop, parts = \
            flat_files.parse_coordinates(seqfeature)
        with self.subTest():
            self.assertEqual(output_start, 2)
        with self.subTest():
            self.assertEqual(output_stop, 20)
        with self.subTest():
            self.assertEqual(parts, 2)

    def test_parse_coordinates_3(self):
        """Verify -1 strand 2-part compound location is parsed correctly."""

        seqfeature = SeqFeature(CompoundLocation([
                        FeatureLocation(
                            ExactPosition(2),
                            ExactPosition(10),
                            strand=-1),
                        FeatureLocation(
                            ExactPosition(8),
                            ExactPosition(20),
                            strand=-1)],
                        'join'),
                        type='CDS',
                        location_operator='join')
        output_start, output_stop, parts = \
            flat_files.parse_coordinates(seqfeature)
        with self.subTest():
            self.assertEqual(output_start, 8)
        with self.subTest():
            self.assertEqual(output_stop, 10)
        with self.subTest():
            self.assertEqual(parts, 2)

    def test_parse_coordinates_4(self):
        """Verify 1 strand 2-part compound location that wraps around
        genome end is parsed correctly."""

        # Wrap-around feature, directly copied from
        # Biopython-parsed ET08 flat file.
        seqfeature = SeqFeature(CompoundLocation([
                        FeatureLocation(
                            ExactPosition(154873),
                            ExactPosition(155445),
                            strand=1),
                        FeatureLocation(
                            ExactPosition(0),
                            ExactPosition(4),
                            strand=1)],
                        'join'),
                        type='CDS',
                        location_operator='join')
        output_start, output_stop, parts = \
            flat_files.parse_coordinates(seqfeature)
        with self.subTest():
            self.assertEqual(output_start, 154873)
        with self.subTest():
            self.assertEqual(output_stop, 4)
        with self.subTest():
            self.assertEqual(parts, 2)

    def test_parse_coordinates_5(self):
        """Verify -1 strand 2-part compound location that wraps around
        genome end is parsed correctly."""

        # Wrap-around feature, directly copied from
        # Biopython-parsed ET08 flat file modified in CLC so that the
        # the wrap-around gene is on the bottom strand and then
        # exported as a GenBank-formatted flat file.
        seqfeature = SeqFeature(CompoundLocation([
                        FeatureLocation(
                            ExactPosition(0),
                            ExactPosition(4),
                            strand=-1),
                        FeatureLocation(
                            ExactPosition(154873),
                            ExactPosition(155445),
                            strand=-1)],
                        'join'),
                        type='CDS',
                        location_operator='join')
        output_start, output_stop, parts = \
            flat_files.parse_coordinates(seqfeature)
        with self.subTest():
            self.assertEqual(output_start, 154873)
        with self.subTest():
            self.assertEqual(output_stop, 4)
        with self.subTest():
            self.assertEqual(parts, 2)

    def test_parse_coordinates_6(self):
        """Verify undefined strand 2-part compound location is not parsed."""

        seqfeature = SeqFeature(CompoundLocation([
                        FeatureLocation(
                            ExactPosition(2),
                            ExactPosition(10),
                            strand=None),
                        FeatureLocation(
                            ExactPosition(8),
                            ExactPosition(20),
                            strand=None)],
                        'join'),
                        type='CDS',
                        location_operator='join')
        output_start, output_stop, parts = \
            flat_files.parse_coordinates(seqfeature)
        with self.subTest():
            self.assertEqual(output_start, -1)
        with self.subTest():
            self.assertEqual(output_stop, -1)
        with self.subTest():
            self.assertEqual(parts, 0)

    def test_parse_coordinates_7(self):
        """Verify 1 strand 3-part compound location is not parsed."""

        seqfeature = SeqFeature(CompoundLocation([
                        FeatureLocation(
                            ExactPosition(2),
                            ExactPosition(10),
                            strand=1),
                        FeatureLocation(
                            ExactPosition(8),
                            ExactPosition(20),
                            strand=1),
                        FeatureLocation(
                            ExactPosition(30),
                            ExactPosition(50),
                            strand=1)],
                        'join'),
                        type='CDS',
                        location_operator='join')
        output_start, output_stop, parts = \
            flat_files.parse_coordinates(seqfeature)
        with self.subTest():
            self.assertEqual(output_start, -1)
        with self.subTest():
            self.assertEqual(output_stop, -1)
        with self.subTest():
            self.assertEqual(parts, 3)

    def test_parse_coordinates_8(self):
        """Verify location of invalid data type is not parsed."""

        seqfeature = SeqFeature(None, type = "CDS", strand = None)
        output_start, output_stop, parts = \
            flat_files.parse_coordinates(seqfeature)
        with self.subTest():
            self.assertEqual(output_start, -1)
        with self.subTest():
            self.assertEqual(output_stop, -1)
        with self.subTest():
            self.assertEqual(parts, 0)

    def test_parse_coordinates_9(self):
        """Verify non-compound location with fuzzy start coordinate
        is parsed correctly."""
        seqfeature = SeqFeature(
                        FeatureLocation(
                            BeforePosition(2),
                            ExactPosition(10)),
                        type = "CDS",
                        strand = 1)
        output_start, output_stop, parts = \
            flat_files.parse_coordinates(seqfeature)
        with self.subTest():
            self.assertEqual(output_start, -1)
        with self.subTest():
            self.assertEqual(output_stop, 10)
        with self.subTest():
            self.assertEqual(parts, 1)

    def test_parse_coordinates_10(self):
        """Verify non-compound location with fuzzy stop coordinate
        is parsed correctly."""
        seqfeature = SeqFeature(
                        FeatureLocation(
                            ExactPosition(2),
                            BeforePosition(10)),
                        type = "CDS",
                        strand = 1)
        output_start, output_stop, parts = \
            flat_files.parse_coordinates(seqfeature)
        with self.subTest():
            self.assertEqual(output_start, 2)
        with self.subTest():
            self.assertEqual(output_stop, -1)
        with self.subTest():
            self.assertEqual(parts, 1)

    def test_parse_coordinates_11(self):
        """Verify 1 strand 2-part compound location with fuzzy start
        coordinate is parsed correctly."""
        seqfeature = SeqFeature(CompoundLocation([
                        FeatureLocation(
                            BeforePosition(2),
                            ExactPosition(10),
                            strand=1),
                        FeatureLocation(
                            ExactPosition(8),
                            ExactPosition(20),
                            strand=1)],
                        "join"),
                        type="CDS",
                        location_operator="join")
        output_start, output_stop, parts = \
            flat_files.parse_coordinates(seqfeature)
        with self.subTest():
            self.assertEqual(output_start, -1)
        with self.subTest():
            self.assertEqual(output_stop, 20)
        with self.subTest():
            self.assertEqual(parts, 2)




    def test_create_seqfeature_dictionary_1(self):
        """Verify feature dictionary is constructed correctly with
        one feature."""
        feature_list = [SeqFeature(type = "CDS")]
        feature_dict = flat_files.create_seqfeature_dictionary(feature_list)
        with self.subTest():
            self.assertEqual(len(feature_dict.keys()), 1)
        with self.subTest():
            self.assertEqual(len(feature_dict["CDS"]), 1)


    def test_create_seqfeature_dictionary_2(self):
        """Verify feature dictionary is constructed correctly with
        no features."""
        feature_list = []
        feature_dict = flat_files.create_seqfeature_dictionary(feature_list)
        with self.subTest():
            self.assertEqual(len(feature_dict.keys()), 0)


    def test_create_seqfeature_dictionary_3(self):
        """Verify feature dictionary is constructed correctly with
        several different feature types."""
        feature_list = [
            SeqFeature(type = "CDS"),
            SeqFeature(type = "CDS"),
            SeqFeature(type = "tRNA"),
            SeqFeature(type = "tmRNA"),
            SeqFeature(type = "other"),
            SeqFeature(type = "gene")]
        feature_dict = flat_files.create_seqfeature_dictionary(feature_list)
        with self.subTest():
            self.assertEqual(len(feature_dict.keys()), 5)
        with self.subTest():
            self.assertEqual(len(feature_dict["CDS"]), 2)
        with self.subTest():
            self.assertEqual(len(feature_dict["tRNA"]), 1)
        with self.subTest():
            self.assertEqual(len(feature_dict["tmRNA"]), 1)
        with self.subTest():
            self.assertEqual(len(feature_dict["other"]), 1)
        with self.subTest():
            self.assertEqual(len(feature_dict["gene"]), 1)




    def test_genome_to_seqrecord_1(self):
        """Verify that genome_to_seqrecord can initialize
        a SeqRecord given the basic conditions"""
        self.gnm.seq = Seq("ATA")
        self.gnm.name = "Trixie"
        self.gnm.id = "Trixie"
        self.gnm.cds_features = []
        self.gnm.host_genus = "Mycobacterium"

        record = flat_files.genome_to_seqrecord(self.gnm)
        self.assertEqual(record.name, "Trixie")
        self.assertEqual(record.features, [])
        self.assertEqual(record.description,\
                "Mycobacterium phage Trixie, complete genome")
        self.assertEqual(record.seq, Seq("ATA"))

    def test_genome_to_seqrecord_2(self):
        """Verify that genome_to_seqrecord can correctly
        populate seqrecord annotations"""

        self.gnm.seq = Seq("ATA")
        self.gnm.date = "2019"
        self.gnm.accession = "gnm12345"
        self.gnm.host_genus = "Mycobacterium"
        self.gnm.id = "Trixie"
        self.gnm.cluster = "A"
        self.gnm.subcluster = "A2"
        self.gnm.annotation_status = "1"
        self.gnm.annotation_author = 1
        self.gnm.retrieve_record = 1

        record = flat_files.genome_to_seqrecord(self.gnm)
        record_comments = record.annotations["comment"]

        self.assertEqual(record.annotations["date"], "2019")
        self.assertEqual(record.annotations["source"],
                         "Mycobacterium phage Trixie")
        self.assertEqual(record_comments[0], "Cluster: A; Subcluster: A2")
        self.assertEqual(record_comments[2],
                         "Annotation Status: 1; Annotation Author: 1")
        self.assertEqual(record_comments[3],
                         "RetrieveRecord: 1")

    def test_genome_to_seqrecord_3(self):
        self.gnm = None
        with self.assertRaises(AssertionError):
            record = flat_files.genome_to_seqrecord(self.gnm)

    def test_genome_to_seqrecord_4(self):
        # self.gnm = genome.Genome()
        self.gnm = ""
        with self.assertRaises(AttributeError):
            record = flat_files.genome_to_seqrecord(self.gnm)



    def test_get_seqrecord_description(self):
        """
        Unittest that tests flat_files.get_seqrecord_description()
        helper function
        """
        self.gnm.seq = Seq("ATA")
        self.gnm.name = "Trixie"
        self.gnm.id = "Trixie"
        self.gnm.cds_features = []
        self.gnm.host_genus = "Mycobacterium"

        description = flat_files.get_seqrecord_description(self.gnm)
        self.assertEqual(description, "Mycobacterium phage Trixie, "
                                      "complete genome")

    def test_get_seqrecord_annotations(self):
        """
        Unittest that tests flat_files.get_seqrecord_annotations()
        helper function
        """

        self.gnm.seq = Seq("ATA")
        self.gnm.date = "2019"
        self.gnm.accession = "gnm12345"
        self.gnm.host_genus = "Mycobacterium"
        self.gnm.id = "Trixie"
        self.gnm.cluster = "A"
        self.gnm.subcluster = "A2"
        self.gnm.annotation_status = "1"
        self.gnm.annotation_author = 1
        self.gnm.retrieve_record = 1

        annotations = flat_files.get_seqrecord_annotations(self.gnm)
        comments = annotations["comment"]

        self.assertEqual(annotations["date"], "2019")
        self.assertEqual(annotations["source"],
                         "Mycobacterium phage Trixie")
        self.assertEqual(comments[0], "Cluster: A; Subcluster: A2")
        self.assertEqual(comments[2],
                         "Annotation Status: 1; Annotation Author: 1")
        self.assertEqual(comments[3], "RetrieveRecord: 1")

    def test_get_seqrecord_annotations_comments(self):
        """
        Unittest that tests flat_files.get_seqrecord_annotations_comments()
        helper function
        """
        self.gnm.seq = Seq("ATA")
        self.gnm.date = "2019"
        self.gnm.accession = "gnm12345"
        self.gnm.host_genus = "Mycobacterium"
        self.gnm.id = "Trixie"
        self.gnm.cluster = "A"
        self.gnm.subcluster = "A2"
        self.gnm.annotation_status = "1"
        self.gnm.annotation_author = 1
        self.gnm.retrieve_record = 1

        comments = flat_files.get_seqrecord_annotations_comments(self.gnm)

        self.assertEqual(comments[0], "Cluster: A; Subcluster: A2")
        self.assertEqual(comments[2],
                         "Annotation Status: 1; Annotation Author: 1")
        self.assertEqual(comments[3], "RetrieveRecord: 1")




class TestFlatFileFunctions2(unittest.TestCase):


    def setUp(self):
        self.qualifier_dict = {"locus_tag": ["SEA_L5_1"],
                               "translation": ["ABCDE"],
                               "transl_table": ["11"],
                               "product": [" unknown "],
                               "function": [" hypothetical protein "],
                               "note": [" gp5 "],
                               "gene": ["2"]}

        self.seqfeature = SeqFeature(FeatureLocation(
                            ExactPosition(2), ExactPosition(10)),
                            type="CDS",
                            strand=1,
                            qualifiers=self.qualifier_dict)


    def test_parse_cds_seqfeature_1(self):
        """Verify CDS features is parsed."""
        cds_ftr = flat_files.parse_cds_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertEqual(cds_ftr.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(cds_ftr._locus_tag_num, "1")
        with self.subTest():
            self.assertEqual(cds_ftr.orientation, "F")
        with self.subTest():
            self.assertEqual(cds_ftr.start, 2)
        with self.subTest():
            self.assertEqual(cds_ftr.stop, 10)
        with self.subTest():
            self.assertEqual(cds_ftr.parts, 1)
        with self.subTest():
            self.assertEqual(cds_ftr.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(cds_ftr.translation, "ABCDE")
        with self.subTest():
            self.assertEqual(cds_ftr.translation_length, 5)
        with self.subTest():
            self.assertEqual(cds_ftr.length, 8)
        with self.subTest():
            self.assertEqual(cds_ftr.translation_table, 11)
        with self.subTest():
            self.assertEqual(cds_ftr.raw_product, "unknown")
        with self.subTest():
            self.assertEqual(cds_ftr.product, "")
        with self.subTest():
            self.assertEqual(cds_ftr.raw_function, "hypothetical protein")
        with self.subTest():
            self.assertEqual(cds_ftr.function, "")
        with self.subTest():
            self.assertEqual(cds_ftr.raw_note, "gp5")
        with self.subTest():
            self.assertEqual(cds_ftr.note, "")
        with self.subTest():
            self.assertEqual(cds_ftr.gene, "2")
        with self.subTest():
            self.assertTrue(isinstance(cds_ftr.seqfeature, SeqFeature))
        with self.subTest():
            self.assertEqual(cds_ftr.name, "2")


    def test_parse_cds_seqfeature_2(self):
        """Verify CDS features is parsed with no locus tag."""
        self.qualifier_dict.pop("locus_tag")
        cds_ftr = flat_files.parse_cds_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertEqual(cds_ftr.start, 2)
        with self.subTest():
            self.assertEqual(cds_ftr.locus_tag, "")
        with self.subTest():
            self.assertEqual(cds_ftr._locus_tag_num, "")
        with self.subTest():
            self.assertEqual(cds_ftr.gene, "2")
        with self.subTest():
            self.assertEqual(cds_ftr.name, "2")


    def test_parse_cds_seqfeature_3(self):
        """Verify CDS features is parsed with no translation."""
        self.qualifier_dict.pop("translation")
        cds_ftr = flat_files.parse_cds_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertEqual(cds_ftr.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(cds_ftr.translation, "")
        with self.subTest():
            self.assertEqual(cds_ftr.translation_length, 0)


    def test_parse_cds_seqfeature_4(self):
        """Verify CDS features is parsed with no translation table."""
        self.qualifier_dict.pop("transl_table")
        cds_ftr = flat_files.parse_cds_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertEqual(cds_ftr.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(cds_ftr.translation_table, 0)


    def test_parse_cds_seqfeature_5(self):
        """Verify CDS features is parsed with no product."""
        self.qualifier_dict.pop("product")
        cds_ftr = flat_files.parse_cds_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertEqual(cds_ftr.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(cds_ftr.raw_product, "")
        with self.subTest():
            self.assertEqual(cds_ftr.product, "")


    def test_parse_cds_seqfeature_6(self):
        """Verify CDS features is parsed with no function."""
        self.qualifier_dict.pop("function")
        cds_ftr = flat_files.parse_cds_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertEqual(cds_ftr.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(cds_ftr.raw_function, "")
        with self.subTest():
            self.assertEqual(cds_ftr.function, "")


    def test_parse_cds_seqfeature_7(self):
        """Verify CDS features is parsed with no note."""
        self.qualifier_dict.pop("note")
        cds_ftr = flat_files.parse_cds_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertEqual(cds_ftr.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(cds_ftr.raw_note, "")
        with self.subTest():
            self.assertEqual(cds_ftr.note, "")


    def test_parse_cds_seqfeature_8(self):
        """Verify CDS features is parsed with no gene."""
        self.qualifier_dict.pop("gene")
        cds_ftr = flat_files.parse_cds_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertEqual(cds_ftr.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(cds_ftr._locus_tag_num, "1")
        with self.subTest():
            self.assertEqual(cds_ftr.gene, "")
        with self.subTest():
            self.assertEqual(cds_ftr.name, "1")


    def test_parse_cds_seqfeature_9(self):
        """Verify CDS features is parsed with 3-part compound location."""
        self.seqfeature = SeqFeature(CompoundLocation([
                        FeatureLocation(
                            ExactPosition(2),
                            ExactPosition(10),
                            strand=1),
                        FeatureLocation(
                            ExactPosition(8),
                            ExactPosition(20),
                            strand=1),
                        FeatureLocation(
                            ExactPosition(30),
                            ExactPosition(50),
                            strand=1)],
                        'join'),
                        type='CDS',
                        location_operator='join',
                        qualifiers = self.qualifier_dict)
        cds_ftr = flat_files.parse_cds_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertEqual(cds_ftr.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(cds_ftr.orientation, "F")
        with self.subTest():
            self.assertEqual(cds_ftr.start, -1)
        with self.subTest():
            self.assertEqual(cds_ftr.stop, -1)
        with self.subTest():
            self.assertEqual(cds_ftr.parts, 3)
        with self.subTest():
            self.assertEqual(cds_ftr.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(cds_ftr.length, 0)


    def test_parse_cds_seqfeature_10(self):
        """Verify CDS features is parsed with fuzzy coordinates."""
        self.seqfeature = SeqFeature(
                            FeatureLocation(
                                BeforePosition(2),
                                ExactPosition(10)),
                            type="CDS",
                            strand=1,
                            qualifiers=self.qualifier_dict)
        cds_ftr = flat_files.parse_cds_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertEqual(cds_ftr.locus_tag, "SEA_L5_1")
        with self.subTest():
            self.assertEqual(cds_ftr.orientation, "F")
        with self.subTest():
            self.assertEqual(cds_ftr.start, -1)
        with self.subTest():
            self.assertEqual(cds_ftr.stop, 10)
        with self.subTest():
            self.assertEqual(cds_ftr.parts, 1)
        with self.subTest():
            self.assertEqual(cds_ftr.coordinate_format, "0_half_open")
        with self.subTest():
            self.assertEqual(cds_ftr.length, 11)




class TestFlatFileFunctions3(unittest.TestCase):


    def setUp(self):
        self.seqfeature1 = SeqFeature(FeatureLocation(
                                      ExactPosition(2), ExactPosition(10)),
                                      type = "CDS",
                                      strand = 1)

        self.seqfeature2 = SeqFeature(FeatureLocation(
                                      ExactPosition(5000), ExactPosition(6000)),
                                      type = "tRNA",
                                      strand = 1)

        self.seqfeature3 = SeqFeature(FeatureLocation(
                                      ExactPosition(1), ExactPosition(11000)),
                                      type = "source",
                                      strand = 1)

        self.seqfeature7 = SeqFeature(FeatureLocation(
                                      ExactPosition(1), ExactPosition(9000)),
                                      type = "source",
                                      strand = 1)

        # Wrap-around feature, directly copied from
        # Biopython-parsed ET08 flat file.
        self.seqfeature4 = SeqFeature(CompoundLocation([
                                    FeatureLocation(
                                        ExactPosition(154873),
                                        ExactPosition(155445),
                                        strand=1),
                                    FeatureLocation(
                                        ExactPosition(0),
                                        ExactPosition(4),
                                        strand=1)],
                                    'join'),
                                    type='CDS',
                                    location_operator='join')

        self.seqfeature5 = SeqFeature(FeatureLocation(
                                      ExactPosition(9), ExactPosition(50)),
                                      type = "CDS",
                                      strand = -1)

        self.seqfeature6 = SeqFeature(FeatureLocation(
                                      ExactPosition(9), ExactPosition(30)),
                                      type = "CDS",
                                      strand = 1)

        self.feature_list = [self.seqfeature1,
                             self.seqfeature2,
                             self.seqfeature3,
                             self.seqfeature4,
                             self.seqfeature5,
                             self.seqfeature6,
                             self.seqfeature7]

        self.reference1 = Reference()
        self.reference1.authors = "Jane"

        self.reference2 = Reference()
        self.reference2.authors = "Doe"

        self.reference3 = Reference()
        self.reference3.authors = "Smith"

        self.refs_list = [self.reference1,
                          self.reference2,
                          self.reference3]

        self.description = "Mycobacterium phage L5 complete genome"
        self.organism = "Gordonia phage KatherineG"
        self.source = "Streptomyces phage phiC31"

        self.date = "23-JAN-2014"

        self.annotation_dict = {"accessions": [" ABC123.1 "],
                                "source": self.source,
                                "organism": self.organism,
                                "references": self.refs_list,
                                "date": self.date}

        self.record = SeqRecord(seq = Seq("atgc"),
                                id = "OPQ123.1",
                                name = "XYZ123",
                                annotations = self.annotation_dict,
                                description = self.description,
                                features = self.feature_list
                                )

        self.filepath = pathlib.Path("/path/to/file/Phage_ZZZ.gb")
        self.exp_date = datetime.strptime(self.date,'%d-%b-%Y')

    def test_parse_genome_data_1(self):
        """Verify retrieved flat file data is parsed correctly."""
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                           gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.name, "KatherineG")
        with self.subTest():
            self.assertEqual(gnm.organism, self.organism)
        with self.subTest():
            self.assertEqual(gnm._organism_name, "KatherineG")
        with self.subTest():
            self.assertEqual(gnm._organism_host_genus, "Gordonia")
        with self.subTest():
            self.assertEqual(gnm.accession, "ABC123")
        with self.subTest():
            self.assertEqual(gnm.description, self.description)
        with self.subTest():
            self.assertEqual(gnm._description_name, "L5")
        with self.subTest():
            self.assertEqual(gnm._description_host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(gnm.source, self.source)
        with self.subTest():
            self.assertEqual(gnm._source_name, "phiC31")
        with self.subTest():
            self.assertEqual(gnm._source_host_genus, "Streptomyces")
        with self.subTest():
            self.assertEqual(gnm.authors, "Jane;Doe;Smith")
        with self.subTest():
            self.assertEqual(gnm.seq, "ATGC")
        with self.subTest():
            self.assertEqual(gnm.length, 4)
        with self.subTest():
            self.assertEqual(gnm.gc, 50.00)
        with self.subTest():
            self.assertEqual(gnm.date, self.exp_date)
        with self.subTest():
            self.assertEqual(len(gnm.cds_features), 4)
        with self.subTest():
            self.assertEqual(len(gnm.source_features), 2)
        with self.subTest():
            self.assertEqual(len(gnm.trna_features), 1)
        with self.subTest():
            self.assertEqual(gnm._cds_features_tally, 4)
        with self.subTest():
            self.assertEqual(gnm._source_features_tally, 2)
        with self.subTest():
            self.assertEqual(gnm._trna_features_tally, 1)
        with self.subTest():
            self.assertEqual(gnm.translation_table, 11)
        with self.subTest():
            self.assertEqual(gnm.id,"KatherineG")
        with self.subTest():
            self.assertEqual(gnm.host_genus, "Gordonia")
        with self.subTest():
            self.assertEqual(gnm.type, "flat_file")
        with self.subTest():
            self.assertEqual(gnm.cds_features[0].genome_id, "KatherineG")
        with self.subTest():
            self.assertEqual(gnm.cds_features[0].id, "KatherineG_CDS_1")
        with self.subTest():
            self.assertEqual(gnm.cds_features[1].id, "KatherineG_CDS_4")
        with self.subTest():
            self.assertEqual(gnm.cds_features[2].id, "KatherineG_CDS_3")
        with self.subTest():
            self.assertEqual(gnm.cds_features[3].id, "KatherineG_CDS_2")

        with self.subTest():
            self.assertEqual(gnm.cds_features[0].start, 2)
        with self.subTest():
            self.assertEqual(gnm.cds_features[0].stop, 10)
        with self.subTest():
            self.assertEqual(gnm.cds_features[0].genome_length, 4)

        with self.subTest():
            self.assertEqual(gnm.cds_features[1].start, 154873)
        with self.subTest():
            self.assertEqual(gnm.cds_features[1].stop, 4)
        with self.subTest():
            self.assertEqual(gnm.cds_features[1].genome_length, 4)

        with self.subTest():
            self.assertEqual(gnm.cds_features[2].start, 9)
        with self.subTest():
            self.assertEqual(gnm.cds_features[2].stop, 50)

        with self.subTest():
            self.assertEqual(gnm.cds_features[3].start, 9)
        with self.subTest():
            self.assertEqual(gnm.cds_features[3].stop, 30)


        with self.subTest():
            self.assertEqual(gnm.source_features[0].genome_id, "KatherineG")
        with self.subTest():
            self.assertEqual(gnm.source_features[0].start, 1)
        with self.subTest():
            self.assertEqual(gnm.source_features[0].stop, 11000)
        with self.subTest():
            self.assertEqual(gnm.source_features[0].id, "KatherineG_SRC_2")
        with self.subTest():
            self.assertEqual(gnm.source_features[1].id, "KatherineG_SRC_1")


    def test_parse_genome_data_2(self):
        """Verify retrieved flat file data is parsed correctly with no
        record name."""
        self.record = SeqRecord(seq = Seq("atgc"),
                                id = "OPQ123.1",
                                annotations = self.annotation_dict,
                                description = self.description,
                                features = self.feature_list
                                )
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                           gnm_type="", genome_id_field="")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.id, "OPQ123.1")
        with self.subTest():
            self.assertEqual(gnm.name, "")
        with self.subTest():
            self.assertEqual(gnm.type, "")

    def test_parse_genome_data_3(self):
        """Verify retrieved flat file data is parsed correctly with no
        record organism."""
        self.annotation_dict.pop("organism")
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                              gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.organism, "")
        with self.subTest():
            self.assertEqual(gnm._organism_name, "")
        with self.subTest():
            self.assertEqual(gnm._organism_host_genus, "")


    def test_parse_genome_data_4(self):
        """Verify retrieved flat file data is parsed correctly with no
        record id."""
        self.record = SeqRecord(seq = Seq("atgc"),
                                name = "XYZ123",
                                annotations = self.annotation_dict,
                                description = self.description,
                                features = self.feature_list,
                                )
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                    gnm_type="flat_file", genome_id_field="")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.id, "")
        with self.subTest():
            self.assertEqual(gnm.name, "XYZ123")


    def test_parse_genome_data_5(self):
        """Verify retrieved flat file data is parsed correctly with no
        accession."""
        self.annotation_dict.pop("accessions")
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                                gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.accession, "")


    def test_parse_genome_data_6(self):
        """Verify retrieved flat file data is parsed correctly with more
        than one accession."""
        self.annotation_dict["accessions"] = [" ABC123.1 ", "TUV456"]
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                                gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.accession, "ABC123")


    def test_parse_genome_data_7(self):
        """Verify retrieved flat file data is parsed correctly with no
        record description."""
        self.record = SeqRecord(seq = Seq("atgc"),
                                id = "OPQ123.1",
                                name = "XYZ123",
                                annotations = self.annotation_dict,
                                features = self.feature_list
                                )

        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                                gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.description, "")
        with self.subTest():
            self.assertEqual(gnm._description_name, "")
        with self.subTest():
            self.assertEqual(gnm._description_host_genus, "")


    def test_parse_genome_data_8(self):
        """Verify retrieved flat file data is parsed correctly with no
        record source."""
        self.annotation_dict.pop("source")
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                                gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.source, "")
        with self.subTest():
            self.assertEqual(gnm._source_name, "")
        with self.subTest():
            self.assertEqual(gnm._source_host_genus, "")


    def test_parse_genome_data_9(self):
        """Verify retrieved flat file data is parsed correctly with no
        references."""
        self.annotation_dict.pop("references")
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                                gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.authors, "")


    def test_parse_genome_data_10(self):
        """Verify retrieved flat file data is parsed correctly with
        references that contain no authors."""
        self.annotation_dict["references"] = [Reference(), Reference()]
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                                gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.authors, "")


    def test_parse_genome_data_11(self):
        """Verify retrieved flat file data is parsed correctly with an
        empty sequence."""
        print(self.record.seq)
        self.record.seq = Seq("")
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                                gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.seq, "")
        with self.subTest():
            self.assertEqual(gnm.length, 0)
        with self.subTest():
            self.assertEqual(gnm.gc, -1)


    def test_parse_genome_data_12(self):
        """Verify retrieved flat file data is parsed correctly with no
        CDS features."""
        self.feature_list = [self.seqfeature2,
                             self.seqfeature3,
                             self.seqfeature7]
        self.record.features = self.feature_list
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                                gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(len(gnm.cds_features), 0)
        with self.subTest():
            self.assertEqual(len(gnm.source_features), 2)
        with self.subTest():
            self.assertEqual(len(gnm.trna_features), 1)
        with self.subTest():
            self.assertEqual(gnm._cds_features_tally, 0)


    def test_parse_genome_data_13(self):
        """Verify retrieved flat file data is parsed correctly with no
        source features."""
        self.feature_list = [self.seqfeature1,
                             self.seqfeature2,
                             self.seqfeature4,
                             self.seqfeature5,
                             self.seqfeature6]
        self.record.features = self.feature_list
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                                gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(len(gnm.cds_features), 4)
        with self.subTest():
            self.assertEqual(len(gnm.source_features), 0)
        with self.subTest():
            self.assertEqual(len(gnm.trna_features), 1)
        with self.subTest():
            self.assertEqual(gnm._source_features_tally, 0)


    def test_parse_genome_data_14(self):
        """Verify retrieved flat file data is parsed correctly with no
        tRNA features."""
        self.feature_list = [self.seqfeature1,
                             self.seqfeature3,
                             self.seqfeature4,
                             self.seqfeature5,
                             self.seqfeature6,
                             self.seqfeature7]
        self.record.features = self.feature_list
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                                gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(len(gnm.cds_features), 4)
        with self.subTest():
            self.assertEqual(len(gnm.source_features), 2)
        with self.subTest():
            self.assertEqual(len(gnm.trna_features), 0)
        with self.subTest():
            self.assertEqual(gnm._trna_features_tally, 0)


    def test_parse_genome_data_15(self):
        """Verify retrieved flat file data is parsed correctly with no
        features."""
        self.record.features = []
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                                gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(len(gnm.cds_features), 0)
        with self.subTest():
            self.assertEqual(len(gnm.source_features), 0)
        with self.subTest():
            self.assertEqual(len(gnm.trna_features), 0)
        with self.subTest():
            self.assertEqual(gnm._cds_features_tally, 0)
        with self.subTest():
            self.assertEqual(gnm._source_features_tally, 0)
        with self.subTest():
            self.assertEqual(gnm._trna_features_tally, 0)


    def test_parse_genome_data_16(self):
        """Verify retrieved flat file data is parsed correctly with no
        filepath provided."""
        gnm = flat_files.parse_genome_data(self.record, gnm_type="flat_file")
        with self.subTest():
            self.assertEqual(gnm.filename, "")
        with self.subTest():
            self.assertEqual(gnm.name, "KatherineG")


    def test_parse_genome_data_17(self):
        """Verify retrieved flat file data is parsed correctly with no
        date provided."""
        self.annotation_dict.pop("date")
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                                gnm_type="flat_file")
        self.exp_date = constants.EMPTY_DATE
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.date, self.exp_date)


    def test_parse_genome_data_18(self):
        """Verify retrieved flat file data is parsed correctly with
        modified translation table."""
        gnm = flat_files.parse_genome_data(self.record, self.filepath,
                                                gnm_type="flat_file",
                                                translation_table=1)
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.translation_table, 1)


    def test_parse_genome_data_19(self):
        """Verify retrieved flat file data is parsed correctly with
        id field and host genus field specified as non-standard fields."""
        gnm = flat_files.parse_genome_data(
                            self.record,
                            self.filepath,
                            gnm_type="flat_file",
                            genome_id_field = "_description_name",
                            host_genus_field = "_description_host_genus")
        with self.subTest():
            self.assertEqual(gnm.filename, "Phage_ZZZ")
        with self.subTest():
            self.assertEqual(gnm.id, "L5")
        with self.subTest():
            self.assertEqual(gnm.host_genus, "Mycobacterium")




class TestFlatFileFunctions4(unittest.TestCase):

    def setUp(self):
        self.string1 = "Mycobacterium phage Trixie"
        self.string2 = "Mycobacterium smegmatis"
        self.string3 = "Gordonia terrae"
        self.qualifier_dict = {"organism": [self.string1],
                               "host": [self.string2],
                               "lab_host": [self.string3]}

        self.seqfeature = SeqFeature(FeatureLocation(
                                ExactPosition(2), ExactPosition(10)),
                                type="source", strand=1,
                                qualifiers=self.qualifier_dict)

    def test_parse_source_seqfeature_1(self):
        """Verify source feature is parsed."""

        src_ftr = flat_files.parse_source_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertIsInstance(src_ftr.seqfeature, SeqFeature)
        with self.subTest():
            self.assertEqual(src_ftr.start, 2)
        with self.subTest():
            self.assertEqual(src_ftr.stop, 10)
        with self.subTest():
            self.assertEqual(src_ftr.organism, self.string1)
        with self.subTest():
            self.assertEqual(src_ftr.host, self.string2)
        with self.subTest():
            self.assertEqual(src_ftr.lab_host, self.string3)
        with self.subTest():
            self.assertEqual(src_ftr._organism_name, "Trixie")
        with self.subTest():
            self.assertEqual(src_ftr._organism_host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(src_ftr._host_host_genus, "Mycobacterium")
        with self.subTest():
            self.assertEqual(src_ftr._lab_host_host_genus, "Gordonia")

    def test_parse_source_seqfeature_2(self):
        """Verify source feature is parsed with no organism qualifier."""

        self.seqfeature.qualifiers.pop("organism")
        src_ftr = flat_files.parse_source_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertEqual(src_ftr.organism, "")
        with self.subTest():
            self.assertEqual(src_ftr._organism_name, "")
        with self.subTest():
            self.assertEqual(src_ftr._organism_host_genus, "")

    def test_parse_source_seqfeature_3(self):
        """Verify source feature is parsed with no host qualifier."""

        self.seqfeature.qualifiers.pop("host")
        src_ftr = flat_files.parse_source_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertEqual(src_ftr.host, "")
        with self.subTest():
            self.assertEqual(src_ftr._host_host_genus, "")

    def test_parse_source_seqfeature_4(self):
        """Verify source feature is parsed with no lab_host qualifier."""

        self.seqfeature.qualifiers.pop("lab_host")
        src_ftr = flat_files.parse_source_seqfeature(self.seqfeature)
        with self.subTest():
            self.assertEqual(src_ftr.lab_host, "")
        with self.subTest():
            self.assertEqual(src_ftr._lab_host_host_genus, "")




if __name__ == '__main__':
    unittest.main()
