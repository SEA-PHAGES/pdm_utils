"""Integration tests for the entire import pipeline."""

import csv
# import time
from datetime import datetime
import unittest
import pymysql
import shutil
from unittest.mock import patch
from pathlib import Path
from pdm_utils import run
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqFeature import ExactPosition, Reference

from pdm_utils.classes import mysqlconnectionhandler as mch
from collections import OrderedDict



# import argparse
# from pdm_utils.pipelines.db_import import import_genome
# from pdm_utils.constants import constants
# from pdm_utils.functions import basic
# from pdm_utils.classes import bundle, genome, ticket, eval
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
# from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
# from Bio.SeqFeature import ExactPosition, Reference
# from unittest.mock import patch
# import getpass

current_date = datetime.today().replace(hour=0, minute=0,
                                        second=0, microsecond=0)

# The following integration tests user the 'pdm_anon' MySQL user.
# It is expected that this user has all privileges for 'test_db' database.
user = "pdm_anon"
pwd = "pdm_anon"
db = "test_db"
unittest_file = Path(__file__)
unittest_dir = unittest_file.parent
schema_file = "test_schema5.sql"
schema_filepath = Path(unittest_dir, "test_files/", schema_file)

# The primary genome used for the following integration tests is Alice.
# Alice is a large genome with:
# 1. a gene that wraps around the genome termini.
# 2. a tail assembly chaperone gene that is a compound feature.
# 3. tRNA and tmRNA genes.
base_flat_file = "test_flat_file_10.gb"
base_flat_file_path = Path(unittest_dir, "test_files/", base_flat_file)


def create_new_db(schema_file, db, user, pwd):
    """Creates a new, empty database."""

    connection = pymysql.connect(host = "localhost",
                                 user = user,
                                 password = pwd,
                                 cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()

    # First, test if a test database already exists within mysql.
    # If there is, delete it so that a fresh test database is installed.
    sql = ("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA "
          f"WHERE SCHEMA_NAME = '{db}'")
    cur.execute(sql)
    result = cur.fetchall()
    if len(result) != 0:
        cur.execute(f"DROP DATABASE {db}")
        connection.commit()

    # Next, create the database within mysql.
    cur.execute(f"CREATE DATABASE {db}")
    connection.commit()
    connection.close()

    # Now import the empty schema from file.
    # Seems like pymysql has trouble with this step, so use subprocess.
    handle = open(schema_filepath, "r")
    command_string = f"mysql -u {user} -p{pwd} {db}"
    command_list = command_string.split(" ")
    proc = subprocess.check_call(command_list, stdin = handle)
    handle.close()




def remove_db(db, user, pwd):
    """Remove the MySQL database created for the test."""
    connection = pymysql.connect(host="localhost",
                                 user=user,
                                 password=pwd,
                                 cursorclass=pymysql.cursors.DictCursor)
    cur = connection.cursor()
    cur.execute(f"DROP DATABASE {db}")
    connection.commit()
    connection.close()


def insert_data_into_phage_table(db, user, pwd, data_dict):
    """Insert data into the phage table."""
    connection = pymysql.connect(host = "localhost",
                                 user = user,
                                 password = pwd,
                                 database = db,
                                 cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()
    sql = (
        "INSERT INTO phage "
        "(PhageID, Accession, Name, "
        "HostStrain, Sequence, SequenceLength, GC, status, "
        "DateLastModified, RetrieveRecord, AnnotationAuthor, "
        "Cluster, Cluster2, Subcluster2) "
        "VALUES ("
        f"'{data_dict['PhageID']}', '{data_dict['Accession']}', "
        f"'{data_dict['Name']}', '{data_dict['HostStrain']}', "
        f"'{data_dict['Sequence']}', {data_dict['SequenceLength']}, "
        f"{data_dict['GC']}, '{data_dict['status']}', "
        f"'{data_dict['DateLastModified']}', "
        f"{data_dict['RetrieveRecord']}, "
        f"{data_dict['AnnotationAuthor']}, '{data_dict['Cluster']}', "
        f"'{data_dict['Cluster2']}', '{data_dict['Subcluster2']}');"
        )
    cur.execute(sql)
    connection.commit()
    connection.close()


def insert_data_into_gene_table(db, user, pwd, data_dict):
    """Insert data into the gene table."""
    connection = pymysql.connect(host = "localhost",
                                 user = user,
                                 password = pwd,
                                 database = db,
                                 cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()
    sql = (
        "INSERT INTO gene "
        "(GeneID, PhageID, Start, Stop, Length, Name, "
        "translation, Orientation, Notes, LocusTag) "
        "VALUES ("
        f"'{data_dict['GeneID']}', '{data_dict['PhageID']}', "
        f"{data_dict['Start']}, {data_dict['Stop']}, "
        f"{data_dict['Length']}, '{data_dict['Name']}', "
        f"'{data_dict['translation']}', '{data_dict['Orientation']}', "
        f"'{data_dict['Notes']}', '{data_dict['LocusTag']}');"
        )
    cur.execute(sql)
    connection.commit()
    connection.close()


phage_table_query = (
    "SELECT "
    "PhageID, Accession, Name, "
    "HostStrain, Sequence, SequenceLength, GC, status, "
    "DateLastModified, RetrieveRecord, AnnotationAuthor, "
    "Cluster, Cluster2, Subcluster2 "
    "FROM phage;"
    )

gene_table_query = (
    "SELECT "
    "GeneID, PhageID, Start, Stop, Length, Name, "
    "translation, Orientation, Notes, LocusTag "
    "FROM gene;"
)



def get_sql_data(db, user, pwd, query):
    """Get data from the database."""
    connection = pymysql.connect(host = "localhost",
                                 user = user,
                                 password = pwd,
                                 database = db,
                                 cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()
    cur.execute(query)
    results = cur.fetchall()
    cur.close()
    connection.close()
    return results


def process_phage_table_data(list_of_sql_results):
    """."""
    x = 0
    while x < len(list_of_sql_results):
        data_dict = list_of_sql_results[x]
        # print(data_dict.keys())
        # input("pause")
        data_dict["Sequence"] = data_dict["Sequence"].decode("utf-8")
        data_dict["SequenceLength"] = int(data_dict["SequenceLength"])
        # data_dict["Notes"] = data_dict["Notes"].decode("utf-8")
        data_dict["GC"] = float(data_dict["GC"])
        x += 1


def process_gene_table_data(list_of_sql_results):
    """."""
    x = 0
    while x < len(list_of_sql_results):
        data_dict = list_of_sql_results[x]
        data_dict["Notes"] = data_dict["Notes"].decode("utf-8")
        data_dict["Length"] = int(data_dict["Length"])
        data_dict["Start"] = int(data_dict["Start"])
        data_dict["Stop"] = int(data_dict["Stop"])
        x += 1

def filter_genome_data(list_of_sql_results, phage_id):
    """."""
    x = 0
    output_dict = {}
    while x < len(list_of_sql_results):
        if list_of_sql_results[x]["PhageID"] == phage_id:
            output_dict = list_of_sql_results[x]
        else:
            pass
        x += 1
    return output_dict



def filter_gene_data(list_of_sql_results, coordinates):
    """."""
    x = 0
    output_dict = {}
    while x < len(list_of_sql_results):
        if (list_of_sql_results[x]["Start"] == coordinates[0] and
                list_of_sql_results[x]["Stop"] == coordinates[1]):
            output_dict = list_of_sql_results[x]
        else:
            pass
        x += 1

    return output_dict



def compare_data(ref_dict, query_dict):
    """Ref dictionary = the expected data.
    Query dictionary = the data retrieved from the database."""
    x = 0
    errors = 0
    for key in ref_dict.keys():
        exp_value = ref_dict[key]
        query_value = query_dict[key]
        if exp_value != query_value:
            errors += 1
            print(key)
            print(exp_value)
            print(query_value)
            input("check value")
    return errors

def get_seq(filepath):
    """Get genome sequence from a flat file so that it can be added
    to a MySQL database for testing."""
    seqrecord = SeqIO.read(filepath, "genbank")
    return seqrecord.seq


l5_phage_table_data = {
    "PhageID": "L5",
    "Accession": "ABC123",
    "Name": "L5",
    "HostStrain": "Mycobacterium",
    "Sequence": "AAAAAAAAAA",
    "SequenceLength": 10,
    "GC": 0,
    "status": "final",
    "DateLastModified": datetime.strptime('1/1/0001', '%m/%d/%Y'),
    "RetrieveRecord": 1,
    "AnnotationAuthor": 1,
    "Cluster": "A2",
    "Cluster2": "A",
    "Subcluster2": "A2"
    }



trixie_phage_table_data = {
    "PhageID": "Trixie",
    "Accession": "BCD456",
    "Name": "Trixie",
    "HostStrain": "Gordonia",
    "Sequence": "GGGGGGGGGGGGGGGGGGGG",
    "SequenceLength": 20,
    "GC": 1,
    "status": "final",
    "DateLastModified": datetime.strptime('1/1/2000', '%m/%d/%Y'),
    "RetrieveRecord": 1,
    "AnnotationAuthor": 1,
    "Cluster": "A3",
    "Cluster2": "A",
    "Subcluster2": "A3"
    }

redrock_phage_table_data = {
    "PhageID": "RedRock",
    "Accession": "BCD456",
    "Name": "RedRock_Draft",
    "HostStrain": "Arthrobacter",
    "Sequence": "CCCCCCCCCCCCCCC",
    "SequenceLength": 15,
    "GC": 1,
    "status": "draft",
    "DateLastModified": datetime.strptime('1/1/2010', '%m/%d/%Y'),
    "RetrieveRecord": 1,
    "AnnotationAuthor": 1,
    "Cluster": "B1",
    "Cluster2": "B",
    "Subcluster2": "B1"
    }


d29_phage_table_data = {
    "PhageID": "D29",
    "Accession": "XYZ123",
    "Name": "D29",
    "HostStrain": "Microbacterium",
    "Sequence": "ATGCATGCATGCATGC",
    "SequenceLength": 16,
    "GC": 0.5,
    "status": "unknown",
    "DateLastModified": datetime.strptime('1/1/0001', '%m/%d/%Y'),
    "RetrieveRecord": 0,
    "AnnotationAuthor": 0,
    "Cluster": "C1",
    "Cluster2": "C",
    "Subcluster2": "C1"
    }






l5_gene_table_data_1 = {
    "GeneID": "L5_0001",
    "PhageID": "L5",
    "Start": 100,
    "Stop": 1100,
    "Length": 1000,
    "Name": "1",
    "translation": "ACTGC",
    "Orientation": "F",
    "Notes": "int",
    "LocusTag": "SEA_L5_0001"
    }



trixie_gene_table_data_1 = {
    "GeneID": "TRIXIE_0001",
    "PhageID": "Trixie",
    "Start": 100,
    "Stop": 1100,
    "Length": 1000,
    "Name": "1",
    "translation": "ACTGC",
    "Orientation": "F",
    "Notes": "int",
    "LocusTag": "SEA_TRIXIE_0001"
    }




def create_min_tkt_dict(old_tkt):

    min_keys = set(["id", "type", "phage_id"])
    new_tkt = {}
    for key in min_keys:
        new_tkt[key] = old_tkt[key]
    return new_tkt

def set_data(old_tkt, keys, value):
    """Set selected keys to certain values."""
    new_tkt = old_tkt.copy()
    for key in keys:
        new_tkt[key] = value
    return new_tkt



l5_ticket_data_complete = {
    "id": 1,
    "type": "add",
    "phage_id": "L5",
    "host_genus": "Mycobacterium",
    "cluster": "A",
    "subcluster": "A2",
    "accession": "ABC123",
    "description_field": "product",
    "annotation_status": "draft",
    "annotation_author": 1,
    "retrieve_record": 1,
    "run_mode": "phagesdb",
    }



alice_ticket_data_complete = {
    "id": 1,
    "type": "add",
    "phage_id": "Alice",
    "host_genus": "Mycobacterium",
    "cluster": "C",
    "subcluster": "C1",
    "accession": "none",
    "description_field": "product",
    "annotation_status": "draft",
    "annotation_author": 1,
    "retrieve_record": 1,
    "run_mode": "pecaan",
    }





l5_ticket_data_min = create_min_tkt_dict(l5_ticket_data_complete)

valid_retrieve_keys = set(["host_genus", "cluster", "subcluster", "accession"])
l5_ticket_data_retrieve = set_data(l5_ticket_data_complete, valid_retrieve_keys, "retrieve")

valid_retain_keys = set(["host_genus", "cluster", "subcluster",
                         "accession", "annotation_author", "retrieve_record"])






def create_import_table(list_of_data_dicts, file_path):
    """Create an import table."""

    headers = set()
    for dict in list_of_data_dicts:
        headers = headers | dict.keys()

    headers_dict = {}
    for header in headers:
        headers_dict[header] = header
    with open(file_path, "w") as file_handle:
        file_writer = csv.DictWriter(file_handle, headers)
        file_writer.writerow(headers_dict)
        for data_dict in list_of_data_dicts:
            file_writer.writerow(data_dict)





def clear_descriptions(record):
    """Remove descriptions from CDS features."""
    x = 0
    while x < len(record.features):
        feature = record.features[x]
        if feature.type == "CDS":
            if "product" in feature.qualifiers.keys():
                feature.qualifiers["product"] = [""]
            if "function" in feature.qualifiers.keys():
                feature.qualifiers["function"] = [""]
            if "note" in feature.qualifiers.keys():
                feature.qualifiers["note"] = [""]
        x += 1





class TestImportGenomeMain1(unittest.TestCase):


    def setUp(self):

        create_new_db(schema_file, db, user, pwd)


        self.base_dir = Path(unittest_dir, "test_wd/test_import")
        self.base_dir.mkdir()

        self.import_table = Path(self.base_dir, "import_table.csv")
        self.genome_folder = Path(self.base_dir, "genome_folder")
        self.genome_folder.mkdir()

        self.output_folder = Path(self.base_dir, "output_folder")
        self.output_folder.mkdir()

        self.log_file = Path(self.output_folder, "import_log.txt")

        self.test_flat_file1 = Path(unittest_dir,
                                    "test_files/test_flat_file_1.gb")

        self.base_flat_file_path = base_flat_file_path

        self.test_import_table1 = Path(unittest_dir,
                                    "test_files/test_import_table_1.csv")


        self.sql_handle = mch.MySQLConnectionHandler()
        self.sql_handle.database = db
        self.sql_handle.username = user
        self.sql_handle.password = pwd


        # Construct minimal Alice genome
        # CDS 252 is a wrap-around compound gene.
        # CDS 124 is a compound feature.
        # CDS 193 is a normal CDS.
        # CDS 139 is a normal CDS.

        # self.alice_cds_data_in_db = {
        #     "GeneID": "",
        #     "PhageID": "",
        #     "Start": "",
        #     "Stop": "",
        #     "Length": "",
        #     "Name": "",
        #     "translation": "",
        #     "Orientation": "",
        #     "Notes": "",
        #     "LocusTag": ""
        #     }

        self.alice_cds_252_translation = (
            "MFDNHPLMSTLITSVRDLNRLANGVIIRNRCEKCDTSAGPDQGV"
            "HFVKTPLGWLYSDPKGKPTTWELFPSDAIHLPVQVIHRVSESEDVSEITENSPHTRKD"
            "DSEATEGAAPSFEPLYPSSHKTERFTVQDTPDGLGAYVFDFGGDAFGAQTCADALAKV"
            "TGKTWYVMHKTVVENTGIYSSMVRHEEESSS")
        self.alice_cds_252_qualifier_dict = OrderedDict(
            [("gene", ["252"]),
             ("locus_tag", ["ALICE_252"]),
             ("note", ["gp252"]),
             ("codon_start", ["1"]),
             ("transl_table", ["11"]),
             ("product", ["hypothetical protein"]),
             ("protein_id", ["AEJ94474.1"]),
             ("translation", [self.alice_cds_252_translation])])
        self.alice_cds_252 = SeqFeature(
                                CompoundLocation(
                                    [FeatureLocation(
                                        ExactPosition(152829),
                                        ExactPosition(153401),
                                        strand=1),
                                    FeatureLocation(
                                        ExactPosition(0),
                                        ExactPosition(4),
                                        strand=1)],
                                    "join"),
                                type="CDS",
                                location_operator="join",
                                qualifiers=self.alice_cds_252_qualifier_dict)

        self.alice_cds_252_coords = (152829, 4)
        self.alice_cds_252_data_in_db = {
            "GeneID": "Alice_CDS_4",
            "PhageID": "Alice",
            "Start": 152829,
            "Stop": 4,
            "Length": 191,
            "Name": "252",
            "translation": self.alice_cds_252_translation,
            "Orientation": "F",
            "Notes": "",
            "LocusTag": "ALICE_252"
            }

        self.alice_cds_124_translation = (
            "MKDEMTTSDVPADPAIDPDLAPPEPRRVVGELVETEPQEHEDPEVTELTDEERSSFVSLLT"
            "CGKHSKKITVMGHPVVIQTLKTGDEMRVGLFTKKYLESQMGFQRAYQVAVCAAGIREIQGK"
            "PLFRELREVTDEDEIFDKNVEAVMELYPIVITQIYQAIMDLEREYAQLAVKLGKTVRLDAS"
            "TELEIRLAYKQGLLTQPSLNRYQRWALRYAIFMDRRLQLQDTEDMLQRQTWYLEPKRYHDL"
            "FLAGAFEPEPIAVAGRDMEEVVDDLDEVDAYFARLEGSQSMSGAQLFAALDEPDEEGWM"
            )
        self.alice_cds_124_qualifier_dict = OrderedDict(
            [("gene", ["124"]),
             ("locus_tag", ["ALICE_124"]),
             ("ribosomal_slippage", [""]),
             ("note", ["gp124"]),
             ("codon_start", ["1"]),
             ("transl_table", ["11"]),
             ("product", ["tail assembly chaperone"]),
             ("protein_id", ["AEJ94379.1"]),
             ("translation", [self.alice_cds_124_translation])])
        self.alice_cds_124 = SeqFeature(
                                CompoundLocation(
                                    [FeatureLocation(
                                        ExactPosition(70374),
                                        ExactPosition(70902),
                                        strand=1),
                                    FeatureLocation(
                                        ExactPosition(70901),
                                        ExactPosition(71285),
                                        strand=1)],
                                    "join"),
                                type="CDS",
                                location_operator="join",
                                qualifiers=self.alice_cds_124_qualifier_dict)

        self.alice_cds_124_coords = (70374, 71285)
        self.alice_cds_124_data_in_db = {
            "GeneID": "Alice_CDS_1",
            "PhageID": "Alice",
            "Start": 70374,
            "Stop": 71285,
            "Length": 303,
            "Name": "124",
            "translation": self.alice_cds_124_translation,
            "Orientation": "F",
            "Notes": "tail assembly chaperone",
            "LocusTag": "ALICE_124"
            }



        self.alice_cds_139_translation = (
            "MKKLVTGGLVAAGITLGLISAPSASAEYMTICPSQVSAVVTANTSCGFADNV"
            "FRGFYRQSGWDPLAYSPATGKVYRMHCAPATTTNWGEAKRCWGVGYGGDLLVVYID"
            )
        self.alice_cds_139_qualifier_dict = OrderedDict(
            [("gene", ["139"]),
             ("locus_tag", ["ALICE_139"]),
             ("note", ["gp139"]),
             ("codon_start", ["1"]),
             ("transl_table", ["11"]),
             ("product", ["hypothetical protein"]),
             ("protein_id", ["AEJ94477.1"]),
             ("translation", [self.alice_cds_139_translation])])
        self.alice_cds_139 = SeqFeature(
                                FeatureLocation(
                                    ExactPosition(88120),
                                    ExactPosition(88447),
                                    strand=-1),
                                type="CDS",
                                qualifiers=self.alice_cds_139_qualifier_dict)

        self.alice_cds_139_coords = (88120, 88447)
        self.alice_cds_139_data_in_db = {
            "GeneID": "Alice_CDS_2",
            "PhageID": "Alice",
            "Start": 88120,
            "Stop": 88447,
            "Length": 108,
            "Name": "139",
            "translation": self.alice_cds_139_translation,
            "Orientation": "R",
            "Notes": "",
            "LocusTag": "ALICE_139"
            }


        self.alice_tmrna_169_qualifier_dict = OrderedDict(
            [('gene', ['169']),
             ('locus_tag', ['ALICE_169'])])
        self.alice_tmrna_169 = SeqFeature(
                                    FeatureLocation(
                                        ExactPosition(95923),
                                        ExactPosition(96358),
                                        strand=1),
                                    type="tmRNA",
                                    qualifiers=self.alice_tmrna_169_qualifier_dict)

        self.alice_trna_170_qualifier_dict = OrderedDict(
            [('gene', ['170']),
             ('locus_tag', ['ALICE_170']),
             ('product', ['tRNA-Gln']),
             ('note', ['tRNA-Gln (ttg)'])])
        self.alice_trna_170 = SeqFeature(
                                FeatureLocation(
                                    ExactPosition(96431),
                                    ExactPosition(96507),
                                    strand=1),
                                type="tRNA",
                                qualifiers=self.alice_trna_170_qualifier_dict)

        self.alice_cds_193_translation = (
            "MGNNRPTTRTLPTGEKATHHPDGRLVLKPKNSLADALSGQLDAQQEASQALTEALV"
            "QTAVTKREAQADGNPVPEDKRVF"
            )
        self.alice_cds_193_qualifier_dict = OrderedDict(
            [("gene", ["193"]),
             ("locus_tag", ["ALICE_193"]),
             ("note", ["gp193"]),
             ("codon_start", ["1"]),
             ("transl_table", ["11"]),
             ("product", ["hypothetical protein"]),
             ("protein_id", ["AEJ94419.1"]),
             ("translation", [self.alice_cds_193_translation])])
        self.alice_cds_193 = SeqFeature(
                                FeatureLocation(
                                    ExactPosition(110297),
                                    ExactPosition(110537),
                                    strand=1),
                                type="CDS",
                                qualifiers=self.alice_cds_193_qualifier_dict)

        self.alice_cds_193_coords = (110297, 110537)
        self.alice_cds_193_data_in_db = {
            "GeneID": "Alice_CDS_3",
            "PhageID": "Alice",
            "Start": 110297,
            "Stop": 110537,
            "Length": 79,
            "Name": "193",
            "translation": self.alice_cds_193_translation,
            "Orientation": "F",
            "Notes": "",
            "LocusTag": "ALICE_193"
            }



        self.alice_source_qualifiers = OrderedDict(
            [("organism", ["Mycobacterium phage Alice"]),
             ("mol_type", ["genomic DNA"]),
             ("isolation_source", ["soil"]),
             ("db_xref", ["taxon:1034128"]),
             ("lab_host", ["Mycobacterium smegmatis mc2 155"]),
             ("country", ["USA: Aledo, TX"]),
             ("lat_lon", ["31.69 N 97.65 W"]),
             ("collection_date", ["01-Oct-2009"]),
             ("collected_by", ["C. Manley"]),
             ("identified_by", ["C. Manley"])])
        self.alice_source = SeqFeature(
                                FeatureLocation(
                                    ExactPosition(0),
                                    ExactPosition(153401),
                                    strand=1),
                                type="source",
                                qualifiers=self.alice_source_qualifiers)
        self.alice_feature_list = [self.alice_source,
                                   self.alice_cds_252, # Wrap around gene
                                   self.alice_cds_124, # Compound gene
                                   self.alice_cds_139, # Bottom strand normal CDS
                                   self.alice_tmrna_169,
                                   self.alice_trna_170,
                                   self.alice_cds_193 # Top strand normal CDS
                                   ]

        self.alice_ref1 = Reference()
        self.alice_ref1.authors = "Hatfull,G.F."
        self.alice_ref2 = Reference()
        self.alice_ref2.authors = (
            "Alferez,G.I., Bryan,W.J., Byington,E.L., Contreras,T.D., "
            "Evans,C.R., Griffin,J.A., Jalal,M.D., Lindsey,C.B., "
            "Manley,C.M., Mitchell,A.M., O'Hara,J.M., Onoh,U.M., "
            "Padilla,E., Penrod,L.C., Regalado,M.S., Reis,K.E., "
            "Ruprecht,A.M., Slater,A.E., Staton,A.C., Tovar,I.G., "
            "Turek,A.J., Utech,N.T., Simon,S.E., Hughes,L.E., "
            "Benjamin,R.C., Serrano,M.G., Lee,V., Hendricks,S.L., "
            "Sheth,N.U., Buck,G.A., Bradley,K.W., Khaja,R., "
            "Lewis,M.F., Barker,L.P., Jordan,T.C., Russell,D.A., "
            "Leuba,K.D., Fritz,M.J., Bowman,C.A., Pope,W.H., "
            "Jacobs-Sera,D., Hendrix,R.W. and Hatfull,G.F.")
        self.alice_ref_list = [self.alice_ref1, self.alice_ref2]
        self.alice_description = "Mycobacterium phage Alice, complete sequence"
        self.alice_organism = "Mycobacterium phage Alice"
        self.alice_source = "Mycobacterium phage Alice"
        self.alice_date = "23-MAR-2018"
        self.alice_accessions = ["JF704092"]
        self.alice_annotation_dict = {
            "accessions": self.alice_accessions,
            "source": self.alice_source,
            "organism": self.alice_organism,
            "references": self.alice_ref_list,
            "date": self.alice_date}
        self.alice_seq = get_seq(base_flat_file_path)
        self.alice_record = SeqRecord(
                                seq = self.alice_seq,
                                id = "JF704092.1",
                                name = "JF704092",
                                annotations = self.alice_annotation_dict,
                                description = self.alice_description,
                                features = self.alice_feature_list
                                )

        self.alice_flat_file = "temp_alice.gb"
        self.alice_flat_file_path = Path(self.genome_folder, "temp_alice.gb")


        self.alice_phage_table_data = {
            "PhageID": "Alice",
            "Accession": "JF704092",
            "Name": "Alice",
            "HostStrain": "Mycobacterium",
            "Sequence": get_seq(base_flat_file_path),
            "SequenceLength": 153401,
            "GC": 64.6808,
            "status": "final",
            "DateLastModified": datetime.strptime('3/23/2018', '%m/%d/%Y'),
            "RetrieveRecord": 1,
            "AnnotationAuthor": 1,
            "Cluster": "C1",
            "Cluster2": "C",
            "Subcluster2": "C1"
            }


    def tearDown(self):
        shutil.rmtree(self.base_dir)
        remove_db(db, user, pwd)



    @patch("pdm_utils.pipelines.db_import.import_genome.setup_sql_handle")
    def test_import_pipeline_1(self, setup_sql_mock):
        """Test pipeline with:
        valid add ticket,
        valid flat file."""

        # Set Alice data for 'draft' genome.
        # for feature in self.alice_record.features:
        #     if "product" in feature.qualifiers.keys():
        #         print("product: ", feature.qualifiers["product"])
        #     if "function" in feature.qualifiers.keys():
        #         print("function: ", feature.qualifiers["function"])
        #     if "note" in feature.qualifiers.keys():
        #         print("note: ", feature.qualifiers["note"])
        clear_descriptions(self.alice_record)

        self.alice_record.annotations["accessions"] = [""]
        self.alice_record.id = ""
        self.alice_record.annotations["organism"] = \
            "Mycobacterium phage Alice_Draft"
        self.alice_record.annotations["source"] = \
            "Mycobacterium phage Alice_Draft"
        self.alice_record.description = \
            "Mycobacterium phage Alice_Draft, complete sequence"

        self.alice_record.features[0].qualifiers["organism"] = [
            "Mycobacterium phage Alice_Draft"]


        self.alice_phage_table_data["Name"] = "Alice_Draft"
        self.alice_phage_table_data["Accession"] = ""
        self.alice_phage_table_data["status"] = "draft"
        self.alice_phage_table_data["DateLastModified"] = current_date

        self.alice_cds_252_data_in_db["LocusTag"] = ""
        self.alice_cds_124_data_in_db["LocusTag"] = ""
        self.alice_cds_139_data_in_db["LocusTag"] = ""
        self.alice_cds_193_data_in_db["LocusTag"] = ""

        self.alice_cds_252_data_in_db["Notes"] = ""
        self.alice_cds_124_data_in_db["Notes"] = ""
        self.alice_cds_139_data_in_db["Notes"] = ""
        self.alice_cds_193_data_in_db["Notes"] = ""
        # self.alice_cds_data_in_db = {
        #     "GeneID": "",
        #     "PhageID": "",
        #     "Start": "",
        #     "Stop": "",
        #     "Length": "",
        #     "Name": "",
        #     "translation": "",
        #     "Orientation": "",
        #     "Notes": "",
        #     "LocusTag": ""
        #     }







        # Export Alice data to flat file.
        SeqIO.write(self.alice_record, self.alice_flat_file_path, "genbank")

        # Set up the database.
        setup_sql_mock.return_value = self.sql_handle
        insert_data_into_phage_table(db, user, pwd, d29_phage_table_data)
        insert_data_into_phage_table(db, user, pwd, redrock_phage_table_data)
        insert_data_into_phage_table(db, user, pwd, trixie_phage_table_data)
        # insert_data_into_gene_table(db, user, pwd, trixie_gene_table_data_1)

        # results = get_sql_data(db, user, pwd, phage_table_query)
        # print(len(results))
        # results2 = get_sql_data(db, user, pwd, gene_table_query)
        # print(len(results2))
        # input("paused")
        # seq = get_seq(self.test_flat_file1)
        # l5_phage_table_data["sequence"] = seq
        # l5_phage_table_data["SequenceLength"] = len(seq)
        # print(l5_phage_table_data["sequence"][:25])
        # print(l5_phage_table_data["SequenceLength"])
        # input("pause")
        # insert_data_into_phage_table(db, user, pwd, l5_phage_table_data)
        # input("pause2")
        # keys = set(["host_genus", "cluster", "subcluster"])
        # l5_ticket_data_retrieve = set_data(l5_ticket_data_complete, keys, "retrieve")


        list_of_data_dicts = [alice_ticket_data_complete]
        create_import_table(list_of_data_dicts, self.import_table)
        # shutil.copy(self.base_flat_file_path, self.genome_folder)
        unparsed_args = ["run.py", "import_dev", db,
                         str(self.genome_folder),
                         str(self.import_table),
                         "-g", "_organism_name",
                         "-p",
                         "-r", "phagesdb",
                         "-d", "product",
                         "-o", str(self.output_folder),
                         "-l", str(self.log_file)
                         ]
        run.main(unparsed_args)

        phage_table_results = get_sql_data(db, user, pwd, phage_table_query)
        process_phage_table_data(phage_table_results)

        output_genome_data = filter_genome_data(phage_table_results, "Alice")
        genome_errors = compare_data(self.alice_phage_table_data,
                                     output_genome_data)

        gene_table_results = get_sql_data(db, user, pwd, gene_table_query)
        process_gene_table_data(gene_table_results)


        cds252_data = filter_gene_data(gene_table_results, self.alice_cds_252_coords)
        cds124_data = filter_gene_data(gene_table_results, self.alice_cds_124_coords)
        cds139_data = filter_gene_data(gene_table_results, self.alice_cds_139_coords)
        cds193_data = filter_gene_data(gene_table_results, self.alice_cds_193_coords)

        cds252_errors = compare_data(self.alice_cds_252_data_in_db, cds252_data)
        cds124_errors = compare_data(self.alice_cds_124_data_in_db, cds124_data)
        cds139_errors = compare_data(self.alice_cds_139_data_in_db, cds139_data)
        cds193_errors = compare_data(self.alice_cds_193_data_in_db, cds193_data)

        for cdsftr in gene_table_results:
            print(cdsftr["GeneID"])
            print(cdsftr["Start"])
            print(cdsftr["Stop"])



        with self.subTest():
            self.assertEqual(len(phage_table_results), 4)
        with self.subTest():
            self.assertEqual(len(gene_table_results), 4)
        with self.subTest():
            self.assertEqual(genome_errors, 0)
        # TODO test for no errors for each CDS feature.

        input("check log")

        # TODO test the exact data in the database by querying and comparing?


        #1. insert alice w/1 CDS
        #2. try to insert alice if already in db
        #3. try to insert alice if another genome has same seq.
        #4. use invalid ticket
        #5. use minimal ticket.
        #6. fail due to genome-level error
        #7. fail due to genome-pair error.
        #8. replace alice genome.
        #9. fail due to CDS error.
        #10. fail due to source error.





if __name__ == '__main__':
    unittest.main()
