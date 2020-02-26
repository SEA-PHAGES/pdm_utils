"""Integration tests for the main find_domains pipeline."""


from datetime import datetime, date
import csv
import getpass
from pathlib import Path
import pymysql
import shutil
import subprocess
import time
import unittest
from unittest.mock import patch
import sqlalchemy
from pdm_utils.constants import constants
from pdm_utils.functions import basic, mysqldb
from pdm_utils.pipelines import find_domains



# The following integration tests user the 'pdm_anon' MySQL user.
# It is expected that this user has all privileges for 'test_db' database.
user = "pdm_anon"
pwd = "pdm_anon"
db = "test_db"
db2 = "Actinobacteriophage"
engine_string1 = f"mysql+pymysql://{user}:{pwd}@localhost/{db}"
engine_string2 = f"mysql+pymysql://{user}:{pwd}@localhost/{db2}"


unittest_file = Path(__file__)
unittest_dir = unittest_file.parent
schema_version = constants.CODE_SCHEMA_VERSION
schema_file = f"test_schema_{schema_version}.sql"
schema_filepath = Path(unittest_dir, "test_files/", schema_file)
version_table_data = {"Version":1, "SchemaVersion":schema_version}


def create_new_db(schema_filepath, db, user, pwd):
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


def insert_data_into_version_table(db, user, pwd, data_dict):
    """Insert data into the version table."""
    connection = pymysql.connect(host = "localhost",
                                 user = user,
                                 password = pwd,
                                 database = db,
                                 cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()
    sql = (
        "INSERT INTO version "
        "(Version, SchemaVersion) "
        "VALUES ("
        f"'{data_dict['Version']}', '{data_dict['SchemaVersion']}');"
        )
    cur.execute(sql)
    connection.commit()
    connection.close()


def insert_data_into_domain_table(db, user, pwd, data_dict):
    """Insert data into the domain table."""
    connection = pymysql.connect(host = "localhost",
                                 user = user,
                                 password = pwd,
                                 database = db,
                                 cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()
    sql = (
        "INSERT INTO domain  "
        "(HitID, DomainID, Name, Description) "
        "VALUES ("
        f"'{data_dict['HitID']}', '{data_dict['DomainID']}',"
        f"'{data_dict['Name']}', '{data_dict['Description']}');"
        )
    cur.execute(sql)
    connection.commit()
    connection.close()




def insert_data_into_gene_domain_table(db, user, pwd, data_dict):
    """Insert data into the gene_domain table."""
    connection = pymysql.connect(host = "localhost",
                                 user = user,
                                 password = pwd,
                                 database = db,
                                 cursorclass = pymysql.cursors.DictCursor)
    cur = connection.cursor()
    sql = (
        "INSERT INTO gene_domain  "
        "(GeneID, HitID, Expect, QueryStart, QueryEnd) "
        "VALUES ("
        f"'{data_dict['GeneID']}', '{data_dict['HitID']}',"
        f"{data_dict['Expect']}, {data_dict['QueryStart']}"
        f"{data_dict['QueryEnd']});"
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
        "Translation, Orientation, Notes, LocusTag) "
        "VALUES ("
        f"'{data_dict['GeneID']}', '{data_dict['PhageID']}', "
        f"{data_dict['Start']}, {data_dict['Stop']}, "
        f"{data_dict['Length']}, '{data_dict['Name']}', "
        f"'{data_dict['Translation']}', '{data_dict['Orientation']}', "
        f"'{data_dict['Notes']}', '{data_dict['LocusTag']}');"
        )
    cur.execute(sql)
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
        "HostGenus, Sequence, Length, GC, Status, "
        "DateLastModified, RetrieveRecord, AnnotationAuthor, "
        "Cluster, Subcluster) "
        "VALUES ("
        f"'{data_dict['PhageID']}', '{data_dict['Accession']}', "
        f"'{data_dict['Name']}', '{data_dict['HostGenus']}', "
        f"'{data_dict['Sequence']}', {data_dict['Length']}, "
        f"{data_dict['GC']}, '{data_dict['Status']}', "
        f"'{data_dict['DateLastModified']}', "
        f"{data_dict['RetrieveRecord']}, "
        f"{data_dict['AnnotationAuthor']}, "
        f"'{data_dict['Cluster']}', '{data_dict['Subcluster']}');"
        )
    cur.execute(sql)
    connection.commit()
    connection.close()





def get_trixie_phage_table_data():
    """Mock phage table data for Trixie."""
    dict = {
        "PhageID": "Trixie",
        "Accession": "BCD456",
        "Name": "Trixie",
        "HostGenus": "Gordonia",
        "Sequence": "GGGGGGGGGGGGGGGGGGGG",
        "Length": 20,
        "GC": 1,
        "Status": "final",
        "DateLastModified": datetime.strptime('1/1/2000', '%m/%d/%Y'),
        "RetrieveRecord": 1,
        "AnnotationAuthor": 1,
        "Cluster": "A",
        "Subcluster": "A3"
        }
    return dict


def get_trixie_gene_table_data_1():
    """Mock gene table data for Trixie."""
    dict = {
        "GeneID": "TRIXIE_0001",
        "PhageID": "Trixie",
        "Start": 100,
        "Stop": 1100,
        "Length": 1000,
        "Name": "1",
        "Translation": "ACTGC",
        "Orientation": "F",
        "Notes": "int",
        "LocusTag": "SEA_TRIXIE_0001"
        }
    return dict

def get_trixie_gene_domain_table_data_1():
    """Mock gene_domain table data for Trixie."""
    dict = {
        "GeneID": "TRIXIE_0001",
        "HitID": "gnl|CDD|334841",
        "Expect": 1.78531e-11,
        "QueryStart": 33,
        "QueryEnd": 115
        }
    return dict

def get_trixie_domain_table_data_1():
    """Mock domain table data for Trixie."""
    dict = {
        "HitID": "gnl|CDD|334841",
        "DomainID": "pfam02195",
        "Name": "ParBc",
        "Description": "ParB-like nuclease domain"
        }
    return dict


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


phage_table_query = (
    "SELECT "
    "PhageID, Accession, Name, "
    "HostGenus, Sequence, Length, GC, Status, "
    "DateLastModified, RetrieveRecord, AnnotationAuthor, "
    "Cluster, Subcluster "
    "FROM phage;")

gene_table_query = (
    "SELECT "
    "GeneID, PhageID, Start, Stop, Length, Name, "
    "Translation, Orientation, Notes, LocusTag, DomainStatus "
    "FROM gene;")

gene_domain_table_query = (
    "SELECT "
    "GeneID, HitID, Expect, QueryStart, QueryEnd "
    "FROM gene_domain;")

domain_table_query = (
    "SELECT "
    "HitID, DomainID, Name, Description "
    "FROM domain;")


def get_domain_insert_statement(data_dict):
    statement = find_domains.INSERT_INTO_DOMAIN.format(
                    data_dict["HitID"], data_dict["DomainID"],
                    data_dict["Name"], data_dict["Description"]
                    )
    return statement

def get_gene_domain_insert_statement(data_dict):
    statement = find_domains.INSERT_INTO_GENE_DOMAIN.format(
                    data_dict["GeneID"], data_dict["HitID"],
                    float(data_dict["Expect"]), int(data_dict["QueryStart"]),
                    int(data_dict["QueryEnd"])
                    )
    return statement

def get_gene_update_statement(data_dict):
    statement = find_domains.UPDATE_GENE.format(
                    data_dict["GeneID"]
                    )
    return statement


# INSERT IGNORE INTO domain (HitID, DomainID, Name, Description) VALUES ("gnl|CDD|334841", "pfam02195", "ParBc", "ParB-like nuclease domain.")
# INSERT IGNORE INTO gene_domain (GeneID, HitID, Expect, QueryStart, QueryEnd) VALUES ("SEA_ABBA_66", "gnl|CDD|334841", 1.78531e-11, 33, 115)

# INSERT_INTO_DOMAIN = """INSERT INTO domain (HitID, DomainID, Name, Description) VALUES ("{}", "{}", "{}", "{}")"""
# INSERT_INTO_GENE_DOMAIN = """INSERT INTO gene_domain (GeneID, HitID, Expect, QueryStart, QueryEnd) VALUES ("{}", "{}", {}, {}, {})"""
# UPDATE_GENE = "UPDATE gene SET DomainStatus = 1 WHERE GeneID = '{}'"


class TestFindDomains1(unittest.TestCase):

    def setUp(self):
        create_new_db(schema_filepath, db, user, pwd)
        insert_data_into_version_table(db, user, pwd, version_table_data)
        insert_data_into_phage_table(db, user, pwd, get_trixie_phage_table_data())
        insert_data_into_gene_table(db, user, pwd, get_trixie_gene_table_data_1())

        self.engine = sqlalchemy.create_engine(engine_string1, echo=False)
        self.connection = self.engine.connect()
        self.trans = self.connection.begin()


    def tearDown(self):
        remove_db(db, user, pwd)
        self.trans.rollback()
        self.engine.dispose()



    def test_execute_statement_1(self):
        """Verify valid data can be inserted into domain table and
        gene_domain table."""
        # INSERT_INTO_DOMAIN
        # INSERT_INTO_GENE_DOMAIN
        # UPDATE_GENE
        domain_data = get_trixie_domain_table_data_1()
        statement1 = get_domain_insert_statement(domain_data)
        find_domains.execute_statement(self.connection, statement1)
        gene_domain_data = get_trixie_gene_domain_table_data_1()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        find_domains.execute_statement(self.connection, statement2)
        self.trans.commit()
        phage_table_results = get_sql_data(db, user, pwd, phage_table_query)
        gene_table_results = get_sql_data(db, user, pwd, gene_table_query)
        gene_domain_table_results = get_sql_data(db, user, pwd, gene_domain_table_query)
        domain_table_results = get_sql_data(db, user, pwd, domain_table_query)
        with self.subTest():
            self.assertEqual(len(phage_table_results), 1)
        with self.subTest():
            self.assertEqual(len(gene_table_results), 1)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 1)
        with self.subTest():
            self.assertEqual(len(gene_domain_table_results), 1)
















###
