"""Helper functions for building mock databases for integration tests."""

from datetime import datetime
from pathlib import Path
import subprocess

import pymysql

from pdm_utils.constants import constants

# Format of the date the script imports into the database.
current_date = datetime.today().replace(hour=0, minute=0,
                                        second=0, microsecond=0)


# It is expected that the 'pdm_anon' MySQL user has all privileges
# for 'pdm_test_*' databases.
USER = "pdm_anon"
PWD = "pdm_anon"
DB = "pdm_test_db"
unittest_file = Path(__file__)
test_dir = unittest_file.parent
test_file_dir = Path(test_dir, "test_files")
SCHEMA_VERSION = constants.CODE_SCHEMA_VERSION
SCHEMA_FILE = f"test_schema_{SCHEMA_VERSION}.sql" # Empty schema
SCHEMA_FILEPATH = Path(test_file_dir, SCHEMA_FILE)
TEST_DB_FILEPATH = Path(test_file_dir, "test_db.sql") # Contains data
VERSION_TABLE_DATA = {"Version": 1, "SchemaVersion": SCHEMA_VERSION}

# Common queries

phage_table_query = "SELECT * FROM phage;"
gene_table_query = "SELECT * FROM gene;"
gene_domain_table_query = "SELECT * FROM gene_domain;"
domain_table_query = "SELECT * FROM domain;"
version_table_query = "SELECT * FROM version;"

# SQLAlchemy setup

def create_engine_string(db=DB, user=USER, pwd=PWD):
    """Generate engine string for SQLAlchemy."""
    engine_string = f"mysql+pymysql://{user}:{pwd}@localhost/{db}"
    return engine_string




# Statement execution

def execute(statement, db=DB, user=USER, pwd=PWD):
    """Execute a statement."""
    # Sometimes a connection is needed without a pre-specified database.
    if db is None:
        database = None
    else:
        database = db
    connection = pymysql.connect(host="localhost",
                                 user=user,
                                 password=pwd,
                                 database=database,
                                 cursorclass=pymysql.cursors.DictCursor)
    cur = connection.cursor()
    cur.execute(statement)
    results = cur.fetchall()
    cur.close()
    connection.commit()
    connection.close()
    return results

def get_data(query, db=DB, user=USER, pwd=PWD):
    """Get data from the database."""
    connection = pymysql.connect(host="localhost",
                                 user=user,
                                 password=pwd,
                                 database=db,
                                 cursorclass=pymysql.cursors.DictCursor)
    cur = connection.cursor()
    cur.execute(query)
    results = cur.fetchall()
    cur.close()
    connection.close()
    return results




# Database construction

def create_empty_test_db(db=DB, user=USER, pwd=PWD):
    """Creates a test database with the current schema version and no data."""
    create_new_db(schema_filepath=SCHEMA_FILEPATH, db=db, user=user, pwd=pwd,
                  version_table_data=VERSION_TABLE_DATA)

def create_filled_test_db(db=DB, user=USER, pwd=PWD):
    """Creates a test database with the current schema version and with data."""
    # No need to add data to version table, since test_db.sql already
    # has data in this table.
    create_new_db(schema_filepath=TEST_DB_FILEPATH, db=db, user=user, pwd=pwd)

def check_if_exists(db=DB, user=USER, pwd=PWD):
    """Checks whether database exists or not."""
    query = ("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA "
             f"WHERE SCHEMA_NAME = '{db}'")
    # No need to connect to a specific database.
    result = get_data(query, db=None, user=user, pwd=pwd)
    if len(result) != 0:
        exists = True
    else:
        exists = False
    return exists

def create_new_db(schema_filepath=None, db=DB, user=USER, pwd=PWD,
                  version_table_data=None):
    """Deletes a database if it exists, then creates a new database."""
    # First, test if a test database already exists within mysql.
    # If there is, delete it so that a fresh test database is installed.
    exists = check_if_exists(db=db, user=user, pwd=pwd)
    if exists:
        remove_db(db=db, user=user, pwd=pwd)

    # Next, create the database within mysql.
    create_db(db=db, user=user, pwd=pwd)

    # Now import the empty schema from file and add version data.
    if schema_filepath is not None:
        install_db(schema_filepath, db=db, user=user, pwd=pwd)
    if version_table_data is not None:
        insert_version_data(version_table_data, db=db, user=user, pwd=pwd)

def create_db(db=DB, user=USER, pwd=PWD):
    """Create a MySQL database for the test."""
    statement = f"CREATE DATABASE {db}"
    execute(statement, db=None, user=user, pwd=pwd)

def remove_db(db=DB, user=USER, pwd=PWD):
    """Remove the MySQL database created for the test."""
    statement = f"DROP DATABASE {db}"
    execute(statement, db=db, user=user, pwd=pwd)

def install_db(schema_filepath, db=DB, user=USER, pwd=PWD):
    """Install data into database from a SQL file."""
    # Seems like pymysql has trouble with this step, so use subprocess.
    handle = open(schema_filepath, "r")
    command_string = f"mysql -u {user} -p{pwd} {db}"
    command_list = command_string.split(" ")
    proc = subprocess.check_call(command_list, stdin=handle)
    handle.close()

def create_schema_file(schema_filepath, db=DB, user=USER, pwd=PWD):
    """Dump empty schema file for a database."""
    # "mysqldump --no-data -u root -p --skip-comments <db_name> > db_schema_before.sql"
    command_string = f"mysqldump --no-data -u {user} -p{pwd} --skip-comments {db}"
    command_list = command_string.split(" ")
    with schema_filepath.open("w") as handle:
        subprocess.check_call(command_list, stdout=handle)



# Inserting data into specific tables

def insert_version_data(data_dict=VERSION_TABLE_DATA, db=DB, user=USER, pwd=PWD):
    """Insert data into the version table."""
    # Unlike other tables, version table data is expected to have one row
    # that reflects the structure of the database. So by default, data_dict
    # has a default value.
    statement = (
        "INSERT INTO version "
        "(Version, SchemaVersion) "
        "VALUES ("
        f"'{data_dict['Version']}', '{data_dict['SchemaVersion']}');"
        )
    execute(statement, db=db, user=user, pwd=pwd)

def insert_domain_data(data_dict, db=DB, user=USER, pwd=PWD):
    """Insert data into the domain table."""
    # No need to insert ID, since that is auto-incremented.
    statement = (
        "INSERT INTO domain  "
        "(HitID, DomainID, Name, Description) "
        "VALUES ("
        f"'{data_dict['HitID']}', '{data_dict['DomainID']}',"
        f"'{data_dict['Name']}', '{data_dict['Description']}');"
        )
    execute(statement, db=db, user=user, pwd=pwd)

def insert_gene_domain_data(data_dict, db=DB, user=USER, pwd=PWD):
    """Insert data into the gene_domain table."""
    # No need to insert ID, since that is auto-incremented.
    statement = (
        "INSERT INTO gene_domain "
        "(GeneID, HitID, Expect, QueryStart, QueryEnd) "
        "VALUES ("
        f"'{data_dict['GeneID']}', '{data_dict['HitID']}', "
        f"{data_dict['Expect']}, {data_dict['QueryStart']}, "
        f"{data_dict['QueryEnd']});"
        )
    execute(statement, db=db, user=user, pwd=pwd)

def insert_phage_data(data_dict, db=DB, user=USER, pwd=PWD):
    """Insert data into the phage table."""

    cluster = data_dict['Cluster']
    if cluster != "NULL":
        cluster = "'{}'".format(cluster)

    subcluster = data_dict['Subcluster']
    if subcluster != "NULL":
        subcluster = "'{}'".format(subcluster)

    notes = data_dict['Notes']
    if notes != "NULL":
        notes = "'{}'".format(notes)

    # Don't encapsulate Cluster or Subcluster with quotes, since they
    # can sometimes be assigned NULL.
    statement = (
        "INSERT INTO phage "
        "(PhageID, Accession, Name, "
        "HostGenus, Sequence, Length, GC, Status, "
        "DateLastModified, RetrieveRecord, AnnotationAuthor, "
        "Cluster, Subcluster, Notes) "
        "VALUES ("
        f"'{data_dict['PhageID']}', '{data_dict['Accession']}', "
        f"'{data_dict['Name']}', '{data_dict['HostGenus']}', "
        f"'{data_dict['Sequence']}', {data_dict['Length']}, "
        f"{data_dict['GC']}, '{data_dict['Status']}', "
        f"'{data_dict['DateLastModified']}', "
        f"{data_dict['RetrieveRecord']}, "
        f"{data_dict['AnnotationAuthor']}, "
        f"{cluster}, {subcluster}, {notes});"
        )
    execute(statement, db=db, user=user, pwd=pwd)

def insert_gene_data(data_dict, db=DB, user=USER, pwd=PWD):
    """Insert data into the gene table."""
    # Since PhamID and DomainStatus can be auto-generated,
    # they may not be in the data_dict.
    if "PhamID" not in data_dict.keys():
        data_dict["PhamID"] = "NULL"
    if "DomainStatus" not in data_dict.keys():
        data_dict["DomainStatus"] = 0

    statement = (
        "INSERT INTO gene "
        "(GeneID, PhageID, Start, Stop, Length, Name, "
        "Translation, Orientation, Notes, LocusTag, "
        "Parts, DomainStatus, PhamID) "
        "VALUES ("
        f"'{data_dict['GeneID']}', '{data_dict['PhageID']}', "
        f"{data_dict['Start']}, {data_dict['Stop']}, "
        f"{data_dict['Length']}, '{data_dict['Name']}', "
        f"'{data_dict['Translation']}', '{data_dict['Orientation']}', "
        f"'{data_dict['Notes']}', '{data_dict['LocusTag']}', "
        f"{data_dict['Parts']}, {data_dict['DomainStatus']}, "
        f"{data_dict['PhamID']});"
        )
    execute(statement, db=db, user=user, pwd=pwd)




# Data processing

def process_phage_table_data(list_of_sql_results):
    """Converts datatype for data retrieved from a few fields
    in the phage table."""
    for x in range(len(list_of_sql_results)):
        data_dict = list_of_sql_results[x]
        data_dict["Sequence"] = data_dict["Sequence"].decode("utf-8")
        data_dict["Length"] = int(data_dict["Length"])
        data_dict["GC"] = float(data_dict["GC"])

        if data_dict["Notes"] is None:
            data_dict["Notes"] = "NULL"
        else:
            data_dict["Notes"] = data_dict["Notes"].decode("utf-8")
        if data_dict["Cluster"] is None:
            data_dict["Cluster"]  = "NULL"
        if data_dict["Subcluster"] is None:
            data_dict["Subcluster"]  = "NULL"

def process_gene_table_data(list_of_sql_results):
    """Converts datatype for data retrieved from a few fields
    in the gene table."""
    for x in range(len(list_of_sql_results)):
        data_dict = list_of_sql_results[x]
        data_dict["Notes"] = data_dict["Notes"].decode("utf-8")
        data_dict["Length"] = int(data_dict["Length"])
        data_dict["Start"] = int(data_dict["Start"])
        data_dict["Stop"] = int(data_dict["Stop"])
        data_dict["Parts"] = int(data_dict["Parts"])

        if data_dict["PhamID"] is None:
            data_dict["PhamID"]  = "NULL"


def filter_genome_data(list_of_sql_results, phage_id):
    """Returns a dictionary of data from the phage table
    based on the indicated PhageID."""
    output_dict = {}
    for x in range(len(list_of_sql_results)):
        if list_of_sql_results[x]["PhageID"] == phage_id:
            output_dict = list_of_sql_results[x]
        else:
            pass
    return output_dict

def filter_gene_data(list_of_sql_results, coordinates):
    """Returns a dictionary of data from the gene table
    based on a tuple of coordinates."""
    output_dict = {}
    for x in range(len(list_of_sql_results)):
        if (list_of_sql_results[x]["Start"] == coordinates[0] and
                list_of_sql_results[x]["Stop"] == coordinates[1]):
            output_dict = list_of_sql_results[x]
        else:
            pass
    return output_dict
