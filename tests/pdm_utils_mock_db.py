"""Helper functions for building mock databases for integration tests."""

from datetime import datetime
from pathlib import Path
import pymysql
import subprocess

# Format of the date the script imports into the database.
current_date = datetime.today().replace(hour=0, minute=0,
                                        second=0, microsecond=0)


# It is expected that the 'pdm_anon' MySQL user has all privileges for 'test_db' database.
user = "pdm_anon"
pwd = "pdm_anon"
db = "test_db"
unittest_file = Path(__file__)
test_dir = unittest_file.parent
test_file_dir = Path(test_dir, "test_files")
schema_version = 8
schema_file = f"test_schema_{schema_version}.sql"
schema_filepath = Path(test_file_dir, schema_file)
version_table_data = {"Version":1, "SchemaVersion":schema_version}

phage_table_query = "SELECT * FROM phage;"
gene_table_query = "SELECT * FROM gene;"
gene_domain_table_query = "SELECT * FROM gene_domain;"
domain_table_query = "SELECT * FROM domain;"
version_table_query = "SELECT * FROM version;"


def execute(statement, db=db, user=user, pwd=pwd):
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

def get_data(query, db=db, user=user, pwd=pwd):
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


def create_new_db(schema_filepath=schema_filepath, db=db, user=user, pwd=pwd):
    """Creates a new, empty database."""
    connection = pymysql.connect(host="localhost",
                                 user=user,
                                 password=pwd,
                                 cursorclass=pymysql.cursors.DictCursor)
    cur = connection.cursor()

    # First, test if a test database already exists within mysql.
    # If there is, delete it so that a fresh test database is installed.
    statement = ("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA "
                 f"WHERE SCHEMA_NAME = '{db}'")

    result = execute(statement, db=None, user=user, pwd=pwd)
    if len(result) != 0:
        remove_db(db=db, user=user, pwd=pwd)

    # Next, create the database within mysql.
    create_db(db=db, user=user, pwd=pwd)

    # Now import the empty schema from file.
    # Seems like pymysql has trouble with this step, so use subprocess.
    install_db(schema_filepath, db=db, user=user, pwd=pwd)


def create_db(db=db, user=user, pwd=pwd):
    """Create a MySQL database for the test."""
    statement = f"CREATE DATABASE {db}"
    execute(statement, db=None, user=user, pwd=pwd)

def remove_db(db=db, user=user, pwd=pwd):
    """Remove the MySQL database created for the test."""
    statement = f"DROP DATABASE {db}"
    execute(statement, db=db, user=user, pwd=pwd)

def install_db(schema_filepath, db=db, user=user, pwd=pwd):
    """Install data into database from a SQL file."""
    handle = open(schema_filepath, "r")
    command_string = f"mysql -u {user} -p{pwd} {db}"
    command_list = command_string.split(" ")
    proc = subprocess.check_call(command_list, stdin=handle)
    handle.close()



def insert_version_data(data_dict, db=db, user=user, pwd=pwd):
    """Insert data into the version table."""
    statement = (
        "INSERT INTO version "
        "(Version, SchemaVersion) "
        "VALUES ("
        f"'{data_dict['Version']}', '{data_dict['SchemaVersion']}');"
        )
    execute(statement, db=db, user=user, pwd=pwd)





def insert_domain_data(data_dict, db=db, user=user, pwd=pwd):
    """Insert data into the domain table."""
    statement = (
        "INSERT INTO domain  "
        "(HitID, DomainID, Name, Description) "
        "VALUES ("
        f"'{data_dict['HitID']}', '{data_dict['DomainID']}',"
        f"'{data_dict['Name']}', '{data_dict['Description']}');"
        )
    execute(statement, db=db, user=user, pwd=pwd)





def insert_gene_domain_data(data_dict, db=db, user=user, pwd=pwd):
    """Insert data into the gene_domain table."""
    statement = (
        "INSERT INTO gene_domain  "
        "(GeneID, HitID, Expect, QueryStart, QueryEnd) "
        "VALUES ("
        f"'{data_dict['GeneID']}', '{data_dict['HitID']}',"
        f"{data_dict['Expect']}, {data_dict['QueryStart']}"
        f"{data_dict['QueryEnd']});"
        )
    execute(statement, db=db, user=user, pwd=pwd)














def insert_phage_data(data_dict, db=db, user=user, pwd=pwd):
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







def insert_gene_data(data_dict, db=db, user=user, pwd=pwd):
    """Insert data into the gene table."""
    statement = (
        "INSERT INTO gene "
        "(GeneID, PhageID, Start, Stop, Length, Name, "
        "Translation, Orientation, Notes, LocusTag, Parts) "
        "VALUES ("
        f"'{data_dict['GeneID']}', '{data_dict['PhageID']}', "
        f"{data_dict['Start']}, {data_dict['Stop']}, "
        f"{data_dict['Length']}, '{data_dict['Name']}', "
        f"'{data_dict['Translation']}', '{data_dict['Orientation']}', "
        f"'{data_dict['Notes']}', '{data_dict['LocusTag']}',"
        f"{data_dict['Parts']});"
        )
    execute(statement, db=db, user=user, pwd=pwd)


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
