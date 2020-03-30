from pathlib import Path
from sqlalchemy import create_engine
import pymysql
import subprocess

#Global constants
USER = "pdm_anon"
PASSWORD = "pdm_anon"
DATABASE = "test_db"
CONNECTION_STRING = f"mysql+pymysql://{USER}:{PASSWORD}@localhost/{DATABASE}"
# TEST_DB_PATH = Path("../tests/integration_heavy/test_files/test_db.sql")
TEST_DB_PATH = Path("../tests/test_files/test_db.sql")
TEMP_DIR = Path.cwd()

def setup_test_db():
    connection = pymysql.connect(host="localhost",
                                 user=USER,
                                 password=PASSWORD,
                                 cursorclass=pymysql.cursors.DictCursor)
    cur = connection.cursor()

    # First, test if a test database already exists within mysql.
    # If there is, delete it so that a fresh test database is installed.
    sql_query = "SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA " + \
               f"WHERE SCHEMA_NAME = '{DATABASE}'"
    cur.execute(sql_query)
    result = cur.fetchall()

    if len(result) != 0:
        cur.execute(f"DROP DATABASE {DATABASE}")
        connection.commit()

    # Next, create the database within mysql.
    cur.execute(f"CREATE DATABASE {DATABASE}")
    connection.commit()
    connection.close()

    # Now import the empty schema from file.
    # Seems like pymysql has trouble with this step, so use subprocess.
    handle = open(TEST_DB_PATH, "r")
    command_string = (f"mysql -u {USER} -p{PASSWORD} {DATABASE} "
                      f"--execute")
    command_list = command_string.split(" ")
    command_list.append(f"SOURCE {TEST_DB_PATH}")
    proc = subprocess.check_call(command_list, stdin=handle)
    handle.close()

def connect_test_db():
    engine = create_engine(CONNECTION_STRING, echo=False)

    return engine

def teardown_test_db():
    connection = pymysql.connect(host="localhost",
                                     user=USER,
                                     password=PASSWORD,
                                     cursorclass=pymysql.cursors.DictCursor)
    cur = connection.cursor()

    cur.execute(f"DROP DATABASE {DATABASE}")
    connection.commit()

    connection.close()
