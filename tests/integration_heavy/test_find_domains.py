"""Integration tests for the main find_domains pipeline."""

import logging
from pathlib import Path
import shutil
import sys
import unittest
from unittest.mock import patch
import sqlalchemy
from pdm_utils import run
from pdm_utils.pipelines import find_domains

# Import helper functions to build mock database and mock flat files
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import test_db_utils
import test_data_utils


# Create the main test directory in which all files will be
# created and managed.
test_root_dir = Path("/tmp", "pdm_utils_tests_find_domains")
if test_root_dir.exists() == True:
    shutil.rmtree(test_root_dir)
test_root_dir.mkdir()

# New folder that will get created/removed for each test.
test_folder = Path(test_root_dir, "output")

# How the output folder and files are named in the find_domains pipeline.
results_folder = Path(find_domains.RESULTS_FOLDER)
results_path = Path(test_folder, results_folder)

# Set up a log file to catch all logging for review.
# Note: this should overwrite logging output file in find_domains pipeline.
# Log file is kept outside of 'output' folder.
log = Path(test_root_dir, "test_log.txt")
log = log.expanduser()
log = log.resolve()
logging.basicConfig(filename=log, filemode="w",level=logging.DEBUG)

# The following integration tests user the 'pdm_anon' MySQL user.
# It is expected that this user has all privileges for 'pdm_test_db' database.
engine_string = test_db_utils.create_engine_string()

# Assumes that output message contains "SQLAlchemy Error..."
error_msg = "SQLAlchemy"


pipeline = "find_domains"
DB = test_db_utils.DB
DB2 = "Actinobacteriophage"

def get_unparsed_args():
    """Returns list of command line arguments to find domains."""
    unparsed_args = ["run.py", pipeline, DB,
                     "-t", str(1),
                     "-o", str(test_folder)]
    return unparsed_args


def get_trixie_gene_table_domain_status_update_data_1():
    """Mock gene table DomainStatus update data for Trixie."""
    dict = {
        "GeneID": "TRIXIE_0001"
        }
    return dict

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

def get_all_domain_data(db):
    """Retrieve all data from domain table."""
    query = "SELECT * FROM domain"
    results = test_db_utils.get_data(query, db=db)
    results2 = []
    for i in range(len(results)):
        data = results[i]
        data["Description"] = data["Description"].decode("utf-8")
        results2.append(data)
    return results2

def statements_to_txns(statement_list):
    """Convert list of statements to list of lists, where each sub-list
    represents a collection of statements meant to be executed as a single
    transaction."""
    new_list = [[stmt] for stmt in statement_list]
    return new_list

def get_all_domain_txns(db):
    """Get all domain data from a database and convert it to a list of
    find_domain INSERT statement transactions."""
    results = get_all_domain_data(db=db)
    statements = [get_domain_insert_statement(result) for result in results]
    txns = statements_to_txns(statements)
    return txns


def get_gene_id_dict(list_of_results):
    """Get a dictionary of data where key = GeneID."""
    dict1 = {}
    for i, dict2 in enumerate(list_of_results):
        key = dict2["GeneID"]
        if key in dict1.keys():
            # list1 = dict1[key]
            # list1.append(list_of_results[i])
            # dict1[key] = list1
            # list1.append(list_of_results[i])
            dict1[key].append(list_of_results[i])
        else:
            dict1[key] = [list_of_results[i]]
    return dict1


# INSERT IGNORE INTO domain (HitID, DomainID, Name, Description) VALUES ("gnl|CDD|334841", "pfam02195", "ParBc", "ParB-like nuclease domain.")
# INSERT IGNORE INTO gene_domain (GeneID, HitID, Expect, QueryStart, QueryEnd) VALUES ("SEA_ABBA_66", "gnl|CDD|334841", 1.78531e-11, 33, 115)

# INSERT_INTO_DOMAIN = """INSERT INTO domain (HitID, DomainID, Name, Description) VALUES ("{}", "{}", "{}", "{}")"""
# INSERT_INTO_GENE_DOMAIN = """INSERT INTO gene_domain (GeneID, HitID, Expect, QueryStart, QueryEnd) VALUES ("{}", "{}", {}, {}, {})"""
# UPDATE_GENE = "UPDATE gene SET DomainStatus = 1 WHERE GeneID = '{}'"



class TestFindDomains10(unittest.TestCase):

    def setUp(self):
        test_db_utils.create_empty_test_db()
        test_db_utils.insert_phage_data(test_data_utils.get_trixie_phage_data())
        self.gene_data_1 = test_data_utils.get_trixie_gene_data()
        self.gene_data_1["GeneID"] = "TRIXIE_0001"
        self.gene_data_2 = test_data_utils.get_trixie_gene_data()
        self.gene_data_2["GeneID"] = "TRIXIE_0002"
        self.gene_data_3 = test_data_utils.get_trixie_gene_data()
        self.gene_data_3["GeneID"] = "TRIXIE_0003"

        test_db_utils.insert_gene_data(self.gene_data_1)

        self.domain_data = test_data_utils.get_trixie_domain_data()
        self.gene_domain_data = test_data_utils.get_trixie_gene_domain_data()

        self.rps_hit_1 = {
            "HitID": "hit_1",
            "DomainID": self.domain_data["DomainID"],
            "Name": self.domain_data["Name"],
            "Description": self.domain_data["Description"],
            "Expect": self.gene_domain_data["Expect"],
            "QueryStart": self.gene_domain_data["QueryStart"],
            "QueryEnd": self.gene_domain_data["QueryEnd"]
            }

        self.rps_hit_2 = {
            "HitID": "hit_2",
            "DomainID": self.domain_data["DomainID"],
            "Name": self.domain_data["Name"],
            "Description": self.domain_data["Description"],
            "Expect": self.gene_domain_data["Expect"],
            "QueryStart": self.gene_domain_data["QueryStart"],
            "QueryEnd": self.gene_domain_data["QueryEnd"]
            }

        self.rps_data = [self.rps_hit_1, self.rps_hit_2]

        translation = "ABCDE"
        self.cds_translation_dict = {
            translation: {self.gene_data_1["GeneID"],
                          self.gene_data_2["GeneID"],
                          self.gene_data_3["GeneID"]}
             }
        self.rps_translation_dict = {translation: self.rps_data}


        self.rps_result_1 = {"Translation": "ABCDE", "Data": self.rps_data}
        self.rps_result_2 = {"Translation": "FGHIJ", "Data": self.rps_data}
        self.rps_result_3 = {"Translation": "MNOPQ", "Data": self.rps_data}
        self.rps_result_4 = {"Translation": "RSTUV", "Data": self.rps_data}

        self.rps_results = [self.rps_result_1, self.rps_result_2,
                            self.rps_result_3, self.rps_result_4]



    def tearDown(self):
        test_db_utils.remove_db()


    def test_construct_gene_update_stmt_1(self):
        """Verify gene table data can be updated."""
        # Use previously validated data.
        stmt = find_domains.construct_gene_update_stmt(self.gene_data_1["GeneID"])
        test_db_utils.execute(stmt)
        results = test_db_utils.get_data(test_db_utils.gene_table_query)
        status = results[0]["DomainStatus"]
        self.assertEqual(status, 1)

    def test_construct_domain_stmt_1(self):
        """Verify domain table data can be inserted."""
        # Use previously validated data.
        hit_id_1 = self.domain_data["HitID"]
        domain_id_1 = self.domain_data["DomainID"]
        name_1 = self.domain_data["Name"]
        desc_1 = self.domain_data["Description"]

        stmt = find_domains.construct_domain_stmt(self.domain_data)
        test_db_utils.execute(stmt)
        results = test_db_utils.get_data(test_db_utils.domain_table_query)

        hit_id_2 = results[0]["HitID"]
        domain_id_2 = results[0]["DomainID"]
        name_2 = results[0]["Name"]
        desc_2 = results[0]["Description"].decode("utf-8")

        with self.subTest():
            self.assertEqual(len(results), 1)
        with self.subTest():
            self.assertEqual(hit_id_1, hit_id_2)
        with self.subTest():
            self.assertEqual(domain_id_1, domain_id_2)
        with self.subTest():
            self.assertEqual(name_1, name_2)
        with self.subTest():
            self.assertEqual(desc_1, desc_2)

    def test_construct_gene_domain_stmt_1(self):
        """Verify gene_domain table data can be inserted."""
        # Use previously validated gene_domain data.
        # Valid domain data needs to be inserted first.
        test_db_utils.insert_domain_data(self.domain_data)

        # Pop gene_id. The dictionary for find_domains doesn't need it as a key.
        gene_id_1 = self.gene_domain_data.pop("GeneID")
        hit_id_1 = self.gene_domain_data["HitID"]
        expect_1 = self.gene_domain_data["Expect"]
        query_start_1 = self.gene_domain_data["QueryStart"]
        query_end_1 = self.gene_domain_data["QueryEnd"]
        stmt = find_domains.construct_gene_domain_stmt(self.gene_domain_data,
                                                       gene_id_1)
        test_db_utils.execute(stmt)
        results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
        gene_id_2 = results[0]["GeneID"]
        hit_id_2 = results[0]["HitID"]
        expect_2 = results[0]["Expect"]
        query_start_2 = results[0]["QueryStart"]
        query_end_2 = results[0]["QueryEnd"]

        with self.subTest():
            self.assertEqual(len(results), 1)
        with self.subTest():
            self.assertEqual(gene_id_1, gene_id_2)
        with self.subTest():
            self.assertEqual(hit_id_1, hit_id_2)
        with self.subTest():
            self.assertEqual(expect_1, expect_2)
        with self.subTest():
            self.assertEqual(query_start_1, query_start_2)
        with self.subTest():
            self.assertEqual(query_end_1, query_end_2)






    def test_construct_sql_txn_1(self):
        """Verify list of SQL statements representing one transaction
        is constructed."""
        gene_id = self.gene_data_1["GeneID"]
        txn = find_domains.construct_sql_txn(gene_id, self.rps_data)
        self.assertEqual(len(txn), 5)


    def test_construct_sql_txns_1(self):
        """Verify list of list of SQL statements representing one transaction
        each are constructed."""
        txns = find_domains.construct_sql_txns(self.cds_translation_dict,
                                              self.rps_translation_dict)
        with self.subTest():
            self.assertEqual(len(txns), 3)
        with self.subTest():
            self.assertEqual(len(txns[0]), 5)
        with self.subTest():
            self.assertEqual(len(txns[1]), 5)
        with self.subTest():
            self.assertEqual(len(txns[2]), 5)




    def test_create_results_dict_1(self):
        """Verify dictionary is constucted correctly."""
        dict = find_domains.create_results_dict(self.rps_results)
        with self.subTest():
            self.assertEqual(len(dict.keys()), 4)
        with self.subTest():
            self.assertEqual(len(dict["ABCDE"]), 2)
        with self.subTest():
            self.assertEqual(len(dict["FGHIJ"]), 2)




class TestFindDomains1(unittest.TestCase):

    def setUp(self):
        test_db_utils.create_empty_test_db()
        test_db_utils.insert_phage_data(test_data_utils.get_trixie_phage_data())
        test_db_utils.insert_gene_data(test_data_utils.get_trixie_gene_data())

        self.engine = sqlalchemy.create_engine(engine_string, echo=False)
        self.connection = self.engine.connect()
        self.trans = self.connection.begin()


    def tearDown(self):
        self.trans.rollback()
        self.engine.dispose()
        test_db_utils.remove_db()


    def test_execute_statement_1(self):
        """Verify valid DomainStatus update data can be inserted
        into gene table."""
        gene_table_results1 = test_db_utils.get_data(test_db_utils.gene_table_query)
        update_data = get_trixie_gene_table_domain_status_update_data_1()
        statement = get_gene_update_statement(update_data)
        results_tup = find_domains.execute_statement(self.connection, statement)
        result = results_tup[0]
        type_error = results_tup[1]
        value_error = results_tup[2]
        msg = results_tup[3]
        self.trans.commit()
        phage_table_results = test_db_utils.get_data(test_db_utils.phage_table_query)
        gene_table_results2 = test_db_utils.get_data(test_db_utils.gene_table_query)
        gene_domain_table_results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        domain_status1 = gene_table_results1[0]["DomainStatus"]
        domain_status2 = gene_table_results2[0]["DomainStatus"]
        with self.subTest():
            self.assertEqual(len(phage_table_results), 1)
        with self.subTest():
            self.assertEqual(len(gene_table_results2), 1)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 0)
        with self.subTest():
            self.assertEqual(len(gene_domain_table_results), 0)
        with self.subTest():
            self.assertEqual(domain_status1, 0)
        with self.subTest():
            self.assertEqual(domain_status2, 1)
        with self.subTest():
            self.assertEqual(result, 0)
        with self.subTest():
            self.assertFalse(type_error)
        with self.subTest():
            self.assertFalse(value_error)
        with self.subTest():
            self.assertIsInstance(msg, str)


    def test_execute_statement_2(self):
        """Verify valid data can be inserted into domain table."""
        domain_data = test_data_utils.get_trixie_domain_data()
        statement = get_domain_insert_statement(domain_data)
        results_tup = find_domains.execute_statement(self.connection, statement)
        result = results_tup[0]
        type_error = results_tup[1]
        value_error = results_tup[2]
        msg = results_tup[3]
        self.trans.commit()
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 1)
        with self.subTest():
            self.assertEqual(result, 0)
        with self.subTest():
            self.assertFalse(type_error)
        with self.subTest():
            self.assertFalse(value_error)


    def test_execute_statement_3(self):
        """Verify valid data can be inserted into domain and gene_domain tables."""
        domain_data = test_data_utils.get_trixie_domain_data()
        statement1 = get_domain_insert_statement(domain_data)
        results_tup1 = find_domains.execute_statement(self.connection, statement1)
        result1 = results_tup1[0]
        type_error1 = results_tup1[1]
        value_error1 = results_tup1[2]
        msg1 = results_tup1[3]
        gene_domain_data = test_data_utils.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        results_tup2 = find_domains.execute_statement(self.connection, statement2)
        result2 = results_tup2[0]
        type_error2 = results_tup2[1]
        value_error2 = results_tup2[2]
        msg2 = results_tup2[3]
        self.trans.commit()
        gene_domain_table_results = test_db_utils.get_data(
                                        test_db_utils.gene_domain_table_query)
        domain_table_results = test_db_utils.get_data(
                                        test_db_utils.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 1)
        with self.subTest():
            self.assertEqual(len(gene_domain_table_results), 1)
        with self.subTest():
            self.assertEqual(result1, 0)
        with self.subTest():
            self.assertEqual(result2, 0)
        with self.subTest():
            self.assertFalse(type_error1)
        with self.subTest():
            self.assertFalse(type_error2)
        with self.subTest():
            self.assertFalse(value_error1)
        with self.subTest():
            self.assertFalse(value_error2)


    def test_execute_statement_4(self):
        """Verify invalid data can NOT be inserted into gene_domain table
        when there is no matching 'HitID' in the domain table."""
        gene_domain_data = test_data_utils.get_trixie_gene_domain_data()
        statement = get_gene_domain_insert_statement(gene_domain_data)
        results_tup = find_domains.execute_statement(self.connection, statement)
        result = results_tup[0]
        type_error = results_tup[1]
        value_error = results_tup[2]
        msg = results_tup[3]
        self.trans.commit()
        gene_domain_table_results = test_db_utils.get_data(
                                        test_db_utils.gene_domain_table_query)
        with self.subTest():
            self.assertEqual(len(gene_domain_table_results), 0)
        with self.subTest():
            self.assertEqual(result, 1)
        with self.subTest():
            self.assertFalse(type_error)
        with self.subTest():
            self.assertFalse(value_error)


    def test_execute_statement_5(self):
        """Verify invalid data can NOT be inserted into domain table when
        there is a duplicated HitID."""
        domain_data = test_data_utils.get_trixie_domain_data()
        statement = get_domain_insert_statement(domain_data)
        results_tup1 = find_domains.execute_statement(self.connection, statement)
        result1 = results_tup1[0]
        type_error1 = results_tup1[1]
        value_error1 = results_tup1[2]
        msg1 = results_tup1[3]
        self.trans.commit()
        domain_table_results1 = test_db_utils.get_data(
                                        test_db_utils.domain_table_query)
        new_trans = self.connection.begin()
        results_tup2 = find_domains.execute_statement(
                                        self.connection, statement)
        result2 = results_tup2[0]
        type_error2 = results_tup2[1]
        value_error2 = results_tup2[2]
        msg2 = results_tup2[3]
        new_trans.commit()
        domain_table_results2 = test_db_utils.get_data(
                                        test_db_utils.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results1), 1)
        with self.subTest():
            self.assertEqual(len(domain_table_results2), 1)
        with self.subTest():
            self.assertEqual(result1, 0)
        with self.subTest():
            self.assertEqual(result2, 0)
        with self.subTest():
            self.assertFalse(type_error1)
        with self.subTest():
            self.assertFalse(type_error2)
        with self.subTest():
            self.assertFalse(value_error1)
        with self.subTest():
            self.assertFalse(value_error2)
        with self.subTest():
            self.assertFalse(error_msg in msg1)
        with self.subTest():
            self.assertTrue(error_msg in msg2)


    def test_execute_statement_6(self):
        """Verify invalid data can NOT be inserted into domain table when
        there is an invalid column name."""
        domain_data = test_data_utils.get_trixie_domain_data()
        statement = get_domain_insert_statement(domain_data)
        statement = statement.replace("Name", "Name_invalid")
        results_tup = find_domains.execute_statement(self.connection, statement)
        result = results_tup[0]
        type_error = results_tup[1]
        value_error = results_tup[2]
        msg = results_tup[3]
        self.trans.commit()
        domain_table_results = test_db_utils.get_data(
                                        test_db_utils.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 0)
        with self.subTest():
            self.assertEqual(result, 1)
        with self.subTest():
            self.assertTrue(error_msg in msg)
        with self.subTest():
            self.assertFalse(type_error)
        with self.subTest():
            self.assertFalse(value_error)


    def test_execute_statement_7(self):
        """Verify invalid data can NOT be inserted due to '%'."""
        domain_data = test_data_utils.get_trixie_domain_data()
        # "Description": "ParB-like nuclease domain"
        description = domain_data["Description"]
        description = description.replace("nuclease domain", "nuclease % domain")
        domain_data["Description"] = description
        statement = get_domain_insert_statement(domain_data)
        results_tup = find_domains.execute_statement(self.connection, statement)
        result = results_tup[0]
        type_error = results_tup[1]
        value_error = results_tup[2]
        msg = results_tup[3]
        self.trans.commit()
        domain_table_results = test_db_utils.get_data(
                                        test_db_utils.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 0)
        with self.subTest():
            self.assertEqual(result, 1)
        with self.subTest():
            self.assertFalse(error_msg in msg)
        with self.subTest():
            self.assertTrue(type_error)
        with self.subTest():
            self.assertFalse(value_error)


    def test_execute_statement_8(self):
        """Verify invalid data can be inserted after '%' is
        replaced with '%%'."""
        domain_data = test_data_utils.get_trixie_domain_data()
        # "Description": "ParB-like nuclease domain"
        description = domain_data["Description"]
        description = description.replace("nuclease domain", "nuclease % domain")
        domain_data["Description"] = description
        statement = get_domain_insert_statement(domain_data)
        statement = statement.replace("%", "%%")
        results_tup = find_domains.execute_statement(self.connection, statement)
        result = results_tup[0]
        type_error = results_tup[1]
        value_error = results_tup[2]
        msg = results_tup[3]
        self.trans.commit()
        domain_table_results = test_db_utils.get_data(
                                            test_db_utils.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 1)
        with self.subTest():
            self.assertEqual(result, 0)
        with self.subTest():
            self.assertFalse(type_error)
        with self.subTest():
            self.assertFalse(value_error)


    def test_execute_statement_9(self):
        """Verify invalid data can NOT be inserted due to '% w'."""
        domain_data = test_data_utils.get_trixie_domain_data()
        # "Description": "ParB-like nuclease domain"
        description = domain_data["Description"]
        description = description.replace("nuclease domain", "nuclease % wdomain")
        domain_data["Description"] = description
        statement = get_domain_insert_statement(domain_data)
        results_tup = find_domains.execute_statement(self.connection, statement)
        result = results_tup[0]
        type_error = results_tup[1]
        value_error = results_tup[2]
        msg = results_tup[3]
        self.trans.commit()
        domain_table_results = test_db_utils.get_data(
                                            test_db_utils.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 0)
        with self.subTest():
            self.assertEqual(result, 1)
        with self.subTest():
            self.assertFalse(error_msg in msg)
        with self.subTest():
            self.assertFalse(type_error)
        with self.subTest():
            self.assertTrue(value_error)


    def test_execute_statement_10(self):
        """Verify invalid data can be inserted after '% w' is
        replaced with '%% w'."""
        domain_data = test_data_utils.get_trixie_domain_data()
        # "Description": "ParB-like nuclease domain"
        description = domain_data["Description"]
        description = description.replace("nuclease domain", "nuclease % wdomain")
        domain_data["Description"] = description
        statement = get_domain_insert_statement(domain_data)
        statement = statement.replace("%", "%%")
        results_tup = find_domains.execute_statement(self.connection, statement)
        result = results_tup[0]
        type_error = results_tup[1]
        value_error = results_tup[2]
        msg = results_tup[3]
        self.trans.commit()
        domain_table_results = test_db_utils.get_data(
                                            test_db_utils.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 1)
        with self.subTest():
            self.assertEqual(result, 0)
        with self.subTest():
            self.assertFalse(type_error)
        with self.subTest():
            self.assertFalse(value_error)




class TestFindDomains2(unittest.TestCase):

    def setUp(self):
        test_db_utils.create_empty_test_db()
        test_db_utils.insert_phage_data(test_data_utils.get_trixie_phage_data())
        test_db_utils.insert_gene_data(test_data_utils.get_trixie_gene_data())

        self.engine = sqlalchemy.create_engine(engine_string, echo=False)
        self.connection = self.engine.connect()


    def tearDown(self):
        self.engine.dispose()
        test_db_utils.remove_db()




    def test_execute_transaction_1(self):
        """Verify function runs with list of zero statements."""
        result = find_domains.execute_transaction(self.connection)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 0)
        with self.subTest():
            self.assertEqual(result, 0)


    def test_execute_transaction_2(self):
        """Verify list of one valid statement can be inserted into domain table."""
        domain_data = test_data_utils.get_trixie_domain_data()
        statement = get_domain_insert_statement(domain_data)
        statements = [statement]
        result = find_domains.execute_transaction(self.connection, statements)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 1)
        with self.subTest():
            self.assertEqual(result, 0)


    def test_execute_transaction_3(self):
        """Verify list of two valid statements can be inserted into
        domain and gene_domain tables."""
        domain_data = test_data_utils.get_trixie_domain_data()
        statement1 = get_domain_insert_statement(domain_data)
        gene_domain_data = test_data_utils.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        statements = [statement1, statement2]
        result = find_domains.execute_transaction(self.connection, statements)
        gene_domain_table_results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 1)
        with self.subTest():
            self.assertEqual(len(gene_domain_table_results), 1)
        with self.subTest():
            self.assertEqual(result, 0)


    def test_execute_transaction_4(self):
        """Verify list of three statements (including one with
        duplicated HitID) are inserted."""
        domain_data1 = test_data_utils.get_trixie_domain_data()
        test_db_utils.insert_domain_data(domain_data1)
        domain_table_results1 = test_db_utils.get_data(test_db_utils.domain_table_query)
        # Duplicate HitID
        statement1 = get_domain_insert_statement(domain_data1)
        # Valid
        gene_domain_data = test_data_utils.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        # Valid
        update_data = get_trixie_gene_table_domain_status_update_data_1()
        statement3 = get_gene_update_statement(update_data)

        statements = [statement1, statement2, statement3]
        result = find_domains.execute_transaction(self.connection, statements)
        gene_table_results = test_db_utils.get_data(test_db_utils.gene_table_query)
        gene_domain_table_results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
        domain_table_results2 = test_db_utils.get_data(test_db_utils.domain_table_query)
        domain_status = gene_table_results[0]["DomainStatus"]
        with self.subTest():
            self.assertEqual(len(domain_table_results1), 1)
        with self.subTest():
            self.assertEqual(len(domain_table_results2), 1)
        with self.subTest():
            self.assertEqual(len(gene_domain_table_results), 1)
        with self.subTest():
            self.assertEqual(result, 0)
        with self.subTest():
            self.assertEqual(domain_status, 1)


    def test_execute_transaction_5(self):
        """Verify list of three valid statements and one invalid statement
        are NOT inserted. All statements rolled back."""
        # Valid
        domain_data1 = test_data_utils.get_trixie_domain_data()
        statement1 = get_domain_insert_statement(domain_data1)
        # Valid
        gene_domain_data = test_data_utils.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        # Invalid
        domain_data2 = test_data_utils.get_trixie_domain_data()
        statement3 = get_domain_insert_statement(domain_data2)
        statement3 = statement3.replace("HitID", "unique_id")
        statement3 = statement3.replace("Name", "Name_invalid")
        # Valid - function should exit before executing this though.
        update_data = get_trixie_gene_table_domain_status_update_data_1()
        statement4 = get_gene_update_statement(update_data)

        statements = [statement1, statement2, statement3, statement4]
        result = find_domains.execute_transaction(self.connection, statements)
        gene_table_results = test_db_utils.get_data(test_db_utils.gene_table_query)
        gene_domain_table_results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        domain_status = gene_table_results[0]["DomainStatus"]
        with self.subTest():
            self.assertEqual(len(domain_table_results), 0)
        with self.subTest():
            self.assertEqual(len(gene_domain_table_results), 0)
        with self.subTest():
            self.assertEqual(result, 1)
        with self.subTest():
            self.assertEqual(domain_status, 0)


    def test_execute_transaction_6(self):
        """Verify list of three valid statements and one invalid statement
        (containing '%') are inserted (since '%' replaced with '%%')."""
        # Valid
        domain_data1 = test_data_utils.get_trixie_domain_data()
        statement1 = get_domain_insert_statement(domain_data1)
        # Valid
        gene_domain_data = test_data_utils.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        # Invalid '%'
        domain_data2 = test_data_utils.get_trixie_domain_data()
        # "Description": "ParB-like nuclease domain"
        description = domain_data2["Description"]
        description = description.replace("nuclease domain", "nuclease % domain")
        domain_data2["Description"] = description
        domain_data2["HitID"] = "unique_id"
        statement3 = get_domain_insert_statement(domain_data2)
        # Valid
        update_data = get_trixie_gene_table_domain_status_update_data_1()
        statement4 = get_gene_update_statement(update_data)

        statements = [statement1, statement2, statement3, statement4]
        result = find_domains.execute_transaction(self.connection, statements)
        gene_table_results = test_db_utils.get_data(test_db_utils.gene_table_query)
        gene_domain_table_results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        domain_status = gene_table_results[0]["DomainStatus"]
        with self.subTest():
            self.assertEqual(len(domain_table_results), 2)
        with self.subTest():
            self.assertEqual(len(gene_domain_table_results), 1)
        with self.subTest():
            self.assertEqual(result, 0)
        with self.subTest():
            self.assertEqual(domain_status, 1)


    def test_execute_transaction_7(self):
        """Verify list of three valid statements and one invalid statement
        (containing '% w') are inserted (since '% w' replaced with '%% w')."""
        # Valid
        domain_data1 = test_data_utils.get_trixie_domain_data()
        statement1 = get_domain_insert_statement(domain_data1)
        # Valid
        gene_domain_data = test_data_utils.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        # Invalid '% w'
        domain_data2 = test_data_utils.get_trixie_domain_data()
        # "Description": "ParB-like nuclease domain"
        description = domain_data2["Description"]
        description = description.replace("nuclease domain", "nuclease % wdomain")
        domain_data2["Description"] = description
        domain_data2["HitID"] = "unique_id"
        statement3 = get_domain_insert_statement(domain_data2)
        # Valid
        update_data = get_trixie_gene_table_domain_status_update_data_1()
        statement4 = get_gene_update_statement(update_data)

        statements = [statement1, statement2, statement3, statement4]
        result = find_domains.execute_transaction(self.connection, statements)
        gene_table_results = test_db_utils.get_data(test_db_utils.gene_table_query)
        gene_domain_table_results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        domain_status = gene_table_results[0]["DomainStatus"]
        with self.subTest():
            self.assertEqual(len(domain_table_results), 2)
        with self.subTest():
            self.assertEqual(len(gene_domain_table_results), 1)
        with self.subTest():
            self.assertEqual(result, 0)
        with self.subTest():
            self.assertEqual(domain_status, 1)


    # TODO this isn't the best way to test the rolllback in the
    # try/except block. The execute_statement is patched, so no changes are
    # actually made to the database. This test just confirms that the except
    # block is entered.
    @patch("pdm_utils.pipelines.find_domains.execute_statement")
    def test_execute_transaction_8(self, es_mock):
        """Verify error inside try/except block is executed."""
        stmt_result1 = 0
        type_error1 = False
        # TODO make sure this is set correctly
        value_error1 = False

        msg1 = "empty"
        mock_result1 = (stmt_result1, type_error1, value_error1, msg1)

        stmt_result2 = 0
        type_error2 = False
        # TODO make sure this is set correctly
        value_error2 = False

        msg2 = 2 # the function expects this to be a string, so this should
                 # break the code and trigger the except block.
        mock_result2 = (stmt_result2, type_error2, value_error2, msg2)
        es_mock.side_effect = [mock_result1, mock_result2]
        # Valid
        domain_data1 = test_data_utils.get_trixie_domain_data()
        statement1 = get_domain_insert_statement(domain_data1)
        # Valid
        gene_domain_data = test_data_utils.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)

        statements = [statement1, statement2]
        result = find_domains.execute_transaction(self.connection, statements)
        gene_domain_table_results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 0)
        with self.subTest():
            self.assertEqual(len(gene_domain_table_results), 0)
        with self.subTest():
            self.assertEqual(result, 1)
        with self.subTest():
            self.assertEqual(es_mock.call_count, 2)




class TestFindDomains3(unittest.TestCase):

    def setUp(self):
        test_db_utils.create_empty_test_db()
        test_db_utils.insert_phage_data(test_data_utils.get_trixie_phage_data())

        cds1 = test_data_utils.get_trixie_gene_data() # GeneID = "TRIXIE_0001"
        cds2 = test_data_utils.get_trixie_gene_data()
        cds3 = test_data_utils.get_trixie_gene_data()
        cds2["GeneID"] = "TRIXIE_0002"
        cds3["GeneID"] = "TRIXIE_0003"
        test_db_utils.insert_gene_data(cds1)
        test_db_utils.insert_gene_data(cds2)
        test_db_utils.insert_gene_data(cds3)
        stmt = "UPDATE gene SET DomainStatus = 0"
        test_db_utils.execute(stmt)
        self.engine = sqlalchemy.create_engine(engine_string, echo=False)

    def tearDown(self):
        test_db_utils.remove_db()
        self.engine.dispose()




    def test_insert_domain_data_1(self):
        """Verify list of domain data can be inserted."""
        # The purpose of this test is to take thousands of rows of domain data
        # from an existing database such as Actinobacteriophage and try to
        # reinsert that data into a test database using the find_domains
        # pipeline. This ensures that as the find_domains pipeline is modified,
        # it is still able to insert all domain data previously encountered.
        # This test could be improved if it is possible to retrieve all
        # domain data from the Conserved Domain Database en masse.
        logging.info("test_insert_domain_data_1")
        txns = get_all_domain_txns(DB2)
        find_domains.insert_domain_data(self.engine, txns)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        rows = len(domain_table_results)
        self.assertEqual(len(txns), rows)


    def test_insert_domain_data_2(self):
        """Verify one transaction of three valid statements can be
        inserted into domain, gene_domain, and gene tables."""
        logging.info("test_insert_domain_data_2")

        domain_data = test_data_utils.get_trixie_domain_data()
        statement1 = get_domain_insert_statement(domain_data)
        # GeneID = "TRIXIE_0001"
        gene_domain_data = test_data_utils.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        update_data = get_trixie_gene_table_domain_status_update_data_1()
        statement3 = get_gene_update_statement(update_data)
        statements = [statement1, statement2, statement3]
        txns = [statements]
        find_domains.insert_domain_data(self.engine, txns)

        gene_domain_table_results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        gene_table_results = test_db_utils.get_data(test_db_utils.gene_table_query)
        gene_table_dict = {}
        for dict in gene_table_results:
            gene_table_dict[dict["GeneID"]] = dict

        domain_status = gene_table_dict["TRIXIE_0001"]["DomainStatus"]
        d_rows = len(domain_table_results)
        gd_rows = len(gene_domain_table_results)

        with self.subTest():
            self.assertEqual(d_rows, 1)
        with self.subTest():
            self.assertEqual(gd_rows, 1)
        with self.subTest():
            self.assertEqual(domain_status, 1)


    def test_insert_domain_data_3(self):
        """Verify one transaction of three valid statements and one
        invalid statement are NOT inserted. All statements rolled back."""
        logging.info("test_insert_domain_data_3")

        # Valid
        domain_data1 = test_data_utils.get_trixie_domain_data()
        statement1 = get_domain_insert_statement(domain_data1)
        # Valid
        gene_domain_data = test_data_utils.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        # Invalid
        domain_data2 = test_data_utils.get_trixie_domain_data()
        statement3 = get_domain_insert_statement(domain_data2)
        statement3 = statement3.replace("HitID", "unique_id")
        statement3 = statement3.replace("Name", "Name_invalid")
        # Valid - function should exit before executing this though.
        update_data = get_trixie_gene_table_domain_status_update_data_1()
        statement4 = get_gene_update_statement(update_data)
        statements = [statement1, statement2, statement3, statement4]
        txns = [statements]
        find_domains.insert_domain_data(self.engine, txns)

        gene_domain_table_results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        gene_table_results = test_db_utils.get_data(test_db_utils.gene_table_query)
        gene_table_dict = {}
        for dict in gene_table_results:
            gene_table_dict[dict["GeneID"]] = dict

        domain_status = gene_table_dict["TRIXIE_0001"]["DomainStatus"]
        d_rows = len(domain_table_results)
        gd_rows = len(gene_domain_table_results)

        with self.subTest():
            self.assertEqual(d_rows, 0)
        with self.subTest():
            self.assertEqual(gd_rows, 0)
        with self.subTest():
            self.assertEqual(domain_status, 0)


    def test_insert_domain_data_4(self):
        """Verify two transactions of valid statements are inserted and
        one transaction of invalid statement is NOT inserted.
        Only invalid transaction statements are rolled back."""
        logging.info("test_insert_domain_data_4")

        # Transaction 1 - all valid
        # "GeneID": "TRIXIE_0001"
        t1_domain_data = test_data_utils.get_trixie_domain_data()
        t1_statement1 = get_domain_insert_statement(t1_domain_data)
        t1_gene_domain_data = test_data_utils.get_trixie_gene_domain_data()
        t1_statement2 = get_gene_domain_insert_statement(t1_gene_domain_data)
        t1_update_data = get_trixie_gene_table_domain_status_update_data_1()
        t1_statement3 = get_gene_update_statement(t1_update_data)
        t1 = [t1_statement1, t1_statement2, t1_statement3]

        # Transaction 2 - invalid
        t2_domain_data2 = test_data_utils.get_trixie_domain_data()
        t2_statement1 = get_domain_insert_statement(t2_domain_data2)
        t2_statement1 = t2_statement1.replace("HitID", "unique_id")
        t2_statement1 = t2_statement1.replace("Name", "Name_invalid")
        t2 = [t2_statement1]

        # Transaction 3 - all valid
        # "GeneID": "TRIXIE_0002"
        t3_gene_domain_data = test_data_utils.get_trixie_gene_domain_data()
        t3_gene_domain_data["GeneID"] = "TRIXIE_0002"
        t3_statement1 = get_gene_domain_insert_statement(t3_gene_domain_data)
        t3_update_data = get_trixie_gene_table_domain_status_update_data_1()
        t3_update_data["GeneID"] = "TRIXIE_0002"
        t3_statement2 = get_gene_update_statement(t3_update_data)
        t3 = [t3_statement1, t3_statement2]

        txns = [t1, t2, t3]
        find_domains.insert_domain_data(self.engine, txns)

        gene_domain_table_results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        gene_table_results = test_db_utils.get_data(test_db_utils.gene_table_query)
        gene_table_dict = {}
        for dict in gene_table_results:
            gene_table_dict[dict["GeneID"]] = dict

        domain_status1 = gene_table_dict["TRIXIE_0001"]["DomainStatus"]
        domain_status2 = gene_table_dict["TRIXIE_0002"]["DomainStatus"]
        d_rows = len(domain_table_results)
        gd_rows = len(gene_domain_table_results)

        with self.subTest():
            self.assertEqual(d_rows, 1)
        with self.subTest():
            self.assertEqual(gd_rows, 2)
        with self.subTest():
            self.assertEqual(domain_status1, 1)
        with self.subTest():
            self.assertEqual(domain_status2, 1)



class TestFindDomains4(unittest.TestCase):

    def setUp(self):
        test_folder.mkdir()
        test_db_utils.create_empty_test_db()
        test_db_utils.insert_phage_data(test_data_utils.get_trixie_phage_data())

        # Translation from Trixie gene 35, contains multiple conserved domains.
        translation = (
            "MQASYISPVDGQRYWGPRNYDNRMDAEAWLASEKRLIDMEQWTPPAERAKKEAANS"
            "ITVEEYTKKWLAERDLAEGTRELYKTHARKRIYPVLGDVAVAEMTPALVRAWWAGM"
            "GKDYPTARRHAYNVLRAVMNTAVEDKLLTENPCRIEQKAAAERDVEALTPEELDIV"
            "AGEVLEHYRVAVYILAWTSLRFGELIELRRKDIDDDGETMMFRVRRGAARVGQKVV"
            "VGNTKTVRSKRPVTVPPHVATMIREHMADRTKMNKGPEALLVTTTQGQRLSKSAFT"
            "RALKKGYRKIGRTDLRVHDLRAVGATYAAQAGATTKELMVRLGHTTPRMAMKYQMA"
            "SEARDAEIAKRMSELAGA")

        cds1 = test_data_utils.get_trixie_gene_data() # GeneID = "TRIXIE_0001"
        cds1["Translation"] = translation

        cds2 = test_data_utils.get_trixie_gene_data()
        cds2["Translation"] = "A"
        cds2["GeneID"] = "TRIXIE_0002"

        cds3 = test_data_utils.get_trixie_gene_data()
        cds3["Translation"] = translation
        cds3["GeneID"] = "TRIXIE_0003"

        test_db_utils.insert_gene_data(cds1)
        test_db_utils.insert_gene_data(cds2)
        test_db_utils.insert_gene_data(cds3)
        stmt = "UPDATE gene SET DomainStatus = 1"
        test_db_utils.execute(stmt)
        self.engine = sqlalchemy.create_engine(engine_string, echo=False)
        self.unparsed_args = get_unparsed_args()


    def tearDown(self):
        shutil.rmtree(test_folder)
        test_db_utils.remove_db()
        self.engine.dispose()




    # @patch("pdm_utils.functions.mysqldb.connect_to_db")
    # def test_main_1(self, ctd_mock):
    #     """Verify one gene is processed."""
    #     logging.info("test_main_1")
    #     ctd_mock.return_value = self.engine
    #
    #     stmt = "UPDATE gene SET DomainStatus = 0 WHERE GeneID = 'TRIXIE_0001'"
    #     test_db_utils.execute(stmt)
    #     run.main(self.unparsed_args)
    #
    #     gene_domain_table_results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
    #     domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
    #     gene_table_results = test_db_utils.get_data(test_db_utils.gene_table_query)
    #
    #     gene_table_dict = get_gene_id_dict(gene_table_results)
    #     domain_status = gene_table_dict["TRIXIE_0001"][0]["DomainStatus"]
    #     d_rows = len(domain_table_results)
    #     gd_rows = len(gene_domain_table_results)
    #
    #     with self.subTest():
    #         self.assertTrue(d_rows > 0)
    #     with self.subTest():
    #         self.assertTrue(gd_rows > 0)
    #     with self.subTest():
    #         self.assertEqual(domain_status, 1)


    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    def test_main_2(self, ctd_mock):
        """Verify no genes are processed if all DomainStatus = 1."""
        logging.info("test_main_2")
        ctd_mock.return_value = self.engine
        run.main(self.unparsed_args)
        gene_domain_table_results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        gene_table_results = test_db_utils.get_data(test_db_utils.gene_table_query)
        d_rows = len(domain_table_results)
        gd_rows = len(gene_domain_table_results)
        with self.subTest():
            self.assertEqual(d_rows, 0)
        with self.subTest():
            self.assertEqual(gd_rows, 0)


    @patch("pdm_utils.functions.mysqldb.connect_to_db")
    def test_main_3(self, ctd_mock):
        """Verify all three genes processed if DomainStatus = 0."""
        logging.info("test_main_3")
        ctd_mock.return_value = self.engine

        stmt = "UPDATE gene SET DomainStatus = 0"
        test_db_utils.execute(stmt)
        run.main(self.unparsed_args)

        gene_domain_table_results = test_db_utils.get_data(test_db_utils.gene_domain_table_query)
        domain_table_results = test_db_utils.get_data(test_db_utils.domain_table_query)
        gene_table_results = test_db_utils.get_data(test_db_utils.gene_table_query)

        gene_table_dict = get_gene_id_dict(gene_table_results)
        gene_domain_table_dict = get_gene_id_dict(gene_domain_table_results)

        domain_status1 = gene_table_dict["TRIXIE_0001"][0]["DomainStatus"]
        domain_status2 = gene_table_dict["TRIXIE_0002"][0]["DomainStatus"]
        domain_status3 = gene_table_dict["TRIXIE_0003"][0]["DomainStatus"]
        d_rows = len(domain_table_results)
        gd_rows1 = len(gene_domain_table_dict["TRIXIE_0001"])
        gd_rows3 = len(gene_domain_table_dict["TRIXIE_0003"])

        with self.subTest():
            self.assertTrue(d_rows > 0)
        with self.subTest():
            self.assertTrue(gd_rows1 > 0)
        with self.subTest():
            self.assertTrue("TRIXIE_0002" not in gene_domain_table_dict.keys())
        with self.subTest():
            self.assertEqual(gd_rows1, gd_rows3)
        with self.subTest():
            self.assertEqual(domain_status1, 1)
        with self.subTest():
            self.assertEqual(domain_status2, 1)
        with self.subTest():
            self.assertEqual(domain_status3, 1)



if __name__ == '__main__':
    unittest.main()
