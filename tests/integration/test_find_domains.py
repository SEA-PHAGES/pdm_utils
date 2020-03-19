"""Integration tests for the main find_domains pipeline."""

from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import sqlalchemy
from pdm_utils.pipelines import find_domains

# Import helper functions to build mock database and mock flat files
unittest_file = Path(__file__)
test_dir = unittest_file.parent.parent
if str(test_dir) not in set(sys.path):
    sys.path.append(str(test_dir))
import pdm_utils_mock_db
import pdm_utils_mock_data

# The following integration tests user the 'pdm_anon' MySQL user.
# It is expected that this user has all privileges for 'test_db' database.
engine_string = pdm_utils_mock_db.create_engine_string()

# Assumes that output message contains "SQLAlchemy Error..."
error_msg = "SQLAlchemy"

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

# INSERT IGNORE INTO domain (HitID, DomainID, Name, Description) VALUES ("gnl|CDD|334841", "pfam02195", "ParBc", "ParB-like nuclease domain.")
# INSERT IGNORE INTO gene_domain (GeneID, HitID, Expect, QueryStart, QueryEnd) VALUES ("SEA_ABBA_66", "gnl|CDD|334841", 1.78531e-11, 33, 115)

# INSERT_INTO_DOMAIN = """INSERT INTO domain (HitID, DomainID, Name, Description) VALUES ("{}", "{}", "{}", "{}")"""
# INSERT_INTO_GENE_DOMAIN = """INSERT INTO gene_domain (GeneID, HitID, Expect, QueryStart, QueryEnd) VALUES ("{}", "{}", {}, {}, {})"""
# UPDATE_GENE = "UPDATE gene SET DomainStatus = 1 WHERE GeneID = '{}'"

class TestFindDomains1(unittest.TestCase):

    def setUp(self):
        pdm_utils_mock_db.create_empty_test_db()
        pdm_utils_mock_db.insert_phage_data(pdm_utils_mock_data.get_trixie_phage_data())
        pdm_utils_mock_db.insert_gene_data(pdm_utils_mock_data.get_trixie_gene_data())

        self.engine = sqlalchemy.create_engine(engine_string, echo=False)
        self.connection = self.engine.connect()
        self.trans = self.connection.begin()


    def tearDown(self):
        pdm_utils_mock_db.remove_db()
        self.trans.rollback()
        self.engine.dispose()




    def test_execute_statement_1(self):
        """Verify valid DomainStatus update data can be inserted
        into gene table."""
        gene_table_results1 = pdm_utils_mock_db.get_data(pdm_utils_mock_db.gene_table_query)
        update_data = get_trixie_gene_table_domain_status_update_data_1()
        statement = get_gene_update_statement(update_data)
        result, type_error, msg = find_domains.execute_statement(
                                    self.connection, statement)
        self.trans.commit()
        phage_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.phage_table_query)
        gene_table_results2 = pdm_utils_mock_db.get_data(pdm_utils_mock_db.gene_table_query)
        gene_domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.gene_domain_table_query)
        domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
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
            self.assertIsInstance(msg, str)


    def test_execute_statement_2(self):
        """Verify valid data can be inserted into domain table."""
        domain_data = pdm_utils_mock_data.get_trixie_domain_data()
        statement = get_domain_insert_statement(domain_data)
        result, type_error, msg = find_domains.execute_statement(
                                    self.connection, statement)
        self.trans.commit()
        domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 1)
        with self.subTest():
            self.assertEqual(result, 0)
        with self.subTest():
            self.assertFalse(type_error)


    def test_execute_statement_3(self):
        """Verify valid data can be inserted into domain and gene_domain tables."""
        domain_data = pdm_utils_mock_data.get_trixie_domain_data()
        statement1 = get_domain_insert_statement(domain_data)
        result1, type_error1, msg1 = find_domains.execute_statement(
                                    self.connection, statement1)

        gene_domain_data = pdm_utils_mock_data.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        result2, type_error2, msg2 = find_domains.execute_statement(
                                    self.connection, statement2)

        self.trans.commit()
        gene_domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.gene_domain_table_query)
        domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
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


    def test_execute_statement_4(self):
        """Verify invalid data can NOT be inserted into gene_domain table
        when there is no matching 'HitID' in the domain table."""
        gene_domain_data = pdm_utils_mock_data.get_trixie_gene_domain_data()
        statement = get_gene_domain_insert_statement(gene_domain_data)
        result, type_error, msg = find_domains.execute_statement(
                                    self.connection, statement)
        self.trans.commit()
        gene_domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.gene_domain_table_query)
        with self.subTest():
            self.assertEqual(len(gene_domain_table_results), 0)
        with self.subTest():
            self.assertEqual(result, 1)
        with self.subTest():
            self.assertFalse(type_error)


    def test_execute_statement_5(self):
        """Verify invalid data can NOT be inserted into domain table when
        there is a duplicated HitID."""
        domain_data = pdm_utils_mock_data.get_trixie_domain_data()
        statement = get_domain_insert_statement(domain_data)
        result1, type_error1, msg1 = find_domains.execute_statement(
                                    self.connection, statement)
        self.trans.commit()
        domain_table_results1 = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
        new_trans = self.connection.begin()
        result2, type_error2, msg2 = find_domains.execute_statement(
                                    self.connection, statement)
        new_trans.commit()
        domain_table_results2 = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
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
            self.assertFalse(error_msg in msg1)
        with self.subTest():
            self.assertTrue(error_msg in msg2)


    def test_execute_statement_6(self):
        """Verify invalid data can NOT be inserted into domain table when
        there is an invalid column name."""
        domain_data = pdm_utils_mock_data.get_trixie_domain_data()
        statement = get_domain_insert_statement(domain_data)
        statement = statement.replace("Name", "Name_invalid")
        result, type_error, msg = find_domains.execute_statement(
                                    self.connection, statement)
        self.trans.commit()
        domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 0)
        with self.subTest():
            self.assertEqual(result, 1)
        with self.subTest():
            self.assertTrue(error_msg in msg)
        with self.subTest():
            self.assertFalse(type_error)


    def test_execute_statement_7(self):
        """Verify invalid data can NOT be inserted due to '%'."""
        domain_data = pdm_utils_mock_data.get_trixie_domain_data()
        # "Description": "ParB-like nuclease domain"
        description = domain_data["Description"]
        description = description.replace("nuclease domain", "nuclease % domain")
        domain_data["Description"] = description
        statement = get_domain_insert_statement(domain_data)
        result, type_error, msg = find_domains.execute_statement(
                                    self.connection, statement)
        self.trans.commit()
        domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 0)
        with self.subTest():
            self.assertEqual(result, 1)
        with self.subTest():
            self.assertFalse(error_msg in msg)
        with self.subTest():
            self.assertTrue(type_error)


    def test_execute_statement_8(self):
        """Verify invalid data can be inserted after '%' is
        replaced with '%%'."""
        domain_data = pdm_utils_mock_data.get_trixie_domain_data()
        # "Description": "ParB-like nuclease domain"
        description = domain_data["Description"]
        description = description.replace("nuclease domain", "nuclease % domain")
        domain_data["Description"] = description
        statement = get_domain_insert_statement(domain_data)
        statement = statement.replace("%", "%%")
        result, type_error, msg = find_domains.execute_statement(
                                    self.connection, statement)
        self.trans.commit()
        domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 1)
        with self.subTest():
            self.assertEqual(result, 0)
        with self.subTest():
            self.assertFalse(type_error)




class TestFindDomains2(unittest.TestCase):

    def setUp(self):
        pdm_utils_mock_db.create_empty_test_db()
        pdm_utils_mock_db.insert_phage_data(pdm_utils_mock_data.get_trixie_phage_data())
        pdm_utils_mock_db.insert_gene_data(pdm_utils_mock_data.get_trixie_gene_data())

        self.engine = sqlalchemy.create_engine(engine_string, echo=False)
        self.connection = self.engine.connect()


    def tearDown(self):
        pdm_utils_mock_db.remove_db()
        self.engine.dispose()




    def test_execute_transaction_1(self):
        """Verify function runs with list of zero statements."""
        result = find_domains.execute_transaction(self.connection)
        domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 0)
        with self.subTest():
            self.assertEqual(result, 0)


    def test_execute_transaction_2(self):
        """Verify list of one valid statement can be inserted into domain table."""
        domain_data = pdm_utils_mock_data.get_trixie_domain_data()
        statement = get_domain_insert_statement(domain_data)
        statements = [statement]
        result = find_domains.execute_transaction(self.connection, statements)
        domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 1)
        with self.subTest():
            self.assertEqual(result, 0)


    def test_execute_transaction_3(self):
        """Verify list of two valid statements can be inserted into
        domain and gene_domain tables"""
        domain_data = pdm_utils_mock_data.get_trixie_domain_data()
        statement1 = get_domain_insert_statement(domain_data)
        gene_domain_data = pdm_utils_mock_data.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        statements = [statement1, statement2]
        result = find_domains.execute_transaction(self.connection, statements)
        gene_domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.gene_domain_table_query)
        domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 1)
        with self.subTest():
            self.assertEqual(len(gene_domain_table_results), 1)
        with self.subTest():
            self.assertEqual(result, 0)


    def test_execute_transaction_4(self):
        """Verify list of three statements (including one with
        duplicated HitID) are inserted."""
        domain_data1 = pdm_utils_mock_data.get_trixie_domain_data()
        pdm_utils_mock_db.insert_domain_data(domain_data1)
        domain_table_results1 = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
        # Duplicate HitID
        statement1 = get_domain_insert_statement(domain_data1)
        # Valid
        gene_domain_data = pdm_utils_mock_data.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        # Valid
        update_data = get_trixie_gene_table_domain_status_update_data_1()
        statement3 = get_gene_update_statement(update_data)

        statements = [statement1, statement2, statement3]
        result = find_domains.execute_transaction(self.connection, statements)
        gene_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.gene_table_query)
        gene_domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.gene_domain_table_query)
        domain_table_results2 = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
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
        domain_data1 = pdm_utils_mock_data.get_trixie_domain_data()
        statement1 = get_domain_insert_statement(domain_data1)
        # Valid
        gene_domain_data = pdm_utils_mock_data.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        # Invalid
        domain_data2 = pdm_utils_mock_data.get_trixie_domain_data()
        statement3 = get_domain_insert_statement(domain_data2)
        statement3 = statement3.replace("HitID", "unique_id")
        statement3 = statement3.replace("Name", "Name_invalid")
        # Valid - function should exit before executing this though.
        update_data = get_trixie_gene_table_domain_status_update_data_1()
        statement4 = get_gene_update_statement(update_data)

        statements = [statement1, statement2, statement3, statement4]
        result = find_domains.execute_transaction(self.connection, statements)
        gene_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.gene_table_query)
        gene_domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.gene_domain_table_query)
        domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
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
        domain_data1 = pdm_utils_mock_data.get_trixie_domain_data()
        statement1 = get_domain_insert_statement(domain_data1)
        # Valid
        gene_domain_data = pdm_utils_mock_data.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)
        # Invalid '%'
        domain_data2 = pdm_utils_mock_data.get_trixie_domain_data()
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
        gene_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.gene_table_query)
        gene_domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.gene_domain_table_query)
        domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
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
    def test_execute_transaction_7(self, es_mock):
        """Verify error inside try/except block is executed."""
        stmt_result1 = 0
        type_error1 = False
        msg1 = "empty"
        mock_result1 = (stmt_result1, type_error1, msg1)

        stmt_result2 = 0
        type_error2 = False
        msg2 = 2 # the function expects this to be a string, so this should
                 # break the code and trigger the except block.
        mock_result2 = (stmt_result2, type_error2, msg2)
        es_mock.side_effect = [mock_result1, mock_result2]
        # Valid
        domain_data1 = pdm_utils_mock_data.get_trixie_domain_data()
        statement1 = get_domain_insert_statement(domain_data1)
        # Valid
        gene_domain_data = pdm_utils_mock_data.get_trixie_gene_domain_data()
        statement2 = get_gene_domain_insert_statement(gene_domain_data)

        statements = [statement1, statement2]
        result = find_domains.execute_transaction(self.connection, statements)
        gene_domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.gene_domain_table_query)
        domain_table_results = pdm_utils_mock_db.get_data(pdm_utils_mock_db.domain_table_query)
        with self.subTest():
            self.assertEqual(len(domain_table_results), 0)
        with self.subTest():
            self.assertEqual(len(gene_domain_table_results), 0)
        with self.subTest():
            self.assertEqual(result, 1)
        with self.subTest():
            self.assertEqual(es_mock.call_count, 2)


if __name__ == '__main__':
    unittest.main()
