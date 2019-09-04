"""Integration tests for the main import pipeline."""



import unittest
from pipelines.db_import import import_main
from constants import constants
from classes import bundle, Genome, GenomePair, Ticket
from classes import MySQLConnectionHandler
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqFeature import ExactPosition, Reference



# The following integration tests user the 'tester' MySQL user.
# It is expected that this user has all privileges for 'test_db' database.
user = "tester"
pwd = "tester"
db = "test_db"


class TestImportMain(unittest.TestCase):


    def setUp(self):
        pass



    # TODO test not complete. Need to unittest evaluate functions first.
    def test_evaluate_flat_file_1(self):
        """."""

        ticket1 = Ticket.GenomeTicket()
        ticket1.id = 1
        ticket1.type = "replace"
        ticket1.phage_id = "Trixie"
        ticket1.host_genus = "Mycobacterium"
        ticket1.cluster = "A"
        ticket1.subcluster = "A2"
        ticket1.annotation_status = "final"
        ticket1.annotation_author = 1
        ticket1.annotation_qc = 1
        ticket1.retrieve_record = 1
        ticket1.description_field = "product"
        ticket1.run_mode = "phagesdb"
        ticket1.accession = "ABC123"

        genome1 = Genome.Genome()
        genome1.type = "flat_file"
        genome1.id = "Trixie"
        genome1.seq = Seq("ATGC")

        bndl = bundle.Bundle()
        bndl.ticket = ticket1
        bndl.genome_dict[genome1.type] = genome1

        sql_handle = MySQLConnectionHandler.MySQLConnectionHandler()
        sql_handle.username = user
        sql_handle.password = pwd
        sql_handle.database = db

        host_set = set(["Mycobacterium", "Gordonia"])
        phage_id_set = set(["L5", "Trixie", "D29"])
        seq_seq = set(["ATGC", "AAAA", "TTTT"])
        cluster_set = set(["A", "B", "C", "D"])
        subcluster_set = set(["A1", "A2", "B1"])






if __name__ == '__main__':
    unittest.main()
