"""Object to provide a formatted filtering query 
for retrieving data from a SQL database."""

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation
from pdm_utils.functions import phamerator
from pdm_utils.classes import genome, cds, mysqlconnectionhandler
from typing import List
from pprint import pprint
import cmd, readline, argparse, os, sys

class Filter:
    def __init__(self, database_name, phage_id_list: List[str] = []):
        """Initializes a Filter object used to filter
        results from a SQL database
        """
        self.database = database_name
        
        self.phage_query = "SELECT Name FROM phage"
        self.gene_query = "SELECT phageID FROM gene"
        
        self.results = []
                
        self.phage_filters = {}
        self.gene_filters = {}

        self.history = []

    def reset(self):
        """Resets created queries and filters 
        for the Filter object"""

        self.phage_query = "SELECT Name FROM phage"
        self.gene_query = "SELECT phageID FROM gene"
        
        self.results = []
                
        self.phage_filters = {}
        self.gene_filters = {}

        self.history = []
         
    def update(self):
        """Resets created queries for the Filter object"""

        
        self.phage_query = "SELECT Name FROM phage " +\
                " and ".join(list(self.phage_filters.values()))

        self.gene_query = "SELECT phageID FROM gene " +\
                " and ".join(list(self.gene_filters.values()))

        self.history.append((self.phage_filters, self.gene_filters))
                    
    def undo(self):
        """Undos last filter option"""

        if self.history:
            dead_filters = self.history.pop()
        
        if self.history:
            last_filters = self.history.pop()

            self.phage_filters = last_filters[0]

            self.gene_filters = last_filters[1]

        else:
            self.reset()

    def establish_connection(self):
        """Creates a mysqlconnectionhandler object 
        and populates its credentials
        """

        sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
        sql_handle.database = self.database
        sql_handle.get_credentials()
        try:
            sql_handle.validate_credentials
        except:
            print("Sql connection to database {} \
                    with username  and password failed".format(database_name))

        return sql_handle 

    def hits(self, sql_handle): 

        self.retrieve_results(sql_handle)
        return len(self.results)

    def retrieve_results(self, sql_handle, visual = False):
       
        self.update()
        if self.phage_filters and self.gene_filters:
            database_phage_results = [] 
            for result in phamerator.retrieve_data\
                    (sql_handle, query = self.phage_query):
                database_phage_results.append(result['Name'])
            database_gene_results = []
            for result in phamerator.retrieve_data\
                    (sql_handle, query = self.gene_query):
                database_gene_results.append(result['PhageID'])
            self.results = set(database_phage_results).\
                    intersection(database_gene_results)

        elif self.gene_filters:

            database_gene_results = []
            for result in phamerator.retrieve_data\
                    (sql_handle, query = self.gene_query):
                database_gene_results.append(result['phageID'])
            self.results = database_gene_results

        elif self.phage_filters:
            
            database_phage_results = []
            for result in phamerator.retrieve_data\
                    (sql_handle, query = self.phage_query):
                database_phage_results.append(result['Name'])
            self.results = database_phage_results

        else:
            database_results = []
            for result in phamerator.retrieve_data\
                    (sql_handle, query = self.phage_query):
                database_results.append(result['Name'])
            self.results = database_results
            
    def accession(self, filter: str):
        self.phage_filters.update({"accession" :\
            "WHERE Accession = '{}'".format(filter)})

    def name(self, filter: str):

        self.phage_filters.update({"name" :\
            "WHERE Name = '{}'".format(filter)})

    def phageID(self, filter: str):
        self.phage_filters.update({"id":\
            "WHERE PhageID = '{}'".format(filter)})

    def cluster(self, filter: str):
        self.phage_filters.update({"cluster":\
            "WHERE Cluster2 = '{}'".format(filter)}) 

    def subcluster(self, filter: str):
        self.phage_filters.update({"subcluster":\
            "WHERE Subcluster2 = '{}'".format(filter)})

    def status(self, filter: str):
        self.phage_filters.update({"status":\
            "WHERE status = '{}'".format(filter)})

    def retrieve_record(self, filter: str):
        self.phage_filters.update({"retrieve":\
            "WHERE RetrieveRecord = '{}'".format(filter)})

    def annotation_qc(self, filter: str):
        self.phage_filters.update({"annotation_qc":\
            "WHERE AnnotationQC = '{}'".format(filter)})

    def annotation_author(self, filter: str):
        self.phage_filters.update({"annotation_author":\
            "WHERE AnnotationAuthor = '{}'".format(filter)})

    def gene_product(self, filter: str):
        self.gene_filters.update({"product":\
            "WHERE Notes = '{}'".format(filter)})

    def gene_id(self, filter: str):
        self.gene_filters.update({"id":\
            "WHERE id  = '{}'".format(filter)})

    def interactive(self, sql_handle = None):
        
        interactive_filter = \
                Interactive_Filter(db_filter = self, sql_handle = sql_handle)
        interactive_filter.cmdloop()

class Interactive_Filter(cmd.Cmd):
       
    def __init__(self, db_filter = None, sql_handle = None):
        super(Interactive_Filter, self).__init__()
        self.filter = db_filter
        self.sql_handle = sql_handle
        self.intro =\
        """---------------Hatfull Helper's Filtering---------------
        Type help or ? to list commands.\n"""
        self.prompt = "(database) (filter)user@localhost: "
        self.data = None

    def preloop(self, db_filter = None, sql_handle = None):
                
        if self.filter == None:
            print("---------------------Database Login---------------------")
            self.filter = Filter(input("MySQL database: "))

        
        if self.sql_handle == None:
            self.sql_handle =\
                    self.filter.establish_connection()
        self.prompt = "({}) (filter){}@localhost: ".format\
                (self.filter.database, self.sql_handle._username)
         
    def do_accession(self, *args):
        """Search for genomes by accession number.
        USAGE: accession
        """
        self.filter.accession(input("Accession: "))
        self.retrieve_hits()

    def do_name(self, *args):
        """Search for genomes by name.
        USAGE: name
        """
        self.filter.name(input("Name: "))
        self.retrieve_hits()

    def do_cluster(self, *args):
        """Search for genomes by cluster.
        USAGE: cluster
        """
        self.filter.cluster(input("Cluster: "))
        self.retrieve_hits()

    def do_subcluster(self, *args):
        """Search for genomes by subcluster.
        USAGE: subcluster
        """
        self.filter.subcluster(input("Subcluster: "))
        self.retrieve_hits()

    def do_status(self, *args):
        """Search for genomes by annotation status.
        USAGE: Status :Status:
        """
        self.filter.status(input("Status: "))
        self.retrieve_hits()

    def do_retrieve_record(self, *args):
        """Search for genomes by retrieve record value.
        USAGE: retrieve_record
        """
        self.filter.retrieve_record(input("Retrieve Record: "))
        self.retrieve_hits()

    def do_annotation_QC(self, *args):
        """Search for genomes by annotation quality control.
        USAGE: annotation_QC
        """
        self.filter.annotation_qc(input("Annotation QC: "))
        self.retrieve_hits()

    def do_annotation_author(self, *args):
        """Search for genomes by annotation authorship.
        USAGE: annotation_author
        """
        self.filter.annotation_author(input("Annotation Author: "))
        self.retrieve_hits()

    def do_gene_product(self, *args):
        """Search for genomes by gene products.
        USAGE: gene_product
        """
        self.filter.gene_product(input("Gene Product: ") )
        self.retrieve_hits()

    def do_gene_ID(self, *args):
        """Search for genomes by gene IDs.
        USAGE: gene_ID
        """
        self.filter.gene_id(input("Gene ID: "))
        self.retrieve_hits()

    def retrieve_hits(self, *args):
        "Function to retrieve the hits for filtering functions"
        hits = self.filter.hits(self.sql_handle)
        if hits <= 0:
            print("\
                    No Database Hits.")
            print("\
                    Reloading last filtering options...\n")
            self.filter.undo()
            self.filter.retrieve_results(self.sql_handle)
        else:
            print("\
                    Database Hits: {}\n".format(hits))

    def do_undo(self, *args):
        """Reverts back to queries generated by previous filters
        USAGE: undo
        """
        self.filter.undo()
        print("\
                Reloaded last filtering options")
        print("\
                Database Hits: {}\n".format(len(self.filter.results)))

    def do_show_results(self, *args):
        """Shows results for current database filtering
        USAGE: show_results
        """
       
        print("\
                Results:")
        for row in range(0, len(self.filter.results), 3):
            pprint(self.filter.results[row:row+3])
        print("\n")

    def do_show_filters(self, *args):
        """Displays current fitlers
        USAGE: shpw_filters
        """
        if self.filter.phage_filters:
            print("\
                    Phage table filters:")
            pprint(self.filter.phage_filters)
            print("\n")

        if self.filter.gene_filters:
            print("\
                    Gene table filters:")
            pprint(self.filter.gene_filters)
            print("\n")

        if not self.filter.phage_filters and not self.filter.gene_filters:
            print("\
                    No current filters applied.")

    def do_reset(self, *args):
        """Resets results history and current filters
        USAGE: reset
        """
        print("\
                Resetting Filters and Results History...\n")
        self.filter.reset()

    def do_clear(self, *args):
        """Clears display terminal
        USAGE: clear
        """

        os.system('cls' if os.name == 'nt' else 'clear')

    def do_return(self, *args):
        """Finish quering and return results
        USAGE: return
        """
        
        print("\
                Exiting Filtering...\n")

        return True
   
    def do_exit(self, *args):
        """Exits program entirely without returning values
        USAGE: exit
        """

        sys.exit(1)

    def postloop(self):

        self.data = self.filter.results

if __name__ == "__main__":
    mainloop = Interactive_Filter() 
    mainloop.cmdloop()
