"""Object to provide a formatted filtering query 
for retrieving data from a SQL database."""

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation
from pdm_utils.functions import phamerator
from pdm_utils.classes import genome, cds, mysqlconnectionhandler
from typing import List
import cmd, readline, argparse, os, sys


#Global file constants


class Filter:
    def __init__(self, database_name, sql_handle, phage_id_list=[]):
        """Initializes a Filter object used to filter
        results from a SQL database
        """
        self.database = database_name
        self.sql_handle = sql_handle

        self.phage_query = "SELECT Name FROM phage"
        self.gene_query = "SELECT phageID FROM gene"
       
        self.results = []
                
        self.phage_filters = {}
        self.gene_filters = {}

        self.history = []
    
        self.phage_attributes = []
        for column in self.sql_handle.execute_query(query="DESCRIBE phage"):
            self.phage_attributes.append(column["Field"].lower())

        self.gene_attributes = []
        for column in self.sql_handle.execute_query(query="DESCRIBE gene"):
            self.gene_attributes.append(column["Field"].lower())

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

        if self.phage_filters: 
            self.phage_query = "SELECT Name FROM phage WHERE " +\
                               " and ".join(list(self.phage_filters.values()))
        else:
            self.phage_query = "SELECT Name FROM phage"

        if self.gene_filters: 
            self.gene_query = "SELECT phageID FROM gene WHERE " +\
                              " and ".join(list(self.gene_filters.values()))
        else:
            self.gene_query = "SELECT phageID FROM gene"

        if self.phage_filters or self.gene_filters: 
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

    def hits(self): 

        self.retrieve_results()
        return len(self.results)

    def retrieve_results(self, verbose = False):
       
        self.update()
        if self.phage_filters:
            phage_results = self.retrieve_phage_results()
        
        if self.gene_filters:
            gene_results = self.retrieve_gene_results()

        if self.phage_filters and self.gene_filters:
            self.results = set(phage_results).intersection(gene_results)
        elif self.phage_filters:
            self.results = phage_results

        elif self.gene_filters:
            self.results = gene_results
            
    def retrieve_phage_results(self):
        """Helper function to retrieve phage table results"""

        database_results = [] 
        for result in phamerator.retrieve_data(
                        self.sql_handle, query= self.phage_query):
            database_results.append(result['Name'])

        return database_results

    def retrieve_gene_results(self):
        """Helper function to retrieve gene table results"""
        database_gene_results = []
        for result in phamerator.retrieve_data(
                        self.sql_handle, query = self.gene_query):
            database_gene_results.append(result['phageID'])

        database_results = []
        for gene_result in database_gene_results:
            if gene_result not in database_results:
                database_results.append(gene_result)

        return database_results
 
    def accession(self, filter_value: str):
        self.phage_filters.update({"accession" :\
                                   "Accession = '{}'".format(filter_value)})

    def name(self, filter_value: str):

        self.phage_filters.update({"name" :\
                                   "Name = '{}'".format(filter_value)})

    def phageID(self, filter_value: str):
        self.phage_filters.update({"id":\
                                   "PhageID = '{}'".format(filter_value)})

    def cluster(self, filter_value: str):
        self.phage_filters.update({"cluster":\
            "Cluster2 = '{}'".format(filter_value)}) 

    def subcluster(self, filter_value: str):
        self.phage_filters.update({"subcluster":\
            "Subcluster2 = '{}'".format(filter_value)})

    def status(self, filter_value: str):
        self.phage_filters.update({"status":\
            "status = '{}'".format(filter_value)})

    def retrieve_record(self, filter_value: str):
        self.phage_filters.update({"retrieve":\
            "RetrieveRecord = '{}'".format(filter_value)})

    def annotation_qc(self, filter_value: str):
        self.phage_filters.update({"annotation_qc":\
            "AnnotationQC = '{}'".format(filter_value)})

    def annotation_author(self, filter_value: str):
        self.phage_filters.update({"annotation_author":\
            "AnnotationAuthor = '{}'".format(filter_value)})

    def gene_product(self, filter_value: str):
        self.gene_filters.update({"product":\
            "Notes = '{}'".format(filter_value)})

    def gene_id(self, filter_value: str):
        self.gene_filters.update({"id":\
            "id  = '{}'".format(filter_value)})

    def add_filter(self, filter, filter_value, gene_selection=False):

        if gene_selection == False: 
            if filter.lower() in self.phage_attributes:
                self.phage_filters.update(
                        {filter: "{} = '{}'".format(filter, filter_value)})
                return True

            elif filter.lower() in self.gene_attributes:
                self.gene_filters.update(
                        {filter: "{} = '{}'".format(filter, filter_value)})
                return True

            else:
                return False

        else:
            if filter.lower() in self.gene_attributes:
                self.gene_filters.update(
                        {filter: "{} = '{}'".format(filter, filter_value)})
                return True

            elif filter.lower() in self.phage_attributes:
                self.phage_filters.update(
                        {filter: "{} = '{}'".format(filter, filter_value)})
                return True

            else:
                return False


    def interactive(self, sql_handle = None):
        
        interactive_filter = Cmd_Filter(db_filter=self, sql_handle=sql_handle)
        interactive_filter.cmdloop()

class Cmd_Filter(cmd.Cmd):
       
    def __init__(self, db_filter=None, sql_handle=None):
        super(Cmd_Filter, self).__init__()
        self.filter = db_filter
        self.sql_handle = sql_handle
        self.intro =\
        """---------------Database Search and Query---------------
        Type help or ? to list commands.\n"""
        self.prompt = "(database) (filter)user@localhost: "
        self.data = None

    def preloop(self):
        
        if self.filter == None:
            if self.sql_handle == None:
                self.sql_handle = \
                        mysqlconnectionhandler.MySQLConnectionHandler()
                print(\
                    "---------------------Database Login---------------------")
                self.sql_handle.database = input("MySQL database: ")
                self.sql_handle.get_credentials()
                try:
                    self.sql_handle.validate_credentials()
                except:
                    print("Unable to create a mysql connection.")
                    exit(1) 
            self.filter = Filter(self.sql_handle.database, self.sql_handle)

        self.prompt = "({}) (export){}@localhost: ".\
                format(self.sql_handle.database, self.sql_handle._username)
            
        

    def do_filter(self, *unparsed_args):
        """
        Formats a MySQL database query for a given attribute
        FILTER OPTIONS: Accession, Name, Cluster
        Subcluster, Annotation_Status, Retrieve_Record,
        Annotation_QC, Annotation_Author, Notes, GeneID
        """
        
        filter = unparsed_args[0]

        if filter == "":
            print(""" Please select a valid attribute:
            FILTER OPTIONS: Accession, Name, Cluster
            Subcluster, Annotation_Status, Retrieve_Record,
            Annotation_QC, Annotation_Author, Notes, GeneID
            """)

        elif not (filter.lower() in self.filter.phage_attributes\
                or filter.lower() in self.filter.gene_attributes):
            
            print("Attribute not in database.\n ")

        else:
            
            filter_value = input("Enter {} value: ".format(filter))
            self.filter.add_filter(filter, filter_value)
            self.retrieve_hits()

    def retrieve_hits(self, *args):
        "Function to retrieve the hits for filtering functions"
        hits = self.filter.hits()
        if hits <= 0:
            print("        No Database Hits.")
            print("        Reloading last filtering options...\n")
            self.filter.undo()
            self.filter.retrieve_results()
        else:
            print("        Database Hits: {}\n".format(hits))

    def do_undo(self, *args):
        """Reverts back to queries generated by previous filters
        USAGE: undo
        """
        self.filter.undo()
        print("        Reloaded last filtering options")
            
        self.do_show_filters()

    def do_show(self,*args):
        """Displays information on current filtering
        SHOW OPTIONS: Results, Filters, Hits"""

        options = ["results", "filters", "hits"]
        option = args[0].lower()

        if option in options:
            if option == "results":
                self.show_results()
            elif option == "filters":
                self.show_filters()
            elif option == "hits":
                self.retrieve_hits()

        else:
            print("""Please input a valid option
            SHOW OPTIONS: Results, Filter,Hits""")

    def show_results(self):
        """Shows results for current database filtering
        """
       
        print("        Results:\n")
        for row in range(0, len(self.filter.results), 3):
            result_row = self.filter.results[row:row+3]
            if len(result_row) == 3:
                print("%-20s %-20s %s" % \
                     (result_row[0], result_row[1], result_row[2]))
            elif len(result_row) == 2:
                print("%-20s %-20s" % \
                     (result_row[0], result_row[1]))

            elif len(result_row) == 1:
                print("%s" % (result_row[0]))

        print("\n")

    def show_filters(self, *args):
        """Displays current fitlers
        """
        if self.filter.phage_filters:
            print("        Phage table filters:")
            for phage_filter in self.filter.phage_filters.values():
                print("    " + phage_filter)
            print("\n")

        if self.filter.gene_filters:
            print("        Gene table filters:")
            for gene_filter in self.filter.gene_filters.values():
                print("    " + gene_filter)
            print("\n")

        if not self.filter.phage_filters and not self.filter.gene_filters:
            print("        No current filters applied.")

    def do_reset(self, *args):
        """Resets results history and current filters
        USAGE: reset
        """
        print("        Resetting Filters and Results History...\n")
        self.filter.reset()

    def do_clear(self, *args):
        """Clears display terminal
        USAGE: clear
        """

        os.system('cls' if os.name == 'nt' else 'clear')
        print(self.intro)

    def do_return(self, *args):
        """Finish quering and return results
        USAGE: return
        """
        
        print("        Exiting Filtering...\n")

        return True
   
    def do_exit(self, *args):
        """Exits program entirely without returning values
        USAGE: exit
        """

        sys.exit(1)

    def postloop(self):

        self.data = self.filter.results

if __name__ == "__main__":
    mainloop = Cmd_Filter() 
    mainloop.cmdloop()
