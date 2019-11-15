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
        #Returns information about the database
        table_dicts = sql_handle.execute_query(
            "SELECT DISTINCT TABLE_NAME FROM INFORMATION_SCHEMA.COLUMNS"
            "WHERE COLUMN_NAME IN ('PhageID') AND "
            "TABLE_SCHEMA='{self.database}'")
        tables = {}
        for data_dict in table_dicts:
            table = data_dict["TABLE_NAME"]
            column_dicts = sql_handle.execute_query(f"SHOW columns IN {table}")
            column_list = []
            for column_dict in column_dicts:
                column_list.append(column_dict["Field"])
            tables.update({table: column_list}) 
     
        if phage_id_list == []:
            self.results = phamerator.create_phage_id_set(self.sql_handle)
        else:
            self.results = phage_id_list
                
        self.history = []
    
    def reset(self):
        """
        Resets created queries and filters 
        for the Filter object
        """
         
    def update(self):
        """
        Resets created queries for the Filter object
        """

    def undo(self):
        """
        Undos last filter option
        """

    def hits(self): 
        """
        Returns length of current results
        """

    def retrieve_results(self, verbose = False):
        """
        Sets results according to current filters
        """

    def add_filter(self, verbose = False):
        """
        Adds a filter to database queries
        """

    def group(self):
        pass

    def sort(self):
        pass 

    def interactive(self, sql_handle = None):
        """
        Function to start interactive filtering 
        """
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
