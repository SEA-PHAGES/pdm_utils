"""Object to provide a formatted filtering query
for retrieving data from a SQL database."""

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation
from pdm_utils.functions import phamerator
from pdm_utils.classes import genome, cds, mysqlconnectionhandler
from typing import List
import cmd, readline, argparse, os, sys, re, string

#Global file constants
group_options = ["cluster", "subcluster",
                 "status",  "host",
                 "author", "record"]


class Filter:
    def __init__(self, database_name, sql_handle, phage_id_list=[]):
        """Initializes a Filter object used to filter
        results from a SQL database
        """
        self.database = database_name
        self.sql_handle = sql_handle
        #Returns information about the database
        table_dicts = sql_handle.execute_query(
            "SELECT DISTINCT TABLE_NAME FROM INFORMATION_SCHEMA.COLUMNS "
            "WHERE COLUMN_NAME IN ('PhageID') AND "
            f"TABLE_SCHEMA='{self.database}'")
        tables = {}
        for data_dict in table_dicts:
            table = data_dict["TABLE_NAME"]
            column_dicts = sql_handle.execute_query(f"SHOW columns IN {table}")
            column_list = []
            for column_dict in column_dicts:
                column_list.append((column_dict["Field"]).lower())
            tables.update({table.lower(): column_list})

        self.tables = tables

        if phage_id_list == []:
            results = list(phamerator.create_phage_id_set(self.sql_handle))
        else:
            results = phage_id_list

        self.history = []
        self.filters = {}
        self.phageIDs = results

    def add_filter(self, table, field, value, verbose=False):
        """
        Adds a filter to database queries
        """
        if field == "PhageID":
            if verbose:
                print("PhageID cannot be a field requested for filtering")
        elif table.lower() in (self.tables).keys():
            if field.lower() in self.tables[table]:
                if table.lower() in self.filters.keys():
                    if field.lower() in (self.filters[table]).keys():
                        (self.filters[table])[field].append(value)
                    else:
                        (self.filters[table])[field] = [value]
                else:
                    self.filters.update({table: {field: [value]}})
            else:
                if verbose:
                    print(
                    f"Field '{field}' requested to be filtered"
                    f" is not in table '{table}'")
        else:
            if verbose:
                print(
                f"Table '{table}' requested to be filtered "
                f"is not in '{self.database}'")

    def refresh(self):
        """
        Refreshes phageIDs
        """
        query = ("SELECT PhageID FROM phage WHERE PhageID IN " + \
                 "('" + "','".join(self.phageIDs) + "')")
        results_dicts = self.sql_handle.execute_query(query)
        results = []
        for result in results_dicts:
            results.append(result["PhageID"])

        self.phageIDs = results

    def update(self, verbose=False, sort=None):
        """
        Updates results list for the Filter object
        """
        query_results = self.phageIDs

        queries = []
        for table in self.filters.keys():
            for field in self.filters[table].keys():
                queries.append(
                    f"{field}='"+ \
                            "','".join((self.filters[table])[field]) + "'")
            query = (f"SELECT PhageID FROM {table} WHERE " +\
                    " and ".join(queries) + \
                    " and PhageID IN " + "('" + \
                            "','".join(query_results) + "')")
            if verbose:
                print(f"Filtering {table} in {self.database} for " +\
                       " and ".join(queries) + "...")
            query_results=[]
            for result in self.sql_handle.execute_query(query):
                query_results.append(result["PhageID"])
            queries = []


        if sort != None and sort.lower() in self.tables["phage"]:
            if verbose:
                print(f"Sorting by '{sort}'...")
            query = ("SELECT PhageID FROM phage "
                     "WHERE PhageID IN " + "('" + \
                            "','".join(query_results) + "') "
                    f"ORDER BY {sort}")
            query_results = []
            for result in self.sql_handle.execute_query(query):
                query_results.append(result["PhageID"])

        elif sort != None:
            if verbose:
                print(f"Sort key '{sort}' not supported.")

        self.phageIDs = query_results
        if len(self.phageIDs) > 0:
            self.history.append(query_results)

    def results(self, verbose=False):
        """
        Sets results according to current filters
        """
        if verbose:
            print("Results:")
            for row in range(0, len(self.phageIDs), 3):
                result_row = self.phageIDs[row:row+3]
                if len(result_row) == 3:
                    print("%-20s %-20s %s" % \
                         (result_row[0], result_row[1], result_row[2]))
                elif len(result_row) == 2:
                    print("%-20s %-20s" % \
                         (result_row[0], result_row[1]))

                elif len(result_row) == 1:
                    print("%s" % (result_row[0]))

        return self.phageIDs

    def undo(self, verbose=False):
        """
        Undos last filter option
        """
        if len(self.history) == 1:
            if verbose:
                print("No results history.")
        else:
            current = self.history.pop()
            self.phageIDs = self.history.pop()
            if verbose:
                print("Returned last valid results.")

    def reset(self, verbose=False):
        """
        Resets created queries and filters
        for the Filter object
        """
        self.phageIDs = list(phamerator.create_phage_id_set(self.sql_handle))
        self.history = []
        if verbose:
            print("Results and results history cleared.")

    def hits(self, verbose=False):
        """
        Returns length of current results
        """
        if verbose:
            print(f"Database hits: {len(self.phageIDs)}")
        return len(self.phageIDs)

    def group(self, group, verbose=False):
        """
        Function that creates a two-dimensional array of
        PhageIDs from the results list of PhageIDs
        separated by a characteristic.
        """
        qual_switch =       {"cluster"    : "Cluster2",
                             "subcluster" : "Subcluster2",
                             "status"     : "Status",
                             "host"       : "HostStrain"}

        qual_int_switch =   {"author"     : self.group_author,
                             "record"     : self.group_record}

        if group.lower() in qual_switch.keys():
            if verbose:
                print(f"Grouping by {group}...")
            return self.group_str_qualitative(qual_switch[group.lower()])

        else:
            if verbose:
                print(f"Group key '{group}' not supported.")
            return {}

    def group_str_qualitative(self, field_name):
        """
        Helper function that groups PhageIDs by phage
        field name into a two-dimensional array
        """
        group_dicts = self.sql_handle.execute_query(
            f"SELECT DISTINCT {field_name} FROM phage")
        group_set = []
        for group_dict in group_dicts:
            group_set.append(group_dict[f"{field_name}"])

        groups = dict.fromkeys(group_set, [])

        phage_list = self.phageIDs
        for group in group_set:
            query = f"SELECT PhageID FROM phage WHERE {field_name}='{group}'" + \
                     " and PhageID IN " + "('" + \
                     "','".join(phage_list) + "')"
            results_dicts = self.sql_handle.execute_query(query)
            if results_dicts != ():
                results_list = []
                for result_dict in results_dicts:
                     results_list.append(result_dict["PhageID"])
                     phage_list.remove(result_dict["PhageID"])
                groups[group] = results_list

        groups.update({"other": phage_list})
        return groups

    def group_author(self):
        """
        Helper function that groups PhageIDs by genome
        AnnotationAuthor value into a two-dimensional array
        """
        author_dicts = self.sql_handle.execute_query(
            "SELECT DISTINCT AnnotationAuthor FROM phage")
        author_set = []

        for author_dict in author_dicts:
            author_set.append(author_dict["AnnotationAuthor"])

        author_groups = dict.fromkeys(author_set, [])
        phage_list = self.phageIDs.copy()
        for author in author_set:
            query = f"SELECT PhageID FROM phage WHERE AnnotationAuthor='{author}'" + \
                     " and PhageID IN " + "('" + \
                     "','".join(phage_list) + "')"
            results_dicts = self.sql_handle.execute_query(query)
            if results_dicts != ():
                results_list = []
                for result_dict in results_dicts:
                     results_list.append(result_dict["PhageID"])
                     phage_list.remove(result_dict["PhageID"])
                author_groups[author] = results_list

        author_groups.update({"other": phage_list})
        return author_groups

    def group_record(self):
        """
        Helper function that groups PhageIDs by genome
        RetrieveRecord value into a two-dimensional array
        """
        record_dicts = self.sql_handle.execute_query(
            "SELECT DISTINCT RetrieveRecord FROM phage")
        record_set = []

        for record_dict in record_dicts:
            record_set.append(record_dict["RetrieveRecord"])

        record_groups = dict.fromkeys(record_set, [])

        phage_list = self.phageIDs.copy()
        for record in record_set:
            query = f"SELECT PhageID FROM phage WHERE RetrieveRecord='{record}'" + \
                     " and PhageID IN " + "('" + \
                     "','".join(self.phageIDs) + "')"
            results_dicts = self.sql_handle.execute_query(query)
            if results_dicts != ():
                results_list = []
                for result_dict in results_dicts:
                     results_list.append(result_dict["PhageID"])
                     phage_list.remove(result_dict["PhageID"])
                record_groups[record] = results_list

        record_groups.update({"other": phage_list})
        return record_groups

    def interactive(self, sql_handle = None):
        """
        Function to start interactive filtering.
        """
        pass

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
    sql_handle = mysqlconnectionhandler.MySQLConnectionHandler()
    sql_handle._username = "root"
    sql_handle._password = "phage"
    sql_handle.database = "Actino_Draft"
    filter = Filter(sql_handle.database, sql_handle)
    filter.add_filter("gene", "Notes", "antirepressor")
    filter.add_filter("phage", "cluster2", "F")
    filter.update()
