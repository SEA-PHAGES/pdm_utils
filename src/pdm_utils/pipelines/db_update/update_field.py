# Import standard python modules
import argparse
import sys
import csv

# Import relevant k_phamerate objects
# from main_scripts.python3.classes.MySQLConnectionHandler import \
#     MySQLConnectionHandler
# from main_scripts.python3.classes.RandomFieldUpdateHandler import \
#     RandomFieldUpdateHandler
from pdm_utils.classes.mysqlconnectionhandler import MySQLConnectionHandler
from pdm_utils.classes.randomfieldupdatehandler import RandomFieldUpdateHandler

# Set up argparse
script_description = """
This is a script intended to handle specific, single-field updates to tables
in a Phamerator database.
"""
parser = argparse.ArgumentParser(description=script_description)
parser.add_argument("database_name", metavar="db", type=str, nargs=1,
                    help="name of the Phamerator database to be phamerated")
parser.add_argument("ticket_file", metavar="tf", type=str, nargs=1,
                    help="path to an update ticket")

# Parse command line arguments.
args = parser.parse_args()
database = args.database_name[0].split("=")[-1]
ticket_file = args.ticket_file[0].split("=")[-1]

# Establish the database connection using the MySQLConnectionHandler object
mysql_handler = MySQLConnectionHandler(database)
mysql_handler.create_connection()
mysql_connection = mysql_handler.connection

# Variables to be used for end summary
tickets_processed = 0
tickets_succeeded = 0
tickets_failed = 0

# Iterate through the tickets and process them sequentially.
f = open(ticket_file, "r")
ticket_reader = csv.reader(f, delimiter=",")
for row in ticket_reader:
    handler = RandomFieldUpdateHandler(mysql_connection)
    handler.table = row[0]        # Which table will be updated?
    handler.field = row[1]        # Which field will be updated?
    handler.value = row[2]        # What value will be put in that field?
    handler.key_name = row[3]   # How will we know which row is the right one?
    handler.key_value = row[4]  # How will we know which row is the right one?
    handler.validate_ticket()   # Make sure all handler attributes are valid
    status = handler.execute_ticket()    # Do what was requested
    if status == 1:
        tickets_processed += 1
        tickets_succeeded += 1
    else:
        tickets_processed += 1
        tickets_failed += 1

print("\nDone iterating through tickets.")
print("{} / {} tickets successfully handled.".format(tickets_succeeded,
                                                     tickets_processed))
if tickets_failed > 0:
    print("{} / {} tickets failed to be handled.".format(tickets_failed,
                                                         tickets_processed))
    # print("\nFailed tickets written to {}")
