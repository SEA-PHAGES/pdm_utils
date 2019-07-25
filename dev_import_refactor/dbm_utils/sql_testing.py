import pymysql
import getpass, subprocess, os
# import paramiko
# paramiko.util.log_to_file("/tmp/paramiko.log")


def get_details():

    user = getpass.getpass(prompt="mysql username: ")
    pwd = getpass.getpass(prompt="mysql password: ")

    return (user, pwd)

def create_database(user, pwd):

    connection = pymysql.connect(host = "localhost",
                                    user = user,
                                    password = pwd,
                                    # db = "test_schema4",
                                    # charset = "utf8mb4",
                                    cursorclass = pymysql.cursors.DictCursor)

    try:
        with connection.cursor() as cursor:
            sql = "CREATE DATABASE test_schema4"
            cursor.execute(sql)
        connection.commit()
    except:
        print("Unable to create database")






def drop_database(user, pwd):

    connection = pymysql.connect(host = "localhost",
                                    user = user,
                                    password = pwd,
                                    # db = "test_schema4",
                                    # charset = "utf8mb4",
                                    cursorclass = pymysql.cursors.DictCursor)

    try:
        with connection.cursor() as cursor:
            sql = "DROP DATABASE test_schema4"
            cursor.execute(sql)
        connection.commit()
    except:
        print("Unable to create database")











def query_database(user, pwd, db = "Actino_Draft"):

    connection = pymysql.connect(host = "localhost",
                                    user = user,
                                    password = pwd,
                                    db = db,
                                    cursorclass = pymysql.cursors.DictCursor)

    with connection.cursor() as cursor:
        sql = "SELECT PhageID FROM phage"
        cursor.execute(sql)

        # Returns a list of all items queried.
        result = cursor.fetchall()


    return result














if __name__ == '__main__':
    print("Get db creds")
    user, pwd = get_details()
    print("Creating temp database")
    create_database(user, pwd)
    print("Dropping temp database")
    drop_database(user, pwd)
    print("Done")
