"""
Script that automates construction and sending of emails to notify
expedited submitters of the import status of their genomes.
"""
import argparse
import csv
import email
import imaplib
import os
import shlex
import shutil
import smtplib
import subprocess as sp
import sys
import pathlib
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

import urllib3

# Email login credentials and message templates
USERNAME = "phamerator.qc@gmail.com"
PASSWORD = "egngzitquqpchdux"
TEMPLATES = {
    "success": "{} has been successfully imported into the Phamerator "
               "database. No further action is needed.\n\nBest,\n\n"
               "A friendly Hatfull lab server",
    "submit": "The following genomes are ready for submission to Genbank:\n"
              "{}\n\nBest,\n\n"
              "A friendly Hatfull lab server",
    "fail": "{} could not be imported into the Phamerator database. See the "
            "attached log file(s) for details.\n\nA Google spreadsheet ("
            "https://docs.google.com/spreadsheets/d/1LU-ueqWIwjiWWjGFVuLlpt"
            "r6ha4zF4re2ScmYhdExPc/edit?usp=sharing) has been created to "
            "help interpret the most common errors. Feel free to reach out "
            "if any remain indecipherable.\n\nBest,\n\n"
            "A friendly Hatfull lab server"}

# URLs
ASN1_LINK = "https://phagesdb.org/data/submissionfiles/"

# Local paths
OUT_DIR = os.getcwd()
ZIP_FILE = os.path.join(OUT_DIR, "submission_files.zip")
ASN1_DIR = os.path.join(OUT_DIR, "submission_files")

# Build HTTPS pool so we can use urllib3 to download content
http_pool = urllib3.PoolManager()

MAINTAINER_EMAILS = {"djs@pitt.edu", "laa89@pitt.edu",
                     "phamerator.qc@gmail.com"}


# Low-level functions
def download_to_file(url, out):
    """
    Downloads the content from the specified url and writes it to the
    specified file.
    :param url: url where content should be downloaded from
    :param out: filepath where content should be written to
    :return:

    This example would save the raw html from the Python organization
    homepage into a file called 'python.html':
    >>> url = "https://www.python.org"
    >>> out = "python.html"
    >>> download_to_file(url, out)
    """
    with open(out, "wb") as fh:
        request = http_pool.request("GET", url, preload_content=False)
        if request.status == 200:
            for chunk in request.stream(512000):    # stream chunks of 500kb
                fh.write(chunk)
        else:
            print("Cannot reach url")
    return


def unzip_file(zip_file):
    """
    Unzips the specified zip file
    :param zip_file: filepath where zip file resides
    :return:
    """
    command = f"unzip -qq {zip_file}"
    sp.check_call(args=shlex.split(command))


def get_submission_files():
    """
    Downloads all the submission files in the queue.
    Unzips the downloaded zip file, then returns a list of the
    submission files' full paths
    :return: asn1s
    """
    download_to_file(ASN1_LINK, ZIP_FILE)
    unzip_file(ZIP_FILE)

    # List the files in ASN1_DIR
    return [os.path.join(ASN1_DIR, x) for x in os.listdir(ASN1_DIR)]


def get_phage_ids_from_table(table_file):
    """
    Parses phage_ids out of import table and returns them.
    :param table_file: the path to the table of import actions
    :return: phage_ids
    """
    phage_ids = list()

    with open(table_file, "r") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            phage_ids.append(row['phage_id'])

    return phage_ids


def get_reply_to(username, password, phage_id):
    """
    Logs into phamerator.qc@gmail.com using the given password and 
searches
    the inbox for an email from no-reply@phagesdb.org with the indicated
    phageid in the subject line.
    :param username: the ticket email account's address
    :param password: the password to log into the ticket email account
    :param phage_id: the phage_id to find an email for
    :return: reply_to
    """
    server = imaplib.IMAP4_SSL("imap.gmail.com")

    # Log into ticket email account
    status, message = server.login(username, password)
    if status != "OK":
        print(f"Could not log into {username} with the given password")
        return None

    # Select the inbox
    status, message = server.select("INBOX")
    if status != "OK":
        print("Could not connect to inbox")
        return None

    # Search the inbox for emails from no-reply@phagesdb.org AND with
    # phage_id in the subject line
    search = f"FROM 'no-reply@phagesdb.org' SUBJECT '{phage_id}'"
    status, message = server.search(None, search)
    if status != "OK":
        print("Could not search the inbox")
        return None

    # Get the index of the first hit
    indexes = message[0].decode("utf-8").split()
    if len(indexes) == 0:
        print("No emails match search criteria")
        return None
    fetch = indexes[0]

    # Fetch the email
    status, message = server.fetch(fetch, "(RFC822)")
    if status != "OK":
        print("Could not fetch the email")
        return None

    # Convert the message to an email.message.Message object
    message = email.message_from_bytes(message[0][1])

    # Return the "Reply-to" field
    return message["Reply-to"]


def add_attachment(message, filepath):
    """
    Adds the given filepath's contents as an attachment to the
    given message.
    :param message: the email to attach a file to
    :param filepath: the file to attach to the email
    :return: the email with the file attached to it
    """
    part = MIMEBase('application', "octet-stream")
    with open(filepath, 'rb') as file:
        part.set_payload(file.read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition',
                    'attachment; filename="{}"'.format(
                        os.path.basename(filepath)))
    message.attach(part)

    return message


def send_email(username, password, to_addr, message):
    """
    Connects to gmail server, sends email, and then closes connection.
    :param username: the ticket email account's address
    :param password: the password to log into the ticket email account
    :param to_addr: primary email recipient
    :param message: email.MIMEMultipart.MIMEMultipart() message
    """
    server = smtplib.SMTP("smtp.gmail.com", 587)
    server.starttls()
    server.login(username, password)
    server.sendmail(username, to_addr, message.as_string())
    server.quit()


# High-level functions
def send_successful_emails(username, password, phage_ids):
    """
    Makes calls to low-level functions to build and send emails
    notifying expedited submitters that their genomes were successfully
    imported into the database.
    :param username: the ticket email account's address
    :param password: the password to log into the ticket email account
    :param phage_ids: phage_id to submission filepath mappings
    successfully imported genomes
    """
    print("Sending successful import emails...")
    print("===================================")

    status = "success"

    # Send submitters their individual emails
    for phage_id in phage_ids.keys():
        # Get submitter's email address
        reply_to = get_reply_to(username, password, phage_id)
        # If email address retrieved, send automated email
        if reply_to is not None and reply_to != "":
            to_addr = list({reply_to} | MAINTAINER_EMAILS)
            send_to_submitter(username, password, to_addr, phage_id, 
status)
            print(f"Automated email sent to {reply_to} for {phage_id}")
        # Else notify that manual email is needed
        else:
            print(f"Manual email required for {phage_id}")

    # Send "ready for submission" email
    to_addr = list(MAINTAINER_EMAILS)
    submit_to_genbank(username, password, to_addr, phage_ids)


def send_failed_emails(username, password, phage_ids):
    """
    Makes calls to low-level functions to build and send emails
    notifying expedited submitters that their genomes failed to be
    imported into the database.
    :param username: the ticket email account's address
    :param password: the password to log into the ticket email account
    :param phage_ids: the phage_ids not successfully imported
    """
    print("Sending unsuccessful import emails...")
    print("=====================================")

    status = "fail"

    # Iterate over unsuccessful phages
    for phage_id in phage_ids.keys():
        # Get submitter's email address
        reply_to = get_reply_to(username, password, phage_id)
        # If email address retrieved, send automated email
        if reply_to is not None and reply_to != "":
            to_addr = list({reply_to} | MAINTAINER_EMAILS)
            attach = phage_ids[phage_id]
            send_to_submitter(username, password, to_addr,
                              phage_id, status, attach)
            print(f"Automated email sent to {reply_to} for {phage_id}")
        # Else notify that manual email is needed
        else:
            print(f"Manual email required for {phage_id}")


def send_to_submitter(username, password, to_addr, phage_id, status,
                      attach=None):
    """
    Builds an email notifying an expedited submitter the import status
    of their genome. Then logs into phamerator.qc@gmail.com using the
    given credentials and sends the email.
    :param username: the ticket email account's address
    :param password: the password to log into the ticket email account
    :param to_addr: the email address(es) to send the email to
    :param phage_id: the phage whose status is being reported
    :param status: the status to report
    :param attach: a list of file attachment(s)
    """
    if attach is None:
        attach = list()

    message = MIMEMultipart()

    # Set up header
    message["From"] = username
    message["Reply-to"] = username
    message["To"] = ','.join(to_addr)

    # If status is success, build successful email
    if status == "success":
        subject = f"{phage_id} successfully imported into Phamerator"
    elif status == "fail":
        subject = f"{phage_id} failed to be imported into Phamerator"
    else:
        print("Unknown status encountered - {}".format(status))
        return

    message["Subject"] = subject

    body = TEMPLATES[status].format(phage_id)
    message.attach(MIMEText(body, 'plain'))

    for filepath in attach:
        add_attachment(message, filepath)

    # Send the email
    send_email(username, password, to_addr, message)


def submit_to_genbank(username, password, to_addr, phage_ids):
    """

    :param username: the ticket email account's address
    :param password: the password to log into the ticket email account
    :param to_addr: the email address(es) to send the email to
    :param phage_ids: the phage_ids and submission files dictionary
    """
    print("Sending submission email...")
    print("===========================")

    status = "submit"

    # Don't bother sending an email if there are no phages ready for submission
    if len(phage_ids) == 0:
        print("No phages ready for submission - aborting")
        return

    phage_string = "\n"
    for phage_id, submission_file in phage_ids.items():
        if submission_file is not None:
            phage_string += f"{phage_id} - (attached)\n"
        else:
            phage_string += f"{phage_id} - (not attached)\n"

    # Build the email
    message = MIMEMultipart()

    # Set up header
    message["From"] = username
    message["Reply-to"] = username
    message["To"] = ','.join(to_addr)

    message["Subject"] = "Phages ready for submission to Genbank"
    body = TEMPLATES[status].format(phage_string)

    # Attach the body
    message.attach(MIMEText(body, 'plain'))

    # Add attachment(s)
    for attach in phage_ids.values():
        for filepath in attach:
            add_attachment(message, filepath)

    # Send the email
    send_email(username, password, to_addr, message)


def parse_args(unparsed_args):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("import_path", type=pathlib.Path,
                        help="path to the import pipeline results")

    args = parser.parse_args(unparsed_args)
    return args 


def main(unparsed_args):
    args = parse_args(unparsed_args)
   
    email_submitters(args.import_path)


def email_submitters(import_path): # success_import_path, failed_import_path):
    success_import_path = import_path.joinpath("success")
    if success_import_path.is_dir():
        # Deal with successful import actions
        import_table = os.path.join(success_import_path, "import_table.csv")
        phage_ids = get_phage_ids_from_table(import_table)
        submission_files = get_submission_files()

        phage_dict = dict()
        for phage_id in phage_ids:
            for submission_file in submission_files:
                if phage_id.lower() in submission_file.lower():
                    attach = phage_dict.get(phage_id, [])
                    attach.append(submission_file)
                    phage_dict[phage_id] = attach
                    break       # Only attach one submission file per genome
            if phage_dict.get(phage_id) is None:
                phage_dict[phage_id] = None

        send_successful_emails(USERNAME, PASSWORD, phage_dict)

    failed_import_path = import_path.joinpath("fail")
    if failed_import_path.is_dir():
        # Deal with failed import actions
        import_table = os.path.join(failed_import_path, "import_table.csv")
        phage_ids = get_phage_ids_from_table(import_table)
        log_dir = os.path.join(failed_import_path, "logs")
        logs = [os.path.join(log_dir, x) for x in os.listdir(log_dir)]

        phage_dict = dict()
        for phage_id in phage_ids:
            for log in logs:
                if phage_id.lower() in log.lower():
                    attach = phage_dict.get(phage_id, [])
                    attach.append(log)
                    phage_dict[phage_id] = attach
            if phage_dict.get(phage_id) is None:
                phage_dict[phage_id] = None

        send_failed_emails(USERNAME, PASSWORD, phage_dict)

    # Cleanup .zip and ASN1 files
    os.remove(ZIP_FILE)
    shutil.rmtree(ASN1_DIR)


if __name__ == "__main__":
    main(sys.argv[1:])

