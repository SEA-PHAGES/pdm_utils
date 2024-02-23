import shlex
from subprocess import Popen, PIPE

PHAMERATOR_STORAGE = "http://databases.hatfull.org"


def get_db(database, config_file):
    command = f"pdm_utils get_db server -db {database} "

    if config_file:
        command += f"-c {config_file} "

    run_command(command)


def get_data(database, config_file, working_dir):
    command = (f"pdm_utils get_data {database} "
               f"-o {working_dir} ")

    if config_file:
        command += f"-c {config_file} "

    run_command(command)


def import_gb(database, config_file, working_dir, data_dir, import_table):
    command = (f"pdm_utils import {database} {data_dir} {import_table} "
               f"-o {working_dir} -p ")

    if config_file:
        command += f"-c {config_file} "

    run_command(command)


def update(database, config_file, update_file):
    command = f"pdm_utils update {database} -v "

    if config_file:
        command += f"-c {config_file} "
    if update_file:
        command += f"-f {update_file} "

    run_command(command)


def phammseqs(database, config_file, threads):
    command = (f"pdm_utils phamerate {database} -t {threads} ")

    if config_file:
        command += f"-c {config_file} "

    run_command(command)


def find_domains(database, config_file, output_dir, cpus):
    command = (f"pdm_utils find_domains {database} --cpus {cpus} "
               f"-o {output_dir} ")

    if config_file:
        command += f"-c {config_file} "

    run_command(command)

def find_transmembrane(database, config_file, output_dir, run_machine):
    command = (f"pdm_utils find_transmembrane {database} --run_machine {run_machine} -b 1000 -Mb 20 -v "
               f"-o {output_dir} ")

    if config_file:
        command += f"-c {config_file} "

    run_command(command)


def convert(database, config_file, version):
    command = (f"pdm_utils convert {database} -n {database}_v{version} "
               f"-s {version} ")

    if config_file:
        command += f"-c {config_file} "

    run_command(command)


def export(database, config_file, working_directory, cpus, name=None,
           phams=False):
    command = (f"pdm_utils export {database} sql -o {working_directory} "
               f"-m {database} ")

    if name:
        command += f"-n {name} "
    if phams:
        command += f"-pho -np {cpus} "

    if config_file:
        command += f"-c {config_file} "

    run_command(command)


def curl_metadata(database, outfile):
    command = (f"wget {PHAMERATOR_STORAGE}/{database}/{outfile.name} "
               f"-O {outfile}")

    run_command(command)


def push(database, key_file, config_file, export_directory):
    command = (f"pdm_utils push -d {export_directory} "
               f"-rd /databases/{database} "
               f"-c {config_file} "
               f"-k {key_file}")

    run_command(command)


def run_command(command):
    """Run the indicated command as a subprocess.

    If `verbose` is True, stdout and stderr will be returned and the
    code block that called this one can decide how to handle it.

    If the command invokes a program that is not globally executable,
    it will raise a FileNotFoundError with a message like "No such
    file or directory: 'xxx'".

    :param command: the command to run
    :type command: str
    :raises: FileNotFoundError
    :return: out, err
    """
    command = shlex.split(command)

    # Reading from stdout/stderr is blocking, no need to call sp.wait()
    with Popen(command, stdout=PIPE, stderr=PIPE, close_fds=True) as process:
        out = process.stdout.read().decode("utf-8")
        err = process.stderr.read().decode("utf-8")

    return out, err
