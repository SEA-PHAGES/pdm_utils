"""Functions for running DeepTMHMM and parsing the output."""

from operator import itemgetter
from itertools import groupby
import pathlib
import shlex, subprocess, os, sys
script_dir = os.getcwd()

def run_deeptmhmm(infile, batch_export_path, machine="local"):
    """Execute a DeepTMHMM job using the input file, and return
    the coordinates of predicted signal or transmembrane domains.

    :param infile: the path to a local FASTA protein sequence file
    :type infile: pathlib.Path
    :
    :return: pred_dict
    """

    print("run_deeptmhmm()")

    if machine == "local":
      predictor_path = os.getenv("DEEPTMHMM_PATH")
      print("predictor_path: ", predictor_path)
      print("batch_export_path: ", batch_export_path)

      if predictor_path and os.path.exists(predictor_path):
        subprocess.run([sys.executable, os.path.join(predictor_path, "predict.py"), "--fasta", infile, "--output-dir", os.path.join(batch_export_path)],
        cwd=predictor_path)
        with open(os.path.join(batch_export_path, "predicted_topologies.3line"), 'r') as f:
          data = f.read()
      else:
        print("Please set the DEEPTMHMM_PATH environment variable.")
        return 

    elif machine == "remote":
        job = _api.cli(args=f"--fasta {infile}")
        outfile = job.get_output_file("/predicted_topologies.3line")
        data = outfile.get_data()
    else:
        raise ValueError(
                    "DeepTMHMM run machine must either be remote or local.")


    return parse_3line(data)


def parse_3line(data):
    """Parse each translation's topology from a 3line FASTA-like file
    returned by DeepTMHMM.

    :param data: raw contents of the 3line file returned by DeepTMHMM
    :type data: str or bytes
    :return: domains
    """
    if isinstance(data, str):
        pass
    elif isinstance(data, bytes):
        data = data.decode("utf-8")     # decode bytes as utf-8 strings
    else:
        raise TypeError(f"unsupported type: {type(data)}; data should be a "
                        f"bytes or str object")

    data = data.split("\n>")                        # separate FASTA entries
    data = [x.rstrip().split("\n") for x in data]   # split to 3lines per entry

    # Add type hint for domains dictionary
    domains: dict[str, list[tuple[str, tuple[int, int]]]] = dict()
    for _, translation, topology in data:
        # Track indices of the signal/membrane-predicted states
        temp_dict = {"S": [], "M": []}
        for i, state in enumerate(topology):
            if state in ("S", "M"):
                temp_dict[state].append(i+1)

        # Interpret the tracked indices to continuous ranges
        domains[translation] = list()
        for r in _find_continuous_ranges(temp_dict["S"]):
            domains[translation].append(("signal", r))
        for r in _find_continuous_ranges(temp_dict["M"]):
            domains[translation].append(("transmembrane", r))

    return domains


def _find_continuous_ranges(data):
    """Identify and return continuous integer ranges in the input data.

    :param data: the input data to scan for continuous ranges
    :type data: list[int]
    """
    # Adapted with little modification from
    # https://stackoverflow.com/questions/2154249
    for k, g in groupby(enumerate(data), lambda x: x[0] - x[1]):
        group = (map(itemgetter(1), g))
        group = list(map(int, group))
        yield group[0], group[-1]
