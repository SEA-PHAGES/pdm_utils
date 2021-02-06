import csv
import os
import textwrap
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

from pdm_utils.classes.fileio import FeatureTableParser
from pdm_utils.functions import multithread

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
TBL_EXCLUDED_QUALIFIERS = ["translation"]
TBL_SPECIAL_QUALIFIERS = ["ribosomal_slippage"]


# READING FUNCTIONS
# -----------------------------------------------------------------------------
def retrieve_data_dict(filepath):
    """Open file and retrieve a dictionary of data.

    :param filepath:
        Path to file containing data and column names.
    :type filepath: Path
    :returns:
        A list of elements, where each element is a dictionary
        representing one row of data.
        Each key is a column name and each value is the data
        stored in that field.
    :rtype: list
    """
    data_dicts = []
    with filepath.open(mode='r') as file:
        file_reader = csv.DictReader(file)
        for dict in file_reader:
            data_dicts.append(dict)
    return data_dicts


def read_feature_table(filehandle):
    """Reads a (five-column) feature table and parses the data into a seqrecord.

    :param filepath: Path to the five-column formatted feature table file.
    :type filepath: Path
    :returns: Returns a Biopython SeqRecord object with the table data.
    :rtype: SeqRecord
    """
    parser = FeatureTableParser(filehandle)

    try:
        record = parser.next()
    except:
        record = None

    return record


def parse_feature_table(filehandle):
    """Takes a (five-column) feature table(s) file handle and parses the data.

    :param filehandle: Handle for a five-column formatted feature table file:
    :returns: Returns a feature table file parser generator.
    :rtype: FeatureTableFileParser
    """
    return FeatureTableParser(filehandle)


def reintroduce_fasta_duplicates(ts_to_gs, filepath):
    """Reads a fasta file and reintroduces (rewrittes) duplicate sequences
    guided by an ungapped translation to sequence-id map

    :param filepath: Path to fasta-formatted multiple sequence file
    :type filepath: pathlib.Path
    :param ts_to_gs: Dictionary mapping unique translations to sequence-ids
    :type ts_to_gs: dict
    """
    with filepath.open(mode="r") as filehandle:
        records = []
        for record in SeqIO.parse(filehandle, "fasta"):
            records.append(record)

    gs_to_ts = {}
    for record in records:
        translation = record.seq
        ungapped_translation = translation.ungap(gap="-")

        geneids = ts_to_gs.get(str(ungapped_translation))
        if geneids is not None:
            for geneid in geneids:
                gs_to_ts[geneid] = str(translation)
        else:
            gs_to_ts[record.id] = str(translation)

    write_fasta(gs_to_ts, filepath)


# WRITING FUNCTIONS
# -----------------------------------------------------------------------------
def export_data_dict(data_dicts, file_path, headers, include_headers=False):
    """Save a dictionary of data to file using specified column headers.

    Ensures the output file contains a specified number of columns,
    and it ensures the column headers are exported as well.

    :param data_dicts:
        list of elements, where each element is a dictionary.
    :type data_dicts: list
    :param file_path: Path to file to export data.
    :type file_path: Path
    :param headers:
        List of strings to define the column order in the file.
        If include_headers is selected, the first row of the file
        will contain each string.
    :type headers: list
    :param include_headers:
        Indicates whether the file should contain a
        row of column names derived from the headers parameter.
    :type include_headers: bool
    """

    headers_dict = {}
    for header in headers:
        headers_dict[header] = header
    # with open(file_path, "w") as file_handle:
    with file_path.open("w") as file_handle:
        file_writer = csv.DictWriter(file_handle, headers)
        if include_headers:
            file_writer.writerow(headers_dict)
        for data_dict in data_dicts:
            file_writer.writerow(data_dict)


def write_fasta(ids_seqs, infile_path, name=None):
    """
    Writes the input genes to the indicated file in FASTA multiple
    sequence format (unaligned).
    :param id_seqs: the ids and sequences to be written to file
    :type genes: dict
    :param infile_path: the path of the file to write the genes to
    :type infile: Path
    :type infile: str
    """
    if isinstance(infile_path, str):
        file_handle = open(infile_path, "w")
    elif isinstance(infile_path, Path):
        file_handle = infile_path.open(mode="w")
    else:
        raise TypeError("File path type not supported.")

    if name is not None:
        file_handle.write(f"#{name}\n")

    for id, seq in ids_seqs.items():
        split_seq = textwrap.wrap(seq, 60)
        wrapped_seq = "\n".join(split_seq)
        file_handle.write(f">{id}\n{wrapped_seq}\n\n")

    file_handle.close()


def write_database(alchemist, version, export_path, db_name=None):
    """Output .sql file from the selected database.

    :param alchemist: A connected and fully built AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param version: Database version information.
    :type version: int
    :param export_path: Path to a valid dir for file creation.
    :type export_path: Path
    """
    if db_name is None:
        db_name = alchemist.database

    sql_path = export_path.joinpath(f"{db_name}.sql")
    os.system(f"mysqldump -u {alchemist.username} -p{alchemist.password} "
              f"--skip-comments {alchemist.database} > {str(sql_path)}")
    version_path = sql_path.with_name(f"{db_name}.version")
    version_path.touch()
    version_path.write_text(f"{version}")


def write_seqrecord(seqrecord, file_path, file_format):
    file_handle = file_path.open(mode="w")

    if isinstance(seqrecord, list):
        for record in seqrecord:
            SeqIO.write(record, file_handle, file_format)
            file_handle.write("\n")
    else:
        SeqIO.write(seqrecord, file_handle, file_format)

    file_handle.close()


def write_feature_table(seqrecord_list, export_path, verbose=False):
    """Outputs files as five_column tab-delimited text files.

    :param seq_record_list: List of populated SeqRecords.
    :type seq_record_list: list[SeqRecord]
    :param export_path: Path to a dir for file creation.
    :type export_path: Path
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    """
    if verbose:
        print("Writing selected data to files...")
    for record in seqrecord_list:
        if verbose:
            print(f"...Writing {record.name}...")
        file_name = f"{record.name}.tbl"
        file_path = export_path.joinpath(file_name)
        file_handle = file_path.open(mode='w')

        accession = record.id
        version = record.annotations.get("sequence_version")
        if version:
            accession = ".".join([accession, version])

        prefix = record.annotations.get("tbl_prefix")
        if not prefix:
            prefix = ""

        file_handle.write(f">Feature {prefix}|{accession}|\n")

        for feature in record.features:
            location = feature.location
            if isinstance(feature.location, CompoundLocation):
                location = feature.location.parts[0]

            if feature.strand == 1:
                start = location.start + 1
                stop = location.end
            elif feature.strand == -1:
                start = location.end
                stop = location.start + 1

            file_handle.write(f"{start}\t{stop}\t{feature.type}\n")

            if isinstance(feature.location, CompoundLocation):
                if len(feature.location.parts) > 1:
                    for location in feature.location.parts[1:]:
                        if feature.strand == 1:
                            start = location.start + 1
                            stop = location.end
                        elif feature.strand == -1:
                            start = location.end
                            stop = location.start + 1

                        file_handle.write(f"{start}\t{stop}\n")

            for key in feature.qualifiers.keys():
                if key in TBL_EXCLUDED_QUALIFIERS:
                    continue
                elif key in TBL_SPECIAL_QUALIFIERS:
                    file_handle.write(f"\t\t\t{key}\n")
                    continue

                qualifier_values = feature.qualifiers[key]
                for value in qualifier_values:
                    file_handle.write(f"\t\t\t{key}\t{value}\n")

        file_handle.write("\n")
        file_handle.close()


# MULTITHREAD WRAPPERS
# -----------------------------------------------------------------------------

def write_seqrecords(seqrecord_list, file_format, export_path,
                     export_name=None, concatenate=False, threads=1,
                     verbose=False):
    """Outputs files with a particuar format from a SeqRecord list.

    :param seq_record_list: List of populated SeqRecords.
    :type seq_record_list: list[SeqRecord]
    :param file_format: Biopython supported file type.
    :type file_format: str
    :param export_path: Path to a dir for file creation.
    :type export_path: Path
    :param concatenate: A boolean to toggle concatenation of SeqRecords.
    :type concaternate: bool
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    """
    record_dictionary = {}
    if concatenate:
        if export_name is None:
            record_dictionary.update({export_path.name: seqrecord_list})
        else:
            record_dictionary.update({export_name: seqrecord_list})
    else:
        for record in seqrecord_list:
            record_dictionary.update({record.name: record})

    if verbose:
        print(f"Writing selected data to files at '{export_path}'...")

    work_items = []
    for record_name, records in record_dictionary.items():
        file_name = f"{record_name}.{file_format}"
        file_path = export_path.joinpath(file_name)

        work_items.append((records, file_path, file_format))
    multithread.multithread(work_items, threads, write_seqrecord,
                            verbose=verbose)
