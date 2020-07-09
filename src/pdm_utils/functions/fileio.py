import csv
import os
import textwrap
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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

def write_fasta(ids_seqs, infile_path):
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

    for id, seq in ids_seqs.items():
        split_seq = textwrap.wrap(seq, 60)
        wrapped_seq = "\n".join(split_seq)
        file_handle.write(f">{id}\n{wrapped_seq}\n\n")

    file_handle.close()

def write_database(alchemist, version, export_path):
    """Output .sql file from the selected database.

    :param alchemist: A connected and fully built AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param version: Database version information.
    :type version: int
    :param export_path: Path to a valid dir for file creation.
    :type export_path: Path
    """
    sql_path = export_path.joinpath(f"{alchemist.database}.sql")
    os.system(f"mysqldump -u {alchemist.username} -p{alchemist.password} "
              f"--skip-comments {alchemist.database} > {str(sql_path)}")
    version_path = sql_path.with_name(f"{alchemist.database}.version")
    version_path.touch()
    version_path.write_text(f"{version}")

def write_seqrecord(seqrecord_list, file_format, export_path, concatenate=False,
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
    if verbose:
        print("Writing selected data to files...")

    record_dictionary = {}
    if concatenate:
        record_dictionary.update({export_path.name:seqrecord_list})
    else:
        for record in seqrecord_list:
            record_dictionary.update({record.name:record})

    for record_name in record_dictionary.keys():
        if verbose:
            print(f"...Writing {record_name}...")
        file_name = f"{record_name}.{file_format}"
        if concatenate:
            file_path = export_path.parent.joinpath(file_name)
            export_path.rmdir()
        else:
            file_path = export_path.joinpath(file_name)

        file_handle = file_path.open(mode='w')
        records = record_dictionary[record_name]
        if isinstance(records, list):
            for record in records:
                SeqIO.write(record, file_handle, file_format)
                file_handle.write("\n")
        else:
            SeqIO.write(record_dictionary[record_name], file_handle, file_format)

        file_handle.close()

def write_five_column_table(seqrecord_list, export_path, verbose=False):
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

        file_handle.write(f">Feature {record.id}\n")

        for feature in record.features[1:]:
            if feature.strand == 1:
                start = feature.location.start
                stop = feature.location.end
            elif feature.strand == -1:
                start = feature.location.end
                stop = feature.location.start

            file_handle.write(f"{start+1}\t{stop}\t{feature.type}\n")
            for key in feature.qualifiers.keys():
                if key in ["translation"]:
                    continue
                file_handle.write(f"\t\t\t{key}\t"
                                  f"{feature.qualifiers[key][0]}\n")

        file_handle.close()

