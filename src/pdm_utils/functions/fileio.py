import csv
import os
import re
import textwrap
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation
from Bio.SeqRecord import SeqRecord

#GLOBAL VARIABLES
#-----------------------------------------------------------------------------
TBL_HEADER_FORMAT = re.compile(">Feature ..\|(\w+)\.(.)\|\n") 
TBL_FEATURE_FORMAT = re.compile("(\d+)\t(\d+)\t(\w+)\n")
TBL_REGIONS_FORMAT = re.compile("(\d+)\t(\d+)\n")
TBL_QUALIFIER_FORMAT = re.compile("\t\t\t(.+)\t(.+)\n")
TBL_OPEN_QUALIFIER_FORMAT = re.compile("\t\t\t(.+)\n")

#READING FUNCTIONS
#-----------------------------------------------------------------------------
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

def read_feature_table(filepath):
    """Reads a (five-column) feature table and parses the data into a seqrecord.

    :param filepath: Path to the five-column formatted feature table file.
    :type filepath: Path
    :returns: Returns a Biopython SeqRecord object with the table data.
    :rtype: list
    """  

    file_handle = filepath.open(mode="r") 

    while True:
        try:
            line = file_handle.readline()
        except:
            return None 
       
        if parse_tbl_data_type(line) == "header":
            break

    tbl_data = [line] 
    while True:
        try:
            line = file_handle.readline()
        except:
            break
        
        data_type = parse_tbl_data_type(line)
        if not (data_type == "header" or data_type is None):
            tbl_data.append(line)
        else:
            break

    record = convert_tbl_data_to_record(tbl_data)   
    return record

#WRITING FUNCTIONS
#-----------------------------------------------------------------------------

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

        file_handle.write(f">Feature gb|{record.id}."
                          f"{record.annotations['sequence_version']}|\n")

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
                if key in ["translation", "ribosomal_slippage"]:
                    if key == "ribosomal_slippage":
                        file_handle.write(f"\t\t\t{key}\n")
                    continue
                file_handle.write(f"\t\t\t{key}\t"
                                  f"{feature.qualifiers[key][0]}\n")
        
        file_handle.write("\n")
        file_handle.close()

#BASE FUNCTIONS
#-----------------------------------------------------------------------------
def convert_tbl_data_to_record(tbl_data):
    """Converts string lines from a five_column feature table to a seqrecord.
    """
    header = tbl_data[0]
    if parse_tbl_data_type(header) != "header":
        raise ValueError("Five column table data does not have a proper header")

    header_split = re.split(TBL_HEADER_FORMAT, header) 
    record = SeqRecord(Seq(''))
    record.id = header_split[1] 
    record.annotations = {"sequence_version" : str(header_split[2])}

    if len(tbl_data) == 1:
        return record

    features = []
    qualifiers = {}
    coordinates = []
    last_data_type = None
    feature_split = []
    for line in tbl_data[1:]:
        data_type = parse_tbl_data_type(line)

        if data_type == "header":
            raise SyntaxError("Found table header within data.")
        elif data_type == "feature":
            if coordinates:
                feature = feature_data_to_seqfeature(coordinates,
                                                        feature_split[3])
                feature.qualifiers = qualifiers
                features.append(feature)
            coordinates = []
            qualifiers = {}

            feature_split = re.split(TBL_FEATURE_FORMAT, line)
            coordinates.append((int(feature_split[1]), int(feature_split[2])))
        elif data_type == "region":
            if not last_data_type in ["feature", "region"]:
                raise SyntaxError(
                        "Five column feature table not formatted correctly.\n"
                        "Orphaned multiple region data found.")

            regions_split = re.split(TBL_REGIONS_FORMAT, line)
            coordinates.append((int(regions_split[1]), int(regions_split[2])))
        elif data_type == "qualifier":
            qualifier_split = re.split(TBL_QUALIFIER_FORMAT, line)
            qualifiers[qualifier_split[1]] = [qualifier_split[2]]
        elif data_type == "open_qualifier":
            open_qualifier_split = re.split(TBL_OPEN_QUALIFIER_FORMAT, line)
            qualifiers[open_qualifier_split[1]] = [""]
        else:
            raise SyntaxError(
                        "Five column feature table not formatted correctly.\n"
                        "Unrecognized data format.")

        last_data_type = data_type

    if coordinates:
        feature = feature_data_to_seqfeature(coordinates,
                                                feature_split[3])
        feature.qualifiers = qualifiers
        features.append(feature)

    record.features = features
    return record

def feature_data_to_seqfeature(coordinates, feature_type):
    locations = []
    for coordinate in coordinates:
        if coordinate[0] <= coordinate[1]:
            start = int(coordinate[0]) - 1
            stop = int(coordinate[1])
            strand = 1
        else:
            start = int(coordinate[1]) - 1
            stop = int(coordinate[0])
            strand = -1
        
        locations.append(FeatureLocation(start, stop, strand=strand))

    if len(locations) == 1:
        feature = SeqFeature(locations[0], type=feature_type)
    else:
        feature = SeqFeature(CompoundLocation(locations), type=feature_type)

    return feature



def parse_tbl_data_type(data_line):
    """Parses a five-column table data line and returns it's structure type.

    :param data_line: Line of data from a five-column feature table.
    :type data_line: str
    :returns: Returns the type of data line the line is from the table.
    :rtype: str
    """
    data_type = None
    if not re.match(TBL_HEADER_FORMAT, data_line) is None:
        data_type = "header"
    elif not re.match(TBL_FEATURE_FORMAT, data_line) is None:
        data_type = "feature"
    elif not re.match(TBL_REGIONS_FORMAT, data_line) is None:
        data_type = "region"
    elif not re.match(TBL_QUALIFIER_FORMAT, data_line) is None:
        data_type = "qualifier"
    elif not re.match(TBL_OPEN_QUALIFIER_FORMAT, data_line) is None:
        data_type = "open_qualifier"
   
    return data_type 
