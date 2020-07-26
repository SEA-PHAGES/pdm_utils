import re

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

class FeatureTableParser:
    def __init__(self, filehandle):
        self.filehandle = filehandle
        self._prev_line = ""

    def __next__(self):
        return self.next()
    
    def next(self):
        if self.filehandle.closed:
            raise IOError("Access to file has been closed unexpectedly.")

        
        if parse_tbl_data_type(self._prev_line) == "header":
            tbl_data = [self._prev_line]
        else:
            tbl_data = []
            for line in self.filehandle:
                self._prev_line = line
                if parse_tbl_data_type(line) != "header":
                    continue
                else:
                    tbl_data = [line]
                    break

            if not tbl_data:
                raise StopIteration("File contains no more Feature Tables")

        for line in self.filehandle:
            self._prev_line = line
            data_type = parse_tbl_data_type(line)
            if not(data_type == "header" or data_type is None):
                tbl_data.append(line)
            else:
                break

        record = convert_tbl_data_to_record(tbl_data)
        return record

    def __iter__(self):
        return self

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

