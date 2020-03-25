"""Represents a collection of data about a tRNA feature that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

import re
from collections import OrderedDict

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Alphabet import IUPAC

from pdm_utils.functions import basic


GENBANK_AMINO_ACIDS = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly",
                       "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser",
                       "Thr", "Trp", "Tyr", "Val", "OTHER"}
SPECIAL_AMINO_ACIDS = {"fMet": "Met", "Ile2": "Ile", "Pyl": "OTHER",
                       "SeC": "OTHER", "Stop": "OTHER", "Sup": "OTHER",
                       "Undet": "OTHER"}


class TrnaFeature:
    def __init__(self):
        """
        Constructor method for a tRNA object.
        """
        # The only attribute a tRNA really needs to have is a sequence
        self.seq = Seq("", IUPAC.ambiguous_dna)
        self.length = 0

        # Information about the tRNA with respect to its parent genome
        self.genome_id = ""     # Identifier for the parent genome
        self.genome_length = -1
        self.start = -1         # Start coord in parent genome, 0-indexed
        self.stop = -1          # Stop coord in parent genome, 0-indexed
        self.coordinate_format = ""
        self.orientation = ""   # "F" or "R", etc. - relative to parent genome

        # MySQL tRNAs will also need these:
        self.id = ""            # Identifier for this gene
        self.name = ""          # tRNA number in the parent genome
        self.amino_acid = ""    # Which aa. is this tRNA likely charged with?
        self.anticodon = ""     # Which anticodon is at bottom of the C-loop?
        self.structure = ""     # What is the predicted secondary structure?

        # Aragorn data
        self.aragorn_data = None

        # tRNAscan-SE data
        self.trnascan_data = None

        # Genbank-formatted flat files should also have these attributes:
        self.seqfeature = None  # BioPython SeqFeature object for this tRNA
        self.locus_tag = ""     # Gene ID comprised of PhageID and Gene name
        self._locus_tag_num = ""
        self.gene = ""          # Gene number parsed from feature
        self.product = ""       # Raw product field
        self.note = ""          # Raw note field

        # Useful for processing data from various sources:
        self.evaluations = list()
        self.type = "tRNA"
        self._start_stop_orient_id = tuple()
        self._end_orient_id = tuple()
        self._start_end_id = tuple()

    # TODO: create base feature class - use this version of the method unless
    #  Travis objects, because it's simpler and better documented.
    def set_locus_tag(self, tag="", delimiter="_", check_value=None):
        """
        Populate the locus tag attribute.
        :param tag: Input locus_tag data
        :type tag: str
        :param delimiter: Value used to split locus_tag data
        :type delimiter: str
        :param check_value: Genome name or other value that will be
        used to parse the locus_tag to identify the feature number
        :type check_value: str
        """
        # Set locus tag
        self.locus_tag = tag

        # If no check_value was given, use self.genome_id
        if check_value is None:
            check_value = self.genome_id

        # Compile a search pattern and split the incoming tag on delimiter
        pattern = re.compile(check_value.lower())
        parts = tag.split(delimiter)

        # Attempt to parse locus tag number from the locus tag
        parsed_value = None
        for i in range(len(parts)):
            # Search each part of the locus tag for check_value
            search_result = pattern.search(parts[i].lower())
            if search_result is not None:
                if i == len(parts) - 1:
                    # If the check_value was found in last part of locus tag,
                    # value is probably everything after check_value
                    parsed_value = parts[i][search_result.end():]
                else:
                    # If the check_value was found anywhere else in locus tag,
                    # value is probably found in the next index
                    parsed_value = parts[i + 1]
                # No need to keep searching if we found a match
                break

        if parsed_value is None:
            # Last-ditch effort to find locus tag number
            parsed_value = parts[-1]

        # If the value we parsed is a digit, populate _locus_tag_num
        if parsed_value.isdigit():
            self._locus_tag_num = parsed_value

    # TODO: create base feature class - fully equivalent to version in Cds,
    #  but documented differently.
    def set_name(self, value=None):
        """
        Set the feature name.
        :param value: Value to use for `name` regardless of `gene`
        and `_locus_tag_num`
        :type value: str
        :return:
        """
        # If a value was given, we'll use this by default
        if value is not None:
            self.name = value
        # If no value was given we'll use the following hierarchy:
        # 1. self.gene
        # 2. self._locus_tag_num
        # 3. leave as default ""
        elif self.gene != "":
            self.name = self.gene
        elif self._locus_tag_num != "":
            self.name = self._locus_tag_num
        else:
            pass

    # TODO: create base feature class - use this version of the method unless
    #  Travis objects, because it's documented more cleanly and avoids using
    #  keyword 'format'.
    def set_orientation(self, value, fmt, capitalize=False):
        """
        Set the orientation based on the indicated format.
        :param value: orientation value
        :type value: int or str
        :param fmt: how orientation should be formatted
        :type fmt: str
        :param capitalize: whether to capitalize the first letter of
        orientation
        :type capitalize: bool
        :return:
        """
        self.orientation = basic.reformat_strand(value, fmt, capitalize)

    # TODO: create base feature class - use this version of the method unless
    #  Travis objects, because it's more Pythonic.
    def get_begin_end(self):
        """
        Accesses feature coordinates in transcription begin-end format.
        :return: (begin, end)
        """
        # Get a copy of the orientation in fr_long format:
        orientation = basic.reformat_strand(self.orientation, "fr_long")

        if orientation == "forward":
            # Rightward transcribed gene
            begin, end = self.start, self.stop
        elif orientation == "reverse":
            # Leftward transcribed gene
            begin, end = self.stop, self.start
        else:
            # Unexpected orientation
            begin = end = -1

        return begin, end                       # tuple format is implicit

    # TODO: create base feature class
    def set_location_id(self):
        """
        Create identifier tuples containing feature location data. For
        this method we only care about gene boundaries and will ignore
        any multi-part elements to the gene.
        :return:
        """
        # Tuple with ALL relevant information
        self._start_stop_orient_id = (self.start, self.stop, self.orientation)

        begin, end = self.get_begin_end()

        # If only end coordinate is used, we need to know what orientation the
        # gene is in - two genes could theoretically be transcribed in opposed
        # directions and have the same stop position
        self._end_orient_id = (end, self.orientation)

        # Or we can have transcription begin/end coordinates, in which case
        # the direction can be gleaned from these values
        self._start_end_id = (begin, end)

    # TODO: create base feature class - use this version of the method unless
    #  Travis objects, because it has slightly improved logic.
    def reformat_start_and_stop(self, fmt):
        """
        Convert existing start and stop coordinates to the indicated
        new format; also updates the coordinate format attribute to
        reflect any change.
        :param fmt: the new desired coordinate format
        :type fmt: str
        :return:
        """
        # Only need to do work if the new coordinate format differs from
        # the existing coordinate format
        if fmt != self.coordinate_format:
            self.start, self.stop = basic.reformat_coordinates(
                self.start, self.stop, self.coordinate_format, fmt)
            self.coordinate_format = fmt

    def set_amino_acid(self, value):
        """
        Set the `amino_acid` attribute.
        :param value: what to use as the amino acid
        :type value: str
        :raise: ValueError
        :return:
        """
        if isinstance(value, str):
            self.amino_acid = value
        else:
            raise ValueError(f"Invalid type '{type(value)}' for amino acid.")

    def set_anticodon(self, value):
        """
        Set the `anticodon` attribute.
        :param value: what to use as the anticodon
        :type value: str
        :raise: ValueError
        :return:
        """
        if isinstance(value, str):
            self.anticodon = value
        else:
            raise ValueError(f"Invalid type '{type(value)}' for anticodon.")

    def set_structure(self, value):
        """
        Set the `structure` attribute.
        :param value: what to use as the structure
        :type value: str
        :raise: ValueError
        :return:
        """
        if isinstance(value, str):
            self.structure = value
        else:
            raise ValueError(f"Invalid type '{type(value)}' for structure.")

    # TODO: create base feature class - fully equivalent to the version in Cds,
    #  but better documented.
    def set_seqfeature(self):
        """
        Create a SeqFeature object with which to populate the
        `seqfeature` attribute.
        :return:
        """
        # SeqFeature coordinates are 0-based half-open
        start, stop = basic.reformat_coordinates(
            self.start, self.stop, self.coordinate_format, "0_half_open")

        # SeqFeature orientation is (-1, 1) instead of ("R", "F")
        strand = basic.reformat_strand(self.orientation, "numeric")

        # Standard genes will have start < stop
        if self.start <= self.stop:
            self.seqfeature = SeqFeature(FeatureLocation(start, stop),
                                         strand=strand, type=self.type)
        # Wrap-around genes will have stop < start
        else:
            self.seqfeature = SeqFeature(CompoundLocation(
                [FeatureLocation(start, self.genome_length),
                 FeatureLocation(0, stop)]), strand=strand, type=self.type)
        # Add feature qualifiers
        self.seqfeature.qualifiers = self.get_qualifiers()

    def get_qualifiers(self):
        """
        Helper function that uses tRNA data to populate the qualifiers
        attribute of `seqfeature`.
        :return: qualifiers OrderedDict()
        """
        qualifiers = OrderedDict()
        qualifiers["gene"] = [self.name]
        if self.locus_tag != "":
            qualifiers["locus_tag"] = [self.locus_tag]
        if self.amino_acid != "":
            if self.amino_acid in SPECIAL_AMINO_ACIDS:
                amino_acid = SPECIAL_AMINO_ACIDS[self.amino_acid]
            else:
                amino_acid = self.amino_acid
            qualifiers["product"] = [f"tRNA-{amino_acid}"]
            if self.anticodon != "":
                qualifiers["note"] = [f"tRNA-{self.amino_acid}"
                                      f"({self.anticodon})"]

        return qualifiers

    # TODO: create base feature class - use this version of the method unless
    #  Travis objects, because the logic is slightly better and it raises an
    #  Error if something is wrong with given seq.
    def set_nucleotide_sequence(self, seq=None, parent_genome_seq=None):
        """
        Set this feature's nucleotide sequence
        :param seq: sequence
        :type seq: str or Seq
        :param parent_genome_seq: parent genome sequence
        :type parent_genome_seq: Seq
        :raise: ValueError
        :return:
        """
        # If seq is given we'll try to use that
        if seq is not None:
            # If seq is a BioPython Seq object, use it directly
            if isinstance(seq, Seq):
                self.seq = seq.upper()
            # If seq is a string, we can coerce it to a BioPython Seq object
            elif isinstance(seq, str):
                self.seq = Seq(seq, IUPAC.ambiguous_dna)
            # If seq is something else entirely, we cannot use it
            else:
                raise ValueError(f"tRNA.seq cannot use type '{type(seq)}'")
        # If instead a parent_genome_seq is given, we'll try to use that
        elif parent_genome_seq is not None and self.seqfeature is not None:
            try:
                self.seq = self.seqfeature.extract(parent_genome_seq)
            # TODO: tighten exception clause, OR let the error pass to a higher
            #  level in the call stack
            except:
                # Leave as default
                pass
        # If neither is given, leave self.seq as default
        else:
            pass

    # TODO: create base feature class - fully equivalent to version in Cds,
    #  but documented differently.
    def set_nucleotide_length(self, use_seq=False):
        """
        Set the nucleotide length of this gene feature.
        :param use_seq: whether to use the Seq feature to calculate
        nucleotide length of this feature
        :type use_seq: bool
        :return:
        """
        if use_seq:
            self.length = len(self.seq)
        else:
            if self.coordinate_format == "0_half_open":
                # Python coordinate format
                self.length = self.stop - self.start
            elif self.coordinate_format == "1_closed":
                # Genbank coordinate format
                self.length = self.stop - self.start + 1
            else:
                # Unknown coordinate format
                self.length = -1


#     #This is an initial attempt at checking the tRNA product description
#     #Ultimately, a regular expression would be better to use
#     #tRNA product example = 'tRNA-Ser (AGC)'
#     #The code block below functions, but it does not fully account for
#     #tRNA-OTHER descriptions, tRNA-Stop descriptions,
#     #and it does not check the accuracy of
#     #the amino acid and anticodon pairing.
#     #The biggest problem is that the expected product and note descriptions
#     #are expected to change after they reach NCBI, so it is not clear
#     #how to best address that issue here, since nothing in the import
#     #table reflects WHERE the annotated genome came from.
#
#     product_error = 0
#
#     #product starts off as lowercase 'trna-ser (agc)'
#     #split1_list = 'trna' and 'ser (agc)'
#     tRNA_product_split1_list = product_field.split('-')
#
#     #If product is missing, an error will already have been thrown.
#     #The product should have a hypthen, so only parse product if it can be
#     #split into two elements.
#     if len(tRNA_product_split1_list) == 2:
#         tRNA_product_split1_prefix = tRNA_product_split1_list[0].strip() #'trna'
#
#         #split2_list = 'ser' and 'agc)'
#         tRNA_product_split2_list = tRNA_product_split1_list[1].split('(')
#         tRNA_product_amino_acid_three = tRNA_product_split2_list[0].strip() #'ser'
#
#         if tRNA_product_amino_acid_three != 'other' and \
#             tRNA_product_amino_acid_three != 'stop' and \
#             len(tRNA_product_amino_acid_three) != 3:
#                 product_error += 1
#
#         #The code block below checks for the presence of an anticodon.
#         #No need to use it currently, since there is so much variability
#         #at the tRNA product field, but retain for future use.
#         # if len(tRNA_product_split2_list) == 2:
#         #
#         #     #split3_list = 'agc' and ''
#         #     tRNA_product_split3_list = tRNA_product_split2_list[1].split(')')
#         #
#         #     #Only check the anticodon if the amino acid is NOT 'other'
#         #     if tRNA_product_amino_acid_three != 'other' and \
#         #         len(tRNA_product_split3_list) == 2:
#         #
#         #         tRNA_product_anticodon = tRNA_product_split3_list[0].strip() #'agc'
#         #         if len(tRNA_product_anticodon) != 3:
#         #             product_error += 1
#         #     else:
#         #         product_error += 1
#         #
#         # else:
#         #     product_error += 1
#     else:
#         product_error += 1
#
#     return product_error

# TODO the code below is pasted from import script.
#  It evaluates tRNA features, and needs to be implemented in the
#  check_tRNA function
# def check_trna_feature(feature):
#
#     #Now that start and stop have been parsed, check if coordinates are
#     fuzzy or not
#     if (tRNA_start.isdigit() and tRNA_stop.isdigit()):
#         tRNA_start = int(tRNA_start)
#         tRNA_stop = int(tRNA_stop)
#
#     else:
#         write_out(output_file,"\nError: tRNA starting at %s has fuzzy
#         coordinates in phage %s."\
#                 % (tRNA_start,phageName))
#         record_errors += 1
#         continue
#
#
#     if len(tRNA_seq) != tRNA_size:
#         write_out(output_file,"\nError: unable to retrieve sequence for
#         tRNA starting at %s in phage %s."\
#                 % (tRNA_start + 1,phageName))
#         record_errors += 1
#         continue
#
#     #Check to see if forward orientation terminal nucleotide is correct =
#     A or C
#     if tRNA_seq[-1] != 'A' and tRNA_seq[-1] != 'C':
#         record_warnings += 1
#         write_out(output_file,"\nWarning: tRNA starting at %s does not appear
#         to have correct terminal nucleotide in %s phage." \
#                 % (tRNA_start + 1,phageName))
#         record_errors += question("\nError: tRNA starting at %s has incorrect
#         terminal nucleotide in %s phage." \
#                 % (tRNA_start + 1,phageName))
#
#     if tRNA_size < 60 or tRNA_size > 100:
#         record_warnings += 1
#         write_out(output_file,"\nWarning: tRNA starting at %s does not appear
#         to be the correct size in %s phage."  \
#                 % (tRNA_start + 1,phageName))
#         record_errors += question("\nError: tRNA starting at %s is incorrect
#         size in %s phage." \
#                 % (tRNA_start + 1,phageName))
#
#         if len(tRNA_product) > 0:
#             if check_tRNA_product(tRNA_product) > 0:
#                 write_out(output_file,"\nError: tRNA starting at %s has
#                 incorrect amino acid or anticodon in %s." \
#                     % (tRNA_start + 1, phageName))
#                 record_errors += 1
#
#         else:
#             write_out(output_file,"\nError: tRNA starting at %s has incorrect
#             product in %s." \
#                 % (tRNA_start + 1, phageName))
#             record_errors += 1

# TODO revamp this code into a function
