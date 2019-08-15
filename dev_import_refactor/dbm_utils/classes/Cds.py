"""Represents a collection of data about a CDS features that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""

from functions import basic
from constants import constants
from classes import Eval
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re




class Cds:
    """Class to hold data about a CDS feature."""

    def __init__(self):

        # The following attributes are common to any CDS.

        # TODO: eventually change how id is computed.
        self.id = "" # Gene ID comprised of PhageID and Gene name
        self.name = "" # Tends to be an integer for SEA-PHAGES.
        self.genome_id = "" # Genome from which CDS feature is derived.
        self.seqfeature = None # Biopython SeqFeature object.
        self.left = "" # Genomic position
        self.right = "" # Genomic position
        self.coordinate_format = "" # Indexing format used for coordinates.
        self.strand = "" #'forward', 'reverse', 'top', 'bottom', etc.
        self.compound_parts = 0 # Number of regions that define the feature

        # TODO either implement the set_wrap() or get rid of this attribute.
        self.wrap = False # Does the feature wrap around the end of the genome?

        self.translation_table = ""
        self.translation = "" # Biopython amino acid Seq object.
        self._translation_length = 0
        self.seq = "" # Biopython nucleotide Seq object.
        self._length = 0


        # The following attributes are common to PhameratorDB.
        self.pham = "" # TODO build method to implement this.
        self.description = "" # Raw gene description
        self.processed_description = "" # Non-generic gene descriptions


        # The following attributes are common to
        # GenBank-formatted flat file records.
        self.locus_tag = "" # Gene ID comprised of PhageID and Gene name
        self._locus_tag_num = ""
        self.gene = "" # Tends to be an integer, but not guaranteed.
        self.product = ""
        self.function = ""
        self.note = ""
        self.processed_product = ""
        self.processed_function = ""
        self.processed_note = ""

        # The following attributes are usefule for processing data
        # from various data sources.
        self.evaluations = []
        self.type = ""
        self._left_right_strand_id = ()
        self._end_strand_id = ()
        self._start_end_id = ()




    def set_locus_tag(self, tag="", delimiter="_", check_value=None):
        """Set locus tag and split tag information."""
        self.locus_tag = tag
        if check_value is None:
            check_value = self.genome_id
        pattern = re.compile(check_value.lower())
        parts = tag.split(delimiter)

        index = 0
        found = False
        while (index < len(parts) and not found):
            search_result = pattern.search(parts[index].lower())
            if search_result != None:
                found = True
            else:
                index += 1
        if found:
            if index == len(parts) - 1:
                value = parts[index][search_result.end():]
            else:
                value = parts[-1]
        else:
            value = parts[-1]
        if value.isdigit():
            self._locus_tag_num = value



    def set_name(self, value=None):
        """Set the feature name.

        Ideally, the name of the CDS will be an integer. This information
        can be stored in multiple fields in the GenBank-formatted flat file.
        The name is first derived from the 'gene' qualifier, then
        the 'locus_tag' qualifier, and finally left empty.
        The 'value' parameter can be used to directly set this attribute
        regardless of the 'gene' and '_locus_tag_num' attributes."""

        # 1. PECAAN Draft:
        #    The 'gene' qualifier should be present and contain an integer.
        # 2. New SEA-PHAGES Final:
        #    The 'gene' qualifier should be present and contain an integer.
        # 3. SEA-PHAGES Final in GenBank:
        #    The 'gene' qualifier may or may not be present, and
        #    it may or may not have an integer.
        #    The 'locus_tag' qualifier may or may not be present,
        #    and may or may not have an integer.
        # 4. Non-SEA-PHAGES in GenBank:
        #    The 'gene' qualifier may or may not be present, and
        #    it may or may not have an integer.
        #    The 'locus_tag' qualifier may or may not be present,
        #    and may or may not have an integer.

        if value is not None:
            self.name = value
        elif self.gene != "":
            self.name = self.gene
        elif self._locus_tag_num != "":
            self.name = self._locus_tag_num
        else:
            self.name = ""


    def set_description(self, value):
        """Set the description and processed_description attributes from the
        selected attribute.
        """

        if value == "product":
            self.description = self.product
            self.processed_description = self.processed_product
        elif value == "function":
            self.description = self.function
            self.processed_description = self.processed_function
        elif value == "note":
            self.description = self.note
            self.processed_description = self.processed_note
        else:
            pass


    def get_start_end(self):
        """Return the coordinates in start-end format."""

        # Ensure format of strand info.
        strand = basic.reformat_strand(self.strand, "fr_long")

        if strand == "forward":
            start = self.left
            end = self.right
        elif strand == "reverse":
            start = self.right
            end = self.left
        else:
            start = -1
            end = -1

        return (start, end)

    def translate_seq(self):
        """Translate the CDS nucleotide sequence.

        Use Biopython to translate the nucleotide sequece.
        The method expects the nucleotide sequence to be a valid CDS
        sequence in which:
          1. it begins with a valid start codon,
          2. it ends with a stop codon,
          3. it contains only one stop codon,
          4. its length is divisible by 3,
          5. it translates non-standard start codons to methionine.
        If these criteria are not met, an empty Seq object is returned.
        """

        try:
            translation = self.seq.translate(table=self.translation_table,
                                             cds=True)
        except:
            translation = Seq("", IUPAC.protein)
        return translation


    def set_translation(self, value=None, translate=False):
        """Set translation and its length."""

        if isinstance(value, Seq):
            self.translation = value.upper()
        elif value is not None:
            try:
                self.translation = Seq(value.upper(), IUPAC.protein)
            except:
                self.translation = Seq("", IUPAC.protein)
        elif translate:
            self.translation = self.translate_seq()
        else:
            self.translation = Seq("", IUPAC.protein)
        self._translation_length = len(self.translation)


    def set_translation_table(self, value):
        """Set translation table integer."""

        try:
            self.translation_table = int(value)
        except:
            self.translation_table = -1


    def set_strand(self, value, format, case=False):
        """Sets strand based on indicated format."""
        self.strand = basic.reformat_strand(value, format, case)


    def set_wrap(self):
        """Determines if the feature wraps around the end of the genome.

        This method assumes that the left and right coordinates are set
        and reflect feature boundaries irrespective of strand."""

        if self.left > self.right:
            self.wrap = True
        else:
            self.wrap = False


    def set_location_id(self):
        """Create a tuple of feature location data.

        For left and right coordinates of the feature, it doesn't matter
        whether the feature is complex with a translational frameshift or not.
        Retrieving the "left" and "right" boundary attributes return the very
        beginning and end of the feature, disregarding the
        inner "join" coordinates.
        If only the feature "end" coordinate is used, strand information is
        required.
        If "start" and "end" coordinates are used instead of "left" and "right"
        coordinates, no strand information is required.
        """
        self._left_right_strand_id = (self.left, self.right, self.strand)

        start, end = self.get_start_end()
        self._end_strand_id = (end, self.strand)
        self._start_end_id = (start, end)


    def reformat_left_and_right(self, new_format):
        """Convert left and right coordinates to new coordinate format.
        This also updates the coordinate format attribute to reflect
        change. However, it does not update start and end attributes.
        """

        new_left, new_right = \
            basic.reformat_coordinates(self.left, \
                                       self.right, \
                                       self.coordinate_format, \
                                       new_format)

        if (new_left != "" and new_right != ""):
            self.left = new_left
            self.right = new_right
            self.coordinate_format = new_format


    def set_nucleotide_length(self, seq=False):
        """Set the length of the nucleotide sequence.

        Nucleotide length can be computed multiple ways.
        If the 'option' parameter is False, the 'left' and 'right'
        coordinates are used to determine the length (and take into
        account the 'coordinate_format' attribute. However, for
        compound features, this value may not be accurate. If the
        'option' parameter is True, the nucleotide sequence is used
        to determine the length. When the sequence reflects the
        entire feature (including compound features), the length
        will be accurate.
        """

        if seq:
            self._length = len(self.seq)
            pass
        else:
            if self.coordinate_format == "0_half_open":
                self._length = \
                    self.right - self.left
            elif self.coordinate_format == "1_closed":
                self._length = \
                    self.right - self.left + 1
            else:
                self._length = -1


    def set_nucleotide_sequence(self, value=None, parent_genome_seq=None):
        """Set the nucleotide sequence of the feature.

        This method can directly set the attribute from a supplied 'value',
        or it can retrieve the sequence from the parent genome using
        Biopython. In this latter case, it relies on a Biopython SeqFeature
        object for the sequence extraction method and coordinates.
        If this object was generated from a Biopython-parsed
        GenBank-formatted flat file, the coordinates are by default
        '0-based half-open', the object contains coordinates for every
        part of the feature (e.g. if it is a compound feature) and
        fuzzy locations. As a result, the retrieved sequence may not
        exactly match the length indicated from the 'left' and 'right'
        coordinates.
        """

        if isinstance(value, Seq):
            self.seq = value.upper()
        elif value is not None:
            try:
                self.seq = Seq(value.upper(), IUPAC.ambiguous_dna)
            except:
                self.seq = Seq("", IUPAC.ambiguous_dna)
        elif parent_genome_seq is not None:
            try:
                self.seq = self.seqfeature.extract(parent_genome_seq)
            except:
                # TODO not sure if this is the proper alphabet to use.
                self.seq = Seq("", IUPAC.ambiguous_dna)
        else:
            self.seq = Seq("", IUPAC.ambiguous_dna)


    def set_seqfeature(self):
        """Set the 'seqfeature' attribute.

        The 'seqfeature' attribute stores a Biopython SeqFeature object,
        which contains methods valuable to extracting sequence data
        relevant to the feature.
        """

        # SeqFeature methods rely on coordinates in 0-based half-open
        # format and strand to be numeric.
        new_left, new_right = \
            basic.reformat_coordinates(self.left, \
                                       self.right, \
                                       self.coordinate_format, \
                                       "0_half_open")

        new_strand = basic.reformat_strand(self.strand, "numeric")

        self.seqfeature = SeqFeature(FeatureLocation(new_left, new_right),
                                     strand=new_strand)







    # Evaluations.

    def check_translation_table(self, check_table=11):
        """Check that the translation table is correct."""

        if self.translation_table == check_table:
            result = "The translation table is correct."
            status = "correct"
        else:
            result = "The translation table is not correct."
            status = "error"
        definition = "Check that the translation table is correct."
        eval = Eval.Eval("CDS0001", definition, result, status)
        self.evaluations.append(eval)


    def check_translation_length(self):
        """Confirm that a translation is present."""
        if self._translation_length < 1:
            result = "A translation is not present."
            status = "error"
        else:
            result = "A translation is present."
            status = "correct"

        definition = "Check that there is a translation present."
        eval = Eval.Eval("CDS0002", definition, result, status)
        self.evaluations.append(eval)


    def check_translation(self):
        """Check that the current and expected translations match."""

        translation = self.translate_seq()
        if self._translation_length < len(translation):
            result = "The translation length is shorter than expected."
            status = "error"
        elif self._translation_length > len(translation):
            result = "The translation length is longer than expected."
            status = "error"
        elif self.translation != translation:
            result = "The translation is different than expected."
            status = "error"
        else:
            result = "The translation is correct."
            status = "correct"
        definition = "Check that the feature contains the expected translation."
        eval = Eval.Eval("CDS0003", definition, result, status)
        self.evaluations.append(eval)


    def check_amino_acids(self, check_set=set()):
        """Check whether all amino acids in the translation are valid."""
        amino_acid_set = set(self.translation)
        amino_acid_error_set = amino_acid_set - check_set

        if len(amino_acid_error_set) > 0:
            result = "There are unexpected amino acids in the translation: " \
                + str(amino_acid_error_set)
            status = "error"
        else:
            result = "There are no unexpected amino acid residues."
            status = "correct"

        definition = "Check validity of amino acid residues."
        eval = Eval.Eval("CDS0004", definition, result, status = status)
        self.evaluations.append(eval)


    def check_strand(self, format="fr_short", case=True):
        """Check if strand is set appropriately."""
        expected_strand = basic.reformat_strand(self.strand,
                                                format=format,
                                                case=case)
        if self.strand == expected_strand:
            result = "The feature strand is correct."
            status = "correct"
        else:
            result = "The feature strand is not correct."
            status = "error"
        definition = "Check if the strand is set appropriately."
        eval = Eval.Eval("CDS0005", definition, result, status)
        self.evaluations.append(eval)


    def check_coordinates(self):
        """Check if coordinates are exact.

        This method assumes that if the coordinates are not exact, they
        have been set to -1 or are not integers.
        """

        if not (str(self.left).isdigit() and str(self.right).isdigit()):
            result = "The feature coordinates are not determined: " \
                + str((self.left, self.right))
            status = "error"
        elif (self.left == -1 or self.right == -1):
            # TODO unit test this elif clause.
            result = "The feature coordinates are not determined: " \
                + str((self.left, self.right))
            status = "error"
        else:
            result = "Feature coordinates are exact."
            status = "correct"

        definition = "Check if the left and right boundary coordinates " + \
                        "are exact or fuzzy."
        eval = Eval.Eval("CDS0006", definition, result, status)
        self.evaluations.append(eval)


    def check_locus_tag_present(self, expect=True):
        """Check if status of locus tag matches expectations."""

        if self.locus_tag != "":
            present = True
        else:
            present = False

        if expect:
            if present:
                result = "The locus_tag qualifier is present."
                status = "correct"
            else:
                result = "The locus_tag qualifier is not present."
                status = "error"
        else:
            if present:
                result = "The locus_tag qualifier is present."
                status = "error"
            else:
                result = "The locus_tag qualifier is not present."
                status = "correct"
        definition = "Check if the locus_tag qualifier status is expected."
        eval = Eval.Eval("CDS0007", definition, result, status)
        self.evaluations.append(eval)


    def check_locus_tag_structure(self, check_value=None, only_typo=False,
                                  prefix_set=set(), caps=True):
        """Check if the locus_tag is structured correctly.

        The 'check_value' parameter provides the genome ID that is expected
        to be present. If None, the 'genome_id' parameter is used by default.
        With the 'only_typo' parameter set to True, a simpler
        structure analysis only checks for whether the genome ID is
        present.
        If a prefix is expected, the 'prefix_set' parameter indicates
        the set of possible prefixes that are expected.
        The 'caps' parameters indicates whether the locus_tag is expected
        to be completely capitalized or not.
        """
        if check_value is None:
            check_value = self.genome_id

        results = []
        if only_typo:
            pattern = re.compile(check_value.lower())
            search_result = pattern.search(self.locus_tag.lower())
            if search_result == None:
                results.append("The genome ID is missing.")
        else:
            # Expected structure: SEA_TRIXIE_20
            parts = self.locus_tag.split("_")
            if caps:
                if self.locus_tag != self.locus_tag.upper():
                    results.append("The capitalization is incorrect.")
            if len(parts) == 3:
                if prefix_set is not None:
                    if parts[0].upper() not in prefix_set:
                        results.append("The prefix is missing.")
                if parts[1].upper() != check_value.upper():
                    results.append("The genome ID is missing.")
                if not parts[2].isdigit():
                    results.append("The feature number is missing.")
            else:
                results.append("The locus_tag does not contain three parts.")
        if len(results) == 0:
            result = "The locus_tag qualifier is structured correctly."
            status = "correct"
        else:
            result = "The locus_tag qualifier is not structured correctly." \
                     + " ".join(results)
            status = "error"
        definition = "Check if the locus_tag qualifier is structured correctly."
        eval = Eval.Eval("CDS", definition, result, status)
        self.evaluations.append(eval)











    # TODO is this needed?
    # TODO implement.
    # TODO unittest
    def check_id_typo(self, check_value=None):
        """Check if the id contains potential typos."""

        if check_value is None:
            check_value = self.id

        pattern = re.compile(check_value.lower())
        search_result = pattern.search(self.id.lower())

        if search_result == None:
            result = "The id has a typo."
            status = "error"
        else:
            result = "The id is correct."
            status = "correct"
        definition = "Check if the id contains a typo."
        eval = Eval.Eval("CDS0009", definition, result, status)
        self.evaluations.append(eval)


    def check_gene_present(self, expect=True):
        """Check if the status of gene matches expectations."""

        if self.gene != "":
            present = True
        else:
            present = False

        if expect:
            if present:
                result = "The gene qualifier is present."
                status = "correct"
            else:
                result = "The gene qualifier is not present."
                status = "error"
        else:
            if present:
                result = "The gene qualifier is present."
                status = "error"
            else:
                result = "The gene qualifier is not present."
                status = "correct"
        definition = "Check if the gene status is expected."
        eval = Eval.Eval("CDS", definition, result, status)
        self.evaluations.append(eval)


    def check_gene_structure(self):
        """Check if the gene qualifier contains an integer."""

        try:
            value = int(self.gene)
        except:
            value = self.gene

        if isinstance(value, int):
            result = "The gene qualifier contains an integer."
            status = "correct"
        else:
            result = "The gene qualifier does not contain an integer."
            status = "error"
        definition = "Check if the gene qualifier contains an integer."
        eval = Eval.Eval("CDS", definition, result, status)
        self.evaluations.append(eval)


    def check_compatible_gene_and_locus_tag(self):
        """Check if the gene and locus_tag attributes contain the same
        gene number."""
        if self.gene == self._locus_tag_num:
            result = "The gene and locus_tag numbers are consistent."
            status = "correct"
        else:
            result = "The gene and locus_tag numbers are not consistent."
            status = "error"
        definition = "Check if the gene and locus_tag numbers are consistent."
        eval = Eval.Eval("CDS", definition, result, status)
        self.evaluations.append(eval)


    def check_description_field(self, attribute="product"):
        """Check if there are CDS descriptions in unexpected fields.

        This method evaluates if the indicated field is empty or generic,
        and other fields contain non-generic data.
        """

        description = ""
        description_set = set()

        if attribute == "product":
            description = self.processed_product
        else:
            if self.processed_product != "":
                description_set.add(self.processed_product)

        if attribute == "function":
            description = self.processed_function
        else:
            if self.processed_function != "":
                description_set.add(self.processed_function)

        if attribute == "note":
            description = self.processed_note
        else:
            if self.processed_note != "":
                description_set.add(self.processed_note)

        if (description == "" and len(description_set) > 0):
            result = "The description is not in the expected field."
            status = "error"
        else:
            result = "The description is in the expected field."
            status = "correct"

        definition = \
            "Check if there is a discrepancy between description fields."
        eval = Eval.Eval("CDS0010", definition, result, status)
        self.evaluations.append(eval)


    def check_generic_data(self, attribute=None):
        """Check if the indicated attribute contains generic data."""

        if attribute == "product":
            original = self.product
            processed = self.processed_product
        elif attribute == "function":
            original = self.function
            processed = self.processed_function
        elif attribute == "note":
            original = self.note
            processed = self.processed_note
        else:
            original = ""
            processed = ""

        if original == processed:
            result = "The '%s' field is correct." % attribute
            status = "correct"
        else:
            result = "The '%s' field is not correct." % attribute
            status = "error"
        definition = "Check if the '%s' field contains generic data." % attribute
        eval = Eval.Eval("CDS0011", definition, result, status)
        self.evaluations.append(eval)





    # TODO implement.
    # TODO unittest.
    def check_valid_description(self, check_set=None, attribute=None):
        """Check if the CDS description in the indicated attribute is valid."""
        pass




###
