"""
A collection of objects to help parse and evaluate ClustalO-generated
Multiple Sequence Alignments (MSAs) and Percent Identity Matrices (PIMs).
"""
from Bio import AlignIO


class InitializationError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class PercentIdentityMatrix:
    """
    Class for handling Clustal Omega percent identity matrices.
    """
    def __init__(self, filepath):
        self.__filepath__ = filepath
        self.__node_count__ = 0
        self.__node_names__ = list()
        self.__matrix__ = list()
        self.__initialized__ = False

    def parse_matrix(self):
        """
        Parses pairwise percent identity matrix.
        :return:
        """
        with open(self.__filepath__, "r") as mat:
            # The first line of the matrix file is the number of genes
            num_nodes = int(mat.readline())
            names = list()
            rows = list()

            # All other lines build the matrix
            for line in mat:
                # Remove trailing whitespace and split on internal whitespace
                line = line.rstrip().split()
                names.append(line[0])
                rows.append([float(x) for x in line[1:]])

        # Sanity check that we parsed the matrix properly
        if num_nodes == len(names) == len(rows):
            self.node_count = num_nodes
            self.node_names = names
            self.matrix = rows
            self.initialized = True

    def check_initialization(self, caller):
        """
        Safe programming feature - raise an exception if a client
        tries to perform operations on uninitialized object.
        :param caller: name of the method that called this one
        """
        if not self.initialized:
            m = f"Cannot call method '{caller}' on uninitialized PercentIdentityMatrix"
            raise InitializationError(m)
        return

    @property
    def node_count(self):
        """
        Returns the number of nodes for this object.
        :return: self.__node_count__
        :rtype: int
        """
        self.check_initialization("node_count")
        return self.__node_count__

    @node_count.setter
    def node_count(self, value):
        """
        Sets the number of nodes for this object using the given
        value.
        :param value: number of nodes in the matrix
        :type value: int
        """
        if not isinstance(value, int):
            raise TypeError(f"Node count must be of type 'int'; '"
                            f"{type(value)}' is not valid.")

        self.__node_count__ = value

    @property
    def node_names(self):
        """
        Returns the node names for this object.
        :return: copy of self.__node_names__
        :rtype: list
        """
        self.check_initialization("node_names")
        return self.__node_names__[:]

    @node_names.setter
    def node_names(self, value):
        """
        Sets the node names for this object using the given value.
        :param value: node names for the genes in this object
        :type value: list
        """
        if not isinstance(value, list):
            raise TypeError(f"Node names must be of type 'list'; '"
                            f"{type(value)}' is not valid.")

        self.__node_names__ = value

    @property
    def matrix(self):
        """
        Returns the matrix for this object.
        :return: copy of self.__matrix__
        :rtype: list
        """
        self.check_initialization("matrix")
        return self.__matrix__[:]

    @matrix.setter
    def matrix(self, value):
        """
        Sets the matrix rows for this object using the given value.
        :param value: rows of pairwise identities for the genes in this
        object
        :type value: list
        """
        if not isinstance(value, list):
            raise TypeError(f"Matrix must be of type 'list'; '"
                            f"{type(value)}' is not valid.")

        self.__matrix__ = value

    @property
    def initialized(self):
        """
        Returns the initialization status for this object.
        :return: self.__initialized__
        :rtype: bool
        """
        return self.__initialized__

    @initialized.setter
    def initialized(self, value):
        """
        Sets the initialization state for this object using the given
        value.
        :param value: is the object initialized?
        :type value: bool
        """
        if not isinstance(value, bool):
            raise TypeError(f"Initialization value must be of type 'bool'; "
                            f"'{type(value)}' is not valid.")

        self.__initialized__ = value

    def get_nearest_neighbors(self, query, threshold):
        """
        Scans the query's row in the percent identity matrix to find
        all target nodes with values >= threshold. Returns the list
        in descending order (highest identity first).
        :param query: the anchoring geneid
        :type query: str
        :param threshold: the percent identity threshold for inclusion
        in the return list
        :type threshold: int, float
        :return: neighbors
        :rtype: list
        """
        self.check_initialization("get_nearest_neighbors")

        # Get the query row (or throw an error if the node name is invalid)
        query_index = self.node_names.index(query)
        query_row = self.get_row(query_index)

        # Iterate over the row - store node names where identity > threshold
        identities = list()
        keep_nodes = list()
        for i in range(len(query_row)):
            # Skip the self-match
            if i == query_index:
                continue
            identity = query_row[i]
            if identity >= threshold:
                identities.append(identity)
                keep_nodes.append(self.node_names[i])

        # Order from highest to lowest identity
        neighbors = list()
        while len(identities) > 0:
            next_index = identities.index(max(identities))
            identities.pop(next_index)
            neighbors.append(keep_nodes.pop(next_index))

        return neighbors

    def get_identity(self, query, target):
        """
        Returns the percent identity between a pair of genes in the
        matrix.
        :param query: row geneid
        :type query: str
        :param target: column geneid
        :type target: str
        :return: identity
        :rtype: float
        """
        q_index = self.node_names.index(query)
        t_index = self.node_names.index(target)
        identity = self.get_row(q_index)[t_index]
        return identity

    def get_distance(self, query, target):
        """
        Returns the distance between a pair of genes in the matrix.
        :param query: row geneid
        :type query: str
        :param target: column geneid
        :type target: str
        :return: distance
        :rtype: float
        """
        distance = 100 - self.get_identity(query, target)
        return distance

    def get_row(self, index):
        """
        Returns a copy of the requested matrix row.
        :param index: index of the row to be returned
        :return: copy of self.__matrix__[index]
        """
        self.check_initialization("get_row")
        return self.__matrix__[index][:]

    def get_centroid(self):
        """
        Calculates the average pairwise identity for every node in
        the matrix; returns the geneid and translation for the node
        with the highest average pairwise identity to the rest of the
        pham.
        :return: centroid
        :rtype: str
        """
        self.check_initialization("get_centroid")

        # Calculate each node's average pairwise identity
        apis = list()
        for i in range(self.node_count):
            geneid = self.node_names[i]
            apis.append(self.get_average_identity(geneid))

        # Find the index of the node with highest average pairwise identity
        centroid_index = apis.index(max(apis))

        centroid = self.node_names[centroid_index]

        return centroid

    def get_average_identity(self, geneid):
        """
        Calculates the average pairwise distance from the query geneid
        to every other member of the pham.
        :param geneid: the anchoring (query) gene
        :return:
        """
        self.check_initialization("get_average_distance")

        index = self.node_names.index(geneid)
        identities = self.get_row(index)
        identities.pop(index)     # Pop the self edge
        api = float(sum(identities))/len(identities)

        return api

    def get_false_positives(self, ident_thresh=30):
        """
        Finds and returns the nodes in the matrix whose best link to
        the pham is below the given threshold.
        :param ident_thresh: percent identity below which a node's
        inclusion in pham is a false positive
        :return: fps
        :rtype: set
        """
        fps = set()

        for i in range(self.node_count):
            identities = self.get_row(i).pop(i)         # Pop the self edge
            if max(identities) < ident_thresh:
                fps.add(self.node_names[i])

        return fps


class MultipleSequenceAlignment:
    """
    Class for handling Clustal Omega multiple sequence alignments.
    """
    def __init__(self, filepath, fmt="clustal"):
        self.__filepath__ = filepath
        self.__fmt__ = fmt
        self.__node_count__ = 0
        self.__node_names__ = list()
        self.__node_lengths__ = list()
        self.__alignment__ = None
        self.__conservation__ = None
        self.__initialized__ = False

    def parse_alignment(self):
        """
        Uses BioPython's AlignIO module to help parse alignment files
        :param fmt: alignment file format
        :type fmt: str
        """
        num_nodes = 0
        names = list()
        lengths = list()

        # Parse the alignment file into an alignment object
        aln = AlignIO.read(self.__filepath__, self.__fmt__)

        # Parse the alignment object
        for record in aln:
            num_nodes += 1
            names.append(record.id)
            lengths.append(len(str(record.seq).replace("-", "")))

        # Sanity check that we parsed the alignment properly
        if num_nodes == len(names) == len(lengths):
            self.node_count = num_nodes
            self.node_names = names
            self.node_lengths = lengths
            self.alignment = aln
            self.initialized = True

        # Populate alignment conservation data
        self.calculate_conservation()

    def check_initialization(self, caller):
        """
        Safe programming feature - raise an exception if a client
        tries to perform operations on uninitialized object.
        :param caller: name of the method that called this one
        """
        if not self.initialized:
            m = f"Cannot call method '{caller}' on uninitialized PercentIdentityMatrix"
            raise InitializationError(m)
        return

    @property
    def node_count(self):
        """
        Returns the number of nodes for this object.
        :return: self.__node_count__
        :rtype: int
        """
        self.check_initialization("node_count")
        return self.__node_count__

    @node_count.setter
    def node_count(self, value):
        """
        Sets the number of nodes for this object using the given
        value.
        :param value: number of nodes in the matrix
        :type value: int
        """
        if not isinstance(value, int):
            raise TypeError(f"Node count must be of type 'int'; '"
                            f"{type(value)}' is not valid.")

        self.__node_count__ = value

    @property
    def node_names(self):
        """
        Returns the node names for this object.
        :return: copy of self.__node_names__
        :rtype: list
        """
        self.check_initialization("node_names")
        return self.__node_names__[:]

    @node_names.setter
    def node_names(self, value):
        """
        Sets the node names for this object using the given value.
        :param value: node names for the genes in this object
        :type value: list
        """
        if not isinstance(value, list):
            raise TypeError(f"Node names must be of type 'list'; '"
                            f"{type(value)}' is not valid.")

        self.__node_names__ = value

    @property
    def node_lengths(self):
        """
        Returns the node lengths for this object
        :return: copy of self.__node_lengths__
        :rtype: list
        """
        self.check_initialization("node_lengths")
        return self.__node_lengths__[:]

    @node_lengths.setter
    def node_lengths(self, value):
        """
        Sets the node lengths for this object using the given value.
        :param value: node lengths for the genes in this object
        :type value: list
        """
        if not isinstance(value, list):
            raise TypeError(f"Node lengths must be of type 'list'; "
                            f"'{type(value)}' is not valid.")

        self.__node_lengths__ = value

    @property
    def alignment(self):
        """
        Returns the BioPython alignment object
        :return: self.__alignment__
        :rtype: BioPython Alignment
        """
        self.check_initialization("alignment")
        return self.__alignment__

    @alignment.setter
    def alignment(self, value):
        """
        Sets the alignment for this object using the given value.
        :param value: alignment for the genes in this object
        :type value: Bio.AlignIO.MultipleSeqAlignment
        """
        if not isinstance(value, AlignIO.MultipleSeqAlignment):
            raise TypeError("Aligment must be of type "
                            "'Bio.AlignIO.MultipleSeqAlignment'")

        self.__alignment__ = value

    @property
    def conservation(self):
        """
        Returns the alignment conservation dictionary
        :return: self.__conservation__
        """
        self.check_initialization("conservation")
        return self.__conservation__

    @conservation.setter
    def conservation(self, value):
        """
        Sets the alignment conservation dictionary with the given value
        :param value: conservation dictionary
        :type value: dict
        """
        if not isinstance(value, dict):
            raise TypeError(f"Conservation dictionary must be of type "
                            f"'dict'; '{type(value)}' is not valid.")

        self.__conservation__ = value

    @property
    def initialized(self):
        """
        Returns the initialization status for this object.
        :return: self.__initialized__
        :rtype: bool
        """
        return self.__initialized__

    @initialized.setter
    def initialized(self, value):
        """
        Sets the initialization state for this object using the given
        value.
        :param value: is the object initialized?
        :type value: bool
        """
        if not isinstance(value, bool):
            raise TypeError(f"Initialization value must be of type 'bool'; "
                            f"'{type(value)}' is not valid.")

        self.__initialized__ = value

    def calculate_conservation(self):
        """
        Uses the consensus string of a Clustal alignment to calculate
        the percent identity and percent similarity down alignment
        columns.
        """
        self.check_initialization("calculate_conservation")

        if self.__fmt__ != "clustal":
            "Conservation only supported for Clustal alignments."
            return

        consensus = str(self.__alignment__.column_annotations)
        stars = consensus.count("*")
        colons = consensus.count(":")
        periods = consensus.count(".")
        similars = stars + colons + periods

        identity = float(stars) / min(self.node_lengths) * 100
        similarity = float(similars) / min(self.node_lengths) * 100
        self.conservation = {"identity": identity, "similarity": similarity}

    def longest_gene(self):
        """
        Returns the geneid and translation of the longest sequence in
        the alignment.
        :return: (geneid, translation)
        """
        self.check_initialization("longest_gene")

        i = self.node_lengths.index(max(self.node_lengths))
        geneid = self.node_names[i]
        translation = self.alignment[geneid].seq.replace("-", "")

        return geneid, translation

    def get_sequence(self, geneid, gaps=True):
        """
        Returns the protein sequence of the record in the alignment
        whose id matches the given geneid
        :param geneid: search key
        :type geneid: str
        :param gaps: leave gaps in the sequence?
        :type gaps: bool
        :return: sequence
        :rtype: str
        """
        self.check_initialization("get_sequence")

        index = self.node_names.index(geneid)
        sequence = str(self.alignment[index].seq)

        if not gaps:
            sequence = sequence.replace("-", "")

        return sequence

    def get_false_positives(self, ident_thresh=30, cover_thresh=50):
        """
        Finds any nodes in the alignment that are false positives to
        defined percent identity and coverage thresholds.
        :param ident_thresh: percent identity cutoff for true positives
        :type ident_thresh: int, float
        :param cover_thresh: percent coverage cutoff for true positives
        :type cover_thresh: int, float
        :return: (identity_fps, coverage_fps)
        """
        self.check_initialization("get_false_positives")

        identity_fps = set()
        coverage_fps = set()

        min_len = min(self.node_lengths)
        max_len = max(self.node_lengths)

        # Identity checks need to be done on 2-member alignments because these
        # lack a percent identity matrix
        if self.node_count == 2:
            coverage = float(min_len) / max_len * 100

            if self.conservation["identity"] < ident_thresh:
                # Because every pham has by definition at least a single valid
                # member, we will say the shorter one is the false positive. If
                # genes are the same length, this is effectively random.
                fp_index = self.node_lengths.index(min_len)
                identity_fps.add(self.node_names[fp_index])
            if coverage < cover_thresh:
                fp_index = self.node_lengths.index(min_len)
                coverage_fps.add(self.node_names[fp_index])
        # No need to check identity for 3+ member phams, as these have a
        # percent identity matrix
        else:
            # Longest gene is anchor
            query_name, query_sequence = self.longest_gene()

            for index in range(self.node_count):
                target_name = self.node_names[index]
                if target_name == query_name:
                    continue
                target_sequence = self.get_sequence(target_name)
                coverage = calculate_coverage(query_sequence, target_sequence)
                if coverage < cover_thresh:
                    coverage_fps.add(target_name)

        return identity_fps, coverage_fps


class FastaMultipleSequence:
    def __init__(self, filepath):
        self.__filepath__ = filepath
        self.__node_count__ = 0
        self.__node_names__ = list()
        self.__sequences__ = list()
        self.__initialized__ = False

    def parse_fasta(self):
        """
        Parses .fasta file
        :return:
        """
        node_count = 0
        names = list()
        sequences = list()

        with open(self.__filepath__, "r") as fasta:
            lines = fasta.readlines()

        for line in lines:
            if line.startswith(">"):
                # header - grab geneid
                node_count += 1
                names.append(line.lstrip(">").rstrip())
            else:
                # sequence - grab it
                sequences.append(line.rstrip())

        if node_count == len(names) == len(sequences):
            self.node_count = node_count
            self.node_names = names
            self.sequences = sequences
            self.initialized = True

    def write_fasta(self):
        """
        Dumps this object to a fasta file
        :return:
        """
        self.check_initialization("write_fasta")

        with open(self.__filepath__, "w") as fasta:
            for i in range(self.node_count):
                header = self.node_names[i]
                sequence = self.sequences[i]
                fasta.write(f">{header}\n{sequence}\n")

    def check_initialization(self, caller):
        """
        Safe programming feature - raise an exception if a client
        tries to perform operations on uninitialized object.
        :param caller: name of the method that called this one
        """
        if not self.initialized:
            m = f"Cannot call method '{caller}' on uninitialized PercentIdentityMatrix"
            raise InitializationError(m)
        return

    @property
    def node_count(self):
        """
        Returns the number of nodes for this object.
        :return: self.__node_count__
        :rtype: int
        """
        self.check_initialization("node_count")
        return self.__node_count__

    @node_count.setter
    def node_count(self, value):
        """
        Sets the number of nodes for this object using the given
        value.
        :param value: number of nodes in the matrix
        :type value: int
        """
        if not isinstance(value, int):
            raise TypeError(f"Node count must be of type 'int'; '"
                            f"{type(value)}' is not valid.")

        self.__node_count__ = value

    @property
    def node_names(self):
        """
        Returns the node names for this object.
        :return: copy of self.__node_names__
        :rtype: list
        """
        self.check_initialization("node_names")
        return self.__node_names__[:]

    @node_names.setter
    def node_names(self, value):
        """
        Sets the node names for this object using the given value.
        :param value: node names for the genes in this object
        :type value: list
        """
        if not isinstance(value, list):
            raise TypeError(f"Node names must be of type 'list'; '"
                            f"{type(value)}' is not valid.")

        self.__node_names__ = value

    @property
    def sequences(self):
        """
        Returns the node names for this object.
        :return: copy of self.__node_names__
        :rtype: list
        """
        self.check_initialization("node_names")
        return self.__sequences__[:]

    @sequences.setter
    def sequences(self, value):
        """
        Sets the node names for this object using the given value.
        :param value: node names for the genes in this object
        :type value: list
        """
        if not isinstance(value, list):
            raise TypeError(f"Node names must be of type 'list'; '"
                            f"{type(value)}' is not valid.")

        self.__sequences__ = value

    @property
    def initialized(self):
        """
        Returns the initialization status for this object.
        :return: self.__initialized__
        :rtype: bool
        """
        return self.__initialized__

    @initialized.setter
    def initialized(self, value):
        """
        Sets the initialization state for this object using the given
        value.
        :param value: is the object initialized?
        :type value: bool
        """
        if not isinstance(value, bool):
            raise TypeError(f"Initialization value must be of type 'bool'; "
                            f"'{type(value)}' is not valid.")

        self.__initialized__ = value


def calculate_coverage(query, target, mode="shared_gapped"):
    """
    Can calculate target's percent coverage of the query in several
    ways. The simplest mode "length" strips gaps from both sequences
    and returns target_len / query_len * 100. The remaining modes will
    calculate coverage of only shared sequence, but differ in their
    handling of internal gaps. Both will trim the alignment length
    down to the first index at either end where neither sequence has
    a gap, and both will count gap-only columns against the alignment
    length. "shared_gapped" will ignore any other gaps in the middle
    of the alignment, while "shared_ungapped" will count both single-
    and double-sequence gaps as internal gaps. Both gapped modes
    calculate coverage as:

        shared_len = (trimmed_stop - trimmed_start + 1) - internal gaps
        coverage = float(shared_len) / ungapped_query_len * 100

    :param query: gapped sequence of query gene
    :param target: gapped sequence of target gene
    :param mode: how to calculate coverage
    :return: coverage
    """
    modes = ["shared_gapped", "shared_ungapped", "length"]
    if mode not in modes:
        print(f"Supported coverage modes: {modes[0]}, {modes[1]}, {modes[2]}")
        raise ValueError(f"Invalid mode for calculate_coverage(): '{mode}'")

    ungapped_query_len = len(query.replace("-", ""))
    ungapped_target_len = len(target.replace("-", ""))

    if mode == "length":
        coverage = float(ungapped_target_len) / ungapped_query_len * 100
    else:
        start = 0
        stop = len(query) - 1
        internal_gaps = 0

        # Trim N-terminus of alignment to the 1st column w/o gap in either seq
        for i in range(len(query)):
            if query[i] != "-" and target[i] != "-":
                break
            else:
                start += 1

        # Trim C-terminus of alignment to the 1st column w/o gap in either seq
        for i in reversed(range(len(query))):
            if query[i] != "-" and target[i] != "-":
                break
            else:
                stop -= 1

        for i in range(start, stop+1, 1):
            # Gap-only columns are always treated as internal gaps
            if query[i] == target[i] == "-":
                internal_gaps += 1
            # Half-gap columns only count if mode is shared_ungapped
            elif (mode == "shared_ungapped") and (query[i] == "-" or target[i] == "-"):
                internal_gaps += 1

        shared_len = (stop - start + 1) - internal_gaps
        coverage = float(shared_len) / ungapped_query_len * 100

    return coverage
