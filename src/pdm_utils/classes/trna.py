"""Represents a collection of data about a tRNA feature that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""






class TrnaFeature:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:

        #Datafields from MySQL database:
        self.type = '' #Feature type: CDS, GenomeBoundary,or tRNA
        self.start = -1 #Position of left boundary of feature, 0-indexed
        self.stop = -1 #Position of right boundary of feature, 0-indexed
        self.orientation = '' #'forward', 'reverse', or 'NA'
        self.length = 0


        #Common to MySQL database
        self.genome_id = ''
        self.id = '' #tRNA ID comprised of PhageID and Gene name
        self.name = ''
        self.notes = ''
        self.search_notes = '' #non-generic descriptions

        #Common to NCBI
        self.locus_tag = '' #Gene ID comprised of PhageID and Gene name
        self.gene = ''
        self.raw_product = ''
        self.product = ''
