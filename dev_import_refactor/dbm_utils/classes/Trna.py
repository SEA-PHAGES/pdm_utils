"""Represents a collection of data about a tRNA feature that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""






class TrnaFeature:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:

        #Datafields from Phamerator database:
        self.type_id = '' #Feature type: CDS, GenomeBoundary,or tRNA
        self.left_boundary = '' #Position of left boundary of feature, 0-indexed
        self.right_boundary = '' #Position of right boundary of feature, 0-indexed
        self.strand = '' #'forward', 'reverse', or 'NA'
        self.length = ''


        #Common to Phamerator
        self.phage_id = ''
        self.gene_id = '' #tRNA ID comprised of PhageID and Gene name
        self.gene_name = ''
        self.notes = ''
        self.search_notes = '' #non-generic descriptions




        #Common to NCBI
        self.locus_tag = '' #Gene ID comprised of PhageID and Gene name
        self.gene_number = ''
        self.product_description = ''
        self.processed_product_description = ''


















###
