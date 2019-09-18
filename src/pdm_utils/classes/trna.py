"""Represents a collection of data about a tRNA feature that are commonly used
to maintain and update SEA-PHAGES phage genomics data.
"""






class TrnaFeature:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:

        #Datafields from Phamerator database:
        self.type = '' #Feature type: CDS, GenomeBoundary,or tRNA
        self.left = -1 #Position of left boundary of feature, 0-indexed
        self.right = -1 #Position of right boundary of feature, 0-indexed
        self.strand = '' #'forward', 'reverse', or 'NA'
        self.length = 0


        #Common to Phamerator
        self.genome_id = ''
        self.id = '' #tRNA ID comprised of PhageID and Gene name
        self.name = ''
        self.notes = ''
        self.search_notes = '' #non-generic descriptions




        #Common to NCBI
        self.locus_tag = '' #Gene ID comprised of PhageID and Gene name
        self.gene = ''
        self.product = ''
        self.processed_product = ''






























#TODO implement this function

# def check_tRNA_product(product_field):
#
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
#
#
#
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
#
#
#     else:
#         product_error += 1
#
#
#     return product_error











#TODO the code below is pasted from import script.
# It evaluates tRNA features, and needs to be implemented in the
# check_tRNA function

# def check_trna_feature(feature):
#
#     #Now that start and stop have been parsed, check if coordinates are fuzzy or not
#     if (tRNA_left.isdigit() and tRNA_right.isdigit()):
#         tRNA_left = int(tRNA_left)
#         tRNA_right = int(tRNA_right)
#
#     else:
#         write_out(output_file,"\nError: tRNA starting at %s has fuzzy coordinates in phage %s."\
#                 % (tRNA_left,phageName))
#         record_errors += 1
#         continue
#
#
#     if len(tRNA_seq) != tRNA_size:
#         write_out(output_file,"\nError: unable to retrieve sequence for tRNA starting at %s in phage %s."\
#                 % (tRNA_left + 1,phageName))
#         record_errors += 1
#         continue
#
#
#     #Check to see if forward strand terminal nucleotide is correct = A or C
#     if tRNA_seq[-1] != 'A' and tRNA_seq[-1] != 'C':
#         record_warnings += 1
#         write_out(output_file,"\nWarning: tRNA starting at %s does not appear to have correct terminal nucleotide in %s phage." \
#                 % (tRNA_left + 1,phageName))
#         record_errors += question("\nError: tRNA starting at %s has incorrect terminal nucleotide in %s phage." \
#                 % (tRNA_left + 1,phageName))
#
#     if tRNA_size < 60 or tRNA_size > 100:
#         record_warnings += 1
#         write_out(output_file,"\nWarning: tRNA starting at %s does not appear to be the correct size in %s phage."  \
#                 % (tRNA_left + 1,phageName))
#         record_errors += question("\nError: tRNA starting at %s is incorrect size in %s phage." \
#                 % (tRNA_left + 1,phageName))
#
#         if len(tRNA_product) > 0:
#             if check_tRNA_product(tRNA_product) > 0:
#                 write_out(output_file,"\nError: tRNA starting at %s has incorrect amino acid or anticodon in %s." \
#                     % (tRNA_left + 1, phageName))
#                 record_errors += 1
#
#         else:
#             write_out(output_file,"\nError: tRNA starting at %s has incorrect product in %s." \
#                 % (tRNA_left + 1, phageName))
#             record_errors += 1

#TODO revamp this code into a function







###
