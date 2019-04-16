"""Represents a structure to directly compare data between two or more genomes."""



class MatchedGenomes:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.ticket = ''
        self.matched_genomes_dict = {}
        self.__list_of_errors = [] #Todo not sure if this is list of messages or error objects
        self.__list_of_sql_queries = [] #Todo not sure if this should be here or not















#Below: old python2 code

class MatchedGenomes:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.__phamerator_genome = ''
        self.__phagesdb_genome = ''
        self.__ncbi_genome = ''

        #Phamerator and NCBI matched data comparison results
        self.__phamerator_ncbi_sequence_mismatch = False
        self.__phamerator_ncbi_sequence_length_mismatch = False
        self.__ncbi_record_header_fields_phage_name_mismatch = False
        self.__ncbi_host_mismatch = False
        self.__phamerator_ncbi_perfect_matched_features = [] #List of MatchedCdsFeature objects
        self.__phamerator_ncbi_imperfect_matched_features = [] #List of MatchedCdsFeature objects
        self.__phamerator_features_unmatched_in_ncbi = [] #List of CdsFeature objects
        self.__ncbi_features_unmatched_in_phamerator = [] #List of CdsFeature objects
        self.__phamerator_ncbi_perfect_matched_features_tally = 0
        self.__phamerator_ncbi_imperfect_matched_features_tally = 0
        self.__phamerator_features_unmatched_in_ncbi_tally = 0
        self.__ncbi_features_unmatched_in_phamerator_tally = 0
        self.__phamerator_ncbi_different_descriptions_tally = 0
        self.__phamerator_ncbi_different_translations_tally = 0
        self.__ph_ncbi_author_error = False





        #Phamerator and phagesdb matched data comparison results
        self.__phamerator_phagesdb_sequence_mismatch = False
        self.__phamerator_phagesdb_sequence_length_mismatch = False
        self.__phamerator_phagesdb_host_mismatch = False
        self.__phamerator_phagesdb_accession_mismatch = False
        self.__phamerator_phagesdb_cluster_subcluster_mismatch = False


        #phagesdb and NCBI matched data comparison results
        self.__phagesdb_ncbi_sequence_mismatch = False
        self.__phagesdb_ncbi_sequence_length_mismatch = False

        #Total errors summary
        self.__contains_errors = False
        self.__total_number_genes_with_errors = 0



    # Define all attribute setters:
    def set_phamerator_genome(self,value):
        self.__phamerator_genome = value
    def set_phagesdb_genome(self,value):
        self.__phagesdb_genome = value
    def set_ncbi_genome(self,value):
        self.__ncbi_genome = value

    def compare_phamerator_ncbi_genomes(self):

        #verify that there is a Phamerator and NCBI genome in the matched genome object
        ph_genome = self.__phamerator_genome
        ncbi_genome = self.__ncbi_genome
        ph_cds_list = ph_genome.get_cds_features()

        if isinstance(ph_genome,PhameratorGenome) and isinstance(ncbi_genome,NcbiGenome):

            if ph_genome.get_sequence() != ncbi_genome.get_sequence():
                self.__phamerator_ncbi_sequence_mismatch = True
            if ph_genome.get_length() != ncbi_genome.get_length():
                self.__phamerator_ncbi_sequence_length_mismatch = True

            #Compare phage names
            pattern1 = re.compile('^' + ph_genome.get_phage_name() + '$')
            pattern2 = re.compile('^' + ph_genome.get_phage_name())

            if find_name(pattern2,ncbi_genome.get_record_description().split(' ')) == 0 or \
                find_name(pattern1,ncbi_genome.get_record_source().split(' ')) == 0 or \
                find_name(pattern1,ncbi_genome.get_record_organism().split(' ')) == 0 or \
                find_name(pattern1,ncbi_genome.get_source_feature_organism().split(' ')) == 0:

                self.__ncbi_record_header_fields_phage_name_mismatch = True

            #Compare host data
            search_host = ph_genome.get_host()
            if search_host == 'Mycobacterium':
                search_host = search_host[:-3]
            pattern3 = re.compile('^' + search_host)

            if (find_name(pattern3,ncbi_genome.get_record_description().split(' ')) == 0 or \
                find_name(pattern3,ncbi_genome.get_record_source().split(' ')) == 0 or \
                find_name(pattern3,ncbi_genome.get_record_organism().split(' ')) == 0 or \
                find_name(pattern3,ncbi_genome.get_source_feature_organism().split(' ')) == 0) or \
                (ncbi_genome.get_source_feature_host() != '' and find_name(pattern3,ncbi_genome.get_source_feature_host().split(' ')) == 0) or \
                (ncbi_genome.get_source_feature_lab_host() != '' and find_name(pattern3,ncbi_genome.get_source_feature_lab_host().split(' ')) == 0):

                self.__ncbi_host_mismatch = True

            #Check author list for errors
            #For genomes with AnnotationAuthor = 1 (Hatfull), Graham is expected
            #to be an author.
            #For genomes with AnnotationAuthor = 0 (non-Hatfull/Genbank), Graham
            #is NOT expected to be an author.
            pattern5 = re.compile('hatfull')
            search_result = pattern5.search(ncbi_genome.get_record_authors().lower())
            if ph_genome.get_annotation_author() == 1 and search_result == None:
                self.__ph_ncbi_author_error = True
            elif ph_genome.get_annotation_author() == 0 and search_result != None:
                self.__ph_ncbi_author_error = True
            else:
                #Any other combination of phamerator and ncbi author can be skipped
                pass


            #Compare CDS features

            #First find all unique start-end-strand cds identifiers for phamerator and ncbi genomes
            ph_start_end_strand_id_set = set()
            ph_start_end_strand_duplicate_id_set = set() #All end_strand ids that are not unique
            for cds in ph_cds_list:
                if cds.get_start_end_strand_id() not in ph_start_end_strand_id_set:
                    ph_start_end_strand_id_set.add(cds.get_start_end_strand_id())
                else:
                    ph_start_end_strand_duplicate_id_set.add(cds.get_start_end_strand_id())
            #Remove the duplicate end_strand ids from the main id_set
            ph_start_end_strand_id_set = ph_start_end_strand_id_set - ph_start_end_strand_duplicate_id_set


            ncbi_cds_list = ncbi_genome.get_cds_features()
            ncbi_start_end_strand_id_set = set()
            ncbi_start_end_strand_duplicate_id_set = set() #All end_strand ids that are not unique
            for cds in ncbi_cds_list:
                if cds.get_start_end_strand_id() not in ncbi_start_end_strand_id_set:
                    ncbi_start_end_strand_id_set.add(cds.get_start_end_strand_id())

                else:
                    ncbi_start_end_strand_duplicate_id_set.add(cds.get_start_end_strand_id())
            #Remove the duplicate end_strand ids from the main id_set
            ncbi_start_end_strand_id_set = ncbi_start_end_strand_id_set - ncbi_start_end_strand_duplicate_id_set

            #Create the perfect matched and unmatched sets
            ph_unmatched_start_end_strand_id_set = ph_start_end_strand_id_set - ncbi_start_end_strand_id_set
            ncbi_unmatched_start_end_strand_id_set = ncbi_start_end_strand_id_set - ph_start_end_strand_id_set
            perfect_matched_cds_id_set = ph_start_end_strand_id_set & ncbi_start_end_strand_id_set


            #From the unmatched sets, created second round of end-strand id sets
            #Phamerator end_strand data
            ph_end_strand_id_set = set()
            ph_end_strand_duplicate_id_set = set() #All end_strand ids that are not unique
            for cds in ph_cds_list:
                if cds.get_start_end_strand_id() in ph_unmatched_start_end_strand_id_set:
                    if cds.get_end_strand_id() not in ph_end_strand_id_set:
                        ph_end_strand_id_set.add(cds.get_end_strand_id())
                    else:
                        ph_end_strand_duplicate_id_set.add(cds.get_end_strand_id())

            #Remove the duplicate end_strand ids from the main id_set
            ph_end_strand_id_set = ph_end_strand_id_set - ph_end_strand_duplicate_id_set


            ncbi_end_strand_id_set = set()
            ncbi_end_strand_duplicate_id_set = set() #All end_strand ids that are not unique
            for cds in ncbi_cds_list:
                if cds.get_start_end_strand_id() in ncbi_unmatched_start_end_strand_id_set:
                    if cds.get_end_strand_id() not in ncbi_end_strand_id_set:
                        ncbi_end_strand_id_set.add(cds.get_end_strand_id())
                    else:
                        ncbi_end_strand_duplicate_id_set.add(cds.get_end_strand_id())
            #Remove the duplicate end_strand ids from the main id_set
            ncbi_end_strand_id_set = ncbi_end_strand_id_set - ncbi_end_strand_duplicate_id_set


            #Create the imperfect matched set
            imperfect_matched_cds_id_set = ph_end_strand_id_set & ncbi_end_strand_id_set


            #Now go back through all cds features and assign to the appropriate dictionary or list
            ph_perfect_matched_cds_dict = {}
            ph_imperfect_matched_cds_dict = {}
            ph_unmatched_cds_list = []


            for cds in ph_cds_list:
                if cds.get_start_end_strand_id() in perfect_matched_cds_id_set:
                    ph_perfect_matched_cds_dict[cds.get_start_end_strand_id()] = cds
                elif cds.get_end_strand_id() in imperfect_matched_cds_id_set:
                    ph_imperfect_matched_cds_dict[cds.get_end_strand_id()] = cds
                else:
                    ph_unmatched_cds_list.append(cds)


            ncbi_perfect_matched_cds_dict = {}
            ncbi_imperfect_matched_cds_dict = {}
            ncbi_unmatched_cds_list = []

            for cds in ncbi_cds_list:
                if cds.get_start_end_strand_id() in perfect_matched_cds_id_set:
                    ncbi_perfect_matched_cds_dict[cds.get_start_end_strand_id()] = cds
                elif cds.get_end_strand_id() in imperfect_matched_cds_id_set:
                    ncbi_imperfect_matched_cds_dict[cds.get_end_strand_id()] = cds
                else:
                    ncbi_unmatched_cds_list.append(cds)


            #Create MatchedCdsFeatures objects
            #Compute matched gene errors and tallies
            #Perfectly matched features
            for start_end_strand_tup in perfect_matched_cds_id_set:

                matched_cds_object = MatchedCdsFeatures()
                matched_cds_object.set_phamerator_feature(ph_perfect_matched_cds_dict[start_end_strand_tup])
                matched_cds_object.set_ncbi_feature(ncbi_perfect_matched_cds_dict[start_end_strand_tup])
                matched_cds_object.compare_phamerator_ncbi_cds_features()

                if matched_cds_object.get_total_errors() > 0:
                    self.__total_number_genes_with_errors += 1



                if matched_cds_object.get_phamerator_ncbi_different_translations():
                    self.__phamerator_ncbi_different_translations_tally += 1
                if matched_cds_object.get_phamerator_ncbi_different_descriptions():
                    self.__phamerator_ncbi_different_descriptions_tally += 1
                self.__phamerator_ncbi_perfect_matched_features.append(matched_cds_object)


            #Imperfectly matched features
            for end_strand_tup in imperfect_matched_cds_id_set:

                matched_cds_object = MatchedCdsFeatures()
                matched_cds_object.set_phamerator_feature(ph_imperfect_matched_cds_dict[end_strand_tup])
                matched_cds_object.set_ncbi_feature(ncbi_imperfect_matched_cds_dict[end_strand_tup])
                matched_cds_object.compare_phamerator_ncbi_cds_features()

                if matched_cds_object.get_total_errors() > 0:
                    self.__total_number_genes_with_errors += 1

                if matched_cds_object.get_phamerator_ncbi_different_descriptions():
                    self.__phamerator_ncbi_different_descriptions_tally += 1
                self.__phamerator_ncbi_imperfect_matched_features.append(matched_cds_object)


            #Compute unmatched error and gene total errors for all unmatched features
            for cds in ph_unmatched_cds_list:
                cds.set_unmatched_error()
                cds.compute_total_cds_errors()
                if cds.get_total_errors() > 0:
                    self.__total_number_genes_with_errors += 1

            for cds in ncbi_unmatched_cds_list:
                cds.set_unmatched_error()
                cds.compute_total_cds_errors()
                if cds.get_total_errors() > 0:
                    self.__total_number_genes_with_errors += 1




            #Set unmatched cds lists
            self.__phamerator_features_unmatched_in_ncbi = ph_unmatched_cds_list
            self.__ncbi_features_unmatched_in_phamerator = ncbi_unmatched_cds_list

            #Now compute the number of features in each category
            self.__phamerator_ncbi_perfect_matched_features_tally = len(self.__phamerator_ncbi_perfect_matched_features)
            self.__phamerator_ncbi_imperfect_matched_features_tally = len(self.__phamerator_ncbi_imperfect_matched_features)
            self.__phamerator_features_unmatched_in_ncbi_tally = len(self.__phamerator_features_unmatched_in_ncbi)
            self.__ncbi_features_unmatched_in_phamerator_tally = len(self.__ncbi_features_unmatched_in_phamerator)


        #If there is no matching NCBI genome, assign all Phamerator genes to Unmatched
        else:

            #Set unmatched cds lists, but do NOT count them in the unmatched tally.
            #The unmatched tally should reflect unmatched genes if there is actually a metching NCBI genome.
            self.__phamerator_features_unmatched_in_ncbi = ph_cds_list

            #Now that all errors have been computed for each gene, compute how many genes have errors
            #If there is no matching NCBI genome, gene error tallies are only computed for the Phamerator genome
            ph_genome.compute_genes_with_errors_tally()
            self.__total_number_genes_with_errors = ph_genome.get_genes_with_errors_tally()




    def compare_phamerator_phagesdb_genomes(self):

        #verify that there is a Phamerator and phagesdb genome in the matched genome object
        ph_genome = self.__phamerator_genome
        pdb_genome = self.__phagesdb_genome

        if isinstance(ph_genome,PhameratorGenome) and isinstance(pdb_genome,PhagesdbGenome):

            if ph_genome.get_sequence() != pdb_genome.get_sequence():
                self.__phamerator_phagesdb_sequence_mismatch = True
            if ph_genome.get_length() != pdb_genome.get_length():
                self.__phamerator_phagesdb_sequence_length_mismatch = True
            if ph_genome.get_accession() != pdb_genome.get_accession():
                self.__phamerator_phagesdb_accession_mismatch = True
            if ph_genome.get_host() != pdb_genome.get_host():
                self.__phamerator_phagesdb_host_mismatch = True
            if ph_genome.get_cluster_subcluster() != pdb_genome.get_cluster() and \
                ph_genome.get_cluster_subcluster() != pdb_genome.get_subcluster():

                self.__phamerator_phagesdb_cluster_subcluster_mismatch = True


    def compare_phagesdb_ncbi_genomes(self):

        #verify that there is a phagesdb and NCBI genome in the matched genome object
        pdb_genome = self.__phagesdb_genome
        ncbi_genome = self.__ncbi_genome

        if isinstance(pdb_genome,PhagesdbGenome) and isinstance(ncbi_genome,NcbiGenome):
            if pdb_genome.get_sequence() != ncbi_genome.get_sequence():
                self.__phagesdb_ncbi_sequence_mismatch = True
            if pdb_genome.get_length() != ncbi_genome.get_length():
                self.__phagesdb_ncbi_sequence_length_mismatch = True



    def compute_total_genome_errors(self):

        ph_genome = self.__phamerator_genome
        ncbi_genome = self.__ncbi_genome
        pdb_genome = self.__phagesdb_genome


        if self.__total_number_genes_with_errors > 0:
            self.__contains_errors = True


        #Phamerator genome
        if ph_genome.get_nucleotide_errors():
            self.__contains_errors = True
        if ph_genome.get_status_description_error():
            self.__contains_errors = True
        if ph_genome.get_status_accession_error():
            self.__contains_errors = True


        #NCBI genome
        if isinstance(ncbi_genome,NcbiGenome):
            if ncbi_genome.get_nucleotide_errors():
                self.__contains_errors = True



        #phagesdb genome
        if isinstance(pdb_genome,PhagesdbGenome):
            if pdb_genome.get_nucleotide_errors():
                self.__contains_errors = True


        #Phamerator-NCBI
        if self.__phamerator_ncbi_sequence_mismatch:
            self.__contains_errors = True
        if self.__phamerator_ncbi_sequence_length_mismatch:
            self.__contains_errors = True
        if self.__ncbi_record_header_fields_phage_name_mismatch:
            self.__contains_errors = True
        if self.__ncbi_host_mismatch:
            self.__contains_errors = True
        if self.__phamerator_ncbi_imperfect_matched_features_tally > 0:
            self.__contains_errors = True
        if self.__phamerator_features_unmatched_in_ncbi_tally > 0:
            self.__contains_errors = True
        if self.__ncbi_features_unmatched_in_phamerator_tally > 0:
            self.__contains_errors = True
        if self.__phamerator_ncbi_different_descriptions_tally > 0:
            self.__contains_errors = True
        if self.__phamerator_ncbi_different_translations_tally > 0:
            self.__contains_errors = True
        if self.__ph_ncbi_author_error:
            self.__contains_errors = True


        #Phamerator-phagesdb
        if self.__phamerator_phagesdb_sequence_mismatch:
            self.__contains_errors = True
        if self.__phamerator_phagesdb_sequence_length_mismatch:
            self.__contains_errors = True
        if self.__phamerator_phagesdb_host_mismatch:
            self.__contains_errors = True
        if self.__phamerator_phagesdb_accession_mismatch:
            self.__contains_errors = True
        if self.__phamerator_phagesdb_cluster_subcluster_mismatch:
            self.__contains_errors = True

        #phagesdb-NCBI
        if self.__phagesdb_ncbi_sequence_mismatch:
            self.__contains_errors = True
        if self.__phagesdb_ncbi_sequence_length_mismatch:
            self.__contains_errors = True


    # Define all attribute getters:
    def get_phamerator_genome(self):
        return self.__phamerator_genome
    def get_phagesdb_genome(self):
        return self.__phagesdb_genome
    def get_ncbi_genome(self):
        return self.__ncbi_genome
    def get_phamerator_ncbi_sequence_mismatch(self):
        return self.__phamerator_ncbi_sequence_mismatch
    def get_phamerator_ncbi_perfect_matched_features(self):
        return self.__phamerator_ncbi_perfect_matched_features
    def get_phamerator_ncbi_imperfect_matched_features(self):
        return self.__phamerator_ncbi_imperfect_matched_features
    def get_phamerator_features_unmatched_in_ncbi(self):
        return self.__phamerator_features_unmatched_in_ncbi
    def get_ncbi_features_unmatched_in_phamerator(self):
        return self.__ncbi_features_unmatched_in_phamerator
    def get_phamerator_ncbi_perfect_matched_features_tally(self):
        return self.__phamerator_ncbi_perfect_matched_features_tally
    def get_phamerator_ncbi_imperfect_matched_features_tally(self):
        return self.__phamerator_ncbi_imperfect_matched_features_tally
    def get_phamerator_features_unmatched_in_ncbi_tally(self):
        return self.__phamerator_features_unmatched_in_ncbi_tally
    def get_ncbi_features_unmatched_in_phamerator_tally(self):
        return self.__ncbi_features_unmatched_in_phamerator_tally
    def get_phamerator_ncbi_sequence_length_mismatch(self):
        return self.__phamerator_ncbi_sequence_length_mismatch
    def get_ncbi_record_header_fields_phage_name_mismatch(self):
        return self.__ncbi_record_header_fields_phage_name_mismatch
    def get_ncbi_host_mismatch(self):
        return self.__ncbi_host_mismatch
    def get_phamerator_ncbi_different_descriptions_tally(self):
        return self.__phamerator_ncbi_different_descriptions_tally
    def get_phamerator_ncbi_different_translation_tally(self):
        return self.__phamerator_ncbi_different_translations_tally

    def get_phagesdb_ncbi_sequence_mismatch(self):
        return self.__phagesdb_ncbi_sequence_mismatch
    def get_phagesdb_ncbi_sequence_length_mismatch(self):
        return self.__phagesdb_ncbi_sequence_length_mismatch

    def get_phamerator_phagesdb_sequence_mismatch(self):
        return self.__phamerator_phagesdb_sequence_mismatch
    def get_phamerator_phagesdb_sequence_length_mismatch(self):
        return self.__phamerator_phagesdb_sequence_length_mismatch
    def get_phamerator_phagesdb_host_mismatch(self):
        return self.__phamerator_phagesdb_host_mismatch
    def get_phamerator_phagesdb_accession_mismatch(self):
        return self.__phamerator_phagesdb_accession_mismatch
    def get_phamerator_phagesdb_cluster_subcluster_mismatch(self):
        return self.__phamerator_phagesdb_cluster_subcluster_mismatch
    def get_total_number_genes_with_errors(self):
        return self.__total_number_genes_with_errors
    def get_contains_errors(self):
        return self.__contains_errors
    def get_ph_ncbi_author_error(self):
        return self.__ph_ncbi_author_error
