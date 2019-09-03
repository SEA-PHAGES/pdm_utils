"""Represents a structure to directly compare data between two or more CDS
features.
"""





#Old python2 code
class CdsPair:

    # Initialize all attributes:
    def __init__(self):
        pass



#Old python2 code
class MatchedCdsFeatures:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.__phamerator_feature = ''
        self.__ncbi_feature = ''

        #Matched data comparison results
        self.__phamerator_ncbi_different_translations = False #True = there are different translations
        self.__phamerator_ncbi_different_start_sites = False #True = there are different start sites
        self.__phamerator_ncbi_different_descriptions = False #True = there are different gene descriptions

        #Total errors summary
        self.__total_errors = 0



    # Define all attribute setters:
    def set_phamerator_feature(self,value):
        self.__phamerator_feature = value
    def set_ncbi_feature(self,value):
        self.__ncbi_feature = value


    def compare_phamerator_ncbi_cds_features(self):

        if self.__phamerator_feature.get_strand() == 'forward':
            if str(self.__phamerator_feature.get_left()) != str(self.__ncbi_feature.get_left()):
                self.__phamerator_ncbi_different_start_sites = True
        elif self.__phamerator_feature.get_strand() == 'reverse':
            if str(self.__phamerator_feature.get_right()) != str(self.__ncbi_feature.get_right()):
                self.__phamerator_ncbi_different_start_sites = True
        else:
            pass


        product_set = set()
        product_set.add(self.__phamerator_feature.get_search_notes())
        product_set.add(self.__ncbi_feature.get_search_product())


        if len(product_set) != 1:
            self.__phamerator_ncbi_different_descriptions = True

        if self.__phamerator_feature.get_translation() != self.__ncbi_feature.get_translation():
            self.__phamerator_ncbi_different_translations = True



        #Compute total errors
        #First add all matched feature errors
        if self.__phamerator_ncbi_different_translations:
            self.__total_errors += 1

        if self.__phamerator_ncbi_different_start_sites:
            self.__total_errors += 1

        if self.__phamerator_ncbi_different_descriptions:
            self.__total_errors += 1

        #Now add all errors from each individual feature
        #You first compute errors for each individual feature.
        #This step is performed here instead of in the mainline code
        #because you need to wait for the feature matching step after the genome matching step
        self.__phamerator_feature.compute_total_cds_errors()
        self.__ncbi_feature.compute_total_cds_errors()
        self.__total_errors += self.__phamerator_feature.get_total_errors()
        self.__total_errors += self.__ncbi_feature.get_total_errors()





    # Define all attribute getters:
    def get_phamerator_feature(self):
        return self.__phamerator_feature
    def get_ncbi_feature(self):
        return self.__ncbi_feature
    def get_phamerator_ncbi_different_start_sites(self):
        return self.__phamerator_ncbi_different_start_sites
    def get_phamerator_ncbi_different_descriptions(self):
        return self.__phamerator_ncbi_different_descriptions
    def get_phamerator_ncbi_different_translations(self):
        return self.__phamerator_ncbi_different_translations
    def get_total_errors(self):
        return self.__total_errors
