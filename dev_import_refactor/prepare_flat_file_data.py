




#Validate folder name.
#TODO complete using previously defined function




#Identify list of files that are Genbank-formatted.
def create_flat_file_list(all_files):


    #TODO refactor?
    #write_out(output_file,"\n\n\n\nAccessing genbank-formatted files for add/replace actions...")


    list_of_failed_genome_files = []
    list_of_valid_genome_files = []

    for filename in all_files:

        valid = validate_flat_file(filename)

        if valid == 1:
            list_of_valid_genome_files.append(filename)
        else:
            list_of_failed_genome_files.append(filename)

            # TODO error handling - file might not have contained a parsable
        # flat file. Test whether the genome object contains data.
        # If it is empty, add filename to list of failed files.

    return list_of_valid_genome_files,list_of_failed_genome_files


# Function iterates through list of files and returns
# a list of GenBank-formatted flat files and a list of file names
# they could not be parsed.


def parse_all_flat_files():
    list_of_flat_file_genomes = []

    for filename in list_of_valid_genome_files:

        filepath = os.path.join(path_to_folder,filename)
        flat_file_genome = parse_flat_file_data(filepath)
        list_of_flat_file_genomes.append(flat_file_genome)

    return list_of_flat_file_genomes


def main(path_to_folder):

    #TODO Validate folder.
    path_to_folder = validate_folder(path_to_folder)

    all_files =  [X for X in os.listdir(path_to_folder) if os.path.isfile(os.path.join(path_to_folder,X))]
    #TODO refactor?
    # write_out(output_file,"\nA total of %s file(s) present in the directory." % len(all_files))


    valid_files,failed_files = create_flat_file_list(all_files)

    flat_file_genomes = parse_all_flat_files(valid_files)

    return flat_file_genomes

###
