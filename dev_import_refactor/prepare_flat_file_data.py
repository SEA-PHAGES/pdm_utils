








def main(path_to_folder):




    #Validate folder name.
    #TODO complete using previously defined function

    #TODO Validate folder.
    path_to_folder = validate_folder(path_to_folder)

    all_files =  [X for X in os.listdir(path_to_folder) if os.path.isfile(os.path.join(path_to_folder,X))]
    #TODO refactor?
    # write_out(output_file,"\nA total of %s file(s) present in the directory." % len(all_files))


    valid_files,failed_files = create_flat_file_list(all_files)

    flat_file_genomes = parse_all_flat_files(valid_files)

    return flat_file_genomes

###
