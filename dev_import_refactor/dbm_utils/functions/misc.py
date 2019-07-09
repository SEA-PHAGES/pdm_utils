"""Misc. functions useful for comparing and processing genomes."""







def match_genomes(list_of_group_objects, genome_dict, key1, key2 = None):
    """Match genome object to another genome object using phage_id.
    The 'key1' parameter provides the type of genome stored in the DataGroup
    genome dictionary to base the match from.
    The 'key2' parameter provides the type of genome to be stored in
    the DataGroup genome dictionary."""

    index = 0
    while index < len(list_of_group_objects):

        group_obj = list_of_group_objects[index]

        try:
            ref_genome = group_obj.genome_dict[key1]
            if ref_genome.phage_id in genome_dict.keys():
                matched_genome = genome_dict[ref_genome.phage_id]

                if key2 is None:
                    group_obj.genome_dict[matched_genome.type] = matched_genome
                else:
                    group_obj.genome_dict[key2] = matched_genome

        except:
            pass

        index += 1












###
