"""Retrieve data from Phamerator and parse into Genome objects

"""






def main(sql_obj):

    retrieved_phamerator_data = prepare_sql_data(sql_obj)

    phamerator_genome_dict = parse_phamerator_genomes(retrieved_phamerator_data)

    phamerator_genome_sets = create_data_sets(phamerator_genome_dict)

    return phamerator_genome_dict,phamerator_genome_sets

###
