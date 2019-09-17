"""Temp notes about update tickets, and potential checks before
implementing update or remove tickets.
"""


# TODO complete function. It should create SQL statements to update data.
# TODO implement.
# TODO unit test.
def check_update_tickets():
    """."""
    pass
    # index = 0
    # while index < len(list_of_update_objects):
    #
    #     bndl = list_of_update_objects[index]
    #
    #     if len(bndl.genome_pairs_dict.keys()) == 0:
    #         # TODO throw an error if there is no matched Phamerator genome?
    #         pass
    #
    #     for key in bndl.genome_pairs_dict.keys():
    #
    #         genome_pair = bndl.genome_pairs_dict[key]
    #
    #         # TODO check for conflicting hosts. It is not common to
    #         # change hosts.
    #         genome_pair.check_xyz()
    #
    #         # TODO check for conflicting status. It is not common to
    #         # change status unless the current phamerator is draft status.
    #         # The only common change is from draft to final.
    #         genome_pair.check_xyz()
    #
    #         # TODO check for conflicting accession.
    #         # It is not common to change from real accession to another
    #         # real accession. But it is common to change from 'none' to
    #         # real accession.
    #         genome_pair.check_xyz()
    #
    #         # TODO check for conflicint authorship.
    #         # It is not common to change authorships.
    #         genome_pair.check_xyz()
    #
    #     index += 1





# TODO complete function. It should create SQL statements to remove data.
# TODO implement.
# TODO unit test.
def check_remove_tickets(list_of_remove_objects, genome_type):
    """."""
    pass
    #
    #
    # index = 0
    # while index < len(list_of_remove_objects):
    #
    #     bndl = list_of_update_objects[index]
    #
    #     try:
    #         gnm = bndl.genomes_dict[genome_type]
    #     except:
    #         # TODO throw an error if there is no matched Phamerator genome?
    #         continue
    #
    #
    #     # TODO list of evaluations for remove ticket.
    #
    #     # TODO check the status of the removing genome.
    #     # It is not common to remove anything but a 'draft' genome.
    #
    #
    #     index += 1
