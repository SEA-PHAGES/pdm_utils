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















    # Below are old unittests before 'update'/'remove' tickets were
    # split to a separate script from 'add'/'replace' tickets.
    #
    #
    # def test_check_remove_ticket_1(self):
    #     """Standard remove ticket."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 0)
    #
    # def test_check_remove_ticket_2(self):
    #     """Primary Phage ID is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.phage_id = "Trixie"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_3(self):
    #     """host_genus is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.host_genus = "Mycobacterium"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_4(self):
    #     """Cluster is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.cluster = "A"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_5(self):
    #     """Subcluster is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.subcluster = "A2"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_6(self):
    #     """Status is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.status = "final"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_7(self):
    #     """Description Field is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.description_field = "Product"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_8(self):
    #     """Accession is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.accession = "ABC123"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_9(self):
    #     """Annotation Author is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.annotation_author = "1"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_10(self):
    #     """Run Mode is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.run_mode = "phagesdb"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_11(self):
    #     """Secondary Phage ID is not present."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.secondary_phage_id = "D29"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    #
    #
    #
    #
    #
    #
    # def test_check_ticket_1(self):
    #     """Check standard update ticket."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.update_ticket.evaluations), 0)
    #
    # def test_check_ticket_2(self):
    #     """Check update ticket with Primary Phage ID not present."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.update_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.update_ticket.evaluations), 1)
    #
    # def test_check_ticket_3(self):
    #     """Check standard add ticket."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 0)
    #
    # def test_check_ticket_4(self):
    #     """Check add ticket with Primary Phage ID already present."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.add_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_ticket_5(self):
    #     """Check standard remove ticket."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 0)
    #
    # def test_check_ticket_6(self):
    #     """Check remove ticket with Primary Phage ID not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.phage_id = "Trixie"
    #     self.remove_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_ticket_7(self):
    #     """Check standard replace ticket."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 0)
    #
    # def test_check_ticket_8(self):
    #     """Check replace ticket when Primary Phage ID is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock","none"])
    #     self.replace_ticket.phage_id = "none"
    #     self.replace_ticket.secondary_phage_id = "none"
    #     self.replace_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_ticket_9(self):
    #     """Check non-standard type of ticket."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.update_ticket.type = "other"
    #     self.update_ticket.check_ticket(phage_id_set)
    #     self.assertEqual(len(self.update_ticket.evaluations), 0)
    #
    ###Above = pasted from test.Ticket.py, since some Ticket methods
    ###were moved to evaluate.py
