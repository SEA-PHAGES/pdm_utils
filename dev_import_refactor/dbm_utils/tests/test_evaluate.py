""" Unit tests for evaluate functions."""



import unittest



class TestEvaluateClass(unittest.TestCase):


    def setUp(self):
        pass










    ###Below = pasted from test.Ticket.py, since some Ticket methods
    ###were moved to evaluate.py
    #
    # def test_check_update_ticket_1(self):
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[0].status, "correct")
    #
    # def test_check_update_ticket_2(self):
    #     """Primary Phage ID not in set."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[0].status, "error")
    #
    # def test_check_update_ticket_3(self):
    #     """Host is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.host = "none"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[1].status, "error")
    #
    # def test_check_update_ticket_4(self):
    #     """Cluster is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.cluster = "none"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[2].status, "error")
    #
    # def test_check_update_ticket_5(self):
    #     """Status is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.status = "none"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[3].status, "error")
    #
    # def test_check_update_ticket_6(self):
    #     """Description Field is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.description_field = "Product"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[4].status, "error")
    #
    # def test_check_update_ticket_7(self):
    #     """Secondary Phage ID is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.secondary_phage_id = "Trixie"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[5].status, "error")
    #
    # def test_check_update_ticket_8(self):
    #     """Annotation Author is 0."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.annotation_author = "0"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[6].status, "correct")
    #
    # def test_check_update_ticket_9(self):
    #     """Annotation Author is not 1 or 0."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.annotation_author = "none"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[6].status, "error")
    #
    # def test_check_update_ticket_10(self):
    #     """Run Mode is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.update_ticket.run_mode = "phagesdb"
    #     self.update_ticket.check_update_ticket(phage_id_set)
    #     self.assertEqual(self.update_ticket.evaluations[7].status, "error")
    #






    #
    # def test_check_add_ticket_1(self):
    #     """Standard add ticket."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 0)
    #
    # def test_check_add_ticket_2(self):
    #     """Primary Phage ID already present."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_3(self):
    #     """Primary Phage ID is none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.primary_phage_id = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_4(self):
    #     """Host is none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.host = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_5(self):
    #     """Cluster is none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.cluster = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_6(self):
    #     """Status is none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.status = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_7(self):
    #     """Status is Final."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.status = "final"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_8(self):
    #     """Description Field is none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.description_field = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_9(self):
    #     """Secondary Phage ID is not none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.secondary_phage_id = "Trixie"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_10(self):
    #     """Annotation Author is 0."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.annotation_author = "0"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 0)
    #
    # def test_check_add_ticket_11(self):
    #     """Annotation Author is not 1 or 0."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.annotation_author = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    # def test_check_add_ticket_12(self):
    #     """Run Mode is none."""
    #     phage_id_set = set(["L5","RedRock"])
    #     self.add_ticket.run_mode = "none"
    #     self.add_ticket.check_add_ticket(phage_id_set)
    #     self.assertEqual(len(self.add_ticket.evaluations), 1)
    #
    #
    #
    #
    #
    #
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
    #     self.remove_ticket.primary_phage_id = "Trixie"
    #     self.remove_ticket.check_remove_ticket(phage_id_set)
    #     self.assertEqual(len(self.remove_ticket.evaluations), 1)
    #
    # def test_check_remove_ticket_3(self):
    #     """Host is not none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.remove_ticket.host = "Mycobacterium"
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
    #
    #
    # def test_check_replace_ticket_1(self):
    #     """Standard replace ticket."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 0)
    #
    # def test_check_replace_ticket_2(self):
    #     """Primary Phage ID is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock","none"])
    #     self.replace_ticket.primary_phage_id = "none"
    #     self.replace_ticket.secondary_phage_id = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_3(self):
    #     """Primary Phage ID not present and different from Secondary Phage ID."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.primary_phage_id = "D29"
    #     self.replace_ticket.secondary_phage_id = "L5"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_4(self):
    #     """Primary Phage ID is present and different from Secondary Phage ID."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.secondary_phage_id = "L5"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 2)
    #
    # def test_check_replace_ticket_5(self):
    #     """Host is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.host = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_6(self):
    #     """Cluster is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.cluster = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_7(self):
    #     """Status is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.status = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_8(self):
    #     """Description Field is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.description_field = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_9(self):
    #     """Secondary Phage ID is not present."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.primary_phage_id = "D29"
    #     self.replace_ticket.secondary_phage_id = "D29"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_10(self):
    #     """Annotation Author is 0."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.annotation_author = "0"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 0)
    #
    # def test_check_replace_ticket_11(self):
    #     """Annotation Author is not 1 or 0."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.annotation_author = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    # def test_check_replace_ticket_12(self):
    #     """Run Mode is none."""
    #     phage_id_set = set(["Trixie","L5","RedRock"])
    #     self.replace_ticket.run_mode = "none"
    #     self.replace_ticket.check_replace_ticket(phage_id_set)
    #     self.assertEqual(len(self.replace_ticket.evaluations), 1)
    #
    #
    #
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
    #     self.remove_ticket.primary_phage_id = "Trixie"
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
    #     self.replace_ticket.primary_phage_id = "none"
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




if __name__ == '__main__':
    unittest.main()
