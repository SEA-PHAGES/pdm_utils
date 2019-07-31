"""Integration tests for the main import pipeline."""



import unittest
from pipelines.database_import import import_main

class TestImportMain(unittest.TestCase):


    def setUp(self):
        pass



    def test_main1_1(self):
        """."""
        ticket_list1 = ["add",
                        "Trixie",
                        "Mycobacterium",
                        "A",
                        "A2",
                        "Final",
                        "Hatfull",
                        "Product",
                        "ABC123",
                        1,
                        1,
                        "phagesdb"]

        ticket_list2 = ["add",
                        "L5",
                        "Mycobacterium",
                        "A",
                        "A2",
                        "Final",
                        "Hatfull",
                        "Product",
                        "ABC123",
                        1,
                        1,
                        "phagesdb"]

        lists_of_tickets = [ticket_list1, ticket_list2]

        list_of_filenames = []
        import_main.main1(lists_of_tickets, list_of_filenames)




if __name__ == '__main__':
    unittest.main()
