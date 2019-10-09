"""Integration tests for misc. functions that interact with webfaction server."""
import unittest
from pdm_utils.constants import constants
from pdm_utils.functions import server


class TestPhameratorFunctions1(unittest.TestCase):

    def setUp(self):
        self.valid_host = constants.DB_HOST
        self.invalid_host = "invalid"


    def test_get_transport_1(self):
        """Verify transport is not setup using a invalid host."""
        transport = server.get_transport(self.invalid_host)
        self.assertIsNone(transport)

    def test_get_transport_2(self):
        """Verify transport is setup using a valid host."""
        transport = server.get_transport(self.valid_host)
        print("before")
        self.assertIsNotNone(transport)
        print("done")

    # def tearDown(self):
    #     """Make sure all connections are closed."""
    #     if transport is not None:
    #         transport.close()




###
