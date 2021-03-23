import unittest
from unittest.mock import Mock
from unittest.mock import patch
from unittest.mock import PropertyMock

from sqlalchemy import MetaData
from sqlalchemy import Table
from sqlalchemy.orm.decl_api import DeclarativeMeta

from pdm_utils.functions import cartography

class TestAutoMapping(unittest.TestCase):
    def setUp(self):
        self.mapper_mock = Mock(spec=DeclarativeMeta)
        self.metadata_mock = Mock(spec=MetaData)

        self.phage_mock = Mock(spec=Table)
        self.gene_mock  = Mock(spec=Table)

        self.classes = {"phage" : self.phage_mock,
                        "gene"  : self.gene_mock}

        type(self.mapper_mock).classes = PropertyMock(
                                                return_value=self.classes)
        type(self.mapper_mock).metadata = PropertyMock(
                                                return_value=self.metadata_mock)

    @patch("pdm_utils.functions.cartography.parsing.translate_table")
    def test_get_map_1(self, translate_table_mock):
        """Verify function structure of get_map().
        """
        translate_table_mock.return_value = "phage"

        map_obj = cartography.get_map(self.mapper_mock, "phage")

        translate_table_mock.assert_called_with(self.metadata_mock, "phage")
        self.assertEqual(map_obj, self.phage_mock)
    

if __name__ == "__main__":
    unittest.main()
