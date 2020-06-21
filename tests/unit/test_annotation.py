import unittest
from unittest.mock import Mock, patch, PropertyMock

from sqlalchemy import Column

from pdm_utils.functions import annotation

class TestAnnotationRetrieval(unittest.TestCase):
    def setUp(self):
        self.mock_alchemist = Mock()

        self.mock_engine = Mock()
        self.mock_execute_obj = Mock()

        self.mock_metadata = Mock()

        self.mock_filter = Mock() 


        self.gene_obj_mock = Mock()
        self.gene_obj_columns = Mock()
        self.tables_dict = {"gene" : self.gene_obj_mock}
 
        self.geneid_column = Column("GeneID")
        self.phageid_column = Column("PhageID")
        self.phamid_column = Column("PhamID")
        self.name_column = Column("Name")
        self.notes_column = Column("Notes")


        type(self.mock_alchemist).engine = PropertyMock(
                                            return_value=self.mock_engine)
        type(self.mock_alchemist).metadata = PropertyMock(
                                            return_value=self.mock_metadata)
        type(self.mock_metadata).tables = PropertyMock(
                                            return_value=self.tables_dict)

        self.mock_engine.execute.return_value = self.mock_execute_obj


        type(self.gene_obj_mock).c = PropertyMock(
                                            return_value=self.gene_obj_columns)

        type(self.gene_obj_columns).GeneID = PropertyMock(
                                            return_value=self.geneid_column)
        type(self.gene_obj_columns).PhageID = PropertyMock(
                                            return_value=self.phageid_column)
        type(self.gene_obj_columns).PhamID = PropertyMock(
                                            return_value=self.phamid_column)
        type(self.gene_obj_columns).Notes = PropertyMock(
                                            return_value=self.notes_column)
        type(self.gene_obj_columns).Name = PropertyMock(
                                            return_value=self.name_column)

    @patch("pdm_utils.functions.annotation.select")
    def test_get_relative_gene_1(self, mock_select):
        """Verify select() function calls of get_relative_gene()."""
        annotation.get_relative_gene(self.mock_alchemist, "Trixie_CDS_2", -1)

        mock_select.assert_any_call([self.geneid_column])
   
    @patch("pdm_utils.functions.annotation.select")
    def test_get_relative_gene_2(self, mock_select):
        """Verify get_relative_gene() raises ValueError at bad GeneID input."""
        with self.assertRaises(ValueError):
            annotation.get_relative_gene(self.mock_alchemist, "BAD_GENE", -1)

    @patch("pdm_utils.functions.annotation.get_relative_gene")
    def test_get_adjacent_genes_1(self, mock_get_relative_gene):
        """Verify function structure of get_adjacent_genes()"""
        annotation.get_adjacent_genes(self.mock_alchemist, "Trixie_CDS_2")

        mock_get_relative_gene.assert_any_call(
                                   self.mock_alchemist, "Trixie_CDS_2", 1)
        mock_get_relative_gene.assert_any_call(
                                   self.mock_alchemist, "Trixie_CDS_2", -1)


if __name__ == "__main__":
    unittest.main()
