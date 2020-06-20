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


    @patch("pdm_utils.functions.annotation.get_adjacent_genes")
    @patch("pdm_utils.functions.annotation.querying.first_column")
    @patch("pdm_utils.functions.annotation.select")
    def test_get_adjacent_phams_1(self, mock_select, mock_first_column,
                                                     mock_get_adjacent_genes):
        """Verify select() call structure of get_adjacent_phams()"""
        mock_select_obj = Mock()
        mock_select_where_obj = Mock()
        mock_select_distinct_obj = Mock()

        mock_select.return_value = mock_select_obj
        mock_select_obj.where.return_value = mock_select_where_obj
        mock_select_where_obj.distinct.return_value = mock_select_distinct_obj

        mock_first_column.return_value = ["Trixie_CDS_2"]
        mock_get_adjacent_genes.return_value = ("Trixie_CDS_1", "Trixie_CDS_3")

        annotation.get_adjacent_phams(self.mock_alchemist, 8675309)

        mock_select.assert_any_call([self.geneid_column])
        mock_select.assert_any_call([self.phamid_column])
        mock_select_obj.where.assert_called()
        mock_select_where_obj.distinct.assert_called()

    @patch("pdm_utils.functions.annotation.get_adjacent_genes")
    @patch("pdm_utils.functions.annotation.querying.first_column")
    @patch("pdm_utils.functions.annotation.select")
    def test_get_adjacent_phams_2(self, mock_select, mock_first_column,
                                                     mock_get_adjacent_genes):
        """Verify first() call structure of get_adjacent_phams()"""
        mock_select_obj = Mock()
        mock_select_where_obj = Mock()
        mock_select_distinct_obj = Mock()

        mock_select.return_value = mock_select_obj
        mock_select_obj.where.return_value = mock_select_where_obj
        mock_select_obj.distinct.return_value = mock_select_distinct_obj
        mock_select_where_obj.distinct.return_value = mock_select_distinct_obj

        mock_first_column.return_value = ["Trixie_CDS_2"]
        mock_get_adjacent_genes.return_value = ("Trixie_CDS_1", "Trixie_CDS_3")

        annotation.get_adjacent_phams(self.mock_alchemist, 8675309)

        mock_first_column.assert_any_call(self.mock_engine, 
                                          mock_select_distinct_obj)
        mock_first_column.assert_any_call(self.mock_engine, 
                                          mock_select_distinct_obj,
                                          in_column=self.geneid_column,
                                          values=["Trixie_CDS_1"])
        mock_first_column.assert_any_call(self.mock_engine, 
                                          mock_select_distinct_obj,
                                          in_column=self.geneid_column,
                                          values=["Trixie_CDS_3"])

    @patch("pdm_utils.functions.annotation.get_adjacent_genes")
    @patch("pdm_utils.functions.annotation.querying.first_column")
    @patch("pdm_utils.functions.annotation.select")  
    def test_get_adjacent_phams_3(self, mock_select, mock_first_column,
                                                     mock_get_adjacent_genes): 
        """Verify get_adjacent_phams() calls get_adjacent_genes()"""
        mock_first_column.return_value = ["Trixie_CDS_2"]
        mock_get_adjacent_genes.return_value = ("Trixie_CDS_1", "Trixie_CDS_3")
         
        annotation.get_adjacent_phams(self.mock_alchemist, 8675309)

        mock_get_adjacent_genes.assert_called_with(self.mock_alchemist, 
                                                                "Trixie_CDS_2")

    @patch("pdm_utils.functions.annotation.func.count")
    @patch("pdm_utils.functions.annotation.querying.first_column")
    @patch("pdm_utils.functions.annotation.select")
    def test_get_count_pham_annotations_1(self, mock_select, mock_first_column,
                                                             mock_count):
        """Verify select() call structure of get_count_pham_annotations()"""
        mock_count_geneid = Mock()
        mock_count.return_value = mock_count_geneid

        mock_first_column.return_value = ["acr".encode("utf-8")]

        annotation.get_count_pham_annotations(self.mock_alchemist, 8675309)

        mock_count.assert_called_with(self.geneid_column)

        mock_select.assert_any_call([self.notes_column])
        mock_select.assert_any_call([mock_count_geneid])

    @patch("pdm_utils.functions.annotation.func.count")
    @patch("pdm_utils.functions.annotation.querying.first_column")
    @patch("pdm_utils.functions.annotation.select")   
    def test_get_count_pham_annotations_2(self, mock_select, mock_first_column,
                                                             mock_count):
        """Verify first() call structure of get_count_pham_annotations()"""
        mock_select_obj = Mock()
        mock_select_where_obj = Mock()
        mock_select_distinct_obj = Mock()

        mock_select.return_value = mock_select_obj
        mock_select_obj.where.return_value = mock_select_where_obj
        mock_select_obj.distinct.return_value = mock_select_distinct_obj
        mock_select_where_obj.distinct.return_value = mock_select_distinct_obj

        mock_count_geneid = Mock()
        mock_count.return_value = mock_count_geneid

        mock_first_column.return_value = ["acr".encode("utf-8")]

        annotation.get_count_pham_annotations(self.mock_alchemist, 8675309)

        mock_count.assert_called_with(self.geneid_column)

        mock_first_column.assert_called_with(self.mock_engine, 
                                             mock_select_distinct_obj)  

if __name__ == "__main__":
    unittest.main()
