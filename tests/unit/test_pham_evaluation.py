import unittest
from decimal import Decimal
from pathlib import Path
from unittest.mock import Mock, patch, PropertyMock

from Bio.Emboss import Applications

from pdm_utils.functions import pham_evaluation


class TestEMBOSSAlign(unittest.TestCase):
    def setUp(self):
        self.mock_emboss_cline = Mock()
        self.mock_emboss_cline.return_value = True, True

        self.mock_query_seq_path = Mock(spec=Path)
        self.mock_target_seq_path = Mock(spec=Path)
        self.mock_outfile_seq_path = Mock(spec=Path)

        self.mock_gapopen = Mock(spec=Decimal)
        self.mock_gapextend = Mock(spec=Decimal)

        self.mock_alignment = Mock()

    @patch("pdm_utils.functions.pham_evaluation.AlignIO.read")
    @patch("pdm_utils.functions.pham_evaluation.create_water_cline")
    @patch("pdm_utils.functions.pham_evaluation.create_stretcher_cline")
    @patch("pdm_utils.functions.pham_evaluation.create_needle_cline")
    def test_pairwise_align_1(self, needle_mock, stretcher_mock, water_mock, 
                                                    alignio_read_mock):
        """Verify pairwise_align() calls EMBOSS-tool methods."""
        needle_mock.return_value = self.mock_emboss_cline
        stretcher_mock.return_value = self.mock_emboss_cline
        water_mock.return_value = self.mock_emboss_cline

        for tool in pham_evaluation.EMBOSS_TOOLS:
            with self.subTest(tool=tool):
                pham_evaluation.pairwise_align(self.mock_query_seq_path,
                                               self.mock_target_seq_path,
                                               self.mock_outfile_seq_path,
                                               tool=tool)

                self.mock_emboss_cline.assert_called()
                alignio_read_mock.assert_called_with(self.mock_outfile_seq_path,
                                                     "emboss")

    @patch("pdm_utils.functions.pham_evaluation.AlignIO.read")
    @patch("pdm_utils.functions.pham_evaluation.create_water_cline")
    @patch("pdm_utils.functions.pham_evaluation.create_stretcher_cline")
    @patch("pdm_utils.functions.pham_evaluation.create_needle_cline")
    def test_pairwise_align_2(self, needle_mock, stretcher_mock, water_mock,
                                                    alignio_read_mock):
        """Verify pairwise_align() calls cline_init with correct parameters."""
        needle_mock.return_value = self.mock_emboss_cline
        stretcher_mock.return_value = self.mock_emboss_cline
        water_mock.return_value = self.mock_emboss_cline

        with self.subTest(tool="needle"):
            pham_evaluation.pairwise_align(self.mock_query_seq_path,
                                           self.mock_target_seq_path,
                                           self.mock_outfile_seq_path,
                                           tool="needle",
                                           gapopen=self.mock_gapopen,
                                           gapextend=self.mock_gapextend)

            needle_mock.assert_called_with(self.mock_query_seq_path,
                                           self.mock_target_seq_path,
                                           self.mock_outfile_seq_path,
                                           self.mock_gapopen, 
                                           self.mock_gapextend)

        with self.subTest(tool="stretcher"):
            pham_evaluation.pairwise_align(self.mock_query_seq_path,
                                           self.mock_target_seq_path,
                                           self.mock_outfile_seq_path,
                                           tool="stretcher",
                                           gapopen=self.mock_gapopen,
                                           gapextend=self.mock_gapextend)

            stretcher_mock.assert_called_with(
                                           self.mock_query_seq_path,
                                           self.mock_target_seq_path,
                                           self.mock_outfile_seq_path,
                                           self.mock_gapopen, 
                                           self.mock_gapextend)
        with self.subTest(tool="water"):
            pham_evaluation.pairwise_align(self.mock_query_seq_path,
                                           self.mock_target_seq_path,
                                           self.mock_outfile_seq_path,
                                           tool="water",
                                           gapopen=self.mock_gapopen,
                                           gapextend=self.mock_gapextend)

            water_mock.assert_called_with( self.mock_query_seq_path,
                                           self.mock_target_seq_path,
                                           self.mock_outfile_seq_path,
                                           self.mock_gapopen, 
                                           self.mock_gapextend)
    
    @patch("pdm_utils.functions.pham_evaluation.AlignIO.read")
    @patch("pdm_utils.functions.pham_evaluation.create_water_cline")
    @patch("pdm_utils.functions.pham_evaluation.create_stretcher_cline")
    @patch("pdm_utils.functions.pham_evaluation.create_needle_cline")
    def test_pairwise_align_3(self, needle_mock, stretcher_mock, water_mock,
                                                    alignio_read_mock):
        """Verify pairwise_align() raises from bad tool input.""" 
        needle_mock.return_value = self.mock_emboss_cline
        stretcher_mock.return_value = self.mock_emboss_cline
        water_mock.return_value = self.mock_emboss_cline

        with self.assertRaises(ValueError):
            pham_evaluation.pairwise_align(self.mock_query_seq_path,
                                           self.mock_target_seq_path,
                                           self.mock_outfile_seq_path,
                                           tool="matcher")

if __name__ == "__main__":
    unittest.main()
