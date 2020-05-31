"""Integration tests for the push_db pipeline."""

from pathlib import Path
import shutil
import unittest
from unittest.mock import patch
from unittest.mock import Mock

from paramiko.sftp_client import SFTPClient
from paramiko.transport import Transport

from pdm_utils import run
from pdm_utils.pipelines import push_db

PIPELINE = "push"

# Create the main test directory in which all files will be
# created and managed. Gets created once for all tests.
test_root_dir = Path("/tmp", "pdm_utils_tests_push_db")
if test_root_dir.exists() == True:
    shutil.rmtree(test_root_dir)
test_root_dir.mkdir()

# Within main test folder, test_run_dir folder will be created/removed
# for each test.
test_run_dir = Path(test_root_dir, "test")
log_file = "log_file.txt"
log_path = Path(test_run_dir, log_file)
upload_dir = Path(test_run_dir, "upload_dir")
upload_file1 = Path(upload_dir, "upload_file1.txt")
upload_file2 = Path(upload_dir, "upload_file2.txt")
upload_file3 = Path(upload_dir, "upload_file3.txt")

def get_unparsed_args(file=None, dir=None, server=None, remote_dir=None):
    """Returns list of command line arguments to convert database."""
    unparsed_args = ["run.py", PIPELINE, "-l", str(log_path)]
    if file is not None:
        unparsed_args.extend(["-f", str(file)])
    if dir is not None:
        unparsed_args.extend(["-d", str(dir)])
    if server is not None:
        unparsed_args.extend(["-s", server])
    if remote_dir is not None:
        unparsed_args.extend(["-rd", remote_dir])
    return unparsed_args


class TestPushDb(unittest.TestCase):

    def setUp(self):
        test_run_dir.mkdir()
        upload_dir.mkdir()
        upload_file1.touch()
        upload_file2.touch()
        self.sftpc_mock = Mock(spec=SFTPClient)
        self.transport_mock = Mock(spec=Transport)

    def tearDown(self):
        shutil.rmtree(test_run_dir)




    @patch("pdm_utils.functions.server.upload_file")
    @patch("pdm_utils.functions.server.setup_sftp_conn")
    @patch("pdm_utils.functions.server.get_transport")
    def test_main_1(self, gt_mock, ssc_mock, uf_mock):
        """Confirm pipeline runs to upload a file and all
        parameters are added."""
        # Note don't actually try to upload data, so patch necessary functions.
        gt_mock.return_value = self.transport_mock
        ssc_mock.return_value = self.sftpc_mock
        uf_mock.return_value = True
        unparsed_args = get_unparsed_args(file=upload_file1, server="invalid",
                                          remote_dir="invalid")
        run.main(unparsed_args)
        with self.subTest():
            gt_mock.assert_called()
        with self.subTest():
            ssc_mock.assert_called()
        with self.subTest():
            uf_mock.assert_called()

    @patch("pdm_utils.functions.server.upload_file")
    @patch("pdm_utils.functions.server.setup_sftp_conn")
    @patch("pdm_utils.functions.server.get_transport")
    def test_main_2(self, gt_mock, ssc_mock, uf_mock):
        """Confirm pipeline runs to upload a directory of files and all
        parameters are added."""
        # Note don't actually try to upload data, so patch necessary functions.
        gt_mock.return_value = self.transport_mock
        ssc_mock.return_value = self.sftpc_mock
        uf_mock.return_value = True
        unparsed_args = get_unparsed_args(dir=upload_dir, server="invalid",
                                          remote_dir="invalid")
        run.main(unparsed_args)
        with self.subTest():
            gt_mock.assert_called()
        with self.subTest():
            ssc_mock.assert_called()
        with self.subTest():
            uf_mock.assert_called()

    @patch("pdm_utils.functions.server.upload_file")
    @patch("pdm_utils.functions.server.setup_sftp_conn")
    @patch("pdm_utils.functions.server.get_transport")
    @patch("sys.exit")
    def test_main_3(self, exit_mock, gt_mock, ssc_mock, uf_mock):
        """Confirm sys.exit is called when server_host is not provided."""
        # Note don't actually try to upload data, so patch necessary functions.
        gt_mock.return_value = self.transport_mock
        ssc_mock.return_value = self.sftpc_mock
        uf_mock.return_value = True
        unparsed_args = get_unparsed_args(file=upload_file1, server=None,
                                          remote_dir="invalid")
        run.main(unparsed_args)
        exit_mock.assert_called()

    @patch("pdm_utils.pipelines.push_db.upload")
    @patch("pdm_utils.functions.server.setup_sftp_conn")
    @patch("pdm_utils.functions.server.get_transport")
    @patch("sys.exit")
    def test_main_4(self, exit_mock, gt_mock, ssc_mock, u_mock):
        """Confirm sys.exit is called when remote_dir is not provided."""
        # Note don't actually try to upload data, so patch necessary functions.
        gt_mock.return_value = self.transport_mock
        ssc_mock.return_value = self.sftpc_mock
        u_mock.return_value = ([], [])
        unparsed_args = get_unparsed_args(file=upload_file1, server="invalid",
                                          remote_dir=None)
        run.main(unparsed_args)
        exit_mock.assert_called()

    @patch("pdm_utils.pipelines.push_db.upload")
    @patch("pdm_utils.functions.server.setup_sftp_conn")
    @patch("pdm_utils.functions.server.get_transport")
    @patch("sys.exit")
    def test_main_5(self, exit_mock, gt_mock, ssc_mock, u_mock):
        """Confirm sys.exit is called when non-existent file is provided."""
        # Note don't actually try to upload data, so patch necessary functions.
        gt_mock.return_value = self.transport_mock
        ssc_mock.return_value = self.sftpc_mock
        u_mock.return_value = ([], [])
        unparsed_args = get_unparsed_args(file=upload_file3, server="invalid",
                                          remote_dir="invalid")
        run.main(unparsed_args)
        exit_mock.assert_called()

    @patch("pdm_utils.functions.server.upload_file")
    @patch("pdm_utils.functions.server.setup_sftp_conn")
    @patch("pdm_utils.functions.server.get_transport")
    def test_main_6(self, gt_mock, ssc_mock, uf_mock):
        """Confirm pipeline does not run completely when no file or directory
        is provided."""
        # Note don't actually try to upload data, so patch necessary functions.
        gt_mock.return_value = self.transport_mock
        ssc_mock.return_value = self.sftpc_mock
        uf_mock.return_value = True
        unparsed_args = get_unparsed_args(server="invalid", remote_dir="invalid")
        run.main(unparsed_args)
        with self.subTest():
            gt_mock.assert_not_called()
        with self.subTest():
            ssc_mock.assert_not_called()
        with self.subTest():
            uf_mock.assert_not_called()

    @patch("pdm_utils.functions.server.upload_file")
    @patch("pdm_utils.functions.server.setup_sftp_conn")
    def test_main_7(self, ssc_mock, uf_mock):
        """Confirm pipeline does not run completely when unable to create
        transport."""
        # Note don't actually try to upload data, so patch necessary functions.
        ssc_mock.return_value = self.sftpc_mock
        uf_mock.return_value = True
        unparsed_args = get_unparsed_args(server="invalid", remote_dir="invalid")
        run.main(unparsed_args)
        with self.subTest():
            ssc_mock.assert_not_called()
        with self.subTest():
            uf_mock.assert_not_called()

    @patch("pdm_utils.functions.server.upload_file")
    @patch("pdm_utils.functions.server.get_transport")
    def test_main_8(self, gt_mock, uf_mock):
        """Confirm pipeline does not run completely when unable to create
        sftp connection."""
        # Note don't actually try to upload data, so patch necessary functions.
        gt_mock.return_value = self.transport_mock
        uf_mock.return_value = True
        unparsed_args = get_unparsed_args(server="invalid", remote_dir="invalid")
        run.main(unparsed_args)
        with self.subTest():
            gt_mock.assert_not_called()
        with self.subTest():
            uf_mock.assert_not_called()

    @patch("pdm_utils.pipelines.push_db.upload")
    @patch("pdm_utils.functions.server.setup_sftp_conn")
    @patch("pdm_utils.functions.server.get_transport")
    def test_main_9(self, gt_mock, ssc_mock, u_mock):
        """Confirm failed files are printed."""
        # Note don't actually try to upload data, so patch necessary functions.
        gt_mock.return_value = self.transport_mock
        ssc_mock.return_value = self.sftpc_mock
        u_mock.return_value = (["file1.txt", "file2.txt"],
                               ["file3.txt", "file4.txt"])
        unparsed_args = get_unparsed_args(file=upload_file1, server="invalid",
                                          remote_dir="invalid")
        run.main(unparsed_args)
        with self.subTest():
            gt_mock.assert_called()
        with self.subTest():
            ssc_mock.assert_called()
        with self.subTest():
            u_mock.assert_called()

if __name__ == '__main__':
    unittest.main()
