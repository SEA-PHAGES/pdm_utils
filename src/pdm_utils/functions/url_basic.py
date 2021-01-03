import re
from collections import deque
import sys
import tempfile

import urllib3
from urllib3.exceptions import (MaxRetryError, PoolError, RequestError)
import yaml
from yaml.scanner import ScannerError

# GLOBAL FUNCTIONS
# -----------------------------------------------------------------------------
POOL_MAX_CONN_ATTEMPTS_MSG = (
                "Maximum connection pool construction attempts reached.\n  "
                "Please check your internet connection.")
RESPONSE_MAX_CONN_ATTEMPTS_MSG = (
                "Maximum request attempts reached.\n\t"
                "Please check the availability at the requested url.")

HREF_NAME_FORMAT = ("""<a href="([-_\w\d./]+)">""")
HREF_DATE_FORMAT = ("""<td align="right">"""
                    """(\d+[-]\d+[-]\d+)\s+\d+[:]\d+\s+</td>""")


# GET FUNCTIONS
# -----------------------------------------------------------------------------
def create_pool(pipeline=False):
    try:
        pool = urllib3.PoolManager()
    except MaxRetryError:
        if pipeline:
            print(POOL_MAX_CONN_ATTEMPTS_MSG)
            sys.exit(1)

        raise PoolError(POOL_MAX_CONN_ATTEMPTS_MSG)

    return pool


def pool_request(url, pool=None,
                 pipeline=False, preload=False, expect_status=200):
    """
    Create a urllib3 PoolManager to access the url link for list of databases
    :returns An HTTPResponse object
    :rtype: urllib3.response.HTTPResponse
    """
    if pool is None:
        pool = create_pool()
    try:
        response = pool.request("GET", url, preload_content=preload)
    except MaxRetryError:
        url_fail_msg = "\n".join([RESPONSE_MAX_CONN_ATTEMPTS_MSG,
                                  "".join(["\t\tRequest failed for ", url])])
        if pipeline:
            print(url_fail_msg)
            sys.exit(1)

        raise RequestError(pool, url, url_fail_msg)

    pool.clear()

    if pipeline:
        if response.status != expect_status:
            print("Received invalid response from server.\n"
                  "Aborting pipeline.")
            sys.exit(1)

    return response


def download_file(file_url, filepath, pool=None, chunk_size=512000,
                  preload=False, verbose=True, expect_status=200):
    """
    Retrieve a file from the server.
    :param file_url:  URL for the desired file:
    :type file_url: str
    :param filepath: Local path where the file can be downloaded to.
    :type filepath: Path
    :returns: If status of the file retrieved from the server is expected.
    :rtype: bool
    """
    request = pool_request(file_url, pool=pool, pipeline=False)
    status = request.status

    if status == expect_status:
        filesize = int(request.getheader("Content-Length"))
        output_file = filepath.open(mode="wb")
        for chunk in request.stream(chunk_size):
            output_file.write(chunk)
            progress = int(request._fp_bytes_read / filesize * 100)
            if verbose:
                print("\r[{}{}] {}%".format("#" * int(progress / 2),
                                            " " * (50 - int(progress / 2)),
                                            progress), end="")
    else:
        if verbose:
            print(" ".join(["ERROR.  HTTP Response", str(status)]))

    if verbose:
        print("")

    return status == expect_status


def spool_file(file_url, pool=None, max_size=512000, chunk_size=512000,
               expect_status=200):
    """
    Retrieve a file from the server.
    :param file_url:  URL for the desired file:
    :type file_url: str
    :param filepath: Local path where the file can be downloaded to.
    :type filepath: Path
    :param expect_status: Status expected from the HTTPresponse object
    :type expect_status: int
    :returns: Returns temporary file handle with the downloaded data
    :rtype: tempfile.SpooledTemporaryFile
    """
    request = pool_request(file_url, pool=pool, pipeline=False)
    status = request.status

    if status == expect_status:
        spooled_file = tempfile.SpooledTemporaryFile(max_size=max_size)
        for chunk in request.stream(chunk_size):
            spooled_file.write(chunk)

    return spooled_file


def get_yaml_data(response):
    """Parses yaml data from a response using the PyYAML library.

    :param response: An HTTPrepsonse object:
    :type response: urllib3.response.HTTPResponse
    :return: Dictionary containing key-value pairs from the YAML file
    :rtype: dict
    """
    try:
        data = yaml.safe_load(response.data)
    except ScannerError:
        data = None

    return data


def get_url_listing_data(response, get_dates=False):
    """
    Get list of file and directory data from the link using a response object
    :param response: PoolManager response for the specified url link
    :type response: urllib3.response.HTTPResponse
    :returns: A list of entry listing data
    :rtype: list[dict]
    """
    listings = list()

    response_data = response.data.decode("utf-8")
    response_data_lines = response_data.split("\n")

    # Get database names
    names_regex = re.compile(HREF_NAME_FORMAT)
    dates_regex = re.compile(HREF_DATE_FORMAT)

    data = []
    for line in response_data_lines:
        name_match = re.match(names_regex, line)
        date_match = re.match(dates_regex, line)

        date = None
        if date_match is not None:
            # Retrieves the first capturing group match to the date regex
            date = date_match[1]

        name = None
        if name_match is not None:
            # Retrieves the first capturing group match to the name regex
            name = name_match[1]

        if name is not None:
            data.append((name, date))

    response.close()

    # creating a list of dictionaries
    for i in range(len(data)):
        db_dict = dict()
        db_dict["num"] = i + 1
        db_dict["name"] = data[i][0]
        db_dict["date"] = data[i][1]

        listings.append(db_dict)

    return listings


# PARSING FUNCTIONIS
# ----------------------------------------------------------------------------
def get_url_listing_names(response, suffix=None, ignore=None):
    """Get a list of file and directory names at the url using a HTTPresponse

    :param response: PoolManager response for the specified url link
    :type response: urllib3.response.HTTPResponse
    :param suffix: Suffix string to filter names for
    :type suffix: str
    :param ignore: Prefix string to filter names against
    :type ignore: str
    :returns: A list of entry listing names
    :rtype: list
    """
    listing = get_url_listing_data(response)

    listing_names = list()

    for data_dict in listing:
        name = data_dict.get("name")

        if name is None:
            continue

        if ignore:
            if name.startswith(ignore):
                continue

        if suffix:
            if len(suffix) >= len(name):
                continue

            if name.endswith(suffix):
                # Removes suffix before appending to names list
                listing_names.append(name[:-(len(suffix))])
        else:
            listing_names.append(name)

    return listing_names


def get_url_listing_dirs(response, ignore="."):
    """Get a list of file and directory names at the url using a HTTPresponse

    :param response: PoolManager response for the specified url link
    :type response: urllib3.response.HTTPResponse
    :param ignore: Prefix string to filter names against
    :type ignore: str
    :returns: A list of entry listing names
    :rtype: list
    """
    return get_url_listing_names(response, suffix="/", ignore=ignore)


def get_url_listing_files(response, file_ext, ignore=None):
    """Get a list of file and directory names at the url using a HTTPresponse

    :param response: PoolManager response for the specified url link
    :type response: urllib3.response.HTTPResponse
    :param file_ext: File extension to filter names for
    :type file_ext: str
    :param ignore: Prefix string to filter names against
    :type ignore: str
    :returns: A list of entry listing names
    :rtype: list
    """
    return get_url_listing_names(response, "".join([".", str(file_ext)]),
                                 ignore=ignore)


class UrlPathGraph():
    """Graph object that manages and stores PathNodes representing
    directories at an url"""
    def __init__(self, url, pool=None, pipeline=False, preload=False,
                 expect_status=200):
        self.nodes = list()

        if pool is None:
            pool = create_pool(pipeline=pipeline)

        self.pool = pool
        self.expect_status = expect_status
        self.pipeline = pipeline

        self.url = url
        self.root = PathNode("", root_path=url, isdir=True)
        self.curr = self.root
        self.nodes.append(self.root)
        self.load_node(self.root)

        if preload:
            self.load_graph()

    def load_graph(self):
        """Loads the structure of a filesystem at the specified url as a graph
        """
        node_stack = deque()

        # Do/while loop that iterates until all directories are explored
        parent = self.root
        while True:
            if not parent.loaded:
                self.load_node(parent)

            for child_name, child_node in parent.children.items():
                # Ignores hidden nodes
                if not child_name.startswith("."):
                    node_stack.append(child_node)

            try:
                parent = node_stack.pop()
            except IndexError:
                break

    def load_node(self, node):
        """Uses a HTTPresponse object to gather data to load a PathNode

        :param node: PathNode to load data and children for
        :type node: PathNode
        """
        response = pool_request(node.get_abs_path(), pool=self.pool,
                                expect_status=self.expect_status,
                                pipeline=self.pipeline)

        node.contents = get_url_listing_names(response)
        for child_dir in get_url_listing_dirs(response):
            child_node = PathNode(child_dir, root_path=self.url,
                                  parent=node)

            self.nodes.append(child_node)

            node.children[child_dir] = child_node

        if len(node.children) > 1:
            node.isdir = True
            if node.parent is not None:
                node.children[".."] = node.parent

        response.close()
        node.loaded = True

    def traverse(self, path):
        """ Traverses a linked list using a filesystem path

        :param path: A filesystem path (node names separated by '/')
        :type path: str
        :return: Returns the path node specified by the filesystem path
        :rtype: PathNode
        """
        curr = self.traverse_relative(path, self.curr)

        if curr is None:
            curr = self.traverse_absolute(path)

        return curr

    def traverse_relative(self, path, node):
        """Traverses a linked list using a filesystem path relative to a node

        :param path: A filesystem path (node names separated by '/')
        :type path: str
        :param node: A PathNode designating the start for the path
        ;type node: PathNode
        :return: Returns the path node specified by the filesystem path
        :rtype: PathNode
        """
        node_stops = path.split("/")

        curr = node
        for node_name in node_stops:
            if not curr.loaded:
                self.load_node(curr)

            curr = curr.children.get(node_name)

            if curr is None:
                return None

        return curr

    def traverse_absolute(self, path):
        """Traverses a linked list using a filesystem path relative to root

        :param path: A filesystem path (node names separated by '/')
        :type path: str
        :return: Returns the path node specified by the filesystem path
        :rtype: PathNode
        """
        node_stops = path.split("/")

        if node_stops[0] == "":
            node_stops.pop(0)

        path = "/".join(node_stops)
        return self.traverse_relative(path, self.root)


class PathNode():
    """Object that stores information about a directory in a file system"""
    def __init__(self, name, root_path="", rel_path=None, abs_path=None,
                 isdir=False, children=None, parent=None, data=None):
        self.name = name
        self.isdir = isdir

        self.root_path = root_path
        self.rel_path = rel_path
        self.abs_path = abs_path

        if children is None:
            children = dict()
        children["."] = self

        self.children = children
        self.parent = parent

        self.data = data
        self.contents = []
        self.loaded = False

    def calc_rel_path(self, end_node=None):
        """Calculates the filesystem path relative to the node reference

        :param end_node: Node to create the relative path from
        :type end_node: PathNode
        :return: Returns the path relative to the node reference
        :rtype: str
        """
        path_nodes = list()
        curr_node = self
        while True:
            path_nodes.insert(0, curr_node.name)

            curr_node = curr_node.parent
            if curr_node == end_node or curr_node is None:
                break

        return "/".join(path_nodes)

    def set_rel_path(self, end_node=None):
        """Sets the filesystem path relative to the node reference

        :param end_node: Node to create the relative path from
        :type end_node: PathNode
        """
        self.rel_path = self.calc_rel_path(end_node=end_node)

    def set_abs_path(self):
        """Sets the filesystem path relative to the root directory
        """
        self.abs_path = "/".join([self.root_path.rstrip("/"),
                                  self.calc_rel_path().lstrip("/")])

    def get_rel_path(self, node=None):
        """Retrieves the filesystem path relative to the deepest node reference

        :return: Returns the relative path
        :rtype: str
        """
        if self.rel_path is None:
            self.set_rel_path()

        return self.rel_path

    def get_abs_path(self):
        """Retrieves the filesystem path relative to the root directory

        :return: Returns the absolute path
        :rtype: str
        """
        if self.abs_path is None:
            self.set_abs_path()

        return self.abs_path

    def print_children(self, max_width=74):
        """Prints the filesystem directory child directories
        """
        names = list()
        max_len = 0
        for child_name, child_node in self.children.items():
            if child_name == ".":
                continue

            if child_node.isdir:
                child_name = "".join([child_name, "/"])

            names.append(child_name)
            if len(child_name) > max_len:
                max_len = len(child_name)

        for i in range(len(names)):
            child = names[i]
            child = "".join([child, " " * (max_len - len(child)), "\t"])
            names[i] = child

        cols = int(max_width / (max_len + (4 - max_len % 4)))

        for row in range(0, len(names), cols):
            child_row = names[row:row+cols]
            print("".join(["\t"] + child_row))

    def print_contents(self, max_width=74):
        """Prints the filesystem directory child contents
        """
        names = list()
        max_len = 0
        for name in self.contents:
            names.append(name)

            if len(name) > max_len:
                max_len = len(name)

        for i in range(len(names)):
            name = names[i]
            name = "".join([name, " " * (max_len - len(name)), "\t"])
            names[i] = name

        cols = int(max_width / (max_len + (4 - max_len % 4)))

        for row in range(0, len(names), cols):
            name_row = names[row:row+cols]
            print("".join(["\t"] + name_row))
