import re
from collections import deque
import sys
import urllib3
from urllib3.exceptions import MaxRetryError
import tempfile

import yaml
from yaml.scanner import ScannerError


MAX_CONN_ATTEMPTS_MSG = ("Maximum connection attempts reached.  "
                         "Please check connection and try again")


# GET FUNCTIONS
# -----------------------------------------------------------------------------
def pool_request(url, pool=None,
                 pipeline=False, preload=False, expect_status=200):
    """
    Create a urllib3 PoolManager to access the url link for list of databases
    :returns An HTTPResponse object
    :rtype: urllib3.response.HTTPResponse
    """
    if pool is None:
        try:
            pool = urllib3.PoolManager()
        except MaxRetryError:
            if pipeline:
                print(MAX_CONN_ATTEMPTS_MSG)
                sys.exit(1)
            else:
                raise MaxRetryError(MAX_CONN_ATTEMPTS_MSG)

    response = pool.request("GET", url, preload_content=preload)

    pool.clear()

    if pipeline:
        if response.status != expect_status:
            print("Received invalid response from server.\n"
                  "Aborting pipeline.")
            sys.exit(1)

    return response


def download_file(file_url, filepath, pool=None, chunk_size=512000,
                  preload=False, verbose=True):
    """
    Retrieve a file from the server.
    :param file_url:  URL for the desired file:
    :type file_url: str
    :param filepath: Local path where the file can be downloaded to.
    :type filepath: Path
    :returns: Status of the file retrieved from the server
    :rtype: bool
    """
    request = pool_request(file_url, pool=pool, pipeline=False)
    status = request.status

    if status == 200:
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

    return status


def spool_file(file_url, pool=None, max_size=512000, chunk_size=512000):
    request = pool_request(file_url, pool=pool, pipeline=False)
    status = request.status

    if status == 200:
        spooled_file = tempfile.SpooledTemporaryFile(max_size=max_size)
        for chunk in request.stream(chunk_size):
            spooled_file.write(chunk)

    return status


def get_yaml_data(response):
    try:
        data = yaml.safe_load(response.data)
    except ScannerError:
        data = None

    return data


def get_url_listing_data(response, get_dates=False):
    """
    Get list of files and directories from the link using a response object
    :param response: PoolManager response for the specified url link
    :type response: urllib3.response.HTTPResponse
    :returns: A list of entry listings
    :rtype: list
    """
    listings = list()

    response_data = response.data.decode("utf-8")
    # print(response_data)

    # Get database names
    names_regex = """<a href="([-_\w\d./]+)">"""
    names_regex = re.compile(names_regex)
    names = names_regex.findall(response_data)

    if get_dates:
        # Get database dates
        dates_regex = ("""<td align="right">"""
                       "(\d+[-]\d+[-]\d+)\s+\d+[:]\d+\s+</td>")
        dates_regex = re.compile(dates_regex)
        dates = dates_regex.findall(response_data)

    response.close()

    # creating a list of dictionaries
    for i in range(len(names)):
        db_dict = dict()
        db_dict["num"] = i+1
        db_dict["name"] = names[i]
        if get_dates:
            db_dict["date"] = dates[i]

        listings.append(db_dict)

    return listings


# PARSING FUNCTIONIS
# ----------------------------------------------------------------------------
def get_url_listing_names(response, suffix=None, ignore=None):
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
                listing_names.append(name[:-(len(suffix))])
        else:
            listing_names.append(name)

    return listing_names


def get_url_listing_dirs(response, ignore="."):
    return get_url_listing_names(response, suffix="/", ignore=ignore)


def get_url_listing_files(response, file_ext):
    return get_url_listing_names(response, "".join([".", str(file_ext)]))


class UrlPathGraph():
    def __init__(self, url, pool=None, pipeline=False, preload=False,
                 expect_status=200):
        self.nodes = list()

        if pool is None:
            pool = urllib3.PoolManager()

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
        node_stack = deque()

        parent = self.root
        while True:
            if not parent.loaded:
                self.load_node(parent)

            for child_name, child_node in parent.children.items():
                if child_name != ".." and child_name != ".":
                    node_stack.append(child_node)

            try:
                parent = node_stack.pop()
            except IndexError:
                break

    def load_node(self, node):
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
        curr = self.traverse_relative(path, self.curr)

        if curr is None:
            curr = self.traverse_absolute(path)

        return curr

    def traverse_relative(self, path, node):
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
        node_stops = path.split("/")

        if node_stops[0] == "":
            node_stops.pop(0)

        path = "/".join(node_stops)
        return self.traverse_relative(path, self.root)


class PathNode():
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

    def set_rel_path(self):
        path_nodes = list()
        curr_node = self
        while True:
            path_nodes.insert(0, curr_node.name)

            curr_node = curr_node.parent
            if curr_node is None:
                break

        self.rel_path = "/".join(path_nodes)

    def set_abs_path(self):
        self.set_rel_path()

        self.abs_path = "/".join([self.root_path.rstrip("/"),
                                  self.rel_path.lstrip("/")])

    def get_rel_path(self):
        if self.rel_path is None:
            self.set_rel_path()

        return self.rel_path

    def get_abs_path(self):
        if self.abs_path is None:
            self.set_abs_path()

        return self.abs_path

    def print_children(self, max_width=74):
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
