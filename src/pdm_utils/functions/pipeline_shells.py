from cmd import Cmd
from pathlib import Path
import shlex
import textwrap

from urllib3.exceptions import RequestError

from pdm_utils.functions import (basic, url_basic)


# SHELL_INTROS
# -----------------------------------------------------------------------------
GET_DB_CMD_INTRO = ("Loaded the pdm_utils get_db interactive environment:\n\n"
                    "Type select [package] to choose a database to install\n"
                    "Use cd [path/to/directory] or ls [path/to/directory] to "
                    "navigate remote contents\n"
                    "Type exit to leave the interactive environment\n\n"
                    "Type help or ? to list more commands.")


# GET_DB_SHELL
# -----------------------------------------------------------------------------
GET_DB_FILE_CMD_COMMANDS = ["desc", "select", "version"]
GET_DB_DIR_CMD_COMMANDS = ["cd"]
GET_DB_PATH_CMD_COMMANDS = (GET_DB_FILE_CMD_COMMANDS +
                            GET_DB_DIR_CMD_COMMANDS + ["ls"])
GET_DB_ARGLESS_COMMANDS = ["exit", "reload"]
GET_DB_SINGLE_ARG_COMMANDS = ["which", "load"]
GET_DB_CMD_COMMANDS = (GET_DB_PATH_CMD_COMMANDS + GET_DB_ARGLESS_COMMANDS +
                       GET_DB_SINGLE_ARG_COMMANDS)


# Inherits from Cmd library:
# Command line interactive faux shell that allows for database selection at
# some url
class GetDBCMD(Cmd):
    def __init__(self, url, name="root", preload=False, pool=None, **kwargs):
        super().__init__(**kwargs)

        if pool is None:
            pool = url_basic.create_pool(pipeline=True)

        self.pool = pool
        self.name = name
        self.preload = preload
        self.load_graph(url)

        self.selected = None
        self.tmp = None
        self.text_wrapper = textwrap.TextWrapper(width=80, initial_indent="\t")

        self.doc_header = ("Commands: (type help [command] for more info)")

    def set_prompt(self):
        prompt = "".join([self.name, "@", self.url])

        if self.dir_graph is not None:
            prompt = "".join([prompt, "::",
                              self.dir_graph.curr.get_rel_path()])

        self.prompt = "".join([prompt, "$ "])

    def load_graph(self, url):
        self.url = url

        self.dir_graph = None
        try:
            self.dir_graph = url_basic.UrlPathGraph(
                                            url, preload=True, pool=self.pool)
        except RequestError as e:
            print(e)
            return

        if self.preload:
            for node in self.dir_graph.nodes:
                if not node.isdir:
                    self.load_yaml_data(node)

    def load_yaml_data(self, node):
        for file_name in node.contents:
            if file_name.endswith(".yml"):
                file_url = "".join([node.get_abs_path(), "/", file_name])

                try:
                    response = url_basic.pool_request(file_url,
                                                      pool=self.dir_graph.pool)
                except RequestError as e:
                    print(e)

                node.data = url_basic.get_yaml_data(response)

    def preloop(self):
        self.set_prompt()

    def precmd(self, line):
        split_args_list = shlex.split(line)
        if split_args_list:
            command = split_args_list[0]

            if command in GET_DB_PATH_CMD_COMMANDS:
                self.tmp = self.parse_path_command(
                                                    split_args_list, command)

            if command in GET_DB_SINGLE_ARG_COMMANDS:
                self.tmp = self.parse_single_arg_command(
                                                    split_args_list, command)
        return line

    def postcmd(self, stop, line):
        self.set_prompt()
        self.tmp = None
        return stop

    # This functionality could be replaced/improved with python's argparse
    # library, however argparse exits on error.  This functionality is
    # toggleable with python3 v3.9, but this is the most recent python3
    # version as of 01/01/2021.  Could be updated as this version becomes more
    # commonly used.
    def parse_path_command(self, split_args_list, command):
        if self.dir_graph is None:
            print("".join(["Request failed.  "
                           "Please check the availability at '", self.url,
                           "'\nLoad a new url with the following:\n"
                           "\tload [url]"]))
            return None

        if len(split_args_list) == 1:
            path = ""
        else:
            path = split_args_list[1]

        if command in GET_DB_FILE_CMD_COMMANDS:
            node_type = "file"
        elif command in GET_DB_DIR_CMD_COMMANDS:
            node_type = "dir"
        else:
            node_type = None

        return self.get_path_node(path, command, node_type)

    def parse_single_arg_command(self, split_args_list, command):
        if len(split_args_list) == 1:
            arg = ""
        else:
            arg = split_args_list[1]

        return arg

    def get_path_node(self, path, command="", node_type=None):
        corrected_path = str(Path(path))
        path_node = self.dir_graph.traverse(corrected_path)

        if path_node is not None:
            if node_type == "file":
                if path_node.isdir:
                    print("".join([command, ": ", path, ": Not a package"]))
                    path_node = None
            elif node_type == "dir":
                if not path_node.isdir:
                    print("".join([command, ": ", path, ": Not a directory"]))
                    path_node = None
        else:
            print("".join([command, ": ", path,
                           ": Not a package or directory"]))

        return path_node

    def do_ls(self, arg):
        ("List directory contents\n"
         "\tUsage: ls [path/to/dir]")
        node = self.tmp
        if node is not None:
            if node.isdir:
                node.print_children()
            else:
                print(node.get_rel_path())

    def do_cd(self, arg):
        ("Change directory\n"
         "\tUsage: cd [path/to/dir]")
        node = self.tmp
        if node is not None:
            self.dir_graph.curr = node

    def do_desc(self, arg):
        ("Describe database package\n"
         "\tUsage: desc [path/to/package]")
        node = self.tmp
        if node is not None:
            if node.data is None:
                self.load_yaml_data(node)

            if node.data is not None:
                lines = []
                if node.data.get("name") is not None:
                    lines.append("Name:")
                    lines += self.text_wrapper.wrap(node.data["name"])

                if node.data.get("date") is not None:
                    lines.append("Date:")
                    lines += self.text_wrapper.wrap(str(node.data["date"]))

                if node.data.get("description") is not None:
                    lines.append("Description:")
                    lines += self.text_wrapper.wrap(node.data["description"])

                for line in lines:
                    print(line)

    def do_version(self, arg):
        ("Version database package\n"
         "\tUsage: version [path/to/package]")
        node = self.tmp
        if node is not None:
            if node.data is None:
                self.load_yaml_data(node)

            if node.data is not None:
                print(": ".join([
                        "Database version", str(node.data.get("version"))]))

    def do_select(self, arg):
        ("Select database package\n"
         "\tUsage: select [path/to/package]")
        node = self.tmp
        if node is not None:
            self.selected = node
            return True

    def do_which(self, arg):
        ("Locate a directory/package\n"
         "\tUsage: which [name]")
        arg = self.tmp
        for node in self.dir_graph.nodes:
            if node.name == arg:
                print(node.get_rel_path())

    def do_clear(self, arg):
        "Bring the command line to the top of screen"
        basic.clear_screen()

    def do_reload(self, arg):
        ("Reload root contents at the current url")
        self.load_graph(self.dir_graph.url)

    def do_load(self, arg):
        ("Reload root contents from a new url\n"
         "\tUsage: load [url]")
        self.load_graph(self.tmp)

    def do_exit(self, arg):
        ("Exits the interactive shell")
        return True
