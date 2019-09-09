print("Hello")
from pdm_utils.constants import constants

def main(parsed_args):
    print(constants.IMPORT_TABLE_SIZE)

    print(parsed_args)

    if parsed_args.database == 'testdb':
        func1()


def func1():
    print("executing function")
