"""Use this script to run all pipelines within the pipelines folder.
It verifies a valid pipeline is selected,
then passes all command line arguments to the main pipeline module.
"""
import argparse

from pdm_utils.pipelines import compare_db
from pdm_utils.pipelines import convert_db
from pdm_utils.pipelines import export_db
from pdm_utils.pipelines import find_domains
from pdm_utils.pipelines import find_transmembrane
from pdm_utils.pipelines import freeze_db
from pdm_utils.pipelines import get_data
from pdm_utils.pipelines import get_db
from pdm_utils.pipelines import get_gb_records
from pdm_utils.pipelines import import_genome
from pdm_utils.pipelines import phamerate
from pdm_utils.pipelines import pham_finder
from pdm_utils.pipelines import push_db
from pdm_utils.pipelines import revise
from pdm_utils.pipelines import pham_review
from pdm_utils.pipelines import update_field

VALID_PIPELINES = {"compare", "convert", "export", "find_domains",
                   "find_transmembrane", "find_phams", "freeze", "get_data",
                   "get_db", "get_gb_records", "import", "phamerate", "push",
                   "revise", "pham_review", "update"}


def main(unparsed_args):
    """Run a pdm_utils pipeline."""
    args = parse_args(unparsed_args)

    if args.pipeline == "compare":
        compare_db.main(unparsed_args)
    elif args.pipeline == "convert":
        convert_db.main(unparsed_args)
    elif args.pipeline == "find_transmembrane":
        find_transmembrane.main(unparsed_args)
    elif args.pipeline == "export":
        export_db.main(unparsed_args)
    elif args.pipeline == "find_domains":
        find_domains.main(unparsed_args)
    elif args.pipeline == "find_phams":
        pham_finder.main(unparsed_args)
    elif args.pipeline == "freeze":
        freeze_db.main(unparsed_args)
    elif args.pipeline == "get_data":
        get_data.main(unparsed_args)
    elif args.pipeline == "get_db":
        get_db.main(unparsed_args)
    elif args.pipeline == "get_gb_records":
        get_gb_records.main(unparsed_args)
    elif args.pipeline == "import":
        import_genome.main(unparsed_args)
    elif args.pipeline == "phamerate":
        phamerate.main(unparsed_args)
    elif args.pipeline == "push":
        push_db.main(unparsed_args)
    elif args.pipeline == "revise":
        revise.main(unparsed_args)
    elif args.pipeline == "pham_review":
        pham_review.main(unparsed_args)
    elif args.pipeline == "update":
        update_field.main(unparsed_args)
    else:
        pass
    print("\n\nPipeline completed")


def parse_args(unparsed_args):
    """
    Use argparse to verify pipeline argument only.

    :param unparsed_args: raw command line args
    :type unparsed_args: list
    """
    run_help = "Command line script to call a pdm_utils pipeline."
    usage = "python3 -m pdm_utils [pipeline]"
    pipeline_help = "Name of the pdm_utils pipeline to run."

    parser = argparse.ArgumentParser(description=run_help, usage=usage)
    parser.add_argument("pipeline", type=str, choices=list(VALID_PIPELINES),
                        help=pipeline_help)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args[1:2])

    return args
