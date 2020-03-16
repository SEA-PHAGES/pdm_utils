"""Use this script to run all pipelines within the pipelines folder.
It verifies a valid pipeline is selected,
then passes all command line arguments to the main pipeline module.
"""
import argparse
import sys

from pdm_utils.pipelines import find_domains
from pdm_utils.pipelines import compare_db
from pdm_utils.pipelines import convert_db
from pdm_utils.pipelines import export_db
from pdm_utils.pipelines import freeze_db
from pdm_utils.pipelines import get_db
from pdm_utils.pipelines import get_data
from pdm_utils.pipelines import get_gb_records
from pdm_utils.pipelines import import_genome
from pdm_utils.pipelines import import_phage
from pdm_utils.pipelines import phamerate
from pdm_utils.pipelines import push_db
from pdm_utils.pipelines import update_field


def main(unparsed_args):
    """Verify a valid pipeline is selected and arguments provided are valid.

    The command line arguments are parsed and performs several basic
    checks on the arguments. Then they are passed to sub-functions to
    specifically validate the arguments based on the selected pipeline.
    """
    RUN_HELP = "Command line script to call a pdm_utils pipeline."
    USAGE = "python3 -m pdm_utils [pipeline]"
    VALID_PIPELINES = {
        "get_db",
        "get_data",
        "get_gb_records",
        "update",
        "import",
        "import_legacy",
        "find_domains",
        "phamerate",
        "export",
        "freeze",
        "compare",
        "convert",
        "push"}
    PIPELINE_HELP = "Name of the pdm_utils pipeline to run."
    pipe_parser = argparse.ArgumentParser(description=RUN_HELP, usage=USAGE)
    pipe_parser.add_argument("pipeline", type=str,
        choices=list(VALID_PIPELINES), help=PIPELINE_HELP)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = pipe_parser.parse_args(unparsed_args[1:2])
    if args.pipeline == "get_data":
        get_data.main(unparsed_args)
    elif args.pipeline == "get_db":
        get_db.main(unparsed_args)
    elif args.pipeline == "update":
        update_field.main(unparsed_args[unparsed_args.index("update") + 1:])
    # Note: import_phage is the legacy import script and will be deprecated.
    # Once import_genome is tested and operational, 'import' will call the
    # 'import_genome' module instead of the 'import_phage' module.
    elif args.pipeline == "import_legacy":
        import_phage.main(unparsed_args)
    elif args.pipeline == "import":
        import_genome.main(unparsed_args)
    elif args.pipeline == "find_domains":
        find_domains.main(unparsed_args[unparsed_args.index("find_domains") + 1:])
    elif args.pipeline == "phamerate":
        phamerate.main(unparsed_args[unparsed_args.index("phamerate") + 1:])
    elif args.pipeline == "export":
        export_db.main(unparsed_args)
    elif args.pipeline == "push":
        push_db.main(unparsed_args)
    elif args.pipeline == "freeze":
        freeze_db.main(unparsed_args)
    elif args.pipeline == "compare":
        compare_db.main(unparsed_args)
    elif args.pipeline == "get_gb_records":
        get_gb_records.main(unparsed_args)
    elif args.pipeline == "convert":
        convert_db.main(unparsed_args)
    else:
        print("Invalid pipeline selection.")
    print("\n\n\nPipeline completed")
