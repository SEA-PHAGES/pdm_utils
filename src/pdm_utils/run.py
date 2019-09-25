"""Use this script to run all pipelines within the pipelines folder.
It verifies a valid pipeline is selected,
then passes all command line arguments to the main pipeline module.
"""
import sys
import argparse
from pdm_utils.pipelines.db_import import import_genome

def main():
    """Verify a valid pipeline is selected and arguments provided are valid.

    The command line arguments are parsed and performs several basic
    checks on the arguments. Then they are passed to sub-functions to
    specifically validate the arguments based on the selected pipeline.
    """
    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    RUN_HELP = "Command line script to call a pdm_utils pipeline."
    VALID_PIPELINES = {
        "retrieve_data",
        "import",
        "cdd",
        "phamerate",
        "export",
        "compare",
        "database_to_file"}
    PIPELINE_HELP = "Name of the pdm_utils pipeline to run."
    pipe_parser = argparse.ArgumentParser(description=RUN_HELP)
    pipe_parser.add_argument("pipeline", type=str,
                        choices=list(VALID_PIPELINES),
                        help=PIPELINE_HELP)
    args = pipe_parser.parse_args(sys.argv[1:2])
    if args.pipeline == "retrieve_data":
        pass
    elif args.pipeline == "import":
        import_genome.run_import(sys.argv)
    elif args.pipeline == "cdd":
        pass
    elif args.pipeline == "phamerate":
        pass
    elif args.pipeline == "export":
        pass
    # TODO eventually 'database_to_file' will be merged into 'export' pipeline.
    elif args.pipeline == "database_to_file":
        pass
    else:
        pass
    print("Pipeline completed")

if __name__ == "__main__":
    main()
