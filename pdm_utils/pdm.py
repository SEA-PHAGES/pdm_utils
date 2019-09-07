"""Use this script to run all pipelines within the pipelines folder.
"""
# This avoids problems wtih import scope.
import sys
# from pipelines.db_import import test
import pipelines

pipeline = sys.argv[1]


if pipeline == 'testrun':
    pipelines.db_import.test.main()

print("Finished")
