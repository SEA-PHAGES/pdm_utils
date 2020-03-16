Code conventions
================

This repo utilizes the following coding conventions:

    1. Use spaces instead of tabs.
    2. Use 'snake_case' for variables and 'PascalCase' for classes.
    3. Docstrings:

        - for all methods and functions.
        - written in reStructuredText for Sphinx autodoc.

    4. Use double (instead of single) quotes for string literals (although not as important in docstrings and comments)
    5. Tests:

        - should be written for all methods and functions.
        - should be constructed for the ``unittest`` module.
        - should have a docstring that briefly states the purpose of the test (although doesn't need to be specifically structured).
        - are split into unit and integration test directories. If the test relies on pure python, it should be stored in the 'unit' directory. These tests run very quickly. If it relies on MySQL, PhagesDB, parsing files, creating files and directories, etc. it should be stored in the 'integration' directory. These tests run more slowly.
