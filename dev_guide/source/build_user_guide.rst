Build Sphinx documentation
==========================

For Sphinx to find the entire ``pdm_utils`` package for autodoc, and to NOT autodoc any other python files in the repo (such as tests), the package directory needs to be stored within a directory that contains no other files. So the primary ``pdm_utils`` package directory is stored within 'src', and the sphinx 'config.py' file records that 'src' is where to start autodoc.

    1. Navigate to the docs directory::

        > cd /path/to/pdm_utils/docs/

    2. To initialize autodoc (this does not need to be run every time)::

        > sphinx-apidoc -o ./source/ ../src/pdm_utils/

    3. Build the docs::

        > make clean
        > sphinx-build -b html ./source ./build

    (make build may also be used instead of sphinx-build, not sure though)

This generates a preview of the html documentation. In order to push the updated documentation to ReadTheDocs:

    1. Merge all source code updates into the master git branch.
    2. Push all changes in the master branch to the GitHub repo.
    3. Login to readthedocs.org.
    4. Choose the 'pdm_utils' project.
    5. Choose 'build version'.
