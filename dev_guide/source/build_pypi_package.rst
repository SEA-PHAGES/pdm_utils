Build the PyPI package
======================

Follow the steps below to push a new version of the ``pdm_utils`` package to PyPI:

    1. Clear contents in the following directories:
        - build/
        - dist/

    2. Remove the following directory:
        - src/pdm_utils.egg-info

    3. Run the test suite and confirm there are no broken tests.

    4. Build the user guide to ensure there are no errors in docstring formatting.

    5. Increment the version number:

        - In the top-level directory, run the following command::

                > python3 ./change_version.py <major, minor, or micro>

        - This updates the version in four separate locations:

            - setup.py (used by PyPI for tracking package versions)
            - src/pdm_utils/__init__.py (can be accessed after installation)
            - docs/conf.py (for Sphinx documentation)

                - version number
                - release number


    .. warning::
        The version number in setup.py must be unique within the ``pdm_utils`` version history in TestPyPI and PyPI databases. Even if the package release is removed from these databases, PyPI stores its version number, and a subsequent package release cannot have the same version number, even if it has been deleted.

    6. By default, the working directory is outside of top-level pdm_utils, but the location can be specified within setup.py using the 'package_dir' and 'packages' parameters. Run the setup script from the working directory::

        > python3 setup.py sdist bdist_wheel


    7. To test the package without uploading to PyPI, install the locally-built package file::

        > pip install /path/to/dist/pdm_utils_XXXX.tar.gz

    8. In a new terminal, open a Python IDE and test the package.

    9. After testing locally, upload the package to TestPyPI::

        > python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*

    10. Open a separate terminal that doesn't utilize the pdm_utils-dev conda environment (for instance, load the pdm_utils-user conda environment), and install the package with pip::

        > python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps pdm_utils

    11. Now upload the package to PyPI::

        > twine upload ./dist/*

    12. Build a new user guide on ReadTheDocs that reflects the new version.

    13. Commit all changes to git and push to GitHub::

        > git commit -m 'incremented package version to <version number>'

    14. Tag the commit with the new package version and push tags to GitHub::

        > git tag -a v<version number> -m 'pdm_utils <version number>'
        > git push origin --tags
