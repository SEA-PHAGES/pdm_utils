.. _install_conda:

Anaconda
========

``pdm_utils`` requires Python >=3.6 and several third-party Python packages:

    - :biopython:`Biopython <>`
    - :networkx:`NetworkX <>`
    - :paramiko:`Paramiko <>`
    - :pymysql:`PyMySQL <>`
    - :pyyaml:`PyYAML <>`
    - :sqlalchemy:`SQLAlchemy <>`
    - :tabulate:`Tabulate <>`

.. note::

    The most recent version of Biopython introduces some problems for the current pdm_utils. We are working on fixing it, but in the meantime, we have edited the conda environment creation command to install the latest non-code-breaking version of Biopython.

Some of them also have Python or non-Python dependencies. Manual installation of these dependencies can be tricky, but the Conda environment manager is a simple, automated alternative. First install Conda, then use Conda to install Python and other dependencies.

Conda is available as part of Anaconda or Miniconda, and complete installation instructions are available in the Conda :conda_docs:`user guide <>`. ``pdm_utils`` has been developed and tested using Conda >=4.5. The directions below highlight installation of Anaconda, but either of these tools is fine since they both install Conda:

#. Install Conda locally through the :anaconda:`Anaconda <>` package.

#. Navigate to the 'Anaconda Distribution' option.

#. Begin the installation:

    - For MacOS: download the Python 3.8 graphical installer and follow the graphical prompts.

    - For Linux:

        A. Download the Python 3.8 x86 Linux installer (e.g. Anaconda3-2019.10-Linux-x86_64.sh) to the Downloads folder.
        B. Open a Terminal window.
        C. Execute the following command::

            > bash ~/Downloads/Anaconda3-2019.10-Linux-x86_64.sh


#. Follow the manufacturer's installation instructions.

    - Accept the license agreement.
    - Install at the default directory.
    - Enter 'yes' when prompted for the installer to run conda init.

#. Optional: execute the following command to prevent Conda from starting up automatically every time a Terminal window is opened::

    > conda config --set auto_activate_base false

#. Close the Terminal window and open a new window.

#. After installing Conda, create an environment to be able to install and use ``pdm_utils`` (the example below creates a Conda environment named 'pdm_utils', but it can be named anything). Enter 'y' when prompted to install all dependencies::

    > conda create --name pdm_utils curl python pip biopython==1.77 networkx paramiko pymysql sqlalchemy tabulate urllib3

#. After the Conda environment is created, it needs to be activated using the following command. The command line prompt will now include '(pdm_utils)', indicating it is operating within this environment::

    > conda activate pdm_utils
    (pdm_utils)>

#. Optional: enter the following command to exit the Conda environment. The default command line prompt will be displayed, and the name of the Conda environment will no longer be displayed::

    (pdm_utils)> conda deactivate
    >


.. note::

    If Conda is used to manage dependencies, the Conda environment must be activated every time you want to use ``pdm_utils``. Otherwise, an error will be encountered.


The 'pdm_utils' Conda environment now contains several required dependencies, and the actual ``pdm_utils`` Python package can be (:ref:`installed <install_pdm_utils_package>`).
