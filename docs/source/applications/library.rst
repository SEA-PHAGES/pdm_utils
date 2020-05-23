Library tutorial
================

``pdm_utils`` provides a library of functions, classes, and methods that leverage SQLAlchemy, a highly-refined, powerful, well-supported 'Database Toolkit for Python'. Tailoring SQLAlchemy tools to the pre-defined pdm_utils MySQL database provide for diverse entry points to the database. This tutorial provides a brief introduction to how the library can be used.

.. toctree::
   :maxdepth: 1

   Direct connection to MySQL <./tutorial/mysql_connection>
   Back-end pdm_utils ORM <./tutorial/back_end_orm>
   Dynamic querying and data filtering <./tutorial/filtering_querying>
   Front-end SQLAlchemy ORM <./tutorial/front_end_orm>




Each section of the tutorial assumes a Python IDE has been initiated within a terminal environment with all dependencies available. In the shell terminal, activate the Conda environment containing the installed ``pdm_utils`` package (if needed) to ensure all dependencies are present. Then open a Python IDE::

    > conda activate pdm_utils
    (pdm_utils)>
    (pdm_utils)> python3
    >>>
