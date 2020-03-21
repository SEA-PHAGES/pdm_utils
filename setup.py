import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pdm_utils",
    version="0.1.1",
    author="Travis Mavrich",
    author_email="trm53@pitt.edu",
    description="MySQL phage genomics database management utilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SEA-PHAGES/pdm_utils",
    packages=setuptools.find_packages(where="src"),
    package_dir={"":"src"},
    install_requires=[
       'biopython>=1.70',
       'pymysql>=0.9.0',
       'paramiko>=2.5.0',
       'tabulate>=0.8.0',
       'sqlalchemy>=1.3.0'
    ],
    project_urls={
        'Documentation': 'https://pdm-utils.readthedocs.io/en/latest/',
        'Source': 'https://github.com/SEA-PHAGES/pdm_utils/',
    },
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research",
        "Topic :: Database",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires='>=3.6'
)
