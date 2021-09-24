import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pdm_utils",
    version="0.9.16",
    author="Christian Gauthier",
    author_email="chg60@pitt.edu",
    description="MySQL phage genomics database management utilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SEA-PHAGES/pdm_utils",
    packages=setuptools.find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        'biopython==1.77',
        'networkx>=2.4',
        'paramiko>=2.7.2',
        'pymysql==0.9.3',
        'pyyaml>=5.3.1',
        'sqlalchemy==1.4.2',
        'tabulate>=0.8.3',
        'urllib3>=1.25.8'
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
