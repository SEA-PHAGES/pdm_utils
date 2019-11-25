import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pdm_utils",
    version="0.0.12",
    author="Travis Mavrich",
    author_email="trm53@pitt.edu",
    description="Phamerator database management utilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SEA-PHAGES/pdm_utils",
    packages=setuptools.find_packages(where="src"),
    package_dir={"":"src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ],
    python_requires='>=3.6',
)
