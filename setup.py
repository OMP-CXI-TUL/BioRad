import setuptools


__version__ = '1.0'

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="BioRad", # Replace with your own username
    version=__version__,
    author="Jakub Riha",
    author_email="jakub.riha@tul.cz",
    description="BioRad SW",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/OMP-CXI-TUL/BioRad",
    packages=setuptools.find_packages(where='src'),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    install_requires=['numpy', 'scipy', 'pyyaml', 'pysqlite3', 'matplotlib'],
    python_requires='>=3.7',
)