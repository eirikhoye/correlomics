Correlomics:uncovering microRNA editing and structural features based on reference and sequencing data

Pipeline for smallRNA sequencing data. 

Processing and quality control of sequencing data using MirTrace. Returns fasta files with collapsed reads.

Generates microRNA reference library from MirGeneDB mature and star microRNA sequences, which are extended with 5 nucleotides at the 5p and 3p end from the genomic loci. 

Generate bowtie reference for each species, bowtie v=1.0.0. 

Map collapsed reads, allowing up to 10 alignments with up to 3 mismatches. Outputs sam files.

Use two separate filtering methods for sam files, to filter for multiple read alignments:
- One selecting for mismatches INSIDE miRNA position 2-18
    - Used for identifying microRNA editing events, including A to I editing.
- One selecting for mismatches OUTSIDE miRNA position 2-18
    - Used for analyzing the preferred read alignment position (0 is canonical start/end).




generates bowtie reference 



Features to add:
- make a sensible file structure
- make it snakemake compatible
- package directory only contain required code, no data
- input data to be specified with YAML file
- use virtualenv not conda (figure out how to make this work with snakemake)



CHATGPT STRUCTURE SUGGESTION:
A common file structure for a Python package is as follows:

mypackage/
    LICENSE
    README.md
    setup.py
    mypackage/
        __init__.py
        module1.py
        module2.py
        ...
    tests/
        __init__.py
        test_module1.py
        test_module2.py
        ...

The top-level directory should contain a LICENSE file, a README.md file, and a setup.py file.
The package's code should be in a directory with the same name as the package.
The package's code should include an __init__.py file, which can be empty or can be used to set up the package.
Any modules or sub-packages should also be included in this directory.
The tests/ directory should contain test files for the package. It should also include an __init__.py file.
It is important to use standard naming conventions for package and module name, also you should follow the PEP8 style guide for writing python code.



TODO
Continue the integration from mirgenedb_editing_events following the flow chart png. 

Start a SnakeMake file and make compatible scripts for each prosessing step.