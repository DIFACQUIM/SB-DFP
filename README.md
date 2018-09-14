# SB-DFP
Statistical-Based Database Fingerprint: Chemical space dependent representation of compound databases

Epigenetic_datasets.tar.gz
  Contains the 28 datasets described in the work in CSV format.

ZINC12_AllClean.tar.gz.part00 to 07
  Contain the following 4 files in CSV format:
    Not_Processed.csv         |     Contains 21 structures not processed by rdkit.
    Repeated_Structures.csv   |     Contains 154 repeated compounds between ZINC 12 AllClean subset and Epigenetic datasets.
    Reference_Dataset.csv     |     Contains 15,403,690 structures used as reference for SB-DFP calculations.
    Decoys.csv                |     Contains 1 million compounds used as decoys in similarity searching runs.
  To uncompress:
    1) Join the 8 parts into a single compressed file with: cat ZINC12_AllClean.tar.gz.part* >> ZINC12_AllClean.tar.gz and
    2) Exctract the files with: tar -xvzf ZINC12_AllClean.tar.gz
    
SB-DFPCalc.ipynb and SB-DFPCalc.py
  The Jupyter notebook with the code employed for DFP and SB-DFP calculations with some examples, and the Python script of such notebook.
  
MACCS.counts and ECFP4.counts
   Files containing the reference counts for the calculation of SB-DFP based on MACCS and ECFP4 respectively.
   These files are needed for the execution of the Python code.
  
DNMT1.csv
  This file is needed for the execution of the Python code (used as example).
  This file is also contained in Epigenetic_datasets.tar.gz
