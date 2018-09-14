
# coding: utf-8

# ## Importing libraries

# In[1]:


from sys import argv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys
from rdkit.Chem import PandasTools
from statsmodels.stats.proportion import proportions_ztest


# ## Definig functions

# In[2]:


def LoadDatasetFromCSV(CSV, ecfp4=True, maccs=True, label="ACTIVE"):
# This function requires a CSV file with the first two columns identified as ID and SMILES, additional columns will be ignored
# The file is loaded as a dataframe and by default both fingerprints, MACCS-166 and ECFP4-2048 are calculated and added as columns
# The last column is added as a label, which by default is "ACTIVE"
    Dataset = pd.read_csv(CSV, usecols=[1])
    if ecfp4 == True:
        ECFP4FP = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(row.SMILES),2,nBits=2048) for index, row in Dataset.iterrows()]
        Dataset["ECFP4FP"] = ECFP4FP
    if maccs == True:
        MACCSFP = [MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(row.SMILES)) for index, row in Dataset.iterrows()]
        Dataset["MACCSFP"] = MACCSFP
    LABEL = [label for index, row in Dataset.iterrows()]
    Dataset["LABEL"] = LABEL
    return Dataset

def DFP_Calc(DF, FP="ECFP4", FORMAT="RDKit"):
# This function requires a dataframe with a column identified as "ECFP4FP" or "MACCSFP" containing the respective fingerprints as RDKit objects
# The input dataframe can be taken from the LoadDatasetFromCSV function
# FP = "ECFP4" or "MACCS" according to the respective DFP
# FORMAT = "RDKit" or "TEXT" according to the output format, RDKit object or TEXT string
    
    if FP == "ECFP4":
        FPSTEXT = [DataStructs.BitVectToText(row.ECFP4FP) for index, row in DF.iterrows()]
        DF_COUNTS = [0 for i in range(len(FPSTEXT[0]))]
        for i in FPSTEXT:
            b = [int(j) for j in i]
            DF_COUNTS = [x + y for x, y in zip(DF_COUNTS, b)]
        DF_PROPORTIONS = [float(x)/DF.shape[0] for x in DF_COUNTS]
        DFP = []
        for i in range(0, len(DF_PROPORTIONS)):
            if DF_PROPORTIONS[i] > 0.5:
                DFP.append(1)
            else:
                DFP.append(0)
        DFP = [str(i) for i in DFP]
        DFP = "".join(DFP)
        DFP_RDKIT = DataStructs.CreateFromBitString(DFP)
    elif FP == "MACCS":
        FPSTEXT = [DataStructs.BitVectToText(row.MACCSFP) for index, row in DF.iterrows()]
        DF_COUNTS = [0 for i in range(len(FPSTEXT[0]))]
        for i in FPSTEXT:
            b = [int(j) for j in i]
            DF_COUNTS = [x + y for x, y in zip(DF_COUNTS, b)]
        DF_PROPORTIONS = [float(x)/DF.shape[0] for x in DF_COUNTS]
        DFP = []
        for i in range(0, len(DF_PROPORTIONS)):
            if DF_PROPORTIONS[i] > 0.5:
                DFP.append(1)
            else:
                DFP.append(0)
        DFP = [str(i) for i in DFP]
        DFP = "".join(DFP)
        DFP_RDKIT = DataStructs.CreateFromBitString(DFP)
    if FORMAT == "RDKit":
        return DFP_RDKIT
    elif FORMAT == "TEXT":
        return DFP

def SBDFP_Calc(DF, FP="ECFP4", FORMAT="RDKit"):
# This function requires a dataframe with a column identified as "ECFP4FP" or "MACCSFP" containing the respective fingerprints as RDKit objects
# The function also requires the files ECFP4.counts or MACCS.counts that contain the "1" Bit counts for the respective fingerprints
# The input dataframe can be taken from the LoadDatasetFromCSV function
# FP = "ECFP4" or "MACCS" according to the respective SB-DFP
# FORMAT = "RDKit" or "TEXT" according to the output format, RDKit object or TEXT string

    if FP == "ECFP4":
        FPSTEXT = [DataStructs.BitVectToText(row.ECFP4FP) for index, row in DF.iterrows()]
        DF_COUNTS = [0 for i in range(len(FPSTEXT[0]))]
        for i in FPSTEXT:
            b = [int(j) for j in i]
            DF_COUNTS = [x + y for x, y in zip(DF_COUNTS, b)]
        REF = open("ECFP4.counts")
        line = REF.readline()
        a = line.split(",")
        REF_COUNTS = [int(x) for x in a]
        SBDFP = []
        for i in range(len(REF_COUNTS)):
            stat, pval = proportions_ztest([REF_COUNTS[i], DF_COUNTS[i]], [15403690,DF.shape[0]], alternative='smaller')
            if pval < 0.01:
                SBDFP.append(1)
            else:
                SBDFP.append(0)
        SBDFP = [str(x) for x in SBDFP]
        SBDFP = "".join(SBDFP)
        SBDFP_RDKIT = DataStructs.CreateFromBitString(SBDFP)
    
    elif FP == "MACCS":
        FPSTEXT = [DataStructs.BitVectToText(row.MACCSFP) for index, row in DF.iterrows()]
        DF_COUNTS = [0 for i in range(len(FPSTEXT[0]))]
        for i in FPSTEXT:
            b = [int(j) for j in i]
            DF_COUNTS = [x + y for x, y in zip(DF_COUNTS, b)]
        REF = open("MACCS.counts")
        line = REF.readline()
        a = line.split(",")
        REF_COUNTS = [int(x) for x in a]
        SBDFP = []
        for i in range(len(REF_COUNTS)):
            stat, pval = proportions_ztest([REF_COUNTS[i], DF_COUNTS[i]], [15403690,DF.shape[0]], alternative='smaller')
            if pval < 0.01:
                SBDFP.append(1)
            else:
                SBDFP.append(0)
        SBDFP = [str(x) for x in SBDFP]
        SBDFP = "".join(SBDFP)
        SBDFP_RDKIT = DataStructs.CreateFromBitString(SBDFP)
    if FORMAT == "RDKit":
        return SBDFP_RDKIT
    elif FORMAT == "TEXT":
        return SBDFP


# ## Examples
# For using the function, the following files are needed: MACCS.counts and ECFP4.counts.
# For the execution of the scrip DNMT1.csv is also needed as example.
# All of them are available in the GitHub repository.

# In[10]:


# Loading a dataset with default parameters
Dataset = LoadDatasetFromCSV("DNMT1.csv")
# Showing the first 5 entries
Dataset.head(5)


# In[17]:


# Calculating DFP/MACCS for the generated dataframe
Dataset_DFP_MACCS = DFP_Calc(Dataset,FP="MACCS",FORMAT="TEXT")
# Showing the DFP
print ("Dataset_DFP_MACCS: "+Dataset_DFP_MACCS)


# In[16]:


# Calculating SB-DFP/ECFP4 for the generated dataframe
Dataset_SBDFP_ECFP4 = SBDFP_Calc(Dataset,FP="ECFP4",FORMAT="TEXT")
# Showing the SB-DFP
print ("Dataset_SBDFP_ECFP4: "+Dataset_SBDFP_ECFP4)

