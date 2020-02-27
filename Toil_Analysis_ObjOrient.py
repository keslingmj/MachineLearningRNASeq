#!/Users/mjk/opt/anaconda3/bin/python3
#### #!/Users/mjk/.virtualenvs/ML/bin/python3
# -*- coding: utf-8 -*-
"""
Toil_Analysis_ObjOrient.py
This script rewrites the same script previously written in .Rmd

Created on Wed Dec 18 15:15:24 2019

@author: mjk
"""

##########
# SET GLOBAL PARAMETERS
##########
RSEM_COUNTS_FILE = "/Users/mjk/Desktop/Tresorit_iOS/projects/RNA-Seq/MachineLearningRNASeq/toilSubsetRSEM382.txt"
TOILGENEANNOT = "/Users/mjk/RNA-Seq_2019/TOIL_Data/gencode.v23.annotation.gene.probemap"
BIOMART_GENE_ATTR_FILE = "/Users/mjk/Desktop/Tresorit_iOS/projects/RNA-Seq/MachineLearningRNASeq/geneAttr.csv"
TCGA_ATTR_FILE = "~/Desktop/Tresorit_iOS/projects/RNA-Seq/data/TCGA_Attributes_full.txt"
WANG_MATRIX = "/Users/mjk/Desktop/Tresorit_iOS/projects/RNA-Seq/data/wangBreastFPKM398_Attrib.txt"
FULL_TOIL_DATA = "/Users/mjk/RNA-Seq_2019/TOIL_Data/TcgaTargetGtex_gene_expected_count"
SEED = 233992812
SEED2 = 1011



###########
## declare packages
###########
import pandas as pd
import re
import numpy as np
import matplotlib as plt
import random
import rpy2.robjects as robjects
from sklearn.model_selection import train_test_split
from rpy2.robjects.packages import importr
from itertools import compress
import copy




# NEEDED FUNCTIONS:
# 1. weightedTrimmedMean
# 2. centerScaleGenes
##########
# 3. pre-process RSEM-Counts File
##########
def processColReorder(df):
    # this code orders the samples by 1. GTEX, 2, TCGA-healthy
    # 3. TCGA-tumor.  Within each of these 3 groups, there is no
    # particular ordering.  Assumes samples are columns
    def colReorder(df, regex, columnReorder):
        hits = [regex.findall(sample) for sample in list(df.columns.values)]
        hits_filter = [m for m in hits if m]
        hits_flatten = [item for sublist in hits_filter for item in sublist]
        columnReorder  = columnReorder + hits_flatten
        return columnReorder
    
    columnReorder = []
    regexes = [re.compile(r'(^GTEX.*)'), re.compile(r'(.*.11$)'),\
               re.compile(r'(.*.01$)')]
    for regex in regexes:
        columnReorder = colReorder(df, regex, columnReorder)
    tmp = df.loc[:, columnReorder]
    return(tmp)


def geneNameMod(toilSubset, TOILGENEANNOT):
    # this code merges gene ID with common gene name,
    # and otherwise returns the toilSubset dataframe unaltered
    
    toilGeneAnnot = pd.read_csv(TOILGENEANNOT, sep="\t", engine="python")
    id2gene = dict(zip(toilGeneAnnot.loc[:,'id'], toilGeneAnnot.loc[:,'gene']))
    toilSubset['gene']=[id2gene[x] for x in toilSubset.index]  # list comp
    

    concatGeneName = toilSubset.loc[:,'gene'] + "-" + toilSubset.index
    toilSubset.rename(index=dict(zip(toilSubset.index, concatGeneName)), inplace=True)
    toilSubset.drop(["gene"], axis=1, inplace=True) 
    return(toilSubset)


# pickRefSample
# filterGenes
    
## NOT FINISHED
def filterGenes(x, cutoff=0.2):    # assumes data in log2-format
    return(x)
    
# addZeroGenes


# 4. Also need to create a structure for holding various dataframes.
class RSeqDF(pd.DataFrame):
    def __init__(self, name, data=None, index=None, columns=None,
                 dtype=None, copy=False):
        pd.DataFrame.__init__(self, data=data, index=index, columns=columns,
                 dtype=dtype, copy=copy)
        self.name = name
        
    def __len__(self):              # returns number of rows in pd.DF
        return self.shape[0]
    
    def __getitem__(self):
        return self

class DF_Set:
    def __init__(self, name, data=None, index=None, columns=None,
                 dtype=None, copy=False):
        objNames = "orig filt zeroGenes nat scaled norm outcome refSmpName \
        refSmp refSmpUnscaled filtScaled".split()
        self._dfs = [RSeqDF('orig', data=data, index=index, \
                 columns=columns, dtype=None, copy=copy)]
        self.nameIdx = dict(zip(objNames, range(11)))
        
    def __getitem__(self, name):                    # looks up DF based on name
        return self._dfs[self.nameIdx[name]]             # name -> posn in list
    
    def filterGenes(self, cutoff=0.2):                     # assumes log2 scale
        ZeroExpGenes = list(compress(self._dfs[0], \
                                     self._dfs[0].quantile(.75, axis=0) < cutoff))
        orig = copy.copy(self._dfs[0])
        self.filt = orig.drop(ZeroExpGenes, axis=1, inplace=False)
        
        self._dfs.append(self.filt)
        self._dfs.append(ZeroExpGenes)
        #return self.filt
    
    def natScale(self):                         # full orig df on natural scale
        self._dfs.append(np.exp2(self._dfs[0]) - 1)
        
        
        
    #def name2Idx(self, name):
    #    return self.nameIdx[name]


##########
# MAIN
##########

# 1. Read in main estimated-counts file
toilSubset = pd.read_csv(RSEM_COUNTS_FILE, sep="\t")
toilSubset.rename(index = dict(zip(toilSubset.index, toilSubset.loc[:,'sample'])),
                  inplace=True)
toilSubset.drop(['sample'], axis=1, inplace=True)

# 2. gene name modification and sort samples
toilSubset = geneNameMod(toilSubset, TOILGENEANNOT)
toilSubset = processColReorder(toilSubset)

# 3. Create Training and Test Sets from toilSubset
# A. transpose the matrix
toilSubsetWide = toilSubset.T

# B. create outcome variable and give it a common index (sample name)

# add method for automatically determining #healthy and #tumor

outcome = pd.DataFrame(pd.concat([pd.Series(np.zeros(185)), pd.Series(np.ones(197))],
                    ignore_index=True))
outcome.rename(index=dict(zip(outcome.index,toilSubsetWide.index)), inplace=True)
all(outcome.index == toilSubsetWide.index)


# C. split df and outcome into training and test sets
toilTrain, toilTest, outcomeTrain, outcomeTest = \
    train_test_split(toilSubsetWide, outcome, test_size=0.25, random_state=SEED)
# but did this maintain the ratio of 'healthy' and 'tumor' samples?

# D. column reorder: GTEX, TCGA-healthy, TCGA-tumor for all 4 df's
toilTrain = (processColReorder(toilTrain.T)).T
outcomeTrain = outcomeTrain.loc[toilTrain.index,:]
toilTest = (processColReorder(toilTest.T)).T
outcomeTest = outcomeTest.loc[toilTest.index,:]


# E. QA test and training sets
all(outcomeTrain.index == toilTrain.index)
all(outcomeTest.index == toilTest.index)
sum(outcomeTrain[0])/len(outcomeTrain.index) # tests even split healthy/tumor

# 4. Create Data Object
# FIRST, GET TO WORK.  THEN CHECK ACCURACY, THEN, WILL WRITE UP PROPERLY
train_orig = RSeqDF('orig', data=toilTrain, index=toilTrain.index, columns=toilTrain.columns)
train_norm = train_orig.transform(lambda x: 2**x - 1)
geneQuantiles = train_orig.quantile(q=0.75, axis='rows') # counter-intuitive
ZeroGenes = geneQuantiles < 0.2
train_filt = train_orig.loc[:, -ZeroGenes]
#edgeR_scalingFactors = edgeR::calcNormFactors(t(train$M_filt),method="TMM")


# 5. Filter Out Zero-Expressed Genes
#r = robjects.r
limma = importr('limma')
edgeR = importr('edgeR')
robjects.r("library(edgeR)")   # working now.
robjects.r("version")

