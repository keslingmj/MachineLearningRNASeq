#!/usr/bin/env python3
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
TOILGENEANNOT = "~/RNA-Seq_2019/TOIL_Data/gencode.v23.annotation.gene.probemap"
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





# NEEDED FUNCTIONS:
# 1. weightedTrimmedMean
# 2. centerScaleGenes
# 3. pre-process RSEM-Counts File
def preProcess(toilSubset, TOILGENEANNOT):
    toilGeneAnnot = pd.read_csv(TOILGENEANNOT, sep="\t")
    id2gene = dict(zip(toilGeneAnnot.loc[:,'id'], toilGeneAnnot.loc[:,'gene']))
    toilSubset['gene']=[id2gene[x] for x in toilSubset.index]  # list comp
    
    def colReorder(toilSubset, regex, columnReorder):
        hits = [regex.findall(sample) for sample in list(toilSubset.columns.values)]
        hits_filter = [m for m in hits if m]
        hits_flatten = [item for sublist in hits_filter for item in sublist]
        columnReorder  = columnReorder + hits_flatten
        return columnReorder
    
    columnReorder = []
    regexes = [re.compile(r'(^GTEX.*)'), re.compile(r'(.*.11$)'),\
               re.compile(r'(.*.01$)')]
    for regex in regexes:
        columnReorder = colReorder(toilSubset, regex, columnReorder)
    
    
    regex = re.compile(r'(^GTEX.*)')
    hits = [regex.findall(sample) for sample in list(toilSubset.columns.values)]
    GTEX_hits = [m for m in hits if m]   
    columnReorder = [item for sublist in GTEX_hits for item in sublist] # nesting is left-to-right
    
    healthy = re.compile(r'(.*.11$)')
    hits = [healthy.findall(sample) for sample in list(toilSubset.columns.values)]
    hits_filter = [m for m in hits if m]
    hits_flatten = [item for sublist in hits_filter for item in sublist]
    columnReorder  = columnReorder + hits_flatten
    #toilSubset.loc[:,GTEX_flatten]  # working!
    
    
    return(id2gene.head())
# pickRefSample
# filterGenes
# addZeroGenes


# 4. Also need to create a structure for holding various dataframes.

##########
# MAIN
##########

# Read in main estimated-counts file
toilSubset = pd.read_csv(RSEM_COUNTS_FILE, sep="\t")
toilSubset.rename(index = dict(zip(toilSubset.index, toilSubset.loc[:,'sample'])),
                  inplace=True)
del toilSubset['sample']

# Dataframe cleanup
preProcess(toilSubset, TOILGENEANNOT)
