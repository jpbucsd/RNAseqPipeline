#it is strange to perform PCA with so few "cells", but this is what was asked of me...

#first use the following command to install in your home directory
#!pip install --user scanpy harmonypy louvain

import scanpy as sc
import anndata as ad
import numpy as np
import scanpy.external as sce
import pandas as pd
import sys
import os

#usage
#python PCA.py -f <file1.genes.results> <file2.genes.results> ... <fileN.genes.results> --oDir <output/directory/of/PCA/chart> --numComps <number of comparissons for PCA, default 2, must be less than sample number. There must be at least 3 files!>
# python PCAtest.py -f /oasis/tscc/scratch/jpburkhardt/RNAseqPipeline/AltSplicing/RNAseqOut/Alpha-Diff.genes.results /oasis/tscc/scratch/jpburkhardt/RNAseqPipeline/AltSplicing/RNAseqOut/Alpha-Prog.genes.results /oasis/tscc/scratch/jpburkhardt/RNAseqPipeline/AltSplicing/RNAseqOut/Hotel-Diff.genes.results /oasis/tscc/scratch/jpburkhardt/RNAseqPipeline/AltSplicing/RNAseqOut/Hotel-Prog.genes.results /oasis/tscc/scratch/jpburkhardt/RNAseqPipeline/AltSplicing/RNAseqOut/Golf-Prog.genes.results /oasis/tscc/scratch/jpburkhardt/RNAseqPipeline/AltSplicing/RNAseqOut/Golf-Diff.genes.results --oDir PCAtest6/PCA --numComps 2

fileFlag=False
ncFlag=False
outFlag=False
oDir=""
nComps=2
files=[]
for arg in sys.argv:
  if arg[0] == '-':
    if arg[1] == '-':
      if arg[2] == 'n':
        if arg[3] == 'u':
          #numComps selected
          fileFlag=False
          ncFlag=True
          outFlag=False
      elif arg[2] == "o":
          #odir sleected
          fileFlag=False
          ncFlag=False
          outFlag=True
    elif arg[1] == 'f':
      #files selected
      fileFlag=True
      ncFlag=False
      outFlag=False
  elif fileFlag:
    files.append(arg)
    print("appending to input for PCA analysis " + str(arg))
  elif ncFlag:
    nComps=int(arg)
    print("setting number of comparissons for PCA analysis to " + str(nComps))
  elif outFlag:
    oDir=arg
    print("setting out directory value for PCA analysis to " + str(oDir))
dataFrames=[]
for file in files:
  file1 = open(file, 'r')
  lines = file1.readlines()

  #the following line creates a dataframe which takes the first column ('gene_id') as the row names
  #column 0 contains gene ids and column 4 contains expected counts of the gene if the hypothesis is correct.
  df = pd.read_csv(file,sep='\t',index_col = 0, header = 0, usecols=[0,4])

  #rename column to be of the dataset
  index=0
  marker=0
  for char in file:
      index+=1
      if char == '/':
          marker=index
  col_name = file[marker:(len(file) - 14)]
  df.rename(columns = {'expected_count':col_name}, inplace = True)
  dataFrames.append(df)


  #print(df.head(10))
dataFrame = dataFrames[0]
for i in range(len(dataFrames) - 1 ):
  dataFrame = pd.merge(dataFrame,dataFrames[i + 1],on=['gene_id'], how='outer').fillna(0)
transposeFrame = dataFrame.T
transposeFrame = transposeFrame.rename(columns={transposeFrame.columns[0]:'datasets'})
print(transposeFrame.head(10))

adata = sc.AnnData(transposeFrame)

sc.pp.pca(adata, n_comps=nComps, zero_center=True, svd_solver='arpack', random_state=0, return_info=False, use_highly_variable=None, dtype='float32', copy=False, chunked=False, chunk_size=None)

#warning, this oDir thing needs to be done to all files with output or else they may attempt to save in TSCC home directory!
imgName=""
if(oDir == ""):
  if not os.path.exists(oDir +  "figures/pca"):
      os.makedirs(oDir +  "figures/pca")
  imgName= "/PCA.png"
else:
  if not os.path.exists(oDir +  "/figures/pca/"):
      os.makedirs(oDir +  "/figures/pca/")
  imgName="/PCA.png"
if not os.path.exists("figures/pca"):
    os.makedirs("figures/pca")
size = 6000/len(files)
sc.pl.pca(adata,save=imgName,size=size,color="datasets",legend_fontsize=10,legend_loc='on data',colorbar_loc=None)
cmd = "mv figures/pca/PCA.png " + oDir + "/figures/pca/"
os.system(cmd)
