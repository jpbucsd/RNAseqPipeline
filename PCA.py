#first use the following command to install in your home directory
#!pip install --user scanpy harmonypy louvain

import scanpy as sc
import anndata as ad
import numpy as np
import scanpy.external as sce
import pandas as pd
import sys
import os

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

#usage
#python PCA.py -f <file1.genes.results> <file2.genes.results> ... <fileN.genes.results> --oDir <output/directory/of/PCA/chart> --numComps <number of comparissons for PCA, default is sample number. There must be at least 2 files!>
fileFlag=False
ncFlag=False
outFlag=False
oDir=""
nComps=0
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
#transposeFrame = transposeFrame.rename(columns={transposeFrame.columns[0]:'datasets'})#this line actually names gene1 datasets instead of replacing geneid
print(transposeFrame.head(10))

#create plot

#shape 0 represents the number of samples while shape 1 represents the number of components, however we must use shape 0 to cut the thousands of samples (genes) as we can only perform PCA with the same or less components as samples.
shrunkFrame = transposeFrame.iloc[:, 10:16]
if nComps == 0 or nComps > transposeFrame.shape[0]:
        nComps = transposeFrame.shape[0]
pca = PCA(nComps)
pcaDF = pd.DataFrame(pca.fit_transform(transposeFrame))
pcaDF.index=transposeFrame.index
pcaDF.columns=['PC%s' % _ for _ in range(len(pcaDF.columns))]
print(pcaDF)
print(pca.explained_variance_ratio_)
ax1 = pcaDF.plot.scatter(x='PC0', y='PC1',cmap='viridis')#to add colors use c= list of colors!
for i, label in enumerate(pcaDF.index):
    ax1.annotate(label, (pcaDF.iloc[:,0][i], pcaDF.iloc[:,1][i]))
ax1.set_xlabel("PC1 (" + str(round(pca.explained_variance_ratio_[0]*100,1)) + "%)")
ax1.set_ylabel("PC2 (" + str(round(pca.explained_variance_ratio_[1]*100,1)) + "%)")
ax1.tick_params(axis='x', labelsize=5)
ax1.tick_params(axis='y', labelsize=5)

#warning, this oDir thing needs to be done to all files with output or else they may attempt to save in TSCC home directory!
imgName=""
if(oDir == ""):
  if not os.path.exists(oDir +  "figures/pca"):
      os.makedirs(oDir +  "figures/pca")
  imgName= 'figures/pca/PCA.png'
else:
  if not os.path.exists(oDir +  "/figures/pca/"):
      os.makedirs(oDir +  "/figures/pca/")
  imgName= oDir + '/figures/pca/PCA.png'
ax1.figure.savefig(imgName)
