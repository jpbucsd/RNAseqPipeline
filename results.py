import pandas as pd
import numpy as np
import sys
import os
from matplotlib import pyplot as plt

#usuage, takes CSV files outputed by DifferentialAnalysis.R as arguments
# results.py -f file1 ... fileN --padjusted 0.05 -- 
#--padjusted or --padj - adjusted pvalue, a padj of 0.5 implies 50% of significant results are false positives. results with Padj above 50% are filtered out
#--log10 or --psig - -log10(pvalue) to determine which results are significant
#--Llog10 to determine which results are worth naming, a log 10 value as threshold for which values to name
#--odir output directory
fileFlag=False
padFlag=False
logFlag=False
llogFlag=False
outFlag=False
oDir=""
padJ=0.5
pval=5
lval=30
files=[]
for arg in sys.argv:
  if arg[0] == '-':
    if arg[1] == '-':
      if arg[2] == 'p':
        if arg[3] == 'a':
          #padj selected
          fileFlag=False
          padFlag=True
          logFlag=False
          llogFlag=False
          outFlag=False
        elif arg[3] == 's':
          #log10 selected
          fileFlag=False
          padFlag=False
          logFlag=True
          llogFlag=False
          outFlag=False
      elif arg[2] == 'l':
          #log10 selected
          fileFlag=False
          padFlag=False
          logFlag=True
          llogFlag=False
          outFlag=False
      elif arg[2] == 'L':
          #log10 selected
          fileFlag=False
          padFlag=False
          logFlag=False
          llogFlag=True
          outFlag=False
      elif arg[2] == "o":
          #odir sleected
          fileFlag=False
          padFlag=False
          logFlag=False
          llogFlag=False
          outFlag=True
    elif arg[1] == 'f':
      #files selected
      fileFlag=True
      padFlag=False
      logFlag=False
      llogFlag=False
      outFlag=False
  elif fileFlag:
    files.append(arg)
  elif padFlag:
    padJ=float(arg)
    print("setting padj value to " + str(padJ))
  elif logFlag:
    pval=float(arg)
    print("setting pval value to " + str(pval))
  elif llogFlag:
    lval=float(arg)
    print("setting label value to " + str(lval))
  elif outFlag:
    oDir=arg
    print("setting odir value to " + str(oDir))
for file in files:
    #a dataframe is made from the deseq csv file. results with no pvalue are cleaved, the first column is named gene ID, and results with padj below < 0.5 are selected for.
    #a dataframe deresultsSig is made from files with a -log10pval above 5, to filter for significant data.
    #using deresultsSig, two dataframes are made, those with positive and negative log2FoldChanges to represent those that are down or up regulated.

    deresults = pd.read_csv(file)
    deresults = deresults[-np.isnan(deresults["pvalue"])]
    deresults = deresults.rename({"Unnamed: 0": "gene_id"}, axis='columns')
    #deresults = deresults[deresults["padj"] < padJ ] #should this step be included now or after graphing?, for now we are moving to after graphing
    
    #deresults = deresults[-1*np.log10(deresults["pvalue"]) >= pval]
    deresultsSig = deresults[-1*np.log10(deresults["pvalue"]) >= pval]
    deresultsNSig = deresults[-1*np.log10(deresults["pvalue"]) < pval]
    deresultsESig = deresults[-1*np.log10(deresults["pvalue"]) >= lval]
    
    fig = plt.figure()
    fig.set_size_inches((10,10))
    volcano = fig.add_subplot(111)
    volcano.scatter(deresultsNSig["log2FoldChange"],-1*np.log10(deresultsNSig["pvalue"]),s=2,alpha=0.5,color="blue")
    volcano.scatter(deresultsSig["log2FoldChange"],-1*np.log10(deresultsSig["pvalue"]),s=2,alpha=0.5,color="red")
    volcano.set_xlabel("Log2FoldChange")
    volcano.set_ylabel("Log10 pvalue")
    volcano.axvline(x=0.0, linestyle="dashed", color="grey",
         linewidth=1)

    #a for loop uses the deresults3 dataframe to annotate the geneIDs of the most significant genes
    for index, row in deresultsESig.iterrows():
        volcano.text(row["log2FoldChange"],-1*np.log10(row["pvalue"]),row["gene_id"], horizontalalignment='left', size=10, color='black')
        print("x: " + str(row["log2FoldChange"]) + ", y: " + str(-1*np.log10(row["pvalue"])) + " , name: " + str(row["gene_id"]))
    
    
    #save the file!
    nMark=0
    for i, char in enumerate(file):
        if char == '/':
            nMark = i
    figName = file[nMark + 1:-4]
    
    if not os.path.exists(oDir):
      os.makedirs(oDir)
    if not os.path.exists(oDir + "/" + figName):
      os.makedirs(oDir  + "/" + figName)
    
    fig.savefig(oDir + "/" + figName + "/" + figName + ".png")
    
    
    #remove values that are likely false positives, I am unsure if this step should be done after graphing or not, for now it will be
    deresultsSig = deresultsSig[deresults["padj"] < padJ ]
    downRegulated = deresultsSig[deresultsSig["log2FoldChange"] < 0]
    upRegulated = deresultsSig[deresultsSig["log2FoldChange"] > 0]
    
    DRlist = ''
    for index, row in downRegulated.iterrows():
        DRlist += row["gene_id"] + '\n'
    fileD = open(oDir + "/" + figName + "/" + "downRegulatedGenes.txt","w")
    fileD.writelines(DRlist)
    fileD.close()

    URlist = ''
    for index, row in upRegulated.iterrows():
        URlist += row["gene_id"] + '\n'
    fileU = open(oDir + "/" + figName + "/" + "upRegulatedGenes.txt","w")
    fileU.writelines(URlist)
    fileU.close()
