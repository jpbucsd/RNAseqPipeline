import pandas as pd
import numpy as np
import sys
import os
from bioinfokit import analys, visuz

#usuage, takes CSV files outputed by DifferentialAnalysis.R as arguments
# python heatmap.py --g genelist.txt (not yet implemented) -f Something_vs_something.csv ... -c (for clustering)
print("running")
fileFlag=False
hierarchicalFlag=False
outFlag=False
gFlag=False
gList=False
files=[]
geneList=""
for arg in sys.argv:
  if arg[0] == '-':
    if arg[1] == '-':
      if arg[2] == "g":
        fileFlag=False
        outFlag=False
        gFlag=True
        gList=True
    elif arg[1] == 'f':
      #files selected
      fileFlag=True
      outFlag=False
      gFlag=False
    elif arg[1] == 'c':
      hierarchicalFlag=True
  elif fileFlag:
    files.append(arg)
  elif outFlag:
    oDir=arg
    print("setting odir value to " + str(oDir))
  elif gFlag:
    geneList=arg
    print("setting gene list value to " + str(geneList))

heat =  pd.read_csv(files[0])

if gList:
    lHeat = pd.DataFrame()

    f = open(geneList, "r")
    for x in f:
      for index, row in heat.iterrows():
        if row["Gene"].strip() == x.strip():
          lHeat.loc[len(lHeat.index)] = row
    heat = lHeat

heat = heat.set_index(heat.columns[0])
print(heat)      
if hierarchicalFlag:
  visuz.gene_exp.hmap(df=heat, dim=(6, 12), tickfont=(6, 4))
  #saves as heatmap.png heatmap_clus.png
else:
  visuz.gene_exp.hmap(df=heat, rowclus=False, colclus=False, dim=(6, 12), tickfont=(6, 4))     
