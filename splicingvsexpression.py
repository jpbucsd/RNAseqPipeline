import pandas as pd
import numpy as np
import sys
import os
from matplotlib import pyplot as plt

#AS_events_sig_spreadsheet
#Prog_vs_Diff.csv

log2 = pd.read_csv("Prog_vs_Diff.csv")
AS = pd.read_csv("AS_events_sig_spreadsheet.csv")

AS['log2FoldChange'] = 0.0

for index, row in AS.iterrows():
  for index2, row2 in log2.iterrows():
    if row["geneSymbol"].strip() == row2[0].strip():
      AS.at[index,'log2FoldChange'] = row2['log2FoldChange']

print(AS)
AS.to_csv('AS_foldChange.csv')

fig = plt.figure()
fig.set_size_inches((10,10))

plt.scatter(AS["IncLevelDifference"], AS['log2FoldChange'],alpha=0.5,color="grey")

ASsmall = AS[AS['log2FoldChange'] < 0.05]
ASrange = ASsmall[ASsmall['log2FoldChange'] > -0.05]

plt.scatter(ASrange["IncLevelDifference"],ASrange["log2FoldChange"],alpha=1,color="black")
plt.axhline(y = 0.05, color = 'black', label = 'upper bound')
plt.axhline(y = -0.05, color = 'black', label = 'lower bound')

plt.axvline(x = 0.05, color = 'grey', label = 'x upper bound')
plt.axvline(x = -0.05, color = 'grey', label = 'x lower bound')

fig.savefig("ASvsFold.png")
