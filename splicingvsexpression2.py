import pandas as pd
import numpy as np
import sys
import os
from matplotlib import pyplot as plt

#AS_events_sig_spreadsheet
#Prog_vs_Diff.csv

log2 = pd.read_csv("Prog_vs_Diff.csv")
AS = pd.read_csv("AS_foldChange.csv")
AS2 = pd.DataFrame();
AS2 = pd.DataFrame(data=AS2, columns=AS.columns)

for index, row in AS.iterrows():
  if row['log2FoldChange'] == 0.0:
    AS.at[index,'log2FoldChange'] = 0.0
  else:
    AS2 = AS2.append(row)

print(AS2)

fig = plt.figure()
fig.set_size_inches((10,10))

plt.scatter(AS2["IncLevelDifference"], AS2['log2FoldChange'],alpha=1,s=5,color="grey")

ASsmall = AS2[AS2['log2FoldChange'] < 1]
ASrange = ASsmall[ASsmall['log2FoldChange'] > -1]

ASL = ASrange.nlargest(10, 'IncLevelDifference')
ASS = ASrange.nsmallest(10, 'IncLevelDifference')

plt.scatter(ASrange["IncLevelDifference"],ASrange["log2FoldChange"],alpha=1,s=5,color="black")

plots = []
for index,row in ASL.iterrows():
  red = plt.scatter(ASL.at[index,"IncLevelDifference"],ASL.at[index,"log2FoldChange"],alpha=1,s=5,color="red",label=ASL.at[index,"geneSymbol"])
  plots.append(red)

for index,row in ASS.iterrows():
  blue = plt.scatter(ASS.at[index,"IncLevelDifference"],ASS.at[index,"log2FoldChange"],alpha=1,s=5,color="blue",label=ASS.at[index,"geneSymbol"])
  plots.append(blue)

plt.axhline(y = 1, color = 'black', label = 'upper bound')
plt.axhline(y = -1, color = 'black', label = 'lower bound')

plt.axvline(x = 0.1, color = 'grey', label = 'x upper bound')
plt.axvline(x = -0.1, color = 'grey', label = 'x lower bound')

plt.xlabel("IncLevelDifference")
plt.ylabel("log2FoldChange")
labels=ASL.geneSymbol.tolist()
labels2=ASS.geneSymbol.tolist()
labels3=labels+labels2
fig.legend(handles=plots)

fig.savefig("AS_DEG.png")
