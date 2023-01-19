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


#make a fixed excel sheet

AS1 = pd.DataFrame();
AS1 = pd.DataFrame(data=AS1,columns=AS.columns)

#we need to make every row in the log2Foldchange the negative version of itself, since we did Prog vs Diff rather than Diff vs Prog
#we also need to check if the 0.0 values are NA or actually 0 for log2fold change, this makes a new dataframe with no NA values

for index, row in AS.iterrows():
  row['log2FoldChange'] = -row['log2FoldChange']
  if row['log2FoldChange'] == 0.0:
    AS.at[index,'log2FoldChange'] = 0.0
    if(log2.loc[log2[0].strip() == row["geneSymbol"].strip()][0]['log2FoldChange'] == 0.0):
      AS1 = AS1.append(row)
      AS2 = AS2.append(row)
  else:
    AS1 = AS1.append(row)
    AS2 = AS2.append(row)

AS1.to_csv('AS_foldChange_Diff_vs_Prog.csv')

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

fig.savefig("AS_DEG.eps",format='eps')
