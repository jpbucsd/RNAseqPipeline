import pandas as pd
import numpy as np
import sys
import os
from matplotlib import pyplot as plt

#AS_events_sig_spreadsheet
#Prog_vs_Diff.csv

log2 = pd.read_csv("Prog_vs_Diff.csv")
AS = pd.read_excel("AS_events_sig_spreadsheet.xlsx")

AS['log2FoldChange'] = 0

for index, row in AS.iterrows():
  for index2, row2 in log2.iterrows():
    if row["geneSymbol"].strip() == row2["gene_id"].strip():
      row['log2FoldChange'] = row2["gene_id"]

display(AS)
