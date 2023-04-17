import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage

csvname = ""
csvflag=False
oflag=False
odir=""

for arg in sys.argv:
  if arg[0] == '-':
    if arg[1] == 'd':
      csvflag=True
    elif arg[1] == 'o':
      oflag=True
  elif csvflag:
    csvname = arg
    csvflag = False
  elif oflag:
    oflag=False
    odir=arg

df = pd.read_csv(csvname)

data = list(df)

hierarchical_cluster = AgglomerativeClustering(n_clusters=2, affinity='euclidean', linkage='ward')
labels = hierarchical_cluster.fit_predict(data)

plt.scatter(x, y, c=labels)
plt.show()
