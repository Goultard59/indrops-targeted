# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 13:59:20 2018

@author: Sophie
"""
import pandas as pd # to deal with dataframe, csv files (groupby function etc...)
import numpy as np # to deal with arrays
import matplotlib.pyplot as plt # to plot graphs
from collections import defaultdict # for Kevin's plots
import os
import argparse

#on verifie si le fichier existe
def isValidFile(parser, arg):
    if not os.path.exists(arg):
        #si il n existe pas afficher une erreur
        parser.error("The file %s does not exist!" % arg)
    else:
        #sinon renvoyer le chemin
        return arg

#on recupere les options utilisateurs 
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--csv", required=True,
                    help="CSV file", type=lambda x: isValidFile(parser, x))
parser.add_argument("-b", "--barcode", required=True,
                    help="barcode file", type=lambda x: isValidFile(parser, x))
parser.add_argument("-f", "--figure", required=True,
                    help="figure file")
#recupere le chemin du fichier de sortie
parser.add_argument("-o", "--outputfile", required=True,
                    help="Output file path")
param = parser.parse_args()

df = pd.read_csv(param.csv)
dfBC = df.groupby(["BC"]).size().reset_index(name="nb_reads")
dfRU = df.groupby(["BC","gene","umi"]).size().reset_index(name="reads_per_umi")
umi = dfRU.groupby(["BC","gene"]).size().reset_index(name="nb_umi")

# FILTER 1 = on the number of reads per BC according to barcode distribution plot AT FO (Keke code / count of the nb of reads/BC)

# Begin by Filter = 0 to have the plots without any filter and then apply a filter according to the plots and re-plot
# ffBC = dfBC with filter on the min nb of reads per BC
    
ffBC = dfBC[(dfBC.nb_reads>=10000)&(dfBC.nb_reads<=100000)]
xfBC = dfBC[(dfBC.nb_reads>=50)&(dfBC.nb_reads<=1000)]
yfBC = dfBC[(dfBC.nb_reads>=10000)&(dfBC.nb_reads<=700000)]

BCRU = dfRU.merge(ffBC, right_on = 'BC', left_on = 'BC')
umiBC = umi.merge(ffBC, right_on = 'BC', left_on = 'BC')

print("number of BC before filtering BC",len(dfBC.index))
print("number of BC after filtering bc 1",len(ffBC.index))
print("number of BC after filtering bc 2",len(xfBC.index))
print("number of BC after filtering bc 3",len(yfBC.index))

#### FILTER 2 = on the nb of reads per umi
#    
ffRU = dfRU[dfRU.reads_per_umi>0]
#    
RUBC=ffRU.merge(dfBC,right_on='BC',left_on='BC')
RUBC = RUBC.groupby(['BC','nb_reads']).size().reset_index(name='nb_umi')
#    
ffumi = ffRU.groupby(["BC","gene"]).size().reset_index(name="nb_umi")
Nb_BC = ffumi.groupby(["BC"]).size().reset_index(name="nb_umi/BC")
#
print("number of BC before filtering umi",len(dfBC.index))
print("number of BC after filtering umi",len(Nb_BC.index))
#

#plot reads/BC 

## Plot keke hist nb of reads per BC after filtering Filter 2 (nb of reads per umi)

## No Filter, use dfBC and drop dfBC.columns[0]
## For FILTER 1 (on nb reads/BC only), use ffBC and drop ffBC.columns[0]
## For FILTER 2 (on nb reads/umi only), use RUBC and drop RUBC.columns[0,2]

  
k2 = dfBC.drop(dfBC.columns[[0]],axis=1)
print(k2)
k2.to_csv(param.outputfile,header=False)
i=0
with open(param.outputfile) as f:
    count_freq = defaultdict(int)
    for line in f:
        count = line.rstrip('\n').split(',')
        count_freq[count[1]] += 1
        i+=1 
    x = np.array(list(count_freq.keys()))
    x = x.astype(int)
    y = np.array(list(count_freq.values()))
    w = x*y

plt.hist(x, np.logspace(0,7,100) , weights=w)
plt.xscale('log')
plt.xlabel(param.barcode)
plt.ylabel('count x #reads/BC')

plt.savefig(param.figure)
    
plt.show()
plt.close()