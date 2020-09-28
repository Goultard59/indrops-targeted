# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 10:46:39 2018

@author: Sophie
"""
import pandas as pd # to deal with dataframe, csv files (groupby function etc...)
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
parser.add_argument("-p", "--polyt", required=True,
                    help="CSV file", type=lambda x: isValidFile(parser, x))
parser.add_argument("--minimum", required=True,
                    help="barcode file", type=int)
parser.add_argument("--maximum", required=True,
                    help="figure file", type=int)
#recupere le chemin du fichier de sortie
parser.add_argument("-o", "--outputfile", required=True,
                    help="Output file path")
parser.add_argument("-m", "--matrice", required=True,
                    help="output matrice file")
param = parser.parse_args()

df = pd.read_csv(param.polyt)
#    
dfBC = df.groupby(["BC"]).size().reset_index(name="nb_reads")
dfRU = df.groupby(["BC","gene","umi"]).size().reset_index(name="reads_per_umi")
umi = dfRU.groupby(["BC","gene"]).size().reset_index(name="nb_umi")
    

### FILTER 1 = on the number of reads per BC according to barcode distribution plot AT FO (Keke code / count of the nb of reads/BC)

### Begin by Filter = 0 to have the plots without any filter and then apply a filter according to the plots and re-plot
#### ffBC = dfBC with filter on the min nb of reads per BC
#    
ffBC = dfBC[(dfBC.nb_reads>=param.minimum) & (dfBC.nb_reads<=param.maximum)]
#
BCRU = dfRU.merge(ffBC, right_on = 'BC', left_on = 'BC')
umiBC = umi.merge(ffBC, right_on = 'BC', left_on = 'BC')

print("number of BC before filtering BC",len(dfBC.index))
print("number of BC after filtering bc",len(ffBC.index))
#    
N1 = dfBC['nb_reads'].sum()
N2 = ffBC['nb_reads'].sum()
print("number of reads before filtering",N1)
print("number of reads after filtering",N2)
##    

# pivot to construc matrice

tab = umiBC.pivot(index="BC", columns="gene", values="nb_umi")
tab.fillna(0.0, inplace=True)
tab.to_csv(param.matrice)
    
tab['tot']=tab.sum(axis=1)
nb_umi=tab['tot'].sum()
print("total number of umi",nb_umi)
print("total number of BC",len(tab.index))
    
gene = BCRU.groupby(["gene"]).size().reset_index(name="nb_umi")
gene['prop']=100*(gene['nb_umi']/gene['nb_umi'].sum())
    
gene.to_csv(param.outputfile)
#   