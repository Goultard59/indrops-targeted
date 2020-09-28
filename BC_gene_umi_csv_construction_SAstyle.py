# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 14:58:40 2018

@author: Sophie
"""

import sys
import os
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator

#check path to file
def isValidFile(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist !" % arg)
    else:
        return arg

###
# argparse
###
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sam", required=True,
                    help="SAM file", type=lambda x: isValidFile(parser, x))
parser.add_argument("-r", "--read2", required=True,
                    help="RD2 file", type=lambda x: isValidFile(parser, x))
parser.add_argument("-f", "--fastqread2", required=True,
                    help="fastq RD2 file", type=lambda x: isValidFile(parser, x))

parser.add_argument("-o", "--outputfile", required=True,
                    help="Output file path")
param = parser.parse_args()

def find_umi(seq):
    ##this function will look for the umi inside (seq) at a given position (i.e before the linker)
    Link="TACGCTACG"
    if Link in seq:
        p = seq.find(Link)#link sequence // p = position of the linker in the read
        return seq[p+15:p+23]
    else:
        return "-"

#initiate reads list
ReadList = {}

###
# aligning part
###
count = 0
count_quality = 0
with open(param.sam, "r") as sm:
    for lb in sm:
        #divise the line
        lb = lb.strip("\n").split("\t")
        # trash reads were gene doesn't align
        if int(lb[4]) < 20 or lb[2] == "*":
            count += 1
            count_quality += 1 
        else :
            #create dictionnary entry and add gene
            count += 1
            key = lb[0]
            ReadList.setdefault(key, [])
            ReadList[lb[0]].append(lb[2])
print("Total reads:",count)
print("reads with gene aligned with MQ<20 is:",count_quality)
print("Part 1 done")

###
# Barcode part
###
count_barcode = 0
with open(param.read2, "r") as rd2:
    for la in rd2:
        #divise the line
        la = la.strip("\n").split(",")
        if la[0] in ReadList.keys():
            #incomplete barcode
            if "*" in la[3]:
                count_barcode += 1
                del ReadList[la[0]]
            else:
                #add barcode to dictionnary
                ReadList[la[0]].append(la[3])
print("BC error:",count_barcode)
print("Part 2 done")

###
# UMI part
###
count_umi = 0
count_survive = 0
with open(param.fastqread2) as in_handle:
    #fastq reader
    for title, seq, qual in FastqGeneralIterator(in_handle):
        title = title.split(" ")
        if title[0] in ReadList.keys():
            umi = find_umi(seq)
            #umi error
            if "-" in umi:
                count_umi += 1
                del ReadList[title[0]]
            else:
                #add umi to dictionnary
                ReadList[title[0]].append(umi)
                count_survive += 1
print("no umi :",count_umi)
print("number of surviving reads :",count_survive)
print("Part 3 done")

#creating the output CSV file
try:
    with open(param.outputfile, "w", newline="") as csv_file:
        csv_file.write("ID,gene,BC,umi\n")
        for key,a_list in sorted(ReadList.items(),key=lambda row:row[0]):
            x = [key,]
            x.extend(a_list)
            csv_file.write(",".join(str(i) for i in x) + "\n")
except IOError:
    print("I/O error")
