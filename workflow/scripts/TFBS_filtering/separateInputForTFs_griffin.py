#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys, os
import argparse
import gzip
import math

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", dest="input", help="Input Bed File")
parser.add_argument("-o","--outdir", dest="outdir", help="Name of the output directory")
parser.add_argument("-v","--verbose", dest="verbose", help="Print updates",default=False)
parser.add_argument("-s","--sites", dest="sites", help="File containing the valid sites")
parser.add_argument("-n","--number", dest="number", help="Number of binding sites per TF")
parser.add_argument("-r","--regions", dest="regions", help="File with the regions for score workflow")
options = parser.parse_args()

try:
    os.mkdir(options.outdir)
    if options.verbose: print ("New directory created.")
except OSError as error:
    if options.verbose: print ("Directory already exists.")

if os.path.exists(options.input):
    if ".gz" in options.input:
        inputFile = gzip.open(options.input,"r")
    else:
        inputFile = open(options.input)
    if options.verbose: print ("File opened.")


number_sites = int(options.number)
chroms = ['chr'+str(m) for m in range(1,23)]
tf_count = 0

data = pd.read_csv(options.input, sep='\t', header=None, names=['Chrom',"Start",'End','TF','peak.count','Strand'])
if options.verbose: print ("Input data read in.")

data = data[~data['TF'].isnull()] # remove entries without a TF name
data['TF'] = data['TF'].replace('/','-') # replace / character
data = data[data['Chrom'].isin(chroms)]

valid_sites = pd.read_csv(options.sites, sep='\t', usecols=['TF_Name','TF_Status'])
valid_sites = valid_sites[valid_sites['TF_Status']!='N'] # exclude unknown binding sites (status N)
valid_sites = valid_sites['TF_Name'].unique()

out_regions = open(options.regions,'wt')
if options.verbose: print ("Config file created.")
out_regions.write("target\tpath\n")

for tf in data['TF'].unique():
    tf_data = data[data['TF']==tf].sort_values(by='peak.count', ascending=False)
    if len(tf_data) >= number_sites and tf in valid_sites:
        top_sites = tf_data.iloc[0:number_sites].reset_index(drop=True)
        #print(tf,len(tf_data))
        out_file = options.outdir + tf + '.' + str(number_sites) + '.bed.gz'
        position = np.floor((top_sites["End"]+top_sites["Start"])/2).astype(int)
        top_sites["Start"] = position
        top_sites["End"] = (position+1)
        top_sites.to_csv(out_file, sep='\t', index=False, compression='gzip', header=False)
        tf_count = tf_count+1
        out_regions.write(tf+"\t"+out_file+"\n")

if options.verbose: print ("Created " + str(tf_count) + " files.")
out_regions.close()
if options.verbose: print ("Script finished.")