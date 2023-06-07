#!/usr/bin/env python

"""
:Author: Lea Burkard
:Date: 01.11.2022
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.signal import savgol_filter

import matplotlib
import seaborn as sns




# variables from snakefile
# = snakemake.input["score"]
#score_refs = snakemake.input["score_ref"]
#score_back = snakemake.input["score_back"]
#score_back_refs = snakemake.input["score_back_ref"]

#sample_ID = snakemake.params["sample"]
#ref_IDs = snakemake.params["ref_IDs"]
#target = snakemake.params["target"]
outfile = snakemake.output[0]
overlay_mode = snakemake.params["overlay_mode"]
smoothing = snakemake.params["smoothing"]
rolling = snakemake.params["rolling"]
background_norm = snakemake.params["background_norm"]
#flank_edge = 500
edge_norm = snakemake.params["edge_norm"]
#bin_size = 1
#mean_norm = True

win_len=snakemake.params["win_len"] #165 #from config
poly=snakemake.params["poly"] #3 #from config



# functions
def calculate_flanking_regions(val: int):
    """Calculates flanking regions for point of interest.

    Args:
        val (int): should be length of value vector

    Raises:
        TypeError: Only integers are allowed

    Returns:
        [iterator]: range of values around center point (e.g. range(-1000,1000))
    """

    if not isinstance(val, int):
        raise TypeError("Only integers are allowed")

    if val % 2 == 0:
        flank = int(val / 2)
        region = range(-flank, flank)
    elif val % 2 == 1:
        flank_l = int(val / 2 - 0.5)
        flank_r = int(val / 2 + 0.5)
        region = range(-flank_l, flank_r)
    return region

def add_sample(path_a: str, path_b: str, overlay_mode:str = "mean",smoothing:bool = False, rolling:bool = False, background_norm:bool = False, window:int = 1000, flank=500, mean_norm:bool = False):
    """Reads .csv file, calculates mean over all rows and divides by trimmed mean.

    Args:
        path_a (str): [path to targets]
        path_b (str): [path to background]
        overlay_mode (str): [Mode of operation, e.g. mean, median]

    Returns:
        [Pandas Object]: [Pandas object containing normalized/processed data]
    """

    if overlay_mode.lower() == "mean":
        if mean_norm:
            sample_a = pd.read_csv(path_a, header=None).iloc[:,1:]#.mean(numeric_only=True)
            mean_value = np.nanmean(sample_a.astype(float).values)
            sample_a = sample_a / mean_value
            sample_a = sample_a.mean(numeric_only=True)
        else:
            sample_a = pd.read_csv(path_a, header=None).iloc[:,1:].mean(numeric_only=True)
        # if background_norm:
        #     sample_b = pd.read_csv(path_b, header=None).iloc[:,1:].mean(numeric_only=True, axis=1)
        #     sample = pd.DataFrame(sample_a / stats.trim_mean(sample_b, 0.1), columns=["value"])
        else:
            sample = pd.DataFrame(sample_a, columns=["value"])
        sample["position"] = calculate_flanking_regions(len(sample))
        sample=sample.set_index("position")
    elif overlay_mode.lower() == "median":
        sample_a = pd.read_csv(path_a, header=None).iloc[:,1:].median(numeric_only=True)
        # if background_norm:
        #     sample_b = pd.read_csv(path_b, header=None).iloc[:,1:].median(numeric_only=True, axis=1)
        #     sample = pd.DataFrame(sample_a / stats.trim_mean(sample_b, 0.1), columns=["value"])
        else:
            sample = pd.DataFrame(sample_a, columns=["value"])
        sample["position"] = calculate_flanking_regions(len(sample))
        sample=sample.set_index("position")
    
    else:
        raise ValueError(f"{overlay_mode} is not a valid keyword.")

    if smoothing:
        if rolling:
            sample = sample.apply(lambda x:savgol_filter(x,window_length=win_len, polyorder=poly)) - sample.apply(lambda x:savgol_filter(x,window_length=win_len, polyorder=poly)).rolling(1000, center=True).median() 
        else:
            sample = sample.apply(lambda x:savgol_filter(x,window_length=win_len, polyorder=poly))
    else:
        if rolling:
            sample = sample - sample.rolling(1000, center=True).median(numeric_only=True)

    return sample



import glob, os.path

c = 0
features = pd.DataFrame()
#TF = "LYL1"

for file in glob.glob("../../cfDNA/results/intermediate/"+run+"/table/GRCh38/target/*_MIDPOINT.csv.gz"):
    sample_id = file.split("--")[1].split(".")[0]
    tf_id = mine_GC.split("--")[0].split("/")[-1]
    
    # sample = pd.DataFrame(add_sample(WPS, WPS_back,overlay_mode,smoothing,rolling,background_norm,mean_norm=mean_norm))
    # edge norm?
    
    if overlay_mode.lower() == "mean":
        if mean_norm:
            sample_a = pd.read_csv(path_a, header=None).iloc[:,1:]#.mean(numeric_only=True)
            mean_value = np.nanmean(sample_a.astype(float).values)
            sample_a = sample_a / mean_value
            sample_a = sample_a.mean(numeric_only=True)
        else:
            sample_a = pd.read_csv(path_a, header=None).iloc[:,1:].mean(numeric_only=True)
        # if background_norm:
        #     sample_b = pd.read_csv(path_b, header=None).iloc[:,1:].mean(numeric_only=True, axis=1)
        #     sample = pd.DataFrame(sample_a / stats.trim_mean(sample_b, 0.1), columns=["value"])
        else:
            sample = pd.DataFrame(sample_a, columns=["value"])
        sample["position"] = calculate_flanking_regions(len(sample))
        sample=sample.set_index("position")
    elif overlay_mode.lower() == "median":
        sample_a = pd.read_csv(path_a, header=None).iloc[:,1:].median(numeric_only=True)
        # if background_norm:
        #     sample_b = pd.read_csv(path_b, header=None).iloc[:,1:].median(numeric_only=True, axis=1)
        #     sample = pd.DataFrame(sample_a / stats.trim_mean(sample_b, 0.1), columns=["value"])
        else:
            sample = pd.DataFrame(sample_a, columns=["value"])
        sample["position"] = calculate_flanking_regions(len(sample))
        sample=sample.set_index("position")

    if smoothing:
        if rolling:
            sample = sample.apply(lambda x:savgol_filter(x,window_length=win_len, polyorder=poly)) - sample.apply(lambda x:savgol_filter(x,window_length=win_len, polyorder=poly)).rolling(1000, center=True).median() 
        else:
            sample = sample.apply(lambda x:savgol_filter(x,window_length=win_len, polyorder=poly))
    else:
        if rolling:
            sample = sample - sample.rolling(1000, center=True).median(numeric_only=True)
    
    if c==0:
        midpoint_coverage = pd.DataFrame(sample)
        midpoint_coverage.columns = midpoint_coverage.columns.astype(str)
        midpoint_coverage.columns.values[-1] = sample_id
    else:
        midpoint_coverage[sample_id] = sample
    c = c+1

print(c,midpoint_coverage)

for col in midpoint_coverage:
    print(col)
    ### First feature: central coverage (central 30 bp)
    data = midpoint_coverage[col].reset_index(drop = False)
    data['position'] = data['position'].astype(int)
    first_feature = data[(data["position"] >= -30) & (data["position"] <= 30)]
    #print(first_feature)

    ### Second feature: mean coverage (central 1000 bp)
    second_feature = data[(data["position"] >= -1000) & (data["position"] <= 1000)]
    #print(second_feature)

    ### Third feature: maximum peak amplitude (using FFT)
    third_feature = data[(data["position"] >= -960) & (data["position"] <= 960)]
    fft_index = 10
    
    for col in first_feature.columns:
        if col != "position":
            #print(col,first_feature[col].mean())
            features.loc[col,"central_coverage_"+tf_id] = first_feature[col].mean(skipna=False)
            #print(col,second_feature[col].mean())
            features.loc[col,"mean_coverage_"+tf_id] = second_feature[col].mean(skipna=False)
            fft_res = np.fft.fft(third_feature[col])
            #print(col,"fft",np.abs(fft_res[fft_index]))
            features.loc[col,"amplitude_"+tf_id] = np.abs(fft_res[fft_index])

#print(features)
features.to_csv(outfile, sep='\t',)