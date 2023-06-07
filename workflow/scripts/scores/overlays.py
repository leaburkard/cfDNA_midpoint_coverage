#!/usr/bin/env python

"""
:Author: Sebastian Röner
:Contact: sebastian.roener@charite.de
:Date: 08.09.2020
"""

import matplotlib
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from scipy.signal import savgol_filter
import numpy as np

matplotlib.use("pdf")


# get variables from snakefile

WPS = snakemake.input["WPS"]
WPS_refs = snakemake.input["WPS_ref"]
COV = snakemake.input["COV"]
COV_refs = snakemake.input["COV_ref"]
WPS_back = snakemake.input["WPS_back"]
WPS_back_refs = snakemake.input["WPS_back_ref"]
COV_back = snakemake.input["COV_back"]
COV_back_refs = snakemake.input["COV_back_ref"]
MP = snakemake.input["MP"]
MP_refs = snakemake.input["MP_ref"]
MP_back = snakemake.input["MP_back"]
MP_back_refs = snakemake.input["MP_back_ref"]

sample_ID = snakemake.params["sample"]
ref_IDs = snakemake.params["ref_IDs"]
target = snakemake.params["target"]
outfile = snakemake.output["plot"]
overlay_mode = snakemake.params["overlay_mode"]
smoothing = snakemake.params["smoothing"]
rolling = snakemake.params["rolling"]
background_norm = snakemake.params["background_norm"]
flank_edge = 500
edge_norm = snakemake.params["edge_norm"]
bin_size = 1
mean_norm = True

win_len=165 #21
poly=3 #2

# def functions


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

    if window % 2 == 0:
        fstart=int(window/2+1)
        fstop=int(-window/2)
    elif window % 2 == 1:
        fstart=int(window/2-0.5+1)
        fstop=int(-window/2+0.5)
    
    if overlay_mode.lower() == "mean":
        if mean_norm:
            sample_a = pd.read_csv(path_a, header=None).iloc[:,1:]#.mean(numeric_only=True)
            mean_value = np.nanmean(sample_a.astype(float).values)
            sample_a = sample_a / mean_value
            sample_a = sample_a.mean(numeric_only=True)
        else:
            sample_a = pd.read_csv(path_a, header=None).iloc[:,1:].mean(numeric_only=True)
        if background_norm:
            sample_b = pd.read_csv(path_b, header=None).iloc[:,1:].mean(numeric_only=True, axis=1)
            sample = pd.DataFrame(sample_a / stats.trim_mean(sample_b, 0.1), columns=["value"])
        else:
            sample = pd.DataFrame(sample_a, columns=["value"])
        sample["position"] = calculate_flanking_regions(len(sample))
        sample=sample.set_index("position")
    elif overlay_mode.lower() == "median":
        sample_a = pd.read_csv(path_a, header=None).iloc[:,1:].median(numeric_only=True)
        if background_norm:
            sample_b = pd.read_csv(path_b, header=None).iloc[:,1:].median(numeric_only=True, axis=1)
            sample = pd.DataFrame(sample_a / stats.trim_mean(sample_b, 0.1), columns=["value"])
        else:
            sample = pd.DataFrame(sample_a, columns=["value"])
        sample["position"] = calculate_flanking_regions(len(sample))
        sample=sample.set_index("position")
    elif overlay_mode.lower() == "confidence":
        sample_a = pd.read_csv(path_a, header=None).iloc[:,1:].T
        if background_norm:
            sample_b = pd.read_csv(path_b, header=None).iloc[:,1:].mean(numeric_only=True, axis=1)
            sample = sample_a / stats.trim_mean(sample_b, 0.1)
        else:
            sample = sample_a
        sample["position"] = calculate_flanking_regions(len(sample))
        sample = sample.set_index("position")
        
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
    
    #sample = sample.iloc[fstart:fstop,:]
        
    if overlay_mode.lower() == "confidence":
        sample = sample.melt(ignore_index=False, var_name="sample_nr")
        
    return sample


# load tables containing position specific scores for all defined target regions
# average over all regions per sample and substract the trimmed mean to normalise

### WPS ###
res_sample = pd.DataFrame(add_sample(WPS, WPS_back,overlay_mode,smoothing,rolling,background_norm,mean_norm=mean_norm))
sys.stderr.write("WPS [%s]: %s %s\n"%(sample_ID,COV, COV_back))

if(edge_norm == True):
    m = pd.concat([res_sample.head(100), res_sample.tail(100)]).mean()
    av_WPS = pd.DataFrame(res_sample/m)
else:
    av_WPS = pd.DataFrame(res_sample)

av_WPS.columns = av_WPS.columns.astype(str)
av_WPS.columns.values[-1] = sample_ID
# feature_data_WPS = av_WPS

for (ref_ID, WPS_ref, WPS_back_ref) in zip(ref_IDs, WPS_refs, WPS_back_refs):
    sys.stderr.write("WPS [%s]: %s %s\n"%(ref_ID, WPS_ref, WPS_back_ref))
    res_ref = add_sample(WPS_ref, WPS_back_ref,overlay_mode,smoothing,rolling,background_norm,mean_norm=mean_norm)["value"]
    if(edge_norm == True):
        m = pd.concat([res_ref.head(100), res_ref.tail(100)]).mean()
        av_WPS[ref_ID] = res_ref/m
    else:
        av_WPS[ref_ID] = res_ref

### COV ###
res_sample = pd.DataFrame(add_sample(COV, COV_back,overlay_mode,smoothing,rolling, background_norm,mean_norm=mean_norm))
sys.stderr.write("COV [%s]: %s %s\n"%(sample_ID, COV, COV_back))

if(edge_norm == True):
    m = pd.concat([res_sample.head(100), res_sample.tail(100)]).mean()
    av_COV = pd.DataFrame(res_sample/m)
else:
    av_COV = pd.DataFrame(res_sample)

av_COV.columns = av_COV.columns.astype(str)
av_COV.columns.values[-1] = sample_ID
# feature_data_COV = av_COV

for (ref_ID, COV_ref, COV_back_ref) in zip(ref_IDs, COV_refs, COV_back_refs):
    sys.stderr.write("COV [%s]: %s %s\n"%(ref_ID, COV_ref, COV_back_ref))
    res_ref = add_sample(COV_ref, COV_back_ref,overlay_mode,smoothing,rolling,background_norm,mean_norm=mean_norm)["value"]
    if(edge_norm == True):
        m = pd.concat([res_ref.head(100), res_ref.tail(100)]).mean()
        av_COV[ref_ID] = res_ref/m
    else:
        av_COV[ref_ID] = res_ref

### MIDPOINT ###
res_sample = pd.DataFrame(add_sample(MP, MP_back,overlay_mode,smoothing,rolling, background_norm,mean_norm=mean_norm))
sys.stderr.write("MP [%s]: %s %s\n"%(sample_ID, MP, MP_back))

if(edge_norm == True):
    m = pd.concat([res_sample.head(100), res_sample.tail(100)]).mean()
    av_MP = pd.DataFrame(res_sample/m)
else:
    av_MP = pd.DataFrame(res_sample)

av_MP.columns = av_MP.columns.astype(str)
av_MP.columns.values[-1] = sample_ID
# feature_data_MP = av_MP

for (ref_ID, MP_ref, MP_back_ref) in zip(ref_IDs, MP_refs, MP_back_refs):
    sys.stderr.write("MP [%s]: %s %s\n"%(ref_ID, MP_ref, MP_back_ref))
    res_ref = add_sample(MP_ref, MP_back_ref,overlay_mode,smoothing,rolling,background_norm,mean_norm=mean_norm)["value"]
    if(edge_norm == True):
        m = pd.concat([res_ref.head(100), res_ref.tail(100)]).mean()
        av_MP[ref_ID] = res_ref/m
    else:
        av_MP[ref_ID] = res_ref

### Bins of midpoints ###
sampleNames = list(av_MP)
indices = av_MP.index.tolist()

position_list = []
data_list = []
count_list = []
for name in sampleNames:
    count_list.append([name,0])

c = 1
for i in indices:
    if c == 1:
        position_list.append(i) #save first position
    for el in count_list:
        el[1] += av_MP.loc[i,el[0]] #add up the midpoints of 15 positions
    if c % bin_size == 0:
        append_list = []
        for el in count_list:
            append_list.append(el[1]/c) #divide by counter
            el[1] = 0
        data_list.append(append_list) #append data
        c = 1
    else:
        c += 1

if c % bin_size != 0:
    append_list = []
    for el in count_list:
        append_list.append(el[1]/c)
        el[1] = 0
    data_list.append(append_list)

av_MP_bins = pd.DataFrame(data_list, index=position_list, columns=sampleNames)

# create line plots and save to a single pdf

if overlay_mode.lower() == "confidence":
    av_WPS_long=av_WPS.reset_index().melt(id_vars=["position", "sample_nr"],  value_name="score", var_name="phenotype").sort_values(by="position")
    av_COV_long=av_COV.reset_index().melt(id_vars=["position", "sample_nr"],  value_name="score", var_name="phenotype").sort_values(by="position")
    av_MP_long=av_MP_bins.reset_index().melt(id_vars=["position", "sample_nr"],  value_name="score", var_name="phenotype").sort_values(by="position")
    Fig_WPS = sns.lineplot(data=av_WPS_long, x="position", y="score", hue="phenotype",)
    plt.suptitle("adjusted WPS: {target} target regions")
    plt.xlabel("Position relative to target site")
    plt.ylabel("normalized WPS")
    plt.close()
    Fig_Cov = sns.lineplot(data=av_COV_long, x="position", y="score", hue="phenotype",)
    plt.suptitle(f"adjusted read coverage: {target} target regions")
    plt.xlabel("Position relative to target site")
    plt.ylabel("normalized read coverage")
    plt.close()
    Fig_MP = sns.lineplot(data=av_MP_long, x="position", y="score", hue="phenotype",)
    plt.suptitle(f"midpoint coverage in {bin_size} bp bins: {target} target regions")
    plt.xlabel("Position relative to target site")
    plt.ylabel("normalized midpoint coverage")
    plt.close()
else:
    Fig_WPS = sns.lineplot(data=av_WPS)
    plt.suptitle(f"adjusted WPS: {target} target regions")
    plt.xlabel("Position relative to target site")
    plt.ylabel("normalized WPS")
    plt.close()
    Fig_Cov = sns.lineplot(data=av_COV)
    plt.suptitle(f"adjusted read coverage: {target} target regions")
    plt.xlabel("Position relative to target site")
    plt.ylabel("normalized read coverage")
    plt.close()
    Fig_MP = sns.lineplot(data=av_MP_bins)
    plt.suptitle(f"midpoint coverage in {bin_size} bp bins: {target} target regions")
    plt.xlabel("Position relative to target site")
    plt.ylabel("normalized midpoint coverage")
    plt.close()

with PdfPages(outfile) as pdf:
    pdf.savefig(Fig_WPS.get_figure())
    pdf.savefig(Fig_Cov.get_figure())
    pdf.savefig(Fig_MP.get_figure())
