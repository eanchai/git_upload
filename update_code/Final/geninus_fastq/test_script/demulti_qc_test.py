
# test IndexHoppin

import os
import pandas as pd
import xml.etree.ElementTree as ET

# set variable path

demux_l1_path = "/Users/duaghk/build/iGeniPipe/testfile/Stats/DemuxSummaryF1L1.txt"
demux_l2_path = "/Users/duaghk/build/iGeniPipe/testfile/Stats/DemuxSummaryF1L2.txt"

samplesheet_path = "/Users/duaghk/build/iGeniPipe/testfile/Stats/SampleSheet_20210924.csv"

sample_df = pd.read_csv(samplesheet_path)
sample_df

# check.

def parse_demux_result(demux_path, i):
    returndf = []
    with open(demux_path) as f:
        line = ""
        while not line.startswith("### Columns"):
            line = f.readline()
        for line in f.readlines():
            idx, count = line.strip().split()
            idx1, idx2 = idx.split("+")
            tmpseries = pd.Series(data=[idx1, idx2, count],
                                  index=["Index1", "Index2", f"Count_{i}"])
            returndf.append(tmpseries)
    returndf = pd.DataFrame(returndf)
    return returndf

def parse_multi_demux(demux_path_list):
    returndf = ""
    for i, path in enumerate(demux_path_list):
        tmpdf = parse_demux_result(path, i+1)
        if type(returndf) == pd.DataFrame:
            returndf = pd.merge(returndf, tmpdf, on=["Index1", "Index2"], how="outer")
        else:
            returndf = tmpdf
    return returndf

demux_path_list = [demux_l1_path, demux_l2_path]
df1 = parse_multi_demux(demux_path_list)

# check index hopping.
sample_df

# how can i check index hoppin?
df2 = df1[df1["Index2"].isin(sample_df["Index2"].tolist()) & df1["Index1"].isin(sample_df["Index"].tolist())]




