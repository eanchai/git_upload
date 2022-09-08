"""
    date: 2022.09.08
    object: 연구용 프로세스 :
        NGS data/fastq 파일의 fastgc 파일에서 각각 종류 분리 및 quality 확인
    module or def: read_file, make_plot function
    author: 
"""
#module
import os
import glob
from selectors import EpollSelector
import pandas as pd
import numpy as np
import matplotlib as mat
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from pathlib import Path
import itertools
mat.rcParams['axes.unicode_minus'] = False

#define functions
def read_file(file_path: Path) -> pd.DataFrame:
#여기에서 연습으로 경로 설정
#fastq_dir = Path("/ess/NGS_data/fastq")
#fastqc_gen = fastq_dir.glob("*/fastqc/*.csv")
    fastq_dir = Path("/home/ubuntu/NGS_data/fastq")
    fastqc_gen=fastq_dir.glob("*/fastqc/*.csv")
    fastqc_gen_list=list(fastqc_gen)
    
    qc_list=[]
    for fn in fastqc_gen_list:
        #fastqc파일에 date추가  
        date= fn.parents[1].name.split("_")[0]
        tmp_df=pd.read_csv(fn)
        tmp_df['date']=date                        
        qc_list.append(tmp_df)
    
    qc_df = pd.concat(qc_list, ignore_index=True)
    qc_df = qc_df.sort_values(by='date')
    
    #SCAN 정보 추가
    qc_df["SCAN"] = [row['SampleID'].rsplit("_",2)[1] for _, row in qc_df.iterrows()]
    qc_df["SCAN"] = qc_df["SCAN"].str.replace('XT',"NaN")
    #TF 정보 저장
    qc_df["Filter"] = [(float(row['FASTQ_Q30_R1']) > 90) & (float(row['FASTQ_Q30_R2']) > 90) & (int(row['FASTQ_TOTAL_BASES(Gb)']) > 24) for _, row in qc_df.iterrows()]
    return qc_df

#plot define
def make_bar_plot(qc_df, target_cols: list, hue_col : str = 'SCAN') -> None:
    target_len = len(target_cols)
    fig_s, ax = plt.subplots(ncols=target_len, figsize=(6*target_len,6))
    for i, col in enumerate(target_cols):                               #enumerate쓰면 index가 같이나오는 tuple나옴
        cutoff = 24 if col == ("FASTQ_TOTAL_BASES(Gb)") else 95         #한줄
        sns.barplot(data=qc_df,
                    x="date",
                    y=col,
                    alpha=0.5,
                    palette='rocket_r',
                    hue=hue_col,
                    ax=ax[i])
        ax[i].axhline(cutoff)
    plt.tight_layout()
    return plt.savefig('Fastqc_barplot.png')

file_dir = Path(".")
qc_df = read_file(file_dir)
make_bar_plot(qc_df, ['FASTQ_Q30_R1', 'FASTQ_Q30_R2', 'FASTQ_TOTAL_BASES(Gb)'])
