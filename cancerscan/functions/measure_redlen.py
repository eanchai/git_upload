"""
    date: 2022.10.25
    object: CancerSCAN R to Py
            MeasuerRealen & Make information, histogram file 
    author: rlo
"""

import os
import logging
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from pytz import timezone
from datetime import datetime, timedelta

class MeasureReadlen:
    def __init__(self, samtools_path: Path, input_bam: str, target_date: str, outdir: str):
        '''
            Intialize
            Intput: 

        '''
        self.samtools_path = samtools_path 
        self.input_bam = input_bam
        self.sample_name = input_bam.split(".")[0]

        #log init
        self.outdir = outdir
        # date time setting.
        input_date = datetime.strptime(target_date, "%Y-%m-%d")
        before_one = input_date - timedelta(days=3)
        self.date = before_one.strftime("%Y-%m-%d")
        # self.date = target_date
        log_dir = os.path.join(outdir, "log")
        os.makedirs(log_dir, exist_ok=True, mode=0o777)
        self.logger = self.make_logger(log_dir, f"make_readlen{self.date}")
        
    @staticmethod
    def make_logger(log_dir, name=None, consoleset=True, streamset=True) -> logging:
        '''
            logger making.
        '''
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)
        loggerformat = "%(asctime)s - %(name)s - %(module)s - %(levelname)s - %(message)s"
        formatter = logging.Formatter(loggerformat)
        # change time zone
        def timetz(*args):
            return datetime.now(tz).timetuple() 
        tz = timezone("Asia/Seoul")
        logging.Formatter.converter = timetz
        if consoleset:
            console = logging.StreamHandler()
            console.setLevel(logging.INFO)
            console.setFormatter(formatter)
            logger.addHandler(console)
        if streamset:
            loggerfile = os.path.join(log_dir, name)
            file_handler = logging.FileHandler(filename=f"{loggerfile}.log")
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        return logger

    def use_samtools(self, view_option: str, thread: int, cut_option: str, output_option: str):
        os.system(f"{self.samtools_path} view {view_option} -@ {thread} {self.input_bam} | cut {cut_option} > {self.sample_name}{output_option}")
        pass

    def carculate_media(self, sample_name_tmp0: Path) -> int:
        tmp0_readlen = np.median(len(f'{sample_name_tmp0}'))
        return tmp0_readlen

    def make_tmp1_df(self, sample_name_tmp1: Path) :
        tmp1_df = pd.read_table(f'{sample_name_tmp1}', header=None, sep="\t")[:-1]
        tmp1_df[1] = tmp1_df[1].abs().astype(int)
        tmp1_df = tmp1_df[~tmp1_df[0].isin(['chrM'])]

        masks = tmp1_df[1] < 1500
        masks = masks.tolist()
        tmp1_df = tmp1_df[masks]
        tmp1_chr = np.unique(tmp1_df[0])

        tmp1_median = np.median(tmp1_df[1])
        return tmp1_df, tmp1_median
    
    def make_hist(self, tmp1_df: pd.DataFrame, tmp1_median: int):
        plt.hist(
            x = tmp1_df[1], 
            bins=np.linspace(tmp1_df[1].min(), tmp1_df[1].max(), 1500),
            color='red'
        )
        # plt.title(f"Distribution of insert size(read length:{read_len}")
        plt.xlabel("Insert size")
        plt.ylabel("Frequency")
        plt.xlim(0,1500)
        plt.axvline(tmp1_median, color='white')
        plt.savefig(f"{self.sample_name}.pdf")
        
    def make_info(self, tmp1_median, tmp0_readlen):
        tmp1_df = pd.DataFrame(tmp1_median, tmp0_readlen)

        f = open(f"{self.sample_name}.info", "w")
        f.write(tmp1_df, sep="\t", quote=False)
        f.close()

    def __call__(self):
        self.logger.info("Make tmp0 file...")
        # self.use_samtools('-f 2', 1, '-f10', '.tmp')
        os.system(f"{self.samtools_path} view -f 2 -@ 1 {self.input_bam} | cut -f10 | head > {self.sample_name}'.tmp0'")

        
        self.logger.info("Carculate readlen...")
        tmp1_readlen = self.carculate_media()

        self.logger.info("Make tmp1 file...")
        # self.use_samtools('-f 2', 1, '-f3,9', '.tmp1')
        os.system(f"{self.samtools_path} view -f 2 -@ 1 {self.input_bam} | cut -f3,9 | head > {self.sample_name}'.tmp0'")

        tmp1_df, tmp1_median = self.make_tmp1_df()

        self.logger.info("Make histogram...")
        self.make_hist(tmp1_df, tmp1_median)
        self.logger.info("Make information file...")
        self.make_info(tmp1_median, tmp1_readlen)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--samtools_path", defualt='/usr/bin/samtools')
    parser.add_argument("--input_bam")
    args = parser.parse_args()


tmp1_df = pd.read_csv('/mnt/c/Users/dlgkr/op/git_upload/cancerscan/sample/CD_22_13185_BD_D_SCN_1.tmp1', sep='\t', header=None)
tmp1_df = tmp1_df[~tmp1_df[0].isin(['chrM'])]
tmp1_chr = np.unique(tmp1_df[0])
type(tmp1_chr)

input_bam = 'CD_22_13185_BD_D_SCN_1.tmp1'
input_bam.split(".")[0]

os.system