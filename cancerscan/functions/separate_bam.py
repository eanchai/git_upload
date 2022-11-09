"""
    date: 2022.11.08
    object: CancerSCAN R to Py
            Bam separation
    author: rlo
"""

import os
import logging
import numpy as np
import pandas as pd
import seaborn as sns
import subprocess as sp
import matplotlib.pyplot as plt
from tqdm import tqdm
from time import time
from pathlib import Path
from pytz import timezone
from datetime import datetime, timedelta

class SeparateBam:
    def __init__(self, input_bam: Path, tmp1_chr: np.ndarray, target_date: str, outdir: str):
        '''
            Intialize
            Intput: Bam file
            Output: each chromosome bam file

        '''
        self.samtools_path = ""
        self.input_bam = input_bam
        self.tmp1_chr = tmp1_chr
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
        self.logger = self.make_logger(log_dir, f"make_separatebam{self.date}")
        
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


    def separate_bam(self, tmp1_chr, sample_name):
        for i in tqdm(range(0, len(tmp1_chr))):
            os.system(f"{self.samtools_path} view -bh -F 12 {self.input_bam} {tmp1_chr[i]} > {sample_name}.{tmp1_chr[i]}")

    def __call__(self):
        self.logger.info("Start BAM separation...")
        self.separate_bam()
        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--samtools_path", defualt='/usr/bin/samtools')
    parser.add_argument("--input_bam")
    args = parser.parse_args()