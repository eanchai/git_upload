"""
    date: 2022.09.14 ~
    object: BAM metric upload
        check BAM file & sql database upload
    author: 
"""

import logging
import pandas as pd
import subprocess as sp
import sqlalchemy as db
from time import time
from pathlib import Path
from pytz import timezone
from datetime import datetime

class BamUpload:
    def __init__(self, data_dir: str, sample_id: str) -> None:
        '''
            Class initialize
        '''
        self.data_dir = Path(data_dir)
        self.sample_id = sample_id
        self.log_dir = self.data_dir.joinpath("log")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self._make_logger(self.log_dir, 'BamUploadTest')
        self.file_ex = ".recal.bam.metrics.txt" #suffix
        self.recal_metric_fn = self.data_dir.joinpath(f"{sample_id}{self.file_ex}")
        # self.bam_dir = Path(f"/home/ubuntu/NGS_data/current_data/{self.file_name + self.file_ex}")

    @staticmethod
    def _make_logger(log_dir: Path, name=None, consoleset=True, streamset=True) -> logging:
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
            loggerfile = log_dir.joinpath(name)
            file_handler = logging.FileHandler(filename=f"{loggerfile}.log")
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        return logger

    def load_file(self) -> pd.DataFrame:
        # sample_id = self.bam_dir.name.split(".")[0]
        with open(self.recal_metric_fn) as f:
            lines = f.readlines()
        header_line = lines[6].replace("\n", '')
        metric_line = lines[7].replace("\n", '')
        file_series = pd.Series(data=metric_line.split("\t"), index=header_line.split("\t"))
        file_series['Sample_ID'] = self.sample_id
        file_dataframe = pd.DataFrame(file_series).T
        return file_dataframe
        
    def write_to_sql(self, data_df: pd.DataFrame):
        db_type='mariadb'
        id='root'
        pw='gw12341234'
        host='db-qc.cg10utv5zpvf.ap-northeast-2.rds.amazonaws.com'
        schema_name='bam_practice'
        url =f"{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}"
        engine = db.create_engine(url)
        conn = engine.connect()
        data_df.to_sql(name='BAM_PRACTICE', con=conn, if_exists='append',index=False)
    
    def __call__(self):
        self.logger.info("Start BAM QC upload pipeline")
        self.logger.info("Load bam QC file...")
        data_df = self.load_file()
        self.logger.info("Sliming data frame...")
        data_df = data_df[
            [
                'Sample_ID',
                'MEAN_TARGET_COVERAGE',
                'PCT_USABLE_BASES_ON_TARGET',
                'PCT_EXC_DUPE',
                'PCT_PF_READS',
                'PCT_EXC_OVERLAP'
            ]
        ]
        self.logger.info("QC data database upload start!")
        try:
            self.write_to_sql(data_df)
        except Exception as e:
            self.logger.info("Error!")
            self.logger.info(e)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", '-d')
    parser.add_argument("--sample_id", '-sid')
    args = parser.parse_args()
    runner = BamUpload(args.data_dir, args.sample_id)
    runner()
    
file_list = Path("/home/ubuntu/NGS_data/bam")
file=list(file_list.glob("*/*.recal.bam.metrics.txt"))

