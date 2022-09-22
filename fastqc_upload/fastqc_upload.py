"""
    date: 2022.09.14 ~
    object: fastqc metric upload
        check fastqc file & sql database upload
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

class FastqcUpload:
    def __init__(self, data_dir: str, sample_id: str) -> None:
        '''
            Class initialize
        '''
        self.data_dir = Path(data_dir)
        self.sample_id = sample_id
        self.log_dir = self.data_dir.joinpath("log")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self._make_logger(self.log_dir, 'FastqcUploadTest')
        self.file_ex = ".csv" #suffix
        common = "fastqc_product_results_"
        self.fastqc_fn = self.data_dir.joinpath(f"{common}{sample_id}{self.file_ex}")
 
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
        file_dataframe = pd.DataFrame(pd.read_csv(self.fastqc_fn))
        return file_dataframe
        
    def write_to_sql(self, data_df: pd.DataFrame):
        db_type='mariadb'
        id='root'
        pw='gw12341234'
        host='db-qc.cg10utv5zpvf.ap-northeast-2.rds.amazonaws.com'
        schema_name='fastqc_practice'
        url =f"{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}"
        engine = db.create_engine(url)
        conn = engine.connect()
        data_df.to_sql(name='FASTQC_PRACTICE', con=conn, if_exists='append',index=False)
    
    def __call__(self):
        self.logger.info("Start fastQC upload pipeline")
        self.logger.info("Load fastQC file...")
        data_df = self.load_file()
        self.logger.info("Sliming data frame...")
        data_df = data_df[
            [
                'SampleID',
                'FASTQ_TOTAL_READ_R1',
                'FASTQ_TOTAL_BASE_R1',
                'FASTQ_Q30_R1',
                'FASTQ_TOTAL_READ_R2',
                'FASTQ_TOTAL_BASE_R2',
                'FASTQ_Q30_R2',
                'FASTQ_TOTAL_BASES(Gb)'
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
    runner = FastqcUpload(args.data_dir, args.sample_id)
    runner()