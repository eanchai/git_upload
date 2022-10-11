"""
    date: 2022.09.30
    object: BAM metric upload
        check BAM file & gc_qc_bi update
    author: 
"""

import logging
import numpy as np
import pandas as pd
import sqlalchemy
import matplotlib as mat
import matplotlib.pyplot as plt
from time import time
from pathlib import Path
from pytz import timezone
from datetime import datetime
import sqlalchemy
from sqlalchemy import create_engine, MetaData, Table, Column, Numeric, Integer, VARCHAR, update

class DBConfig:
    def __init__(self):
        self.db_type='mariadb'
        self.id='root'
        self.pw='gw12341234'
        self.host='db-qc.cg10utv5zpvf.ap-northeast-2.rds.amazonaws.com'
        self.schema_name='geninus'

class BamUpload(DBConfig):
    def __init__(self, data_dir: str, sample_id: str) -> None:
        '''
            Class initialize
        '''
        DBConfig.__init__(self)
        self.sample_id = sample_id
        self.data_dir = Path(data_dir)
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

    #bam 파일에서 필요한 부분 추출 -> data_df로 저장
    def load_file(self) -> pd.DataFrame:
        # sample_id = self.bam_dir.name.split(".")[0]
        with open(self.recal_metric_fn) as f:
            lines = f.readlines()
        header_line = lines[6].replace("\n", '')
        metric_line = lines[7].replace("\n", '')
        file_series = pd.Series(data=metric_line.split("\t"), index=header_line.split("\t"))
        file_series['Sample_ID'] = self.sample_id
        data_df = pd.DataFrame(file_series).T
        
        data_df  = data_df[
            [
                'Sample_ID',
                'MEAN_TARGET_COVERAGE',
                'PCT_USABLE_BASES_ON_TARGET',
                'PCT_EXC_DUPE',
                'PCT_PF_READS',
                'PCT_EXC_OVERLAP'
            ]
        ]
        return data_df
    
    #db에 연결
    def connect_db(self) -> sqlalchemy.engine:
        url = f'{self.db_type}+pymysql://{self.id}:{self.pw}@{self.host}/{self.schema_name}'
        engine = create_engine(url)
        conn = engine.connect()
        
        meta = MetaData(bind=engine)
        MetaData.reflect(meta)
        metadata = sqlalchemy.MetaData()
        table = sqlalchemy.Table('gc_bamqc', metadata, autoload=True, autoload_with=engine)
        return conn, table
    
    ##data_df에서 bam_df로 변환
    def revise_df(self, data_df: pd.DataFrame) -> pd.DataFrame:
        bam_df = data_df.rename(columns={
                'SampleID':'SAMPLE_ID',
                'MEAN_TARGET_COVERAGE':'BAM_MEAN_DEPTH',
                'PCT_USABLE_BASES_ON_TARGET':'BAM_ON_TARGET_RATE'
                })
        
        return bam_df

    #sql에 업로드 
    def write_to_sql(self, table: sqlalchemy.Table, conn: sqlalchemy.engine, bam_df: pd.DataFrame):
        update_table = bam_df
        # update_table.insert(0, 'SAMPLE_ID', sample_id)
        # update_table.insert(1, 'FASTQ_TYPE', read_id)
        # update_table = update_table.rename(columns={'Base':"BASE"})
        update_list = [row.to_dict() for _, row in update_table.iterrows()]
        query = sqlalchemy.update(table).values(update_list)
        result_proxy = conn.execute(query, update_list)
        result_proxy.close()
        
        #bam_df.to_sql(name='qc_bamqc', con=conn, if_exists='replace', index=False)
             
    def __call__(self):
        self.logger.info("Start BAM QC upload pipeline")
        self.logger.info("Start fastqc_overall update pipeline")
        self.logger.info("Connect db and parsing the data...")
        
        #DB connection initialize
        table, conn = self.connect_db()
        
        self.logger.info("Load bam QC file...")
        data_df = self.load_file()
        
        self.logger.info("Revise the Bam dataframe...")
        bam_df = self.revise_df(data_df)
        
        self.logger.info("QC data database upload start!")
        try:
            # query = self.write_to_sql(table, conn, bam_df)
            # bam_df.to_sql(name='gc_qc_bi', con=conn, if_exists='replace', index=False, query = query)
            self.write_to_sql(table, conn, bam_df)
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

