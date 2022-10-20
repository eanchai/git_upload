"""
    date: 2022.10.12 ~
    object: BAM metric upload
        check BAM file & gc_bamqc upload
    author: 
"""

import logging
import pandas as pd
import sqlalchemy
from time import time
from pathlib import Path
from pytz import timezone
from datetime import datetime
import sqlalchemy
from sqlalchemy import create_engine, MetaData, Table, Column, Numeric, Integer, VARCHAR, update

class DBConfig:
    def __init__(self):
        self.db_type='mariadb'
        self.id='eckim'
        self.pw='rladmsco12!'
        self.host='geninus-maria-211117-dev.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com'
        self.schema_name='gh2'

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
        self.logger = self._make_logger(self.log_dir, 'BamUpload')
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
    
    #db에 연결
    def connect_db(self) -> sqlalchemy.engine:
        url = f'{self.db_type}+pymysql://{self.id}:{self.pw}@{self.host}/{self.schema_name}'
        engine = create_engine(url)
        conn = engine.connect()
        return conn
    
    # split sequence date and fc_id
    #이거 해야하는데 현재 테스트 파일 올릴때는 이걸 지금 올릴 수가 없음
    # def _split_fc_dir_to_date_and_id(self) -> str | str:
    #     tokens = self.fc_dir.split("_")
    #     seq_date = tokens[0]
    #     fc_id = tokens[-1]
    #     fc_id = fc_id[1:]
    #     return seq_date, fc_id

    #bam파일에서 원하는 부분 추출해서 bam_df로 저장
    def load_file(self) -> pd.DataFrame:
        with open(self.recal_metric_fn) as f:
            lines = f.readlines()
        header_line = lines[6].replace("\n", '')
        metric_line = lines[7].replace("\n", '')
        file_series = pd.Series(data=metric_line.split("\t"), index=header_line.split("\t"))
        bam_df = pd.DataFrame(file_series).T
        return bam_df
    
    def parse_bam_df(self, bam_df: pd.DataFrame) -> pd.DataFrame:
        # insert sample id
        # seq_date, fc_id = self._split_fc_dir_to_date_and_id()
        seq_date = "date"
        fc_id = 'flow_cell_id'
        bam_df.insert(0, 'SEQ_DATE', seq_date)
        bam_df.insert(1, 'FC_ID', fc_id)
        bam_df.insert(2,'SAMPLE_ID', self.sample_id)
        bam_df['CREATE_USER'] = self.id
        bam_df['CREATE_DATE'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        return bam_df
    
    def write_to_sql(self, conn: sqlalchemy.engine, bam_df: pd.DataFrame):
        bam_df.to_sql(name='gc_rsc_bamqc', con=conn, if_exists='append', index=False)
    
    def __call__(self):
        self.logger.info("Start BAM QC upload pipeline")
        self.logger.info("Start fastqc_overall update pipeline")
        self.logger.info("Connect db and parsing the data...")
        
        #DB connection
        conn = self.connect_db()
        
        self.logger.info("Load bam QC file and Revise...")
        bam_df = self.load_file()
        parsed_bam_df = self.parse_bam_df(bam_df)
        
        self.logger.info("Upload to database_gc_bamqc")
        try:
            self.write_to_sql(conn, parsed_bam_df)
        except Exception as e:
            self.logger.info("Error in upload sql")
            self.logger.info(e)
            
            
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", '-d')
    parser.add_argument("--sample_id", '-sid')
    args = parser.parse_args()
    runner = BamUpload(args.data_dir, args.sample_id)
    runner()      