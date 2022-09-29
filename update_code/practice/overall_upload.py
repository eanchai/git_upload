"""
    date: 2022.09.28
    object: compare expr_line output & update gc_qc_bi
    author: 
"""

import logging
import numpy as np
import pandas as pd
import sqlalchemy as db
import matplotlib as mat
import matplotlib.pyplot as plt
from time import time
from pathlib import Path
from pytz import timezone
from datetime import datetime
import sqlalchemy
from sqlalchemy import create_engine, MetaData, Table, Column, Numeric, Integer, VARCHAR, update

#db connect용
class DBConfig:
    def __init__(self):
        self.db_type = 'mariadb'
        self.id = 'yckim'
        self.pw =  'rladmsco12!'
        self.host = 'geninus-maria-211117.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com'
        self.schema_name = 'gh2'
        
class OverallUpdate:
    def __init__(self, data_dir: str, sample_id: str):
        '''
        class initaialize
        '''
        DBConfig.__init__(self)
        self.file_ex = ".csv"
        self.sample_id = sample_id
        self.data_dir = Path(data_dir)
        common = "fastqc_product_results_"
        self.log_dir = self.data_dir.joinpath("log")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self._make_logger(self.log_dir, 'OverallUploadTest')
        self.qc_df = self.data_dir.joinpath(f"{common}{sample_id}{self.file_ex}")
        
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
    
    def connect_db(self) -> sqlalchemy.engine:
        url = f'{self.db_type}+pymysql://{self.id}:{self.pw}@{self.host}/{self.schema_name}'
        engine = create_engine(url)
        conn = engine.connect()
        return conn
    
    #threshold: tb_expr_se_line_df connect
    def parse_db_to_df(self, conn, query: str) -> pd.DataFrame:
        tb_expr_seq_line_df = pd.read_sql(query, conn)
        tb_expr_seq_line_df = tb_expr_seq_line_df.rename(columns={'sample_id':'SampleID'})
        return tb_expr_seq_line_df
    
    #call에서 query 날려서 common부분 filter 추가해서 merge생성
    def filter_merged_qc_df(self, qc_df: pd.DataFrame, tb_expr_seq_line_df: pd.DataFrame) -> pd.DataFrame:
        qc_df['SampleID'] = [f"'{x}'" for x in qc_df['SampleID']] 
        qc_df['SampleID'] = [x.replace("'","") for x in qc_df['SampleID']]
        qc_merged = pd.merge(qc_df, tb_expr_seq_line_df, how='left', on='SampleID')
        qc_merged = qc_merged.rename(columns={'SampleID':'SAMPLE_ID'})
        qc_merged['FILTER'] = ['PASS' if row['FASTQ_TOTAL_BASES(Gb)'] > row['data_output'] else 'FAIL' for _, row in qc_merged.iterrows()]
        return qc_merged
    
    #merge에서 db table처럼 table 수정, filter_merged_qc_df값을 여기 아래 넣으면,, 될걸
    def revise_df(self, qc_merged: pd.DataFrame) -> pd.DataFrame:
        #원래 connect_db랑 연결할 수 있는데 test에 올리니까..
        #revise_df(self, conn, qc_merged)해서 쓰면 될듯..?
        db_type = 'mariadb'
        id = 'root'
        pw =  'gw12341234'
        host = 'db-qc.cg10utv5zpvf.ap-northeast-2.rds.amazonaws.com'
        schema_name = 'geninus'
        url = f'{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}'
        engine = create_engine(url)
        meta = MetaData(bind=engine)
        MetaData.reflect(meta)
        con = engine.connect()
        
        mytable = Table('gc_qc_bi', meta)
        meta.create_all(engine)
        table_cols = [x.name for x in mytable.columns]
        not_in_table_cols = [x for x in table_cols if x not in qc_merged.columns]
        # not_in_table_cols = [x for x in not_in_table_cols if x not in ['IDX', 'CREATE_DATE', 'LAST_UPDATE_DATE', 'CREATE_USER', 'LAST_UPDATE_USER']]
        for col in not_in_table_cols:
            qc_merged[col] = None
        qc_merged['LAST_UPDATE_DATE'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        qc_merged = qc_merged[table_cols]
        
        def update_sql(self, conn, qc_merged: pd.DataFrame):
            qc_merged.to_sql(name='gc_qc_bi', con = conn, if_exists='replace', index=False)

def __call__(self, qc_path: Path):