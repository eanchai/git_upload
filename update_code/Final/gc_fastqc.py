"""
    date: 2022.10.11
    object: compare expr_line output & update new Table(테이블 명)
    author: 
"""

import logging
import sqlalchemy
import numpy as np
import pandas as pd
import sqlalchemy as db
import matplotlib as mat
import matplotlib.pyplot as plt
from time import time
from pathlib import Path
from pytz import timezone
from datetime import datetime
from sqlalchemy import create_engine, MetaData, Table, Column, Numeric, Integer, VARCHAR, update

##db connect_ gh2 
# class DBConfig:
#     def __init__(self):
#         self.db_type = 'mariadb'
#         self.id = 'yckim'
#         self.pw =  'rladmsco12!'
#         self.host = 'geninus-maria-211117.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com'
#         self.schema_name = 'gh2'

#db connect_geninus dump _ 연습용 서버
class DBConfig:
    def __init__(self):
        self.db_type = 'mariadb'
        self.id = 'root'
        self.pw =  'gw12341234'
        self.host = 'db-qc.cg10utv5zpvf.ap-northeast-2.rds.amazonaws.com'
        self.schema_name = 'geninus'
        
class FastqcUpdate(DBConfig):
    def __init__(self, data_dir: str, sample_id: str):
        '''
        Class initialize
        '''
        DBConfig.__init__(self)
        self.file_ex = ".csv"
        self.sample_id = sample_id
        self.data_dir = Path(data_dir)
        common = "fastqc_product_results_"
        #logger 정의
        self.log_dir = self.data_dir.joinpath("log")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self._make_logger(self.log_dir, 'OverallUploadTest')
        #qc_df 미리 여기서 불러오기
        self.qc_df = self.data_dir.joinpath(f"{common}{self.sample_id}{self.file_ex}")
        
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
    
    #db connect 정의
    def connect_db(self) -> sqlalchemy.engine:
        url = f'{self.db_type}+pymysql://{self.id}:{self.pw}@{self.host}/{self.schema_name}'
        engine = create_engine(url)
        conn = engine.connect()
        return conn
    
    #thresholds는 tb_expr_seq_line_df에 존재, output이랑 비교하기 위해서 connect
    def parse_db_to_df(self, conn, query: str) -> pd.DataFrame:
        tb_expr_seq_line_df = pd.read_sql(query, conn)
        tb_expr_seq_line_df = tb_expr_seq_line_df.rename(columns={'sample_id':'SampleID'})
        return tb_expr_seq_line_df
    
    #call에서 query 날려서 common merge & filter
    def filter_merged_qc_df(self, qc_df: pd.DataFrame, tb_expr_seq_line_df: pd.DataFrame) -> pd.DataFrame:
        qc_df['SampleID'] = [f"'{x}'" for x in qc_df['SampleID']] 
        qc_df['SampleID'] = [x.replace("'","") for x in qc_df['SampleID']]
        qc_merged = pd.merge(qc_df, tb_expr_seq_line_df, how='left', on='SampleID')
        
        qc_merged_df = qc_merged.rename(columns={'SampleID':'SAMPLE_ID'})
        qc_merged_df['FASTQ_OVERALL'] = ['PASS' if row['FASTQ_TOTAL_BASES(Gb)'] > row['data_output'] else 'FAIL' for _, row in qc_merged.iterrows()]
        qc_merged_df['FASTQ_TOTAL_READ'] = qc_merged_df['FASTQ_TOTAL_READ_R1'] + qc_merged_df['FASTQ_TOTAL_READ_R2']
        return qc_merged_df
    
    #db table처럼 table 수정, call에서 qc_merged return 되는 filter_merged_qc_df 값 넣기
    def revise_df(self, qc_merged_df: pd.DataFrame) -> pd.DataFrame:
        db_type = 'mariadb'
        id = 'root'
        pw =  'gw12341234'
        host = 'db-qc.cg10utv5zpvf.ap-northeast-2.rds.amazonaws.com'
        schema_name = 'geninus'
        url = f'{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}'
        engine = create_engine(url)
        meta = MetaData(bind=engine)
        MetaData.reflect(meta)
        mytable = Table('gc_qc_bi', meta)
        meta.create_all(engine)
        
        table_cols = [x.name for x in mytable.columns]
        not_in_table_cols = [x for x in table_cols if x not in qc_merged_df.columns]
        # not_in_table_cols = [x for x in not_in_table_cols if x not in ['IDX', 'CREATE_DATE', 'LAST_UPDATE_DATE', 'CREATE_USER', 'LAST_UPDATE_USER']]
        for col in not_in_table_cols:
            qc_merged_df[col] = None
        qc_merged_df['LAST_UPDATE_DATE'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        qc_upload_df = qc_merged_df[table_cols]
        return qc_upload_df
    
    #gh2 업로드 시 이렇게 쓰면 될듯
    #def update_sql(self, conn, qc_merged: pd.DataFrame):
    #    qc_merged.to_sql(name='gc_qc_bi', con = conn, if_exists='replace', index=False)
    def update_sql(self, qc_upload_df: pd.DataFrame):
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
        
        qc_upload_df.to_sql(name='gc_qc_bi', con = con, if_exists='append', index=False)
        
    def __call__(self):
        self.logger.info("Start fastqc_overall update pipeline")
        self.logger.info("Connect db and parsing the data...")
        
        #DB connection initialize
        conn = self.connect_db()

        # load QC df
        qc_df = pd.read_csv(self.qc_df)
        qc_df['SampleID'] = [f"'{x}'" for x in qc_df['SampleID']] 
        
        # load tb_expr_seq_line from database
        query = f'''
            select sample_id, data_output from tb_expr_seq_line
            where sample_id in ({",".join(qc_df['SampleID'].tolist())})
            '''
        tb_expr_seq_line_df = self.parse_db_to_df(conn, query)
    
        self.logger.info("Filtering Data with threshold...")
        qc_merged_df = self.filter_merged_qc_df(qc_df, tb_expr_seq_line_df)
    
        self.logger.info("Revise the table form...")
        qc_upload_df = self.revise_df(qc_merged_df)
    
        self.logger.info("Update test database...")
        try:
            self.update_sql(qc_upload_df)
        except Exception as e:
            self.logger.info("Error in update_sql")
            self.logger.info(e)
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", '-d')
    parser.add_argument("--sample_id", '-sid')
    args = parser.parse_args()
    runner = FastqcUpdate(args.data_dir, args.sample_id)
    runner()