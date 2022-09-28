"""
    date: 2022.09.14 ~
    object: seq_line_output 비교 및 업데이트 
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

#tb에서 fastqc 파일과 동일한걸 뽑아서 긁기
#dir, file name 설정해서 변수에 넣기 
class OverallUpdate:
    def __init__(self, data_dir: str, sample_id: str):
        '''
        class initaialize
        '''
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
        
    def compare_output(self) -> pd.DataFrame:
        #db connect
        db_type = 'mariadb'
        id = 'yckim'
        pw =  'rladmsco12!'
        host = 'geninus-maria-211117.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com'
        schema_name = 'gh2'
        url = f'{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}'
        engine = create_engine(url)
        conn = engine.connect()
        #data merge
        qc_df = pd.read_csv(self.qc_df)
        qc_df['SampleID'] = [f"'{x}'" for x in qc_df['SampleID']] 
        qc_df
        query = f'''
        select sample_id, data_output from tb_expr_seq_line
        where sample_id in ({",".join(qc_df['SampleID'].tolist())})
        '''
        tb_expr_seq_line = pd.read_sql(query, conn)
        qc_df['SampleID'] = [x.replace("'","") for x in qc_df['SampleID']]
        tb_expr_seq_line = tb_expr_seq_line.rename(columns={'sample_id':'SampleID'})
        qc_merged = pd.merge(qc_df, tb_expr_seq_line, how='left', on='SampleID')
        qc_merged = qc_merged.rename(columns={'SampleID':'SAMPLE_ID'})
        qc_merged['FILTER'] = ['PASS' if row['FASTQ_TOTAL_BASES(Gb)'] > row['data_output'] else 'FAIL' for _, row in qc_merged.iterrows()]
        return qc_merged
    
    def revise_df(self):
        ## connect server
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
        qc_merged = self.compare_output()
        table_cols = [x.name for x in mytable.columns]
        not_in_table_cols = [x for x in table_cols if x not in qc_merged.columns]
        # not_in_table_cols = [x for x in not_in_table_cols if x not in ['IDX', 'CREATE_DATE', 'LAST_UPDATE_DATE', 'CREATE_USER', 'LAST_UPDATE_USER']]
        for col in not_in_table_cols:
            qc_merged[col] = None
        qc_merged['LAST_UPDATE_DATE'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        
        qc_merged = qc_merged[table_cols]
        print(qc_merged)
    
    def update_sql(self):
        '''
        update는 원래 gh2, gc_qc_bi에 해야함
        '''
        ## connect server
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
        qc_merged = self.revise_df()
        
        qc_merged.to_sql('gc_qc_bi', con = con, if_exists='replace', index=False)
        ## query
        # for _, row in self.compare_output().iterrows():
        #     query = mytable.insert().values(
        #         SAMPLE_ID = {row['SampleID']},
        #         FASTQ_TOTAL_READ = {row['FASTQ_TOTAL_READ_R1'] + row['FASTQ_TOTAL_READ_R2']},
        #         FASTQ_Q30 = None,
        #         FASTQ_OVERALL = {row['FILTER']},
        #         FASTQ_TOTAL_READ_R1 = {row['FASTQ_TOTAL_READ_R1']},
        #         FASTQ_GC_CONTENTS_R1 = None ,
        #         FASTQ_OVERALL_R1 = None,
        #         FASTQ_Q30_R1 = {row['FASTQ_Q30_R1']},
        #         FASTQ_TOTAL_READ_R2 = {row['FASTQ_TOTAL_READ_R2']},
        #         FASTQ_GC_CONTENTS_R2 = None,
        #         FASTQ_OVERALL_R2 = None,
        #         FASTQ_Q30_R2 = {row['FASTQ_Q30_R2']},
        #         BAM_MEAN_DEPTH = None,
        #         BAM_CAPTURE_EFFICIENCY = None,
        #         BAM_ON_TARGET_RATE = None,
        #         BAM_DUP_RATE = None,
        #         BAM_PR_SCORE = None,
        #         BAM_UNIFORM = None,
        #         BAM_TUMOR_VOL = None,
        #         BAM_OVERALL_QUALITY = None,
        #         BAM_GENDER_ESTIMATED = None,
        #         ANAL_PF_QC = None,
        #         ANAL_PF_COMM = None,
        #         PCT_TARGET_04X_MEAN = None,
        #         ALIGN_RATE = None,
        #         UNIQUE_READ = None,
        #         UNIQUE_MEAN_DEPTH = None,
        #         BAM_MEAN_INSERT_SIZE = None,
        #         COEF = None,
        #         MEDIAN_EXON_COVERAGE = None
        #         # LAST_UPDATE_DATE = None,
        #         # CREATE_USER = '',
        #         # LAST_UPDATE_USER = 'root'
        # )

                     
    def __call__(self):
        self.logger.info("Start fastqc_overall update pipeline")
        self.logger.info("load fastqc file, Compare output..")
        self.compare_output()
        self.logger.info("Revise database...")
        self.revise_df()
        self.logger.info("Update database...")
        try:
            self.update_sql()
        except Exception as e:
            self.logger.info("Error!")
            self.logger.info(e)
            
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", '-d')
    parser.add_argument("--sample_id", '-sid')
    args = parser.parse_args()
    runner = OverallUpdate(args.data_dir, args.sample_id)
    runner()
    

