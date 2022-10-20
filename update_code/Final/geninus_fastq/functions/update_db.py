#!env /usr/bin/python3

'''

Date: 2021.11.11
Authors: duaghk

FASTQ information update.

'''

# import library

import logging
import sqlalchemy
import pandas as pd
from pathlib import Path
from datetime import datetime
from sqlalchemy import create_engine
from .configs import DBConfig, ToolConfig, SequencingConfig


class UpdateFastQC(DBConfig, ToolConfig, SequencingConfig):
    def __init__(self, date):
        DBConfig.__init__(self)
        ToolConfig.__init__(self)
        SequencingConfig.__init__(self)
        self.log_dir = self.log_dir.joinpath("UpdateQC")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self.make_logger(f"UpdateQC_{date}")
        
    def make_logger(self, name=None, consoleset=True, streamset=True) -> logging:
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
            loggerfile = self.log_dir.joinpath(f"{name}.log")
            file_handler = logging.FileHandler(filename=loggerfile)
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        return logger
    
    #dict가 pickle 파일 -> 이게 fastqc_dict를 받아야함
    def _split_fc_dir_to_date_and_id(self, fastqc_result_path: str) -> str :
        tokens = fastqc_result_path.split("_")
        seq_date = tokens[0]
        fc_id = tokens[-1]
        fc_id = fc_id[1:]
        return seq_date, fc_id

    def parse_df(self, sample_id: str, read_type: str, target_df: pd.DataFrame, fastqc_result_path: str) -> pd.DataFrame:
        seq_date, fc_id = self._split_fc_dir_to_date_and_id(fastqc_result_path)        
        target_df.insert(0, 'SEQ_DATE', seq_date)
        target_df.insert(1, 'FC_ID', fc_id)
        target_df.insert(2, 'SAMPLE_ID', sample_id)
        target_df.insert(3, 'TYPE', read_type)
        return target_df

    def update_fastqc_dict(self, qc_dict: dict, fastqc_result_path: str) -> None:
        target_key_dict = {'bsq':"gc_rsc_qbsq", 'bsc':"gc_rsc_qbsc", 'sge':"gc_rsc_qsge"}
        for sample_id, read_type in qc_dict.items():
            for key, table_name in target_key_dict.items():
                target_df = qc_dict[sample_id][read_type][key]
                target_df = self.parse_df(sample_id, read_type, target_df, fastqc_result_path)
                target_df.to_sql(name=table_name, con=conn, if_exists='append', index=False)

    def connect_db(self) -> sqlalchemy.engine:
        url = f'{self.db_type}+pymysql://{self.id}:{self.pw}@{self.db_address}/{self.db_name}'
        engine = create_engine(url)
        conn = engine.connect()
        return conn
    
    #csv로 나오는게 fastqc -> 이게 product_table을 받아야함
    def parse_db_to_df(self, conn, query: str) -> pd.DataFrame:
        # db_type = 'mariadb'
        # id = 'yckim'
        # pw =  'rladmsco12!'
        # host = 'geninus-maria-211117.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com'
        # schema_name = 'gh2'
        # url = f'{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}'
        # engine = create_engine(url)
        # con = engine.connect()
        tb_expr_seq_line_df = pd.read_sql(query, conn)
        tb_expr_seq_line_df = tb_expr_seq_line_df.rename(columns={'sample_id':'SampleID'})
        return tb_expr_seq_line_df
    
    def merged_qc_expr(self, qc_df: pd.DataFrame, tb_expr_seq_line_df: pd.DataFrame) -> pd.DataFrame:
        qc_merged = pd.merge(qc_df, tb_expr_seq_line_df, how='left', on='SampleID')
        qc_merged_df = qc_merged.rename(columns={'SampleID':'SAMPLE_ID'})
        qc_merged_df = qc_merged_df.rename(columns={'data_output':'DATA_OUTPUT'})
        qc_merged_df['FASTQ_OVERALL'] = ['PASS' if row['FASTQ_TOTAL_BASES(Gb)'] > row['data_output'] else 'FAIL' for _, row in qc_merged.iterrows()]
        qc_merged_df['FASTQ_TOTAL_READ'] = qc_merged_df['FASTQ_TOTAL_READ_R1'] + qc_merged_df['FASTQ_TOTAL_READ_R2']
        qc_merged_df['CREATE_DATE'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        return qc_merged_df
    
    def add_flowcell(self, qc_merged_df: pd.DataFrame, fastqc_result_path: str) -> pd.DataFrame:
        seq_date, fc_id = self._split_fc_dir_to_date_and_id(fastqc_result_path)
        qc_merged_df.insert(0,'SEQ_DATE', seq_date)        
        qc_merged_df.insert(1,'FC_ID', fc_id)
        return qc_merged_df

    def update_sql(self, conn: create_engine, qc_fastqc_update_df: pd.DataFrame):
        qc_fastqc_update_df.to_sql(name="gc_rsc_fastqc", con=conn, if_exists='append', index = False)

    def update_fastqc_table(self, conn:create_engine, qc_fastqc_update_df: pd.DataFrame) -> None:
        try:
            self.update_sql(conn, qc_fastqc_update_df)
        except Exception as e:
            self.logger.info("Error in update sql")                       
            self.logger.info(e)
        return

    def __call__(self, fastqc_result_path, qc_dict, qc_table) -> None:
        self.update_fastqc_dict(qc_dict, fastqc_result_path)

        qc_df = qc_table.copy()
        qc_df['SampleID'] = [f"'{x}'" for x in qc_table['SampleID']] 

        query = f'''
            select sample_id, data_output from tb_expr_seq_line
            where sample_id in ({",".join(qc_df['SampleID'].tolist())})
            '''
            
        tb_expr_seq_line_df = self.parse_db_to_df(query)
        self.logger.info("Filtering Data with threshold...")
        qc_merged_df = self.merged_qc_expr(qc_table, tb_expr_seq_line_df)
        
        self.logger.info("Add flow cell ID on db...")
        qc_fastqc_update_df = self.add_flowcell(qc_merged_df, fastqc_result_path)
        self.update_fastqc_table(qc_fastqc_update_df)
        pass
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--date")
    parser.add_argument("--fastqc_dir")
    args = parser.parse_args()
    
    dbupdatater = UpdateFastQC(args.date)
    dbupdatater(Path(args.fastqc_dir))




