"""
    date: 2022.10.13
    object: pickle module 사용, db 확인
    author: 
"""

import pickle 
import logging
import sqlalchemy
import numpy as np
import pandas as pd
from typing import Dict
from pathlib import Path
from pytz import timezone
from datetime import datetime
from sqlalchemy import create_engine

##db connect_ gh2 
# class DBConfig:
#     def __init__(self):
#         self.db_type = 'mariadb'
#         self.id = 'yckim'
#         self.pw =  'rladmsco12!'
#         self.host = 'geninus-maria-211117.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com'
#         self.schema_name = 'gh2'

#db connect용
class DBConfig:
    def __init__(self):
        self.db_type='mariadb'
        self.id='eckim'
        self.pw='rladmsco12!'
        self.host='geninus-maria-211117-dev.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com'
        self.schema_name='gh2'
        
class PickleUpload(DBConfig):
    def __init__(self, data_dir: str, file_name: str, info_name: str) -> None:
        '''
            Class initialize
        '''
        DBConfig.__init__(self)
        self.data_dir = Path(data_dir)
        self.file_name = file_name
        self.info_name = info_name
        self.log_dir = self.data_dir.joinpath("log")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self._make_logger(self.log_dir, 'PickleUpload')
        
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
    
    def load_pickle(self) -> dict:
        with open(f'{self.data_dir}/{self.file_name}','rb') as file:
            data = pickle.load(file)
        return data
    
    def connect_db(self) -> sqlalchemy.engine:
        url = f'{self.db_type}+pymysql://{self.id}:{self.pw}@{self.host}/{self.schema_name}'
        engine = create_engine(url)
        conn = engine.connect()
        return conn
    
    def pickle_table(self, info_name: str, data: dict) -> dict:
        for sample_id in data.keys():
            for read_id in data[sample_id].keys():
                update_table = data[sample_id][read_id][info_name]
                update_table.insert(1, 'SEQ_DATE', 'date')
                update_table.insert(2, 'FC_ID', 'fc_id')
                update_table.insert(3, 'SAMPLE_ID', sample_id)
                update_table.insert(5, 'IDX', 0)
                update_table.insert(4, 'FASTQ_TYPE', read_id)
                update_table['CREATE_DATE'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                update_table = update_table.rename(columns={'Base':"BASE"})
        return update_table
    
    def _replace_space(self, update_table: pd.DataFrame) -> pd.DataFrame :
        col_list = []
        for col in update_table.columns.tolist():
            new_col = col.replace(" ","_")
            col_list.append(new_col)
        update_table.columns = col_list
        return update_table
    
    def upload_to_sql(self, info_name: str, update_table: pd.DataFrame, conn: sqlalchemy.engine):
        update_table.to_sql(f'gc_rsc_p{info_name}', con=conn, if_exists='append', index = False)
               
    def __call__(self):
        self.logger.info("Start pickle file upload to gc_qc_pbsc")
        self.logger.info("Load pickle data and parsing the data...")
        
        data = self.load_pickle()
        info_name = self.info_name
        conn= self.connect_db()
        
        self.logger.info("Revise the database...")
        update_table = self.pickle_table(info_name, data)
        
        self.logger.info("Replace space to under bar...")
        update_table = self._replace_space(update_table)
        
        self.logger.info(f"Upload to database_gc_p{info_name}...")
        self.upload_to_sql(info_name, update_table, conn)
            
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", '-d')
    parser.add_argument("--file_name", '-fd')
    parser.add_argument("--info_name","-in")
    args = parser.parse_args()
    runner = PickleUpload(args.data_dir, args.file_name, args.info_name)
    runner()
