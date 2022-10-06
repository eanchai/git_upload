"""
    date: 2022.10.05
    object: pickle module 사용, db 확인
    author: 
"""

import pickle 
import logging
import sqlalchemy
import numpy as np
import pandas as pd
import matplotlib as mat
import matplotlib.pyplot as plt
from typing import Dict
from pathlib import Path
from pytz import timezone
from datetime import datetime
from sqlalchemy import create_engine, Table, Column

#db connect용
class DBConfig:
    def __init__(self):
        self.db_type='mariadb'
        self.id='root'
        self.pw='gw12341234'
        self.host='db-qc.cg10utv5zpvf.ap-northeast-2.rds.amazonaws.com'
        self.schema_name='geninus'

class PickleUpload(DBConfig):
    def __init__(self, data_dir: str, file_name: str) -> None:
        '''
            Class initialize
        '''
        DBConfig.__init__(self)
        self.data_dir = Path(data_dir)
        self.file_name = file_name
        self.log_dir = self.data_dir.joinpath("log")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self._make_logger(self.log_dir, 'PickleUploadTest')
           
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
        metadata = sqlalchemy.MetaData()
        table = sqlalchemy.Table('gc_chart_pbsc', metadata, autoload=True, autoload_with=engine)
        return conn, table
    
    def revise_bsc_to_sql(self, data: dict, table: sqlalchemy.Table, conn: sqlalchemy.engine) -> dict:
        info_name = 'bsc'
        for sample_id in data.keys():
            for read_id in data[sample_id].keys():
                update_table = data[sample_id][read_id]['bsc']
                update_table.insert(0, 'SAMPLE_ID', sample_id)
                update_table.insert(1, 'IDX', 0)
                update_table.insert(2, 'FASTQ_TYPE', read_id)
                update_table.insert(8, 'CREATE_DATE', datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                
                update_table = update_table.rename(columns={'Base':"BASE"})
                update_list = [row.to_dict() for _, row in update_table.iterrows()]
                query = sqlalchemy.update(table).values(update_list)
                result_proxy = conn.execute(query)
                # try:
                #     query = sqlalchemy.insert(table).values(update_list)
                #     result_proxy = conn.execute(query)
                # except ValueError:
                #     query = sqlalchemy.update(table).values(update_list)
                #     result_proxy = conn.execute(query)
                result_proxy.close()
                    
    def __call__(self):
        self.logger.info("Start pickle file upload to gc_qc_pbsc")
        self.logger.info("Load pickle data and parsing the data...")
        
        data = self.load_pickle()
        conn, table = self.connect_db()
        
        self.logger.info("Revise the database and update...")
        self.revise_bsc_to_sql(data, table, conn)
            
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", '-d')
    parser.add_argument("--file_name", '-fd')
    args = parser.parse_args()
    runner = PickleUpload(args.data_dir, args.file_name)
    runner()
