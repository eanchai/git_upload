"""
    date: 2022.10.05
    object: pickle module 사용, db 확인
    author: 
"""

import pickle 
import logging
from turtle import st
from typing import Dict
import numpy as np
import pandas as pd
import sqlalchemy as db
import matplotlib as mat
import matplotlib.pyplot as plt
from pathlib import Path
from pytz import timezone
from datetime import datetime
from sqlalchemy import create_engine, Table, Column

class PickleUpload():
    def __init__(self, data_dir: str, file_name: str) -> None:
        '''
            Class initialize
        '''
        self.data_dir = data_dir
        self.file_name = file_name
        
        
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
    
    def load_pickle(self, data_dir: str, file_name: str) -> dict:
        with open(f'{data_dir}+{file_name}','rb') as file:
            data = pickle.load(file)
        return data
    
    def append_id(self, data: dict) -> dict:
        sample_id = [sample_id for sample_id in list(data.keys())]
        for sample in data.keys():
        #    print(sample)
            for read_id in data[sample].keys():
                for info_name, info in data[sample][read_id].items():
                    print(info_name)
                    print(info)
                    
