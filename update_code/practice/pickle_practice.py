"""
    date: 2022.09.14 ~
    object: pickle module 사용, db 확인
    author: 
"""

import logging
import numpy as np
import pickle 
import pandas as pd
import sqlalchemy as db
import matplotlib as mat
import matplotlib.pyplot as plt
from pathlib import Path
from sqlalchemy import create_engine, Table, Column

#data 불러오기
db_type='mariadb'
id='yckim'
pw="rladmsco12!"
host='geninus-maria-211117.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com'
schema_name='gh2'
url =f"{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}"

engine = create_engine(url)
conn = engine.connect()
query = 'select SAMPLE_ID, BASE, FASTQ_TYPE, G, A, T, C from gc_chart_pbsc order by SAMPLE_ID;'
data = pd.read_sql(query, conn)

#db update
with open("data.pickle", "wb") as f:
    pickle.dump(data, f)
    
with open("data.pickle","rb") as fi:
    test = pickle.load(fi)
 
[t['SAMPLE_ID'] for _, t in test.key()]