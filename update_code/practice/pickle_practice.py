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

#data 불러오기, 이건 sqlalchemy core를 사용하는 방식
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

#pickle 파일을 sqlalchemy를 사용해서 처리하는 방식, pickle파일 내부 확인
import pickle
import sqlalchemy as db

#engine연결
db_type='mariadb'
id='yckim'
pw="rladmsco12!"
host='geninus-maria-211117.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com'
schema_name='gh2'
url =f"{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}"
engine = create_engine(url)
conn = engine.connect()
metadata = db.MetaData()
table = db.Table('gc_chart_pbsc', metadata, autoload=True, autoload_with=engine)

values_list = [info for info in data[sample_id]['R1']['bsc'].items()]

query = db.insert('gc_chart_pbsc')
result_proxy = conn.execute(query, values_list)
result_proxy.close()
    

#pickle파일 내부 확인
 
with open('/home/ubuntu/NGS_data/pickle/fastqc_result.pickle', 'rb') as file:
    data = pickle.load(file)

sample_id = list(data.keys())[0]

data[sample_id].keys()

data[sample_id]['R1'].keys()

data[sample_id]['R1']['bsq']
data[sample_id]['R1']['bsc']
data[sample_id]['R1']['sge']
data[sample_id]['R1']['Total_bases']
data[sample_id]['R1']['Total_reads']
data[sample_id]['R1']['Q30']

sample_id = [sample_id for sample_id in list(data.keys())]

for sample in data.keys():
    for read_id in data[sample].keys():
        for info_name, info in data[sample][read_id].items():
            print(info_name)
            print(info)
            print()

import sqlalchemy
from sqlalchemy import insert, update
df = data[sample_id]['R1']['bsc']
df.update('gc_chart_pbsc', conn)



list(table.columns)


info_name = 'bsc'
for sample_id in data.keys():
    for read_id in data[sample].keys():
        update_table = data[sample_id][read_id]['bsc']
        update_table.insert(0, 'SAMPLE_ID', sample_id)
        update_table.insert(1, 'FASTQ_TYPE', read_id)
        update_table = update_table.rename(columns={'Base':"BASE"})
        update_list = [row.to_dict() for _, row in update_table.iterrows()]
        result_proxy = conn.execute(query, update_list)
        result_proxy.close()
        
# import library
import pickle
import sqlalchemy
from sqlalchemy import insert, update, create_engine, 
from datetime import datetime

# connect db
# engine연결
db_type='mariadb'
id='yckim'
pw="rladmsco12!"
host='geninus-maria-211117.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com'
schema_name='gh2'


db_type='mariadb'
id='root'
pw='gw12341234'
host='db-qc.cg10utv5zpvf.ap-northeast-2.rds.amazonaws.com'
schema_name='geninus'



url =f"{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}"
engine = create_engine(url)
conn = engine.connect()
metadata = sqlalchemy.MetaData()
table = sqlalchemy.Table('gc_chart_pbsc', metadata, autoload=True, autoload_with=engine)



# load data
 
with open('/home/ubuntu/NGS_data/pickle/fastqc_result.pickle', 'rb') as file:
    data = pickle.load(file)



# update data

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
        




