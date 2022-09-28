"""
    date: 2022.09.23
    object: check base
    author: 
"""

import logging
import numpy as np
import pandas as pd
import seaborn as sns
import subprocess as sp
import sqlalchemy as db
import matplotlib as mat
import matplotlib.pyplot as plt
from time import time
from pathlib import Path
from pytz import timezone
from datetime import datetime
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

#일단 하나 샘플만 빼서 진행, 샘플 넣는걸로 def 짜면 될듯
data_19 = data[data.SAMPLE_ID == 'CD_19_08402_DE_CS']
data_base = data_19[['G','A','T','C']]
fig, ax = plt.subplots(figsize=(12,6))
ax.set(ylim = (0, 100))
for _, base in enumerate(data_base):
    sns.lineplot(data = data_19,
                 x = "BASE",
                 y = base,
                 label = base)
    plt.xticks(rotation = 45)


