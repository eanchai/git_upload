"""
    date: 2022.09.14 ~
    object: AWS database와 연결
        NGS data/fastq 파일의 fastgc 파일에서 각각 종류 분리 및 quality 확인
    author: 
"""

from pkgutil import ImpImporter
from re import T


import sqlalchemy as db
from pathlib import Path
import pandas as pd

#database에 연결해서 가져오기 - data에 저장
db_type='mariadb'
id='admin'
pw="Geninus1!"
host='db-practice.cg10utv5zpvf.ap-northeast-2.rds.amazonaws.com'
schema_name='geninus'
url =f"{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}"

engine = db.create_engine(url)
conn = engine.connect()             #database 연결

query = 'select * from gc_qc_bi order by create_date desc;'
data = pd.read_sql(query, conn)
data.sql.name.contains("FASTQ")

data = pd.read_csv(Path())

data = data.rename(columns={name_a:name:b, ~~})

data.to_sql(conn, 'gh2')

#database 가져오기
engine = db.create_engine('mysql+pymysql://admin:Geninus1!@db-practice.cg10utv5zpvf.ap-northeast-2.rds.amazonaws.com/geninus')
connection = engine.connect()
metadata = db.MetaData()
table = db.Table('gc_qc_bi', metadata, autoload = True, autoload_with= engine)


conn.close()