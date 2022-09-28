from sqlalchemy import create_engine, Table, Column
from pathlib import Path
import pandas as pd

db_type='mariadb'
id='admin'
pw="Geninus1!"
host='db-practice.cg10utv5zpvf.ap-northeast-2.rds.amazonaws.com'
schema_name='geninus'
url =f"{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}"

engine = create_engine(
    url
)


conn = engine.connect()

query = 'select * from gc_qc_sample order by create_date desc limit 100;'
data = pd.read_sql(query, conn)


data = pd.read_csv(Path())

data = data.rename(columns={name_a:name:b, ~~})

data.to_sql(conn, 'gh2')