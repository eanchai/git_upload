
from asyncio.windows_events import NULL
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

db_type = 'mariadb'
id = 'yckim'
pw =  'rladmsco12!'
host = 'geninus-maria-211117.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com'
schema_name = 'gh2'
url = f'{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}'
engine = create_engine(url)
conn = engine.connect()
        #data merge
        
qc_df = pd.read_csv("/home/ubuntu/NGS_data/fastqc/fastqc_product_results_dual_10.csv")
qc_df['SampleID'] = [f"'{x}'" for x in qc_df['SampleID']] 
query = f'''
select sample_id, data_output from tb_expr_seq_line
where sample_id in ({",".join(qc_df['SampleID'].tolist())})
'''
tb_expr_seq_line = pd.read_sql(query, conn)
qc_df['SampleID'] = [x.replace("'","") for x in qc_df['SampleID']] 
qc_merged = pd.merge(qc_df, tb_expr_seq_line, how='left', left_on='SampleID', right_on='sample_id')
qc_merged['FILTER'] = ['PASS' if row['FASTQ_TOTAL_BASES(Gb)'] > row['data_output'] else 'FAIL' for _, row in qc_merged.iterrows()]
qc_merged.drop(columns='sample_id')

        db_type = 'mariadb'
        id = 'root'
        pw =  'gw12341234'
        host = 'db-qc.cg10utv5zpvf.ap-northeast-2.rds.amazonaws.com'
        schema_name = 'geninus'
        url = f'{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}'
        engine = create_engine(url)
        meta = MetaData(bind=engine)
        MetaData.reflect(meta)
        mytable = Table('gc_qc_bi', meta)
        meta.create_all(engine)
        ## query
        for _, row in qc_merged.drop(columns='sample_id').iterrows():
             query = mytable.insert().values(SAMPLE_ID = {row['SampleID']},
                FASTQ_TOTAL_READ = {row['FASTQ_TOTAL_READ_R1'] + row['FASTQ_TOTAL_READ_R2']},
                FASTQ_Q30 = None,
                FASTQ_OVERALL = {row['FILTER']},
                FASTQ_TOTAL_READ_R1 = {row['FASTQ_TOTAL_READ_R1']},
                FASTQ_GC_CONTENTS_R1 = None ,
                FASTQ_OVERALL_R1 = None,
                FASTQ_Q30_R1 = {row['FASTQ_Q30_R1']},
                FASTQ_TOTAL_READ_R2 = {row['FASTQ_TOTAL_READ_R2']},
                FASTQ_GC_CONTENTS_R2 = None,
                FASTQ_OVERALL_R2 = None,
                FASTQ_Q30_R2 = {row['FASTQ_Q30_R2']},
                BAM_MEAN_DEPTH = None,
                BAM_CAPTURE_EFFICIENCY = None,
                BAM_ON_TARGET_RATE = None,
                BAM_DUP_RATE = None,
                BAM_PR_SCORE = None,
                BAM_UNIFORM = None,
                BAM_TUMOR_VOL = None,
                BAM_OVERALL_QUALITY = None,
                BAM_GENDER_ESTIMATED = None,
                ANAL_PF_QC = None,
                ANAL_PF_COMM = None,
                PCT_TARGET_04X_MEAN = None,
                ALIGN_RATE = None,
                UNIQUE_READ = None,
                UNIQUE_MEAN_DEPTH = None,
                BAM_MEAN_INSERT_SIZE = None,
                COEF = None,
                MEDIAN_EXON_COVERAGE = None
                # LAST_UPDATE_DATE = None,
                # CREATE_USER = '',
                # LAST_UPDATE_USER = 'root'
             )

table_cols = [x.name for x in mytable.columns]
not_in_table_cols = [x for x in table_cols if x not in qc_merged.columns]
for col in not_in_table_cols:
        qc_merged[col] = None

qc_merged = qc_merged[table_cols]
qc_merged = qc_merged.reindex(columns=table_cols)
qc_merged.to_sql()


datetime.now().strftime("%Y-%m-%d %H:%M:%S")





d = list(mytable.columns)
d[0].name





for _, row in qc_merged.drop(columns='sample_id').iterrows():
        gc_qc_bidf = pd.DataFrame({
                'SAMPLE_ID' : [{row['SampleID']}],
                'FASTQ_TOTAL_READ' : {row['FASTQ_TOTAL_READ_R1'] + row['FASTQ_TOTAL_READ_R2']},
                'FASTQ_Q30' : None,
                'FASTQ_OVERALL' : {row['FILTER']},
                'FASTQ_TOTAL_READ_R1' : {row['FASTQ_TOTAL_READ_R1']},
                'FASTQ_GC_CONTENTS_R1' : None ,
                'FASTQ_OVERALL_R1' : None,
                'FASTQ_Q30_R1' : {row['FASTQ_Q30_R1']},
                'FASTQ_TOTAL_READ_R2' : {row['FASTQ_TOTAL_READ_R2']},
                'FASTQ_GC_CONTENTS_R2' : None,
                'FASTQ_OVERALL_R2' : None,
                'FASTQ_Q30_R2' : {row['FASTQ_Q30_R2']},
                'BAM_MEAN_DEPTH' : None,
                'BAM_CAPTURE_EFFICIENCY ': None,
                'BAM_ON_TARGET_RATE' : None,
                'BAM_DUP_RATE' : None,
                'BAM_PR_SCORE' : None,
                'BAM_UNIFORM' : None,
                'BAM_TUMOR_VOL' : None,
                'BAM_OVERALL_QUALITY' : None,
                'BAM_GENDER_ESTIMATED' : None,
                'ANAL_PF_QC' : None,
                'ANAL_PF_COMM' : None,
                'PCT_TARGET_04X_MEAN' : None,
                'ALIGN_RATE' : None,
                'UNIQUE_READ' : None,
                'UNIQUE_MEAN_DEPTH' : None,
                'BAM_MEAN_INSERT_SIZE' : None,
                'COEF' : None,
                'MEDIAN_EXON_COVERAGE' : None
                # LAST_UPDATE_DATE : None,
                # CREATE_USER : '',
                # LAST_UPDATE_USER : 'root'
        })
        
        
sql_rows = []

sql_rows = '({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
 {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})'.format()