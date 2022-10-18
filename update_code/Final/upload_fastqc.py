
import pandas as pd
from pathlib import Path
from datetime import datetime
from sqlalchemy import create_engine

class UploadFastQC:
    
    #dict가 pickle 파일
    def _split_fc_dir_to_date_and_id(self) -> str | str:
        tokens = self.fc_dir.split("_")
        seq_date = tokens[0]
        fc_id = tokens[-1]
        fc_id = fc_id[1:]
        return seq_date, fc_id

    def parse_df(self, sample_id: str, read_type: str, target_df: pd.DataFrame) -> pd.DataFrame:
        seq_date, fc_id = self._split_fc_dir_to_date_and_id()        
        target_df.insert(0, 'SEQ_DATE', seq_date)
        target_df.insert(1, 'FC_ID', fc_id)
        target_df.insert(2, 'SAMPLE_ID', sample_id)
        target_df.insert(3, 'TYPE', read_type)
        return target_df

    def upload_fastqc_dict(self, qc_dict: dict) -> None:
        target_key_dict = {'bsq':"gc_qbsq", 'bsc':"gc_qbsc", 'sge':"gc_qsge"}
        for sample_id, read_type in qc_dict.items():
            for key, table_name in target_key_dict.items():
                target_df = qc_dict[sample_id][read_type][key]
                target_df = self.parse_df(sample_id, read_type, target_df)
                target_df.to_sql(name=table_name, con=conn, if_exists='append', index=False)

    #csv로 나오는게 fastqc
    def parse_db_to_df(self, query: str) -> pd.DataFrame:
        db_type = 'mariadb'
        id = 'yckim'
        pw =  'rladmsco12!'
        host = 'geninus-maria-211117.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com'
        schema_name = 'gh2'
        url = f'{db_type}+pymysql://{id}:{pw}@{host}/{schema_name}'
        engine = create_engine(url)
        con = engine.connect()
        tb_expr_seq_line_df = pd.read_sql(query, con)
        tb_expr_seq_line_df = tb_expr_seq_line_df.rename(columns={'sample_id':'SampleID'})
        return tb_expr_seq_line_df
    
    def merged_qc_expr(self, qc_df: pd.DataFrame, tb_expr_seq_line_df: pd.DataFrame) -> pd.DataFrame:
        # qc_df['SampleID'] = [x.replace("'","") for x in qc_df['SampleID']]
        qc_merged = pd.merge(qc_df, tb_expr_seq_line_df, how='left', on='SampleID')
        qc_merged_df = qc_merged.rename(columns={'SampleID':'SAMPLE_ID'})
        qc_merged_df = qc_merged_df.rename(columns={'data_output':'DATA_OUTPUT'})
        qc_merged_df['FASTQ_OVERALL'] = ['PASS' if row['FASTQ_TOTAL_BASES(Gb)'] > row['data_output'] else 'FAIL' for _, row in qc_merged.iterrows()]
        qc_merged_df['FASTQ_TOTAL_READ'] = qc_merged_df['FASTQ_TOTAL_READ_R1'] + qc_merged_df['FASTQ_TOTAL_READ_R2']
        qc_merged_df['CREATE_DATE'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        return qc_merged_df
    
    def add_flowcell(self, qc_merged_df: pd.DataFrame) -> pd.DataFrame:
        #flow_cell = str(self.data_dir).split('/')[4].rsplit('_')[3]
        # flow_cell_id = data_dir.parent.name.split("_")[-1][1:]
        seq_date, fc_id = self._split_fc_dir_to_date_and_id()
        qc_merged_df.insert(0,'SEQ_DATE', seq_date)        
        qc_merged_df.insert(1,'FC_ID', fc_id)
        return qc_merged_df

    def upload_sql(self, conn: create_engine, qc_fastqc_upload_df: pd.DataFrame):
        qc_fastqc_upload_df.to_sql(name="qc_fastqc", con=conn, if_exists='append', index = False)

    def upload_fastqc_table(self, conn:create_engine, qc_fastqc_upload_df: pd.DataFrame) -> None:
        try:
            self.upload_sql(conn, qc_fastqc_upload_df)
        except Exception as e:
            self.logger.info("Error in upload sql")                       
            self.logger.info(e)
        return

    def __call__(self, qc_dict, qc_table) -> None:
        self.upload_fastqc_dict(qc_dict)

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
        qc_fastqc_upload_df = self.add_flowcell(qc_merged_df)
        self.upload_fastqc_table(qc_fastqc_upload_df)
        pass


