'''
    date: 2022.10.12 ~
    object: compare expr_line & update new Table(gc_fastqc)
    author:
'''
    
import logging
import sqlalchemy
import pandas as pd
# import sqlalchemy as db
# from time import time
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

#db connect_geninus dump _ 연습용 서버
class DBConfig:
    def __init__(self):
        self.db_type = 'mariadb'
        self.id = 'root'
        self.pw =  'gw12341234'
        self.host = 'db-qc.cg10utv5zpvf.ap-northeast-2.rds.amazonaws.com'
        self.schema_name = 'geninus'
        
class FastQCUpload(DBConfig):
    def __init__(self, data_dir: str, type_len: str):
        '''
        Class initialize
        '''
        DBConfig.__init__(self)
        self.file_format = ".csv"
        self.type_len = type_len
        self.data_dir = Path(data_dir)
        common = "fastqc_product_results_"
        #fastqc파일 qc_df에 저장
        self.qc_df = self.data_dir.joinpath(f"{common}{self.type_len}{self.file_format}")
        #logger파일
        self.log_dir = self.data_dir.joinpath("log")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self._make_logger(self.log_dir, 'FastQC_Upload')
         
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
        
    #database connect 
    def connect_db(self) -> sqlalchemy.engine:
        url = f'{self.db_type}+pymysql://{self.id}:{self.pw}@{self.host}/{self.schema_name}'
        engine = create_engine(url)
        conn = engine.connect()
        return conn
    
    #tb_expr_seq_line 불러와서 output 비교
    # 원래 아래 코드로 gh2에서 돌려야하는데 지금 geninus폴더에 expr이 없음
    # def parse_db_to_df(self, conn, query: str) -> pd.DataFrame:
    #     tb_expr_seq_line_df = pd.read_sql(query, conn)
    #     tb_expr_seq_line_df = tb_expr_seq_line_df.rename(columns={'type_len':'SampleID'})
    #     return tb_expr_seq_line_df
    
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
        qc_df['SampleID'] = [x.replace("'","") for x in qc_df['SampleID']]
        qc_merged = pd.merge(qc_df, tb_expr_seq_line_df, how='left', on='SampleID')
        qc_merged_df = qc_merged.rename(columns={'SampleID':'SAMPLE_ID'})
        qc_merged_df = qc_merged_df.rename(columns={'data_output':'DATA_OUTPUT'})
        qc_merged_df['FASTQ_OVERALL'] = ['PASS' if row['FASTQ_TOTAL_BASES(Gb)'] > row['data_output'] else 'FAIL' for _, row in qc_merged.iterrows()]
        qc_merged_df['FASTQ_TOTAL_READ'] = qc_merged_df['FASTQ_TOTAL_READ_R1'] + qc_merged_df['FASTQ_TOTAL_READ_R2']
        qc_merged_df['CREATE_DATE'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        return qc_merged_df
    
    def add_flowcell(self, data_dir: Path, qc_merged_df: pd.DataFrame) -> pd.DataFrame:
        #flow_cell = str(self.data_dir).split('/')[4].rsplit('_')[3]
        flow_cell_id = data_dir.parent.name.split("_")[-1][1:]
        qc_merged_df.insert(1,'FLOW_CELL_ID',flow_cell_id)
        return qc_merged_df
    
    def upload_sql(self, conn: sqlalchemy.engine, qc_fastqc_upload_df: pd.DataFrame):
        qc_fastqc_upload_df.to_sql(name="qc_fastqc", con=conn, if_exists='append', index = False)
    
    def __call__(self):
        self.logger.info("Start fastqc_overall update pipeline")
        self.logger.info("Connect db...")
        conn = self.connect_db()
        data_dir = self.data_dir
        qc_df = pd.read_csv(self.qc_df)
        qc_df['SampleID'] = [f"'{x}'" for x in qc_df['SampleID']] 

        # load QC df 
        # # load tb_expr_seq_line from database
        self.logger.info("Parsing the data...")
        query = f'''
            select sample_id, data_output from tb_expr_seq_line
            where sample_id in ({",".join(qc_df['SampleID'].tolist())})
            '''
        tb_expr_seq_line_df = self.parse_db_to_df(query)
        
        self.logger.info("Filtering Data with threshold...")
        qc_merged_df = self.merged_qc_expr(qc_df, tb_expr_seq_line_df)
        
        self.logger.info("Add flow cell ID on db...")
        qc_fastqc_upload_df = self.add_flowcell(data_dir, qc_merged_df)
        
        self.logger.info("Upload to Database_gc_fastqc...")
        try:
            self.upload_sql(conn, qc_fastqc_upload_df)
        except Exception as e:
            self.logger.info("Error in upload sql")                       
            self.logger.info(e)
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", '-d')
    parser.add_argument("--type_len", '-sid')
    args = parser.parse_args()
    runner = FastQCUpload(args.data_dir, args.type_len)
    runner()