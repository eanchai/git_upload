#!env /usr/bin/python3

'''

Date: 2021.11.11
Authors: duaghk

Parse DB and parse sample information by index-type

'''

# import library
#!env /usr/bin/python3

'''

Date: 2021.11.11
Authors: duaghk

Parse DB and parse sample information by index-type

'''

# import library

import os
import pandas as pd
from sqlalchemy import create_engine, MetaData, Table, select

class ParseDB:
    '''
        Connect GeniHub database and make samplesheet
        Input: 
            Configs for db info.
            date for find only after date.
            db_version: Genihub1 or Genihub2 version.
        Returns:
            fc_dir with samplesheet.
            dictionary format?
            dict:
                fc_id:
                    fc_dir
                    samplesheet
        
    '''
    def __init__(self, configs: dict, date: str, db_version: str, logger=None):
        self.configs = configs
        self.db_info = configs[db_version]
        connect_query = (
            f"{self.db_info['db_type']}+pymysql://"
            f"{self.db_info['id']}:"
            f"{self.db_info['pw']}@"
            f"{self.db_info['db_address']}/"
            f"{self.db_info['db_name']}"
        )
        self.engine = create_engine(connect_query)
        self.conn = self.engine.connect()
        self.date = date
        self.logger = logger
    
    def get_fcid(self) -> list:
        '''
            Find fc_id by expr_dt
            Returns:
                expr_df: pd.DataFrame.
                    Columns(gh2):
                        expr_dt, run_id, equips, equip_side, fc_id
        '''
        meta = MetaData()
        expr_table = Table(self.db_info["flowcell_table"], meta, autoload=True, autoload_with=self.engine)
        expr_columns = [expr_table.c[x] for x in expr_table.c.keys() if x in self.db_info["flowcell_table_columns"]]
        query = select(expr_columns).where(expr_table.c['expr_dt'] >= self.date)
        expr_df = pd.read_sql_query(query, self.conn)
        # drop fc_id == None
        expr_df = expr_df[~expr_df["fc_id"].isna()]
        # self.logger.info(f"{len(expr_df)} flowcell id found.")
        return expr_df
    
    def find_fc_dir(self, fc_id) -> str:
        # sequencer dir is in configs['sequencer_dir']
        for dir in self.configs['sequencer_dir']:
            tgt_dir = [x for x in os.listdir(dir) if x.endswith(fc_id)]
            if len(tgt_dir):
                yield tgt_dir[0]

    def get_sample_id(self, run_id: str = None, fc_id: str = None) -> list:
        meta = MetaData()
        sample_table = Table(self.db_info["sample_table"], meta, autoload=True, autoload_with=self.engine)
        sample_columns = [sample_table.c[x] for x in sample_table.c.keys() if x in self.db_info["sample_table_columns"]]
        if run_id:
            query = select(sample_columns).where(sample_table.c["run_id"] == run_id)
        elif fc_id:
            query = select(sample_columns).where(sample_table.c["fc_id"] == fc_id)
        sample_df = pd.read_sql_query(query, self.conn)
        sample_df = sample_df[~sample_df["sample_id"].isna()]
        # if run_id:
        #     self.logger.info(f"Found {len(sample_df)} samples in run_id: {run_id}")
        # elif fc_id:
        #     self.logger.info(f"Found {len(sample_df)} samples in fc_id: {fc_id}")
        sample_list = sample_df["sample_id"].tolist()
        return sample_list

    def get_sample_service_code(self, sample_list: list) -> pd.DataFrame:
        meta = MetaData()
        service_table = Table(self.db_info["service_table"], meta, autoload=True, autoload_with=self.engine)
        service_columns = [service_table.c[x] for x in service_table.c.keys() if x in self.db_info["service_table_columns"]]
        query = select(service_columns, service_table.c["sample_id"]).in_(sample_list)
        service_df = pd.read_sql_query(query, self.conn)
        # now join.
        return service_df

    def get_bulk_sample_index(self, sample_list: list) -> pd.DataFrame:
        meta = MetaData()
        index_table = Table(self.db_info["bulk_index_table"], meta, autoload=True, autoload_with=self.engine)
        index_columns = [index_table.c[x] for x in index_table.c.keys() if x in self.db_info["bulk_index_table_columns"]]
        query = select(index_columns, index_table.c["sample_id"].in_(sample_list))
        # if run_id:
        #     query = query.where(index_table.c["run_id"] == run_id)
        # elif fc_id:
        #     query = query.where(index_table.c["fc_id"] == fc_id)
        index_df = pd.read_sql_query(query, self.conn)
        # filtering only need columns.
        index_df = index_df[["sample_id", "i5_index", "i7_index"]]
        return index_df

    def get_sc_sample_index(self, sample_list: list) -> pd.DataFrame:
        meta = MetaData()
        index_table = Table(self.db_info["sc_index_table"], meta, autoload=True, autoload_with=self.engine)
        index_columns = [index_table.c[x] for x in index_table.c.keys() if x in self.db_info["sc_index_table_columns"]]
        query = select(index_columns, index_table.c["sample_id"].in_(sample_list))
        # if run_id:
        #     query = query.where(index_table.c["run_id"] == run_id)
        # elif fc_id:
        #     query = query.where(index_table.c["fc_id"] == fc_id)
        index_df = pd.read_sql_query(query, self.conn)
        # filtering only need columns.
        index_df["Index"] = index_df.apply(lambda x: f"{x['index_plate_cd']}-{x['index-no_cd']}")
        index_df = index_df[["sample_id", "index_plate_cd", "index_no_cd"]]
        return index_df


    def get_single_fc_result(self, expr_rows: pd.Series) -> dict:
        returndict = {}
        # returndict["fc_dir"] = list(self.find_fc_dir(expr_rows["fc_id"]))
        sample_list = self.get_sample_id(expr_rows["run_id"])
        index_df = self.get_sample_index(sample_list)
        sample_df = self.get_sample_service_code(index_df)
        returndict["sample_df"] = sample_df
        return returndict

    # need to add status update query.

    def __call__(self):
        expr_df = self.get_fcid()
        returndict = {row["run_id"]: self.get_single_fc_result(row) for _, row in expr_df.iterrows()}
        return returndict


import logging
from pytz import timezone
from datetime import datetime, timedelta
from configs import configs

def make_logger(log_dir, name=None, consoleset=True, streamset=True) -> logging:
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
        loggerfile = os.path.join(log_dir, name)
        file_handler = logging.FileHandler(filename=f"{loggerfile}.log")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    return logger
log_dir = os.path.join("./", "log")
os.makedirs(log_dir, exist_ok=True)
logger = make_logger(log_dir, "ParseDB")
curr_dir = '2021-11-29'
curr_date = datetime.now(timezone("Asia/Seoul"))
before_7 = curr_date - timedelta(days=30)
curr_date = curr_date.strftime("%Y-%m-%d")
before_7 = before_7.strftime("%Y-%m-%d")
db_parser = ParseDB(configs, before_7, "gh2", logger)

fc_id = db_parser.get_fcid()
fc_id

# get target fc_id
test_fc_id = "21-0001"
samplelist = db_parser.get_sample_id(run_id = test_fc_id)
samplelist

index_df = db_parser.get_bulk_sample_index(samplelist)
index_df

index_df2 = db_parser.get_sc_sample_index(samplelist)
index_df2


meta = MetaData()
index_table = Table(db_parser.db_info["sc_index_table"], meta, autoload=True, autoload_with=db_parser.engine)
index_columns = [index_table.c[x] for x in index_table.c.keys() if x in db_parser.db_info["sc_index_table_columns"]]
query = select(index_columns, index_table.c["sample_id"].in_(samplelist))
# if run_id:
#     query = query.where(index_table.c["run_id"] == run_id)
# elif fc_id:
#     query = query.where(index_table.c["fc_id"] == fc_id)
index_df = pd.read_sql_query(query, db_parser.conn)
# filtering only need columns.
index_df
index_df["Index"] = [f"{rows['index_plate_cd']}-{rows['index_no_cd']}" for _, rows in index_df.iterrows()]

index_df = index_df[["sample_id", "index_plate_cd", "index_no_cd"]]




configs

