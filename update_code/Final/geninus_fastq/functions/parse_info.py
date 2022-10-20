
'''
    Date: 2022-07-20
    Author: duaghk
    
    Sequencing information parser

'''

# import library
import logging
import pandas as pd
from pathlib import Path
from pytz import timezone
from datetime import datetime
from sqlalchemy import create_engine, MetaData, Table, select

# custum config
from .configs import DBConfig, SequencingConfig

class ParseInfo(DBConfig, SequencingConfig):
    def __init__(self, date) -> None:
        DBConfig.__init__(self)
        SequencingConfig.__init__(self)
        connect_query = (
            f'{self.db_type}+pymysql://'
            f'{self.id}:'
            f'{self.pw}@'
            f'{self.db_address}/'
            f'{self.db_name}'
        )
        self.engine = create_engine(connect_query)
        self.conn = self.engine.connect()
        self.date = date
        self.log_dir = self.log_dir.joinpath('ParseInfo')
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self.make_logger(f"ParseInfo_{date}")
        self.logger.info(f'Sequencing Inforamtion parsing start. Target date: {self.date}')
    
    def make_logger(self, name=None, consoleset=True, streamset=True) -> logging:
        '''
            logger making.
        '''
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)
        loggerformat = '%(asctime)s - %(name)s - %(module)s - %(levelname)s - %(message)s'
        formatter = logging.Formatter(loggerformat)
        # change time zone
        def timetz(*args):
            return datetime.now(tz).timetuple() 
        tz = timezone('Asia/Seoul')
        logging.Formatter.converter = timetz
        if consoleset:
            console = logging.StreamHandler()
            console.setLevel(logging.INFO)
            console.setFormatter(formatter)
            logger.addHandler(console)
        if streamset:
            loggerfile = self.log_dir.joinpath(f'{name}.log')
            file_handler = logging.FileHandler(filename=loggerfile)
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        return logger
    
    # DataBase Parsing.
    def get_fcid(self) -> pd.DataFrame:
        '''
            Find fc_id by expr_dt
            Returns:
                expr_df: pd.DataFrame.
                    Columns(gh2):
                        expr_dt, run_id, equips, equip_side, fc_id
        '''
        meta = MetaData()
        expr_table = Table(self.fc_table, meta, autoload=True, autoload_with=self.engine)
        expr_columns = [expr_table.c[x] for x in expr_table.c.keys() if x in self.fc_table_columns]
        query = select(expr_columns).filter(expr_table.c['expr_dt'].like(f'{self.date}%'))
        expr_df = pd.read_sql_query(query, self.conn)
        # drop fc_id == None
        expr_df = expr_df[~expr_df['fc_id'].isna()]
        return expr_df
    
    def find_fc_dir(self, fc_id: str) -> str:
        '''
            Find flow cell dir for run bcl2fastq or ranger.
            sequencer directory are prefixed.
            Input:
                fc_id: str. flow cell id like HXX3R8XRG
            Returns:
                fc_dir: str, full path of flowcell dir.
        '''
        # sequencer dir is in configs['sequencer_dir']
        tgt_dir = [list(x.glob(f'*{fc_id}')) for x in self.sequence_dir_list]
        tgt_dir = [x for y in tgt_dir for x in y]
        return tgt_dir[0] if len(tgt_dir) else Path()

    def get_sample_list(self, run_id: str = None, fc_id: str = None) -> list:
        '''
            Get sample id using run_id.
            Bulk and Single cell are not splitted. use only one value.
            Input:
                run_id: str, like 21-0248
                fc_id: str, like HXX3R8XRG
            Returns:
                sample_id_list: list. 
        '''
        meta = MetaData()
        sample_table = Table(self.sample_table, meta, autoload=True, autoload_with=self.engine)
        sample_columns = [sample_table.c[x] for x in sample_table.c.keys() if x in self.sample_table_columns]
        if run_id:
            query = select(sample_columns).where(sample_table.c['run_id'] == run_id)
        elif fc_id:
            query = select(sample_columns).where(sample_table.c['fc_id'] == fc_id)
        sample_df = pd.read_sql_query(query, self.conn)
        sample_df = sample_df[~sample_df['sample_id'].isna()]
        sample_list = sample_df['sample_id'].tolist()
        return sample_list

    def get_sample_service_code(self, sample_list: list) -> pd.DataFrame:
        meta = MetaData()
        service_table = Table(self.service_table, meta, autoload=True, autoload_with=self.engine)
        service_columns = [service_table.c[x] for x in service_table.c.keys() if x in self.service_table_columns]
        query = select(service_columns, service_table.c['sample_id'].in_(sample_list))
        service_df = pd.read_sql_query(query, self.conn)
        return service_df

    def get_bulk_sample_idx(self, sample_list: list) -> pd.DataFrame:
        meta = MetaData()
        index_table = Table(self.bulk_idx_table, meta, autoload=True, autoload_with=self.engine)
        index_columns = [index_table.c[x] for x in index_table.c.keys() if x in self.bulk_idx_table_columns]
        query = select(index_columns, index_table.c['sample_id'].in_(sample_list))
        index_df = pd.read_sql_query(query, self.conn)
        index_df = index_df[self.bulk_idx_table_columns]
        return index_df

    def get_sc_sample_idx(self, sample_list: list) -> pd.DataFrame:
        meta = MetaData()
        index_table = Table(self.sc_idx_table, meta, autoload=True, autoload_with=self.engine)
        index_columns = [index_table.c[x] for x in index_table.c.keys() if x in self.sc_idx_table_columns]
        query = select(index_columns, index_table.c['sample_id'].in_(sample_list))
        index_df = pd.read_sql_query(query, self.conn)
        # filtering only need columns.
        index_df['Index'] = [f"{rows['index_plate_cd']}-{rows['index_no_cd']}" for _, rows in index_df.iterrows()]
        index_df = index_df[['sample_id', 'Index']]
        return index_df

    def get_fc_run_info(self, sequencer_fc_dir: Path) -> dict:
        '''
        get_flowcell_runinfo 

        Get flow cell run read information.
        e.g. read1:151 cycle, read2: 10 cycle, etc.

        Args:
            fc_dir (Path): fc_dir path. absolute path recommended

        Returns:
            dict: _description_
        '''

        # import library
        import xml.etree.ElementTree as ET  
        runinfo_path = sequencer_fc_dir.joinpath('RunInfo.xml')
        if runinfo_path.exists():
            # make xml tree.
            tree = ET.parse(runinfo_path)
            root = tree.getroot()
            # get total read length
            read_list = list(root.iter('Read'))
            read_dict = {x.attrib['Number']: ['Y' if x.attrib['IsIndexedRead'] == 'N' else 'I', x.attrib['NumCycles']] for x in read_list}
        else:
            read_dict = {}
        return read_dict

    def get_single_fc_result(self, expr_rows: pd.Series) -> dict:
        return_dict = {}
        return_dict['sequencer_fc_dir'] = self.find_fc_dir(expr_rows['fc_id'])
        sample_list = self.get_sample_list(run_id=expr_rows['run_id'])
        service_df = self.get_sample_service_code(sample_list)
        bulk_index_df = self.get_bulk_sample_idx(sample_list)
        sc_index_df = self.get_sc_sample_idx(sample_list)
        # merging with service code.
        bulk_df = pd.merge(bulk_index_df, service_df, on='sample_id', how='inner')
        sc_df = pd.merge(sc_index_df, service_df, on='sample_id', how='inner')
        return_dict['bulk_sample_df'] = bulk_df
        return_dict['sc_sample_df'] = sc_df
        return_dict['runinfo'] = self.get_fc_run_info(return_dict['sequencer_fc_dir'])
        return return_dict

    def __call__(self):
        expr_df = self.get_fcid()
        return_dict = {row['run_id']: self.get_single_fc_result(row) for _, row in expr_df.iterrows()}
        self.conn.close()
        return return_dict


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--date", help='Wanted date. format like 2022-01-01')
    args = parser.parse_args()
    db_parser = ParseInfo(args.date)
    db_parsed = db_parser()
    # print(db_parsed)
    for key, values in db_parsed.items():
        if values['sequencer_fc_dir'].name != '':
            save_file_dir = Path(db_parser.fastq_dir, values['sequencer_fc_dir'].name, 'samplesheet')
            save_file_dir.mkdir(parents=True, exist_ok=True)
            sequencer_info_file_fn = save_file_dir.joinpath(f"{key}_sequencer_info.txt")
            bulk_sample_info_file_fn = save_file_dir.joinpath(f"{key}_bulk_sample_info.csv")
            sc_sample_info_file_fn = save_file_dir.joinpath(f"{key}_sc_sample_info.csv")
            with open(sequencer_info_file_fn, 'w') as f:
                f.write(f"flowcell_dir: {values['sequencer_fc_dir']}\n")
                runinfo_line = [f"{v[0]}{v[1]}" for _, v in values['runinfo'].items()]
                f.write(f"RunInfo: {','.join(runinfo_line)}\n")
            values['bulk_sample_df'].to_csv(bulk_sample_info_file_fn, mode='a', header=True, index=False)
            values['sc_sample_df'].to_csv(sc_sample_info_file_fn, mode='a', header=True, index=False)
        else:
            pass
    






