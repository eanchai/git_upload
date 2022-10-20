#!env /usr/bin/python3

'''

Date: 2021.11.11
Authors: duaghk

Config file. for setting default arguments.

'''

from pathlib import Path

class DBConfig:
    def __init__(self) -> None:
        self.db_type = 'mariadb'
        self.id = 'pipeline'
        self.pw = 'vk2vmfkd!s'
        self.db_address = 'geninus-maria-211117.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com'
        self.db_name = 'gh2'
        self.fc_table = 'tb_expr_seq_header'
        self.fc_table_columns = ['expr_dt', 'run_id', 'equip', 'equip_side', 'fc_id'] 
        self.sample_table = 'tb_expr_seq_line'
        self.sample_table_columns = ['run_id', 'sample_id']
        self.bulk_idx_table = 'tb_expr_pre_pcr'
        self.bulk_idx_table_columns = ['sample_id', 'i5_index', 'i7_index']
        self.sc_idx_table = 'tb_expr_sc_lib'
        self.sc_idx_table_columns = ['sample_id', 'index_plate_cd', 'index_no_cd']
        self.service_table = 'tb_order_line'
        self.service_table_columns = ['sample_id', 'service_cd']

class ToolConfig:
    def __init__(self, tool_dir: str = '/ess/research_operation/bin') -> None:
        # set tool path
        self.tool_dir = Path(tool_dir)
        self.bcl2fastq = self.tool_dir.joinpath('bcl2fastq_v2.20', 'bin', 'bcl2fastq')
        self.cellranger = self.tool_dir.joinpath('cellranger-6.1.2', 'bin', 'cellranger')
        self.spaceranger = self.tool_dir.joinpath('spaceranger-1.3.1', 'bin', 'spaceranger')
        self.fastqc = self.tool_dir.joinpath('FastQC', 'fastqc')
        # set tool threads
        self.bcl_threads = 16
        self.fastqc_threads = 2

class SequencingConfig:
    def __init__(self) -> None:
        self.main_dir = Path('/ess/NGS_data')
        self.sequence_dir_list = [self.main_dir.joinpath('NovaSeq_data'), self.main_dir.joinpath('NextSeq_data')]
        self.fastq_dir = self.main_dir.joinpath('fastq')
        self.log_dir = self.fastq_dir.joinpath('log')

class MailingConfig:
    def __init__(self) -> None:
        self.login_id = 'research-operation@kr-geninus.com'
        self.login_pw = 'ROuser!2341234'
        self.from_mail = 'research-operation@kr-geninus.com'
        self.to_mail_list = [
            # 'duaghk@kr-geninus.com',
            # 'development@kr-geninus.com',
            # 'research-operation@kr-geninus.com',
            'genome_center@kr-geninus.com'
            ]

class UploadConfig:
    def __init__(self) -> None:
        self.profile = 'fastq_uploader'
        self.upload_bucket="geninus-ro-fastq"
        pass