#!env /usr/bin/python3

'''

Date: 2021.11.11
Authors: duaghk

Config file. for setting default arguments.

'''
configs = {
    # set db configs.
    # "gh1": {
    #     "db_type": "mysql",
    #     "flowcell_table": "",
    #     "flowcell_table_columns": [],
    #     "sample_table": "",
    #     "sample_table_columns": [],
    #     "index_table": "",
    #     "index_table_columns": []
    # },
    "gh2": {
        "id": "pipeline",
        "pw": "vk2vmfkd!s",
        "db_address": "geninus-maria-211117.cobyqiuirug6.ap-northeast-2.rds.amazonaws.com",
        "db_name": "gh2", 
        "db_type": "mariadb",
        "flowcell_table": "tb_expr_seq_header",
        "flowcell_table_columns": ["expr_dt", "run_id", "equip", "equip_side", "fc_id"],
        "sample_table": "tb_expr_seq_line",
        "sample_table_columns": ['run_id', 'sample_id'],
        "bulk_index_table": "tb_expr_pre_pcr",
        "bulk_index_table_columns": ["sample_id", "i5_index", "i7_index"],
        "sc_index_table": "tb_expr_sc_lib",
        "sc_index_table_columns": ["sample_id", "index_plate_cd", "index_no_cd"],
        "service_table": "tb_order_line",
        "service_table_columns": ["sample_id", "service_cd"]
    },
    # tool path setting.
    "bcl2fastq": "/ess/research_operation/bin/bcl2fastq_v2.20/bin/bcl2fastq",
    "cellranger": '/ess/research_operation/bin/cellranger-6.1.2/bin/cellranger',
    "spaceranger": '/ess/research_operation/bin/spaceranger-1.3.1/bin/spaceranger',
    "fastqc": '/ess/research_operation/bin/FastQC/fastqc',
    # sequencer target dir setting.
    "sequencer_dir": [
        # '/data2/NGS_data/NextSeq_data', '/data2/NGS_data/NovaSeq_data/Runs/MyRun',
        # '/data6/NGS_data/NextSeq_data', '/data6/NGS_data/NovaSeq_data/Runs', 
        '/ess/NGS_data/NextSeq_data', '/ess/NGS_data/NovaSeq_data'
        ],
    # set threads
    "bcl_threads": 16,
    "fastqc_threads": 2,
    # set mailing configs.
    "mailing": {
        "login_id": 'research_operation@kr-geninus.com',
        "login_pw": "ROuser!2341234",
        "from_mail": 'research_operation@kr-geninus.com',
        "to_mail": [
            "development@kr-geninus.com",
            "research-operation@kr-geninus.com"
            "mj.park@kr-geninus.com"
            ]
    }
}
