#!env /usr/bin/python3

'''

Date: 2021.11.11
Authors: duaghk

Main script for integrate fastq making.

'''

# import library

import sys
import logging
from time import time, sleep
from pytz import timezone
from datetime import datetime, timedelta
from pathlib import Path
import subprocess as sp

from .functions import ParseInfo, MakeSampleSheet, Bcl2fastqRunner, RangerRunner, FastQCRunner, ParseFastQC, QueueController, SendEmail, update_db, UploadFastq, ToolConfig, SequencingConfig


class MakeFastq(ToolConfig, SequencingConfig):
    '''
        Module for BCL2FASTQ running.
    '''
    def __init__(self, target_date: str, overwrite: bool = False):
        '''
            Initialize
            Input:
                outdir: Main dir for run process.
                    e.g.: /data6/NGS_data/fastq <- here
                sge: Using queue instead of process
            Config file:
                value:
                    bcl2fastq: path
                    cellranger: path
                    spaceranger: path
                    threads: int
        '''
        # date time setting.
        ToolConfig.__init__(self)
        SequencingConfig.__init__(self)
        input_date = datetime.strptime(target_date, "%Y-%m-%d")
        # before_one = input_date - timedelta(days=3)
        self.date = input_date.strftime("%Y-%m-%d")
        # self.date = target_date
        self.log_dir = self.log_dir.joinpath("GeninusFastq")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self.make_logger(target_date)
        self.overwrite = overwrite

    def make_logger(self, name=None, consoleset=True, streamset=True) -> logging:
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
            loggerfile = self.log_dir.joinpath(f"{name}.log")
            file_handler = logging.FileHandler(filename=loggerfile)
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        return logger

    def _mixed_check(self, samplesheet_list):
        idx_type_list = list(set([x.name.split("_")[1] for x in samplesheet_list]))
        return len(idx_type_list) > 1
    
    def _make_dir(self, fastq_fc_output_dir):
        '''
            Make dir for running.
            Need directories:
                fastq_dir, script_dir, samplesheet_dir, fastqc_dir, log_dir.
        '''
        # dir_name_list = ["fastq", "script", "samplesheet", "fastqc", "log"]
        # Not use for loop. Because need to return pinpoint directory.

        # make fastq dir
        fastq_dir = fastq_fc_output_dir.joinpath("fastq")
        fastq_dir.mkdir(parents=True, exist_ok=True, mode=0o777)
        # make samplesheet dir
        samplesheet_dir = fastq_fc_output_dir.joinpath("samplesheet")
        samplesheet_dir.mkdir(parents=True, exist_ok=True, mode=0o777)
        # make fastqc dir
        fastqc_dir = fastq_fc_output_dir.joinpath("fastqc")
        fastqc_dir.mkdir(parents=True, exist_ok=True, mode=0o777)
        # make script dir
        script_dir = fastq_fc_output_dir.joinpath("script")
        script_dir.mkdir(parents=True, exist_ok=True, mode=0o777)

        return fastq_dir, samplesheet_dir, fastqc_dir, script_dir

    @staticmethod
    def _check_and_sleep_until_RTA_complete(sequencer_fc_dir):
        while True:
            filelist = [x.name for x in sequencer_fc_dir.glob("*Complete.txt")]
            if "RTAComplete.txt" in filelist and "CopyComplete.txt" in filelist:
                break
            else:
                sleep(300)
        pass

    def run_single_fc(self, info: dict):
        '''
        run_single_fc 

        Single flow cell info batch run.

        Args:
            info (dict): Information for flowcell id.
                Keys:
                    fc_dir
                    bulk_sample_df
                    single_sample_df
        '''
        # # check if fc_dir exists.
        # if not info['fc_dir']:
        #     self.logger.info(f"{run_id} has no FlowCell directory!")
        #     sys.exit()
        # # if not len(info['bulk_sample_df']) and not len(info['sc_sample_df']):
        # if not len(info['bulk_sample_df']):
        #     self.logger.info(f"{run_id} has no Bulk Samplesheet!")
        #     sys.exit()
        # make outdir path
        sequencer_fc_dir = info['sequencer_fc_dir']
        self._check_and_sleep_until_RTA_complete(sequencer_fc_dir)
        fc_name = sequencer_fc_dir.name
        fastq_fc_output_dir = self.fastq_dir.joinpath(fc_name)
        fastq_fc_output_dir.mkdir(parents=True, exist_ok=True, mode=0o777)
        # make directories.
        fastq_dir, samplesheet_dir, fastqc_dir, script_dir = self._make_dir(fastq_fc_output_dir)
        # make logger for flowcell.
        # fc_logger = self.make_logger(log_dir, name=fc_name)
        # make sample sheet
        self.logger.info("Make SampleSheet...")
        samplesheet_list = self.make_ss(
            samplesheet_dir=samplesheet_dir,
            bulk_sample_df=info['bulk_sample_df'], 
            sc_sample_df=info['sc_sample_df'],
            )
        mixed = self._mixed_check(samplesheet_list)  # dual,single mixed sequencing check.
        # get bcl2fastq commands about each samplesheet.
        fastq_run_dict = {}
        for samplesheet_path in samplesheet_list:
            if samplesheet_path.name.split("_")[1] in ["chromium", "visium"]:
                # cmd, fastq_outdir = self.ranger(fc_dir, fastq_dir, samplesheet_path, fc_logger)
                pass
            else:
                print(info['runinfo'])
                bases_mask = ",".join([f"{v[0]}{v[1]}" for _,v in info['runinfo'].items()])
                cmd, fastq_fc_output_dir = self.bcl2fastq(sequencer_fc_dir, bases_mask, samplesheet_path, mixed)
                idx_type = fastq_fc_output_dir.name
            fastq_run_dict[idx_type] = {'cmd':cmd, 'outdir': fastq_fc_output_dir}
        # overwrite check.
        if self.overwrite:
            for k, v in fastq_run_dict.items():
                fastq_fc_output_dir = v['outdir']
                if fastq_fc_output_dir.exists():
                    sp.run(f"rm -rf {v['outdir']}", shell=True)
                else:
                    pass
        else:
            for k, v in fastq_run_dict.items():
                fastq_fc_output_dir = v['outdir']
                if fastq_fc_output_dir.exists():
                    fastq_list = [x for x in fastq_fc_output_dir.glob("*.fastq.gz")]
                    if len(fastq_list):
                        fastq_run_dict[k]['cmd']= 'echo ""'

        # need to generate fastq script dir.
        fastq_script_dir = script_dir.joinpath("fastq")
        fastq_script_dir.mkdir(parents=True, exist_ok=True, mode=0o777)
        self.queuecontroller(fastq_run_dict, fastq_script_dir, self.bcl_threads)
        # need to remove Undetermined.
        for idx_type, v in fastq_run_dict.items():
            fastq_fc_output_dir = v['outdir']
            undetermined_r1 = fastq_fc_output_dir.joinpath("Undetermined_S0_R1_001.fastq.gz")
            undetermined_r2 = fastq_fc_output_dir.joinpath("Undetermined_S0_R2_001.fastq.gz")
            undetermined_r3 = fastq_fc_output_dir.joinpath("Undetermined_S0_R3_001.fastq.gz")
            if undetermined_r1.exists(): sp.run(f"rm -rf {undetermined_r1}", shell=True)
            if undetermined_r2.exists(): sp.run(f"rm -rf {undetermined_r2}", shell=True)
            if undetermined_r3.exists(): sp.run(f"rm -rf {undetermined_r3}", shell=True)

        # QC check.
        # need to generate fastqc script dir.
        fastqc_script_dir = script_dir.joinpath("fastqc")
        fastqc_script_dir.mkdir(parents=True, exist_ok=True, mode=0o777)
        fastqc_cmd_info_dict = {}
        for idx_type, v in fastq_run_dict.items():
            tmpoutdir = fastqc_dir.joinpath(idx_type)
            tmpoutdir.mkdir(parents=True, exist_ok=True)
            tmpdict = self.qcrunner(
                v['outdir'], 
                tmpoutdir,
                idx_type
                )
            fastqc_cmd_info_dict.update(tmpdict)
        # fastqc_cmd_info_dict = self.qcrunner(fastq_dir, fastqc_dir, fc_logger)
        self.queuecontroller(fastqc_cmd_info_dict, fastqc_script_dir, self.fastqc_threads)
        # QC parsing.
        fastqc_result_path_list = []
        for samplesheet_prefix in fastqc_dir.glob("*"):
            if samplesheet_prefix.is_dir():
                qc_dict , qc_table = self.qcparser(fastqc_dir.joinpath(samplesheet_prefix))
                fastqc_result_path = fastqc_dir.joinpath(f'fastqc_product_results_{samplesheet_prefix.name}.csv')
                
                # qc_table.to_csv(fastqc_result_path, header=True, index=False)
                # Need to add Upload data.
                self.updater(fastqc_result_path, qc_dict, qc_table)
                
                fastqc_result_path_list.append(fastqc_result_path)
            
            else:
                continue
        
        # Need to add mailing.
        headline = f'{self.date} {fc_name} BCL2FASTQ completed.'
        main_text = f"""
        {self.date}  {fc_name}'s BCL2FASTQ completed. \n
        Please check your samples in attached files.
        """

        self.mailer(headline, main_text, fastqc_result_path_list)
        return fc_name

    def _initialize_module(self):
        try:
            self.parse_db = ParseInfo(self.date)
            self.logger.info("Database parse module initialized.")
            self.make_ss = MakeSampleSheet(self.date)
            self.logger.info("Samplesheet parse module initialized.")
            self.bcl2fastq = Bcl2fastqRunner(self.date)
            self.logger.info("Bulk seq fastq make module initialized.")
            # init ranger
            # self.ranger = RangerRunner(
            #     self.configs["cellranger"],
            #     self.configs["spaceranger"],
            #     self.configs["bcl_threads"],
            #     )
            # self.logger.info("Single cell fastq make module initialized.")
            # init check_quality
            self.qcrunner = FastQCRunner(self.date)
            self.logger.info("Quality check module initialized.")
            # init queue controller.
            self.queuecontroller = QueueController(self.date)
            self.logger.info("Queue control module initalized.")
            # init Fastq output parser.
            self.qcparser = ParseFastQC(self.date)
            self.logger.info("FASTQC result parse module initalized")
            # init mailing.
            self.mailer = SendEmail(self.date)
            self.logger.info("Send email module initialized.")
            self.uploader = UploadFastq(self.date)
            self.logger.info("Fastq upload module initialized.")
            
            
            self.updater = update_db.UpdateFastQC(self.date)
            self.logger.info("Update FastQC to DataBase.")
            
            
             
        except Exception as e:
            self.logger.info("Module initialization error. Check the below error log.")
            self.logger.info(e)
            sys.exit()

    def __call__(self):
        # init class
        self._initialize_module()
        # init db parser.
        fc_data_dict = self.parse_db()
        fc_name_list = []
        for run_id, info in fc_data_dict.items():
            if not info['sequencer_fc_dir']:
                self.logger.info(f"{run_id} has no FlowCell directory!")
                continue
            # if not len(info['bulk_sample_df']) and not len(info['sc_sample_df']):
            if not len(info['bulk_sample_df']):
                self.logger.info(f"{run_id} has no Bulk Samplesheet!")
                continue
            fc_name = self.run_single_fc(info)
            fc_name_list.append(fc_name)
        for fc_name in fc_name_list:
            self.uploader(fc_name)
        if len(fc_name_list):
            headline = f"{self.date} AWS bucket uploaded."
            fc_list_string = '\n'.join(fc_name_list)
            main_text = f"""{self.date} FASTQ uploaded.
            Uploaded flowcell directory:
            {fc_list_string}"""
            self.mailer(headline, main_text)
        else:
            pass
        pass

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--target_date", required=True,
        help='Sequencing batch date. Module will run target day.')
    parser.add_argument("--overwrite", default=False, action='store_true')
    args = parser.parse_args()
    mf = MakeFastq(
        target_date=args.target_date,
        overwrite=args.overwrite)
    exit(mf())


