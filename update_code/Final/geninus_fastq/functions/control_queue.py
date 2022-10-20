#!env /usr/bin/python3

'''

Date: 2021.12.12
Authors: duaghk

Queue control module.

'''
import logging
from pytz import timezone
from datetime import datetime
import time
import subprocess as sp
from pathlib import Path
from .configs import SequencingConfig

class QueueController(SequencingConfig):
    '''
        Queue controller for using subprocess.
    '''
    def __init__(self, date) -> None:
        SequencingConfig.__init__(self)
        self.log_dir = self.log_dir.joinpath("QueueControl")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self.make_logger(f"QueueControl_{date}")

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

    @staticmethod
    def make_queue_script(prefix: str, cmd: str, script_dir: Path) -> Path:
        '''
            Make script for queue jobs.
            Input:
                idx_type: str. [dual, single, visium, chromium]
                idx_len: int, [4, 8, 10, -]. -: SingleCell.
                cmd: Process run commands with done file generation.
                script_dir: script save dir. need to end with script purpose prefix.
            Returns:
                script_path: str. full script path.
                    e.g. /data/NGS_data/fastq/211206_A01192_0279_AHHWCYDRXY/scripts/bcl2fastq/dual_10.sh
                    e.g. /data/NGS_data/fastq/211206_A01192_0279_AHHWCYDRXY/scripts/fastqc/dual_10.sh
        '''
        prefix_dir = script_dir.joinpath(prefix)
        prefix_dir.mkdir(parents=True, exist_ok=True, mode=0o777)
        script_path = prefix_dir.joinpath(f"{prefix}.sh")
        with open(script_path, 'w') as f:
            f.write(cmd)
        return script_path


    @staticmethod
    def _submit_queue(threads: int, queue_script: str) -> str:
        '''
            functions for submit job script
            Input:
                threads: int, threads for job.
                queue_script: str, queue_script path.
            Returns:
                job_id: str, submitted job id.
        '''
        queue_cmd = f"qsub -wd {queue_script.parent} -q op.q -q proddev.q -pe make {threads} {queue_script}"
        # get job id.
        job_id = sp.check_output(queue_cmd, shell=True, universal_newlines=True)
        job_id = job_id.split()[2]
        return job_id
    
    def submit_multi_queue(self, threads: int, queue_script_list: list) -> list:
        job_id_list = [self._submit_queue(threads, x) for x in queue_script_list]
        return job_id_list
    
    @staticmethod
    def check_queue(job_id_list: list) -> None:
        while True:
            qstat_out = sp.check_output("qstat", shell=True, universal_newlines=True)
            # 3 break point.
            if qstat_out == '':
                break
            qstat_out = qstat_out.split("\n")[2:-1]
            qstat_list = [x.split()[0] for x in qstat_out]
            running_job = [x for x in job_id_list if x in qstat_list]
            if len(running_job):
                time.sleep(60)
            else:
                break
    
    @staticmethod
    def check_done_file(outdir_list: list) -> list:
        done_file_list = [x for outdir in outdir_list for x in outdir.glob("*.done")]
        return done_file_list

    def __call__(self, cmd_outdir_dict:dict, outdir: str, threads: int) -> bool:
        '''
            Input:
                cmddict:
                    Key: prefix. like: dual_8, dual_10, fastqc, etc.
                    Values:
                        cmd: str, commands.
                        outdir: str, outputdir.
                outdir: script output directory.
        '''
        # script dir need to set in __main__ script.
        queue_script_list = [self.make_queue_script(prefix, cmd_dir["cmd"], outdir) for prefix, cmd_dir in cmd_outdir_dict.items()]
        self.logger.info(f"{len(queue_script_list)} job script generated.")
        joblist = self.submit_multi_queue(threads, queue_script_list)
        self.logger.info(f"{len(joblist)} Queued!.")
        self.logger.info(f"Job id list: {', '.join(joblist)}")
        self.logger.info("Wait until queued job completed.")
        self.check_queue(joblist)
        self.logger.info("Job all finished.")
        outdir_list = [cmd_dir["outdir"] for cmd_dir in cmd_outdir_dict.values()]
        outdir_list = list(set(outdir_list))
        done_file_list = self.check_done_file(outdir_list)
        if len(done_file_list) == len(queue_script_list):
            self.logger.info("All job normally terminated.")
        else:
            self.logger.info("Queue process not normal termination!")
        # return done_file_list


if __name__ == "__main__":
    import pickle
    import argparse
    # set parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--cmd_info_pickle", default=None)
    parser.add_argument("--date", type=int)
    parser.add_argument("--threads", default=2, type=int)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    with open(args.cmd_info_pickle, 'rb') as f:
        queue_cmd_info_dict = pickle.load(f)
    q_ctrler = QueueController(args.date)
    done_list = q_ctrler(queue_cmd_info_dict, args.outdir, args.threads)
    for donefile in done_list:
        print(donefile)




        
