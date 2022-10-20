#!env /usr/bin/python3

'''

Date: 2022.07.25
Authors: duaghk

Upload AWS scripts.

'''

import sys
import logging
import subprocess as sp
from pytz import timezone
from pathlib import Path
from datetime import datetime
from .configs import SequencingConfig, UploadConfig

class UploadFastq(SequencingConfig, UploadConfig):
    def __init__(self, date) -> None:
        SequencingConfig.__init__(self)
        UploadConfig.__init__(self)
        self.log_dir = self.log_dir.joinpath("UploadFastq")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self.make_logger(f"UploadFastq_{date}")

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

    def run_subprocess(self, cmd):
        self.logger.info(f"Command: {cmd}")
        proc = sp.run(cmd, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, universal_newlines=True)
        if proc.returncode:
            self.logger.info(proc.stderr)
            sys.exit()
        self.logger.info(f"Stderr: {proc.stderr}")
        self.logger.info(f"Stdout: {proc.stdout}")

    def __call__(self, fc_dir: str) -> None:
        upload_dir_path = self.fastq_dir.joinpath(fc_dir)
        upload_bucket_path = f"s3://{self.upload_bucket}/{fc_dir}"
        cmd = f"aws s3 sync --no-progress {upload_dir_path} {upload_bucket_path}"
        self.run_subprocess(cmd)








