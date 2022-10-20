#!env /usr/bin/python3

'''

Date: 2021.11.11
Authors: duaghk

Test scripts.

'''

# import library
import os
import pytest
import logging
# import testing modules.
from functions import parse_db, make_samplesheet, bcl_to_fastq, ranger


class TestMakeSS:
    def __init__(self, input_ss, outdir):
        self.input_ss = input_ss
        self.outdir = outdir
        pass
    def make_logger(self, log_dir, name=None, consoleset=True, streamset=True):
        '''
            logger making.
        '''
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)
        loggerformat = "%(asctime)s - %(name)s - %(module)s - %(levelname)s - %(message)s"
        formatter = logging.Formatter(loggerformat)
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

    def test_init(self):
        pass

