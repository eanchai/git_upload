#!env /usr/bin/python3

'''

Date: 2021.12.12
Authors: duaghk

change fastq.gz name and run fastqc.

'''

import re
import logging
import subprocess as sp
from pathlib import Path
from pytz import timezone
from datetime import datetime
from .configs import ToolConfig, SequencingConfig

class FastQCRunner(ToolConfig, SequencingConfig):
    def __init__(self, date):
        '''
            Input:
                fastqc: fastqc tool path.
                threads: for running fastqc. fastqc use 2 threads default in one fastq.gz.
                logger: logging threads.
        '''
        ToolConfig.__init__(self)
        SequencingConfig.__init__(self)
        self.log_dir = self.log_dir.joinpath("RunFastQC")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self.make_logger(f"RunFastQC_{date}")
        pass

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
    

    def find_fastq(self, fastq_dir: str) -> list:
        '''
            Find fastq in input fastqdir.
        '''
        returnlist = [x for x in fastq_dir.glob("*.fastq.gz") if "Undetermined" not in str(x)]
        self.logger.info(f"{len(returnlist)} fastq.gz Found!")
        return returnlist

    def _change_fastqname(self, fastq_path: str, idx_type: str):
        '''
            Change raw bcl2fastq name.
            Input: fastq full path.
            Functions:
                change fastq name raw to ${sample_id}_R[1-2].fastq.gz
            Return:
                None.
        '''
        dirname = fastq_path.parent
        basename = fastq_path.name
        sample_id = re.split(r"_S[0-9]+_R[1-3]_[0-9]{3}.fastq.gz", basename)[0]

        if sample_id == basename:
            new_path = fastq_path
        else:
            if idx_type == 'ls_4':
                read_prefix = re.split(sample_id, basename)[1]
                read_prefix = read_prefix.split("_")[2]  # Cuz, read index always like _S001_R1_001.fastq.gz
                if read_prefix == "R2":
                    read_prefix = "idx"
                elif read_prefix == 'R3':
                    read_prefix = 'R2'
                else:
                    pass
                new_path = dirname.joinpath(f"{sample_id}.{read_prefix}.fastq.gz")
                mv_cmd = f"mv {fastq_path} {new_path}"
                self.logger.info(f"Fastq rename commands: {mv_cmd}")
                sp.run(mv_cmd, shell=True)
            else:
                read_prefix = re.split(sample_id, basename)[1]
                read_prefix = read_prefix.split("_")[2]  # Cuz, read index always like _S001_R1_001.fastq.gz
                new_path = dirname.joinpath(f"{sample_id}.{read_prefix}.fastq.gz")
                mv_cmd = f"mv {fastq_path} {new_path}"
                self.logger.info(f"Fastq rename commands: {mv_cmd}")
                sp.run(mv_cmd, shell=True)
        return new_path
    
    def change_fastq_multi(self, fastq_list: list, idx_type: str):
        new_fastq_list = [self._change_fastqname(x, idx_type) for x in fastq_list]
        return new_fastq_list

    def make_fastqc_commands(self, fastq_path: Path, qc_outdir: Path):
        fastq_name = fastq_path.name.rsplit(".", 2)[0]
        donefile = qc_outdir.joinpath(f"{fastq_name}.done")
        if donefile.exists():
            cmd = 'echo ""'
        else:
            cmd = (
                f"{self.fastqc} -o {qc_outdir} -t {self.fastqc_threads} --extract {fastq_path}"
            )
        # add done file generator.
            cmd += f" && touch {qc_outdir.joinpath(fastq_path.name.rsplit('.', 2)[0] + '.done')}"
        return cmd
    
    def __call__(self, fastq_dir: Path, outdir: Path, idx_type: str):
        fastq_list = self.find_fastq(fastq_dir)
        new_fastq_list = self.change_fastq_multi(fastq_list, idx_type)
        new_fastq_list = [x for x in new_fastq_list if str(x).split(".")[1] != 'idx']
        # need to make prefix.
        # fastqc_commands = [self.make_fastqc_commands(x, qc_outdir) for x in new_fastq_list]
        fastqc_commands = {
            x.name.rsplit(".", 2)[0]: {
                "cmd":self.make_fastqc_commands(x, outdir),
                "outdir":outdir
                } for x in new_fastq_list
            }
        return fastqc_commands

if __name__ == "__main__":
    import argparse
    # set parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--date")
    parser.add_argument("--fastq_dir", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    qc_runner = FastQCRunner(args.date)
    qc_command_info_dict = qc_runner(Path(args.fastq_dir), Path(args.outdir), Path(args.fastq_dir).name)
    for k,v in qc_command_info_dict.items():
        print(f"{k}'s FASTQC commands: {v['cmd']}, FASTQC output directory: {v['outdir']}")
        with open(v['outdir'].joinpath(f"{k}_fastqc.sh"), 'w') as f:
            f.write(v['cmd'])
