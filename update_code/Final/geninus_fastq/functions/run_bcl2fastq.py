#!env /usr/bin/python3

'''

Date: 2021.11.11
Authors: duaghk

Run BCL2FASTQ script.

'''

# import library
import logging
from datetime import datetime
from pytz import timezone
from pathlib import Path

# import config.
from .configs import ToolConfig, SequencingConfig

class Bcl2fastqRunner(ToolConfig, SequencingConfig):
    '''
        Run BCL2FASTQ modules.
    '''
    def __init__(self, date):
        ToolConfig.__init__(self)
        SequencingConfig.__init__(self)
        # check if Sequencing was dual or single only sample.
        # thread settings. l:p:w = 1:2:1
        self.l_threads = int((self.bcl_threads//2))
        self.p_threads = int(self.bcl_threads)
        self.w_threads = self.l_threads
        self.log_dir = self.log_dir.joinpath("RunBcl2Fastq")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self.make_logger(f"RunBcl2Fastq_{date}")
        self.logger.info(f'Sequencing Inforamtion parsing start. Target date: {date}')
    # define functions
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
    def parse_bases_mask(bases_mask: str, idx_type: str) -> str:
        mask_splitted = bases_mask.split(',')
        if idx_type == "ls_4":
            mask_splitted[1] = 'I4Y4'
            if len(mask_splitted) == 4:
                mask_splitted[2] = 'N*'
        return_str = ",".join(mask_splitted)
        return return_str

    def run_bcl2fastq(self, sequencer_fc_dir: str, bases_mask:str, samplesheet_path: str,
                      barcode_mismatch:bool = False, mask_short:bool = False, 
                      minimum_trimmed:int = False, ignore_missing_bcls:bool = False) -> str:
        '''
            BCL2FASTQ command maker.
            Input:
                bcldir: path, BCL directory path.
                bases_mask: str, base_mask for masking options. e.g. Y*,I4Y4N*,N*,Y*
                samplesheet_path: path, bcl2fastq input samplesheet path.
                outdir: path, Fastq writing dir. Path must have fc_dir.
                Other options: bcl2fastq options. if None, not added.
            Returns:
                BCL2FASTQ running commands with done file make commands.
        '''
        # get bcl dir.
        fc_name = sequencer_fc_dir.name
        idx_type = samplesheet_path.name.split(".")[0].split("_", 1)[-1]
        bases_mask = self.parse_bases_mask(bases_mask, idx_type)
        fastq_fc_output_dir = self.fastq_dir.joinpath(fc_name, 'fastq', idx_type)
        fastq_fc_output_dir.mkdir(parents=True, exist_ok=True)
        cmd = (
            f"{self.bcl2fastq} "
            f"--runfolder-dir {sequencer_fc_dir} "
            f"--use-bases-mask {bases_mask} "
            f"--loading-threads {self.l_threads} "
            f"--processing-threads {self.p_threads} "
            f"--writing-threads {self.w_threads} "
            f"--sample-sheet {samplesheet_path} "
            "--no-lane-splitting "
            f"--output-dir {fastq_fc_output_dir} "
        )
        if barcode_mismatch: 
            cmd += f"--barcode-mismatches 0 "
        if mask_short:
            cmd += f"--mask-short-adapter-reads 0 "
        if minimum_trimmed:
            cmd += f"--minimum-trimmed-read-length 0 "
        if ignore_missing_bcls:
            cmd += "--ignore-missing-bcls"
        # add done file generator.
        cmd += f" && touch {fastq_fc_output_dir.joinpath(f'{idx_type}.done')}"
        
        return cmd, fastq_fc_output_dir

    def __call__(self, sequencer_fc_dir, bases_mask, samplesheet_path: str, mixed_seq: bool = False) -> list:
        '''
            One samplesheet, one commands, one prefix.
        '''
        self.logger.info("BCL2FASTQ running start.")
        # get RunInfo information. 
        sampleinfo = Path(samplesheet_path).name  # split samplesheet path for sample information.
        print(sampleinfo)
        _, idx_type, idx_len = sampleinfo.split(".")[0].split("_")
        self.logger.info(f"Index type: {idx_type} / Index length: {idx_len}")
        if idx_type == "ls":
            bcl_commands, fastq_fc_output_dir = self.run_bcl2fastq(sequencer_fc_dir, bases_mask, samplesheet_path, 
                                              barcode_mismatch=True, mask_short=True, minimum_trimmed=True)
                                            #   ignore_missing_bcls=True)
        else:
            if mixed_seq and idx_type == "single":
                bcl_commands, fastq_fc_output_dir = self.run_bcl2fastq(sequencer_fc_dir, bases_mask, samplesheet_path, barcode_mismatch=True)
            else:
                bcl_commands, fastq_fc_output_dir = self.run_bcl2fastq(sequencer_fc_dir, bases_mask, samplesheet_path,)
            self.logger.info(f"Commands: {bcl_commands}")

        self.logger.info(f"BCL2FASTQ running command generated.")
        return bcl_commands, fastq_fc_output_dir


if __name__ == "__main__":
    # import more library
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--date", required=True)
    parser.add_argument("--sequencer_fc_dir", type=str, required=True)
    parser.add_argument("--bases_mask", type=str, required=True)
    parser.add_argument("--mixed_seq", default = "False", action='store_true')
    parser.add_argument("--samplesheet_path", required=True)

    args = parser.parse_args()
    run_cls = Bcl2fastqRunner(args.date)
    run_cmd, fastq_fc_output_dir = run_cls(Path(args.sequencer_fc_dir), args.bases_mask, Path(args.samplesheet_path), args.mixed_seq)
    # os.system(f"qsh -cwd -q all.q -pe make {args.threads} {run_cmd}")
    print(f"Fastq Outdir: {fastq_fc_output_dir}")
    print(f"Command: {run_cmd}")
