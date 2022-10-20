#!env /usr/bin/python3

'''

Date: 2021.11.11
Authors: duaghk

Run CellRanger mkfastq

'''

import os

class RangerRunner:
    def __init__(self, cellranger: str , spaceranger: str, threads: int):
        self.cellranger = cellranger
        self.spaceranger = spaceranger
        self.threads = threads

    def run_ranger(self, rundir: str, outdir: str, samplesheet_path: str) -> str:
        # set outdir 
        samplesheet_info = os.path.basename(samplesheet_path)
        samplesheet_info = os.path.splitext(samplesheet_info)[0]
        _, idx_type, _ = samplesheet_info.split("_")
        # Make index type-wise outdir
        fc_outdir = os.path.join(outdir, idx_type)
        os.makedirs(fc_outdir, exist_ok=True)

        if idx_type == "visium":
           cmd = f"{self.spaceranger} mkfastq "
        else:
            cmd = f"{self.cellranger} mkfastq "
        cmd += (
            f"--run={rundir} "
            f"--samplesheet={samplesheet_path} "
            f"--output-dir={fc_outdir} "
            # f"--localcores={self.threads} "
            # f"--localmems={self.threads*4} "
        )
        return cmd, fc_outdir
    
    def __call__(self, rundir: str, outdir: str, samplesheet_path: str, logger=None) -> str:
        ranger_command, fc_outdir = self.run_ranger(rundir, outdir, samplesheet_path)
        logger.info(f"SingleCell cmd: {ranger_command}")
        return ranger_command, fc_outdir 
