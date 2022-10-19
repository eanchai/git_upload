#!env /usr/bin/python3

'''

Date: 2021.12.12
Authors: duaghk

Parse fastqc module.

'''

import logging
import pandas as pd
from pathlib import Path
from pytz import timezone
from datetime import datetime
from .configs import SequencingConfig


class ParseFastQC(SequencingConfig):
    def __init__(self, date):
        SequencingConfig.__init__(self)
        self.log_dir = self.log_dir.joinpath("ParseQC")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self.make_logger(f"ParseQC_{date}")

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
    def _line_splitter(self, fastqc_output: str) -> list:
        '''
            fastqc_data.txt file line splitter.
            fastqc_data.txt are separated by '>>END_MODULE\n'.
            So, blocked line split need.
        '''
        returnlist = []
        with open(fastqc_output) as f:
            lines = f.readlines()
        while True:
            try:
                endindex = lines.index(">>END_MODULE\n")
                returnlist.append(lines[:endindex])
                lines = lines[endindex+1:]
            except ValueError:
                break
        return returnlist

    def _make_table(self, lines: list) -> pd.DataFrame:
        '''
            Fastqc output line to dataframe.
        '''
        returndf = []
        for line in lines:
            # print(line)
            if line.startswith(">>"):
                continue
            if line.startswith("#"):
                header = line.strip().split("\t")
                header[0] = header[0].strip("#")
            else:
                tmp_values = line.strip().split("\t")
                tmpseries = pd.Series(data = tmp_values, index = header)
                returndf.append(tmpseries)

        returndf = pd.DataFrame(returndf)
        for col in returndf.columns:
            try:
                returndf[col] = returndf[col].astype("float")
            except:
                pass
        return returndf
    
    def _get_gc_percent(self, basic_lines: list) -> float:
        gc_line = [x for x in basic_lines if x.startswith("%\GC")][0]
        gc_percent = gc_line.strip().split("\t")[1]
        return gc_percent

    def _get_base_numbers(self, seq_len_table: pd.DataFrame) -> float:
        try:
            seq_len_table["base_count"] = seq_len_table["Length"] * seq_len_table["Count"]
        except:
            seq_len_table["base_count"] = 0
        return seq_len_table["base_count"].sum()
        
    def _parse_fastqc(self, fastqc_output: Path) -> dict:
        '''
            Parse fastqc output.
            lines explain.
                lines[0]: Basic statistics.
                lines[1]: Per base sequence quality.
                lines[2]: Per tile sequence quality. (ignore)
                lines[3]: Per sequence quality score.
                lines[4]: Per base sequence content.
                lines[5]: Per sequence GC content. (ignore)
                lines[6]: Per base N content. (ignore)
                lines[7]: Sequence Length Distribution.
                lines[8]: Sequence Duplication Levels. (ignore)
                lines[9]: Sequence Overrepresented sequence (ignore)
                lines[10]: Adapter Content. (ignore)
        '''
        # init value
        returndict = {}
        lines = self._line_splitter(fastqc_output)
        abnormality = 0 if len(lines) == 11 else -1
        bsq = self._make_table(lines[1])
        sqc = self._make_table(lines[3+abnormality])
        bsc = self._make_table(lines[4+abnormality])
        sgc = self._make_table(lines[5+abnormality])
        seq_len_table = self._make_table(lines[7+abnormality])
        product_base_size = self._get_base_numbers(seq_len_table)
        product_read_size = seq_len_table["Count"].sum()
        q30 = sqc["Count"][sqc["Quality"] >= 30].sum()/sqc["Count"].sum() * 100
        returndict["bsq"] = bsq
        returndict["bsc"] = bsc
        returndict['sge'] = sgc
        returndict["Total_bases"] = product_base_size
        returndict["Total_reads"] = product_read_size
        returndict["Q30"] = round(q30, 3)
        return returndict
    
    @staticmethod
    def _get_dir_list(fastqc_dir: Path) -> list:
        sample_fastqc_dirs = [x for x in fastqc_dir.glob("*") if x.is_dir()]
        return sample_fastqc_dirs

    def parse_multi_fastqc(self, fastqc_dir: Path) -> dict:
        returndict = {}
        fastqc_dir_list = self._get_dir_list(fastqc_dir)
        for sample_qc_dir in fastqc_dir_list:
            prefix = sample_qc_dir.name
            sample_id, read_id = prefix.split(".")
            read_id = read_id.split("_")[0]
            if sample_id not in returndict.keys():
                returndict[sample_id] = {}
            returndict[sample_id][read_id] = self._parse_fastqc(sample_qc_dir.joinpath("fastqc_data.txt"))
        return returndict

    @staticmethod
    def _make_sample_fastqc_series(sample_id: str, dict_items: dict) -> pd.Series:
        data_list = [sample_id]
        index_list = [
            "SampleID",
            "FASTQ_TOTAL_BASE_R1", "FASTQ_TOTAL_READ_R1", "FASTQ_Q30_R1", 
            "FASTQ_TOTAL_BASE_R2", "FASTQ_TOTAL_READ_R2", "FASTQ_Q30_R2"
            ]
        # check whether R1 info exists
        if "R1" in dict_items.keys():
            r1_list = [dict_items["R1"]["Total_bases"], dict_items["R1"]["Total_reads"], dict_items["R1"]["Q30"]]
        else:
            r1_list = ["", "", ""]
        # check whether R2 info exists
        if "R2" in dict_items.keys():
            r2_list = [dict_items["R2"]["Total_bases"], dict_items["R2"]["Total_reads"], dict_items["R2"]["Q30"]]
        else:
            r2_list = ["", "", ""]
        data_list = data_list + r1_list + r2_list
        returnseries = pd.Series(data=data_list, index=index_list)
        return returnseries
    
    @staticmethod
    def _calculate_total_base(r1_base: int or str, r2_base: int or str) -> int:
        int_list = [int(x) for x in [r1_base, r2_base] if x != ""]
        return sum(int_list)
    
    def make_product_result_table(self, result_dict: dict):
        result_table = [self._make_sample_fastqc_series(sample_id, v) for sample_id, v in result_dict.items()]
        result_table = pd.DataFrame(result_table)
        result_table["FASTQ_TOTAL_BASES(Gb)"] = [
            self._calculate_total_base(rows["FASTQ_TOTAL_BASE_R1"], rows["FASTQ_TOTAL_BASE_R2"])
            for _, rows in result_table.iterrows()]
        giga = 10**9
        result_table["FASTQ_TOTAL_BASES(Gb)"] = [round(x/giga, 2) for x in result_table["FASTQ_TOTAL_BASES(Gb)"]]
        result_table = result_table.sort_values(by='SampleID')
        # Integer change.
        for col in ["FASTQ_TOTAL_BASE_R1", "FASTQ_TOTAL_READ_R1", "FASTQ_TOTAL_BASE_R2", "FASTQ_TOTAL_READ_R2"]:
            try:
                result_table[col] = result_table[col].astype(int)
            except:
                pass
        return result_table

    def __call__(self, fastqc_dir: Path):
        fastqc_dict = self.parse_multi_fastqc(fastqc_dir)
        self.logger.info(f"{len(fastqc_dict.keys())} numbers FASTQC result parsed.")
        product_table = self.make_product_result_table(fastqc_dict)
        return fastqc_dict, product_table

if __name__ == "__main__":
    import logging
    from pytz import timezone
    from datetime import datetime

    import argparse
    import pickle
    parser = argparse.ArgumentParser()
    parser.add_argument("--date")
    parser.add_argument("--fastqc_dir")
    parser.add_argument("--outdir")
    args = parser.parse_args()

    qcparser = ParseFastQC(args.date)
    parsed_dict, parsed_table = qcparser(Path(args.fastqc_dir))
    with open(Path(args.outdir, "fastqc_result.pickle"), 'wb') as f:
        pickle.dump(parsed_dict, f)
    parsed_table.to_csv(Path(args.outdir, "fastqc_product_result.csv"), header=True, index=False)

