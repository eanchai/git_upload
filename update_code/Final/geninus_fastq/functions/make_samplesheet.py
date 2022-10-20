#!env /usr/bin/python3

'''

Date: 2021.11.23
Authors: duaghk

Make samplesheet by index types.

'''

# import library
import logging
import pandas as pd
from pathlib import Path
from pytz import timezone
from datetime import datetime
from .configs import SequencingConfig

class MakeSampleSheet(SequencingConfig):
    '''
        SampleSheet making scripts.
        Input:
            output-dir with fc_dir, (e.g. /data6/NGS_data/fastq/211203_A01192_0278_BHHWFCDRXY)
            logger: for logging.
            bulk_sample_df:
                DB parsed samplesheet: pd.DataFrame format.
                Columns:
                    sample_id, i5_index, i7_index
            sc_sample_df:
                DB parsed samplesheet: pd.DataFrame format.
                Columns:
                    sample_id, Index(SI id like SI-TT-A1)
        Output:
            SampleSheet paths.
            Save format: csv
            Save file:
                format:SampleSheet_{index_type}_{index_length}.csv
            Save SampleSheet:
                SingleCell:
                    Lane('*' masked), Sample, Index
                BulkSeq:
                    with header.
                    Lane, SampleID, sample_name, sample_plate, sample_well, Index2, Index, Sampleproject, Description(optional)
    '''
    def __init__(self, date):
        '''
            Initlize.
            Input: 
                outdir: str, must be full path with fc_dir.
                logger: logging module. for logging.

        '''
        # for bulk header.
        SequencingConfig.__init__(self)
        self.header = (
            '[Data],,,,,,,,\n,,,,,,,,\n'
        )
        self.log_dir = self.log_dir.joinpath("MakeSampleSheet")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self.make_logger(f"MakeSampleSheet_{date}")

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
    # @staticmethod
    # def _make_sc_df_index(sc_sample_df:pd.DataFrame) -> pd.DataFrame:
    #     '''
    #         make single cell index.
    #         Input:
    #             sc_sample_df: pd.DataFrame
    #                 Columns: ["sample_id", "index_plate_cd", "index_no_cd", "service_cd"]
    #         Returns:
    #             pd.DataFrame.
    #                 Columns: ["sample_id", "Index"]
    #     '''
    #     sc_sample_df["Index"] = [f"{rows['index_plate_cd']}-{rows['index_no_cd']}" for _, rows in sc_sample_df.iterrows()]
    #     # column reordering.
    #     sc_sample_df = sc_sample_df[["sample_id", "Index", "service_cd"]]
    #     return sc_sample_df
    
    @staticmethod
    def _rename_columns(sample_df: pd.DataFrame, singlecell:bool = False) -> pd.DataFrame:
        '''
            rename columns for sample_df
            Input:
                sample_df: pd.DataFrame.
                    Columns:
                        singlecell: ["sample_id", "Index", "service_cd"]
                        bulk: ["sample_df", "i5_index", "i7_index", "service_cd"]
            Returns:
                pd.DataFrame.
                    Columns:
                        singlecell: ["Sample", "Index", "Sampleproject"]
                        bulk: ["SampleID", "Index2", "Index", "Sampleproject"]
        '''
        if singlecell:
            rename_columns = {"sample_id": "Sample", "service_cd": "Sampleproject"}
        else:
            rename_columns = {"sample_id": "SampleID", "i5_index": "Index2", 
                              "i7_index": "Index", "service_cd": "Sampleproject"}
        sample_df = sample_df.rename(columns=rename_columns)
        return sample_df
    
    @staticmethod
    def _add_additional_columns(sample_df: pd.DataFrame, singlecell:bool = False) -> pd.DataFrame:
        '''
            Add samplesheet default columns for requirement of bcl2fastq.
            SingleCell sample ddataframe only need Lane.
            bulk sample dataframe need many.
        '''
        if singlecell:
            sample_df.insert(0, "Lane", "*")
        else:
            sample_df.insert(0, "Lane", "")
            sample_df.insert(2, "sample_name", "")
            sample_df.insert(3, "sample_plate", "")
            sample_df.insert(4, "sample_well", "")
        return sample_df

    @staticmethod
    def _make_mask(masking_series, target_val) -> list:
        '''
            For service code masking.
            service code are 3 letters.
            first are S or J.
        '''
        masking = pd.Series([True if x[1:3] == target_val else False for x in masking_series])
        # masking = [True if x[1:3] == target_val else False for x in masking_series]
        return masking

    def _tag_index_type(self, samplesheet:pd.DataFrame, singlecell:bool = False) -> pd.DataFrame:
        '''
            Add index type for split data.
        '''
        # init value
        samplesheet["idx_type"] = ""
        if singlecell:
            samplesheet.loc[self._make_mask(samplesheet["Sampleproject"], "SV"), "idx_type"] = "visium"
            samplesheet.loc[~self._make_mask(samplesheet["Sampleproject"], "SV"), "idx_type"] = "chromium"
        else:
            samplesheet.loc[self._make_mask(samplesheet["Sampleproject"], "LN").tolist(), "idx_type"] = "LS"
            samplesheet.loc[
                (~self._make_mask(samplesheet["Sampleproject"], "LN")) & (samplesheet["Index2"].isna()),
                "idx_type"] = "single"
            samplesheet.loc[
                (~self._make_mask(samplesheet["Sampleproject"], "LN")) & (~samplesheet["Index2"].isna()),
                "idx_type"] = "dual"
        return samplesheet

    @staticmethod
    def _arrange_ls_index(samplesheet:pd.DataFrame) -> pd.DataFrame:
        '''
            LiquidSCAN need to arrange index.
            forward 4 sequence: Sample Index
            backward 4 sequence: UMI-like
        '''
        samplesheet["Index"] = [row["Index"][:4] 
                                if row["idx_type"] =="LS" 
                                else row["Index"] 
                                for _, row in samplesheet.iterrows()]
        return samplesheet

    @staticmethod
    def _tag_index_len(samplesheet:pd.DataFrame) -> pd.DataFrame:
        '''
            Tag samplesheet index length for groupby saving.
            Index column is in samplesheet.
        '''
        samplesheet["idx_len"] = [len(row["Index"]) 
                                  for _, row in samplesheet.iterrows()]
        return samplesheet
    
    @staticmethod
    def _drop_columns(samplesheet: pd.DataFrame, columns: list) -> pd.DataFrame:
        samplesheet = samplesheet.drop(columns=columns)
        return samplesheet
    
    def _add_bulk_header(self, filepath: str):
        with open(filepath, 'a') as f:
            f.write(self.header)

    def parse_sc_sample_df(self, sc_sample_df: pd.DataFrame) -> pd.DataFrame:
        '''
            Parsing singlce cell sample dataframe.
            Input:
                sc_sample_df: pd.DataFrame. single cell index information.
                    Columsn: ["sample_id", "Index", "service_cd"]
            Returns:
                returndict: named SampleSheet dictionary.
                    chromium_ss: pd.DataFrame. chromium dataframe for running cellranger
                        Columns: ["SampleID", "Index"]
                    visium_ss: pd.DataFrame. visium dataframe for running spaceranger
                        Columns: ["SampleID", "Index"]
        '''
        returndict = {}
        # sc_sample_df = self._make_sc_df_index(sc_sample_df)
        sc_sample_df = self._rename_columns(sc_sample_df, singlecell=True)
        sc_sample_df = self._add_additional_columns(sc_sample_df, singlecell=True)
        sc_sample_df = self._tag_index_type(sc_sample_df, singlecell=True)
        # split.
        chromium_ss = sc_sample_df[sc_sample_df["idx_type"].isin(["chromium"])]
        chromium_ss = self._drop_columns(chromium_ss, ["idx_type", "Sampleproject"])
        visium_ss = sc_sample_df[sc_sample_df["idx_type"].isin(["visium"])]
        visium_ss = self._drop_columns(visium_ss, ["idx_type", "Sampleproject"])
        if len(chromium_ss):
            returndict["chromium_-"] = chromium_ss
        if len(visium_ss):
            returndict["visium_-"] = visium_ss
        return returndict
    
    def parse_bulk_sample_df(self, bulk_sample_df: pd.DataFrame) -> pd.DataFrame:
        '''
            Parsing bulk sample dataframe.
            Input:
                bulk_sample_df: pd.DataFrame. bulk sequencing index information.
                    Columsn: ["sample_id", "i5_index, "i7_index", "service_cd"]
            Returns:
                returndict: named SampleSheet dictionary.
                    All dataframes' Columns are samme.
                    Columns: ["Lane", "SampleID", "sample_name", "sample_plate", "sample_well", "Index", "Index2", "Sampleproject"]
                    ls_ss: pd.DataFrame. LS indexed dataframe for running LiquidSCAN bcl2fastq.
                    single_8: pd.DataFrame. Index only dataframe. indexs' length are 8-base.
                    dual_8: pd.DataFrame. Index2, Index dataframe. indexs' length are 8-base.
                    dual_10: pd.DataFrame. Index2, Index dataframe. indexs' length are 10-base.
        '''
        returndict = {}
        bulk_sample_df = bulk_sample_df[~bulk_sample_df['sample_id'].isin(['sample_id'])]
        bulk_sample_df = self._rename_columns(bulk_sample_df)
        bulk_sample_df = self._add_additional_columns(bulk_sample_df)
        # print(bulk_sample_df)
        bulk_sample_df = self._tag_index_type(bulk_sample_df)
        bulk_sample_df = self._tag_index_len(bulk_sample_df)
        bulk_sample_df = self._arrange_ls_index(bulk_sample_df)
        
        # split.
        ls_ss = bulk_sample_df[bulk_sample_df["idx_type"].isin(["LS"])]
        single_8 = bulk_sample_df[bulk_sample_df["idx_type"].isin(["single"])]
        dual_8 = bulk_sample_df[(bulk_sample_df["idx_type"].isin(["dual"])) & (bulk_sample_df['idx_len'].isin([8]))]
        dual_10 = bulk_sample_df[(bulk_sample_df["idx_type"].isin(["dual"])) & (bulk_sample_df['idx_len'].isin([10]))]
        # drop.
        ls_ss = self._drop_columns(ls_ss, ["idx_type", "idx_len"])
        single_8 = self._drop_columns(single_8, ["idx_type", "idx_len"])
        dual_8 = self._drop_columns(dual_8, ["idx_type", "idx_len"])
        dual_10 = self._drop_columns(dual_10, ["idx_type", "idx_len"])
        if len(ls_ss): 
            returndict["ls_4"] = ls_ss
        if len(single_8): 
            returndict["single_8"] = single_8
        if len(dual_8): 
            returndict["dual_8"] = dual_8
        if len(dual_10):
            returndict["dual_10"] = dual_10
        return returndict

    def __call__(self, samplesheet_dir:Path, bulk_sample_df: pd.DataFrame = None, sc_sample_df: pd.DataFrame = None) -> list:
        sample_fn_list = []
        # make outdir. 
        if type(bulk_sample_df) == pd.DataFrame and len(bulk_sample_df):
            self.logger.info(f"Bulk sequencing sample data inserted. Sample numbers: {len(bulk_sample_df)}")
            bulk_df_dict = self.parse_bulk_sample_df(bulk_sample_df)
        else:
            bulk_df_dict = None
        if type(sc_sample_df) == pd.DataFrame and len(sc_sample_df):
            single_df_dict = self.parse_sc_sample_df(sc_sample_df)
            self.logger.info(f"SingleCell sequencing sample data inserted. Sample numbers: {len(sc_sample_df)}")
        else:
            single_df_dict = None
        # save data with index length
        if bulk_df_dict:
             for k, v in bulk_df_dict.items():
                self.logger.info(f"{len(v)} {k} SampleSheet saving...")
                filename = samplesheet_dir.joinpath(f"SampleSheet_{k}.csv")
                if filename.exists():
                    filename.unlink()
                self._add_bulk_header(filename)
                v.drop_duplicates().to_csv(filename, mode='a', header=True, index=False)
                sample_fn_list.append(filename)
        if single_df_dict:
            for k, v in single_df_dict.items():
                self.logger.info(f"{len(v)} {k} SampleSheet saving...")
                filename = samplesheet_dir.joinpath(f"SampleSheet_{k}.csv")
                if filename.exists():
                    filename.unlink()
                v.drop_duplicates().to_csv(filename, header=True, index=False)
                sample_fn_list.append(filename)
        return sample_fn_list


if __name__ == "__main__":
    import logging
    from pytz import timezone
    from datetime import datetime

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--date", help='date.')
    parser.add_argument("--samplesheet_dir", help="fastq_fc_output_dir's samplesheet dir.")
    parser.add_argument("--bulk_sample_info_path", default=None)
    parser.add_argument("--sc_sample_info_path", default=None)
    args = parser.parse_args()

    if args.bulk_sample_info_path:
        if args.bulk_sample_info_path.endswith("xlsx"):
            bulk_ss = pd.read_excel(args.bulk_sample_info_path)
        else:
            bulk_ss = pd.read_csv(args.bulk_sample_info_path)
    else:
        bulk_ss = None
    if args.sc_sample_info_path:
        if args.sc_sample_info_path.endswith("xlsx"):
            sc_ss = pd.read_excel(args.sc_sample_info_path)
        else:
            sc_ss = pd.read_csv(args.sc_sample_info_path)
    else:
        sc_ss = None
    make_ss = MakeSampleSheet(args.date)
    sample_fn_list = make_ss(Path(args.samplesheet_dir), bulk_ss, sc_ss)
    print(sample_fn_list)
    
    
    