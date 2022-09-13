import pandas as pd
from pathlib import Path

# define functions.
def read_metrics(metric_path: Path) -> pd.Series:
    sample_id = metric_path.name.split(".")[0]
    with open(metric_path) as f:
        lines = f.readlines()
    header_line = lines[6]
    metric_line = lines[7]
    return_series = pd.Series(data=metric_line.split("\t"), index=header_line.split("\t"))
    return_series['Sample ID'] = sample_id
    return return_series


def make_metric_df(data_dir: Path) -> pd.DataFrame:
    data_list = list(data_dir.glob("*recal.bam.metrics.txt"))
    metrics = [read_metrics(x) for x in data_list]
    metrics_df = pd.DataFrame(metrics)
    for col in metrics_df.columns:
        try:
            metrics_df[col] = metrics_df[col].astype(float)
        except:
            pass
    return metrics_df

#check MEAN, MEDIAN, --COVERAGE CAP 200비교
def coverage(metric_path: Path) -> pd.Series:
    with open(metric_path) as f:
        lines = f.readlines()
    coverage_line = lines[1].splot(" --")[6]
    coverage = float(coverage_line.split(" ")[1])
    
#위의 def 안썼음,,,어어,,
with open("CD_22_14544_BD_D_SRW_1.recal.bam.metrics.txt") as f:
    files=f.readlines()
coverage=float(files[1].split(" --")[6].split(" ")[1])


#define_function 
data_dir = Path('.')
metrics_df =make_metric_df(data_dir)
metrics_df["Coverage_Filter"] = ['Warning' if (row['MEAN_TARGET_COVERAGE'] < coverage) & (row['MEDIAN_TARGET_COVERAGE'] < coverage) else 'PASS'  for _, row in metrics_df.iterrows()]
metrics_df.to_csv("metric_summary.csv", header=True, index=False)
