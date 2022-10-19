import os
import pandas as pd

fn = "/Users/duaghk/build/iGeniPipe/testfile/qc_test/CD_21_15384_XT_CF.R2_fastqc/fastqc_data.txt"


# valuelist[0] # default info
# valuelist[1] # Per base sequence quality\tpass, Base mena median lower quartile upper quartile 10th percentile 90th percentile
# valuelist[2] # Per tile sequence quality
# valuelist[3] # Per sequence quality score, Quality count
# valuelist[4] # Per base sequence content
# valuelist[5] # Per sequenece GC content, GC content Count
# valuelist[6] # Per base N content, Base N-Count
# valuelist[7] # Sequence Length Distribution, Length Count, (151, read count)
# valuelist[8] # Sequence Duplication Levels
# valuelist[9] # Overrepresented sequence
# valuelist[10] # Adapter Content


def _line_splitter(fastqc_output):
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

def _make_table(lines):
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

lines = _line_splitter(fn)
bsq = _make_table(lines[1])
sqc = _make_table(lines[3])
bsc = _make_table(lines[4])
sgc = _make_table(lines[5])
seq_len_table = _make_table(lines[7])

seq_len_table
seq_len_table["base_count"] = seq_len_table["Length"] * seq_len_table["Count"]
seq_len_table["base_count"].sum()

# test
lines.index(">>END_MODULE\n")
lines.count(">>END_MODULE\n")

lines[0]

gc_line = [x for x in lines[0] if x.startswith("%GC")]
gc_line

sqc
sqc["Count"][sqc["Quality"] >= 30].sum()/sqc["Count"].sum() * 100

tgt_dir = "/Users/duaghk/build"
targetdictlist = [os.path.join(tgt_dir, x) for x in os.listdir(tgt_dir) if os.path.isdir(os.path.join(tgt_dir, x))]

dirs = [entry.path for entry in os.scandir(tgt_dir) if entry.is_dir()]
dirs

next(os.walk(tgt_dir))[1]
