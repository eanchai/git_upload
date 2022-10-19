
# test Stat json file.

import os
import pandas as pd
 
print(os.getcwd())

statpath = "./testfile/Stats"
print(os.listdir(statpath))

tgtfile = os.path.join(statpath, "DemultiplexingStats.xml")

def get_read_info(barcoderoot, sample_id):
    barcode_seq = barcoderoot.attrib['name']
    barcodecount = 0
    perfectbarcodecount = 0
    for lane in barcoderoot:
        for count in lane:
            if count.tag == "BarcodeCount":
                barcodecount += int(count.text)
            elif count.tag == "PerfectBarcodeCount":
                perfectbarcodecount += int(count.text)
    returnseries = pd.Series(
        data=[sample_id, barcode_seq, barcodecount, perfectbarcodecount],
        index=["SampleID", "Barcode", "BarcodeCount", "PerfectBarcodeCount"]
        )
    return returnseries



def get_demultiplexing_info(runinfo_path):
    import xml.etree.ElementTree as ET  
    tree = ET.parse(runinfo_path)
    root = tree.getroot()
    readlist = list(root.iter("Sample"))
    returndf = [get_read_info(y, x.attrib['name']) for x in readlist for y in x if y.attrib['name']!='all']
    returndf = pd.DataFrame(returndf)
    return returndf

runstat = get_demultiplexing_info(tgtfile)
runstat


