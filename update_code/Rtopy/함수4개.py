"""
    date: 2022.09.26
    object: function, make dataframe for report
    author: 
"""
import os
import logging
import numpy as np
import pandas as pd
import subprocess as sp
import sqlalchemy as db
from time import time
from pathlib import Path
from pytz import timezone
from datetime import datetime

# Generate patient's barcode by insering either 1("het") or 2("hom") to correct position(rs number)
def make_patient_barcode(Gene_Bcode, ptn, Gene_refATGC, logdir, ptnID):

    ptn_Bcode = pd.DataFrame(index=range(1), columns=Gene_Bcode.columns)
    ptn_Bcode.iloc[0] = 0
    
    for i in range(2, len(ptn_Bcode.columns)):
        rs = ptn_Bcode.columns[i]
        homhet = ptn[ptn.rsid==rs].homhet.values[0]

        if type(homhet)==float: # NaN인 경우
            reverse = 0
            score = 100 ### Change it back to 2 when ref information for all is.na 
        elif type(homhet)!=float:

            # Is the reference reversed? 1 = reversed; 0 = not reversed
            if ptn[ptn.rsid==rs].Ref.values[0] == Gene_refATGC.iloc[:,i].values[0]:
                reverse = 0 # If NA, assume that it's not reversed.
            else:
                reverse = 1

            if (homhet == "hom") & (reverse == 1):
                score = 0
            elif (homhet == "hom") & (reverse == 0):
                score = 2
            elif (homhet == "het"):
                score = 1
            elif (homhet == "_") & (reverse == 1):
                score = 2
            elif (homhet == "_") & (reverse == 0):
                score = 0
            else:
                score = 0
        ptn_Bcode.iloc[0,i] = score

    ## CYP2D6: NA > 90% -> *5   //   NA < 90% -> Wobble(100) -> 0
    hundred_count = 0
    for i in range(len(ptn_Bcode)):
        if ptn_Bcode.iloc[0,i] == 100:
            hundred_count = hundred_count+1
            
    if ( hundred_count > (len(ptn_Bcode)*0.9) ):
        ptn_Bcode.iloc[0] = 100
    elif ( hundred_count <= (len(ptn_Bcode)*0.9) ):
        for i in range(len(ptn_Bcode.columns)):
            if ptn_Bcode.iloc[0,i]==100:
                ptn_Bcode.iloc[0,i] = 0

    ptn['two_nuc'] = '-'
    for i in range(len(ptn)):
        if ptn.homhet.iloc[i]=='hom':
            ptn.two_nuc.iloc[i]=ptn.Alt.iloc[i]+ptn.Alt.iloc[i]
        elif ptn.homhet.iloc[i]=='het':
            ptn.two_nuc.iloc[i]=ptn.Ref.iloc[i]+ptn.Alt.iloc[i]
        else:
            ptn.two_nuc.iloc[i]=ptn.Ref.iloc[i]+ptn.Ref.iloc[i]
    rsnames=ptn_Bcode.columns

    bins2 = pd.DataFrame(index=range(1), columns=range(len(rsnames)))

    for i in range(len(rsnames)):
        if len(ptn[ptn.rsid==rsnames[i]]['two_nuc'].values)>0:
            bin2 = ptn[ptn.rsid==rsnames[i]]['two_nuc'].values[0]
        else:
            bin2 = '.'
        bins2.iloc[0,i] = bin2

    bins2.columns = ptn_Bcode.columns
    ptn_Bcode_check = pd.concat([ptn_Bcode, bins2])
    
    os.chdir(logdir)
    ptn_Bcode_check.to_csv(f"{logdir}/_LOG_ptn_Bcode_{ptn_Bcode_check.columns[0]}.txt", sep='\t', index=False)

    return(ptn_Bcode)

# find_patient_star creates patient's barcode
def find_patient_star(ptn_Bcode, Gene_Bcode, outdir, ptnID, priority):

    # Find the best match between patient's barcode and gene barcode and return best prediction of patient's star number
    report = []

    for i in range(len(Gene_Bcode)):

        # Subtract to compare patient's barcode with each row of gene barcode
        p=ptn_Bcode.iloc[0,2:]
        g=Gene_Bcode.iloc[i,2:]

        # Change all the matched position with "0 / FALSE" and mismatch with "1 / TRUE"
        # See if every position had a match
        if all(p.values==g.values):
            report.append(Gene_Bcode.iloc[i,0])

    if len(report)==0 :
        report = return_best_match(ptn_Bcode, Gene_Bcode, outdir, ptnID, priority)
    return(report)

def return_best_match(ptn_Bcode, Gene_Bcode, outdir, ptnID, priority):

    dat = Gene_Bcode
    dath = dat.iloc[:,0:2]
    datc = dat.iloc[:,2:] ## star(dath), barcode(datc) split
    ptn = ptn_Bcode
    vec = ptn.iloc[:, 2:] ## ptn matrix (row number = datc)
    vec.mat = pd.DataFrame(np.repeat(vec.values, len(datc), axis = 0))
    vec.mat.columns = datc.columns

    res = abs(datc - vec.mat)
    sums = res.sum(axis=1)
    idx = np.where(sums == min(sums))
    close_dat = datc.iloc[idx].reset_index(drop=True)
    delta_h = dath.iloc[idx]
    dat_idx = dat.iloc[idx]

    common = find_highest_priority(delta_h.iloc[:,0], priority).values[0]
    
    delta = (close_dat - vec).reset_index(drop=True)
    fin_delta = pd.concat([delta_h, delta], axis=0)
    
    fin_log = pd.concat([fin_delta, dat_idx, ptn], axis=0).reset_index(drop=True) ## prn) write the log (delta / best.ref / ptn)
    os.chdir(outdir)
    gname = fin_log.columns[0]
    fin_log.to_csv(f"{outdir}/{ptnID}_LOG_{gname}_bestMatch.txt", sep='\t', index=False)

    return([fin_log.iloc[0,0]]) ## return string

# Compare other possible star numbers that have the same barcode and return one with the highest allele frequency(highest priority)
def find_highest_priority(gene_ptn_star, gene_allele_priority2):
    
    bins=[] 

    if len(gene_ptn_star) == 0:
        return("NA")
    elif len(gene_ptn_star) == 1:
        return(gene_ptn_star)
    else:
        
        for i in range(len(gene_ptn_star)): 
            vec = gene_ptn_star[i].split('_')
            score1 = np.where(gene_allele_priority2 == re.sub("[A-Z]", "", vec[0]))[0]+1
            score2 = np.where(gene_allele_priority2 == re.sub("[A-Z]", "", vec[1]))[0]+1
            sum = score1 + score2
            bins.append(sum)

        gene_ptn_star_final = [gene_ptn_star[bins.index(min(bins))]]
        return(gene_ptn_star_final)
