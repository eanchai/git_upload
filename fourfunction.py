"""
    date: 2022.09.26
    object: function, make dataframe for report
    author: 
"""
import re
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

class fourfunction:
    def __init__(self, Gene_Bcode: str, ptn: str, Gene_refATGC: str, logdir: str, ptnID: str, pritority: str, outdir: str, gene_ptn_star, gene_allele_priority2) -> pd.DataFrame:
        '''
        Class initialize
        '''
        self.Gene_Bcode = Gene_Bcode #B allele
        self.ptn = ptn  #protien
        self.Gene_refATGC = Gene_refATGC #reference 
        self.logdir = logdir
        self.ptnID = ptnID
        self.pritority = pritority #우선순위
        self.outdir = outdir
        self.gene_ptn_star = gene_ptn_star
        self.gene_allele_priority2 = gene_allele_priority2
        
# Generate patient's barcode by insering either 1("het") or 2("hom") to correct position(rs number)
    def make_patient_barcode(self) -> pd.DataFrame:
        ptn_Bcode = pd.DataFrame(index=range(1), columns=self.Gene_Bcode.columns)
        ptn_Bcode.iloc[0] = 0
        ptn_Bcode_col = ptn_Bcode.columns
        
        for ptn_length in range(2, len(ptn_Bcode_col)):
            rs = ptn_Bcode_col[ptn_length]
            homhet = self.ptn[self.ptn.rsid == rs].homhet.values[0]

            if type(homhet) == float: # NaN인 경우
                reverse = 0
                score = 100 ### Change it back to 2 when ref information for all is.na 
            else:
                # Is the reference reversed? 1 = reversed; 0 = not reversed
                homhet_dict = {
                    'hom':1,
                    'het':0,
                    '_':-1
                }
                score = reverse - homhet_dict[homhet] if reverse else 1 + homhet_dict[homhet]
            ptn_Bcode.iloc[0,ptn_length] = score

    ## CYP2D6: NA > 90% -> *5   //   NA < 90% -> Wobble(100) -> 0
        hundred_count = 0
        for ptn_Bcode_length in range(len(ptn_Bcode)):
            if ptn_Bcode.iloc[0,ptn_Bcode_length] == 100:
                hundred_count = hundred_count+1
            
        if (hundred_count > (len(ptn_Bcode)*0.9)):
            ptn_Bcode.iloc[0] = 100
        elif (hundred_count <= (len(ptn_Bcode)*0.9)):
            for i in range(len(ptn_Bcode_col)):
                if ptn_Bcode.iloc[0,i]==100:
                    ptn_Bcode.iloc[0,i] = 0

        self.ptn['two_nuc'] = '-'
        for ptn_length in range(len(self.ptn)):
            if self.ptn.homhet.iloc[ptn_length] == 'hom':
                self.ptn.two_nuc.iloc[ptn_length] = self.ptn.Alt.iloc[ptn_length] + self.ptn.Alt.iloc[ptn_length]
            elif self.ptn.homhet.iloc[ptn_length] == 'het':
                self.ptn.two_nuc.iloc[ptn_length] = self.ptn.Ref.iloc[ptn_length] + self.ptn.Alt.iloc[ptn_length]
            else:
                self.ptn.two_nuc.iloc[ptn_length] = self.ptn.Ref.iloc[ptn_length] + self.ptn.Ref.iloc[ptn_length]
        
        bins2 = pd.DataFrame(index=range(1), columns=range(len(ptn_Bcode_col)))
        for i in range(len(ptn_Bcode_col)):
            if len(self.ptn[self.ptn.rsid == ptn_Bcode_col[i]]['two_nuc'].values) > 0:
                bin2 = self.ptn[self.ptn.rsid == ptn_Bcode_col[i]]['two_nuc'].values[0]
            else:
                bin2 = '.'
            bins2.iloc[0,i] = bin2

        bins2.columns = ptn_Bcode.columns
        ptn_Bcode_check = pd.concat([ptn_Bcode, bins2])
    
        os.chdir(self.logdir)
        ptn_Bcode_check.to_csv(f"{self.logdir}/_LOG_ptn_Bcode_{ptn_Bcode_check.columns[0]}.txt", sep='\t', index=False)

        return(ptn_Bcode)

# find_patient_star creates patient's barcode
    def find_patient_star(self):

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
            report = self.return_best_match(ptn_Bcode, Gene_Bcode, outdir, ptnID, priority)
        return(report)

    def return_best_match(self):

        dat = self.Gene_Bcode
        dath = dat.iloc[:,0:2]
        datc = dat.iloc[:,2:] ## star(dath), barcode(datc) split
        ptn = self.ptn_Bcode
        vec = ptn.iloc[:, 2:] ## ptn matrix (row number = datc)
        vec.mat = pd.DataFrame(np.repeat(vec.values, len(datc), axis = 0))
        vec.mat.columns = datc.columns

        res = abs(datc - vec.mat)
        sums = res.sum(axis=1)
        idx = np.where(sums == min(sums))
        close_dat = datc.iloc[idx].reset_index(drop=True)
        delta_h = dath.iloc[idx]
        dat_idx = dat.iloc[idx]

        common = self.find_highest_priority(delta_h.iloc[:,0], priority).values[0]
    
        delta = (close_dat - vec).reset_index(drop=True)
        fin_delta = pd.concat([delta_h, delta], axis=0)
    
        fin_log = pd.concat([fin_delta, dat_idx, ptn], axis=0).reset_index(drop=True) ## prn) write the log (delta / best.ref / ptn)
        os.chdir(outdir)
        gname = fin_log.columns[0]
        fin_log.to_csv(f"{outdir}/{ptnID}_LOG_{gname}_bestMatch.txt", sep='\t', index=False)

        return([fin_log.iloc[0,0]]) ## return string

# Compare other possible star numbers that have the same barcode and return one with the highest allele frequency(highest priority)
    def find_highest_priority(self):
    
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