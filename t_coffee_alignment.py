#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 11:50:37 2022

@author: Saubhik
"""
import os, glob, re, subprocess, pickle
import pandas as pd

def t_coffee_alignment(directory,output_file_name,space_rename=1):
    """
    This function does t-coffee alignment with following parameters

    Parameters
    ----------
    directory : str
        Give a directory path where the fasta files for the allignment exists.
    output_excel_file_name : str
        Give the name of output excel file without extension.
    space_rename: bool
        True if fasta file names contain space
        False if fasta file names does not contain space

    Returns
    -------
    An excel file in the input directory with t_coffee alignment scores.

    """
    
    os.chdir(directory)
    
    #This part is only needed to remove  space from your fasta file name 
    if bool(space_rename)== True:
        files=[file for file in glob.glob("*.fasta")]
        for gene_file in files:
            new_name=gene_file.replace(' ','_')
            gene_file=os.rename(gene_file,new_name)
    
    files=glob.glob("*.fasta")
    all_scores=[]
    failed={}
    for gene_file in files:
        result=subprocess.run(['t_coffee',gene_file], capture_output=True, text=True)
        error=result.stderr
    
        pattern="SCORE=\d*"
        pattern_1="Nseq=\d*"
        try:
            with open(gene_file.replace('.fasta','.aln'),'r') as f:
                for idx, i in enumerate(f.readlines()):
                    if idx==0:
                        print(i)
                        break
        except Exception as e:
            failed[gene_file.replace('.fasta',"")]=error
            continue
        try:
            os.remove(gene_file.replace('.fasta','.html'))
            os.remove(gene_file.replace('.fasta','.dnd'))     
        except Exception as e:
            failed[gene_file.replace('.fasta',"")]=error
            continue
        score=re.findall(pattern, i)[0]
        score=int(score.split('=')[1])
        Nseq=re.findall(pattern_1, i)[0]
        Nseq=int(Nseq.split('=')[1])
        all_scores.append([gene_file.replace('.fasta',""),score,Nseq])
    
    with open(directory+'//'+output_file_name+'failed.pickle', 'wb') as handle:
        pickle.dump(failed, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    all_scores=pd.DataFrame(all_scores,columns=['gene_name','alignment_score','number_seq'])
    
    with open(directory+'//'+output_file_name+'.pickle', 'wb') as handle:
        pickle.dump(all_scores, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    excel={}
    for idx, name in enumerate(all_scores.gene_name):
        if "CDS" in name:
            name=name.rstrip("_nuc_CDS")
            if name not in list(excel.keys()):
                excel[name]=[]
            score=list(all_scores.alignment_score)[idx]
            number_seq=list(all_scores.number_seq)[idx]
            excel[name]=excel[name]+[score,number_seq]
        if "ISR" in name:
            name=name.rstrip("_nuc_ISR")
            if name not in list(excel.keys()):
                excel[name]=[]
            score=list(all_scores.alignment_score)[idx]
            number_seq=list(all_scores.number_seq)[idx]
            excel[name]=excel[name]+[score,number_seq]
        if "Rem_UTR" in name:
            name=name.rstrip("_nuc_Rem_UTR")
            if name not in list(excel.keys()):
                excel[name]=[]
            score=list(all_scores.alignment_score)[idx]
            number_seq=list(all_scores.number_seq)[idx]
            excel[name]=excel[name]+[score,number_seq]
        if "UTR" in name and "Rem_UTR" not in name :
            name=name.rstrip("_nuc_UTR")
            if name not in list(excel.keys()):
                excel[name]=[]
            score=list(all_scores.alignment_score)[idx]
            number_seq=list(all_scores.number_seq)[idx]
            excel[name]=excel[name]+[score,number_seq]

    columns=['cds_score','cds_nseq','isr_score','isr_nseq','rem_utr_score','rem_utr_nseq','3-utr_score','3-utr_nseq']
    excel_df=pd.DataFrame.from_dict(excel, orient='index' )
    excel_df.columns=columns
    excel_df.to_excel(directory+'//'+output_file_name+'.xlsx')

    
    return failed
    
