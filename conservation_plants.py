# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 19:50:18 2022

@author: Saubhik
"""

import os,glob,re, pickle
import numpy as np
import pandas as pd
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
Entrez.email = 'saubhiksom@iisc.ac.in'
#folder paths
path=r"nuc_alignments\\"
protein_path=r"prot_alignments\\"
if os.isdir(path) == False:
    os.makedirs(path)
if os.isdir(protein_path) == False:
    os.makedirs(protein_path)
#clearing fasta files from folder paths
files=glob.glob(path+"*.fasta")
for file in files:
    os.remove(file)

files=glob.glob(protein_path+"*.fasta")
for file in files:
    os.remove(file)

#commandline paths    
clustalw_exe=r"Add clustalw2 path"
blastn=r"NCBI/blast-BLAST_VERSION+/bin/blastn"
failed={}

#re patterns used
pattern="NM_.*|XM_.*|AT.*|NR_.*|XR_.*"
pattern_2= '(NM_[^ ]*|XM_[^ ]*|NR_[^ ]*|XR_[^ ]*)'
pattern_1="\d+"
pattern_for_AT="AT\S*"

#loading the annontation file
annotation= np.loadtxt('annotation_file.txt',delimiter="\t", dtype=object, skiprows=1)
#getting the gene_names
readthrough_df=pd.read_excel('Table S2 List of RT-positive genes.xlsx', skiprows=2)
positive_transript=[]
for i in readthrough_df['Transcript id']:
    positive_transript+= i.split(";")
while "" in positive_transript: positive_transript.remove('')

genes={}
for gene_id, gene_symbol in zip(readthrough_df['Gene id'], readthrough_df['Gene symbol']):
    if pd.notnull(gene_symbol)== False:
        gene_symbol=gene_id
    genes[gene_symbol]=gene_id

#databases used for local BLASTn
database=['a_lyrata_rna_db','brassica_rna_db','camellina_rna_db','caspella_rna_db',
          'eutrema_rna_db','raphanus_rna_db']

#getting CDS, ISR, Rem_UTR and 3'UTR sequences of all homologues from NCBI online

Homologues={}
for keys in genes:
    Homologues_for_sp={}
    print(keys)
    counter=0
    values=genes[keys]
    for fasta in SeqIO.parse('transcriptome.fasta', 'fasta'):
        header=fasta.description
        if values in header:
            sequence= str(fasta.seq)
            counter+=1
            for index, annot in enumerate(annotation[0:,0]):
                
                if annot == re.findall(pattern_for_AT, header)[0] and annot in positive_transript:
                    start=int(annotation[index,1])
                    stop1=int(annotation[index,2])
                    stop2=int(annotation[index,3])
                    cds= sequence[start-1:stop1]
                    isr=sequence[stop1:stop2]
                    rem_utr=sequence[stop2:]
                    all_3utr=sequence[stop1:]
                    assert len(isr)+len(rem_utr)==len(all_3utr), "position error"
            #adding the arabidopsis th seq retrieved by SRA
            Homologues_for_sp[header]={'CDS':cds,"ISR":isr,"Rem_UTR":rem_utr,"UTR":all_3utr}
            ###########
            fname=path+keys+'_'+str(counter)+".fasta"
            with open(fname,'w') as f:
                f.write('>'+header+'\n'+sequence+'\n')
            for sp in database:
                blastn_cline = NcbiblastnCommandline(blastn,query=fname,db=r'NCBI/blast-BLAST_VERSION+/bin/%s'%sp, outfmt=5,evalue=0.0001, out=path+"test.XML", max_target_seqs=1)
                try:
                    blastn_cline()
                except Exception as e1:
                    # print(e1)
                    failed[keys]=[sp,e1]
                    continue
                with open(path+"test.XML") as result_handle:
                    E_VALUE_THRESH=0.01
                    # blast_records = NCBIXML.parse(result_handle)
                    blast_record = NCBIXML.read(result_handle)
                    homolo_gene=homolo_seq=None
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            query_cover=hsp.align_length / blast_record.query_length*100
                            if hsp.expect < E_VALUE_THRESH and query_cover>50 :
                                # print("sequence:", alignment.title)
                                homolo_gene=re.findall(pattern, alignment.title)[0]
                                entrez_id=re.findall(pattern_2, alignment.title)[0]
                                # homolo_seq=mapping[sp][homolo_gene]
                
                os.remove(path+"test.XML")
                with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=entrez_id) as handle:
                    seq_record = SeqIO.read(handle, "gb") # using "gb" as an alias for "genbank"
               
                for features in seq_record.features:
                    if features.type == "CDS":
                      location=str(features.location)
                      CDS_location_ATG= int(re.findall(pattern_1, location)[0])
                      CDS_location_STOP=int(re.findall(pattern_1, location)[1])
                      gene_name=seq_record.description
                      CDS_seq=str(seq_record.seq[CDS_location_ATG:CDS_location_STOP])
                      UTR_seq=str(seq_record.seq[CDS_location_STOP:])
                      ISR_protein=str(Seq.translate(seq_record.seq[CDS_location_STOP:], to_stop=True))
                      ISR_seq=str(seq_record.seq[CDS_location_STOP:CDS_location_STOP+len(ISR_protein)*3+3])
                      Rem_UTR_seq=str(seq_record.seq[CDS_location_STOP+len(ISR_protein)*3+3:])
                      assert len(ISR_seq)+len(Rem_UTR_seq)==len(UTR_seq), "position error"
                      Homologues_for_sp[str(seq_record.id)+" "+gene_name.replace("PREDICTED: ","")]={"CDS":CDS_seq,
                                                                                                     "ISR":ISR_seq,
                                                                                                     'Rem_UTR':Rem_UTR_seq,
                                                                                                     "UTR":UTR_seq}
                      
                      break
            Homologues[keys+'_'+str(counter)]=Homologues_for_sp
            Homologues_for_sp={}
            os.remove(fname)    

with open('Homologues.pickle','wb') as handle:
    pickle.dump(Homologues, handle, protocol=pickle.HIGHEST_PROTOCOL)
Homologues=pickle.load(open('Homologues.pickle','rb'))

#writing the nucleotide fasta files for allignment
for gene in Homologues:
    for homolo_sp in (Homologues[gene]):
        for element in (Homologues[gene][homolo_sp]):
            nuc_seq_now=Homologues[gene][homolo_sp][element]
            if nuc_seq_now !='':
                with open(path+gene+"_nuc_"+element+'.fasta', 'a') as f:
                    f.write('>'+element+'_'+homolo_sp+'\n'+nuc_seq_now+'\n')

#writing the protein fasta files for allignment
for gene in Homologues:
    for homolo_sp in (Homologues[gene]):
        for element in (Homologues[gene][homolo_sp]):
            prot_seq_now= str(Seq.translate(Homologues[gene][homolo_sp][element],stop_symbol='X'))
            if prot_seq_now != 'X' and prot_seq_now != '' : # and protein_seq_now.count('X')!=len(protein_seq_now)
                with open(protein_path+gene+"_prot_"+element+'.fasta', 'a') as f:
                    f.write('>'+element+'_'+homolo_sp+'\n'+prot_seq_now+'\n')

                    
#clustalw nucleotide allignment (optional)
"""
files=glob.glob(path+"*.fasta")

for file in files:
            
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile= file, outfile=(file.replace(".fasta",".aln")),outorder= "INPUT")
    try:
        clustalw_cline()
    except Exception as e:
        # print(e)
        failed[file]=e
        # os.remove(file)
        continue
    # align = AlignIO.read(r'alignments/%s' %(gene_file.replace(".fasta",".aln")), "clustal")
    try:
        #os.remove(file)
        os.remove(file.replace(".fasta",".dnd"))   
    except Exception:
        continue

#clustalw protein allignment 
files=glob.glob(protein_path+"*.fasta")

for file in files:
            
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile= file, outfile=(file.replace(".fasta",".aln")),outorder= "INPUT")
    try:
        clustalw_cline()
    except Exception as e2:
        # print(e2)
        failed[file]= e2
        # os.remove(file)
        continue
    # align = AlignIO.read(r'alignments/%s' %(gene_file.replace(".fasta",".aln")), "clustal")
    try:
        #os.remove(file)
        os.remove(file.replace(".fasta",".dnd"))   
    except Exception:
        continue

""" 
