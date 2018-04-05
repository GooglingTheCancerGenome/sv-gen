#!/usr/bin/env python

import os
import argparse as ap
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

parser = ap.ArgumentParser(description='Provide directory structure.')
parser.add_argument('-rd','--RUNDIR',dest='run_dir',default=os.path.join(os.getcwd(),'run_test'))
parser.add_argument('-gf','--GFASTA',dest='gen_fasta',default='/Users/tschafers/Documents/reference/chr17_human.fasta')

args = parser.parse_args()
os.makedirs(args.run_dir)

gen_temp_dir = os.path.join(args.run_dir,'genomes/tmp')
gen_tumor_dir =  os.path.join(args.run_dir,'genomes/Tumor')
gen_normal_dir = os.path.join(args.run_dir,'genomes/Normal')

## Make subdirecotires
os.makedirs(gen_temp_dir)
os.makedirs(gen_tumor_dir)
os.makedirs(gen_normal_dir)

record = SeqIO.index(args.gen_fasta, "fasta")
seq_id_A_N = record['17'].id+'A_N'
seq_id_B_N = record['17'].id+'B_N'
seq_id_A_T = record['17'].id+'A_T'
seq_id_N = record['17'].id+'_N'


r_a = SeqRecord(record['17'].seq, id=seq_id_A_N)
r_b = SeqRecord(record['17'].seq, id=seq_id_B_N)
r_t = SeqRecord(record['17'].seq, id=seq_id_A_T)
r_N = SeqRecord(r_a.seq+r_b.seq,  id=seq_id_N)

SeqIO.write(r_a,gen_temp_dir+'/chr17A_N.fasta','fasta')
SeqIO.write(r_b,gen_temp_dir+'/chr17B_N.fasta','fasta')
SeqIO.write(r_t,gen_temp_dir+'/chr17A_T.fasta','fasta')
SeqIO.write(r_N,gen_normal_dir+'/chr17_N.fasta','fasta')









