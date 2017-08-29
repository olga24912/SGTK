from Bio import SeqIO
import seaborn as sns
import warnings
import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels as sm
import matplotlib
import matplotlib.pyplot as plt
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

break_len=100
seq_records=[]

def splitRecord(record):
    global seq_records
    print(record)
    lens=[]
    slen = len(record)
    print(slen)
    cur = 0
    loc = 500
    scale = 25264235
    while cur + 500 < slen:
        l = scale*st.powerlognorm.rvs(c=53, s=4, size=1)[0]
        l += 500
        l = int(l)
        if (cur + l <= slen):
            cur += l
            lens.append(l)
            seq_records.append(SeqRecord(Seq(str(record.seq[cur - l: cur])),
                           id=(record.id + str(cur-l)), name=(record.id + str(cur-l)),
                                     description=""))
        
        cur += break_len

    

for seq_record in SeqIO.parse("/home/olga/bio-project/data/art_data/ref.fasta", "fasta"):
    splitRecord(seq_record)

random.shuffle(seq_records)
SeqIO.write(seq_records, "/home/olga/bio-project/data/art_data/gen_contigs.fasta", "fasta")

