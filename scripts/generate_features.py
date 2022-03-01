import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


length, mw, ar = [],[],[]
ii, pi, nc, gr = [],[],[],[]
ids, fam = [],[]


for record in SeqIO.parse('AMPSphere_v.2021-03.faa', 'fasta'):
    print(f'... processing record {record.id}')
    ptn = ProteinAnalysis(str(record.seq))
    ids.append(record.id)
    fam.append(record.description.split(' | ')[1])
    length.append(len(record.seq))
    mw.append(ptn.molecular_weight())
    ar.append(ptn.aromaticity())
    ii.append(ptn.instability_index())
    pi.append(ptn.isoelectric_point())
    nc.append(ptn.charge_at_pH(pH=7.0))
    gr.append(ptn.gravy())

print(f'... creating data frame')
df = pd.DataFrame(np.asarray([ids, fam, length, mw,
                              ar, ii, pi, nc, gr]).T,
                  columns=['id', 'family', 'length',
                           'molecular_weight', 'aromaticity',
                           'instability_index', 'isoelectric_point',
                           'charge', 'gravy'])

print(f'... exporting data')
df.to_csv('features_plot_ampsphere.tsv', sep='\t', header=True, index=None)
