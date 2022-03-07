import pandas as pd

# The files needed to generate the resource "AMPSphere_ProGenomes2.tsv.gz"
# are available in:
# ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/resource_table/pgenomes.tar.gz
data = pd.read_table('AMPSphere_ProGenomes2.tsv.gz',
                     sep='\t',
                     header='infer')

# The files needed to generate the resource "pgenomes_samples.tsv"
# can be downloaded from:
# wget -O pgenomes_samples.tsv http://progenomes.embl.de/data/proGenomes2.1_specI_clustering.tab
ref = pd.read_table('pgenomes_samples/pgenomes_samples.tsv',
                    sep='\t',
                    header=None,
                    names=['specI', 'genome'])

df1 = data.merge(on='genome', right=ref)
df2 = data[~data.AMP.isin(df1.AMP)]
df = pd.concat([df1, df2]).sort_values(by='AMP')

df.to_csv('pgenomes_AMP_specI.tsv',
          sep='\t',
          header=True,
          index=None)

del data, ref, df1, df2

# The files needed to generate the resource "complete_amps_associated_taxonomy.tsv.gz"
# are available in:
# ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/resource_table/create_complete_amps_associated_taxonomy.py
df2 = pd.read_table('complete_amps_associated_taxonomy.tsv.gz',
                    sep='\t',
                    header='infer')

df.rename({'GMSC10': 'gmsc', 'AMP': 'amp', 'genome': 'sample', 'species': 'source'}, axis=1, inplace=True)
df2.rename({'name': 'source'}, axis=1, inplace=True)

df = df[['gmsc', 'amp', 'sample', 'source', 'specI']]
df2 = df2[['gmsc', 'amp', 'sample', 'source', 'specI']]

df['is_metagenomic'] = [False for x in range(len(df))]
df2['is_metagenomic'] = [True for x in range(len(df2))]

gmsc_genes = pd.concat([df, df2])
gmsc_genes = gmsc_genes.fillna('*')
gmsc_genes = gmsc_genes.sort_values(by='gmsc').reset_index(drop=True)

gmsc_genes.to_csv('complete_gmsc_pgenomes_metag.tsv',
                  sep='\t',
                  header=True,
                  index=None)

del df, df2

# The files needed to generate the resource "metadata_hreadable/reduced_metadata.tsv"
# and the table are available in:
# ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/resource_table/metadata_hreadable.tar.gz
metadata = pd.read_table('metadata_hreadable/reduced_metadata.tsv',
                         sep='\t',
                         header='infer')

metadata.rename({'sample_accession': 'sample'}, axis=1, inplace=True)

df = gmsc_genes.merge(on='sample', right=metadata)
df2 = gmsc_genes[~gmsc_genes.gmsc.isin(df.gmsc)]

gdf = pd.concat([df, df2])
gdf = gdf.sort_values(by='gmsc')
gdf = gdf.reset_index(drop=True)

gdf.to_csv('gmsc_amp_genes_envohr_source.tsv',
           sep='\t',
           header=True,
           index=None)
