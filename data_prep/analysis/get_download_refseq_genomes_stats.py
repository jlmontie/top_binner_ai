import sys

import pandas as pd
import plotly.express as px
from tqdm import tqdm

from idbd_bio_utils import NcbiTaxonomy
sys.path.insert(0, '../../scripts')
from taxa_tools import get_class, get_genome_size

nodes = '/data/analysis_group1/jmontgomery/ai_binner/refseq/ncbi_tax_dumps_2012-02-02/nodes.dmp'
names = '/data/analysis_group1/jmontgomery/ai_binner/refseq/ncbi_tax_dumps_2012-02-02/names.dmp'
merged = '/data/analysis_group1/jmontgomery/ai_binner/refseq/ncbi_tax_dumps_2012-02-02/merged.dmp'
ncbi = NcbiTaxonomy(nodes=nodes, names=names, merged=merged)
df = pd.read_csv('../data/ncbi_download_metadata.txt', sep='\t')
df['org_class'] = df['species_taxid'].apply(get_class, ncbi=ncbi)
df['abs_path'] = df['local_filename'].str[16:]
# print(f"Calculating genome sizes")
# tqdm.pandas()
# df['genome_size'] = df['abs_path'].progress_apply(get_genome_size)
df_no_dups = df.drop_duplicates(subset=['assembly_accession'])
# df_no_dups.to_csv('../data/ncbi_download_metadata_with_stats.txt', index=False, sep='\t')

stats_file = open('../data/ncbi_download_overall_stats.txt', 'w')
org_stats = df['org_class'].value_counts()
stats_file.write(f"Organism class stats:\n")
stats_file.write(f"{org_stats}\n\n")
# print(f"Organism class stats:\n{org_stats}")
org_stats = org_stats.to_frame().reset_index()
org_stats.columns = ['organism', 'total references']
species_stats = df['species_taxid'].value_counts().to_frame().reset_index()
species_stats.columns = ['taxid', 'total references']
stats_file.write(f"Species with single reference:\n{len(species_stats[species_stats['total references'] == 1])}\n")
stats_file.write(f"Species with multiple references:\n{len(species_stats[species_stats['total references'] > 1])}\n")
stats_file.write(f"Missing species taxids:\n{(df['species_taxid'].isna().sum())}\n")
stats_file.write(f"Duplicate accessions:\n{df.duplicated(subset=['assembly_accession']).sum()}")

# fig1 = px.bar(org_stats, y='total references', x='organism', log_y=True)
# fig1.write_html('../plots/references_per_organism.html')

# fig2 = px.scatter(species_stats.sort_values(['total references'], ascending=False),
#                   y='total references', x='taxid', log_y=True)
# fig2.update_xaxes(showticklabels=False)
# fig2.write_html('../plots/references_per_taxa.html')

