import glob
import gzip
import os
from collections import defaultdict

import pysam
from tqdm import tqdm
import pandas as pd

from idbd_bio_utils import NcbiTaxonomy

dump_dir = '/data/analysis_group1/jmontgomery/ai_binner/refseq/ncbi_tax_dumps_2012-02-02/'
merged = os.path.join(dump_dir, 'merged.dmp')
nodes = os.path.join(dump_dir, 'nodes.dmp')
names = os.path.join(dump_dir, 'names.dmp')
ncbi = NcbiTaxonomy(merged=merged, nodes=nodes, names=names)
pangenome_dir = '/data/analysis_group1/jmontgomery/ai_binner/pangenomes/fastas'
fasta_paths_file_dir = '/home/jmontgomery/ai_binner/data_prep/fasta_paths_dashing_filtered'

fasta_paths_file_ls = os.listdir(fasta_paths_file_dir)
genome_stats = defaultdict(list)
for fasta_paths_file in tqdm(fasta_paths_file_ls):
    taxid = fasta_paths_file.split('_')[0]
    org_name = ncbi.get_name(int(taxid))
    with open(os.path.join(fasta_paths_file_dir, fasta_paths_file)) as file:
        ref_genome = file.readline().strip()
        genomes = 1
        for line in file:
            genomes += 1
    with pysam.FastxFile(ref_genome, 'rt') as ref_h:
        original_genome_len = 0
        for entry in ref_h:
            original_genome_len += len(entry.sequence)
    pangenome = glob.glob(os.path.join(pangenome_dir, f"{taxid}_pangenome.fa"))[0]
    with pysam.FastxFile(pangenome) as pan_h:
        pangenome_len = 0
        for entry in pan_h:
            pangenome_len += len(entry.sequence)
    genome_stats['taxid'].append(taxid)
    genome_stats['Organism'].append(org_name)
    genome_stats['Original genome size'].append(original_genome_len)
    genome_stats['Pangenome size'].append(pangenome_len)
    genome_stats['Pangenome expansion factor'].append(pangenome_len / original_genome_len)
    genome_stats['taxid reference genomes'].append(genomes)

genome_stats_df = pd.DataFrame(genome_stats)
genome_stats_df.to_csv('pangenome_compression_stats_.csv', index=False)
