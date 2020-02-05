# Data Prep Steps
Below is the sequence of steps performed to prepare refseq genomes for a
simulated read database and for word embedding. This work was performed on
asimov and all absolute paths refer to this server.

## Download refseq fastas from NCBI
Download all refseq genomes from NCBI for bacteria, fungi, viral, protozoa.
Genomes downloaded with
[ncbi-genome-download](https://pypi.org/project/ncbi-genome-download/).
All but contigs downloaded (i.e., complete genomes, chromosomes, and scaffolds).

**command**: see `commands/download_refseq.sh`

**output**:
1. downloaded refseq genomes saved on asimov in
`/data/analysis_group1/jmontgomery/ai_binner/refseq/genomes`
2. medata describing downloaded fastas - `data/ncbi_download_metadata.txt`

**analysis**: `analysis/get_download_refseq_genomes_stats.py`.

This script analyzes the metadata file `data/ncbi_download_metadata.txt` to
determine the organism class and genome size for each fasta. The metdata file
was rewritten with additional columns to the file below.

**analysis output**:

1. modified metadata table - `data/ncbi_download_metadata_with_stats.txt`
2. summary of download stats - `data/ncbi_download_overall_stats.txt`
3. plots of genomes per reference and species taxid
    - `plots/references_per_organism.html`
    - `plots/references_per_taxa.html`

## Pangenomes
Genomes are reduced by creating a pangenome for all bacterial species that have
more than one genome fasta.

**command**: see `commands/build_pangenomes.sh`

**output**: pangenomes saved on asimov in
`/data/analysis_group1/jmontgomery/ai_binner/pangenomes/fastas`.

**analysis**: `analysis/get_pangenome_stats.py` and
`AnalyzePangenomeStats.ipynb`

This analysis first runs `get_pangenome_stats.py`. This outputs a file which is
consumed by the `AnalyzePangenomeStats.ipynb` notebook. (This analysis
was separated int two stages because of the lengthy compute time of the first
script.) The notebook shows the compression effect of pangenome creation.

**analysis output**:
1. File with paths to pangenomes and single-reference fastas. This file is used
for reference in building the simulated read database and the word embedding
corpuses.
2. `pangenome_compression_stats_.csv`. This is the stats file consumed by the
notebook with stats on the compression effect.

## Genome clustering
K-means clustering reduces the genomes to 1,200 of the most diverse and
representative genomes. Clustering is performed on a jaccard index matrix with
1,200 nodes and the nodes are selected for the reduced data set.

**command**: 


## Simulated read database

## Word embedding
Further data preparation (i.e. corpus preparation) are covered further in the
`word_embedding` directory. See README in that location for more details.