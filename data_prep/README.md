# Data Prep Steps
The steps followed to prepare the refseq genome fasta files for the models.
The code for each subheading is contained within a directory of the same name.

## Download refseq fastas from NCBI
Download all refseq genomes from NCBI for bacteria, fungi, viral, protozoa.
Genomes downloaded with
[ncbi-genome-download](https://pypi.org/project/ncbi-genome-download/).
All but contigs downloaded (i.e., complete genomes, chromosomes, and scaffolds).

**command**: see `commands/download_refseq.sh`

**output**: `data/ncbi_download_metadata.txt`.

**analysis**: `data/ncbi_download_metadata.txt` was analyzed with
`analysis/get_download_refseq_genomes_stats.py` to determine organism class
for each fasta. The metdata file was rewritten with additional columns
to `data/ncbi_download_metadata_with_stats.txt`.

## Pangenomes
Genomes are reduced by creating a pangenome for all bacterial species that have
more than one genome fasta

