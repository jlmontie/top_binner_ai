# Top-level Binner Using AI
This is development for a top-level binner that will classify sequencing reads
using AI algorithms.

## Installation
No installation is required. The python packages used in this project are
managed with pipenv and may be installed in a pipenv virtual environment as
follows.
1. Clone the repository
2. `cd top_binner_ai`
3. `pipenv install`

## Repo Directory Structure
### Data Prep
This directory contains the code preparing data for various algorithms,
including word embedding and neural net training. The README in this directory
provides more details. Here is a summary of the development process.

1. Download refseq fastas from NCBI
2. Create pangenomes for bacterial species with more than one reference fasta.
3. Select a representative subset of genomes using clustering to further reduce
model training data.
4. Create simulated read database for model training.
