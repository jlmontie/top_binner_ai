OUTPUT=/data/analysis_group1/jmontgomery/ai_binner/refseq/genomes/
META=/home/jmontgomery/top_binner_ai/data_prep/data/ncbi_download_metadata.txt

ncbi-genome-download bacteria,fungi,viral,protozoa -F fasta -l complete,chromosome,scaffold -o $OUTPUT --flat-output -p 24 -r 2 -m $META -v

# nohup bash download_refseq.sh > download_refseq.out 2>&1