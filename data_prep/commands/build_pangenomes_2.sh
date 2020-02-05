SCRIPT=/home/jmontgomery/top_binner_ai/data_prep/pangenomes/build_pangenome.py
METADATA_FILE=/home/jmontgomery/top_binner_ai/data_prep/data/ncbi_download_metadata_with_stats.txt
PANGENOME_DIR=/data/analysis_group1/jmontgomery/ai_binner/pangenomes/
THREADS=72

python $SCRIPT $METADATA_FILE $PANGENOME_DIR $THREADS --skip_inprogress

# nohup bash build_pangenomes_2.sh > build_pangenomes_2.out 2>&1