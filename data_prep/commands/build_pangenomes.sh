SCRIPT=/home/jmontgomery/ai_binner/data_prep_clean_start/pangenomes/build_pangenome.py
METADATA_FILE=/home/jmontgomery/ai_binner/data_prep_clean_start/data/ncbi_download_metadata_with_stats.txt
PANGENOME_DIR=/data/analysis_group1/jmontgomery/ai_binner/pangenomes/
THREADS=24

python $SCRIPT $METADATA_FILE $PANGENOME_DIR $THREADS

# nohup bash build_pangenomes.sh > build_pangenomes.out 2>&1