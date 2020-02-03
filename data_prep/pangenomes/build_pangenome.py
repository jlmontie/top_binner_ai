import glob
import os
import subprocess as sp
import textwrap
import time
import argparse
from functools import partial
from multiprocessing import Pool

import pandas as pd
import pysam


class pangenome:

    def __init__(self, fasta_list, pangenome_dir, taxid):
        self._start_time = time.time()
        self.fasta_list = fasta_list
        self._taxid = str(taxid)
        self._create_directory_structure(pangenome_dir)

    def _create_directory_structure(self, pangenome_dir):
        if not os.path.isdir(pangenome_dir):
            try:
                os.mkdir(pangenome_dir)
            except:
                pass
        self._tmp_dir = os.path.join(pangenome_dir, 'tmp')
        if not os.path.isdir(self._tmp_dir):
            try:
                os.mkdir(self._tmp_dir)
            except:
                pass
        self._parent_dir = pangenome_dir
        self._nucmer_prefix = os.path.join(self._tmp_dir, self._taxid)
        self._delta = self._nucmer_prefix + '.delta'
        self._done = self._nucmer_prefix + '.done'
        self._pangenomes_dir = os.path.join(self._parent_dir, 'pangenomes')
        self._pangenome = os.path.join(self._parent_dir, 'pangenomes', self._taxid + '_pangenome.fa')
        if not os.path.isdir(self._pangenomes_dir):
            try:
                os.mkdir(self._pangenomes_dir)
            except:
                pass
        self._query = os.path.join(self._tmp_dir, self._taxid + '_query.fa')

    def _make_query_multifasta(self, query_fasta, cover_end=False):
        with pysam.FastxFile(query_fasta) as fh, open(self._query, 'w') as query_fh:
            for entry in fh:
                sequence = entry.sequence
                i = 0
                if len(sequence) < self._chunk_size:  # if sequence is too short to get coverage
                    return
                if cover_end:
                    while i * self._chunk_size < len(sequence):
                        fragment = sequence[i * self._chunk_size:i * self._chunk_size + self._chunk_size]
                        query_fh.write(f">{entry.name}.{i}\n{fragment}\n")
                        i += 1
                else:
                    while (i + 1) * self._chunk_size < len(sequence):
                        fragment = sequence[i * self._chunk_size:i * self._chunk_size + self._chunk_size]
                        query_fh.write(f">{entry.name}.{i}\n{fragment}\n")
                        i += 1

    def _copy_reference(self, reference_path):
        if reference_path.endswith('.gz'):
            sp.call(f"install -C -m 777 -o jmontgomery {reference_path} {self._pangenomes_dir}", shell=True)
            sp.call(f"mv {os.path.join(self._pangenomes_dir, os.path.basename(reference_path))} {self._pangenome + '.gz'}", shell=True)
            sp.call(f"gunzip -f {self._pangenome + '.gz'}", shell=True)
        else:
            sp.call(f"install -C -m 777 -o jmontgomery {reference_path} {self._pangenomes_dir}", shell=True)
            sp.call(f"mv {os.path.basename(reference_path)} {self._pangenomes_dir}", shell=True)

    def _prep_reference_and_query(self):
        reference_path = self.fasta_list[0]
        self._copy_reference(reference_path)
        with pysam.FastxFile(reference_path) as ref_fh:
            for entry in ref_fh:
                self._reference_name = entry.name
                break
        self._query_list = self.fasta_list[1:]

    def _nucmer_alignment(self):
        sp.call(f"nucmer -p {self._nucmer_prefix} {self._pangenome} {self._query}", shell=True)

    def _parse_delta_file(self):
        if not os.path.isfile(self._delta):
            print(f'No delta file created for genome {os.path.basename(self._current_query)}. Skipping.')
            return 1
        with open(self._delta) as delta_file:
            mismatching_queries = []
            for line in delta_file:
                if self._reference_name in line:
                    query_name = line.split(' ')[1]
                    data = next(delta_file).split(' ')
                    identity_pct = 1 - (100 * int(data[4]) / self._chunk_size)
                    length_pct = 100 * ((1 + abs(int(data[2]) - int(data[3]))) / self._chunk_size)
                    if identity_pct < self._max_identity_pct and length_pct < self._max_length_pct:
                        mismatching_queries.append(query_name)
        self._mismatching_queries = mismatching_queries
        return 0

    def _update_pangenome(self):
        with pysam.FastxFile(self._query) as fh, open(self._pangenome, 'a') as pan_fh:
            for entry in fh:
                if entry.name in self._mismatching_queries:
                    pan_fh.write(textwrap.fill(entry.sequence, 80))
                    pan_fh.write('\n')

    def _cleanup(self):
        if os.path.isfile(self._delta):
            os.remove(self._delta)
        if os.path.isfile(self._query):
            os.remove(self._query)
        with open(self._done, 'w') as done_file:
            done_file.write(f"Finished processing {len(self._query_list)} fastas in {time.time() - self._start_time}s")

    def create_pangenome(self, max_identity_pct=90, max_length_pct=70,
            chunk_size=500):
        self._chunk_size = chunk_size
        self._max_identity_pct = max_identity_pct
        self._max_length_pct = max_length_pct
        self._prep_reference_and_query()
        for idx, query in enumerate(self._query_list):
            self._current_query = query
            print(f"\nTaxid {self._taxid}, comparing reference to {os.path.basename(query)} ({idx + 1}/{len(self._query_list)}).\n")
            self._make_query_multifasta(query)
            self._nucmer_alignment()
            status = self._parse_delta_file()
            if status == 1:
                continue
            print(f"\nTaxid {self._taxid}, {len(self._mismatching_queries)} chunks added to pangenome\n")
            self._update_pangenome()
            print(f"\nTaxid {self._taxid}, query {os.path.basename(query)} complete. Taxid {100 * (idx + 1) / len(self._query_list):0.2f}% complete.\n")
        self._cleanup()


def build_pangenome(group, pangenome_dir=None, skip_existing=True):
    taxid = group[0]
    group_df = group[1]
    group_df = group_df.sort_values(['genome_size'], ascending=False)
    fasta_list = group_df['abs_path'].to_list()
    if not all(group_df['org_class'] == 'bacteria'):
        return
    if not len(fasta_list) > 1:
        return
    if skip_existing:
        done_file = os.path.join(pangenome_dir, 'tmp', str(taxid)) + '.done'
        if os.path.isfile(done_file):
            print(f"Taxid {taxid} complete. Skipping.")
            return
    pan = pangenome(fasta_list, pangenome_dir, taxid)
    pan.create_pangenome(max_identity_pct=90, max_length_pct=70,
        chunk_size=500)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('Build pangenomes from'
        'ncbi_download_metadata file.'))
    parser.add_argument('metadata_file', type=str,
        help=('Path to ncbi_download_metadata file. Must include abs_path and'
            'genome_size columns.'))
    parser.add_argument('pangenome_dir', type=str, help=('Directory to save'
        'pangenomes and temporary files.'))
    parser.add_argument('threads', type=int,
        help='Parallel processing threads.')
    args = parser.parse_args()
    metadata = pd.read_csv(args.metadata_file, sep='\t')
    meta_group = metadata.groupby('species_taxid')
    p = Pool(args.threads)
    p.map(partial(build_pangenome, pangenome_dir=args.pangenome_dir), meta_group)
    p.close()
    p.join()
