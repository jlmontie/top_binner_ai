import argparse
import json
import multiprocessing as mp
import os
import subprocess as sp

import numpy as np
from tqdm import tqdm
from scipy.spatial.distance import cosine
from sklearn.cluster import KMeans


def triangle_to_matrix(triangle_filepath):
    with open(triangle_filepath) as tr_file:
        tr_file.readline()
        data_ls = []
        for line in tr_file:
            data_ls.append([float(item) if item != '-' else 0
                            for item in line.strip().split('\t')[1:]])
        data = np.array(data_ls)
        data = data.T + data
        np.fill_diagonal(data, 1)
    return data


def run_dashing(paths_file, tmp_dir):
    if tmp_dir is None:
        tmp_dir = 'tmp'
    distances_file = os.path.join(tmp_dir, 'distances.txt')
    sizes_file = os.path.join(tmp_dir, 'sizes.txt')
    sp.call(f"dashing cmp -k 25 -p 48 -F {paths_file} -O {distances_file} -o {sizes_file}", shell=True)
    jaccard_matrix = triangle_to_matrix(distances_file)
    return jaccard_matrix


def get_centroids(jaccard_matrix, max_nodes):
    if jaccard_matrix.shape[0] < max_nodes:
        n_clusters = jaccard_matrix.shape[0]
    else:
        n_clusters = max_nodes
    neigh = KMeans(n_clusters=n_clusters)
    neigh.fit(jaccard_matrix)
    centroid_samples_idx = []
    for centroid_idx in range(neigh.cluster_centers_.shape[0]):
        centroid_position = neigh.cluster_centers_[centroid_idx, :]
        cos_dist_ls = []
        for jaccard_matrix_idx in range(jaccard_matrix.shape[0]):
            sample_jaccard = jaccard_matrix[jaccard_matrix_idx, :]
            cos_dist_ls.append(cosine(centroid_position, sample_jaccard))
        closest_sample_idx = np.argmin(cos_dist_ls)
        centroid_samples_idx.append(closest_sample_idx)
    return centroid_samples_idx


def reduce_genus(paths_file, corpus_paths_filepath, max_nodes, tmp_dir):
    with open(paths_file) as infile:
        paths_ls = []
        for line in infile:
            paths_ls.append(line.strip())
    jaccard_matrix = run_dashing(paths_file, tmp_dir)
    centroid_samples_idx = get_centroids(jaccard_matrix, max_nodes)
    centroid_samples_paths = np.array(paths_ls)[centroid_samples_idx]
    with open(corpus_paths_filepath, 'w') as outfile:
        for path in centroid_samples_paths:
            outfile.write(f"{path}\n")


def main(paths_file, output, max_nodes, tmp_dir):
    reduce_genus(paths_file, output, max_nodes, tmp_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reduce genomes within a taxid by K-means clustering of jaccard similarities.')
    parser.add_argument('paths_file', type=str, help='Path to text file with fasta paths.')
    parser.add_argument('output', type=str, help='Path to file where reduced corpus genome paths will be written.')
    parser.add_argument('max_nodes', type=int, help='(Optional) Number of clustering nodes or total genomes to return.')
    parser.add_argument('-t', '--temp_dir', type=str, help='(Optional) Directory to write temporary files. Default creates tmp in current working directory.')

    args = parser.parse_args()
    main(args.paths_file, args.output, args.max_nodes, args.temp_dir)
