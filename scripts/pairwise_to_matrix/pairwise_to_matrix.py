import argparse
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import time
import os
from scipy.sparse import csr_matrix

class EfficientSimilarityMatrix:
    @staticmethod
    def load_csv_to_sparse_matrix(csv_file_path, similarity_column):
        print("Loading CSV file...")
        start_time = time.time()

        dask_df = dd.read_csv(csv_file_path, usecols=['query_md5', 'match_md5', similarity_column])
        dask_df = dask_df.categorize(columns=['query_md5', 'match_md5'])
        dask_df['query_code'] = dask_df['query_md5'].cat.codes
        dask_df['match_code'] = dask_df['match_md5'].cat.codes

        with ProgressBar():
            matrix_data = dask_df.compute()

        rows = matrix_data['query_code'].values
        cols = matrix_data['match_code'].values
        data = matrix_data[similarity_column].values

        num_items = max(rows.max(), cols.max()) + 1
        sparse_matrix = coo_matrix((data, (rows, cols)), shape=(num_items, num_items))

        elapsed_time = time.time() - start_time
        print(f"CSV file loaded in {elapsed_time:.2f} seconds.")

        return sparse_matrix, matrix_data['query_md5'].cat.categories

    @staticmethod
    def save_sparse_matrix_to_hdf5(sparse_matrix, md5_categories, hdf5_file_path):
        print(f"Saving data to HDF5 file: {hdf5_file_path}")
        start_time = time.time()

        with h5py.File(hdf5_file_path, 'w') as f:
            f.create_dataset('data', data=sparse_matrix.data)
            f.create_dataset('indices', data=sparse_matrix.col)
            f.create_dataset('indptr', data=sparse_matrix.row)
            md5_array = np.array(md5_categories.tolist(), dtype=h5py.string_dtype())
            f.create_dataset('md5_mapping', data=md5_array)

        elapsed_time = time.time() - start_time
        print(f"HDF5 file saved in {elapsed_time:.2f} seconds.")

    @staticmethod
    def load_hdf5_to_sparse_matrix(hdf5_file_path):
        print(f"Loading HDF5 file: {hdf5_file_path}")
        start_time = time.time()

        with h5py.File(hdf5_file_path, 'r') as f:
            data = f['data'][:]
            indices = f['indices'][:]
            indptr = f['indptr'][:]
            md5_mapping = f['md5_mapping'][:]

        sparse_matrix = coo_matrix((data, (indptr, indices)))

        elapsed_time = time.time() - start_time
        print(f"HDF5 file loaded in {elapsed_time:.2f} seconds.")
        return sparse_matrix, md5_mapping

    @staticmethod
    def write_sparse_matrix_to_tsv(sparse_matrix, md5_mapping, tsv_file_path, batch_size=1000):
        print(f"Writing data to TSV file: {tsv_file_path}")
        start_time = time.time()

        sparse_matrix_csr = csr_matrix(sparse_matrix)

        with open(tsv_file_path, 'w') as f:
            header = '\t'.join(md5_mapping.astype(str)) + '\n'
            f.write(header)

            num_rows = sparse_matrix_csr.shape[0]
            for i in range(0, num_rows, batch_size):
                end = min(i + batch_size, num_rows)
                rows = sparse_matrix_csr[i:end].toarray()
                for row in rows:
                    row_str = '\t'.join(map(str, row)) + '\n'
                    f.write(row_str)

        elapsed_time = time.time() - start_time
        print(f"TSV file written in {elapsed_time:.2f} seconds.")


def file_exists(filepath):
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"The file '{filepath}' does not exist.")

def main():
    parser = argparse.ArgumentParser(description="Efficient Similarity Matrix Processing")
    subparsers = parser.add_subparsers(dest="command", required=True)

    parser_hdf5 = subparsers.add_parser("to_hdf5")
    parser_hdf5.add_argument("-p", "--pairwise-csv", type=str, required=True, help="Path to the pairwise CSV file")
    parser_hdf5.add_argument("-m", "--metric", type=str, required=True, choices=['containment', 'max_containment', 'jaccard', 'intersect_hashes'], help="Similarity metric column in CSV")
    parser_hdf5.add_argument("-o", "--output-prefix", type=str, required=True, help="Output prefix for HDF5 file")

    parser_tsv = subparsers.add_parser("to_tsv")
    parser_tsv.add_argument("--hdf5", type=str, required=True, help="Path to the HDF5 file")
    parser_tsv.add_argument("--output_tsv", type=str, required=True, help="Output path for TSV file")

    args = parser.parse_args()

    if args.command == "to_hdf5":
        file_exists(args.pairwise_csv)
        matrix, md5_categories = EfficientSimilarityMatrix.load_csv_to_sparse_matrix(args.pairwise_csv, args.metric)
        EfficientSimilarityMatrix.save_sparse_matrix_to_hdf5(matrix, md5_categories, args.output_prefix + ".hdf5")

    elif args.command == "to_tsv":
        file_exists(args.hdf5)
        matrix, md5_mapping = EfficientSimilarityMatrix.load_hdf5_to_sparse_matrix(args.hdf5)
        EfficientSimilarityMatrix.write_sparse_matrix_to_tsv(matrix, md5_mapping, args.output_tsv)

if __name__ == "__main__":
    main()
