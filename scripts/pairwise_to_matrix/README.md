# Pairwise Command Output to Sparse/Dense Matrix

## Description
Efficient processing of large pairwise comparison CSV files to HDF5 and exporting HDF5 files containing sparse matrices to TSV.


## Installation

- Libraries: `pandas`, `dask`, `numpy`, `scipy`, `h5py`

```
pip install pandas dask numpy scipy h5py
```

## Usage
**Convert CSV to HDF5:**

```
python pairwise_to_matrix.py to_hdf5 -p <path_to_csv> -m <metric> -o <output_prefix>
```

- `<path_to_csv>`: Path to the CSV file.
- `<metric>`: Similarity metric (choose from `containment`, `max_containment`, `jaccard`, `intersect_hashes`).
- `<output_prefix>`: Prefix for the output HDF5 file.

**Convert HDF5 to TSV:**
```
python pairwise_to_matrix.py to_tsv --hdf5 <path_to_hdf5> --output_tsv <output_tsv>
```

- `<path_to_hdf5>`: Path to the HDF5 file.
- `<output_tsv>`: Path for the output TSV file.
