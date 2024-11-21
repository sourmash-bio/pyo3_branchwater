#! /usr/bin/env python
import sys
import argparse
import csv
import numpy

from sourmash.logging import debug_literal, error, notify, print_results
from sourmash import sourmash_args

def load_labelinfo_csv(filename):
    "Load file output by 'sourmash compare --labels-to'"
    with sourmash_args.FileInputCSV(filename) as r:
        labelinfo = list(r)

    labelinfo.sort(key=lambda row: int(row["sort_order"]))
    return labelinfo


def main():
    p = argparse.ArgumentParser()
    p.add_argument('mat', help='numpy similarity matrix')
    p.add_argument('labels_from', help='CSV file output by --labels-to')
    p.add_argument('-o', '--output-pairwise-csv', required=True,
                   help="save CSV to this file")
    p.add_argument('-u', '--use-column', default='jaccard')
    args = p.parse_args()

    with open(args.mat, "rb") as fp:
        mat = numpy.load(fp)
    notify(f"...got {mat.shape[0]} x {mat.shape[1]} matrix.", *mat.shape)

    labelinfo = load_labelinfo_csv(args.labels_from)

    n_rows = 0
    with open(args.output_pairwise_csv, "w", newline="") as outfp:
        w = csv.writer(outfp)
        w.writerow(["query_name", "match_name", args.use_column])
        for i in range(mat.shape[0]):
            for j in range(i + 1):
                val = mat[i][j]
                if val > 0.0:
                    w.writerow([labelinfo[j]["label"],
                                labelinfo[i]["label"],
                                val])
                    n_rows += 1

    notify(f"wrote {n_rows} rows to '{args.output_pairwise_csv}'")


if __name__ == '__main__':
    sys.exit(main())
