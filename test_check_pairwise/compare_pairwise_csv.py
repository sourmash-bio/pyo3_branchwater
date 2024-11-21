#! /usr/bin/env python
import sys
import argparse
import csv
import numpy
from collections import defaultdict

from sourmash.logging import debug_literal, error, notify, print_results
from sourmash import sourmash_args

def main():
    p = argparse.ArgumentParser()
    p.add_argument('csv1')
    p.add_argument('csv2')
    p.add_argument('-u', '--use-column', default='jaccard')
    p.add_argument('-R', '--round-to', default=10, type=int)
    args = p.parse_args()

    notify(f"loading '{args.csv1}'")
    with sourmash_args.FileInputCSV(args.csv1) as r:
        rows1 = list(r)
    notify(f"loaded {len(rows1)} rows.")

    notify(f"loading '{args.csv2}'")
    with sourmash_args.FileInputCSV(args.csv2) as r:
        rows2 = list(r)
    notify(f"loaded {len(rows2)} rows.")

    queries1 = set(( row["query_name"] for row in rows1 ))
    matches1 = set(( row["match_name"] for row in rows1 ))

    queries2 = set(( row["query_name"] for row in rows2 ))
    matches2 = set(( row["match_name"] for row in rows2 ))

    assert queries1 == queries2, queries1 ^ queries2
    assert matches1 == matches2, matches1 ^ matches2

    colname = args.use_column

    d1 = defaultdict(dict)
    for row in rows1:
        q = row["query_name"]
        m = row["match_name"]
        value = float(row[colname])

        d1[q][m] = value

    #import pprint
    #pprint.pprint(d1)

    fail = False
    num_matching = 0
    num_missing = 0
    num_ident = 0
    num_nomatch = 0
    for row in rows2:
        q = row["query_name"]
        m = row["match_name"]
        val2 = float(row[colname])

        val1 = d1[q].get(m, -1)

        if val1 == -1:
            num_missing += 1
            fail = True
        elif round(val1, args.round_to) != round(val2, args.round_to):
            notify(f"MISMATCH: '{q}' vs '{m}', {val1} != {val2}")
            fail = True
            num_nomatch += 1
        else:
            if val1 == 1.0: # and val2 == 1.0
                num_ident += 1
            else:
                num_matching += 1

    print(f"values at 100%: {num_ident}")
    print(f"matching: {num_matching}")
    print(f"non-matching: {num_nomatch}")
    print(f"num missing: {num_missing}")

    if fail:
        notify("differences detected! failing!")
        sys.exit(-1)

if __name__ == '__main__':
    sys.exit(main())
