# fastgather, fastmultigather, and manysearch - an introduction

This repository implements four sourmash plugins, `fastgather`, `fastmultigather`, `multisearch`, and `manysearch`. These plugins make use of multithreading in Rust to provide very fast implementations of `search` and `gather`. With large databases, these commands can be hundreds to thousands of times faster, and 10-50x lower memory. 

The main *drawback* to these plugin commands is that their inputs and outputs are not as rich as the native sourmash commands. In particular, this means that input databases need to be prepared differently, and the output may be most useful as a prefilter in conjunction with regular sourmash commands.

## Preparing the database

All four commands use
_text files containing lists of signature files_, or "fromfiles", for the search database, and `multisearch`, `manysearch` and `fastmultigather` use "fromfiles" for queries, too.

(Yes, this plugin will eventually be upgraded to support zip files; keep an eye on [sourmash#2230](https://github.com/sourmash-bio/sourmash/pull/2230).)

To prepare a fromfile from a database, first you need to split the database into individual files:
```
mkdir gtdb-reps-rs214-k21/
cd gtdb-reps-rs214-k21/
sourmash sig split -k 21 /group/ctbrowngrp/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k21.zip -E .sig.gz
cd ..
```

and then build a "fromfile":
```
find gtdb-reps-rs214-k21/ -name "*.sig.gz" -type f > list.gtdb-reps-rs214-k21.txt
```

## Running the commands

### Running `manysearch`

The `manysearch` command finds overlaps between one or more query genomes, and one or more subject (meta)genomes. It is the core command we use for searching petabase-scale databases.

`manysearch` takes two file lists as input, and outputs a CSV:
```
sourmash scripts manysearch query-list.txt podar-ref-list.txt -o results.csv
```

To run it, you need to provide two "fromfiles" containing lists of paths to signature files (`.sig` or `.sig.gz`). If you create a fromfile as above with GTDB reps, you can generate a query fromfile like so:

```
head -10 list.gtdb-reps-rs214-k21.txt > list.query.txt
```
and then run `manysearch` like so:

```
sourmash scripts manysearch list.query.txt list.gtdb-rs214-k21.txt  -o query.x.gtdb-reps.csv -k 21 --cores 4
```

The results file here, `query.x.gtdb-reps.csv`, will have 8 columns: `query` and `query_md5`, `match` and `match_md5`, and `containment`, `jaccard`, `max_containment`, and `intersect_hashes`.

### Running `multisearch`

The `multisearch` command compares one or more query genomes, and one or more subject genomes. It differs from `manysearch` by loading everything into memory.

`manysearch` takes two file lists as input, and outputs a CSV:
```
sourmash scripts manysearch query-list.txt podar-ref-list.txt -o results.csv
```

To run it, you need to provide two "fromfiles" containing lists of paths to signature files (`.sig` or `.sig.gz`). If you create a fromfile as above with GTDB reps, you can generate a query fromfile like so:

```
head -10 list.gtdb-reps-rs214-k21.txt > list.query.txt
```
and then run `manysearch` like so:

```
sourmash scripts manysearch list.query.txt list.gtdb-rs214-k21.txt  -o query.x.gtdb-reps.csv -k 21 --cores 4
```

The results file here, `query.x.gtdb-reps.csv`, will have 8 columns: `query` and `query_md5`, `match` and `match_md5`, and `containment`, `jaccard`, `max_containment`, and `intersect_hashes`.

### Running `fastgather`

The `fastgather` command is a much faster version of `sourmash gather``.

`fastgather` takes a query metagenome and a file list as the database, and outputs a CSV:
```
sourmash scripts fastgather query.sig.gz podar-ref-list.txt -o results.csv --cores 4
```

#### Using `fastgather` to create a picklist for `sourmash gather`

One handy use case for `fastgather` is to create a picklist that can be used by `sourmash gather`. This makes full use of the speed of `fastgather` while producing a complete set of `gather` outputs.

For example, if `list.gtdb-rs214-k21.txt` contains the paths to all GTDB RS214 genomes in `sig.gz` files, as above, then the following command will do a complete gather against GTDB:

```
sourmash scripts fastgather SRR606249.trim.sig.gz \
    list.gtdb-rs214-k21.txt -o SRR606249.fastgather.csv -k 21
```

This CSV file can then be used as a picklist for `sourmash gather` like so:

```
sourmash gather SRR606249.trim.sig.gz /group/ctbrowngrp/sourmash-db/gtdb-rs214/gtdb-rs214-k21.zip \
    --picklist SRR606249.fastgather.csv:match_name:ident \
    -o SRR606249.gather.csv
```

Here the picklist should be used on a sourmash collection that contains a manifest - this will prevent sourmash from loading any sketches other than the ones in the fastgather CSV file. We recommend using zip file databases - manifests are produced automatically when `-o filename.zip` is used with `sketch dna`, and they also be prepared with `sourmash sig cat`. (If you are using a GTDB database, as above, then you already have a manifest!)

#### Example of picklist usage

A complete example Snakefile implementing the above workflow is available [in the 2023-swine-usda](https://github.com/ctb/2023-swine-usda/blob/main/Snakefile) repository. Note, it is slightly out of date at the moment!

### Running `fastmultigather`

`fastmultigather` takes a file list of query metagenomes and a file list for the database, and outputs many CSVs:
```
sourmash scripts fastmultigather query-list.txt podar-ref-lists.txt --cores 4
```

The main advantage that `fastmultigather` has over `fastgather` is that you only load the database files once, which can be a significant time savings for large databases!

#### Output files for `fastmultigather`

`fastmultigather` will output two CSV files for each query, a `prefetch` file containing all overlapping matches between that query and the database, and a `gather` file containing the minimum metagenome cover for that query in the database.

The prefetch CSV will be named `{basename}.prefetch.csv`, and the gather CSV will be named `{basename}.gather.csv`.  Here, `{basename}` is the filename, stripped of its path.

**Warning:** At the moment, if two different queries have the same `{basename}`, the CSVs for one of the queries will be overwritten by the other query. The behavior here is undefined in practice, because of multithreading: we don't know what queries will be executed when or files will be written first.

## Notes on concurrency and efficiency

Each command does things slightly differently, with implications for CPU and disk load. You can measure threading efficiency with `/usr/bin/time -v` on Linux systems, and disk load by number of complaints received when running.

(The below info is for fromfile lists. If you are using mastiff indexes, very different performance parameters apply. We will update here as we benchmark and improve!)

`manysearch` loads all the queries at the beginning, and then loads one database sketch from disk per thread. The compute-per-database-sketch is dominated by I/O. So your number of threads should be chosen with care for disk load. We typically limit it to `-c 32` for shared disks.

`multisearch` loads all the queries and database sketches once, at the beginning, and then uses multithreading to search across all matching sequences. For large databases it is extremely efficient at using all available cores. So 128 threads or more should work fine!

Like `multisearch`, `fastgather` loads everything at the beginning, and then uses multithreading to search across all matching sequences. For large databases it is extremely efficient at using all available cores. So 128 threads or more should work fine!

`fastmultigather` loads the entire database once, and then loads one query from disk per thread. The compute-per-query can be significant, though, so multithreading efficiency here is less dependent on I/O and the disk is less likely to be saturated with many threads. We suggest limiting threads to between 32 and 64 to decrease shared disk load.
