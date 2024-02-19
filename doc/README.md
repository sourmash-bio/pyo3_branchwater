# The branchwater plugin for sourmash

| command | functionality | docs |
| -------- | -------- | -------- |
| `manysketch` | Rapidly build sketches for many input files     | [link](#Running-manysketch)     |
| `fastgather` | Multithreaded `gather` of **one** metagenome against a database| [link](#Running-fastgather)
| `fastmultigather` | Multithreaded `gather` of **multiple** metagenomes against a database | [link](#Running-fastmultigather)
| `manysearch` | Multithreaded containment search for many queries in many large metagenomes | [link](#Running-manysearch)
| `multisearch` | Multithreaded comparison of multiple sketches, in memory | [link](#Running-multisearch)
| `pairwise` | Multithreaded pairwise comparison of multiple sketches, in memory | [link](#Running-multisearch)

This repository implements multithreaded plugins for [sourmash](https://sourmash.readthedocs.io/) that provide very fast implementations of `sketch`, `search`, and `gather`. These commands are typically hundreds to thousands of times faster, and 10-50x lower memory, than the current sourmash code. For example, a `gather` of SRR606249 with sourmash v4.8.6 against GTDB rs214 takes 40 minutes and 14 GB of RAM, while `fastgather` with 64 cores takes only 2 minutes and 2 GB of RAM.

The main *drawback* to these plugin commands is that their inputs and outputs are not as rich as the native sourmash commands. This means that your input files may need to be prepared differently, and the output may in some cases be most useful as a prefilter in conjunction with regular sourmash commands - see the instructions below for using `fastgather` to create picklists for sourmash.

## Input file formats

sourmash supports a variety of different storage formats for sketches (see [sourmash docs](https://sourmash.readthedocs.io/en/latest/command-line.html#choosing-signature-output-formats)), and the branchwater plugin works with some (but not all) of them. Branchwater _also_ supports an additional database type, a RocksDB-based inverted index, that is not yet supported by sourmash (through v4.8.6).

**As of v0.9.0, we recommend using zip files or manifest CSVs whenever you need to provide multiple sketches.**

| command | query input | database format |
| -------- | -------- | -------- |
| `manysketch`     | CSV with input fasta/fastq paths (details below)    | _produces_ Zip database |
| `gather`     | Single metagenome in sig, zip, manifest CSV, or fromfile     | Zip, manifest CSV, or fromfile |
| `fastmultigather` | Multiple metagenomes in sig, zip, manifest CSV, or fromfile | Zip, manifest CSV, fromfile, or rocksdb index |
| `manysearch` | Multiple genomes in sig, zip, manifest CSV, or fromfile | Zip, manifest CSV, fromfile, or rocksdb index |
| `multisearch` | Multiple sketches in sig, zip, manifest CSV, or fromfile | Multiple sketches in sig, zip, manifest CSV, or fromfile |
| `pairwise` | Multiple sketches in sig, zip, manifest CSV, or fromfile | N/A 

### Using zipfiles

When working with large collections of small sketches such as genomes, we suggest using zipfiles as produced by sourmash (e.g. using `sourmash sig cat` or `manysketch`). Zip files have a few nice features:

* sketches are compressed in zip files;
* zip files can contain many sketches, including incompatible types (e.g. multiple k-mer sizes);
* zip files contain "manifests" listing their contents;
* subsets of zip files can be efficiently selected and loaded depending on what is needed;
* in particular, _single_ sketches can be loaded on demand, supporting lower memory requirements for certain kinds of searches.

For all these reasons, zip files are the most efficient and effective basic storage type for sketches in sourmash, and as of the branchwater plugin v0.9.0, they are fully supported!

You can create zipfiles with sourmash like so:
```
sourmash sig cat <list of sketches> -o sigs.zip
```

### Using manifests instead of zip files - why and when?

There are various places where we recommend using manifests instead of zip files. Why?

Well, first, if you are using a zip file created by sourmash, you are already using a manifest! And you will get all of the benefits described above!
 
But if you want to use a collection of multiple very large metagenomes (as search targets in `manysearch`, or as queries in `fastmultigather`), then standalone manifests might be a good solution for you.

This is for two specific reasons:
* first, metagenome sketches are often extremely large (100s of MBs to GBs), and it is not ideal to zip many large sketches into a single zip file;
* second, both `manysearch` and `fastmultigather` take a single argument that specifies collections of metagenomes which need to be loaded on demand, because they cannot fit into memory;

so the question becomes, how do you provide collections of large metagenomes to `manysearch` and `fastmultigather` in a single filename?

And the answer is: manifests. Manifests are a sourmash filetype that contains information about sketches without containing the actual sketch content, and they can be used as "catalogs" of sketch content.

The branchwater plugin supports manifest CSVs.  These can be created from lists of sketches by using `sourmash sig collect` or `sourmash sig manifest`; for example,
```
sourmash sig manifest <from file> -o manifest.csv
```
will create a manifest CSV from a list of sketches.

### Using RocksDB inverted indexes

The branchwater plugin also supports a database type that is not yet supported by sourmash: inverted indexes stored in a RocksDB database. These indexes provide fast and low-memory lookups when searching very large datasets, and are used for the branchwater petabase scale search hosted at [branchwater.sourmash.bio](https://branchwater.sourmash.bio). 

Some commands - `fastmultigather` and `manysearch` - support using these RocksDB-based inverted indexes. They can be created by running `sourmash scripts index`.

### Using "fromfiles"

**Note: We no longer recommend using "fromfiles". Use zip files or manifests instead.**

You can make a fromfile by listing a collection of .sig.gz files like so:
```
find /path/to/directory/ -name "*.sig.gz" -type f > directory.txt
```

When using a fromfile for search, we load all signatures into memory at the start in order to generate a manifest. To avoid memory issues, the signatures are not kept in memory, but instead re-loaded as described below for each command (see: Notes on concurrency and efficiency). This makes using fromfiles less efficient than `zip` files or manifests (as of v0.9.0).

## Running the commands

### Running `manysketch`

The `manysketch` command sketches one or more FASTA/FASTQ files into a zipped sourmash signature collection (`zip`). `manysketch` uses one thread per input file, so it can (very) efficiently sketch many files at once; and, because sequence file parsing is entirely implemented in Rust, it is much, _much_ faster than `sourmash sketch` for large FASTQ files.

To run `manysketch`, you need to build a text file list of FASTA/FASTQ files, with one on each line (`manysketch.csv`, below).  A simple way to do this for a directory is this command snippet:
```
echo name,genome_filename,protein_filename > manysketch.csv
for i in *.fa.gz
do
echo $i,$i,
done >> manysketch.csv
```

You can then run:

```
sourmash scripts manysketch manysketch.csv -o fa.zip
```
The output will be written to `fa.zip`

You can check if all signatures were written properly with
```
sourmash sig summarize fa.zip
```

To modify sketching parameters, use `--param-str` or `-p` and provide valid param string(s)
```
sourmash scripts manysketch fa.csv -o fa.zip -p k=21,k=31,k=51,scaled=1000,abund -p protein,k=10,scaled=200
```
See [the sourmash sketch docs](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-sketch-make-sourmash-signatures-from-sequence-data) for more information on param strings.

### Running `multisearch`

The `multisearch` command compares one or more query genomes, and one or more subject genomes. It differs from `manysearch` by loading all genomes into memory.

`multisearch` takes two input collections and outputs a CSV:
```
sourmash scripts multisearch query.sig.gz database.zip -o results.csv
```


The results file `results.csv`, will have 8 columns: `query` and `query_md5`, `match` and `match_md5`, and `containment`, `jaccard`, `max_containment`, and `intersect_hashes`.

The `pairwise` command does the same comparisons as `multisearch` but takes
only a single collection of sketches, for which it calculates all the pairwise comparisons. Since the comparisons are symmetric, it approximately
twice as fast as `multisearch`.

### Running `fastgather`

The `fastgather` command is a much faster version of `sourmash gather`.

`fastgather` takes a single query metagenome and a database, and outputs a CSV:
```
sourmash scripts fastgather query.sig.gz database.zip -o results.csv --cores 4
```

#### Using `fastgather` to create a picklist for `sourmash gather`

One handy use case for `fastgather` is to create a picklist that can be used by `sourmash gather`. This makes full use of the speed of `fastgather` while producing a complete set of `gather` outputs.

For example, if you run a complete `gather` against GTDB rs214,

```
sourmash scripts fastgather SRR606249.trim.sig.gz \
    gdtb-rs214-k21.zip -o SRR606249.fastgather.csv -k 21
```

The resulting CSV file can then be used as a picklist for `sourmash gather` like so:

```
sourmash gather SRR606249.trim.sig.gz gtdb-rs214-k21.zip \
    --picklist SRR606249.fastgather.csv:match_name:ident \
    -o SRR606249.gather.csv
```

#### Example of picklist usage

A complete example Snakefile implementing the above workflow is available [in the sourmash-slainte Snakefile](https://github.com/dib-lab/sourmash-slainte/blob/main/Snakefile).

### Running `fastmultigather`

`fastmultigather` takes a collection of query metagenomes and a collection of sketches as a database, and outputs many CSVs:
```
sourmash scripts fastmultigather queries.manifest.csv database.zip --cores 4
```
We suggest using a manifest CSV for the queries.

The main advantage that `fastmultigather` has over running `fastgather` on multiple queries is that you only load the database files once with `fastmultigather`, which can be a significant time savings for large databases!

#### Output files for `fastmultigather`

`fastmultigather` will output two CSV files for each query, a `prefetch` file containing all overlapping matches between that query and the database, and a `gather` file containing the minimum metagenome cover for that query in the database.

The prefetch CSV will be named `{signame}.prefetch.csv`, and the gather CSV will be named `{signame}.gather.csv`.  Here, `{signame}` is the name of your sourmash signature.

**Warning:** At the moment, if two different queries have the same `{signame}`, the CSVs for one of the queries will be overwritten by the other query. The behavior here is undefined in practice, because of multithreading: we don't know what queries will be executed when or files will be written first.

### Running `manysearch`

The `manysearch` command compares one or more collections of query sketches, and one or more collections of subject sketches. It is the core command we use for searching petabase-scale databases of metagenomes for contained genomes.

`manysearch` takes two collections as input, and outputs a CSV:
```
sourmash scripts manysearch queries.zip metagenomes.manifest.csv -o results.csv
```
We suggest using a manifest CSV for the metagenome collection.

The results file here, `query.x.gtdb-reps.csv`, will have 8 columns: `query` and `query_md5`, `match` and `match_md5`, and `containment`, `jaccard`, `max_containment`, and `intersect_hashes`.

## Notes on concurrency and efficiency

Each command does things slightly differently, with implications for CPU and disk load. You can measure threading efficiency with `/usr/bin/time -v` on Linux systems, and disk load by number of complaints received when running.

`manysketch` loads one sequence file from disk per thread and sketches it using all signature params simultaneously.

`manysearch` loads all the queries at the beginning, and then loads one database sketch from disk per thread. The compute-per-database-sketch is dominated by I/O. So your number of threads should be chosen with care for disk load. We typically limit it to `-c 32` for shared disks. We suggest using a manifest CSV file for the database sketches.

`multisearch` loads all the queries and database sketches once, at the beginning, and then uses multithreading to search across all matching sequences. For large databases it is extremely efficient at using all available cores. So 128 threads or more should work fine! Zipfiles and manifests should work well.

`pairwise` acts just like `multisearch`, but only loads one file (and then does all comparisons between all pairs within that file).

Like `multisearch` and `pairwise`, `fastgather` loads everything at the beginning, and then uses multithreading to search across all matching sequences. For large databases it is extremely efficient at using all available cores. So 128 threads or more should work fine! We suggest using zipfile or manifests for the database.

`fastmultigather` loads the entire database once, and then loads one query from disk per thread. The compute-per-query can be significant, though, so multithreading efficiency here is less dependent on I/O and the disk is less likely to be saturated with many threads. We suggest limiting threads to between 32 and 64 to decrease shared disk load.

## Appendix 1 - `index` to create a low-memory index

The command `sourmash scripts index` makes an on-disk inverted index
for low memory fast search. Indexing takes a while, but then search
takes fewer resources.

Currently only `fastmultigather` and `manysearch` can use this kind of index.

We suggest using the extension `.rocksdb` for these databases, as we
use [RocksDB](https://rocksdb.org/) for the underlying database storage
mechanism.
