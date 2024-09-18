# The branchwater plugin for sourmash

| command | functionality | docs |
| -------- | -------- | -------- |
| `manysketch` | Rapidly build sketches for many input files     | [link](#Running-manysketch)     |
| `fastgather` | Multithreaded `gather` of **one** metagenome against a database| [link](#Running-fastgather)
| `fastmultigather` | Multithreaded `gather` of **multiple** metagenomes against a database | [link](#Running-fastmultigather)
| `manysearch` | Multithreaded containment search for many queries in many large metagenomes | [link](#Running-manysearch)
| `multisearch` | Multithreaded comparison of multiple sketches, in memory | [link](#Running-multisearch-and-pairwise)
| `pairwise` | Multithreaded pairwise comparison of multiple sketches, in memory | [link](#Running-multisearch-and-pairwise)
| `cluster` | cluster sequences based on similarity data from `pairwise` or `multisearch` | [link](#Running-cluster)
| `index` | build a RocksDB inverted index for efficient containment queries | [link](#Running-index)
| `sigcat` | concatenate signatures to a zipfile | [link](#Running-sigcat) |

This repository implements multithreaded plugins for
[sourmash](https://sourmash.readthedocs.io/) that provide very fast
implementations of `sketch`, `search`, and `gather`. These commands
are typically hundreds to thousands of times faster, and 10-50x lower
memory, than the current sourmash code. For example, a `gather` of
SRR606249 with sourmash v4.8.6 against GTDB rs214 takes 40 minutes and
14 GB of RAM, while `fastgather` with 64 cores takes only 2 minutes
and 2 GB of RAM.

The main *drawback* to these plugin commands is that their inputs are
more restricted than the native sourmash commands, and in some cases
their outputs are in a different format. This means that your input
files may need to be prepared differently, and the output may need to
be processed differently.  The plugin commands are also a bit less
user friendly, because (for now) we're more focused on speed than
polish and user experience.

**Note:** As of v0.9.5, the outputs of `fastgather` and `fastmultigather` almost completely match the output of `sourmash gather`; see below for details.

## Input file formats

sourmash supports a variety of different storage formats for sketches (see [sourmash docs](https://sourmash.readthedocs.io/en/latest/command-line.html#choosing-signature-output-formats)), and the branchwater plugin works with some (but not all) of them. Branchwater _also_ supports an additional database type, a RocksDB-based inverted index, that is not (yet) supported natively by sourmash (through v4.8.11).

<!-- **As of v0.9.0, we recommend using zip files or manifest CSVs whenever you need to provide multiple sketches.** CTB -->

| command | command input | database format |
| -------- | -------- | -------- |
| `manysketch`     | CSV with input fasta/fastq paths (details below)    | _produces_ Zip database |
| `gather`     | Single metagenome in sig, zip, or pathlist     | Zip or pathlist |
| `fastmultigather` | Multiple metagenomes in sig, zip, or pathlist | Zip, pathlist, or rocksdb index |
| `manysearch` | Multiple genomes in sig, zip, or pathlist | Zip, pathlist, or rocksdb index |
| `multisearch` | Multiple sketches in sig, zip, or pathlist | Multiple sketches in sig, zip, or pathlist |
| `pairwise` | Multiple sketches in sig, zip, or pathlist | N/A |
| `cluster`| Output from `pairwise` or `multisearch`| N/A |
| `index` | Multiple sketches in sig, zip, or pathlist | N/A |

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

<!--
CTB

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
sourmash sig manifest <pathlist> -o manifest.csv
```
will create a manifest CSV from a list of sketches.
-->

### Using RocksDB inverted indexes

The branchwater plugin also supports a database type that is not yet supported by sourmash: inverted indexes stored in a RocksDB database. These indexes provide fast and low-memory lookups when searching very large datasets, and are used for the branchwater petabase scale search hosted at [branchwater.sourmash.bio](https://branchwater.sourmash.bio). 

Some commands - `fastmultigather` and `manysearch` - support using these RocksDB-based inverted indexes. They can be created by running `sourmash scripts index`. See [the `index` documentation, below](#Running-index).

### Using "pathlists"

<!-- **Note: We no longer recommend using "pathlists". Use zip files or manifests instead.** CTB -->

You can make a pathlist by listing a collection of .sig.gz files like so:
```
find /path/to/directory/ -name "*.sig.gz" -type f > directory.txt
```

When using a pathlist for search, we load all signatures into memory at the start in order to generate a manifest. To avoid memory issues, the signatures are not kept in memory, but instead re-loaded as described below for each command (see: Notes on concurrency and efficiency). This makes using pathlists less efficient than `zip` files (as of v0.9.0).


<!-- or manifests... CTB -->

## Running the commands

### Running `manysketch`

The `manysketch` command sketches a list of FASTA/FASTQ files into a
zipped sourmash signature collection (`zip`). `manysketch` uses one
thread per input file, so it can (very) efficiently sketch many files
at once; and, because sequence file parsing is entirely implemented in
Rust, it is much, _much_ faster than `sourmash sketch` for large FASTQ
files. However, it does not currently support translation,
i.e. protein signature generation from DNA FASTA.

`manysketch` can be used to sketch collections of genomes, pairs of R1/R2
files, and more flexible collections based on prefix - which approach is
used depends on the columns provided to `manysketch` in the input CSV file.

#### Specifying input FASTA

To run `manysketch`, you need to build a text file list of FASTA/FASTQ files (see `manysketch.csv` example, below).

The following formats are accepted:
- 3 columns: `name,genome_filename,protein_filename`
  >`genome_filename` entries are considered DNA FASTA, `protein_filename` entries are considered protein FASTA.
- 3 columns: `name,read1,read2`
  > All entries considered DNA FASTA, and both `read1` and `read2` files are used as input for a single sketch with name `name`.
- 4 columns: `name,input_moltype,prefix,exclude`
  > This filetype uses `glob` to find files that match `prefix` but do not match `exclude`. As such, `*` are ok in the `prefix` and `exclude` columns. Since we are dealing with "prefixes" here, we automatically search with `*` on the end of the `prefix` entry.

A simple way to build a manysketch input file for a directory is this command snippet:
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

#### singleton sketching

`manysketch` also supports building independent sketches for each record in a FASTA file (`--singleton`).

You can run:

```
sourmash scripts manysketch manysketch.csv -o fa.zip --singleton
```
The output will be written to `fa.zip`

You can check if all signatures were written properly with
```
sourmash sig summarize fa.zip
```
The number of sketches per parameter combination should equal the total number of records in all input FASTA.
The `name` column will not be used. Instead, each sketch will be named from the FASTA record name.

### Running `multisearch` and `pairwise`

The `multisearch` command compares one or more query genomes, and one or more subject genomes. It differs from `manysearch` by loading all genomes into memory.

`multisearch` takes two input collections and outputs a CSV:
```
sourmash scripts multisearch query.sig.gz database.zip -o results.csv
```


The results file `results.csv`, will have 8 columns: `query` and `query_md5`, `match` and `match_md5`, and `containment`, `jaccard`, `max_containment`, and `intersect_hashes`.

The `pairwise` command does the same comparisons as `multisearch` but takes
only a single collection of sketches, for which it calculates all the pairwise comparisons. Since the comparisons are symmetric, it is approximately
twice as fast as `multisearch`.

The `-t/--threshold` for `multisearch` and `pairwise` applies to the
containment of query-in-target and defaults to 0.01.  To report 
_any_ overlap between two sketches, set the threshold to 0.

### Running `fastgather`

The `fastgather` command is a much faster version of `sourmash gather`.

`fastgather` takes a single query metagenome and a database, and outputs a CSV:
```
sourmash scripts fastgather query.sig.gz database.zip -o results.csv --cores 4
```

As of v0.9.5, `fastgather` outputs the same columns as `sourmash gather`, with only a few exceptions:
* `match_name` is output instead of `name`;
* `match_md5` is output instead of `md5`;
* `match_filename` is output instead of `filename`, and the value is different;
* `potential_false_negative` is not present in `fastgather` output;

### Running `fastmultigather`

`fastmultigather` takes a collection of query metagenomes and a collection of sketches as a database, and outputs many CSVs:
```
sourmash scripts fastmultigather queries.manifest.csv database.zip --cores 4 --save-matches
```
<!-- We suggest using a manifest CSV for the queries. CTB -->

The main advantage that `fastmultigather` has over running `fastgather` on multiple queries is that `fastmultigather` only needs to load the database once for all queries, unlike with `fastgather`; this can be a significant time savings for large databases!

#### Output files for `fastmultigather`

On a database of sketches (but not on RocksDB indexes) `fastmultigather` will output two CSV files for each query, a `prefetch` file containing all overlapping matches between that query and the database, and a `gather` file containing the minimum metagenome cover for that query in the database.

The prefetch CSV will be named `{signame}.prefetch.csv`, and the gather CSV will be named `{signame}.gather.csv`.  Here, `{signame}` is the name of your sourmash signature.

`--save-matches` is an optional flag that will save the matched hashes for each query in a separate sourmash signature `{signame}.matches.sig`. This can be useful for debugging or for further analysis.

When searching against a RocksDB index, `fastmultigather` will output a single file containing all gather results, specified with `-o/--output`. No prefetch results will be output.

`fastmultigather` gather CSVs provide the same columns as `fastgather`, above.

**Warning:** At the moment, if two different queries have the same `{signame}`, the CSVs for one of the queries will be overwritten by the other query. The behavior here is undefined in practice, because of multithreading: we don't know what queries will be executed when or files will be written first.

### Running `manysearch`

The `manysearch` command compares one or more collections of query sketches, and one or more collections of subject sketches. It is the core command we use for searching petabase-scale databases of metagenomes for contained genomes.

`manysearch` takes two collections as input, and outputs a CSV:
```
sourmash scripts manysearch queries.zip metagenomes.manifest.csv -o results.csv
```
<!-- We suggest using a manifest CSV for the metagenome collection. -->

The results file here, `query.x.gtdb-reps.csv`, will have the
following columns: `query`, `query_md5`, `match_name`, `match_md5`,
`containment`, `jaccard`, `max_containment`, `intersect_hashes`,
`query_containment_ani`.

If you run `manysearch` _without_ using a RocksDB database (that is,
against regular sketches), the results file will also have the
following columns: , `match_containment_ani`,
`average_containment_ani`, and `max_containment_ani`.

Finally, if using sketches that have abundance information, the
results file will also contain the following columns: `average_abund`,
`median_abund`, `std_abund`, `n_weighted_found`, and `total_weighted_hashes`.

See
[the prefetch CSV output column documentation](https://sourmash.readthedocs.io/en/latest/classifying-signatures.html#appendix-e-prefetch-csv-output-columns)
for information on these various columns.

The `-t/--threshold` for `manysearch` applies to the containment of
query-in-database (e.g. genome in metagenome) and defaults to 0.01.
To report _any_ overlap between two sketches, set the threshold to 0.
(This will produce many, many results when searching a collection of
metagenomes!)

By default, `manysearch` will display the contents of the CSV file in a
human-readable format. This can be disabled with `-N/--no-pretty-print`
when executing large searches.

### Running `cluster`

The `cluster` command conducts graph-based clustering via the sequence similarity measures in `pairwise` or `multisearch` outputs. It is a new command and we are exploring its utility.

`cluster` takes the csv output of `pairwise` or `multisearch` input, and outputs two CSVs:

1. `-o`, `--output` will contain the names of the clusters and the `ident` of each sequence included in the cluster (e.g. `Component_1, name1;name2`)

```
cluster,nodes
Component_1,name1;name2;name3
Component_2,name4
```

2. `--cluster-sizes` will contain information on cluster size, with a counts for the number of clusters of that size. For the two clusters above, the counts would look like this:

```
cluster_size,count
3,1
1,1
```

`cluster` takes a `--similarity_column` argument to specify which of the similarity columns, with the following choices: `containment`, `max_containment`, `jaccard`, `average_containment_ani`, `maximum_containment_ani`. All values should be input as fractions (e.g. 0.9 for 90%)

### Running `index`

The `index` subcommand creates a RocksDB inverted index that can be
used as a database for `manysearch` (containment queries into
mixtures) and `fastmultigather` (mixture decomposition against a
database of genomes).

RocksDB inverted indexes support fast, low-latency, and low-memory
queries, in exchange for a time-intensive indexing step and extra disk
space. They are also used for the
[branchwater.sourmash.bio](https://branchwater.sourmash.bio/)
real-time SRA metagenome query.

Note that RocksDB indexes do not support abundance estimation for
containment queries.

To run `index`, provide it with multiple sketches in sig, zip, or
pathlist format, and specify the desired output directory; we suggest
using the `.rocksdb` extension for RocksDB databases, e.g. `-o
gtdb-rs214-k31.rocksdb`.

By default, as of v0.9.7, `index` will store a copy of the sketches
along with the inverted index.  This will substantially increase the
disk space required for large databases.  You can provide an optional
`--no-internal-storage` to `index` to store them externally, which
reduces the disk space needed for the index.  Read below for technical
details!

#### Internal vs external storage of sketches in a RocksDB index

(The below applies to v0.9.7 and later of the plugin; for v0.9.6 and
before, only external storage was implemented.)

RocksDB indexes support containment queries (a la the
[branchwater application](https://github.com/sourmash-bio/branchwater)),
as well as `gather`-style mixture decomposition (see
[Irber et al., 2022](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2)).
For this plugin, the `manysearch` command supports a RocksDB index for
the database for containment queries, and `fastmultigather` can use a
RocksDB index for the database of genomes.

RocksDB indexes contain references to the sketches used to construct
the index. If `--internal-storage` is set (which is the default), a
copy of the sketches is stored within the RocksDB database directory;
if `--no-internal-storage` is provided, then the references point to
the original source sketches used to construct the database, wherever
they reside on your disk.

The sketches *are not used* by `manysearch`, but *are used* by
`fastmultigather`: with v0.9.6 and later, you'll get an error if you
run `fastmultigather` against a RocksDB index where the sketches
cannot be loaded.

What this means is therefore a bit complicated, but boils down to
the following two approaches:

1. The safest thing to do is build a RocksDB index and use internal
storage (the default). This will consume more disk space but your
RocksDB database will always be usable for both `manysearch` and
`fastmultigather`, as well as the branchwater app.
2. If you want to avoid storing duplicate copies of your sketches,
then specify `--no-internal-storage` and provide a stable absolute
path to the source sketches.  This will again support both
`manysearch` and `fastmultigather`, as well as the branchwater app.
If the source sketches later become unavailable, `fastmultigather`
will stop working (although `manysearch` and the branchwater app
should be fine).

You should (for the moment) avoid specifying relative paths to the
sketches when running `index`.  Follow
[sourmash_branchwater_plugin#415](https://github.com/sourmash-bio/sourmash_plugin_branchwater/issues/415)
if better support for relative paths is of interest!

#### Links and more materials

Note that RocksDB indexes are implemented in the core
[sourmash Rust library](https://crates.io/crates/sourmash), and
used in downstream software packages (this plugin, and
[the branchwater application code](https://github.com/sourmash-bio/branchwater)).  The above documentation applies to sourmash core v0.15.0.

### Running `sigcat`

The `sigcat` command combines signatures into a single sourmash Zipfile. It's equivalent `sourmash` command is `sourmash sig cat`.

For example:
```
sourmash scripts sigcat file1.sig file2.sig -o all.zip
```
will combine all signatures in `file1.sig` and `file2.sig` and put them in the file `all.zip`.

`sigcat` can be used to select out specific signatures as well, by `ksize`, `moltype`, and `scaled`, where signatures will be downsampled to the desired scaled.

## Notes on concurrency and efficiency

Each command does things somewhat differently, with implications for
CPU and disk load; moreover, the performance will vary depending on
the database format. You can measure threading efficiency with
`/usr/bin/time -v` on Linux systems, and disk load by number of
complaints received when running.  Your best bet is to
[just ask the team for suggestions](https://github.com/dib-lab/sourmash/issues),
but you can lightly skim the docs below and play with some small data
sets, too!

---

`manysketch` loads one sequence file from disk per thread and sketches
it using all signature params simultaneously.

`manysearch` loads all the queries at the beginning, and then loads
one database sketch from disk per thread. The
compute-per-database-sketch is dominated by I/O. So your number of
threads should be chosen with care for disk load. We typically limit
it to `-c 32` for shared disks. <!-- We suggest using a manifest CSV
file for the database sketches. CTB -->

`multisearch` loads all the queries and database sketches once, at the
beginning, and then uses multithreading to search across all matching
sequences. For large databases it is extremely efficient at using all
available cores. So 128 threads or more should work fine!
Zipfiles <!-- and manifests CTB --> should work well.

`pairwise` acts just like `multisearch`, but only loads one file (and
then does all comparisons between all pairs within that file).

Like `multisearch` and `pairwise`, `fastgather` loads everything at
the beginning, and then uses multithreading to search across all
matching sequences. For large databases it is extremely efficient at
using all available cores. So 128 threads or more should work fine! We
suggest using zipfiles <!-- or manifests CTB --> for the database.

`fastmultigather` loads the entire database once, and then loads one
query from disk per thread. The compute-per-query can be significant,
though, so multithreading efficiency here is less dependent on I/O and
the disk is less likely to be saturated with many threads. We suggest
limiting threads to between 32 and 64 to decrease shared disk load.

`cluster` loads the entire file multithreaded, and then populates the
graph sequentially.

## Appendix 1 - `index` to create a low-memory index

The command `sourmash scripts index` makes an on-disk inverted index
for low memory fast search. Indexing takes a while, but then search
takes fewer resources.

Currently only `fastmultigather` and `manysearch` can use this kind of index.

`fastmultigather` with this index produces a complete set of `sourmash gather` columns.

We suggest using the extension `.rocksdb` for these databases, as we
use [RocksDB](https://rocksdb.org/) for the underlying database storage
mechanism.
