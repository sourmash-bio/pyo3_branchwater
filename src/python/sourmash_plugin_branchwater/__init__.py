#! /usr/bin/env python
import sys
import argparse
from sourmash.plugins import CommandLinePlugin
from sourmash.logging import notify
import os
import importlib.metadata

from . import sourmash_plugin_branchwater
from . import prettyprint

__version__ = importlib.metadata.version("sourmash_plugin_branchwater")


def print_version():
    notify(
        f"=> sourmash_plugin_branchwater {__version__}; cite Irber et al., doi: 10.1101/2022.11.02.514947\n"
    )


def get_max_cores():
    try:
        if "SLURM_CPUS_ON_NODE" in os.environ:
            return int(os.environ["SLURM_CPUS_ON_NODE"])
        elif "SLURM_JOB_CPUS_PER_NODE" in os.environ:
            cpus_per_node_str = os.environ["SLURM_JOB_CPUS_PER_NODE"]
            return int(cpus_per_node_str.split("x")[0])
        else:
            return os.cpu_count()
    except Exception:
        return os.cpu_count()


def set_thread_pool(user_cores):
    avail_threads = get_max_cores()
    num_threads = min(avail_threads, user_cores) if user_cores else avail_threads
    if user_cores and user_cores > avail_threads:
        notify(
            f"warning: only {avail_threads} threads available, using {avail_threads}"
        )
    actual_rayon_cores = sourmash_plugin_branchwater.set_global_thread_pool(num_threads)
    return actual_rayon_cores


class Branchwater_Manysearch(CommandLinePlugin):
    command = "manysearch"
    description = "search many metagenomes for contained genomes"

    def __init__(self, p):
        super().__init__(p)
        p.add_argument("query_paths", help="input file of sketches")
        p.add_argument("against_paths", help="input file of sketches")
        p.add_argument(
            "-o", "--output", required=True, help="CSV output file for matches"
        )
        p.add_argument(
            "-t",
            "--threshold",
            default=0.01,
            type=float,
            help="containment threshold for reporting matches (default: 0.01)",
        )
        p.add_argument(
            "-k",
            "--ksize",
            default=31,
            type=int,
            help="k-mer size at which to select sketches",
        )
        p.add_argument(
            "-s",
            "--scaled",
            default=None,
            type=int,
            help="scaled factor at which to do comparisons",
        )
        p.add_argument(
            "-m",
            "--moltype",
            default="DNA",
            choices=["DNA", "protein", "dayhoff", "hp", "skipm1n3", "skipm2n3"],
            help="molecule type: DNA, protein, dayhoff, hp, or skipmer (skipm1n3 or skipm2n3); default DNA",
        )
        p.add_argument(
            "-c",
            "--cores",
            default=0,
            type=int,
            help="number of cores to use (default is all available)",
        )
        p.add_argument(
            "-P",
            "--pretty-print",
            action="store_true",
            default=True,
            help="display results after search finishes (default: True)",
        )
        p.add_argument(
            "-N",
            "--no-pretty-print",
            action="store_false",
            dest="pretty_print",
            help="do not display results (e.g. for large output)",
        )
        p.add_argument(
            "--ignore-abundance",
            action="store_true",
            help="do not do expensive abundance calculations",
        )
        p.add_argument(
            "-A",
            "--output-all-comparisons",
            action="store_true",
            help="ignore threshold and output all comparisons; against a RocksDB database, this will only output comparisons with some overlap",
        )

    def main(self, args):
        print_version()
        notify(
            f"ksize: {args.ksize} / scaled: {args.scaled} / moltype: {args.moltype} / threshold: {args.threshold}"
        )
        num_threads = set_thread_pool(args.cores)

        notify(
            f"searching all sketches in '{args.query_paths}' against '{args.against_paths}' using {num_threads} threads"
        )

        super().main(args)
        status = sourmash_plugin_branchwater.do_manysearch(
            args.query_paths,
            args.against_paths,
            args.threshold,
            args.ksize,
            args.scaled,
            args.moltype,
            args.output,
            args.ignore_abundance,
            args.output_all_comparisons,
        )
        if status == 0:
            notify(f"...manysearch is done! results in '{args.output}'")

            if args.pretty_print:
                prettyprint.pretty_print_manysearch(args.output)
        return status


class Branchwater_Fastgather(CommandLinePlugin):
    command = "fastgather"
    description = "massively parallel sketch gather"

    def __init__(self, p):
        super().__init__(p)
        p.add_argument("query_sig", help="metagenome sketch")
        p.add_argument("against_paths", help="input file of sketches")
        p.add_argument(
            "-o",
            "--output-gather",
            required=True,
            help="save gather output (minimum metagenome cover) to this file",
        )
        p.add_argument(
            "--output-prefetch", help="save prefetch output (all overlaps) to this file"
        )
        p.add_argument(
            "-t",
            "--threshold-bp",
            default=50000,
            type=float,
            help="threshold in estimated base pairs, for reporting matches (default: 50kb)",
        )
        p.add_argument(
            "-k",
            "--ksize",
            default=31,
            type=int,
            help="k-mer size at which to do comparisons (default: 31)",
        )
        p.add_argument(
            "-s",
            "--scaled",
            default=None,
            type=int,
            help="scaled factor at which to do comparisons (default: determined from query)",
        )
        p.add_argument(
            "-m",
            "--moltype",
            default="DNA",
            choices=["DNA", "protein", "dayhoff", "hp", "skipm1n3", "skipm2n3"],
            help="molecule type: DNA, protein, dayhoff, hp, or skipmer (skipm1n3 or skipm2n3); default DNA",
        )
        p.add_argument(
            "-c",
            "--cores",
            default=0,
            type=int,
            help="number of cores to use (default is all available)",
        )

    def main(self, args):
        print_version()
        notify(
            f"ksize: {args.ksize} / scaled: {args.scaled} / moltype: {args.moltype} / threshold bp: {args.threshold_bp}"
        )

        num_threads = set_thread_pool(args.cores)

        notify(
            f"gathering all sketches in '{args.query_sig}' against '{args.against_paths}' using {num_threads} threads"
        )
        super().main(args)
        status = sourmash_plugin_branchwater.do_fastgather(
            args.query_sig,
            args.against_paths,
            int(args.threshold_bp),
            args.ksize,
            args.scaled,
            args.moltype,
            args.output_gather,
            args.output_prefetch,
        )
        if status == 0:
            notify(f"...fastgather is done! gather results in '{args.output_gather}'")
            if args.output_prefetch:
                notify(f"prefetch results in '{args.output_prefetch}'")
        return status


class Branchwater_Fastmultigather(CommandLinePlugin):
    command = "fastmultigather"
    description = "massively parallel sketch multigather"

    def __init__(self, p):
        super().__init__(p)
        p.add_argument("query_paths", help="input file of sketches to query")
        p.add_argument(
            "against_paths",
            help="input file of sketches to search against \
                       OR a branchwater indexed database generated with 'sourmash scripts index'",
        )
        p.add_argument(
            "-t",
            "--threshold-bp",
            default=50000,
            type=float,
            help="threshold in estimated base pairs, for reporting matches (default: 50kb)",
        )
        p.add_argument(
            "-k",
            "--ksize",
            default=31,
            type=int,
            help="k-mer size at which to do comparisons (default: 31)",
        )
        p.add_argument(
            "-s",
            "--scaled",
            default=None,
            type=int,
            help="scaled factor at which to do comparisons (default: determined from query collection)",
        )
        p.add_argument(
            "-m",
            "--moltype",
            default="DNA",
            choices=["DNA", "protein", "dayhoff", "hp", "skipm1n3", "skipm2n3"],
            help="molecule type: DNA, protein, dayhoff, hp, or skipmer (skipm1n3 or skipm2n3); default DNA",
        )
        p.add_argument(
            "-c",
            "--cores",
            default=0,
            type=int,
            help="number of cores to use (default is all available)",
        )
        p.add_argument(
            "-o",
            "--output",
            help="CSV output file containing gather matches",
        )
        p.add_argument(
            "--create-empty-results",
            "--create-empty-prefetch-results",
            action="store_true",
            default=False,
            help="create empty prefetch results file for each query, even if no matches (non-RockSDB only)",
        )
        p.add_argument(
            "--save-matches",
            action="store_true",
            default=False,
            help="save matched hashes for every input to a signature",
        )

    def main(self, args):
        print_version()
        notify(
            f"ksize: {args.ksize} / scaled: {args.scaled} / moltype: {args.moltype} / threshold bp: {args.threshold_bp} / save matches: {args.save_matches}"
        )

        num_threads = set_thread_pool(args.cores)

        notify(
            f"gathering all sketches in '{args.query_paths}' against '{args.against_paths}' using {num_threads} threads"
        )
        super().main(args)
        status = sourmash_plugin_branchwater.do_fastmultigather(
            args.query_paths,
            args.against_paths,
            int(args.threshold_bp),
            args.ksize,
            args.scaled,
            args.moltype,
            args.output,
            args.save_matches,
            args.create_empty_results,
        )
        if status == 0:
            notify(f"...fastmultigather is done!")
        return status


class Branchwater_Index(CommandLinePlugin):
    command = "index"
    description = "Build Branchwater RevIndex"

    def __init__(self, p):
        super().__init__(p)
        p.add_argument("siglist", help="input file of sketches")
        p.add_argument(
            "-o", "--output", required=True, help="output file for the index"
        )
        p.add_argument(
            "-k",
            "--ksize",
            default=31,
            type=int,
            help="k-mer size at which to select sketches",
        )
        p.add_argument(
            "-s",
            "--scaled",
            default=None,
            type=int,
            help="scaled factor at which to select sketches (default: chosen as max from collection)",
        )
        p.add_argument(
            "-m",
            "--moltype",
            default="DNA",
            choices=["DNA", "protein", "dayhoff", "hp", "skipm1n3", "skipm2n3"],
            help="molecule type: DNA, protein, dayhoff, hp, or skipmer (skipm1n3 or skipm2n3); default DNA",
        )
        p.add_argument(
            "-c",
            "--cores",
            default=0,
            type=int,
            help="number of cores to use (default is all available)",
        )
        p.add_argument(
            "--internal-storage",
            default=True,
            action="store_true",
            help="build indexes that contain sketches and are relocatable (default: True)",
        )
        p.add_argument(
            "--no-internal-storage",
            "--no-store-sketches",
            action="store_false",
            help="do not store sketches in the index; index may not be relocatable (default: False)",
            dest="internal_storage",
        )

    def main(self, args):
        notify(
            f"ksize: {args.ksize} / scaled: {args.scaled} / moltype: {args.moltype} "
        )

        num_threads = set_thread_pool(args.cores)

        notify(f"indexing all sketches in '{args.siglist}'")

        super().main(args)
        status = sourmash_plugin_branchwater.do_index(
            args.siglist,
            args.ksize,
            args.scaled,
            args.moltype,
            args.output,
            False,  # colors - currently must be false?
            args.internal_storage,
        )
        if status == 0:
            notify(f"...index is done! results in '{args.output}'")
        return status


class Branchwater_Check(CommandLinePlugin):
    command = "check"
    description = "Check Branchwater RevIndex"

    def __init__(self, p):
        super().__init__(p)
        p.add_argument("index", help="RocksDB index file created with 'index'")
        p.add_argument("--quick", action="store_true")
        p.add_argument(
            "--writable",
            "--upgrade",
            action="store_true",
            help="open database in read-write mode to upgrade the internal format if needed",
        )

    def main(self, args):
        notify(f"checking index '{args.index}'")
        super().main(args)
        status = sourmash_plugin_branchwater.do_check(
            args.index, args.quick, args.writable
        )
        if status == 0:
            notify(f"...index is ok!")
        return status


class Branchwater_Multisearch(CommandLinePlugin):
    command = "multisearch"
    description = "massively parallel in-memory sketch search"

    def __init__(self, p):
        super().__init__(p)
        p.add_argument("query_paths", help="input file of sketches")
        p.add_argument("against_paths", help="input file of sketches")
        p.add_argument(
            "-o", "--output", required=True, help="CSV output file for matches"
        )
        p.add_argument(
            "-t",
            "--threshold",
            default=0.01,
            type=float,
            help="containment threshold for reporting matches (default: 0.01)",
        )
        p.add_argument(
            "-k",
            "--ksize",
            default=31,
            type=int,
            help="k-mer size at which to select sketches",
        )
        p.add_argument(
            "-s",
            "--scaled",
            default=None,
            type=int,
            help="scaled factor at which to do comparisons (default: determined from query collection)",
        )
        p.add_argument(
            "-m",
            "--moltype",
            default="DNA",
            choices=["DNA", "protein", "dayhoff", "hp", "skipm1n3", "skipm2n3"],
            help="molecule type: DNA, protein, dayhoff, hp, or skipmer (skipm1n3 or skipm2n3); default DNA",
        )
        p.add_argument(
            "-c",
            "--cores",
            default=0,
            type=int,
            help="number of cores to use (default is all available)",
        )
        p.add_argument(
            "-a", "--ani", action="store_true", help="estimate ANI from containment"
        )
        p.add_argument(
            "-p",
            "--prob-significant-overlap",
            action="store_true",
            help="estimate probability of overlap for significance ranking of search results, of the specific query and match, given all queries and possible matches",
        )
        p.add_argument(
            "-A",
            "--output-all-comparisons",
            action="store_true",
            help="ignore threshold and output all comparisons",
        )

    def main(self, args):
        print_version()
        notify(
            f"ksize: {args.ksize} / scaled: {args.scaled} / moltype: {args.moltype} / threshold: {args.threshold}"
        )

        num_threads = set_thread_pool(args.cores)

        notify(
            f"searching all sketches in '{args.query_paths}' against '{args.against_paths}' using {num_threads} threads"
        )
        notify(
            f"estimate ani? {args.ani} / estimate probability of overlap? {args.prob_significant_overlap}"
        )

        super().main(args)
        status = sourmash_plugin_branchwater.do_multisearch(
            args.query_paths,
            args.against_paths,
            args.threshold,
            args.ksize,
            args.scaled,
            args.moltype,
            args.ani,
            args.prob_significant_overlap,
            args.output_all_comparisons,
            args.output,
        )
        if status == 0:
            notify(f"...multisearch is done! results in '{args.output}'")
        return status


class Branchwater_Pairwise(CommandLinePlugin):
    command = "pairwise"
    description = "massively parallel in-memory pairwise comparisons"

    def __init__(self, p):
        super().__init__(p)
        p.add_argument("sig_paths", help="input file of sketches")
        p.add_argument(
            "-o", "--output", required=True, help="CSV output file for matches"
        )
        p.add_argument(
            "-t",
            "--threshold",
            default=0.01,
            type=float,
            help="containment threshold for reporting matches",
        )
        p.add_argument(
            "-k",
            "--ksize",
            default=31,
            type=int,
            help="k-mer size at which to select sketches",
        )
        p.add_argument(
            "-s",
            "--scaled",
            default=None,
            type=int,
            help="scaled factor at which to do comparisons",
        )
        p.add_argument(
            "-m",
            "--moltype",
            default="DNA",
            choices=["DNA", "protein", "dayhoff", "hp", "skipm1n3", "skipm2n3"],
            help="molecule type: DNA, protein, dayhoff, hp, or skipmer (skipm1n3 or skipm2n3); default DNA",
        )
        p.add_argument(
            "-c",
            "--cores",
            default=0,
            type=int,
            help="number of cores to use (default is all available)",
        )
        p.add_argument(
            "-a", "--ani", action="store_true", help="estimate ANI from containment"
        )
        p.add_argument(
            "--write-all",
            "--write-self-comparisons",
            action="store_true",
            help="write self comparisons for all sketches",
        )
        p.add_argument(
            "-A",
            "--output-all-comparisons",
            action="store_true",
            help="ignore threshold and output all comparisons",
        )

    def main(self, args):
        print_version()
        notify(
            f"ksize: {args.ksize} / scaled: {args.scaled} / moltype: {args.moltype} / threshold: {args.threshold}"
        )

        num_threads = set_thread_pool(args.cores)

        notify(
            f"pairwise-comparing all sketches in '{args.sig_paths}' using {num_threads} threads"
        )

        super().main(args)
        status = sourmash_plugin_branchwater.do_pairwise(
            args.sig_paths,
            args.threshold,
            args.ksize,
            args.scaled,
            args.moltype,
            args.ani,
            args.write_all,
            args.output_all_comparisons,
            args.output,
        )
        if status == 0:
            notify(f"...pairwise is done! results in '{args.output}'")
        return status


class Branchwater_SingleSketch(CommandLinePlugin):
    command = "singlesketch"
    description = "sketch a single sample"

    def __init__(self, p):
        super().__init__(p)
        p.add_argument(
            "input_filenames", help="input file(s); use '-' for stdin", nargs="+"
        )
        p.add_argument(
            "-o",
            "--output",
            required=True,
            help="output file for the signature or - for stdout",
        )
        p.add_argument(
            "-p",
            "--param-string",
            action="append",
            type=str,
            default=[],
            help="parameter string for sketching (default: k=31,scaled=1000)",
        )
        p.add_argument(
            "-n",
            "--name",
            help="optional name for the signature, default is the basename of input path",
        )
        p.add_argument(
            "-I",
            "--input-moltype",
            "--input-molecule-type",
            choices=["DNA", "dna", "protein"],
            default="DNA",
            help="molecule type of input sequence (DNA or protein)",
        )

    def main(self, args):
        print_version()
        if not args.param_string:
            args.param_string = ["k=31,scaled=1000,dna"]

        # Check and append 'dna' if no moltype is found in a param string
        moltypes = ["dna", "protein", "dayhoff", "hp", "skipm1n3", "skipm2n3"]
        updated_param_strings = []

        for param in args.param_string:
            # If no moltype is found in the current param string, append 'dna'
            if not any(mt in param for mt in moltypes):
                param += ",dna"
            updated_param_strings.append(param)

        notify(f"params: {updated_param_strings}")

        # Convert the list of param strings to a single string
        args.param_string = "_".join(updated_param_strings).lower()

        # If --name is not provided, default to input_filename, but if the source file is -, set name to empty string
        signature_name = (
            args.name
            if args.name
            else (
                os.path.basename(args.input_filenames[0])
                if args.input_filenames[0] != "-"
                else ""
            )
        )

        notify(
            f"sketching {len(args.input_filenames)} files ({args.input_moltype}) with params '{args.param_string}' and name '{signature_name}' using a single thread"
        )

        super().main(args)
        status = sourmash_plugin_branchwater.do_singlesketch(
            args.input_filenames,
            args.input_moltype,
            args.param_string,
            args.output,
            signature_name,
        )  # Pass the name to Rust
        if status == 0:
            notify(f"...singlesketch is done! results in '{args.output}'")
        return status


class Branchwater_Manysketch(CommandLinePlugin):
    command = "manysketch"
    description = "massively parallel sketching"

    def __init__(self, p):
        super().__init__(p)
        p.add_argument(
            "fromfile_csv",
            help="a csv file containing paths to FASTA files. \
                        Columns must be: 'name,genome_filename,protein_filename' or 'name,read1,read2'",
        )
        p.add_argument(
            "-o", "--output", required=True, help="output zip file for the signatures"
        )
        p.add_argument(
            "-p",
            "--param-string",
            action="append",
            type=str,
            default=[],
            help="parameter string for sketching (default: k=31,scaled=1000)",
        )
        p.add_argument(
            "-c",
            "--cores",
            default=0,
            type=int,
            help="number of cores to use (default is all available)",
        )
        p.add_argument(
            "-s",
            "--singleton",
            action="store_true",
            help="build one sketch per FASTA record, i.e. multiple sketches per FASTA file",
        )
        p.add_argument(
            "-f",
            "--force",
            action="store_true",
            help="allow use of individual FASTA files in more than more sketch",
        )

    def main(self, args):
        print_version()
        if not args.param_string:
            args.param_string = ["dna,k=31,scaled=1000"]
        notify(f"params: {args.param_string}")

        # convert to a single string for easier rust handling
        args.param_string = "_".join(args.param_string)
        # lowercase the param string
        args.param_string = args.param_string.lower()

        num_threads = set_thread_pool(args.cores)

        notify(
            f"sketching all files in '{args.fromfile_csv}' using {num_threads} threads"
        )

        super().main(args)
        status = sourmash_plugin_branchwater.do_manysketch(
            args.fromfile_csv,
            args.param_string,
            args.output,
            args.singleton,
            args.force,
        )
        if status == 0:
            notify(f"...manysketch is done! results in '{args.output}'")
        return status


class Branchwater_Cluster(CommandLinePlugin):
    command = "cluster"
    description = 'cluster from "pairwise" or "multisearch" results'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument(
            "pairwise_csv",
            help="a csv file containing similarity information. \
                        Currently, only a branchwater 'pairwise' or 'multisearch' file will work",
        )
        p.add_argument(
            "-o", "--output", required=True, help="output csv file for the clusters"
        )
        p.add_argument(
            "--cluster-sizes",
            default=None,
            help="output file for the cluster size histogram",
        )
        p.add_argument(
            "--similarity-column",
            type=str,
            default="average_containment_ani",
            choices=[
                "containment",
                "max_containment",
                "jaccard",
                "average_containment_ani",
                "max_containment_ani",
            ],
            help="column to use as similarity measure",
        )
        p.add_argument(
            "-t",
            "--threshold",
            type=float,
            default=0.95,
            help="similarity threshold for clustering. Default: 95%% ANI (0.95)",
        )
        p.add_argument(
            "-c",
            "--cores",
            default=0,
            type=int,
            help="number of cores to use (default is all available)",
        )

    def main(self, args):
        print_version()

        num_threads = set_thread_pool(args.cores)

        notify(
            f"generating clusters for comparisons in '{args.pairwise_csv}' using {num_threads} threads"
        )

        super().main(args)
        status = sourmash_plugin_branchwater.do_cluster(
            args.pairwise_csv,
            args.output,
            args.similarity_column,
            args.threshold,
            args.cluster_sizes,
        )
        if status == 0:
            notify(f"...clustering is done! results in '{args.output}'")
            notify(f"                       cluster counts in '{args.cluster_sizes}'")
        return status


class Branchwater_Fastagather(CommandLinePlugin):
    command = "fastagather"
    description = "massively parallel gather directly from FASTA"

    def __init__(self, p):
        super().__init__(p)
        p.add_argument("query_fa", help="FASTA file")
        p.add_argument("against_paths", help="database file of sketches")
        p.add_argument(
            "-I",
            "--input-moltype",
            "--input-molecule-type",
            choices=["DNA", "dna", "protein"],
            default="DNA",
            help="molecule type of input sequence (DNA or protein)",
        )
        p.add_argument(
            "-o",
            "--output-gather",
            required=True,
            help="save gather output (minimum metagenome cover) to this file",
        )
        p.add_argument(
            "-t",
            "--threshold-hashes",
            default=4,
            type=float,
            help="threshold in hashes, for reporting matches (default: 4)",
        )
        p.add_argument(
            "-k",
            "--ksize",
            default=31,
            type=int,
            help="k-mer size at which to do comparisons (default: 31)",
        )
        p.add_argument(
            "-s",
            "--scaled",
            default=1000,
            type=int,
            help="scaled factor at which to do comparisons (default: 10)",
        )
        p.add_argument(
            "-m",
            "--moltype",
            default="DNA",
            choices=["DNA", "protein", "dayhoff", "hp"],
            help="molecule type for search (DNA, protein, dayhoff, hp, skipm1n3, or skipm2n3; default DNA)",
        )
        p.add_argument(
            "-c",
            "--cores",
            default=0,
            type=int,
            help="number of cores to use (default is all available)",
        )

    def main(self, args):
        print_version()
        notify(
            f"ksize: {args.ksize} / scaled: {args.scaled} / moltype: {args.moltype} / threshold hashes: {args.threshold_hashes}"
        )

        num_threads = set_thread_pool(args.cores)

        notify(
            f"gathering all sketches in '{args.query_fa}' against '{args.against_paths}' using {num_threads} threads"
        )
        super().main(args)
        status = sourmash_plugin_branchwater.do_fastagather(
            args.query_fa,
            args.against_paths,
            args.input_moltype,
            int(args.threshold_hashes),
            args.ksize,
            args.scaled,
            args.moltype,
            args.output_gather,
        )
        if status == 0:
            notify(f"...fastgather is done! gather results in '{args.output_gather}'")
        return status
