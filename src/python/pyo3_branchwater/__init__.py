#! /usr/bin/env python
import sys
import argparse
from sourmash.plugins import CommandLinePlugin
from sourmash.logging import notify
import os

from . import pyo3_branchwater


def get_max_cores():
    try:
        if 'SLURM_CPUS_ON_NODE' in os.environ:
            return int(os.environ['SLURM_CPUS_ON_NODE'])
        elif 'SLURM_JOB_CPUS_PER_NODE' in os.environ:
            cpus_per_node_str = os.environ['SLURM_JOB_CPUS_PER_NODE']
            return int(cpus_per_node_str.split('x')[0])
        else:
            return os.cpu_count()
    except Exception:
        return os.cpu_count()


def set_thread_pool(user_cores):
    avail_threads = get_max_cores()
    num_threads = min(avail_threads, user_cores) if user_cores else avail_threads
    if user_cores and user_cores > avail_threads:
        notify(f"warning: only {avail_threads} threads available, using {avail_threads}")
    actual_rayon_cores = pyo3_branchwater.set_global_thread_pool(num_threads)
    return actual_rayon_cores


class Branchwater_Manysearch(CommandLinePlugin):
    command = 'manysearch'
    description = 'massively parallel sketch search'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('query_paths',
                       help="a text file containing paths to .sig/.sig.gz files")
        p.add_argument('against_paths',
                       help="a text file containing paths to .sig/.sig.gz files")
        p.add_argument('-o', '--output', required=True,
                       help='CSV output file for matches')
        p.add_argument('-t', '--threshold', default=0.01, type=float,
                       help='containment threshold for reporting matches')
        p.add_argument('-k', '--ksize', default=31, type=int,
                       help='k-mer size at which to select sketches')
        p.add_argument('-s', '--scaled', default=1000, type=int,
                       help='scaled factor at which to do comparisons')
        p.add_argument('-c', '--cores', default=0, type=int,
                       help='number of cores to use (default is all available)')

    def main(self, args):
        notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold: {args.threshold}")

        num_threads = set_thread_pool(args.cores)

        notify(f"searching all sketches in '{args.query_paths}' against '{args.against_paths}' using {num_threads} threads")

        super().main(args)
        status = pyo3_branchwater.do_manysearch(args.query_paths,
                                                args.against_paths,
                                                args.threshold,
                                                args.ksize,
                                                args.scaled,
                                                args.output)
        if status == 0:
            notify(f"...manysearch is done! results in '{args.output}'")
        return status


class Branchwater_Fastgather(CommandLinePlugin):
    command = 'fastgather'
    description = 'massively parallel sketch gather'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('query_sig', help="metagenome sketch")
        p.add_argument('against_paths', help="a text file containing paths to .sig/.sig.gz files")
        p.add_argument('-o', '--output-gather', required=True,
                       help="save gather output (minimum metagenome cover) to this file")
        p.add_argument('--output-prefetch',
                       help="save prefetch output (all overlaps) to this file")
        p.add_argument('-t', '--threshold-bp', default=50000, type=float,
                       help='threshold in estimated base pairs, for reporting matches (default: 50kb)')
        p.add_argument('-k', '--ksize', default=31, type=int,
                       help='k-mer size at which to do comparisons (default: 31)')
        p.add_argument('-s', '--scaled', default=1000, type=int,
                       help='scaled factor at which to do comparisons (default: 1000)')
        p.add_argument('-c', '--cores', default=0, type=int,
                        help='number of cores to use (default is all available)')


    def main(self, args):
        notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold bp: {args.threshold_bp}")

        num_threads = set_thread_pool(args.cores)


        notify(f"gathering all sketches in '{args.query_sig}' against '{args.against_paths}' using {num_threads} threads")
        super().main(args)
        status = pyo3_branchwater.do_countergather(args.query_sig,
                                                   args.against_paths,
                                                   int(args.threshold_bp),
                                                   args.ksize,
                                                   args.scaled,
                                                   args.output_gather,
                                                   args.output_prefetch)
        if status == 0:
            notify(f"...fastgather is done! gather results in '{args.output_gather}'")
            if args.output_prefetch:
                notify(f"prefetch results in '{args.output_prefetch}'")
        return status


class Branchwater_Fastmultigather(CommandLinePlugin):
    command = 'fastmultigather'
    description = 'massively parallel sketch multigather'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('query_paths', help="a text file containing paths to .sig/.sig.gz files to query")
        p.add_argument('against_paths', help="a text file containing paths to .sig/.sig.gz files to search against")
        p.add_argument('-t', '--threshold-bp', default=50000, type=float,
                       help='threshold in estimated base pairs, for reporting matches (default: 50kb)')
        p.add_argument('-k', '--ksize', default=31, type=int,
                       help='k-mer size at which to do comparisons (default: 31)')
        p.add_argument('-s', '--scaled', default=1000, type=int,
                       help='scaled factor at which to do comparisons (default: 1000)')
        p.add_argument('-c', '--cores', default=0, type=int,
                       help='number of cores to use (default is all available)')


    def main(self, args):
        notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold bp: {args.threshold_bp}")

        num_threads = set_thread_pool(args.cores)

        notify(f"gathering all sketches in '{args.query_paths}' against '{args.against_paths}' using {num_threads} threads")
        super().main(args)
        status = pyo3_branchwater.do_multigather(args.query_paths,
                                                 args.against_paths,
                                                 int(args.threshold_bp),
                                                 args.ksize,
                                                 args.scaled)
        if status == 0:
            notify(f"...fastmultigather is done!")
        return status


class Branchwater_Index(CommandLinePlugin):
    command = 'index'
    description = 'Build Branchwater RevIndex'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('siglist',
                       help="a text file containing paths to .sig/.sig.gz files")
        p.add_argument('-o', '--output', required=True,
                       help='output file for the index')
        p.add_argument('-t', '--threshold', default=0.01, type=float,
                       help='containment threshold for reporting matches')
        p.add_argument('-k', '--ksize', default=31, type=int,
                       help='k-mer size at which to select sketches')
        p.add_argument('-s', '--scaled', default=1000, type=int,
                       help='scaled factor at which to do comparisons')
        p.add_argument('--save-paths', action='store_true',
                       help='save paths to signatures into index. Default: save full sig into index')
        p.add_argument('-c', '--cores', default=0, type=int,
                       help='number of cores to use (default is all available)')

    def main(self, args):
        notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold: {args.threshold}")

        num_threads = set_thread_pool(args.cores)

        notify(f"indexing all sketches in '{args.siglist}'")

        sig_paths = None
        if args.siglist.endswith('.zip'):
            notify(f"trying to load from zip file '{args.siglist}'")
            from sourmash.index import ZipFileLinearIndex
            zipidx = ZipFileLinearIndex.load(args.siglist, use_manifest=True)
            zip_path = os.path.abspath(args.siglist)
            allfiles = zipidx.storage._filenames()
            sig_paths = [os.path.join(zip_path, x) for x in allfiles if x.endswith('.sig') or x.endswith('.sig.gz')]
            # get paths from zipindex; pass to pyo3_branchwater

        super().main(args)
        status = pyo3_branchwater.do_index(args.siglist,
                                                args.ksize,
                                                args.scaled,
                                                args.threshold,
                                                args.output,
                                                args.save_paths,
                                                False, # colors - currently must be false?
                                                sig_paths)
        if status == 0:
            notify(f"...index is done! results in '{args.output}'")
        return status

class Branchwater_Check(CommandLinePlugin):
    command = 'check'
    description = 'Check Branchwater RevIndex'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('index',
                       help='index file')
        p.add_argument('--quick', action='store_true')

    def main(self, args):
        notify(f"checking index '{args.index}'")
        super().main(args)
        status = pyo3_branchwater.do_check(args.index, args.quick)
        if status == 0:
            notify(f"...index is ok!")
        return status

# class Branchwater_Update(CommandLinePlugin):
#     command = 'update'
#     description = 'Update Branchwater RevIndex'

#     def __init__(self, p):
#         super().__init__(p)
#         p.add_argument('siglist',
#                        help="a text file containing paths to .sig/.sig.gz files")
#         p.add_argument('-o', '--output', required=True,
#                        help='output file for the index')
#         p.add_argument('-t', '--threshold', default=0.01, type=float,
#                        help='containment threshold for reporting matches')
#         p.add_argument('-k', '--ksize', default=31, type=int,
#                        help='k-mer size at which to select sketches')
#         p.add_argument('-s', '--scaled', default=1000, type=int,
#                        help='scaled factor at which to do comparisons')
#         p.add_argument('--save-paths', action='store_true',
#                         help='save paths to signatures into index. Default: save full sig into index')
#        p.add_argument('-c', '--cores', default=0, type=int,
#                         help='number of cores to use (default is all available)')

#     def main(self, args):
#         notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold: {args.threshold}")
#         num_threads = set_thread_pool(args.cores)
#         notify(f"updating index with all sketches in '{args.siglist}'")
#         super().main(args)
#         status = pyo3_branchwater.do_update(args.siglist,
#                                                 args.ksize,
#                                                 args.scaled,
#                                                 args.threshold,
#                                                 args.output,
#                                                 args.save_paths,
#                                                 False) # colors - currently must be false
#         if status == 0:
#             notify(f"...index is done! results in '{args.output}'")
#         return status

class Branchwater_Search(CommandLinePlugin):
    command = 'search'
    description = 'Search against Branchwater RevIndex'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('query_paths',
                       help="a text file containing paths to .sig/.sig.gz files")
        p.add_argument('index',
                       help="an index prepared with pyo3_branchwater's sourmash scripts index")
        p.add_argument('-o', '--output',
                       help='CSV output file for matches')
        p.add_argument('-t', '--threshold-bp', default=50000, type=float,
                       help='threshold in estimated base pairs, for reporting matches (default: 50kb)')
        p.add_argument('--containment-threshold', default=0.01, type=float,
                       help='containment threshold for reporting matches')
        p.add_argument('-k', '--ksize', default=31, type=int,
                       help='k-mer size at which to select sketches')
        p.add_argument('-s', '--scaled', default=1000, type=int,
                       help='scaled factor at which to do comparisons')
        p.add_argument('-c', '--cores', default=0, type=int,
                       help='number of cores to use (default is all available)')

    def main(self, args):
        notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold_bp: {args.threshold_bp}")

        num_threads = set_thread_pool(args.cores)

        notify(f"searching all sketches in '{args.query_paths}' against '{args.index}' using {num_threads} threads")
        super().main(args)
        status = pyo3_branchwater.do_mastiffmanysearch(args.query_paths,
                                              args.index,
                                              args.ksize,
                                              args.scaled,
                                              int(args.threshold_bp),
                                              args.containment_threshold,
                                              args.output)
        if status == 0:
            if args.output:
                notify(f"...search is done! results in '{args.output}'")
            else:
                notify(f"...search is done!")
        return status

class Branchwater_Gather(CommandLinePlugin):
    command = 'gather'
    description = 'Gather against Branchwater RevIndex'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('query_paths',
                       help="a text file containing paths to .sig/.sig.gz files")
        p.add_argument('index',
                       help="an index prepared with pyo3_branchwater's sourmash scripts index")
        p.add_argument('-o', '--output',
                       help='CSV output file for matches')
        p.add_argument('-t', '--threshold-bp', default=50000, type=float,
                       help='threshold in estimated base pairs, for reporting matches (default: 50kb)')
        p.add_argument('-k', '--ksize', default=31, type=int,
                       help='k-mer size at which to select sketches')
        p.add_argument('-s', '--scaled', default=1000, type=int,
                       help='scaled factor at which to do comparisons')
        p.add_argument('-c', '--cores', default=0, type=int,
                       help='number of cores to use (default is all available)')

    def main(self, args):
        notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold_bp: {args.threshold_bp}")

        num_threads = set_thread_pool(args.cores)

        notify(f"gathering all sketches in '{args.query_paths}' against '{args.index}' using {num_threads} threads")
        super().main(args)
        status = pyo3_branchwater.do_mastiffmanygather(args.query_paths,
                                              args.index,
                                              args.ksize,
                                              args.scaled,
                                              int(args.threshold_bp),
                                              args.output)
        if status == 0:
            if args.output:
                notify(f"...gather is done! results in '{args.output}'")
            else:
                notify(f"...gather is done!")
        return status
