import os
import pandas

import pytest

import sourmash_plugin_branchwater as branch
from . import sourmash_tst_utils as utils
from .sourmash_tst_utils import get_test_data, make_file_list


def test_basic():
    sigfile = get_test_data('SRR606249.sig.gz')
    res = branch.api.api_load_collection(sigfile, 31, 100_000, 'DNA')
    assert res.location == sigfile
    assert len(res) == 1


def test_fail():
    # try to load a (nonexistent) collection
    sigfile = get_test_data('XXX_SRR606249.sig.gz')

    with pytest.raises(RuntimeError):
        res = branch.api.api_load_collection(sigfile, 31, 100_000, 'DNA')
    # @CTB should do something better here than RuntimeError ;)


def test_load_rocksdb(runtmp):
    # test basic index!
    siglist = runtmp.output("db-sigs.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(siglist, [sig2, sig47, sig63])

    output = runtmp.output("db.rocksdb")

    runtmp.sourmash("scripts", "index", siglist, "-o", output)
    assert os.path.exists(output)

    db = branch.api.BranchRevIndex(output)

    assert db.ksize() == 31
    assert db.moltype() == "DNA"
    assert db.min_max_scaled() == (1000, 1000) # @CTB guaranteed to be only one
    assert len(db) == 3


def test_load_rocksdb_get_selection(runtmp):
    # test basic index!
    siglist = runtmp.output("db-sigs.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(siglist, [sig2, sig47, sig63])

    output = runtmp.output("db.rocksdb")

    runtmp.sourmash("scripts", "index", siglist, "-o", output)
    assert os.path.exists(output)

    db = branch.api.BranchRevIndex(output)

    selection = db.selection()

    assert selection.ksize() == 31
    assert selection.moltype() == "DNA"
    assert selection.scaled() == 1000


def test_fastmultigather_rocksdb(runtmp):
    siglist = runtmp.output("db-sigs.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query = get_test_data('SRR606249.sig.gz')

    make_file_list(siglist, [sig2, sig47, sig63])

    output = runtmp.output("db.rocksdb")

    runtmp.sourmash("scripts", "index", siglist, "-o", output)
    assert os.path.exists(output)

    db = branch.api.BranchRevIndex(output)

    query_coll = branch.api.api_load_collection(query, 31, 100_000, 'DNA')

    csv_out = runtmp.output("xxx.csv")
    status = db.fastmultigather_against(query_coll,
                                        db.selection(),
                                        0,
                                        csv_out)
    print(f"status: {status}")
    df = pandas.read_csv(csv_out)
    assert len(df) == 3
    keys = set(df.keys())
    assert {
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "f_match",
        "intersect_bp",
    }.issubset(keys)

    print(df.to_markdown())


def test_basic_collection_load():
    sigfile = get_test_data('SRR606249.sig.gz')
    res = branch.api.api_load_collection(sigfile, 31, 100_000, 'DNA')
