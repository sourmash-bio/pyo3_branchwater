"""
Test 'sourmash scripts fastmultigather'
"""

import os
import pytest
import pandas
import shutil

import sourmash
from . import sourmash_tst_utils as utils
from .sourmash_tst_utils import (
    get_test_data,
    make_file_list,
    zip_siglist,
    index_siglist,
)


def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "fastmultigather")

    assert "usage:  fastmultigather" in runtmp.last_result.err


def test_simple(runtmp, zip_against):
    # test basic execution!
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output("against.zip"))

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_list,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        runtmp.output("SRR606249.gather.csv"),
        in_directory=runtmp.output(""),
    )

    print(os.listdir(runtmp.output("")))

    g_output = runtmp.output("SRR606249.gather.csv")
    p_output = runtmp.output("SRR606249.prefetch.csv")
    assert os.path.exists(p_output)

    # check prefetch output (only non-indexed gather)
    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
    }

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    print(df)
    assert len(df) == 3
    keys = set(df.keys())
    assert {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
        "gather_result_rank",
    }.issubset(keys)


def test_simple_list_of_zips(runtmp):
    # test basic execution!
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.sig.zip")
    sig47 = get_test_data("47.sig.zip")
    sig63 = get_test_data("63.sig.zip")

    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_list,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        runtmp.output("SRR606249.gather.csv"),
        in_dir=runtmp.output(""),
    )

    print(os.listdir(runtmp.output("")))

    g_output = runtmp.output("SRR606249.gather.csv")
    p_output = runtmp.output("SRR606249.prefetch.csv")
    assert os.path.exists(p_output)

    # check prefetch output (only non-indexed gather)
    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
    }

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    print(df)
    assert len(df) == 3
    keys = set(df.keys())
    assert {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
        "gather_result_rank",
    }.issubset(keys)


def test_simple_space_in_signame(runtmp):
    # test basic execution!
    query = get_test_data("SRR606249.sig.gz")
    renamed_query = runtmp.output("in.zip")
    name = "my-favorite-signame has spaces"
    # rename signature
    runtmp.sourmash("sig", "rename", query, name, "-o", renamed_query)

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    against_list = runtmp.output("against.txt")

    make_file_list(against_list, [sig2, sig47, sig63])

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        renamed_query,
        against_list,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        runtmp.output("my-favorite-signame.gather.csv"),
        in_directory=runtmp.output(""),
    )

    print(os.listdir(runtmp.output("")))

    g_output = runtmp.output("my-favorite-signame.gather.csv")
    p_output = runtmp.output("my-favorite-signame.prefetch.csv")
    assert os.path.exists(p_output)
    assert os.path.exists(g_output)


def test_simple_zip_query(runtmp):
    # test basic execution!
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_list,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        runtmp.output("SRR606249.gather.csv"),
        in_directory=runtmp.output(""),
    )

    print(os.listdir(runtmp.output("")))

    g_output = runtmp.output("SRR606249.gather.csv")
    p_output = runtmp.output("SRR606249.prefetch.csv")

    # check prefetch output (only non-indexed gather)
    assert os.path.exists(p_output)
    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
    }

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
        "gather_result_rank",
    }.issubset(keys)


def test_simple_read_manifests(runtmp):
    # test basic execution!
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    against_list = runtmp.output("against.txt")
    against_mf = runtmp.output("against.csv")
    query_mf = runtmp.output("query.csv")

    make_file_list(against_list, [sig2, sig47, sig63])

    runtmp.sourmash("sig", "manifest", query, "-o", query_mf)
    runtmp.sourmash("sig", "manifest", against_list, "-o", against_mf)

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_mf,
        against_list,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        runtmp.output("SRR606249.gather.csv"),
        in_directory=runtmp.output(""),
    )

    print(os.listdir(runtmp.output("")))

    g_output = runtmp.output("SRR606249.gather.csv")
    p_output = runtmp.output("SRR606249.prefetch.csv")

    # check prefetch output (only non-indexed gather)
    assert os.path.exists(p_output)
    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
    }

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
        "gather_result_rank",
    }.issubset(keys)


def test_simple_indexed(runtmp, zip_query, toggle_internal_storage):
    # test basic execution!
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    g_output = runtmp.output("out.csv")
    against_db = index_siglist(
        runtmp,
        against_list,
        runtmp.output("test.rocksdb"),
        toggle_internal_storage=toggle_internal_storage,
    )
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_db,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        g_output,
    )

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    expected_keys = {
        "match_name",
        "query_filename",
        "query_n_hashes",
        "match_filename",
        "f_match_orig",
        "query_bp",
        "query_abundance",
        "match_containment_ani",
        "intersect_bp",
        "total_weighted_hashes",
        "n_unique_weighted_found",
        "query_name",
        "gather_result_rank",
        "moltype",
        "query_containment_ani",
        "sum_weighted_found",
        "f_orig_query",
        "ksize",
        "max_containment_ani",
        "std_abund",
        "scaled",
        "average_containment_ani",
        "f_match",
        "f_unique_to_query",
        "average_abund",
        "unique_intersect_bp",
        "median_abund",
        "query_md5",
        "match_md5",
        "remaining_bp",
        "f_unique_weighted",
    }
    assert keys == expected_keys


def test_simple_indexed_query_manifest(runtmp, toggle_internal_storage):
    # test basic execution!
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_mf = runtmp.output("query.csv")
    against_list = runtmp.output("against.txt")

    make_file_list(against_list, [sig2, sig47, sig63])
    runtmp.sourmash("sig", "manifest", query, "-o", query_mf)

    g_output = runtmp.output("out.csv")
    against_db = index_siglist(
        runtmp,
        against_list,
        runtmp.output("db"),
        toggle_internal_storage=toggle_internal_storage,
    )
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_mf,
        against_db,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        g_output,
    )

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    expected_keys = {
        "match_name",
        "query_filename",
        "query_n_hashes",
        "match_filename",
        "f_match_orig",
        "query_bp",
        "query_abundance",
        "match_containment_ani",
        "intersect_bp",
        "total_weighted_hashes",
        "n_unique_weighted_found",
        "query_name",
        "gather_result_rank",
        "moltype",
        "query_containment_ani",
        "sum_weighted_found",
        "f_orig_query",
        "ksize",
        "max_containment_ani",
        "std_abund",
        "scaled",
        "average_containment_ani",
        "f_match",
        "f_unique_to_query",
        "average_abund",
        "unique_intersect_bp",
        "median_abund",
        "query_md5",
        "match_md5",
        "remaining_bp",
        "f_unique_weighted",
    }
    assert keys == expected_keys


def test_missing_querylist(runtmp, capfd, indexed, zip_query, toggle_internal_storage):
    # test missing querylist
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    if zip_query:
        query_list = runtmp.output("query.zip")
    # do not make query_list!
    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(
            runtmp,
            against_list,
            runtmp.output("db"),
            toggle_internal_storage=toggle_internal_storage,
        )

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts", "fastmultigather", query_list, against_list, "-s", "100000"
        )

    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: No such file or directory" in captured.err


def test_sig_query(runtmp, capfd, indexed):
    # sig file is now fine as a query
    query = get_test_data("SRR606249.sig.gz")

    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output("out.csv")
    output_params = ["-o", g_output]

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))
    else:
        p_output = runtmp.output("SRR606249.prefetch.csv")

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query,
        against_list,
        "-s",
        "100000",
        *output_params,
    )

    captured = capfd.readouterr()
    print(captured.err)
    if not indexed:
        # check prefetch output (only non-indexed gather)
        assert os.path.exists(p_output)
        df = pandas.read_csv(p_output)
        assert len(df) == 3
        keys = set(df.keys())
        assert {
            "query_filename",
            "query_name",
            "query_md5",
            "match_name",
            "match_md5",
            "intersect_bp",
        }.issubset(keys)

    # check gather output (both)
    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    if indexed:
        assert {
            "query_name",
            "query_md5",
            "match_name",
            "match_md5",
            "f_match",
            "intersect_bp",
        }.issubset(keys)
    else:
        assert {
            "query_filename",
            "query_name",
            "query_md5",
            "match_name",
            "match_md5",
            "gather_result_rank",
            "intersect_bp",
        }.issubset(keys)


def test_missing_query(runtmp, capfd, indexed):
    # test missing query
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, "no-exist"])
    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))

    runtmp.sourmash(
        "scripts", "fastmultigather", query_list, against_list, "-s", "100000"
    )

    captured = capfd.readouterr()
    print(captured.err)
    assert "WARNING: could not load sketches from path 'no-exist'" in captured.err
    assert "WARNING: 1 query paths failed to load. See error messages above."


def test_nomatch_query(runtmp, capfd, indexed, zip_query):
    # test nomatch file in querylist
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")
    badsig1 = get_test_data("1.fa.k21.sig.gz")

    make_file_list(query_list, [sig2, badsig1])
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))
    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))

    runtmp.sourmash(
        "scripts", "fastmultigather", query_list, against_list, "-s", "100000"
    )

    captured = capfd.readouterr()
    print(captured.err)
    assert "WARNING: skipped 1 query paths - no compatible signatures." in captured.err


def test_missing_against(runtmp, capfd, zip_against):
    # test missing against
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])

    if zip_against:
        against_list = runtmp.output("against.zip")
    # do not make against_list

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts", "fastmultigather", query_list, against_list, "-s", "100000"
        )

    captured = capfd.readouterr()
    print(captured.err)

    assert "Error: No such file or directory" in captured.err


def test_sig_against(runtmp, capfd):
    # against file can be a sig now
    query = get_test_data("SRR606249.sig.gz")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")

    g_output = runtmp.output("SRR606249.gather.csv")
    p_output = runtmp.output("SRR606249.prefetch.csv")
    runtmp.sourmash(
        "scripts", "fastmultigather", query, sig2, "-s", "100000", "-o", g_output
    )

    captured = capfd.readouterr()
    print(captured.err)

    # check prefetch output (only non-indexed gather)
    assert os.path.exists(p_output)
    df = pandas.read_csv(p_output)
    assert len(df) == 1
    keys = set(df.keys())
    assert {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
    }.issubset(keys)

    # check gather output
    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 1
    keys = set(df.keys())
    assert {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
        "gather_result_rank",
    }.issubset(keys)


def test_bad_against(runtmp, capfd):
    # test bad 'against' file - in this case, one containing a nonexistent file
    query = get_test_data("SRR606249.sig.gz")
    query_list = runtmp.output("query.txt")
    make_file_list(query_list, [query])

    against_list = runtmp.output("against.txt")
    sig2 = get_test_data("2.fa.sig.gz")
    make_file_list(against_list, [sig2, "no exist"])

    # should succeed, but with error output.
    runtmp.sourmash(
        "scripts", "fastmultigather", query_list, against_list, "-s", "100000"
    )

    captured = capfd.readouterr()
    print(captured.err)

    assert "WARNING: could not load sketches from path 'no exist'" in captured.err
    assert (
        "WARNING: 1 search paths failed to load. See error messages above."
        in captured.err
    )


def test_empty_against(runtmp, capfd):
    # test bad 'against' file - in this case, an empty one
    query = get_test_data("SRR606249.sig.gz")
    query_list = runtmp.output("query.txt")
    make_file_list(query_list, [query])

    against_list = runtmp.output("against.txt")
    make_file_list(against_list, [])

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts", "fastmultigather", query_list, against_list, "-s", "100000"
        )

    captured = capfd.readouterr()
    print(captured.err)

    assert "Sketch loading error: No such file or directory" in captured.err
    assert "No search signatures loaded, exiting." in captured.err


def test_nomatch_in_against(runtmp, capfd, zip_against):
    # test an against file that has a non-matching ksize sig in it
    query = get_test_data("SRR606249.sig.gz")
    query_list = runtmp.output("query.txt")
    make_file_list(query_list, [query])

    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig1 = get_test_data("1.fa.k21.sig.gz")
    make_file_list(against_list, [sig2, sig1])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output("against.zip"))

    runtmp.sourmash(
        "scripts", "fastmultigather", query_list, against_list, "-s", "100000"
    )

    captured = capfd.readouterr()
    print(captured.err)

    assert "WARNING: skipped 1 search paths - no compatible signatures." in captured.err


def test_md5(runtmp, zip_query):
    # test correct md5s present in output
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_list,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        runtmp.output("SRR606249.gather.csv"),
        in_directory=runtmp.output(""),
    )

    print(os.listdir(runtmp.output("")))

    g_output = runtmp.output("SRR606249.gather.csv")
    p_output = runtmp.output("SRR606249.prefetch.csv")

    # check prefetch output (only non-indexed gather)
    assert os.path.exists(p_output)
    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
    }

    md5s = set(df["match_md5"])
    for against_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
            assert ss.md5sum() in md5s

    # check gather output (mostly same for indexed vs non-indexed version)
    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
    }.issubset(keys)

    md5s = set(df["match_md5"])
    for against_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
            assert ss.md5sum() in md5s


def test_md5_indexed(runtmp, zip_query):
    # test correct md5s present in output
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    g_output = runtmp.output("out.csv")
    against_list = index_siglist(runtmp, against_list, runtmp.output("db"))
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_list,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        g_output,
    )

    # check gather output (mostly same for indexed vs non-indexed version)
    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    expected_keys = {
        "match_name",
        "query_filename",
        "query_n_hashes",
        "match_filename",
        "f_match_orig",
        "query_bp",
        "query_abundance",
        "match_containment_ani",
        "intersect_bp",
        "total_weighted_hashes",
        "n_unique_weighted_found",
        "query_name",
        "gather_result_rank",
        "moltype",
        "query_containment_ani",
        "sum_weighted_found",
        "f_orig_query",
        "ksize",
        "max_containment_ani",
        "std_abund",
        "scaled",
        "average_containment_ani",
        "f_match",
        "f_unique_to_query",
        "average_abund",
        "unique_intersect_bp",
        "median_abund",
        "query_md5",
        "match_md5",
        "remaining_bp",
        "f_unique_weighted",
    }
    assert keys == expected_keys

    md5s = set(df["match_md5"])
    for against_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
            assert ss.md5sum() in md5s


def test_csv_columns_vs_sourmash_prefetch(runtmp, zip_query, zip_against):
    # the column names should be strict subsets of sourmash prefetch cols
    query = get_test_data("SRR606249.sig.gz")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    make_file_list(query_list, [query])
    against_list = runtmp.output("against.txt")
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))
    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output("against.zip"))

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_list,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        runtmp.output("SRR606249.gather.csv"),
        in_directory=runtmp.output(""),
    )

    g_output = runtmp.output("SRR606249.gather.csv")
    p_output = runtmp.output("SRR606249.prefetch.csv")

    assert os.path.exists(p_output)
    assert os.path.exists(g_output)
    # now run sourmash prefetch
    sp_output = runtmp.output("sourmash-prefetch.csv")
    runtmp.sourmash(
        "prefetch", query, against_list, "-o", sp_output, "--scaled", "100000"
    )

    gather_df = pandas.read_csv(g_output)
    g_keys = set(gather_df.keys())
    assert {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "gather_result_rank",
        "intersect_bp",
    }.issubset(g_keys)
    g_keys.remove("gather_result_rank")  # 'rank' is not in sourmash prefetch!

    sourmash_prefetch_df = pandas.read_csv(sp_output)
    sp_keys = set(sourmash_prefetch_df.keys())
    print(g_keys - sp_keys)
    diff_keys = g_keys - sp_keys
    assert diff_keys == set(
        [
            "remaining_bp",
            "f_match_orig",
            "f_unique_weighted",
            "average_abund",
            "unique_intersect_bp",
            "std_abund",
            "sum_weighted_found",
            "median_abund",
            "n_unique_weighted_found",
            "f_unique_to_query",
            "f_orig_query",
            "total_weighted_hashes",
            "f_match",
        ]
    )


def test_csv_columns_vs_sourmash_gather_fullresults(runtmp):
    # the column names should be identical to sourmash gather cols
    query = get_test_data("SRR606249.sig.gz")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    make_file_list(query_list, [query])
    against_list = runtmp.output("against.txt")
    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output("SRR606249.gather.csv")
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_list,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        runtmp.output("SRR606249.gather.csv"),
    )

    assert os.path.exists(g_output)
    # now run sourmash gather
    sg_output = runtmp.output(".csv")
    runtmp.sourmash(
        "gather", query, against_list, "-o", sg_output, "--scaled", "100000"
    )

    gather_df = pandas.read_csv(g_output)
    g_keys = set(gather_df.keys())
    expected_keys = {
        "match_name",
        "query_filename",
        "query_n_hashes",
        "match_filename",
        "f_match_orig",
        "query_bp",
        "query_abundance",
        "match_containment_ani",
        "intersect_bp",
        "total_weighted_hashes",
        "n_unique_weighted_found",
        "query_name",
        "gather_result_rank",
        "moltype",
        "query_containment_ani",
        "sum_weighted_found",
        "f_orig_query",
        "ksize",
        "max_containment_ani",
        "std_abund",
        "scaled",
        "average_containment_ani",
        "f_match",
        "f_unique_to_query",
        "average_abund",
        "unique_intersect_bp",
        "median_abund",
        "query_md5",
        "match_md5",
        "remaining_bp",
        "f_unique_weighted",
    }
    assert g_keys == expected_keys

    sourmash_gather_df = pandas.read_csv(sg_output)
    sg_keys = set(sourmash_gather_df.keys())
    print(sg_keys)
    modified_keys = ["match_md5", "match_name", "match_filename"]
    sg_keys.update(
        modified_keys
    )  # fastmultigather is more explicit (match_md5 instead of md5, etc)
    print("g_keys - sg_keys:", g_keys - sg_keys)
    assert not g_keys - sg_keys, g_keys - sg_keys


def test_csv_columns_vs_sourmash_gather_indexed(runtmp):
    # the column names should be identical to sourmash gather cols
    query = get_test_data("SRR606249.sig.gz")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    make_file_list(query_list, [query])
    against_list = runtmp.output("against.txt")
    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output("out.csv")
    against_db = index_siglist(runtmp, against_list, runtmp.output("db"))
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_db,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        g_output,
    )

    assert os.path.exists(g_output)
    # now run sourmash gather
    sg_output = runtmp.output("sourmash-gather.csv")
    runtmp.sourmash(
        "gather", query, against_list, "-o", sg_output, "--scaled", "100000"
    )

    gather_df = pandas.read_csv(g_output)
    g_keys = set(gather_df.keys())
    expected_keys = {
        "match_name",
        "query_filename",
        "query_n_hashes",
        "match_filename",
        "f_match_orig",
        "query_bp",
        "query_abundance",
        "match_containment_ani",
        "intersect_bp",
        "total_weighted_hashes",
        "n_unique_weighted_found",
        "query_name",
        "gather_result_rank",
        "moltype",
        "query_containment_ani",
        "sum_weighted_found",
        "f_orig_query",
        "ksize",
        "max_containment_ani",
        "std_abund",
        "scaled",
        "average_containment_ani",
        "f_match",
        "f_unique_to_query",
        "average_abund",
        "unique_intersect_bp",
        "median_abund",
        "query_md5",
        "match_md5",
        "remaining_bp",
        "f_unique_weighted",
    }
    assert g_keys == expected_keys

    sourmash_gather_df = pandas.read_csv(sg_output)
    sg_keys = set(sourmash_gather_df.keys())
    print(sg_keys)
    modified_keys = ["match_md5", "match_name", "match_filename"]
    sg_keys.update(
        modified_keys
    )  # fastmultigather is more explicit (match_md5 instead of md5, etc)
    print("g_keys - sg_keys:", g_keys - sg_keys)
    assert not g_keys - sg_keys, g_keys - sg_keys


def test_simple_protein(runtmp):
    # test basic protein execution
    sigs = get_test_data("protein.zip")

    sig_names = ["GCA_001593935", "GCA_001593925"]

    gather_out = runtmp.output("xxx.csv")
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        sigs,
        sigs,
        "-s",
        "100",
        "--moltype",
        "protein",
        "-k",
        "19",
        "-o",
        gather_out,
    )

    all_df = pandas.read_csv(gather_out)
    for qsig in sig_names:
        p_output = runtmp.output(os.path.join(qsig + ".prefetch.csv"))
        assert os.path.exists(p_output)

        df = all_df[all_df["match_name"] == qsig]

        assert len(df) == 1
        keys = set(df.keys())
        assert {
            "query_filename",
            "query_name",
            "query_md5",
            "match_name",
            "match_md5",
            "intersect_bp",
            "gather_result_rank",
        }.issubset(keys)
        print(df)
        # since we're just matching to identical sigs, the md5s should be the same
        assert set(df["query_md5"]) == set(df["match_md5"])


def test_simple_dayhoff(runtmp):
    # test basic protein execution
    sigs = get_test_data("dayhoff.zip")

    sig_names = ["GCA_001593935", "GCA_001593925"]

    gather_out = runtmp.output("xxx.csv")
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        sigs,
        sigs,
        "-s",
        "100",
        "--moltype",
        "dayhoff",
        "-k",
        "19",
        "-o",
        gather_out,
    )

    all_df = pandas.read_csv(gather_out)
    for qsig in sig_names:
        p_output = runtmp.output(os.path.join(qsig + ".prefetch.csv"))
        assert os.path.exists(p_output)

        df = all_df[all_df["match_name"] == qsig]

        assert len(df) == 1
        keys = set(df.keys())
        assert {
            "query_filename",
            "query_name",
            "query_md5",
            "match_name",
            "match_md5",
            "intersect_bp",
            "gather_result_rank",
        }.issubset(keys)
        print(df)
        # since we're just matching to identical sigs, the md5s should be the same
        assert set(df["query_md5"]) == set(df["match_md5"])


def test_simple_hp(runtmp):
    # test basic protein execution
    sigs = get_test_data("hp.zip")

    sig_names = ["GCA_001593935", "GCA_001593925"]

    gather_out = runtmp.output("xxx.csv")
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        sigs,
        sigs,
        "-s",
        "100",
        "--moltype",
        "hp",
        "-k",
        "19",
        "-o",
        gather_out,
    )

    all_df = pandas.read_csv(gather_out)
    for qsig in sig_names:
        p_output = runtmp.output(os.path.join(qsig + ".prefetch.csv"))
        assert os.path.exists(p_output)

        df = all_df[all_df["match_name"] == qsig]
        assert len(df) == 1
        keys = set(df.keys())
        assert {
            "query_filename",
            "query_name",
            "query_md5",
            "match_name",
            "match_md5",
            "intersect_bp",
            "gather_result_rank",
        }.issubset(keys)
        print(df)
        # since we're just matching to identical sigs, the md5s should be the same
        assert set(df["query_md5"]) == set(df["match_md5"])


def test_simple_protein_indexed(runtmp):
    # test basic protein execution
    sigs = get_test_data("protein.zip")

    sigs_db = index_siglist(
        runtmp, sigs, runtmp.output("db"), ksize=19, moltype="protein", scaled=100
    )

    out_csv = runtmp.output("out.csv")
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        sigs,
        sigs_db,
        "-s",
        "100",
        "--moltype",
        "protein",
        "-k",
        "19",
        "-o",
        out_csv,
    )

    assert os.path.exists(out_csv)

    df = pandas.read_csv(out_csv)
    assert len(df) == 2
    keys = set(df.keys())
    expected_keys = {
        "match_name",
        "query_filename",
        "query_n_hashes",
        "match_filename",
        "f_match_orig",
        "query_bp",
        "query_abundance",
        "match_containment_ani",
        "intersect_bp",
        "total_weighted_hashes",
        "n_unique_weighted_found",
        "query_name",
        "gather_result_rank",
        "moltype",
        "query_containment_ani",
        "sum_weighted_found",
        "f_orig_query",
        "ksize",
        "max_containment_ani",
        "std_abund",
        "scaled",
        "average_containment_ani",
        "f_match",
        "f_unique_to_query",
        "average_abund",
        "unique_intersect_bp",
        "median_abund",
        "query_md5",
        "match_md5",
        "remaining_bp",
        "f_unique_weighted",
    }
    assert keys == expected_keys
    print(df)
    # since we're just matching to identical sigs, the md5s should be the same
    assert df["query_md5"][0] == df["match_md5"][0]
    assert df["query_md5"][1] == df["match_md5"][1]


def test_simple_dayhoff_indexed(runtmp):
    # test basic protein execution
    sigs = get_test_data("dayhoff.zip")

    sigs_db = index_siglist(
        runtmp, sigs, runtmp.output("db"), ksize=19, moltype="dayhoff", scaled=100
    )

    out_csv = runtmp.output("out.csv")
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        sigs,
        sigs_db,
        "-s",
        "100",
        "--moltype",
        "dayhoff",
        "-k",
        "19",
        "-o",
        out_csv,
    )

    assert os.path.exists(out_csv)

    df = pandas.read_csv(out_csv)
    assert len(df) == 2
    keys = set(df.keys())
    expected_keys = {
        "match_name",
        "query_filename",
        "query_n_hashes",
        "match_filename",
        "f_match_orig",
        "query_bp",
        "query_abundance",
        "match_containment_ani",
        "intersect_bp",
        "total_weighted_hashes",
        "n_unique_weighted_found",
        "query_name",
        "gather_result_rank",
        "moltype",
        "query_containment_ani",
        "sum_weighted_found",
        "f_orig_query",
        "ksize",
        "max_containment_ani",
        "std_abund",
        "scaled",
        "average_containment_ani",
        "f_match",
        "f_unique_to_query",
        "average_abund",
        "unique_intersect_bp",
        "median_abund",
        "query_md5",
        "match_md5",
        "remaining_bp",
        "f_unique_weighted",
    }
    assert keys == expected_keys
    print(df)
    # since we're just matching to identical sigs, the md5s should be the same
    assert df["query_md5"][0] == df["match_md5"][0]
    assert df["query_md5"][1] == df["match_md5"][1]


def test_simple_hp_indexed(runtmp):
    # test basic protein execution
    sigs = get_test_data("hp.zip")

    sigs_db = index_siglist(
        runtmp, sigs, runtmp.output("db"), ksize=19, moltype="hp", scaled=100
    )

    out_csv = runtmp.output("out.csv")
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        sigs,
        sigs_db,
        "-s",
        "100",
        "--moltype",
        "hp",
        "-k",
        "19",
        "-o",
        out_csv,
    )

    assert os.path.exists(out_csv)

    df = pandas.read_csv(out_csv)
    assert len(df) == 2
    keys = set(df.keys())
    expected_keys = {
        "match_name",
        "query_filename",
        "query_n_hashes",
        "match_filename",
        "f_match_orig",
        "query_bp",
        "query_abundance",
        "match_containment_ani",
        "intersect_bp",
        "total_weighted_hashes",
        "n_unique_weighted_found",
        "query_name",
        "gather_result_rank",
        "moltype",
        "query_containment_ani",
        "sum_weighted_found",
        "f_orig_query",
        "ksize",
        "max_containment_ani",
        "std_abund",
        "scaled",
        "average_containment_ani",
        "f_match",
        "f_unique_to_query",
        "average_abund",
        "unique_intersect_bp",
        "median_abund",
        "query_md5",
        "match_md5",
        "remaining_bp",
        "f_unique_weighted",
    }
    assert keys == expected_keys
    print(df)
    # since we're just matching to identical sigs, the md5s should be the same
    assert df["query_md5"][0] == df["match_md5"][0]
    assert df["query_md5"][1] == df["match_md5"][1]


def test_indexed_full_output(runtmp):
    # test correct md5s present in output
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output("out.csv")
    against_db = index_siglist(runtmp, against_list, runtmp.output("rocksdb"))
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_db,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        g_output,
    )

    # check full gather output
    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    expected_keys = {
        "match_name",
        "query_filename",
        "query_n_hashes",
        "match_filename",
        "f_match_orig",
        "query_bp",
        "query_abundance",
        "match_containment_ani",
        "intersect_bp",
        "total_weighted_hashes",
        "n_unique_weighted_found",
        "query_name",
        "gather_result_rank",
        "moltype",
        "query_containment_ani",
        "sum_weighted_found",
        "f_orig_query",
        "ksize",
        "max_containment_ani",
        "std_abund",
        "scaled",
        "average_containment_ani",
        "f_match",
        "f_unique_to_query",
        "average_abund",
        "unique_intersect_bp",
        "median_abund",
        "query_md5",
        "match_md5",
        "remaining_bp",
        "f_unique_weighted",
    }
    assert keys == expected_keys
    results = df.values.tolist()

    # check a few columns
    average_ani = set(df["average_containment_ani"])
    avg_ani = set([round(x, 4) for x in average_ani])
    assert avg_ani == {0.9221, 0.9306, 0.9316}

    f_unique_weighted = set(df["f_unique_weighted"])
    f_unique_weighted = set([round(x, 4) for x in f_unique_weighted])
    assert f_unique_weighted == {0.0063, 0.002, 0.0062}

    unique_intersect_bp = set(df["unique_intersect_bp"])
    unique_intersect_bp = set([round(x, 4) for x in unique_intersect_bp])
    assert unique_intersect_bp == {4400000, 1800000, 2200000}


def test_nonindexed_full_vs_sourmash_gather(runtmp):
    query = get_test_data("SRR606249.sig.gz")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    make_file_list(query_list, [query])
    against_list = runtmp.output("against.txt")
    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output("SRR606249.gather.csv")
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_list,
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        g_output,
    )

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert os.path.exists(g_output)
    # now run sourmash gather
    sg_output = runtmp.output(".csv")
    runtmp.sourmash(
        "gather", query, against_list, "-o", sg_output, "--scaled", "100000"
    )

    gather_df = pandas.read_csv(g_output)
    g_keys = set(gather_df.keys())

    sourmash_gather_df = pandas.read_csv(sg_output)
    sg_keys = set(sourmash_gather_df.keys())
    print(sg_keys)
    modified_keys = ["match_md5", "match_name", "match_filename"]
    sg_keys.update(
        modified_keys
    )  # fastmultigather is more explicit (match_md5 instead of md5, etc)
    print("g_keys - sg_keys:", g_keys - sg_keys)
    assert not g_keys - sg_keys, g_keys - sg_keys

    for index, row in sourmash_gather_df.iterrows():
        print(row.to_dict())

    fmg_intersect_bp = set(gather_df["intersect_bp"])
    g_intersect_bp = set(sourmash_gather_df["intersect_bp"])
    assert fmg_intersect_bp == g_intersect_bp == set([4400000, 4100000, 2200000])

    fmg_f_orig_query = set([round(x, 4) for x in gather_df["f_orig_query"]])
    g_f_orig_query = set([round(x, 4) for x in sourmash_gather_df["f_orig_query"]])
    assert fmg_f_orig_query == g_f_orig_query == set([0.0098, 0.0105, 0.0052])

    fmg_f_match = set([round(x, 4) for x in gather_df["f_match"]])
    g_f_match = set([round(x, 4) for x in sourmash_gather_df["f_match"]])
    assert fmg_f_match == g_f_match == set([0.439, 1.0])

    fmg_f_unique_to_query = set(
        [round(x, 3) for x in gather_df["f_unique_to_query"]]
    )  # rounding to 4 --> slightly different!
    g_f_unique_to_query = set(
        [round(x, 3) for x in sourmash_gather_df["f_unique_to_query"]]
    )
    assert fmg_f_unique_to_query == g_f_unique_to_query == set([0.004, 0.01, 0.005])

    fmg_f_unique_weighted = set([round(x, 4) for x in gather_df["f_unique_weighted"]])
    g_f_unique_weighted = set(
        [round(x, 4) for x in sourmash_gather_df["f_unique_weighted"]]
    )
    assert fmg_f_unique_weighted == g_f_unique_weighted == set([0.0063, 0.002, 0.0062])

    fmg_average_abund = set([round(x, 4) for x in gather_df["average_abund"]])
    g_average_abund = set([round(x, 4) for x in sourmash_gather_df["average_abund"]])
    assert fmg_average_abund == g_average_abund == set([8.2222, 10.3864, 21.0455])

    fmg_median_abund = set([round(x, 4) for x in gather_df["median_abund"]])
    g_median_abund = set([round(x, 4) for x in sourmash_gather_df["median_abund"]])
    assert fmg_median_abund == g_median_abund == set([8.0, 10.5, 21.5])

    fmg_std_abund = set([round(x, 4) for x in gather_df["std_abund"]])
    g_std_abund = set([round(x, 4) for x in sourmash_gather_df["std_abund"]])
    assert fmg_std_abund == g_std_abund == set([3.172, 5.6446, 6.9322])

    g_match_filename_basename = [
        os.path.basename(filename) for filename in sourmash_gather_df["filename"]
    ]
    fmg_match_filename_basename = [
        os.path.basename(filename) for filename in gather_df["match_filename"]
    ]
    assert all(
        [
            x in fmg_match_filename_basename
            for x in ["2.fa.sig.gz", "63.fa.sig.gz", "47.fa.sig.gz"]
        ]
    )
    assert fmg_match_filename_basename == g_match_filename_basename

    assert list(sourmash_gather_df["name"]) == list(gather_df["match_name"])
    assert list(sourmash_gather_df["md5"]) == list(gather_df["match_md5"])

    fmg_f_match_orig = set([round(x, 4) for x in gather_df["f_match_orig"]])
    g_f_match_orig = set([round(x, 4) for x in sourmash_gather_df["f_match_orig"]])
    assert fmg_f_match_orig == g_f_match_orig == set([1.0])

    fmg_unique_intersect_bp = set(gather_df["unique_intersect_bp"])
    g_unique_intersect_bp = set(sourmash_gather_df["unique_intersect_bp"])
    assert (
        fmg_unique_intersect_bp
        == g_unique_intersect_bp
        == set([4400000, 1800000, 2200000])
    )

    fmg_gather_result_rank = set(gather_df["gather_result_rank"])
    g_gather_result_rank = set(sourmash_gather_df["gather_result_rank"])
    assert fmg_gather_result_rank == g_gather_result_rank == set([0, 1, 2])

    fmg_remaining_bp = list(gather_df["remaining_bp"])
    assert fmg_remaining_bp == [415600000, 413400000, 411600000]
    ### Gather remaining bp does not match, but I think this one is right?
    # g_remaining_bp = list(sourmash_gather_df['remaining_bp'])
    # print("gather remaining bp: ", g_remaining_bp) #{4000000, 0, 1800000}
    # assert fmg_remaining_bp == g_remaining_bp == set([])

    fmg_query_containment_ani = set(
        [round(x, 4) for x in gather_df["query_containment_ani"]]
    )
    g_query_containment_ani = set(
        [round(x, 4) for x in sourmash_gather_df["query_containment_ani"]]
    )
    assert fmg_query_containment_ani == {0.8442, 0.8613, 0.8632}
    # gather cANI are nans here -- perhaps b/c sketches too small
    # assert fmg_query_containment_ani == g_query_containment_ani == set([0.8632, 0.8444, 0.8391])
    print("fmg qcANI: ", fmg_query_containment_ani)
    print("g_qcANI: ", g_query_containment_ani)

    fmg_n_unique_weighted_found = set(gather_df["n_unique_weighted_found"])
    g_n_unique_weighted_found = set(sourmash_gather_df["n_unique_weighted_found"])
    assert (
        fmg_n_unique_weighted_found == g_n_unique_weighted_found == set([457, 148, 463])
    )

    fmg_sum_weighted_found = set(gather_df["sum_weighted_found"])
    g_sum_weighted_found = set(sourmash_gather_df["sum_weighted_found"])
    assert fmg_sum_weighted_found == g_sum_weighted_found == set([920, 457, 1068])

    fmg_total_weighted_hashes = set(gather_df["total_weighted_hashes"])
    g_total_weighted_hashes = set(sourmash_gather_df["total_weighted_hashes"])
    assert fmg_total_weighted_hashes == g_total_weighted_hashes == set([73489])


def test_rocksdb_gather_against_index_with_sigs(runtmp, capfd):
    # fastmultigather should succeed if indexed sigs are stored internally.
    query = get_test_data("SRR606249.sig.gz")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")
    shutil.copyfile(sig2, runtmp.output("2.fa.sig.gz"))
    shutil.copyfile(sig47, runtmp.output("47.fa.sig.gz"))
    shutil.copyfile(sig63, runtmp.output("63.fa.sig.gz"))

    query_list = runtmp.output("query.txt")
    make_file_list(query_list, [query])
    against_list = runtmp.output("against.txt")
    make_file_list(against_list, ["2.fa.sig.gz", "47.fa.sig.gz", "63.fa.sig.gz"])

    # index! note: '--internal-storage' defaults to True
    runtmp.sourmash("scripts", "index", against_list, "-o", "subdir/against.rocksdb")

    # remove the external storage out from under the rocksdb.
    os.unlink(runtmp.output("2.fa.sig.gz"))
    os.unlink(runtmp.output("47.fa.sig.gz"))
    os.unlink(runtmp.output("63.fa.sig.gz"))

    g_output = runtmp.output("zzz.csv")

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        "subdir/against.rocksdb",
        "-s",
        "100000",
        "-t",
        "0",
        "-o",
        g_output,
        in_location=runtmp.output(""),
    )


def test_rocksdb_no_internal_storage_gather_fails(runtmp, capfd):
    # force gather to fail b/c we make an index with no internal sketches
    query = get_test_data("SRR606249.sig.gz")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")
    shutil.copyfile(sig2, runtmp.output("2.fa.sig.gz"))
    shutil.copyfile(sig47, runtmp.output("47.fa.sig.gz"))
    shutil.copyfile(sig63, runtmp.output("63.fa.sig.gz"))

    query_list = runtmp.output("query.txt")
    make_file_list(query_list, [query])
    against_list = runtmp.output("against.txt")
    make_file_list(against_list, ["2.fa.sig.gz", "47.fa.sig.gz", "63.fa.sig.gz"])

    runtmp.sourmash(
        "scripts",
        "index",
        against_list,
        "--no-internal-storage",
        "-o",
        "subdir/against.rocksdb",
    )

    # remove the external storage out from under the rocksdb.
    # this will make gather fail.
    os.unlink(runtmp.output("2.fa.sig.gz"))
    os.unlink(runtmp.output("47.fa.sig.gz"))
    os.unlink(runtmp.output("63.fa.sig.gz"))

    g_output = runtmp.output("zzz.csv")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts",
            "fastmultigather",
            query_list,
            "subdir/against.rocksdb",
            "-s",
            "100000",
            "-t",
            "0",
            "-o",
            g_output,
            in_location=runtmp.output(""),
        )

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    captured = capfd.readouterr()
    print(captured.err)

    assert "Error gathering matches:" in captured.err
    assert "1 failed gathers. See error messages above." in captured.err


def test_save_matches(runtmp):
    # test basic execution!
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_list,
        "-s",
        "100000",
        "-t",
        "0",
        "--save-matches",
        "-o",
        runtmp.output("SRR606249.gather.csv"),
        in_directory=runtmp.output(""),
    )

    print(os.listdir(runtmp.output("")))

    g_output = runtmp.output("SRR606249.gather.csv")
    p_output = runtmp.output("SRR606249.prefetch.csv")
    m_output = runtmp.output("SRR606249.matches.sig")

    assert os.path.exists(g_output)
    assert os.path.exists(p_output)
    assert os.path.exists(m_output)

    # check prefetch output (only non-indexed gather)
    df = pandas.read_csv(p_output)

    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
    }

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)

    assert len(df) == 3
    keys = set(df.keys())
    assert {
        "query_filename",
        "query_name",
        "query_md5",
        "match_name",
        "match_md5",
        "intersect_bp",
        "gather_result_rank",
    }.issubset(keys)

    # can't test against prefetch because matched k-mers can overlap
    match_ss = list(sourmash.load_file_as_signatures(m_output, ksize=31))[0]
    match_mh = match_ss.minhash
    matches_sig_len = len(match_mh)

    # right size?
    assert sum(df["intersect_bp"]) >= matches_sig_len * 100_000

    # containment?
    mg_ss = list(sourmash.load_file_as_signatures(query, ksize=31))[0]
    assert match_mh.contained_by(mg_ss.minhash) == 1.0
    assert mg_ss.minhash.contained_by(match_mh) < 1


def test_create_empty_prefetch_results(runtmp):
    # sig2 has 0 hashes in common with 47 and 63
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    make_file_list(query_list, [sig2])
    make_file_list(against_list, [sig47, sig63])

    gather_out = runtmp.output("SRR606249.gather.csv")
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_list,
        "-s",
        "100000",
        "-t",
        "0",
        "--create-empty-results",
        "-o",
        gather_out,
        in_directory=runtmp.output(""),
    )

    print(os.listdir(runtmp.output("")))

    p_output = runtmp.output("CP001071.1.prefetch.csv")
    assert os.path.exists(p_output)


def test_simple_against_scaled(runtmp, zip_against):
    # we shouldn't automatically downsample query
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    downsampled_sigs = runtmp.output("ds.sig.zip")
    runtmp.sourmash(
        "sig",
        "downsample",
        "--scaled",
        "120_000",
        sig2,
        sig47,
        sig63,
        "-o",
        downsampled_sigs,
    )

    query_list = runtmp.output("query.txt")
    make_file_list(query_list, [query])

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts",
            "fastmultigather",
            query_list,
            downsampled_sigs,
            "-t",
            "0",
            in_directory=runtmp.output(""),
        )


def test_simple_query_scaled(runtmp):
    # test basic execution w/automatic scaled selection based on query
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_list,
        "-t",
        "0",
        "-o",
        runtmp.output("SRR606249.gather.csv"),
        in_directory=runtmp.output(""),
    )

    print(os.listdir(runtmp.output("")))

    g_output = runtmp.output("SRR606249.gather.csv")
    assert os.path.exists(g_output)


def test_exit_no_against(runtmp, indexed):
    # test that it exits properly when nothing to search
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(
            runtmp,
            against_list,
            runtmp.output("db"),
        )

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts",
            "fastmultigather",
            query_list,
            against_list,
            "-s",
            "1000",
            in_directory=runtmp.output(""),
        )


def test_simple_query_scaled_indexed(runtmp):
    # test basic execution w/automatic scaled selection based on query
    # (on a rocksdb)
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])
    against_list = index_siglist(
        runtmp, against_list, runtmp.output("against.rocksdb"), scaled=1000
    )

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts",
            "fastmultigather",
            query_list,
            against_list,
            "-o",
            "foo.csv",
            "-t",
            "0",
            in_directory=runtmp.output(""),
        )


def test_equal_matches(runtmp, indexed):
    # check that equal matches get returned from fastmultigather
    # NOTE: the use of a bunch of bottom hashes in the artifical sketches
    # below makes some of the numbers weird if you downsample etc. So
    # be careful!
    base = sourmash.MinHash(scaled=1, ksize=31, n=0)

    a = base.copy_and_clear()
    b = base.copy_and_clear()
    c = base.copy_and_clear()

    a.add_many(range(0, 1000))
    b.add_many(range(1000, 2000))
    c.add_many(range(0, 2000))

    ss = sourmash.SourmashSignature(a, name="g_a")
    sourmash.save_signatures([ss], open(runtmp.output("a.sig"), "wb"))
    ss = sourmash.SourmashSignature(b, name="g_b")
    sourmash.save_signatures([ss], open(runtmp.output("b.sig"), "wb"))
    ss = sourmash.SourmashSignature(c, name="g_mg")
    sourmash.save_signatures([ss], open(runtmp.output("mg.sig"), "wb"))

    against_list = runtmp.output("combined.sig.zip")
    runtmp.sourmash("sig", "cat", "a.sig", "b.sig", "-o", against_list)

    outfile = runtmp.output("g_mg.gather.csv")
    if indexed:
        against_list = index_siglist(
            runtmp,
            against_list,
            runtmp.output("db"),
        )

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        "mg.sig",
        against_list,
        "--threshold-bp=0",
        "-o",
        outfile,
    )

    df = pandas.read_csv(runtmp.output(outfile))
    assert len(df) == 2
    assert set(df["intersect_bp"]) == {1000}


def test_explicit_scaled(runtmp, indexed):
    # check that an explicit downsampling with -s is respected.
    query = get_test_data("SRR606249.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])
    against_list = zip_siglist(runtmp, against_list, runtmp.output("against.zip"))

    outfile = runtmp.output("SRR606249.gather.csv")
    if indexed:
        against_list = index_siglist(
            runtmp,
            against_list,
            runtmp.output("db"),
        )

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        query_list,
        against_list,
        "-s",
        "150000",
        "-t",
        "0",
        "-o",
        outfile,
        in_directory=runtmp.output(""),
    )

    print(os.listdir(runtmp.output("")))

    assert os.path.exists(outfile)
    df = pandas.read_csv(outfile)
    print(df)
    assert len(df) == 3
    assert set(list(df["scaled"])) == {150_000}
    assert round(df["f_unique_to_query"].sum(), 6) == round(0.01836514223, 6)


def test_rocksdb_v0_9_5(runtmp):
    # there was a RevIndex format change between this plugin v0.9.5 and
    # v0.9.12; test that databases can be opened etc.

    sig2 = get_test_data("2.fa.sig.gz")

    rocksdb_dir = get_test_data("rocksdb/podar-ref-subset.branch0_9_5.rocksdb")
    rocksdb_zip = get_test_data("rocksdb/podar-ref-subset.sig.zip")

    target_rocksdb = runtmp.output("podar-ref-subset.branch0_9_5.rocksdb")
    target_zip = runtmp.output("podar-ref-subset.sig.zip")
    shutil.copytree(rocksdb_dir, runtmp.output(target_rocksdb))
    shutil.copyfile(rocksdb_zip, target_zip)

    # fails b/c Column family not found: storage
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts",
            "fastmultigather",
            sig2,
            "podar-ref-subset.branch0_9_5.rocksdb",
            "-o",
            "out.csv",
            "-s",
            "100_000",
            in_dir=runtmp.output(""),
        )

    # upgrade with 'check'
    runtmp.sourmash(
        "scripts",
        "check",
        "--upgrade",
        "podar-ref-subset.branch0_9_5.rocksdb",
        in_dir=runtmp.output(""),
    )

    # should now work...
    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        sig2,
        "podar-ref-subset.branch0_9_5.rocksdb",
        "-o",
        "out.csv",
        "-s",
        "100_000",
        in_dir=runtmp.output(""),
    )


def test_rocksdb_v0_9_13_internal(runtmp):
    # test that databases created with v0.9.13 w/internal storage can be
    # opened/searched.
    sig2 = get_test_data("2.fa.sig.gz")

    rocksdb_dir = get_test_data(
        "rocksdb/podar-ref-subset.branch0_9_13.internal.rocksdb"
    )

    # no need to copy zip file: internal storage
    target_rocksdb = runtmp.output("podar-ref-subset.branch0_9_13.internal.rocksdb")
    shutil.copytree(rocksdb_dir, runtmp.output(target_rocksdb))

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        sig2,
        "podar-ref-subset.branch0_9_13.internal.rocksdb",
        "-o",
        "out.csv",
        "-s",
        "100_000",
        in_dir=runtmp.output(""),
    )

    assert os.path.exists(runtmp.output("out.csv"))


def test_rocksdb_v0_9_13_external(runtmp):
    # test that databases created with v0.9.13 w/xternal storage can be
    # opened/searched.
    sig2 = get_test_data("2.fa.sig.gz")

    rocksdb_dir = get_test_data(
        "rocksdb/podar-ref-subset.branch0_9_13.external.rocksdb"
    )
    rocksdb_zip = get_test_data("rocksdb/podar-ref-subset.sig.zip")

    # note: copy zipfile b/c external storage.
    target_rocksdb = runtmp.output("podar-ref-subset.branch0_9_13.external.rocksdb")
    target_zip = runtmp.output("podar-ref-subset.sig.zip")
    shutil.copytree(rocksdb_dir, runtmp.output(target_rocksdb))
    shutil.copyfile(rocksdb_zip, target_zip)

    runtmp.sourmash(
        "scripts",
        "fastmultigather",
        sig2,
        "podar-ref-subset.branch0_9_13.external.rocksdb",
        "-o",
        "out.csv",
        "-s",
        "100_000",
        in_dir=runtmp.output(""),
    )

    assert os.path.exists(runtmp.output("out.csv"))
