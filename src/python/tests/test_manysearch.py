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
        runtmp.sourmash("scripts", "manysearch")

    assert "usage:  manysearch" in runtmp.last_result.err


def test_simple(runtmp, zip_query, zip_against):
    # test basic execution!
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))
    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output("against.zip"))

    runtmp.sourmash(
        "scripts", "manysearch", query_list, against_list, "-o", output, "-t", "0.01"
    )
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 5

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        # identical?
        if row["match_name"] == row["query_name"]:
            assert row["query_md5"] == row["match_md5"], row
            assert float(row["containment"] == 1.0)
            assert float(row["jaccard"] == 1.0)
            assert float(row["max_containment"] == 1.0)
            assert float(row["query_containment_ani"] == 1.0)
            assert float(row["match_containment_ani"] == 1.0)
            assert float(row["average_containment_ani"] == 1.0)
            assert float(row["max_containment_ani"] == 1.0)

        else:
            # confirm hand-checked numbers
            q = row["query_name"].split()[0]
            m = row["match_name"].split()[0]
            cont = float(row["containment"])
            jaccard = float(row["jaccard"])
            maxcont = float(row["max_containment"])
            intersect_hashes = int(row["intersect_hashes"])
            query_ani = float(row["query_containment_ani"])
            match_ani = float(row["match_containment_ani"])
            average_ani = float(row["average_containment_ani"])
            max_ani = float(row["max_containment_ani"])
            jaccard = round(jaccard, 4)
            cont = round(cont, 4)
            maxcont = round(maxcont, 4)
            query_ani = round(query_ani, 4)
            match_ani = round(match_ani, 4)
            average_ani = round(average_ani, 4)
            max_ani = round(max_ani, 4)
            print(
                q,
                m,
                f"{jaccard:.04}",
                f"{cont:.04}",
                f"{maxcont:.04}",
                f"{query_ani:.04}",
                f"{match_ani:.04}",
                f"{average_ani:.04}",
                f"{max_ani:.04}",
            )

            if q == "NC_011665.1" and m == "NC_009661.1":
                assert jaccard == 0.3207
                assert cont == 0.4828
                assert maxcont == 0.4885
                assert intersect_hashes == 2529
                assert query_ani == 0.9768
                assert match_ani == 0.9772
                assert average_ani == 0.977
                assert max_ani == 0.9772

            if q == "NC_009661.1" and m == "NC_011665.1":
                assert jaccard == 0.3207
                assert cont == 0.4885
                assert maxcont == 0.4885
                assert intersect_hashes == 2529
                assert query_ani == 0.9772
                assert match_ani == 0.9768
                assert average_ani == 0.977
                assert max_ani == 0.9772


def test_simple_output_all(runtmp, zip_query, zip_against):
    # test basic execution!
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))
    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output("against.zip"))

    runtmp.sourmash(
        "scripts",
        "manysearch",
        query_list,
        against_list,
        "-o",
        output,
        "-t",
        "0.01",
        "-A",
    )
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 9

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        # identical?
        if row["match_name"] == row["query_name"]:
            assert row["query_md5"] == row["match_md5"], row
            assert float(row["containment"] == 1.0)
            assert float(row["jaccard"] == 1.0)
            assert float(row["max_containment"] == 1.0)
            assert float(row["query_containment_ani"] == 1.0)
            assert float(row["match_containment_ani"] == 1.0)
            assert float(row["average_containment_ani"] == 1.0)
            assert float(row["max_containment_ani"] == 1.0)

        else:
            # confirm hand-checked numbers
            q = row["query_name"].split()[0]
            m = row["match_name"].split()[0]
            cont = float(row["containment"])
            jaccard = float(row["jaccard"])
            maxcont = float(row["max_containment"])
            intersect_hashes = int(row["intersect_hashes"])
            query_ani = float(row["query_containment_ani"])
            match_ani = float(row["match_containment_ani"])
            average_ani = float(row["average_containment_ani"])
            max_ani = float(row["max_containment_ani"])
            jaccard = round(jaccard, 4)
            cont = round(cont, 4)
            maxcont = round(maxcont, 4)
            query_ani = round(query_ani, 4)
            match_ani = round(match_ani, 4)
            average_ani = round(average_ani, 4)
            max_ani = round(max_ani, 4)
            print(
                q,
                m,
                f"{jaccard:.04}",
                f"{cont:.04}",
                f"{maxcont:.04}",
                f"{query_ani:.04}",
                f"{match_ani:.04}",
                f"{average_ani:.04}",
                f"{max_ani:.04}",
            )

            if q == "NC_011665.1" and m == "NC_009661.1":
                assert jaccard == 0.3207
                assert cont == 0.4828
                assert maxcont == 0.4885
                assert intersect_hashes == 2529
                assert query_ani == 0.9768
                assert match_ani == 0.9772
                assert average_ani == 0.977
                assert max_ani == 0.9772
            elif q == "NC_009661.1" and m == "NC_011665.1":
                assert jaccard == 0.3207
                assert cont == 0.4885
                assert maxcont == 0.4885
                assert intersect_hashes == 2529
                assert query_ani == 0.9772
                assert match_ani == 0.9768
                assert average_ani == 0.977
                assert max_ani == 0.9772
            else:
                assert jaccard == 0
                assert cont == 0
                assert maxcont == 0
                assert intersect_hashes == 0
                assert query_ani == 0
                assert match_ani == 0
                assert average_ani == 0
                assert max_ani == 0


def test_simple_abund(runtmp):
    # test with abund sig
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")
    query_list = runtmp.output("query.txt")
    make_file_list(query_list, [sig2, sig47, sig63])

    against = get_test_data("SRR606249.sig.gz")

    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts",
        "manysearch",
        query_list,
        against,
        "-o",
        output,
        "--scaled",
        "100000",
        "-k",
        "31",
        "-t",
        "0.01",
    )

    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 3

    dd = df.to_dict(orient="index")
    dd = list(sorted(dd.values(), key=lambda x: x["query_name"]))
    print(dd)

    row = dd[0]
    query_name = row["query_name"].split()[0]
    average_abund = round(float(row["average_abund"]), 4)
    median_abund = round(float(row["median_abund"]), 4)
    std_abund = round(float(row["std_abund"]), 4)
    n_weighted_found = int(row["n_weighted_found"])
    total_weighted_hashes = int(row["total_weighted_hashes"])

    assert query_name == "CP001071.1"
    assert average_abund == round(21.045454545454500, 4)
    assert median_abund == 21.5
    assert std_abund == round(5.644605411181010, 4)
    assert n_weighted_found == 463
    assert total_weighted_hashes == 73489

    row = dd[1]
    query_name = row["query_name"].split()[0]
    average_abund = round(float(row["average_abund"]), 4)
    median_abund = round(float(row["median_abund"]), 4)
    std_abund = round(float(row["std_abund"]), 4)
    n_weighted_found = int(row["n_weighted_found"])
    total_weighted_hashes = int(row["total_weighted_hashes"])

    assert query_name == "NC_009661.1"
    assert average_abund == round(11.365853658536600, 4)
    assert median_abund == 11.0
    assert std_abund == round(4.976805212676670, 4)
    assert n_weighted_found == 466
    assert total_weighted_hashes == 73489

    row = dd[2]
    query_name = row["query_name"].split()[0]
    average_abund = round(float(row["average_abund"]), 4)
    median_abund = round(float(row["median_abund"]), 4)
    std_abund = round(float(row["std_abund"]), 4)
    n_weighted_found = int(row["n_weighted_found"])
    total_weighted_hashes = int(row["total_weighted_hashes"])

    assert query_name == "NC_011665.1"
    assert average_abund == round(10.386363636363600, 4)
    assert median_abund == 10.5
    assert std_abund == round(6.932190750047300, 4)
    assert n_weighted_found == 457
    assert total_weighted_hashes == 73489


def test_simple_indexed(runtmp, zip_query, indexed_query):
    # test basic execution!
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    if indexed_query:
        query_list = index_siglist(runtmp, query_list, runtmp.output("query_db"))

    output = runtmp.output("out.csv")

    against_list = index_siglist(runtmp, against_list, runtmp.output("db"))

    print("query_list is:", query_list)
    runtmp.sourmash(
        "scripts", "manysearch", query_list, against_list, "-o", output, "-t", "0.01"
    )
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 5

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        # identical?
        if row["match_name"] == row["query_name"]:
            assert float(row["containment"] == 1.0)
            assert float(row["query_containment_ani"] == 1.0)
        else:
            # confirm hand-checked numbers
            q = row["query_name"].split()[0]
            m = row["match_name"].split()[0]
            cont = float(row["containment"])
            intersect_hashes = int(row["intersect_hashes"])
            query_ani = float(row["query_containment_ani"])
            cont = round(cont, 4)
            query_ani = round(query_ani, 4)
            print(q, m, f"{cont:.04}", f"{query_ani:.04}")

            if q == "NC_011665.1" and m == "NC_009661.1":
                assert cont == 0.4828
                assert intersect_hashes == 2529
                assert query_ani == 0.9768

            if q == "NC_009661.1" and m == "NC_011665.1":
                assert cont == 0.4885
                assert intersect_hashes == 2529
                assert query_ani == 0.9772


def test_simple_list_of_zips(runtmp):
    # test basic execution!
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.sig.zip")
    sig47 = get_test_data("47.sig.zip")
    sig63 = get_test_data("63.sig.zip")

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts", "manysearch", query_list, against_list, "-o", output, "-t", "0.01"
    )
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 5

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        # identical?
        if row["match_name"] == row["query_name"]:
            assert float(row["containment"] == 1.0)
            assert float(row["query_containment_ani"] == 1.0)
        else:
            # confirm hand-checked numbers
            q = row["query_name"].split()[0]
            m = row["match_name"].split()[0]
            cont = float(row["containment"])
            intersect_hashes = int(row["intersect_hashes"])
            query_ani = float(row["query_containment_ani"])
            cont = round(cont, 4)
            query_ani = round(query_ani, 4)
            print(q, m, f"{cont:.04}", f"{query_ani:.04}")

            if q == "NC_011665.1" and m == "NC_009661.1":
                assert cont == 0.4828
                assert intersect_hashes == 2529
                assert query_ani == 0.9768

            if q == "NC_009661.1" and m == "NC_011665.1":
                assert cont == 0.4885
                assert intersect_hashes == 2529
                assert query_ani == 0.9772


def test_simple_with_cores(runtmp, capfd, indexed, zip_query):
    # test basic execution with -c argument (that it runs, at least!)
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts",
        "manysearch",
        query_list,
        against_list,
        "-o",
        output,
        "-c",
        "4",
        "-t",
        "0.01",
    )
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 5

    result = runtmp.last_result
    print(result.err)
    assert " using 4 threads" in result.err


def test_simple_threshold(runtmp, indexed, zip_query):
    # test with a simple threshold => only 3 results
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts", "manysearch", query_list, against_list, "-o", output, "-t", "0.5"
    )
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 3


def test_simple_scaled(runtmp, indexed, zip_query):
    # test with a different (explicitly specified) scaled
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts", "manysearch", query_list, against_list, "-o", output, "-s", "10_000"
    )
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 5
    assert set(list(df["scaled"])) == {10000}


def test_simple_scaled_fail(runtmp, capfd, indexed, zip_query):
    # test with a different scaled, that fails
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")
    against = get_test_data("SRR606249.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [against])

    if indexed:
        against_list = index_siglist(
            runtmp, against_list, runtmp.output("db"), scaled=100_000
        )

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    output = runtmp.output("out.csv")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts",
            "manysearch",
            query_list,
            against_list,
            "-o",
            output,
            "-s",
            "10_000",
        )

    captured = capfd.readouterr()
    print(captured.err)


def test_simple_manifest(runtmp, indexed):
    # test with a simple threshold => only 3 results
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    query_mf = runtmp.output("qmf.csv")
    against_mf = runtmp.output("amf.csv")

    runtmp.sourmash("sig", "manifest", query_list, "-o", query_mf)
    runtmp.sourmash("sig", "manifest", against_list, "-o", against_mf)

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))
    else:
        against_list = against_mf

    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts", "manysearch", query_mf, against_list, "-o", output, "-t", "0.5"
    )
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 3


def test_missing_query(runtmp, capfd, indexed, zip_query):
    # test with a missing query list
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    # make_file_list(query_list, [sig2, sig47, sig63]) # don't make query
    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))

    if zip_query:
        query_list = runtmp.output("query.zip")

    output = runtmp.output("out.csv")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "manysearch", query_list, against_list, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)

    assert "Error: No such file or directory" in captured.err


def test_sig_query(runtmp, capfd, indexed):
    # test with a single sig query (a .sig.gz file)
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))

    output = runtmp.output("out.csv")

    runtmp.sourmash("scripts", "manysearch", sig2, against_list, "-o", output)


def test_bad_query_2(runtmp, capfd, indexed):
    # test with a bad query list (a missing file)
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")
    make_file_list(query_list, [sig2, "no-exist"])
    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))
    output = runtmp.output("out.csv")

    runtmp.sourmash("scripts", "manysearch", query_list, against_list, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)

    assert "WARNING: could not load sketches from path 'no-exist'" in captured.err
    assert (
        "WARNING: 1 query paths failed to load. See error messages above."
        in captured.err
    )


def test_missing_against(runtmp, capfd, indexed):
    # test with a missing against list
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])
    # do not create against_list

    output = runtmp.output("out.csv")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "manysearch", query_list, against_list, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)

    assert "Error: No such file or directory" in captured.err


def test_nomatch_against(runtmp, capfd):
    # nonmatching against file (SRR606249 has scaled=100_000)
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")
    nomatch_sketch = get_test_data("SRR606249.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [nomatch_sketch])

    output = runtmp.output("out.csv")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts",
            "manysearch",
            query_list,
            against_list,
            "-o",
            output,
            "-s",
            "1000",
        )

    captured = capfd.readouterr()
    assert "No search signatures loaded, exiting." in captured.err


def test_bad_against(runtmp, capfd):
    # test with a bad against list (a missing file)
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")
    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, "no-exist"])

    output = runtmp.output("out.csv")

    runtmp.sourmash("scripts", "manysearch", query_list, against_list, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)

    assert "WARNING: could not load sketches from path 'no-exist'" in captured.err
    assert (
        "WARNING: 1 search paths failed to load. See error messages above."
        in captured.err
    )


def test_empty_query(runtmp, indexed, capfd):
    # test with an empty query list
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [])
    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))

    output = runtmp.output("out.csv")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "manysearch", query_list, against_list, "-o", output)

    print(runtmp.last_result.err)
    captured = capfd.readouterr()
    print(captured.err)
    assert "No query signatures loaded, exiting." in captured.err


def test_nomatch_query(runtmp, capfd, indexed, zip_query):
    # test a non-matching (diff ksize) in query; do we get warning message?
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig1 = get_test_data("1.fa.k21.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63, sig1])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")
    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))

    runtmp.sourmash("scripts", "manysearch", query_list, against_list, "-o", output)
    assert os.path.exists(output)

    captured = capfd.readouterr()
    print(captured.err)

    assert "WARNING: skipped 1 query paths - no compatible signatures." in captured.err


def test_load_only_one_bug(runtmp, capfd, indexed, zip_against):
    # check that we behave properly when presented with multiple against
    # sketches
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig1_k31 = get_test_data("1.fa.k31.sig.gz")

    # note: this was created as a 3-sketch-in-one-signature directly
    # via sourmash sketch dna -p k=21,k=31,k=51.
    sig1_all = get_test_data("1.combined.sig.gz")

    make_file_list(query_list, [sig1_k31])
    make_file_list(against_list, [sig1_all])

    output = runtmp.output("out.csv")

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output("against.zip"))
    elif indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))

    runtmp.sourmash("scripts", "manysearch", query_list, against_list, "-o", output)
    assert os.path.exists(output)

    captured = capfd.readouterr()
    print(captured.err)

    assert not "WARNING: skipped 1 paths - no compatible signatures." in captured.err
    assert not "WARNING: no compatible sketches in path " in captured.err


def test_load_only_one_bug_as_query(runtmp, capfd, indexed, zip_query):
    # check that we behave properly when presented with multiple query
    # sketches in one file, with only one matching.
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig1_k31 = get_test_data("1.fa.k31.sig.gz")

    # note: this was created as a 3-sketch-in-one-signature directly
    # via sourmash sketch dna -p k=21,k=31,k=51.
    sig1_all = get_test_data("1.combined.sig.gz")

    make_file_list(query_list, [sig1_all])
    make_file_list(against_list, [sig1_k31])

    output = runtmp.output("out.csv")

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))
    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    runtmp.sourmash("scripts", "manysearch", query_list, against_list, "-o", output)

    assert os.path.exists(output)

    captured = capfd.readouterr()
    print(captured.err)
    print(runtmp.last_result.out)

    assert not "WARNING: skipped 1 paths - no compatible signatures." in captured.err
    assert not "WARNING: no compatible sketches in path " in captured.err


def test_md5(runtmp, indexed, zip_query):
    # test that md5s match what was in the original files, not downsampled etc.
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output("db"))

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    runtmp.sourmash(
        "scripts", "manysearch", query_list, against_list, "-o", output, "-t", "0.01"
    )
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 5

    md5s = list(df["query_md5"])
    print(md5s)

    for query_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(query_file, ksize=31):
            assert ss.md5sum() in md5s

    if not indexed:  # indexed search cannot produce match_md5
        md5s = list(df["match_md5"])
        print(md5s)

        for against_file in (sig2, sig47, sig63):
            for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
                assert ss.md5sum() in md5s


def test_simple_protein(runtmp):
    # test basic execution with proteins
    protsigs = get_test_data("protein.zip")
    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts",
        "manysearch",
        protsigs,
        protsigs,
        "-k",
        "19",
        "-s",
        "100",
        "--moltype",
        "protein",
        "-o",
        output,
        "-t",
        "0.01",
    )

    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 4

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        print(row)
        # identical?
        if row["match_name"] == row["query_name"]:
            assert row["query_md5"] == row["match_md5"], row
            assert float(row["containment"] == 1.0)
            assert float(row["jaccard"] == 1.0)
            assert float(row["max_containment"] == 1.0)
            assert float(row["query_containment_ani"]) == 1.0
            assert float(row["match_containment_ani"]) == 1.0
            assert float(row["average_containment_ani"]) == 1.0
            assert float(row["max_containment_ani"]) == 1.0
        else:
            # confirm hand-checked numbers
            q = row["query_name"].split()[0]
            m = row["match_name"].split()[0]
            cont = float(row["containment"])
            jaccard = float(row["jaccard"])
            maxcont = float(row["max_containment"])
            intersect_hashes = int(row["intersect_hashes"])
            query_ani = float(row["query_containment_ani"])
            match_ani = float(row["match_containment_ani"])
            average_ani = float(row["average_containment_ani"])
            max_ani = float(row["max_containment_ani"])

            jaccard = round(jaccard, 4)
            cont = round(cont, 4)
            maxcont = round(maxcont, 4)
            query_ani = round(query_ani, 4)
            match_ani = round(match_ani, 4)
            average_ani = round(average_ani, 4)
            max_ani = round(max_ani, 4)
            print(
                q,
                m,
                f"{jaccard:.04}",
                f"{cont:.04}",
                f"{maxcont:.04}",
                intersect_hashes,
                f"{query_ani:.04}",
                f"{match_ani:.04}",
                f"{average_ani:.04}",
                f"{max_ani:.04}",
            )

            if q == "GCA_001593925" and m == "GCA_001593935":
                assert jaccard == 0.0434
                assert cont == 0.1003
                assert maxcont == 0.1003
                assert intersect_hashes == 342
                assert query_ani == 0.9605
                assert match_ani == 0.9547
                assert average_ani == 0.9576
                assert max_ani == 0.9605

            if q == "GCA_001593935" and m == "GCA_001593925":
                assert jaccard == 0.0434
                assert cont == 0.0712
                assert maxcont == 0.1003
                assert intersect_hashes == 342
                assert query_ani == 0.9547
                assert match_ani == 0.9605
                assert average_ani == 0.9576
                assert max_ani == 0.9605


def test_simple_protein_indexed(runtmp):
    # test basic execution with proteins
    protsigs = get_test_data("protein.zip")
    output = runtmp.output("out.csv")

    protsigs_db = index_siglist(
        runtmp, protsigs, runtmp.output("db"), ksize=19, moltype="protein", scaled=100
    )

    runtmp.sourmash(
        "scripts",
        "manysearch",
        protsigs,
        protsigs_db,
        "-k",
        "19",
        "-s",
        "100",
        "--moltype",
        "protein",
        "-o",
        output,
        "-t",
        "0.01",
    )

    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 4

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        print(row)
        # identical?
        if row["match_name"] == row["query_name"]:
            assert float(row["containment"] == 1.0)
            assert float(row["query_containment_ani"] == 1.0)
        else:
            # confirm hand-checked numbers
            q = row["query_name"].split()[0]
            m = row["match_name"].split()[0]
            cont = float(row["containment"])
            query_ani = float(row["query_containment_ani"])
            intersect_hashes = int(row["intersect_hashes"])

            cont = round(cont, 4)
            query_ani = round(query_ani, 4)
            print(q, m, f"{cont:.04}", intersect_hashes, f"{query_ani:.04}")

            if q == "GCA_001593925" and m == "GCA_001593935":
                assert cont == 0.1003
                assert intersect_hashes == 342
                assert query_ani == 0.9605

            if q == "GCA_001593935" and m == "GCA_001593925":
                assert cont == 0.0712
                assert intersect_hashes == 342
                assert query_ani == 0.9547


def test_simple_dayhoff(runtmp):
    # test basic execution with dayhoff
    protsigs = get_test_data("dayhoff.zip")
    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts",
        "manysearch",
        protsigs,
        protsigs,
        "-k",
        "19",
        "-s",
        "100",
        "--moltype",
        "dayhoff",
        "-o",
        output,
        "-t",
        "0.01",
    )

    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 4

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        print(row)
        # identical?
        if row["match_name"] == row["query_name"]:
            assert row["query_md5"] == row["match_md5"], row
            assert float(row["containment"] == 1.0)
            assert float(row["jaccard"] == 1.0)
            assert float(row["max_containment"] == 1.0)
            assert float(row["query_containment_ani"]) == 1.0
            assert float(row["match_containment_ani"]) == 1.0
            assert float(row["average_containment_ani"]) == 1.0
            assert float(row["max_containment_ani"]) == 1.0

        else:
            # confirm hand-checked numbers
            q = row["query_name"].split()[0]
            m = row["match_name"].split()[0]
            cont = float(row["containment"])
            jaccard = float(row["jaccard"])
            maxcont = float(row["max_containment"])
            intersect_hashes = int(row["intersect_hashes"])
            query_ani = float(row["query_containment_ani"])
            match_ani = float(row["match_containment_ani"])
            average_ani = float(row["average_containment_ani"])
            max_ani = float(row["max_containment_ani"])

            jaccard = round(jaccard, 4)
            cont = round(cont, 4)
            maxcont = round(maxcont, 4)
            query_ani = round(query_ani, 4)
            match_ani = round(match_ani, 4)
            average_ani = round(average_ani, 4)
            max_ani = round(max_ani, 4)
            print(
                q,
                m,
                f"{jaccard:.04}",
                f"{cont:.04}",
                f"{maxcont:.04}",
                intersect_hashes,
                f"{query_ani:.04}",
                f"{match_ani:.04}",
                f"{average_ani:.04}",
                f"{max_ani:.04}",
            )

            if q == "GCA_001593925" and m == "GCA_001593935":
                assert jaccard == 0.1326
                assert cont == 0.2815
                assert maxcont == 0.2815
                assert intersect_hashes == 930
                assert query_ani == 0.978
                assert match_ani == 0.9722
                assert average_ani == 0.9751
                assert max_ani == 0.978

            if q == "GCA_001593935" and m == "GCA_001593925":
                assert jaccard == 0.1326
                assert cont == 0.2004
                assert maxcont == 0.2815
                assert intersect_hashes == 930
                assert query_ani == 0.9722
                assert match_ani == 0.978
                assert average_ani == 0.9751
                assert max_ani == 0.978


def test_simple_dayhoff_indexed(runtmp):
    # test indexed execution with dayhoff
    protsigs = get_test_data("dayhoff.zip")
    output = runtmp.output("out.csv")

    protsigs_db = index_siglist(
        runtmp, protsigs, runtmp.output("db"), ksize=19, moltype="dayhoff", scaled=100
    )

    runtmp.sourmash(
        "scripts",
        "manysearch",
        protsigs,
        protsigs_db,
        "-k",
        "19",
        "-s",
        "100",
        "--moltype",
        "dayhoff",
        "-o",
        output,
        "-t",
        "0.01",
    )

    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 4

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        print(row)
        # identical?
        if row["match_name"] == row["query_name"]:
            assert float(row["containment"] == 1.0)
            assert float(row["query_containment_ani"] == 1.0)
        else:
            # confirm hand-checked numbers
            q = row["query_name"].split()[0]
            m = row["match_name"].split()[0]
            cont = float(row["containment"])
            query_ani = float(row["query_containment_ani"])
            intersect_hashes = int(row["intersect_hashes"])

            cont = round(cont, 4)
            query_ani = round(query_ani, 4)
            print(q, m, f"{cont:.04}", intersect_hashes, f"{query_ani:.04}")

            if q == "GCA_001593925" and m == "GCA_001593935":
                assert cont == 0.2815
                assert intersect_hashes == 930
                assert query_ani == 0.978

            if q == "GCA_001593935" and m == "GCA_001593925":
                assert cont == 0.2004
                assert intersect_hashes == 930
                assert query_ani == 0.9722


def test_simple_hp(runtmp):
    # test basic execution with hp
    protsigs = get_test_data("hp.zip")
    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts",
        "manysearch",
        protsigs,
        protsigs,
        "-k",
        "19",
        "-s",
        "100",
        "--moltype",
        "hp",
        "-o",
        output,
        "-t",
        "0.01",
    )

    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 4

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        print(row)
        # identical?
        if row["match_name"] == row["query_name"]:
            assert row["query_md5"] == row["match_md5"], row
            assert float(row["containment"] == 1.0)
            assert float(row["jaccard"] == 1.0)
            assert float(row["max_containment"] == 1.0)
            assert float(row["query_containment_ani"]) == 1.0
            assert float(row["match_containment_ani"]) == 1.0
            assert float(row["average_containment_ani"]) == 1.0
            assert float(row["max_containment_ani"]) == 1.0
        else:
            # confirm hand-checked numbers
            q = row["query_name"].split()[0]
            m = row["match_name"].split()[0]
            cont = float(row["containment"])
            jaccard = float(row["jaccard"])
            maxcont = float(row["max_containment"])
            intersect_hashes = int(row["intersect_hashes"])
            query_ani = float(row["query_containment_ani"])
            match_ani = float(row["match_containment_ani"])
            average_ani = float(row["average_containment_ani"])
            max_ani = float(row["max_containment_ani"])

            jaccard = round(jaccard, 4)
            cont = round(cont, 4)
            maxcont = round(maxcont, 4)
            query_ani = round(query_ani, 4)
            match_ani = round(match_ani, 4)
            average_ani = round(average_ani, 4)
            max_ani = round(max_ani, 4)
            print(
                q,
                m,
                f"{jaccard:.04}",
                f"{cont:.04}",
                f"{maxcont:.04}",
                intersect_hashes,
                f"{query_ani:.04}",
                f"{match_ani:.04}",
                f"{average_ani:.04}",
                f"{max_ani:.04}",
            )

            if q == "GCA_001593925" and m == "GCA_001593935":
                assert jaccard == 0.4983
                assert cont == 0.747
                assert maxcont == 0.747
                assert intersect_hashes == 1724
                assert query_ani == 0.9949
                assert match_ani == 0.9911
                assert average_ani == 0.993
                assert max_ani == 0.9949

            if q == "GCA_001593935" and m == "GCA_001593925":
                assert jaccard == 0.4983
                assert cont == 0.5994
                assert maxcont == 0.747
                assert intersect_hashes == 1724
                assert query_ani == 0.9911
                assert match_ani == 0.9949
                assert average_ani == 0.993
                assert max_ani == 0.9949


def test_simple_hp_indexed(runtmp):
    # test indexed execution with hp, indexed
    protsigs = get_test_data("hp.zip")
    output = runtmp.output("out.csv")

    protsigs_db = index_siglist(
        runtmp, protsigs, runtmp.output("db"), ksize=19, moltype="hp", scaled=100
    )

    runtmp.sourmash(
        "scripts",
        "manysearch",
        protsigs,
        protsigs_db,
        "-k",
        "19",
        "-s",
        "100",
        "--moltype",
        "hp",
        "-o",
        output,
        "-t",
        "0.01",
    )

    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 4

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        print(row)
        # identical?
        if row["match_name"] == row["query_name"]:
            assert float(row["containment"] == 1.0)
            assert float(row["query_containment_ani"]) == 1.0

        else:
            # confirm hand-checked numbers
            q = row["query_name"].split()[0]
            m = row["match_name"].split()[0]
            cont = float(row["containment"])
            intersect_hashes = int(row["intersect_hashes"])
            query_ani = float(row["query_containment_ani"])

            cont = round(cont, 4)
            query_ani = round(query_ani, 4)
            print(q, m, f"{cont:.04}", intersect_hashes, f"{query_ani:.04}")

            if q == "GCA_001593925" and m == "GCA_001593935":
                assert cont == 0.747
                assert intersect_hashes == 1724
                assert query_ani == 0.9949

            if q == "GCA_001593935" and m == "GCA_001593925":
                assert cont == 0.5994
                assert intersect_hashes == 1724
                assert query_ani == 0.9911


def test_pretty_print(runtmp):
    # test pretty-printing of output
    query = get_test_data("hmp-queries.sig.zip")
    against = get_test_data("hmp-against.sig.zip")

    outcsv = runtmp.output("xxx.csv")

    runtmp.sourmash("scripts", "manysearch", query, against, "-o", outcsv)
    print(runtmp.last_result.out)

    expected = """\
query             p_genome avg_abund   p_metag   metagenome name
--------          -------- ---------   -------   ---------------
B. fragilis I1345   96.7%     7.3      27.5%     CD136
B. fragilis I1345   96.7%     7.5      22.6%     CD237
F. prausnitzii      58.4%    25.3      30.7%     CD136
"""
    for line in expected.splitlines():
        assert expected in runtmp.last_result.out, expected


def test_no_pretty_print(runtmp):
    # test turning off pretty-printing of output
    query = get_test_data("hmp-queries.sig.zip")
    against = get_test_data("hmp-against.sig.zip")

    outcsv = runtmp.output("xxx.csv")

    runtmp.sourmash("scripts", "manysearch", query, against, "-o", outcsv, "-N")
    print(runtmp.last_result.out)

    # if this fails in the future, it might be because the order of the
    # output gets shuffled by multithreading. consider refactoring to
    # do line by line?
    expected = "p_genome"
    assert expected not in runtmp.last_result.out


def test_bug_550(runtmp):
    # check a bug where a manifest made from a .sig file causes problems
    # due to a problem with the way Signature::name() behaved in sourmash
    # before r0.18.0.
    # see https://github.com/sourmash-bio/sourmash_plugin_branchwater/issues/550
    fa_file = get_test_data("short.fa")
    sig_out = runtmp.output("short.sig")
    mf_out = runtmp.output("short.mf.csv")
    csv_out = runtmp.output("out.csv")

    runtmp.sourmash("sketch", "dna", fa_file, "-o", sig_out)
    runtmp.sourmash("sig", "collect", "-F", "csv", sig_out, "-o", mf_out)
    runtmp.sourmash("scripts", "manysearch", mf_out, mf_out, "-o", csv_out)

    assert os.path.exists(csv_out)


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
            "manysearch",
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
        "manysearch",
        sig2,
        "podar-ref-subset.branch0_9_5.rocksdb",
        "-o",
        "out.csv",
        "-s",
        "100_000",
        in_dir=runtmp.output(""),
    )
