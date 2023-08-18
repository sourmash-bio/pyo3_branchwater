import os
import pytest
import pandas

import sourmash_tst_utils as utils


def get_test_data(filename):
    thisdir = os.path.dirname(__file__)
    return os.path.join(thisdir, 'test-data', filename)


def make_file_list(filename, paths):
    with open(filename, 'wt') as fp:
        fp.write("\n".join(paths))
        fp.write("\n")


def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather')

    assert 'usage:  fastgather' in runtmp.last_result.err


def test_simple(runtmp):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '-s', '100000')
    assert os.path.exists(g_output)

    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_file', 'match', 'match_md5sum', 'rank', 'overlap'}


def test_simple_with_prefetch(runtmp):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')
    assert os.path.exists(g_output)
    assert os.path.exists(p_output)

    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_file', 'match', 'match_md5sum', 'rank', 'overlap'}

    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_file', 'match', 'match_md5sum', 'overlap'}


def test_missing_query(runtmp, capfd):
    # test missing query
    query = runtmp.output('no-such-file')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', query, against_list,
                        '-o', g_output, '--output-prefetch', p_output,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: No such file or directory ' in captured.err


def test_bad_query(runtmp, capfd):
    # test non-sig query
    query = runtmp.output('no-such-file')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    # since 'query' needs to be a sig, this breaks it.
    make_file_list(query, [sig2])

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', query, against_list,
                        '-o', g_output, '--output-prefetch', p_output,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: expected value at line 1' in captured.err


def test_missing_against(runtmp, capfd):
    # test missing against
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    #make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', query, against_list,
                        '-o', g_output, '--output-prefetch', p_output,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: No such file or directory ' in captured.err


def test_bad_against(runtmp, capfd):
    # test bad 'against' file - in this case, use a .sig.gz file.
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', query, sig2,
                        '-o', g_output, '--output-prefetch', p_output,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: invalid line in fromfile ' in captured.err