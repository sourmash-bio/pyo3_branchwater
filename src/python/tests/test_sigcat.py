import os
import pytest
import sourmash
from sourmash import index

from . import sourmash_tst_utils as utils

def get_test_data(filename):
    thisdir = os.path.dirname(__file__)
    return os.path.join(thisdir, 'test-data', filename)

def make_file_list(filename, paths):
    with open(filename, 'wt') as fp:
        fp.write("\n".join(paths))
        fp.write("\n")

def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'sigcat')

    assert 'usage:  sigcat' in runtmp.last_result.err

def zip_siglist(runtmp, siglist, db):
    runtmp.sourmash('sig', 'cat', siglist,
                    '-o', db)
    return db


def test_simple_sigcat(runtmp, capfd):
    # test basic execution!
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2])
    make_file_list(against_list, [sig47, sig63])

    output = runtmp.output('out.zip')

    runtmp.sourmash('scripts', 'sigcat', query_list, against_list,
                    '-o', output)
    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    
    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 3
    captured = capfd.readouterr()
    print(captured.out)
    assert f"Concatenated 3 signatures into '{output}'." in captured.out


def test_simple_sigcat_with_zip(runtmp, capfd):
    # test basic execution with a zip involved
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2])
    make_file_list(against_list, [sig47, sig63])

    output = runtmp.output('out.zip')
    # zip the against_list
    against_list = zip_siglist(runtmp, against_list, runtmp.output('db.zip'))

    runtmp.sourmash('scripts', 'sigcat', query_list, against_list,
                    '-o', output)
    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    
    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 3
    captured = capfd.readouterr()
    print(captured.out)
    assert f"Concatenated 3 signatures into '{output}'." in captured.out


def test_output_manifest(runtmp, capfd):
    # test basic manifest-generating functionality.

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'sigcat', sig2, sig47, sig63, '-o', output)

    loader = sourmash.load_file_as_index(output)

    rows = []
    siglist = []
    for (sig, loc) in loader._signatures_with_internal():
        row = index.CollectionManifest.make_manifest_row(sig, loc)
        rows.append(row)
        siglist.append(sig)

    manifest = index.CollectionManifest(rows)

    assert len(manifest) == len(rows)
    assert len(manifest) == 3

    md5_list = [ row['md5'] for row in manifest.rows ]
    assert 'f3a90d4e5528864a5bcc8434b0d0c3b1' in md5_list
    assert '09a08691ce52952152f0e866a59f6261' in md5_list
    assert '38729c6374925585db28916b82a6f513' in md5_list

    for sig in siglist:
        assert sig in manifest
        assert sig.minhash.ksize == 31
        assert sig.minhash.moltype == 'DNA'
        assert sig.minhash.scaled == 1000
    captured = capfd.readouterr()
    print(captured.out)
    assert f"Concatenated 3 signatures into '{output}'." in captured.out


def test_sigcat_missing_sigfile(runtmp, capfd):
    # test missing sig file
    sigf = runtmp.output('sig_pathlist.txt') #define but don't create

    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')
    output = runtmp.output('out.zip')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'sigcat', sig47, sigf, sig63, '-o', output)


    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: No such file or directory:" in captured.err


def test_sigcat_output_not_zipfile(runtmp, capfd):
    # test output not zipfile
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')
    output = runtmp.output('out.sig')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'sigcat', sig47, sig63, '-o', output)


    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: Output file must end with '.zip'" in captured.err


def test_sigcat_incompatible_sigfiles(runtmp, capfd):
    # test incompatible sig files

    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')
    output = runtmp.output('out.zip')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'sigcat', sig47, \
                        sig63, '-o', output, '-k', '9')


    captured = capfd.readouterr()
    print(captured.err)
    assert "WARNING: skipped 1 analysis paths - no compatible signatures." in captured.err
    assert "Error: No analysis signatures loaded, exiting" in captured.err


def test_sigcat_incompatible_sigfiles_force(runtmp, capfd):
    # test incompatible sig files

    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')
    sighp = get_test_data('hp.zip')
    output = runtmp.output('out.zip')

    # with pytest.raises(utils.SourmashCommandFailed):
    runtmp.sourmash('scripts', 'sigcat', sig47, sig63, sighp, \
                    '-o', output, '-m', 'hp', '--force')


    captured = capfd.readouterr()
    print(captured.err)
    assert "WARNING: skipped 1 analysis paths - no compatible signatures." in captured.err
    assert "Error: No analysis signatures loaded, exiting" not in captured.err
    print(captured.out)
    assert f"Concatenated 2 signatures into '{output}'." in captured.out


def test_sigcat_incompatible_sigfiles_force_fail(runtmp, capfd):
    # test incompatible sig files -- not recoverable

    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')
    output = runtmp.output('out.zip')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'sigcat', sig47, sig63, \
                    '-o', output, '-m', 'hp', '--force')


    captured = capfd.readouterr()
    print(captured.err)
    assert "WARNING: skipped 1 analysis paths - no compatible signatures." in captured.err
    assert "Error: No analysis signatures loaded, exiting" not in captured.err
    assert "Error: No signatures to concatenate, exiting" in captured.err


def test_sigcat_bad_input_type(runtmp, capfd):
    # test using bad input (fasta file instead of sig)
    sig47 = get_test_data('47.fa.sig.gz')
    fa1 = get_test_data('short.fa')
    output = runtmp.output('out.zip')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'sigcat', sig47, fa1, '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: No analysis signatures loaded, exiting." in captured.err


def test_sigcat_multk_multsc_multmol(runtmp, capfd):
    # test cat with multiple types of sigs
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig1 = get_test_data('1.fa.k21.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig_hp = get_test_data('hp.zip')
    sig_sc100k = get_test_data('SRR606249.sig.gz')

    output = runtmp.output('out.zip')

    runtmp.sourmash('scripts', 'sigcat', sig1, sig2, sig_hp, sig_sc100k, '-o', output)
    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 5
    captured = capfd.readouterr()
    print(captured.out)
    assert f"Concatenated 5 signatures into '{output}'." in captured.out


def test_sigcat_multk_multsc_multmol_selectk_fail(runtmp, capfd):
    # test cat with multiple types of sigs -- select on ksize
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig1 = get_test_data('1.fa.k21.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig_hp = get_test_data('hp.zip')
    sig_sc100k = get_test_data('SRR606249.sig.gz')

    output = runtmp.output('out.zip')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'sigcat', sig1, sig2, sig_hp, sig_sc100k, '-k', '21','-o', output)

    assert not runtmp.last_result.out # stdout should be empty

    captured = capfd.readouterr()
    print(captured.err)
    assert "WARNING: skipped 1 analysis paths - no compatible signatures." in captured.err
    assert "Error: No analysis signatures loaded, exiting" in captured.err


def test_sigcat_multk_multsc_multmol_selectk_force(runtmp, capfd):
    # test cat with multiple types of sigs -- select on ksize
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig1 = get_test_data('1.fa.k21.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig_hp = get_test_data('hp.zip')
    sig_sc100k = get_test_data('SRR606249.sig.gz')

    output = runtmp.output('out.zip')

    runtmp.sourmash('scripts', 'sigcat', sig1, sig2, sig_hp, sig_sc100k, '-k', '21','-o', output, '--force')
    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 1
    captured = capfd.readouterr()
    print(captured.out)
    assert f"Concatenated 1 signatures into '{output}'." in captured.out


def test_sigcat_multk_multsc_multmol_selectmoltype_fail(runtmp, capfd):
    # test cat with multiple types of sigs -- select on moltype
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig1 = get_test_data('1.fa.k21.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig_hp = get_test_data('hp.zip')
    sig_sc100k = get_test_data('SRR606249.sig.gz')

    output = runtmp.output('out.zip')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'sigcat', sig1, sig2, sig_hp, sig_sc100k, '-m', 'hp','-o', output)

    assert not runtmp.last_result.out # stdout should be empty

    captured = capfd.readouterr()
    print(captured.err)
    assert "WARNING: skipped 1 analysis paths - no compatible signatures." in captured.err
    assert "Error: No analysis signatures loaded, exiting" in captured.err


def test_sigcat_multk_multsc_multmol_selectmoltype_force(runtmp, capfd):
    # test cat with multiple types of sigs -- select on moltype
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig1 = get_test_data('1.fa.k21.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig_hp = get_test_data('hp.zip')
    sig_sc100k = get_test_data('SRR606249.sig.gz')

    output = runtmp.output('out.zip')

    runtmp.sourmash('scripts', 'sigcat', sig1, sig2, sig_hp, sig_sc100k, '-m', 'hp','-o', output, '--force')
    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 2
    captured = capfd.readouterr()
    print(captured.out)
    assert f"Concatenated 2 signatures into '{output}'." in captured.out


def test_sigcat_multk_multsc_multmol_selectscaled_fail(runtmp, capfd):
    # test cat with multiple types of sigs -- select on scaled
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig1 = get_test_data('1.fa.k21.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig_hp = get_test_data('hp.zip')
    sig_sc100k = get_test_data('SRR606249.sig.gz')

    output = runtmp.output('out.zip')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'sigcat', sig1, sig2, sig_hp, sig_sc100k, '--scaled', '100','-o', output)

    assert not runtmp.last_result.out # stdout should be empty

    captured = capfd.readouterr()
    print(captured.err)
    assert "WARNING: skipped 1 analysis paths - no compatible signatures." in captured.err
    assert "Error: No analysis signatures loaded, exiting" in captured.err


def test_sigcat_multk_multsc_multmol_selectscaled_force(runtmp, capfd):
    # test cat with multiple types of sigs -- select on scaled
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig1 = get_test_data('1.fa.k21.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig_hp = get_test_data('hp.zip')
    sig_sc100k = get_test_data('SRR606249.sig.gz')

    output = runtmp.output('out.zip')

    runtmp.sourmash('scripts', 'sigcat', sig1, sig2, sig_hp, sig_sc100k, '--scaled', '150','-o', output, '--force')
    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 2
    for sig in sigs:
        assert sig.minhash.scaled == 150
    captured = capfd.readouterr()
    print(captured.out)
    assert f"Concatenated 2 signatures into '{output}'." in captured.out
