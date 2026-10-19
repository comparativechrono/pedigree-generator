import subprocess
import sys

import pandas as pd
import pytest

from PedGen import gen_genotype_export, pedigree_generator
from PedGen.cli import main


def read_ped(path):
    return pd.read_csv(path, sep='\t', header=None,
                       names=['FamilyID', 'IndividualID', 'PaternalID', 'MaternalID', 'Sex', 'Phenotype'])


def test_writes_a_ped_file_named_after_the_family(tmp_path, capsys):
    assert main(['--mode', 'AR', '-g', '4', '-c', '3', '-f', 'FAM2', '-s', '7', '-d', str(tmp_path)]) == 0
    ped = read_ped(tmp_path / 'FAM2.ped')
    assert set(ped['FamilyID']) == {'FAM2'} and len(ped) > 4
    assert set(ped['Sex']) <= {1, 2} and set(ped['Phenotype']) <= {1, 2}
    out = capsys.readouterr().out
    assert out.startswith('Family FAM2: ') and 'wrote' in out


def test_the_same_seed_gives_the_same_family(tmp_path):
    for name in ('a', 'b'):
        main(['-s', '3', '-g', '3', '-d', str(tmp_path / name)])
    assert (tmp_path / 'a' / 'FAM1.ped').read_text() == (tmp_path / 'b' / 'FAM1.ped').read_text()


def test_as_python_does(tmp_path):
    # the command writes what pedigree_generator and gen_PED_export write
    main(['--mode', 'AD', '-g', '3', '-c', '2', '-s', '11', '--partners', '0.5', '--ancestors', '0.1',
          '--allele-frequency', '0.05', '--unaffected-partners', '-d', str(tmp_path / 'command')])
    df = pedigree_generator(max_children=2, FamilyID='FAM1', mode='AD', generation_count=3, SpouseLikelihood=0.5,
                            AffectedSpouse=False, BackpropLikelihood=0.1, alt_freq=0.05, seed=11)
    from PedGen import gen_PED_export
    path = gen_PED_export(df, str(tmp_path / 'python'))
    assert (tmp_path / 'command' / 'FAM1.ped').read_text() == open(path).read()


def test_cousins_and_min_affected_as_python_does(tmp_path):
    assert main(['--mode', 'AR', '-g', '5', '-s', '3', '--cousins', '1', '--min-affected', '2',
                 '-d', str(tmp_path / 'command')]) == 0
    assert (read_ped(tmp_path / 'command' / 'FAM1.ped')['Phenotype'] == 2).sum() >= 2
    df = pedigree_generator(max_children=3, FamilyID='FAM1', mode='AR', generation_count=5, seed=3,
                            CousinLikelihood=1.0, min_affected=2)
    from PedGen import gen_PED_export
    path = gen_PED_export(df, str(tmp_path / 'python'))
    assert (tmp_path / 'command' / 'FAM1.ped').read_text() == open(path).read()


def test_no_family_with_enough_affected(tmp_path, capsys):
    assert main(['-g', '1', '--min-affected', '2', '-d', str(tmp_path)]) == 1
    assert 'pedgen: none of 1,000 pedigrees had 2 affected people' in capsys.readouterr().err
    assert list(tmp_path.iterdir()) == []


def test_genotypes(tmp_path):
    table = tmp_path / 'genotypes.tsv'
    main(['-s', '5', '-d', str(tmp_path), '--genotypes', str(table)])
    df = pedigree_generator(max_children=3, FamilyID='FAM1', mode='AD', generation_count=3, seed=5)
    header, row = [line.split('\t') for line in table.read_text().splitlines()]
    assert header[:3] == ['variant', 'chrom', 'af'] and row[:3] == ['causal', 'autosome', '0.001']
    assert header[3:] == [str(i) for i in df['IndividualID']] and row[3:] == [str(g) for g in df['Genotype']]
    # without a frequency there is no af column
    gen_genotype_export(df, tmp_path / 'plain.tsv', variant='v1')
    assert (tmp_path / 'plain.tsv').read_text().startswith('variant\tchrom\t1\t')


def test_pednet_reads_the_genotypes(tmp_path):
    pednet = pytest.importorskip('pednet')
    main(['--mode', 'AR', '-g', '3', '-s', '2', '-d', str(tmp_path), '--genotypes', str(tmp_path / 'g.tsv')])
    family = pednet.read_ped(str(tmp_path / 'FAM1.ped'))['FAM1']
    variants = pednet.read_genotype_table(str(tmp_path / 'g.tsv'))
    score = pednet.rank_variants(family, variants)[0]
    assert score.id == 'causal' and score.log10_bf == score.log10_bf


@pytest.mark.parametrize('family', ['FAM 2', ''])
def test_family_with_a_space(tmp_path, family, capsys):
    assert main(['-f', family, '-d', str(tmp_path)]) == 1
    assert 'must be a word with no spaces' in capsys.readouterr().err
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('args', [['--partners', '1.5'], ['-g', '0'], ['--mode', 'XL'], ['-c', 'two'],
                                  ['--cousins', '2'], ['--min-affected', '0']])
def test_bad_arguments(args):
    with pytest.raises(SystemExit) as stop:
        main(args)
    assert stop.value.code == 2


def test_matplotlib_is_imported_only_to_draw():
    code = 'import sys, PedGen; print("matplotlib" in sys.modules)'
    out = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, check=True).stdout
    assert out.strip() == 'False'
