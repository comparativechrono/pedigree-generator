import builtins
import os

import matplotlib
matplotlib.use('Agg')

import pytest

from PedGen import construct_pedigree_graph, gen_PED_export, pedigree_generator, plot_pedigree_tree
from PedGen.core import generation_layout, generations, inheritance_weights


def hardy_weinberg(q):
    p = 1 - q
    return {0: p * p, 1: 2 * p * q, 2: q * q}


@pytest.mark.parametrize('q', [0.02, 0.1, 0.3, 0.5])
def test_ancestors_keep_hardy_weinberg_proportions(q):
    # A person drawn from the population, then their parents drawn by the reverse weights, must give parents
    # who are themselves in Hardy-Weinberg proportions.
    weights = inheritance_weights(1 - q, q)
    hw = hardy_weinberg(q)
    father = {0: 0.0, 1: 0.0, 2: 0.0}
    mother = {0: 0.0, 1: 0.0, 2: 0.0}
    for child, p_child in hw.items():
        pairs = weights['reverse_genotypes'][child]
        w = weights['reverse_weights'][child]
        for (f, m), wi in zip(pairs, w):
            father[f] += p_child * wi / sum(w)
            mother[m] += p_child * wi / sum(w)
    for g in (0, 1, 2):
        assert father[g] == pytest.approx(hw[g])
        assert mother[g] == pytest.approx(hw[g])


def test_every_way_of_having_a_child_is_listed():
    weights = inheritance_weights(0.9, 0.1)
    # two carriers can have a child with no alternate allele
    assert (1, 1) in weights['reverse_genotypes'][0]
    for child in (0, 1, 2):
        for pair in weights['reverse_genotypes'][child]:
            assert child in weights['forward_genotypes'][pair]


def test_forward_weights_follow_mendel():
    weights = inheritance_weights(0.9, 0.1)
    assert weights['forward_genotypes'][(1, 1)] == [2, 1, 0]
    assert weights['forward_weights'][(1, 1)] == [0.25, 0.5, 0.25]
    assert weights['forward_genotypes'][(2, 0)] == [1]
    assert weights['forward_genotypes'][(0, 0)] == [0]


@pytest.mark.parametrize('mode', ['AD', 'AR'])
def test_generated_pedigrees_are_mendelian(mode):
    weights = inheritance_weights(0.8, 0.2)
    for seed in range(40):
        df = pedigree_generator(max_children=3, FamilyID='F', mode=mode, generation_count=4,
                                BackpropLikelihood=0.5, alt_freq=0.2, seed=seed).set_index('IndividualID')
        assert df.loc[1, 'Phenotype'] == 2
        for person, row in df.iterrows():
            g = int(row['Genotype'])
            affected = g >= 1 if mode == 'AD' else g == 2
            assert row['Phenotype'] == (2 if affected else 1)
            father, mother = int(row['PaternalID']), int(row['MaternalID'])
            assert (father == 0) == (mother == 0)
            if father:
                assert df.loc[father, 'Sex'] == 1 and df.loc[mother, 'Sex'] == 2
                pair = (int(df.loc[father, 'Genotype']), int(df.loc[mother, 'Genotype']))
                assert g in weights['forward_genotypes'][pair]


def test_a_seed_makes_the_pedigree_reproducible():
    a = pedigree_generator(3, 'F', 'AD', 4, seed=7)
    b = pedigree_generator(3, 'F', 'AD', 4, seed=7)
    assert a.equals(b)
    others = [pedigree_generator(3, 'F', 'AD', 4, seed=s) for s in range(8, 14)]
    assert any(not a.equals(o) for o in others)


def test_ped_export(tmp_path):
    df = pedigree_generator(3, 'FAM7', 'AR', 3, seed=1)
    path = gen_PED_export(df, str(tmp_path / 'out'))
    assert os.path.basename(path) == 'FAM7.ped'
    lines = open(path).read().splitlines()
    assert len(lines) == len(df)
    for line in lines:
        cols = line.split('\t')
        assert len(cols) == 6 and cols[0] == 'FAM7'


def test_graph_reads_text_ids():
    import pandas as pd
    df = pd.DataFrame({'FamilyID': ['F'] * 3, 'IndividualID': ['a', 'b', 'c'], 'PaternalID': ['0', '0', 'a'],
                       'MaternalID': ['0', '0', 'b'], 'Sex': [1, 2, 1], 'Phenotype': [2, 1, 2]})
    G = construct_pedigree_graph(df)
    assert set(G.nodes) == {'a', 'b', 'c'}
    assert set(G.edges) == {('a', 'c'), ('b', 'c')}


def test_graph_reads_decimal_and_missing_parents():
    # a parent column with a missing value is read as decimals, so its zeros are 0.0
    import pandas as pd
    df = pd.DataFrame({'FamilyID': ['F'] * 4, 'IndividualID': [1, 2, 3, 4], 'PaternalID': [0.0, 0.0, 1.0, float('nan')],
                       'MaternalID': [0.0, 0.0, 2.0, 0.0], 'Sex': [1, 2, 1, 2], 'Phenotype': [2, 1, 2, 1]})
    G = construct_pedigree_graph(df)
    assert set(G.nodes) == {1, 2, 3, 4}
    assert set(G.edges) == {(1, 3), (2, 3)}


def test_generation_layout_puts_partners_in_one_row():
    # with the ancestors of people who married in drawn too
    for seed in range(30):
        df = pedigree_generator(3, 'F', 'AD', 4, BackpropLikelihood=0.5, seed=seed)
        G = construct_pedigree_graph(df)
        pos = generation_layout(G)
        for child in G.nodes:
            parents = list(G.predecessors(child))
            if len(parents) == 2:
                assert pos[parents[0]][1] == pytest.approx(pos[parents[1]][1])
                assert pos[child][1] == pytest.approx(pos[parents[0]][1] - 1)


def test_generations_when_people_of_two_generations_have_children():
    # 4 has a child with his niece 6: whatever the order of the rows, each child is one generation below its lower parent
    import pandas as pd
    rows = [(1, 0, 0, 1), (2, 0, 0, 2), (3, 1, 2, 1), (4, 1, 2, 1), (5, 0, 0, 2), (6, 3, 5, 2), (7, 4, 6, 1)]
    for order in (rows, rows[::-1], rows[4:] + rows[:4]):
        df = pd.DataFrame([('F', i, f, m, sex, 1) for i, f, m, sex in order],
                          columns=['FamilyID', 'IndividualID', 'PaternalID', 'MaternalID', 'Sex', 'Phenotype'])
        assert generations(construct_pedigree_graph(df)) == {1: 0, 2: 0, 3: 1, 4: 1, 5: 1, 6: 2, 7: 3}


def test_plot_without_pygraphviz(monkeypatch):
    real_import = builtins.__import__

    def no_pygraphviz(name, *args, **kwargs):
        if name in ('pygraphviz', 'networkx.drawing.nx_agraph'):
            raise ImportError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', no_pygraphviz)
    df = pedigree_generator(3, 'F', 'AD', 3, seed=2)
    fig = plot_pedigree_tree(df)
    assert fig is not None
