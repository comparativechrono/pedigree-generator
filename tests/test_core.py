import builtins
import json
import os

import matplotlib
matplotlib.use('Agg')

import pytest

from PedGen import construct_pedigree_graph, gen_GA4GH_export, gen_PED_export, pedigree_generator, plot_pedigree_tree
from PedGen.core import affected_founder_weights, generation_layout, generations, inheritance_weights


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


def test_the_allele_is_rare_unless_said_otherwise():
    import inspect
    assert inspect.signature(pedigree_generator).parameters['alt_freq'].default == 0.001


def test_an_affected_founder_is_as_people_with_the_allele_are():
    # a dominant founder is homozygous in the share of affected people who are, q^2 / (q^2 + 2pq): a third at q = 0.5,
    # and one in 2,000 at the default (it used to be a fifth, whatever the frequency)
    het, hom = affected_founder_weights(0.5)
    assert hom / (het + hom) == pytest.approx(1 / 3)
    het, hom = affected_founder_weights(0.001)
    assert hom / (het + hom) == pytest.approx(0.001 / 1.999)
    assert affected_founder_weights(0) == (1, 0) and affected_founder_weights(1) == (0, 1)
    homozygous = sum(pedigree_generator(1, 'F', 'AD', 1, alt_freq=0.5, seed=s).loc[0, 'Genotype'] == 2 for s in range(900))
    assert 250 < homozygous < 350
    assert pedigree_generator(1, 'F', 'AD', 1, alt_freq=0, seed=1).loc[0, 'Genotype'] == 1


# what main made before the allele was rare by default: a recessive family at 0.1 is the same, draw for draw, as no
# random number is drawn for cousins when CousinLikelihood is 0
SEED_3_AT_0_1 = [[1, 0, 0, 1, 2, 2], [2, 0, 0, 2, 1, 0], [3, 1, 2, 2, 1, 1], [4, 0, 0, 1, 1, 0], [5, 4, 3, 2, 1, 0],
                 [6, 7, 8, 1, 1, 0], [7, 0, 0, 1, 1, 0], [8, 0, 0, 2, 1, 0], [9, 6, 5, 1, 1, 0], [10, 6, 5, 2, 1, 0],
                 [11, 6, 5, 1, 1, 0], [12, 1, 2, 1, 1, 1], [13, 0, 0, 2, 1, 0], [14, 12, 13, 2, 1, 0],
                 [15, 0, 0, 1, 1, 0], [16, 15, 14, 2, 1, 0], [17, 15, 14, 1, 1, 0], [18, 12, 13, 1, 1, 1],
                 [19, 0, 0, 2, 1, 1], [20, 18, 19, 2, 1, 0], [21, 18, 19, 2, 1, 0], [22, 18, 19, 1, 1, 1]]


def test_a_recessive_family_at_the_old_frequency_is_as_before():
    df = pedigree_generator(3, 'F', 'AR', 4, alt_freq=0.1, seed=3)
    rows = df[['IndividualID', 'PaternalID', 'MaternalID', 'Sex', 'Phenotype', 'Genotype']].astype(int).values.tolist()
    assert rows == SEED_3_AT_0_1


def ancestors(df, person):
    parents = df.loc[person, ['PaternalID', 'MaternalID']]
    found = set()
    for parent in (int(x) for x in parents if int(x)):
        found |= {parent} | ancestors(df, parent)
    return found


def generation_of(df, person):
    # the founder's generation is 0, and their descendants' is one below their parent's in the family
    if person == 1:
        return 0
    father, mother = int(df.loc[person, 'PaternalID']), int(df.loc[person, 'MaternalID'])
    if not father:
        return None
    above = [g for g in (generation_of(df, father), generation_of(df, mother)) if g is not None]
    return min(above) + 1 if above else None


def test_cousins_marry():
    # with CousinLikelihood=1, a descendant of the founder marries a cousin whenever one is free: a man and a woman of
    # the same generation who share an ancestor, and are not brother and sister. Each has no other partner, and the
    # pedigree stays Mendelian.
    weights = inheritance_weights(0.999, 0.001)
    cousins = set()
    for seed in range(30):
        df = pedigree_generator(3, 'F', 'AR', 5, seed=seed, CousinLikelihood=1.0).set_index('IndividualID')
        partners = {}
        for person, row in df.iterrows():
            father, mother = int(row['PaternalID']), int(row['MaternalID'])
            if not father:
                continue
            assert df.loc[father, 'Sex'] == 1 and df.loc[mother, 'Sex'] == 2
            assert partners.setdefault(father, mother) == mother and partners.setdefault(mother, father) == father
            pair = (int(df.loc[father, 'Genotype']), int(df.loc[mother, 'Genotype']))
            assert int(row['Genotype']) in weights['forward_genotypes'][pair]
            if ancestors(df, father) & ancestors(df, mother):
                cousins.add((seed, father, mother))
                assert generation_of(df, father) == generation_of(df, mother) >= 2
                assert (df.loc[father, 'PaternalID'], df.loc[father, 'MaternalID']) != \
                       (df.loc[mother, 'PaternalID'], df.loc[mother, 'MaternalID'])
    assert len(cousins) >= 5
    # and with CousinLikelihood=0, no one does
    for seed in range(30):
        df = pedigree_generator(3, 'F', 'AR', 5, seed=seed).set_index('IndividualID')
        assert all(not (ancestors(df, int(r['PaternalID'])) & ancestors(df, int(r['MaternalID'])))
                   for _, r in df.iterrows() if int(r['PaternalID']))


def test_min_affected():
    df = pedigree_generator(3, 'F', 'AR', 5, seed=1, CousinLikelihood=1.0, min_affected=2)
    assert (df['Phenotype'] == 2).sum() >= 2
    assert df.equals(pedigree_generator(3, 'F', 'AR', 5, seed=1, CousinLikelihood=1.0, min_affected=2))
    # one generation is only the founder
    with pytest.raises(ValueError, match='none of 1,000 pedigrees had 2 affected people'):
        pedigree_generator(3, 'F', 'AD', 1, seed=1, min_affected=2)


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


def test_ga4gh_export(tmp_path):
    # the GA4GH Pedigree has everyone, with their sexes and parents, and the condition: had by the affected, and
    # marked as excluded for the unaffected, so that the PED file's phenotypes come back from it
    from pedigree_ga4gh import from_ga4gh_pedigree, validate
    df = pedigree_generator(3, 'FAM3', 'AR', 4, seed=12, CousinLikelihood=1, min_affected=2)
    path = gen_GA4GH_export(df, str(tmp_path / 'out'))
    assert path == str(tmp_path / 'out' / 'FAM3_ga4gh-pedigree.json')
    doc = json.load(open(path))
    assert doc['id'] == 'FAM3' and len(doc['individuals']) == len(df)
    assert not [f for f in validate(doc) if f.severity == 'error']
    kindred, notes = from_ga4gh_pedigree(doc)
    assert not notes
    for row in df.itertuples(index=False):
        person = kindred.member(str(row.IndividualID))
        father, mother, _ = kindred.genetic_parents(person.id)
        assert (father or '0', mother or '0') == (str(row.PaternalID), str(row.MaternalID))
        assert person.sex == {1: 'MALE', 2: 'FEMALE'}[row.Sex]
        assert [(t.id, t.label, t.excluded) for t in person.terms] == [('MONDO:0000001', 'disease', row.Phenotype == 1)]
    assert 'autosomal recessive inheritance over 4 generations, allele frequency 0.001, seed 12' in doc['narrative']


def test_ga4gh_export_of_a_phenotype(tmp_path):
    # an HPO term is a phenotype, not a disease
    df = pedigree_generator(2, 'FAM4', 'AD', 3, seed=5)
    doc = json.load(open(gen_GA4GH_export(df, str(tmp_path), condition=('hp:0000365', 'Hearing impairment'))))
    features = [f for i in doc['individuals'] for f in i['phenotypicFeatures']]
    assert len(features) == len(df) and not any('diseases' in i for i in doc['individuals'])
    assert all(f['type'] == {'id': 'HP:0000365', 'label': 'Hearing impairment'} for f in features)
    assert sum(not f.get('excluded') for f in features) == (df['Phenotype'] == 2).sum()


@pytest.mark.parametrize('condition', [('Huntington', 'Huntington disease'), ('XYZ:1', 'unknown ontology'),
                                       ('MONDO:0007739', ' '), 'MONDO:0007739'])
def test_ga4gh_export_refuses_a_condition_it_cannot_write(tmp_path, condition):
    # a GA4GH file can only name a condition by a term of an ontology it can declare, with a label
    df = pedigree_generator(2, 'FAM5', 'AD', 2, seed=1)
    with pytest.raises(ValueError):
        gen_GA4GH_export(df, str(tmp_path), condition=condition)
    assert not list(tmp_path.iterdir())


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
