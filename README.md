# Pedigree Generation Package

PedGen is a Python package for generating genetically realistic pedigrees. Generated pedigrees can be drawn as directed acyclic graphs and exported as PED files.

Pedigree generation and PED file export are also available in a [Streamlit web app](https://pedigree-generatorgit-kzmzvjzpxmt9yyvnryeq7h.streamlit.app).

## Installation

All the Python dependencies are installed with the package:
```bash
pip install PedGen
```

### Drawing with Graphviz (optional)

`plot_pedigree_tree` draws with Graphviz's `dot` layout when PyGraphviz is installed, and otherwise draws one row per generation. PyGraphviz needs Graphviz itself, which has to be installed first:

**Linux**
```bash
sudo apt install graphviz libgraphviz-dev
```

**macOS (Homebrew)**
```bash
brew install graphviz
```

Then install PedGen with PyGraphviz:
```bash
pip install "PedGen[graphviz]"
```

## Usage
The main function, `pedigree_generator`, builds a realistic pedigree by Mendelian inheritance, with alternate allele frequencies in Hardy-Weinberg equilibrium. The pedigree is returned as a pandas DataFrame laid out like a PED file.

### Pedigree Generation:

**`pedigree_generator(max_children, FamilyID, mode, generation_count, SpouseLikelihood = 0.6, AffectedSpouse = True, BackpropLikelihood = 0.25, alt_freq = 0.001, seed = None, CousinLikelihood = 0.0, min_affected = 1)`** returns a DataFrame with the family relationships, sex, phenotype and genotype of everyone in the generated pedigree.

#### **Parameters**
- **max_children** (int): the maximum number of children for any couple in the pedigree
- **FamilyID** (string): the family ID given to the generated pedigree
- **mode** ['AD', 'AR']: the mode of inheritance used to generate the pedigree (AD = autosomal dominant, AR = autosomal recessive)
- **generation_count** (int): the number of generations in the pedigree
- **SpouseLikelihood** (float, 0-1): the chance that a partner is generated for each child as the pedigree is built forwards, so that the child has children of their own (default: 0.6). In every generation but the last, the first child always has a partner, so that the pedigree reaches `generation_count` generations
- **AffectedSpouse** (bool): whether the generated partners can carry the alternate allele (default: True)
- **BackpropLikelihood** (float, 0-1): the chance that a partner's parents are generated as the pedigree is built backwards, and their parents in turn (default: 0.25)
- **alt_freq** (float, 0-1): the frequency of the alternate allele in the population, which partners and their ancestors are drawn from (default: 0.001, a rare allele, as a rare disease's is). It also sets how likely an affected founder is to have two copies in a dominant pedigree: in the share of affected people who do, so almost never for a rare allele
- **seed** (int): makes the pedigree reproducible, so that the same seed and settings always give the same pedigree (default: None, which uses Python's shared random number generator)
- **CousinLikelihood** (float, 0-1): the chance that a descendant of the founder who has a partner marries a cousin rather than someone from outside the family, when a cousin of the other sex has no partner yet (default: 0). The pedigree is built one branch at a time, so the cousins are those in the branches built already
- **min_affected** (int): makes pedigrees until one has at least this many affected people, as families come to notice through their affected members (default: 1). After 1,000 pedigrees without one, it raises a `ValueError` that says what to change

#### A rare recessive condition

With a rare allele, someone who marries into the family almost never carries it, so a recessive condition seldom turns up again in a family unless relatives marry. That is how it is with rare recessive diseases, and why they are commoner in the children of cousins. So for a recessive pedigree with more than the founder affected, give a chance that cousins marry, four or more generations (the founder's grandchildren are the first cousins), and `min_affected`:

```python
df = pedigree_generator(max_children=3, FamilyID="AR1", mode="AR", generation_count=5,
                        CousinLikelihood=1.0, min_affected=2, seed=1)
```

At an allele frequency of 0.001, with cousins marrying whenever they can, about one family in 35 of five generations has a second person affected, and one in 70 of four generations; when they don't, almost none do. So `min_affected` makes 35 or 70 pedigrees, on average, before it finds one, which takes about a second. With a commoner allele, such as 0.02 (about the frequency of the alleles that cause cystic fibrosis in northern Europe, where one person in 25 carries one), a recessive condition also turns up again through partners who carry it.

The returned DataFrame has one row per person, with the columns FamilyID, IndividualID, PaternalID, MaternalID, Sex (1 male, 2 female), Phenotype (1 unaffected, 2 affected) and Genotype (the number of alternate alleles: 0, 1 or 2).

### Additional Functionality:
**`gen_PED_export(df, output_dir='.')`** writes the generated pedigree as a PED file named after its family ID (for example `FAM1.ped`) and returns the file's path.
- **df** (DataFrame): a pedigree from `pedigree_generator`
- **output_dir** (string): the folder to write the PED file in, made if it doesn't exist (default: the working directory)

**`gen_genotype_export(df, path, alt_freq=None, variant='causal')`** writes the Genotype column, each person's number of copies of the alternate allele, as a table that [PedNet](https://github.com/comparativechrono/PedNet) reads as a candidate variant (`pednet rank FAM1.ped genotypes.tsv`), and returns the path.
- **df** (DataFrame): a pedigree from `pedigree_generator`
- **path** (string): the file to write
- **alt_freq** (float, optional): the alternate allele's frequency, as given to `pedigree_generator`, which the table then records
- **variant** (string): the name of the variant in the table (default: `causal`)

**`plot_pedigree_tree(df, title="Pedigree (Tree Layout)")`** draws the generated pedigree as a directed graph, with squares for males and circles for females, and affected people in red.
- **df** (DataFrame): a pedigree from `pedigree_generator`
- **title** (string): the title to show above the drawing

**`construct_pedigree_graph(df)`** builds a `networkx.DiGraph` of the pedigree, with an edge from each parent to each child, for further analysis.

### From the command line

Installing PedGen also installs the `pedgen` command, which generates a pedigree and writes it as a PED file named after the family:

```bash
pedgen --mode AR --generations 4 --family FAM2 --seed 7 -d results     # writes results/FAM2.ped
pedgen --seed 1 --genotypes FAM1-genotypes.tsv                          # FAM1.ped, and each person's copies of the allele
pedgen --mode AR --generations 5 --cousins 1 --min-affected 2 --seed 1  # a recessive condition in children of cousins
```

Its options are `pedigree_generator`'s parameters: `--mode`, `--generations`, `--max-children`, `--family`, `--seed`, `--partners` (SpouseLikelihood), `--ancestors` (BackpropLikelihood), `--allele-frequency` (alt_freq), `--unaffected-partners` (AffectedSpouse = False), `--cousins` (CousinLikelihood) and `--min-affected` (min_affected), with `-d` for the folder to write in and `--genotypes` for the table that `gen_genotype_export` writes. `pedgen --help` lists them, with their defaults.

PedGen imports matplotlib only to draw, so generating and exporting pedigrees doesn't need it, and is quicker to start.

## Tests

```bash
pip install -e ".[test]"
pytest
```

## License

[MIT](https://choosealicense.com/licenses/mit/)
