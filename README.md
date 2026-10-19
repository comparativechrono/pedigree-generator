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

**`pedigree_generator(max_children, FamilyID, mode, generation_count, SpouseLikelihood = 0.6, AffectedSpouse = True, BackpropLikelihood = 0.25, alt_freq = 0.1, seed = None)`** returns a DataFrame with the family relationships, sex, phenotype and genotype of everyone in the generated pedigree.

#### **Parameters**
- **max_children** (int): the maximum number of children for any couple in the pedigree
- **FamilyID** (string): the family ID given to the generated pedigree
- **mode** ['AD', 'AR']: the mode of inheritance used to generate the pedigree (AD = autosomal dominant, AR = autosomal recessive)
- **generation_count** (int): the number of generations in the pedigree
- **SpouseLikelihood** (float, 0-1): the chance that a partner is generated for each child as the pedigree is built forwards, so that the child has children of their own (default: 0.6). The first child always has a partner, so that the pedigree reaches `generation_count` generations
- **AffectedSpouse** (bool): whether the generated partners can carry the alternate allele (default: True)
- **BackpropLikelihood** (float, 0-1): the chance that a partner's parents are generated as the pedigree is built backwards, and their parents in turn (default: 0.25)
- **alt_freq** (float, 0-1): the frequency of the alternate allele in the population, which partners and their ancestors are drawn from (default: 0.1)
- **seed** (int): makes the pedigree reproducible, so that the same seed and settings always give the same pedigree (default: None, which uses Python's shared random number generator)

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
```

Its options are `pedigree_generator`'s parameters: `--mode`, `--generations`, `--max-children`, `--family`, `--seed`, `--partners` (SpouseLikelihood), `--ancestors` (BackpropLikelihood), `--allele-frequency` (alt_freq) and `--unaffected-partners` (AffectedSpouse = False), with `-d` for the folder to write in and `--genotypes` for the table that `gen_genotype_export` writes. `pedgen --help` lists them, with their defaults.

PedGen imports matplotlib only to draw, so generating and exporting pedigrees doesn't need it, and is quicker to start.

## Tests

```bash
pip install -e ".[test]"
pytest
```

## License

[MIT](https://choosealicense.com/licenses/mit/)
