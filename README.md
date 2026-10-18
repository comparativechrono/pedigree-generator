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

**`plot_pedigree_tree(df, title="Pedigree (Tree Layout)")`** draws the generated pedigree as a directed graph, with squares for males and circles for females, and affected people in red.
- **df** (DataFrame): a pedigree from `pedigree_generator`
- **title** (string): the title to show above the drawing

**`construct_pedigree_graph(df)`** builds a `networkx.DiGraph` of the pedigree, with an edge from each parent to each child, for further analysis.

## Tests

```bash
pip install -e ".[test]"
pytest
```

## License

[MIT](https://choosealicense.com/licenses/mit/)
