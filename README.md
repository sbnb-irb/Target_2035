# Target_2035
**PocketVec** is a strategy based on the concept of *similar pockets bind similar molecules*, that allows for the generation of 128-bit binding site descriptors in which each bit encodes the docking rank of a particular molecule. It was shown that similar pockets had similar rankings for the same sets of molecules, leading to more similar descriptors than pockets that are not similar. This makes PocketVec a powerful tool for assessing pocket similarity across the whole human pocketome. 

For further information on the PocketVec methodology, we recommend reading the original publication:

        Comajuncosa-Creus, A., Jorba, G., Barril, X., & Aloy, P. (2024). Comprehensive detection and characterization of human druggable pockets through binding site descriptors. Nature Communications, 15(1), 7917..

Here, we first present an update of the PocketVec descriptors, which approximately doubles the number of described pockets acrosss all four datasets. We then present a novel strategy that leverages PocketVec descriptors to systematically cluster pockets with similar ligand-binding profiles, first at the Pfam domain level and subsequently across different Pfam domains. From the generated clusters, we select 200 representative pockets using two alternative prioritization strategies designed to identify pockets with distinct ligand-binding profiles while maximizing coverage of the generated pocket space. 

Further methodological and conceptual details are provided in the associated publication:

        **Add publication**

This repository contains the updated PocketVec descriptors and the results presented in the manuscript. It also contains the required code to define Pfam-level pockets and select their representatives, identify ligand-associated pockets, and perform the pocket clustering and representative-selection procedures.

# Repository structure

```
Target_2035/
├── Code/
│   ├── Notebooks/
│   │   ├── 1_Clustering.ipynb
│   │   └── 2_Selection.ipynb
│   └──utils/
│       └── utils.py
├── Examples/
├── Results/
├── environment.yml
└── README.md
```

- `Code/` contains the notebooks and functions implementing the clustering and selection of pockets.
- `Examples/` contains the data required to run the examples provided in the notebooks.
- `Results/` contains the results generated for the publication. See `Results/README.md` for further information.
- `environment.yml` contains the Conda environment with the required dependencies to run the provided code.

# Installation

1. Clone this repository to your local folder:

        git clone https://github.com/sbnb-irb/Target_2035.git

2. Create and activate a conda environment with all the requirements:

        conda env create --name target_2035 --file=environment.yml
        conda activate target_2035



# Provided code
The analysis is divided into two notebooks:

- `1_Clustering.ipynb`: exemplifies, for a selected Pfam domain, how to define its Pfam-level pockets, obtain the pocket representative and identify if the pockets are associated to ligands. It also demonstartes the clustering of the representative pcokets. 
- `2_Selection.ipynb`: exemplifies, for a given cluster, how to select the representative pockets according to the two selection strategies. 

Both notebooks use functions defined in `utils.py` and data provided in the `Examples/` and `Results/` directories.s