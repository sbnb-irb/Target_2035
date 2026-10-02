# Target_2035
This repository contains the code necessary to:

 1. Define Pfam-level pockets and obtain its representative. 
 2. Identify ligand-associated pockets. 
 3. Cluster a reduced dataset of pockets.
 4. Select representatives from a cluster following two distinct selection strategies: novel- and experimentally-based.

 as well as data relevant to the updated PocketVec descriptors and the clustering and selection pipeline. 

 The repository is organized in the following directories:

- `Code`: contains the necessary notebooks and functions to apply the presented pipeline. 
- `Results`: contain the papers results. For more information check `Results/README.md`
- `Examples`: contains the data necessary for the examples shown in `Code`.

# Installation

1. Clone this repository to your local folder:

        git clone https://github.com/sbnb-irb/Target_2035.git

2. Create and activate a conda environment with all the requirements:

        conda env create --name target_2035 --file=environment.yml
        conda activate target_2035


