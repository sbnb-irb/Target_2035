# Results README

This repository contains the data generated and used in the analysis described in the manuscript. The files are organized into four distinct directories: **PocketVec**, **Descriptors**, **ProSPECCTs**, and **Selection**.


## 1. PocketVec

The `PocketVec` directory contains the protein-, domain-, and pocket-level data used to generate the updated PocketVec datasets.

| File | Description |
|---|---|
| `Human_Proteome_Uniprot.tsv` | For each reviewed human protein downloaded from UniProt: UniProt ID, entry name, protein name, gene name, organism, organism ID, database, review status and protein length. |
| `Uniprot_to_Pfam.tsv` | For each reviewed human protein with a Pfam domain or family: number of Pfam domains/families and the number of repeats of each Pfam domain/family. |
| `PDB_LIG_domain.tsv` | For each domain (defined by protein, Pfam domain, and start/end positions): reference PDB structure, structure resolution, number of structures, number of identified pockets, number of aligned ligands, and pocket information. |
| `PDB_PD_domain.tsv` | For each domain (defined by protein, Pfam domain, and start/end positions): detected pocket identifier, reference PDB structure, druggability score, Prank rescore, and pocket buriedness. |
| `AF2_LIG_domain.tsv` | For each domain (defined by protein, Pfam domain, and start/end positions): PDB-LIG centroid used to define the pocket, pocket buriedness, and minimum pLDDT. |
| `AF2_PD_domain.tsv` | For each domain (defined by protein, Pfam domain, and start/end positions): detected pocket identifier, druggability score, Prank rescore, pocket buriedness, and minimum pLDDT. |
| `PDB_LIG_pocket.tsv` | For each pocket located in a domain: reference PDB structure, pocket identifier, and pocket-centroid coordinates (`X`, `Y`, `Z`). |
| `PDB_PD_pocket.tsv` | For each pocket located in a domain: reference PDB structure, druggability score, Prank rescore, pocket buriedness, and pocket-centroid coordinates (`X`, `Y`, `Z`). |
| `AF2_LIG_pocket.tsv` | For each pocket located in a domain: pocket buriedness, minimum pLDDT, and pocket-centroid coordinates (`X`, `Y`, `Z`). |
| `AF2_PD_pocket.tsv` | For each pocket located in a domain: druggability score, Prank rescore, pocket buriedness, minimum pLDDT, and pocket-centroid coordinates (`X`, `Y`, `Z`). |

## 2. Descriptors

The `Descriptors` directory contains the updated PocketVec descriptors generated for each dataset (PDB-LIG, PDB-PD, AF2-LIG and AF2-PD).

| File | Description |
|---|---|
| `descriptors_PDB_LIG_80.pkl` | PocketVec descriptors for the PDB-LIG dataset. The dictionary uses pocket identifiers as keys and the corresponding PocketVec descriptors as values. Outlier descriptors are removed. |
| `descriptors_PDB_PD_80.pkl` | PocketVec descriptors for the PDB-PD dataset. The dictionary uses pocket identifiers as keys and the corresponding PocketVec descriptors as values. Outlier descriptors are removed. |
| `descriptors_AF2_LIG_80.pkl` | PocketVec descriptors for the AF2-LIG dataset. The dictionary uses pocket identifiers as keys and the corresponding PocketVec descriptors as values. Outlier descriptors are removed. |
| `descriptors_AF2_PD_80.pkl` | PocketVec descriptors for the AF2-PD dataset. The dictionary uses pocket identifiers as keys and the corresponding PocketVec descriptors as values. Outlier descriptors are removed. |

## 3. ProSPECCTs

The `ProSPECCTs` directory contains the descriptors generated for the different ProSPECCTs datasets.

| File | Description |
|---|---|
| `D1.pkl` | PocketVec descriptors for ProSPECCTs dataset D1. |
| `D1.2.pkl` | PocketVec descriptors for ProSPECCTs dataset D1.2. |
| `D2.pkl` | PocketVec descriptors for ProSPECCTs dataset D2. |
| `D3.pkl` | PocketVec descriptors for ProSPECCTs dataset D3. |
| `D4.pkl` | PocketVec descriptors for ProSPECCTs dataset D4. |
| `D5.pkl` | PocketVec descriptors for ProSPECCTs dataset D5. |
| `D5.2.pkl` | PocketVec descriptors for ProSPECCTs dataset D5.2. |
| `D7.pkl` | PocketVec descriptors for ProSPECCTs dataset D7. |
| `AF2_D1.pkl` | PocketVec descriptors generated for ProSPECCTs dataset D1 using AlphaFold structures. |

## 4. Selection

The `Selection` directory contains the files used for pocket clustering and selection, including the defined Pfam-level pockets, ligand-associated pockets, reduced and full clustered datasets, as well as the the final selected representatives for both selection strategies

| File | Description |
|---|---|
| `PDB_LIG_Pfam_pockets.tsv` | For each identified Pfam-level pocket in the PDB-LIG dataset (defined by Pfam domain and pocket identifier): list of the included pockets. |
| `AF2_PD_Pfam_pockets.tsv` | For each identified Pfam-level pocket in the AF2-PD dataset (defined by Pfam domain and pocket identifier): list of the included pockets. |
| `representative_pockets.tsv` | For each identified Pfam-level pocket in the AF2-PD dataset (defined by Pfam domain and pocket identifier): list of included pockets, chosen representative pocket, overall mean cosine distance, representative mean cosine distance, and cluster assigned to the representative. Cosine-distance values are generated only when more than one pocket is included. |
| `AF2_PD_ligand_information.tsv` | For each pocket in the AF2-PD dataset: whether the pocket is ligand-associated (`True`/`False`), best docking score (i.e. lowest), pocket buriedness, and Prank score. |
| `full_dataset.tsv` | For each pocket in the AF2-PD dataset: assigned cluster, whether the corresponding Pfam-level pocket, Pfam domain, or protein has a ligand; whether the protein or domain has a structure in the PDB; and whether the domain is being structurally determined by the SGC. |
| `novel_selection.tsv` | For each selected pocket obtained using the novel selection strategy: cluster it belongs to, cluster rank, and Mahalanobis distance. |
| `experimental_selection.tsv` | For each selected pocket obtained using the experimental selection strategy: cluster it belongs to, cluster rank, Mahalanobis distance, and structure origin (`SGC`, `PDB`, or `AF2`). |

## 5. Abbreviations

**AF2**, AlphaFold 2; **PDB**, Protein Data Bank; **Pfam**, protein family database; **SGC**, Structural Genomics Consortium; **pLDDT**, predicted Local Distance Difference Test.
