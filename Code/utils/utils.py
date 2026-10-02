from openbabel import pybel
import random 
from pathlib import Path
import subprocess 
import re
import pandas as pd
import numpy as np
from Bio.PDB import *
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import cdist
from scipy.spatial.distance import mahalanobis
import shutil

def sd_to_pdb(input_file, output_file):
    """
    Convert a centroid SD file to a PDB.

    Args:
        input_file (str | Path): Path to input structure file (SD format)
        output_file (str | Path): Path to output structure file (PDB format)
    """
    try:
        mol = next(pybel.readfile("sdf", str(input_file)))
        mol.write("pdb", str(output_file), overwrite=True)
    except Exception as e:
        print(f"[FAILED] {input_file}: {e}")

def convert_pockets(structure_path, domains):
    """
    Convert all pocket .sd files of each domain to PDB.

    Args:
        structure_path (str | Path): Folder of the Pfam family (<outpath>/<pfam>)
        domains (list): Names of the domain subfolders to process
    """
    for domain in domains:
        print(domain)
        for pocket in sorted((Path(structure_path) / domain).glob("*.sd")):
            sd_to_pdb(pocket, pocket.with_suffix(".pdb"))

def select_reference(domains, random_seed = 42):
    """
    From a list of possible domains, select the reference structure.
    
    Args:
        domains (list): List of the different domains to choose from. 
        random_seed (int): Seed to use for the random selection. 
    """

    random.seed(random_seed)
    ref_dom = random.sample(domains, 1)[0]
    return ref_dom

def find_tmalign():
    """
    Locates the TMalign executable installed in the active conda environment.
    """

    for name in ("TMalign", "tmalign"):
        found = shutil.which(name)
        if found:
            return found

    bin_dir = Path(sys.executable).parent
    for name in ("TMalign", "tmalign"):
        if (bin_dir / name).is_file():
            return str(bin_dir / name)

    raise FileNotFoundError(
        "TMalign not found. Install it with `conda install -c bioconda tmalign` "
        "(it is included in environment.yml) and make sure the notebook uses that environment."
    )

def run_TM_align(path_to_reference, path_to_structure, matrix_file, tmalign=None):
    """
    Align a structure onto the reference with TMalign, writting the matrix file and returning the RMSD.

    Args: 
        path_to_reference (str | Path) : Path of the reference structure. 
        path_to_structure (str | Path) : Path to the structure we want to superimpose. 
        matrix_file (str | Path) : Path to where we want to write the superposition matrix.
        tmalign (str | Path, optional): Path to a TMalign executable. By default, the one from the active environment is used.
    
    """
    tmalign = str(tmalign) if tmalign else find_tmalign()

    result = subprocess.run(
        [tmalign, str(path_to_structure), str(path_to_reference), "-m", str(matrix_file)],
        capture_output=True, text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(result.stderr)

    match = re.search(r"RMSD=\s*([\d.]+)", result.stdout)

    if match is None:
        raise RuntimeError(f"Could not find RMSD in TMalign output for {path_to_structure}")
    return float(match.group(1))

def read_mtx(file):
    """
    Read a TMalign matrix file. Returns translation t (3,) and rotation u (3, 3).

    Args: 
        file (str | Path) : Path of the superposition matrix. 

    """
    f = open(file, "r").readlines()[2:5]
    f = np.array([i.split()[1:] for i in f], dtype='float32')
    t = np.array(f[:,0])
    u = np.array([f[0][1:], f[1][1:], f[2][1:]])
    return t, u

def superpose_file(input_file, output_file, matrix_file):
    """
    Apply a TMalign superposition (translation t, rotation u) to all atoms of a PDB file.

    Args:
        input_file (str | Path) : Path of the structure to be superimposed. 
        output_file (str | Path) : Path to save the superposed structure. 
        matrix_file (str | Path) : Path to save the superposed structure. 

    """

    t, u = read_mtx(matrix_file)
    structure = PDBParser(QUIET=True).get_structure("st", str(input_file))
    
    for atom in structure.get_atoms():
        x0, y0, z0 = atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]
        x = t[0] + u[0][0]*x0 + u[0][1]*y0 + u[0][2]*z0
        y = t[1] + u[1][0]*x0 + u[1][1]*y0 + u[1][2]*z0
        z = t[2] + u[2][0]*x0 + u[2][1]*y0 + u[2][2]*z0
        atom.coord = np.array([x, y, z])
        
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_file))   

def superpose_structures(pfam_path, domains, reference):
    """
    Superpose all structures and centroids of each domain onto the reference structure.
    
    Args:
        pfam_path (str | Path): Folder of the Pfam family (<outpath>/<pfam>)
        domains (list): List of the different domains. 
        reference (str): Reference structure. 
    """

    path_to_reference = pfam_path / reference / f"{reference}.pdb"
    rmsd_values = {}

    for domain in domains:
        domain_path = pfam_path / domain
        matrix_file = domain_path / "alignment_matrix.txt"

        try:
            rmsd_values[domain] = run_TM_align(
                path_to_reference, domain_path / f"{domain}.pdb", matrix_file
            )

        except RuntimeError as e:
            print(f"[FAILED] {domain}: {e}")
            rmsd_values[domain] = None
            continue

        for pdb in sorted(domain_path.glob("*.pdb")):
            if pdb.stem.endswith("_aligned"):
                continue
            superpose_file(pdb, pdb.with_name(f"{pdb.stem}_aligned.pdb"), matrix_file)

    pd.Series(rmsd_values, name="RMSD").rename_axis("domain").to_csv(pfam_path / "RMSD.tsv", sep="\t")

def read_centroid(pdb_file):
    """
    Return the (x, y, z) coordinates of the single atom in a centroid PDB file.
    
    Args:
        pdb_file (str | Path): Folder of the Pfam family (<outpath>/<pfam>)

    """
    
    structure = PDBParser(QUIET=True).get_structure("st", str(pdb_file))
    atoms = list(structure.get_atoms())
    if len(atoms) != 1:
        raise ValueError(f"{pdb_file} has {len(atoms)} atoms, expected 1")
    return atoms[0].coord

def cluster_centroids(pfam_path, metric="euclidean", distance = 10, linkage = "complete"):
    """
    Cluster the superposed pocket centroids of a Pfam family by their coordinates.

    Args:
        pfam_path (str | Path): Folder of the Pfam family (<outpath>/<pfam>)
        metric (str): Distance metric between centroids
        distance (float): Distance threshold (Å) above which clusters are not merged
        linkage (str): Linkage criterion ("complete", "average", "single", "ward")

    """

    pfam = pfam_path.name
    coordinates = []
    
    for pocket in sorted(pfam_path.glob("*/*fpocket*_aligned.pdb")):
        if pocket.stem.endswith("_aligned"):
            x, y, z = read_centroid(pocket)
            pocket_name = pocket.stem.removesuffix("_aligned")
            coordinates.append([pocket_name, x, y, z])

    if not coordinates:
        raise ValueError(f"No aligned centroids found in {pfam_path}")

    coordinates = pd.DataFrame(coordinates, columns=["Descriptor", "X", "Y", "Z"])

    if len(coordinates) == 1:  
        coordinates["Pocket"] = 0
    else:
        clustering = AgglomerativeClustering(
            metric=metric, distance_threshold=distance, n_clusters=None, linkage=linkage
        )
        coordinates["Pocket"] = clustering.fit_predict(coordinates[["X", "Y", "Z"]])

    out = (
        coordinates.groupby("Pocket")["Descriptor"].agg(",".join)
        .rename("Included pockets")
        .reset_index()
    )
    out.insert(0, "Pfam domain", pfam)
    out.to_csv(pfam_path / f"{pfam}_pockets.tsv", index=False, sep="\t")

def find_representative(pockets, descriptors):
    """
    Compute the cosine distance between pockets of the same cluster and select the one with the lowest cosine distance. 

    Args:
        pockets (list): Names of the pockets in the cluster
        descriptors (dict): Pocket : PocketVec descriptor

    """

    if len(pockets) == 1:
        return [pockets[0],  np.nan, np.nan]

    X = np.array([descriptors[p] for p in pockets])
    D = cdist(X, X, metric="cosine")
    mean_dist = D.sum(axis=1) / (len(pockets) - 1) 
    best = mean_dist.argmin()

    return [pockets[best], round(mean_dist.mean(), 3), round(mean_dist[best], 3)]
    
def cluster_centroids_lig(pfam_path, metric="euclidean", distance = 10, linkage = "complete"):
    """
    Cluster the superposed pocket centroids of a Pfam domain by their coordinates.

    Args:
        pfam_path (str | Path): Folder of the Pfam family (<outpath>/<pfam>)
        metric (str): Distance metric between centroids
        distance (float): Distance threshold (Å) above which clusters are not merged
        linkage (str): Linkage criterion ("complete", "average", "single", "ward")
    """

    pfam = pfam_path.name
    coordinates = []
        
    for pocket in sorted(pfam_path.glob("*CTR*.pdb")):
        x, y, z = read_centroid(pocket)
        pocket_name = pocket.stem
        coordinates.append([pocket_name, "LIG", x, y, z])

    for pocket in sorted(pfam_path.glob("*fpocket*.pdb")):
        x, y, z = read_centroid(pocket)
        pocket_name = pocket.stem
        coordinates.append([pocket_name, "PD", x, y, z])

    
    if not coordinates:
        raise ValueError(f"No aligned centroids found in {pfam_path}")

    coordinates = pd.DataFrame(coordinates, columns=["Pocket", "Method", "X", "Y", "Z"])

    clustering = AgglomerativeClustering(
        metric=metric, distance_threshold=distance, n_clusters=None, linkage=linkage
    )
    coordinates["Cluster"] = clustering.fit_predict(coordinates[["X", "Y", "Z"]])

    lig_clusters = set(coordinates.loc[coordinates["Method"] == "LIG", "Cluster"])

    coordinates["Ligand-associated"] = (
        coordinates["Cluster"].isin(lig_clusters)
        .astype("boolean")                           
        .where(coordinates["Method"] == "PD", pd.NA) 
    )
    
    coordinates.to_csv(pfam_path / f"{pfam}_ligand.tsv", index=False, sep="\t")

def cluster_pockets(pockets, descriptors, outfile, metric="cosine", distance = 0.22, linkage="complete"):
    """
    Cluster PocketVec descriptors.

    Args:
        pockets (list): Names of the pockets to cluster. 
        descriptors (dict): Pocket : PocketVec descriptor
        outfile (str): File to save the results. 
        metric (str): Distance metric between centroids
        distance (float): Distance threshold above which clusters are not merged
        linkage (str): Linkage criterion ("complete", "average", "single", "ward")
    """

    reduced_pockets = [i for i in descriptors.keys() if i in pockets]
    reduced_descriptors = [descriptors[i] for i in reduced_pockets]
    labels = AgglomerativeClustering(metric=metric, distance_threshold = distance, n_clusters=None, linkage=linkage).fit_predict(reduced_descriptors)
    clusters = [(i,j) for i,j in zip(reduced_pockets, labels)]
    clusters_df = pd.DataFrame(clusters, columns= ["Pocket", "Cluster"]).sort_values(by="Pocket")
    clusters_df.to_csv(outfile, sep="\t", index= False)
    
def compute_buriedness_threshold(lig_metrics, percentage=99):
    """
    Given a reference sample, compute the value for which the required percentage is included.

    Args:
        lig_metrics (Series): metrics of the reference sample. 
        percentage (int): percentage of the included sample. 
    """
    lig = np.array(lig_metrics)
    threshold = np.percentile(lig, 100 - percentage)
    return threshold

def apply_novelty_filter(pockets, metrics, threshold):
    """
    Given a set of pockets, filter out those that do not follow the novelty startegy's requirements. 
    Args:
        pockets (Dataframe): pockets to filter. 
        metrics (Dataframe): Dataframe containing relevant information for the filtering. 
        threshold (int): pre-computed buriedness threshold. 
    """

    pockets = pockets.copy()
    pockets["Best docking score"] = pockets["Pocket"].map(
        metrics.set_index("Pocket")["Best docking score"]
    )

    pockets["Buriedness"] = pockets["Pocket"].map(
        metrics.set_index("Pocket")["Buriedness"]
    )
    pockets["Prank score"] = pockets["Pocket"].map(
            metrics.set_index("Pocket")["Prank score"]
        )
    mask = (
        ~pockets["Pfam domain has ligand"]
        & ~pockets["Protein has ligand"]
        & ~pockets["Protein has PDB"]
        & (pockets["Buriedness"] >= threshold)
    )
    return pockets[mask]
    
def compute_mahalanobis(pockets, lig_data, metrics = ["Best docking score" ,"Buriedness", "Prank score"]):
    """
    Given a set of pockets and a reference sample, compute the Mahalanobis distance. 
    Args:
        pockets (Dataframe): pockets we want to compute the Mahalanobis distance for. 
        lig_data (Dataframe): reference sample. 
        metrics (list): metrics used to compute the Mahalanobis distance. 
    """

    ref_data = lig_data[metrics].values
    pocket_data = pockets[metrics].values
    
    mu = ref_data.mean(axis=0)
    sigma = np.cov(ref_data.T)

    sigma_inv = np.linalg.inv(sigma)
        
    distances = np.array([
        mahalanobis(x, mu, sigma_inv) 
        for x in pocket_data
    ])

    pockets["Mahalanobis"] = [i for i in distances]
    
def choose_novel_representatives(filtered_pockets, outfile, num_select = 2):
    """
    Given a set of pockets, select the required number according to the novel selection strategy. 
    Args:
        filtered_pockets (Dataframe): pockets from which we want to select. 
        outfile (str): file to save the selection results. 
        num_select (int): number of pockets to be selected. 
    """

    pockets = filtered_pockets.copy()
    pockets["Rank"] = pockets["Mahalanobis"].rank(ascending=True)
    cluster = re.fullmatch(r"cluster_(\d+)_novel_selection", outfile.stem).group(1)

    if len(pockets) <= num_select:
        print("Not enough pockets")
        return

    to_write = []
    for rank in range(num_select):
        selection = pockets.loc[pockets["Rank"] == rank + 1, "Pocket"].to_list()[0]
        mahalanobis_val = pockets.loc[pockets["Pocket"] == selection, "Mahalanobis"].to_list()[0]
        to_write.append([cluster, selection, round(mahalanobis_val, 4), rank + 1])

    to_write_df = pd.DataFrame(to_write, columns=["Cluster", "Pocket", "Mahalanobis distance", "Rank"])
    to_write_df.to_csv(outfile, sep="\t", index=False)

def apply_experimental_filter(pockets, metrics, threshold):
    """
    Given a set of pockets, filter out those that do not follow the experimental startegy's requirements. 
    Args:
        pockets (Dataframe): pockets to filter. 
        metrics (Dataframe): Dataframe containing relevant information for the filtering. 
        threshold (int): pre-computed buriedness threshold. 
    """

    pockets = pockets.copy()
    pockets["Best docking score"] = pockets["Pocket"].map(
        metrics.set_index("Pocket")["Best docking score"]
    )

    pockets["Buriedness"] = pockets["Pocket"].map(
        metrics.set_index("Pocket")["Buriedness"]
    )
    pockets["Prank score"] = pockets["Pocket"].map(
            metrics.set_index("Pocket")["Prank score"]
        )
    mask = (
        ~pockets["Pfam domain has ligand"]
        & ~pockets["Protein has ligand"]
        & (pockets["Buriedness"] >= threshold)
    )
    return pockets[mask]

def select_pocket(data, num_select, to_write, type_data, cluster):
    """
    Given a subset of pockets, select the required number of them according to their Mahalanobis distance. 
    Args:
        data (Dataframe): pockets from which we want to select. 
        num_select (int): number of pockets to be selected. 
        to_write(list): contains the already selected pockets. 
        type_data (str): origin of the pocket (SGC, PDB or AF2).
        cluster (str): cluster number. 
    """

    to_select = range(0, min(len(data), num_select))
    num_select -= min(len(data), num_select)

    for rank in to_select:
        selection = data.loc[data["Rank"] == rank + 1, "Pocket"].to_list()[0]
        mahalanobis_val = data.loc[data["Pocket"] == selection, "Mahalanobis"].to_list()[0]
        to_write.append([cluster, selection, round(mahalanobis_val, 4), len(to_write) + 1, type_data])

    return num_select

def choose_experimental_representatives(filtered_pockets, outfile, num_select = 2):
    """
    Given a set of pockets, select the required number according to the experimental selection strategy. 
    Args:
        filtered_pockets (Dataframe): pockets from which we want to select. 
        outfile (str): file to save the selection results. 
        num_select (int): number of pockets to be selected. 
    """

    pockets = filtered_pockets.copy()
    pockets["Rank"] = pockets["Mahalanobis"].rank(ascending=True)
    cluster = re.fullmatch(r"cluster_(\d+)_experimental_selection", outfile.stem).group(1)
    
    if len(pockets) <= num_select:
        print("Not enough pockets")
        return
    
    sgc_filtered = pockets[pockets["Domain has SGC"] == True].copy()
    to_write = []
    if len(sgc_filtered) > 0:
        sgc_filtered["Rank"] = sgc_filtered["Mahalanobis"].rank(ascending=True)
        num_select = select_pocket(sgc_filtered, num_select, to_write, "SGC", cluster)

    if num_select != 0:
        pdb_filtered = pockets[(pockets["Domain has PDB"] == True) & (pockets["Domain has SGC"] == False)].copy()
        if len(pdb_filtered) > 0:
            pdb_filtered["Rank"] = pdb_filtered["Mahalanobis"].rank(ascending=True)
            num_select = select_pocket(pdb_filtered, num_select, to_write, "PDB", cluster)
                        
    if num_select != 0:
        af2_filtered = pockets[(pockets["Protein has PDB"] == False) & (pockets["Domain has SGC"] == False)].copy()
        if len(af2_filtered) > 0:
            af2_filtered["Rank"] = af2_filtered["Mahalanobis"].rank(ascending=True)
            num_select = select_pocket(af2_filtered, num_select, to_write, "AF2", cluster)
                    
    to_write_df = pd.DataFrame(to_write, columns=["Cluster", "Pocket", "Mahalanobis distance", "Rank", "Origin"])
    to_write_df.to_csv(outfile, sep="\t", index=False)




















