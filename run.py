import os
import sys
import argparse
import json
from pathlib import Path

import pandas as pd
import scanpy as sc

#repo_dir = Path(__file__).parent
#sys.path.insert(0, str(repo_dir))



# Parse command line arguments
parser = argparse.ArgumentParser(description="Scanpy single-cell analysis pipeline")
parser.add_argument(
    "--output_dir", "-o", type=str, required=True, help="Output directory"
)
parser.add_argument("--name", "-n", type=str, required=True, help="Dataset name")
parser.add_argument(
    "--data.h5ad", dest="data_h5ad", type=str, required=True, help="Input h5ad file"
)
parser.add_argument(
    "--data.clusters_truth_num",
    dest="clusters_truth_num_path",
    type=str,
    required=True,
    help="input true number of clusters path",
)
parser.add_argument(
    "--method_name", type=str, choices=["random"], help="Method to run"
)
parser.add_argument(
    "--d_cluster",
    type=int,
    help="delta number of clusters with respect to true number of clusters",
)
parser.add_argument(
    "--n_comp",
    type=int,
    default=50,
    help="number of PCA components to use for KNN graph construction",
)
parser.add_argument(
    "--n_neig",
    type=int,
    default=15,
    help="number of neighbors to use for KNN graph construction",
)
parser.add_argument(
    "--n_hvg",
    type=int,
    default=1000,
    help="number of highly variable genes to use",
)

parser.add_argument("--seed", type=int, default=28, help="seed number")

# only have manual filtering
parser.add_argument(
    "--filter",
    type=str,
    choices=["manual"],
    default="manual",
    help="filtering strategy (manual uses suggested cutoffs; auto uses package QC)",
)
args, _ = parser.parse_known_args()

# time object to store time involved (in seconds) in each step
timings: dict[str, float | None] = {
    "gpu_load": None,
    "find_mit_gene": None,
    "filter": None,
    "normalization": None,
    "hvg": None,
    "scaling": None,
    "pca": None,
    "t_sne": None,
    "umap": None,
    "louvain": None,
    "leiden": None,
}

clustering_info: dict[str, dict[str, float | int | None]] = {
    "resolutions": {
        "louvain": None,
        "leiden": None,
    },
    "num_runs": {
        "louvain": None,
        "leiden": None,
    },
}

sc.settings.verbosity = 3

# data ####
adata = sc.read_h5ad(args.data_h5ad)
adata.var_names_make_unique()
with open(args.clusters_truth_num_path, "r") as f:
    clusters_truth_num = int(f.read())
n_cluster = clusters_truth_num + args.d_cluster

rng = np.random.default_rng(arg.seed)
louvain = rng.integers(0, n_clusters, size=adata.n_obs)
leiden = rng.integers(0, n_clusters, size=adata.n_obs)

# Save timings as JSON
with open(os.path.join(args.output_dir, f"{args.name}.timings.json"), "w") as f:
    json.dump(timings, f, indent=2)

# Save clustering metadata as JSON
with open(os.path.join(args.output_dir, f"{args.name}.clustering_info.json"), "w") as f:
    json.dump(clustering_info, f, indent=2)

# Save cluster assignments as TSV
adata.obs[["louvain", "leiden"]].reset_index(names="cell_id").to_csv(
    os.path.join(args.output_dir, f"{args.name}.clusters.tsv"),
    sep="\t",
    index=False,
)

hvg_list = adata.var.names()
hvg_list.to_csv(
    os.path.join(args.output_dir, f"{args.name}.hvgs.tsv"),
    sep="\t",
    index=False,
    header=False,
)

# Save PCA coordinates as TSV
pca_df = pd.DataFrame(
    adata.obsm["X_pca"],
    index=adata.obs_names,
    columns=[f"PC{i + 1}" for i in range(adata.obsm["X_pca"].shape[1])],
)
pca_df.to_csv(
    os.path.join(args.output_dir, f"{args.name}.pca.tsv"),
    sep="\t",
    index=True,
    index_label="cell_id",
)
