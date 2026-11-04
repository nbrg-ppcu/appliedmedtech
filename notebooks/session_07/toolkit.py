"""Toolkit for Session 7: a bone cell atlas (dimension reduction and clustering).

You do not need to read this file to do the lab: every function is shown in use in the notebook, and
`help(function_name)` prints its description. It is here so that nothing is a black box.

How to read this file:
  - Big blocks of '#' lines separate the sections (SETTINGS, DATA, VISUALIZATION, CLUSTERING, ...).
  - In each section the functions the notebook uses come first; below a thin '# ---------- internal ...'
    line come the internal helpers they need.
  - Names starting with _ (like _scatter) are internal helpers: the notebook never calls them.
  - Every function gets its data through its arguments - nothing is hidden in global variables.
    `cells` is always the table from load_bone_cells(); `data` is any table of numbers with one row
    per cell (the gene values cells[GENES], or coordinates).
  - Every docstring ends with the internal helpers that function calls, or, for an internal helper,
    the functions that use it.

The functions come in families. Functions that work together are listed together:

  DATA                load_bone_cells(), reveal_expert_labels()
  LOOKING AT GENES    plot_two_genes()
  DIMENSION REDUCTION compute_pca(), compute_pca_loadings(), compute_distance_matrix(), compute_pcoa(),
                      compute_nmds(), compute_umap()
  VISUALIZATION       plot_dimred() - the one plotting function for every dimension reduction -
                      and plot_explained_variance()
  CLUSTERING          build_cluster_tree() -> plot_cluster_tree() / cut_cluster_tree();
                      cluster_hdbscan(); cluster_leiden()
  DESCRIBING CLUSTERS find_marker_genes(), cluster_composition(), cluster_quality(), compare_clusterings()
  SCORING CLUSTERS    mean_silhouette(), plot_silhouette_by_cluster_count() (no labels needed);
                      score_against_experts(), count_cell_types_per_donor() (with the expert labels)

Dimension reduction and clustering are separate: clustering works on whatever `data` you give it. In
the lab it gets the 62 gene values; a map from dimension reduction is only used to look at the clusters.
"""
from __future__ import annotations

import inspect
import random
import warnings
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import HDBSCAN
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.neighbors import NearestNeighbors

warnings.filterwarnings("ignore")

__all__ = [
    "GENES", "DISTANCES", "MARKER_GENES",
    "load_bone_cells", "reveal_expert_labels", "plot_two_genes",
    "compute_pca", "compute_pca_loadings", "compute_distance_matrix", "compute_pcoa", "compute_nmds", "compute_umap",
    "plot_dimred", "plot_explained_variance",
    "build_cluster_tree", "plot_cluster_tree", "cut_cluster_tree", "cluster_hdbscan", "cluster_leiden",
    "find_marker_genes", "cluster_composition", "cluster_quality", "compare_clusterings",
    "mean_silhouette", "plot_silhouette_by_cluster_count", "score_against_experts", "count_cell_types_per_donor",
]

Labels = Sequence | pd.Series            # one label per cell (a clustering)


####################################################################################################
# SETTINGS - THE GENES, THE DISTANCES AND THE MARKER TABLE
####################################################################################################
GENES: list[str] = [
    "PTPRC", "PECAM1", "CDH5", "PDGFRB", "COL1A1", "RUNX2", "SP7", "LEPR", "CXCL12", "SOX9", "ALPL", "BGLAP",
    "IBSP", "SPP1", "SOST", "DMP1", "PHEX", "ADIPOQ", "LPL", "COL2A1", "ACAN", "COL10A1", "GJA5", "SEMA3G",
    "STAB2", "FLT4", "TYROBP", "CSF1R", "LYZ", "CD14", "FCN1", "C1QA", "CD163", "S100A8", "S100A9", "CSF3R",
    "ACP5", "CTSK", "MMP9", "DCSTAMP", "CD79A", "MS4A1", "CD3E", "CD3D", "IL7R", "CD34", "KIT", "HBB", "GYPA",
    "ACTB", "B2M", "MALAT1", "FOS", "JUN", "HSPA1A", "DUSP1", "XIST", "RPS4Y1", "MKI67", "TOP2A", "CDKN2A",
    "IL1B"]
"""The 62 gene columns of bone_cells.csv, in file order. cells[GENES] is the table of gene values."""

DISTANCES: tuple[str, ...] = ("euclidean", "correlation", "braycurtis", "aitchison")
"""The distances compute_pcoa(), compute_nmds() and compute_distance_matrix() understand."""

MARKER_GENES: dict[str, str] = {
    "osteoprogenitor": "RUNX2, ALPL, COL1A1; SP7 low",
    "osteoblast": "RUNX2, SP7, BGLAP, IBSP, SPP1, COL1A1 high",
    "osteocyte": "SOST, DMP1, PHEX, BGLAP",
    "Adipo-CAR": "LEPR, CXCL12, ADIPOQ, LPL",
    "Osteo-CAR": "LEPR, CXCL12, ALPL, SPP1",
    "chondrocyte": "SOX9, COL2A1, ACAN, COL10A1",
    "arterial EC": "PECAM1, CDH5, GJA5, SEMA3G",
    "sinusoidal EC": "PECAM1, CDH5, STAB2, FLT4",
    "classical monocyte": "PTPRC, TYROBP, CD14, LYZ, FCN1",
    "macrophage": "PTPRC, TYROBP, C1QA, CD163",
    "neutrophil": "PTPRC, S100A8, S100A9, CSF3R",
    "osteoclast": "PTPRC, TYROBP, CSF1R, ACP5, CTSK, MMP9, DCSTAMP",
    "B cell": "PTPRC, CD79A, MS4A1",
    "T cell": "PTPRC, CD3E, CD3D, IL7R",
    "HSPC": "PTPRC, CD34, KIT; MKI67",
    "erythroid progenitor": "PTPRC (low), HBB, GYPA; MKI67",
}
"""The marker table of the data dictionary: cell type -> the genes to look for. Names match bone_labels.csv."""


####################################################################################################
# DATA - LOADING THE CELLS AND THE EXPERT LABELS
# load_bone_cells() first; reveal_expert_labels() only in part F.
####################################################################################################
def load_bone_cells(data_dir: str) -> pd.DataFrame:
    """Read bone_cells.csv from `data_dir` (a folder path or a web address).

    Returns the table `cells`: one row per cell, 13 metadata columns, then the 62 genes. Most other
    functions take this table (or a part of it) as their first argument.

    Internal helpers it calls: _folder().
    """
    return pd.read_csv(_folder(data_dir) + "bone_cells.csv")


def reveal_expert_labels(data_dir: str, cells: pd.DataFrame) -> pd.DataFrame:
    """Read the expert labels (bone_labels.csv), put in the same row order as `cells`.

    Only for part F. Returns the table `expert_labels` with the columns cell_id, compartment, lineage and
    cell_type; score_against_experts() and count_cell_types_per_donor() take it as an argument.

    Internal helpers it calls: _folder().
    """
    expert_labels = pd.read_csv(_folder(data_dir) + "bone_labels.csv").set_index("cell_id")
    return expert_labels.loc[cells["cell_id"]].reset_index()


# ---------- internal helpers of this section ---------------------------------------------------------
def _folder(data_dir: str) -> str:
    """The folder path or web address with a '/' at the end.

    Internal helper - the notebook never calls it directly. Used by: load_bone_cells(),
        reveal_expert_labels().
    """
    return data_dir if data_dir.endswith("/") else data_dir + "/"


def _values_per_cell(cells: pd.DataFrame | None, color: str | Labels) -> pd.Series:
    """Turn a column name or gene name (looked up in `cells`) or a list of labels (one per cell) into one
    value per cell. Used by every plotting and describing function, so they all accept the same kinds of
    `color` / `clusters` argument.

    Internal helper - the notebook never calls it directly. Used by: cluster_composition(),
        cluster_quality(), compare_clusterings(), find_marker_genes(), mean_silhouette(),
        plot_dimred(), plot_two_genes(), score_against_experts().
    """
    if isinstance(color, str):
        if cells is None:
            raise ValueError(f"To colour by '{color}', also give the table: cells=cells.")
        if color not in cells.columns:
            raise KeyError(f"'{color}' is neither a metadata column nor a gene of the table.")
        return cells[color].reset_index(drop=True)
    return pd.Series(np.asarray(color))


####################################################################################################
# LOOKING AT GENES
####################################################################################################
def plot_two_genes(cells: pd.DataFrame, gene_x: str, gene_y: str, color: str | Labels | None = None) -> None:
    """Two genes of `cells` against each other, one dot per cell, coloured like plot_dimred().

    Internal helpers it calls: _colour_title(), _scatter(), _values_per_cell().
    """
    ax = plt.figure(figsize=(6.5, 5.5)).gca()
    values = None if color is None else _values_per_cell(cells, color)
    _scatter(cells[[gene_x, gene_y]].values, [gene_x, gene_y], values, color, _colour_title(color), ax)
    plt.show()


####################################################################################################
# DIMENSION REDUCTION - COMPUTING NEW COORDINATES
# Each compute_...() returns a table of coordinates, one row per cell, that plot_dimred() can draw.
####################################################################################################
def compute_pca(cells: pd.DataFrame, n_components: int = 20, genes: Sequence[str] = GENES) -> pd.DataFrame:
    """Principal component analysis (PCA) of the genes in `cells`.

    PCA finds new axes (principal components, PCs), each a weighted mix of the genes: PC1 is the
    direction along which the cells differ most, PC2 the most of what is left, and so on. The genes are
    used as they are, because they already share one log scale.

    Returns `pca_coordinates`: one row per cell, the columns PC1 ... PC<n_components>. The fitted model
    travels with the table (in `.attrs`), so plot_dimred(method="pca", ...), plot_explained_variance() and
    compute_pca_loadings() can show the explained variance and the gene weights.
    """
    model = PCA(n_components=n_components, random_state=0).fit(cells[list(genes)])
    pca_coordinates = pd.DataFrame(model.transform(cells[list(genes)]),
                                   columns=[f"PC{i + 1}" for i in range(n_components)])
    pca_coordinates.attrs.update(method="pca", model=model, genes=list(genes))
    return pca_coordinates


def compute_pca_loadings(pca_coordinates: pd.DataFrame, component: int = 1) -> pd.DataFrame:
    """The weight of every gene on one principal component, sorted from the largest to the most negative.

    pca_coordinates  the output of compute_pca().
    component        which component: 1 for PC1, 2 for PC2, ...
    Genes at the top are high at one end of the component, genes at the bottom at the other end; genes
    near 0 do not matter for it. Look at the two ends with .head() and .tail().
    Also used by plot_dimred(method="pca", show_loadings=True).
    """
    model = pca_coordinates.attrs["model"]
    weights = pd.Series(model.components_[component - 1], index=pca_coordinates.attrs["genes"])
    weights = weights.sort_values(ascending=False)
    return pd.DataFrame({"gene": weights.index, f"weight on PC{component}": weights.round(3).values})


def compute_distance_matrix(cells: pd.DataFrame, distance: str = "euclidean",
                            genes: Sequence[str] = GENES) -> np.ndarray:
    """The distance between every pair of cells: a 1,000 x 1,000 table. Used by compute_pcoa() and compute_nmds().

    distance:
      "euclidean"   straight-line distance over the genes (the log values as they are).
      "correlation" 1 - the correlation of two cells' gene profiles: only the pattern counts, not the level.
      "braycurtis"  the share of the total amount that differs, gene by gene. It is computed on the
                    amounts (counts per 10,000, i.e. the log values turned back), as in microbiome work.
      "aitchison"   the distance for compositions (shares of a whole). It is the Euclidean distance after
                    centring each cell's log values on that cell's mean, so it compares the ratios
                    between genes, not the totals. (Pseudocount 1, already in the log1p values.)
    """
    values = cells[list(genes)].values
    if distance == "braycurtis":
        return squareform(pdist(np.expm1(values), metric="braycurtis"))
    if distance == "aitchison":
        centred_log_ratios = values - values.mean(axis=1, keepdims=True)
        return squareform(pdist(centred_log_ratios, metric="euclidean"))
    if distance in ("euclidean", "correlation"):
        return squareform(pdist(values, metric=distance))
    raise ValueError(f"distance must be one of {DISTANCES}, not '{distance}'.")


def compute_pcoa(cells: pd.DataFrame, distance: str = "euclidean", n_axes: int = 10,
                 genes: Sequence[str] = GENES) -> pd.DataFrame:
    """Principal coordinate analysis (PCoA), also called classical multidimensional scaling.

    You give it the cells and a distance (see compute_distance_matrix() for the choices), and it finds
    coordinates whose distances match the distances between the cells as well as possible. With the
    Euclidean distance it gives the same map as PCA, up to the sign of the axes.

    Returns a table with the columns PCo1 ... PCo<n_axes>. The distance and the share of the distance
    information the first two axes keep travel with the table, for plot_dimred(method="pcoa", ...).
    """
    D = compute_distance_matrix(cells=cells, distance=distance, genes=genes)
    n = len(D)
    centring = np.eye(n) - np.ones((n, n)) / n
    inner_products = -0.5 * centring @ (D ** 2) @ centring
    eigenvalues, eigenvectors = np.linalg.eigh(inner_products)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    coordinates = eigenvectors[:, :n_axes] * np.sqrt(np.maximum(eigenvalues[:n_axes], 0))
    kept = eigenvalues[:2].sum() / eigenvalues[eigenvalues > 0].sum()
    pcoa_coordinates = pd.DataFrame(coordinates, columns=[f"PCo{i + 1}" for i in range(n_axes)])
    pcoa_coordinates.attrs.update(method="pcoa", distance=distance, kept=float(kept))
    return pcoa_coordinates


def compute_nmds(cells: pd.DataFrame, distance: str = "euclidean", seed: int = 0,
                 genes: Sequence[str] = GENES) -> pd.DataFrame:
    """Non-metric multidimensional scaling (NMDS) in 2 dimensions.

    Like PCoA it starts from the distances between the cells (see compute_distance_matrix() for
    `distance`), but it keeps only their ORDER: which cell is closer to which. It improves a starting map
    step by step; `seed` sets the random parts. Takes about half a minute.

    Returns a table with the columns NMDS1 and NMDS2. The stress (how badly the map breaks the order of
    the distances: 0 perfect, below 0.2 usable - scikit-learn's normalised stress) travels with the table,
    for plot_dimred(method="nmds", ...). Works with old and new scikit-learn versions.
    """
    D = compute_distance_matrix(cells=cells, distance=distance, genes=genes)
    parameters = inspect.signature(MDS.__init__).parameters
    settings = dict(n_components=2, n_init=1, max_iter=300, random_state=seed)
    if "metric_mds" in parameters:                 # scikit-learn 1.8 and newer
        model = MDS(metric_mds=False, metric="precomputed", init="classical_mds", **settings)
        coordinates = model.fit_transform(D)
    else:                                          # older scikit-learn
        model = MDS(metric=False, dissimilarity="precomputed", **settings)
        start = compute_pcoa(cells=cells, distance=distance, n_axes=2, genes=genes).values
        coordinates = model.fit_transform(D, init=start)
    nmds_coordinates = pd.DataFrame(coordinates, columns=["NMDS1", "NMDS2"])
    nmds_coordinates.attrs.update(method="nmds", distance=distance, stress=float(model.stress_))
    return nmds_coordinates


def compute_umap(input_coordinates: pd.DataFrame, n_components: int = 2, n_neighbors: int = 15,
                 min_dist: float = 0.3, seed: int = 0) -> pd.DataFrame:
    """A UMAP map of the cells (uniform manifold approximation and projection).

    input_coordinates  what UMAP works on: usually pca_coordinates, as is standard for single-cell data
                       (the gene values cells[GENES] work too).
    n_components       how many UMAP axes: 2 for a map you can look at.
    n_neighbors        how many nearest cells count as neighbours: small shows fine detail, large the
                       overall shape.
    min_dist           how tightly points may pack.
    seed               the random start.
    UMAP links each cell to its nearest neighbours and lays them out so that those neighbourhoods stay
    together. Groups of dots on the map are called islands. Takes 10 to 20 seconds.

    Returns a table with the columns UMAP1, UMAP2, ..., with its settings attached for plot_dimred().
    """
    import umap

    model = umap.UMAP(n_components=n_components, n_neighbors=n_neighbors, min_dist=min_dist, random_state=seed)
    umap_coordinates = pd.DataFrame(model.fit_transform(np.asarray(input_coordinates)),
                                    columns=[f"UMAP{i + 1}" for i in range(n_components)])
    umap_coordinates.attrs.update(method="umap", n_neighbors=n_neighbors, min_dist=min_dist, seed=seed)
    return umap_coordinates


# ---------- internal: method name -> compute function, used by plot_dimred() -------------------------
_COMPUTE = {"pca": compute_pca, "pcoa": compute_pcoa, "nmds": compute_nmds, "umap": compute_umap}


####################################################################################################
# VISUALIZATION - PLOTTING A DIMENSION REDUCTION
# plot_dimred() is the one function the notebook calls. It picks the internal plotter for the method
# (_plot_pca, _plot_pcoa, _plot_nmds, _plot_umap); every plotter draws with _scatter().
####################################################################################################
def plot_dimred(method: str, coordinates: pd.DataFrame | None = None, cells: pd.DataFrame | None = None,
                color: str | Labels | list[str] | None = None, x_axis: str | None = None,
                y_axis: str | None = None, show_loadings: bool = False, n_columns: int = 3,
                **method_settings) -> None:
    """Plot a dimension reduction: the one plotting function for PCA, PCoA, NMDS and UMAP.

    method           "pca", "pcoa", "nmds" or "umap": the method that made (or should make) the coordinates.
    coordinates      the coordinates you made earlier, e.g. pca_coordinates = compute_pca(cells=cells).
                     If you leave it out, new coordinates are computed now with `method_settings`.
    cells            the table from load_bone_cells(): needed to colour by a column or gene name.
    color            what colours the dots: a metadata column ("facility") or gene ("PTPRC") of `cells`,
                     a clustering (one label per cell) - or a LIST of names, which draws one panel each.
    x_axis, y_axis   which two axes to draw (PCA and PCoA only), e.g. x_axis="PC3", y_axis="PC4".
                     Default: the first two.
    show_loadings    PCA only: also draw the genes with the largest weights on the two components shown
                     (computed with compute_pca_loadings(), drawn by _plot_loadings()).
    n_columns        panels per row.
    method_settings  settings for computing new coordinates when `coordinates` is left out, e.g.
                     plot_dimred(method="nmds", cells=cells, color="PTPRC", distance="braycurtis") or
                     plot_dimred(method="umap", cells=cells, color="PTPRC", input_coordinates=pca_coordinates,
                                 n_neighbors=50).

    How it works: plot_dimred() looks up `method` and calls the matching private plotter - _plot_pca(),
    _plot_pcoa(), _plot_nmds() or _plot_umap() - which adds what is special about that method (the % of
    variation, the distance and stress, the UMAP settings) and draws with _scatter().
    It only draws. To keep coordinates for later, compute them first with a compute_...() function.

    Internal helpers it calls: _colour_title(), _plot_loadings(), _plot_nmds(), _plot_pca(),
        _plot_pcoa(), _plot_umap(), _values_per_cell().
    """
    method = method.lower()
    if method not in _PLOT:
        raise ValueError(f"method must be one of {list(_PLOT)}, not '{method}'.")
    if coordinates is None:
        if method != "umap":
            method_settings.setdefault("cells", cells)
        coordinates = _COMPUTE[method](**method_settings)
    elif method_settings:
        print(f"Note: coordinates were given, so these settings were not used: {method_settings}")
    many_colourings = isinstance(color, list) and len(color) > 0 and all(isinstance(c, str) for c in color) \
        and len(color) != len(coordinates)
    colourings = color if many_colourings else [color]
    panels = [("scatter", c) for c in colourings]
    if show_loadings and method == "pca":
        panels += [("loadings", int(axis[2:])) for axis in (x_axis or "PC1", y_axis or "PC2")]
    n_columns = min(n_columns, len(panels))
    n_rows = int(np.ceil(len(panels) / n_columns))
    fig, axes = plt.subplots(n_rows, n_columns, figsize=(5.4 * n_columns, 4.5 * n_rows), squeeze=False)
    for ax, (kind, what) in zip(axes.flat, panels):
        if kind == "scatter":
            values = None if what is None else _values_per_cell(cells, what)
            _PLOT[method](coordinates, values, _colour_title(what), x_axis, y_axis, ax)
        else:
            _plot_loadings(coordinates, what, ax)
    for ax in axes.flat[len(panels):]:
        ax.axis("off")
    plt.tight_layout()
    plt.show()


def plot_explained_variance(pca_coordinates: pd.DataFrame) -> None:
    """Bar chart of the share of the variation each principal component carries (a 'scree plot').

    pca_coordinates is the output of compute_pca(). Part of the PCA family:
    compute_pca() -> plot_explained_variance() / compute_pca_loadings() / plot_dimred(method="pca").
    """
    shares = pca_coordinates.attrs["model"].explained_variance_ratio_ * 100
    plt.figure(figsize=(7, 3.5))
    plt.bar(range(1, len(shares) + 1), shares, color="0.45")
    for i, share in enumerate(shares[:5]):
        plt.text(i + 1, share + 0.3, f"{share:.0f}%", ha="center", fontsize=8)
    plt.xlabel("principal component")
    plt.ylabel("% of the variation")
    plt.xticks(range(1, len(shares) + 1))
    plt.show()


# ---------- internal plotters: one per method, called by plot_dimred() -------------------------------
def _plot_pca(coordinates: pd.DataFrame, values: pd.Series | None, label: str, x_axis: str | None,
              y_axis: str | None, ax) -> None:
    """Private plotter of plot_dimred(method="pca"): any two PCs, with their % of the variation on the axes.

    Internal helper - the notebook never calls it directly. Used by: plot_dimred().
    """
    x_axis, y_axis = x_axis or "PC1", y_axis or "PC2"
    shares = coordinates.attrs["model"].explained_variance_ratio_
    names = [f"{axis} ({shares[int(axis[2:]) - 1]:.0%} of the variation)" for axis in (x_axis, y_axis)]
    _scatter(coordinates[[x_axis, y_axis]].values, names, values, label, f"PCA · {label}", ax)


def _plot_loadings(pca_coordinates: pd.DataFrame, component: int, ax, n_genes_each_end: int = 6) -> None:
    """Private plotter of plot_dimred(method="pca", show_loadings=True): the genes at both ends of one PC.

    Uses compute_pca_loadings(); bars to the right are genes high at the positive end of the component.

    Internal helper - the notebook never calls it directly. Used by: plot_dimred().
    """
    loadings = compute_pca_loadings(pca_coordinates=pca_coordinates, component=component)
    ends = pd.concat([loadings.tail(n_genes_each_end), loadings.head(n_genes_each_end)])
    weights = ends.iloc[:, 1].values
    ax.barh(range(len(ends)), weights, color=["#7b6fb0" if w < 0 else "#2f6355" for w in weights])
    ax.set_yticks(range(len(ends)))
    ax.set_yticklabels(ends["gene"], fontsize=8)
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_xlabel(f"weight on PC{component}")
    ax.set_title(f"the genes of PC{component}")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _plot_pcoa(coordinates: pd.DataFrame, values: pd.Series | None, label: str, x_axis: str | None,
               y_axis: str | None, ax) -> None:
    """Private plotter of plot_dimred(method="pcoa"): two principal coordinates, the distance in the title.

    Internal helper - the notebook never calls it directly. Used by: plot_dimred().
    """
    x_axis, y_axis = x_axis or "PCo1", y_axis or "PCo2"
    distance = coordinates.attrs.get("distance", "?")
    kept = coordinates.attrs.get("kept")
    first_two = kept is not None and (x_axis, y_axis) == ("PCo1", "PCo2")
    extra = f", first 2 axes keep {kept:.0%}" if first_two else ""
    _scatter(coordinates[[x_axis, y_axis]].values, [x_axis, y_axis], values, label,
             f"PCoA ({distance}{extra}) · {label}", ax)


def _plot_nmds(coordinates: pd.DataFrame, values: pd.Series | None, label: str, x_axis: str | None,
               y_axis: str | None, ax) -> None:
    """Private plotter of plot_dimred(method="nmds"): the 2D map, with the distance and the stress in the title.

    Internal helper - the notebook never calls it directly. Used by: plot_dimred().
    """
    distance = coordinates.attrs.get("distance", "?")
    stress = coordinates.attrs.get("stress")
    stress_text = f", stress {stress:.2f}" if stress is not None else ""
    _scatter(coordinates[["NMDS1", "NMDS2"]].values, ["NMDS1", "NMDS2"], values, label,
             f"NMDS ({distance}{stress_text}) · {label}", ax)


def _plot_umap(coordinates: pd.DataFrame, values: pd.Series | None, label: str, x_axis: str | None,
               y_axis: str | None, ax) -> None:
    """Private plotter of plot_dimred(method="umap"): the first two UMAP axes, its settings in the title.

    The axes carry no ticks: positions and distances on a UMAP map have no meaning of their own.

    Internal helper - the notebook never calls it directly. Used by: plot_dimred().
    """
    settings = coordinates.attrs
    settings_text = f"n_neighbors={settings.get('n_neighbors', '?')}, seed={settings.get('seed', '?')}"
    _scatter(coordinates[["UMAP1", "UMAP2"]].values, ["UMAP1", "UMAP2"], values, label,
             f"UMAP ({settings_text}) · {label}", ax, dot_size=6)
    ax.set_xticks([])
    ax.set_yticks([])


# ---------- internal: method name -> plotter, used by plot_dimred() ----------------------------------
_PLOT = {"pca": _plot_pca, "pcoa": _plot_pcoa, "nmds": _plot_nmds, "umap": _plot_umap}


# ---------- internal drawing helpers: used by every plot ---------------------------------------------
def _scatter(xy: np.ndarray, axis_names: Sequence[str], values: pd.Series | None, color, title: str,
             ax: plt.Axes, dot_size: float = 8) -> None:
    """Draw one scatter plot of the cells.

    Colouring rules: a gene or a numeric metadata column (given by name) gets a colour scale with a colour
    bar; anything else (facility, donor, clusters) gets one colour per group and a legend. Groups called
    -1 or "noise" are drawn in light grey, behind the others.

    Internal helper - the notebook never calls it directly. Used by: _plot_nmds(), _plot_pca(),
        _plot_pcoa(), _plot_umap(), plot_two_genes().
    """
    if values is None:
        ax.scatter(xy[:, 0], xy[:, 1], s=dot_size, color="0.5")
    elif isinstance(color, str) and pd.api.types.is_numeric_dtype(values):
        points = ax.scatter(xy[:, 0], xy[:, 1], s=dot_size, c=values, cmap="viridis")
        plt.colorbar(points, ax=ax, shrink=0.8)
    else:
        groups = sorted(values.unique(), key=_legend_order)
        for group, colour in zip(groups, _palette(len(groups))):
            in_group = (values == group).values
            is_noise = str(group) in ("-1", "noise")
            ax.scatter(xy[in_group, 0], xy[in_group, 1], s=dot_size, label="noise" if is_noise else group,
                       color="0.82" if is_noise else colour, zorder=1 if is_noise else 2)
        ax.legend(fontsize=7, markerscale=2, frameon=False, bbox_to_anchor=(1.0, 1.0), loc="upper left",
                  ncol=2 if len(groups) > 16 else 1)
    ax.set_xlabel(axis_names[0])
    ax.set_ylabel(axis_names[1])
    ax.set_title(title)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _palette(n_groups: int) -> list:
    """n distinct colours without greys (grey is kept for noise).

    Internal helper - the notebook never calls it directly. Used by: _scatter().
    """
    if n_groups <= 9:
        return [plt.get_cmap("tab10")(i) for i in range(10) if i != 7][:n_groups]
    base = [plt.get_cmap("tab20")(i) for i in range(20) if i not in (14, 15)]
    base += [plt.get_cmap("tab20b")(i) for i in range(20)]
    return [base[i % len(base)] for i in range(n_groups)]


def _legend_order(group) -> tuple:
    """Sort order for legend entries: noise last, numbers in numeric order, text alphabetically.

    Internal helper - the notebook never calls it directly. Used by: _scatter(),
        find_marker_genes().
    """
    is_noise = str(group) in ("-1", "noise")
    try:
        return (is_noise, 0, float(group), "")
    except (TypeError, ValueError):
        return (is_noise, 1, 0.0, str(group))


def _colour_title(color) -> str:
    """The title text for a colouring: the column or gene name, 'clusters', or nothing.

    Internal helper - the notebook never calls it directly. Used by: plot_dimred(),
        plot_two_genes().
    """
    if color is None:
        return ""
    return color if isinstance(color, str) else "clusters"


####################################################################################################
# CLUSTERING - HIERARCHICAL TREE
# build_cluster_tree() -> plot_cluster_tree() to look at it -> cut_cluster_tree() to get clusters.
# `data` is any table of numbers with one row per cell; in the lab: the gene values cells[GENES].
####################################################################################################
def build_cluster_tree(data: pd.DataFrame | np.ndarray, linkage_method: str = "ward") -> np.ndarray:
    """Hierarchical clustering: the full merge history of the cells (a 'linkage matrix').

    data            one row per cell, e.g. gene_values = cells[GENES].
    linkage_method  how the distance between two GROUPS of cells is measured:
      "ward"     merge the pair whose merge makes the spread inside the groups grow least (compact groups)
      "average"  the mean distance over all pairs of cells
      "complete" the farthest pair of cells
      "single"   the closest pair of cells (tends to make chains)
    Starting from single cells, it merges the two closest groups again and again until one group is left.
    Returns `cluster_tree`, for plot_cluster_tree() and cut_cluster_tree().
    """
    return linkage(np.asarray(data), method=linkage_method)


def plot_cluster_tree(cluster_tree: np.ndarray, n_clusters: int | None = None, title: str = "") -> None:
    """Draw the top 40 branches of a tree from build_cluster_tree(); with n_clusters, a dashed line where
    it falls into that many clusters.

    The number under a branch is how many cells it holds. A long vertical line means two groups that are
    far apart.
    """
    n_clusters = n_clusters if n_clusters and n_clusters > 1 else None    # a cut needs at least 2 clusters
    threshold = cluster_tree[-(n_clusters - 1), 2] if n_clusters else 0
    plt.figure(figsize=(11, 4.5))
    dendrogram(cluster_tree, truncate_mode="lastp", p=40, color_threshold=threshold,
               above_threshold_color="0.6", show_leaf_counts=True, leaf_font_size=8)
    if n_clusters:
        cut_height = (cluster_tree[-(n_clusters - 1), 2] + cluster_tree[-n_clusters, 2]) / 2
        plt.axhline(cut_height, ls="--", color="k", lw=1)
    plt.ylabel("merge height (distance)")
    plt.xlabel("branch (number of cells in brackets)")
    plt.title(title or (f"cut into {n_clusters} clusters" if n_clusters else ""))
    plt.show()


def cut_cluster_tree(cluster_tree: np.ndarray, n_clusters: int) -> pd.Series:
    """Cut a tree from build_cluster_tree() into n_clusters clusters: one cluster number (1 ... n) per cell,
    in the same order as the rows the tree was built from."""
    return pd.Series(fcluster(cluster_tree, n_clusters, criterion="maxclust"), name=f"tree cut into {n_clusters}")


####################################################################################################
# CLUSTERING - DENSITY (HDBSCAN) AND GRAPH (LEIDEN)
# Both take `data` (one row per cell, in the lab the gene values) and return one cluster label per cell.
####################################################################################################
def cluster_hdbscan(data: pd.DataFrame | np.ndarray, min_cluster_size: int = 15, min_samples: int = 5,
                    selection_mode: str = "eom") -> pd.Series:
    """Density-based clustering with HDBSCAN: clusters are crowded regions, the rest is noise (-1).

    data              one row per cell, e.g. gene_values = cells[GENES].
    min_cluster_size  the smallest group it will call a cluster.
    min_samples       how cautious it is about calling a region dense.
    selection_mode    "eom" (excess of mass, the default) keeps the groups that stay stable over the widest
                      range of densities; "leaf" keeps the smallest, finest groups (and more noise).
    You do not choose the number of clusters. Returns one label per cell; -1 means noise.
    """
    labels = HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples,
                     cluster_selection_method=selection_mode).fit_predict(np.asarray(data))
    print(f"{labels.max() + 1} clusters, {np.mean(labels == -1):.0%} of cells are noise")
    return pd.Series(labels, name="hdbscan")


def cluster_leiden(data: pd.DataFrame | np.ndarray, resolution: float = 1.0, n_neighbors: int = 15,
                   seed: int = 0) -> pd.Series:
    """Graph-based clustering with Leiden, the default of scanpy (Seurat uses the closely related Louvain).

    data          one row per cell, e.g. gene_values = cells[GENES].
    resolution    how finely to split: higher values give more, smaller clusters.
    n_neighbors   every cell is linked to its n_neighbors nearest other cells; the links make a graph.
    seed          the random start.
    Leiden splits the graph into communities: groups with many links inside and few between. Every cell
    gets a cluster (no noise). Returns one cluster number per cell.
    """
    import igraph

    values = np.asarray(data)
    neighbours = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(values).kneighbors(values)[1]
    edges = sorted({(min(i, j), max(i, j)) for i, row in enumerate(neighbours) for j in row[1:]})
    graph = igraph.Graph(n=len(values), edges=edges)
    igraph.set_random_number_generator(random.Random(seed))
    communities = graph.community_leiden(objective_function="modularity", resolution=resolution,
                                         n_iterations=-1)
    print(f"{max(communities.membership) + 1} clusters")
    return pd.Series(communities.membership, name="leiden")


####################################################################################################
# DESCRIBING CLUSTERS
# Each takes the table `cells` and a clustering (one label per cell) and returns a table.
####################################################################################################
def find_marker_genes(cells: pd.DataFrame, clusters: Labels, n_top_genes: int = 5) -> pd.DataFrame:
    """For each cluster: its size and the genes that are much higher in it than in all other cells.

    The number after each gene is the mean in the cluster minus the mean elsewhere (log scale). Compare
    them with MARKER_GENES (the marker table of the data dictionary) to name the cluster.

    Internal helpers it calls: _legend_order(), _values_per_cell().
    """
    gene_values = cells[GENES].reset_index(drop=True)
    labels = _values_per_cell(cells, clusters)
    rows = []
    for cluster in sorted(labels.unique(), key=_legend_order):
        in_cluster = (labels == cluster).values
        difference = gene_values[in_cluster].mean() - gene_values[~in_cluster].mean()
        top = difference.nlargest(n_top_genes)
        rows.append({"cluster": cluster, "cells": int(in_cluster.sum()),
                     "top genes (difference)": ", ".join(f"{gene} ({d:.1f})" for gene, d in top.items())})
    return pd.DataFrame(rows).set_index("cluster")


def cluster_composition(cells: pd.DataFrame, clusters: Labels, by_column: str = "facility") -> pd.DataFrame:
    """For each cluster: the % of its cells in each group of the column `by_column` of `cells` (facility,
    donor_id, ...), and its size.

    A real cell type should be found in every facility; a cluster from one place only needs a closer look.

    Internal helpers it calls: _values_per_cell().
    """
    labels = _values_per_cell(cells, clusters).rename("cluster")
    groups = cells[by_column].reset_index(drop=True)
    table = pd.crosstab(labels, groups, normalize="index").mul(100).round(0)
    table["cells"] = labels.value_counts().reindex(table.index)
    return table


def cluster_quality(cells: pd.DataFrame, clusters: Labels) -> pd.DataFrame:
    """For each cluster: the median n_counts (RNA molecules per cell), n_genes and pct_mito, and its size.

    Doublets usually have more n_counts than single cells; damaged or stressed cells have a high pct_mito.

    Internal helpers it calls: _values_per_cell().
    """
    labels = _values_per_cell(cells, clusters).rename("cluster")
    quality = cells[["n_counts", "n_genes", "pct_mito"]].reset_index(drop=True)
    table = quality.groupby(labels).median().round(1)
    table["cells"] = labels.value_counts().reindex(table.index)
    return table


def compare_clusterings(row_clusters: Labels, column_clusters: Labels) -> pd.DataFrame:
    """Cross-table of two clusterings of the same cells: one row per cluster of `row_clusters`, one column
    per cluster of `column_clusters`; each number counts the cells the two share.

    Internal helpers it calls: _values_per_cell().
    """
    return pd.crosstab(_values_per_cell(None, row_clusters).rename("rows"),
                       _values_per_cell(None, column_clusters).rename("columns"))


####################################################################################################
# SCORING CLUSTERS
# mean_silhouette() and plot_silhouette_by_cluster_count() need no labels;
# score_against_experts() and count_cell_types_per_donor() need the table from reveal_expert_labels().
####################################################################################################
def mean_silhouette(data: pd.DataFrame | np.ndarray, clusters: Labels) -> float:
    """The mean silhouette of a clustering, measured on `data` (noise, label -1, is left out).

    data      one row per cell: the same table the clustering was made from, e.g. gene_values.
    For one cell: a = mean distance to the other cells of its own cluster, b = mean distance to the cells
    of the nearest other cluster, silhouette = (b - a) / max(a, b). Near 1: deep inside its cluster;
    near 0: on a border; negative: it would fit better next door. Needs no labels.

    Internal helpers it calls: _values_per_cell().
    """
    labels = _values_per_cell(None, clusters).values
    kept = labels != -1
    return round(float(silhouette_score(np.asarray(data)[kept], labels[kept])), 3)


def plot_silhouette_by_cluster_count(data: pd.DataFrame | np.ndarray, cluster_tree: np.ndarray,
                                     cluster_counts: Sequence[int] = range(2, 21)) -> None:
    """The mean silhouette of `cluster_tree` (built from `data`) cut into n clusters, for every n in
    `cluster_counts`; plotted. Uses cut_cluster_tree() and mean_silhouette() for each n."""
    scores = [mean_silhouette(data=data, clusters=cut_cluster_tree(cluster_tree=cluster_tree, n_clusters=n))
              for n in cluster_counts]
    plt.figure(figsize=(7, 3.5))
    plt.plot(list(cluster_counts), scores, marker="o", color="0.3")
    plt.xlabel("number of clusters")
    plt.ylabel("mean silhouette")
    plt.xticks(list(cluster_counts))
    plt.show()


def score_against_experts(clusters: Labels, expert_labels: pd.DataFrame, level: str = "cell_type") -> pd.Series:
    """Adjusted Rand index (ARI) of a clustering against the expert labels.

    expert_labels  the table from reveal_expert_labels().
    level          "compartment", "lineage" or "cell_type".
    1 = the same grouping, 0 = no better than chance; the names of the clusters do not matter, only which
    cells are together. Returns three numbers: the ARI over all cells (noise counted as one more cluster),
    the ARI over the cells not called noise, and how many cells that second score used.

    Internal helpers it calls: _values_per_cell().
    """
    expert = expert_labels[level].values
    labels = _values_per_cell(None, clusters).values
    kept = labels != -1
    return pd.Series({"ARI, all cells": round(adjusted_rand_score(expert, labels), 3),
                      "ARI, without noise": round(adjusted_rand_score(expert[kept], labels[kept]), 3),
                      "cells kept": int(kept.sum())}, dtype=object)


def count_cell_types_per_donor(cells: pd.DataFrame, expert_labels: pd.DataFrame, cell_types: Sequence[str],
                               level: str = "cell_type") -> pd.DataFrame:
    """For each donor in `cells`: facility, condition and the number of cells of each listed expert cell type.

    expert_labels  the table from reveal_expert_labels(), in the same row order as `cells`.
    """
    counts = pd.crosstab(cells["donor_id"].values, expert_labels[level].values)[list(cell_types)]
    counts.index.name = "donor_id"
    return cells.groupby("donor_id")[["facility", "condition"]].first().join(counts)
