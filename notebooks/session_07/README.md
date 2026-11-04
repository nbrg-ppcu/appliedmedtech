# Session 7: A bone cell atlas

Unsupervised learning: dimension reduction and clustering of single bone and marrow cells.

[Student notebook](S7_bone_atlas.ipynb) · [Open in Colab](https://colab.research.google.com/github/nbrg-ppcu/appliedmedtech/blob/main/notebooks/session_07/S7_bone_atlas.ipynb)

Open the notebook in Colab, save a copy in Drive, and run the setup cells in order. The first cell installs the packages; the data load directly from GitHub. A Python CPU runtime is sufficient.

| File | What it is |
|---|---|
| [S7_lab_handout.md](S7_lab_handout.md) | the situation, the files, how the session runs: **read first** |
| [bone_DICTIONARY.md](data/bone_DICTIONARY.md) | what every column is, who the donors of each facility are, and the marker table: **read before the notebook** |
| [S7_bone_atlas.ipynb](S7_bone_atlas.ipynb) | the notebook: 11 tasks, each new tool first shown in a worked example |
| [toolkit.py](toolkit.py) | the toolkit the notebook uses (the notebook downloads it in Colab); `help(function_name)` describes each function |
| [bone_cells.csv](data/bone_cells.csv) | 1,000 cells from 16 donors: metadata and 62 genes, no cell types |
| [bone_labels.csv](data/bone_labels.csv) | the expert labels, for the last part only |

**Software:** Python 3 with `numpy`, `pandas`, `scipy`, `scikit-learn`, `umap-learn`, `igraph` and `matplotlib`. NMDS takes about half a minute and each UMAP map 5 to 15 seconds; everything else takes seconds.

The data are synthetic: built to behave like real bone single-cell data, with no real person in them.
