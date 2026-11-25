# Session 10: Protein language models

ESM-2: amino-acid tokens, last-layer embeddings, pseudo-log-likelihood and mutation scores.

[Student notebook](ProteinLM_practice.ipynb) · [Open in Colab](https://colab.research.google.com/github/nbrg-ppcu/appliedmedtech/blob/main/notebooks/session_10/ProteinLM_practice.ipynb)

Open the notebook in Colab, save a copy in Drive, and run the setup cells in order. The first code cell installs the packages; the model and data load directly from Hugging Face. Select a Python 3 GPU runtime; a T4 is sufficient. CPU execution is supported but slower.

| File | What it is |
|---|---|
| [S10_lab_handout.md](S10_lab_handout.md) | the situation, the data, how the session runs: **read first** |
| [protein_DICTIONARY.md](data/protein_DICTIONARY.md) | what every dataset column means, the test flag, and the arrays created in the notebook: **read before the notebook** |
| [ProteinLM_practice.ipynb](ProteinLM_practice.ipynb) | the notebook: 12 exercises, with worked examples and all helper code inside it |
| [Course-ProteinLoc](https://huggingface.co/datasets/neuralbioinfo/Course-ProteinLoc) | the protein localization dataset; downloaded by the notebook |

**Software:** Python 3 with `torch`, `transformers`, `datasets`, `numpy`, `pandas`, `matplotlib`, `seaborn`, `umap-learn` and `tqdm`. The setup cell installs the additional packages; Colab supplies PyTorch and the core scientific Python packages. No manual uploads, Drive mount, or Hugging Face token are required.

The notebook uses 1,000 proteins by default. For a shorter CPU run, set `N_PROTEINS = 100`; use at least 51 for the UMAP comparison exercise. Exercise cells are intentionally incomplete, but the supplied examples run in order. Mutation scores describe model preferences; they are not measurements of protein function.
