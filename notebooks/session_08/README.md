# Session 8: Protein language models

[Student notebook](ProteinLM_practice.ipynb) · [Open this branch in Colab](https://colab.research.google.com/github/nbrg-ppcu/appliedmedtech/blob/session-08-protein-language-models/notebooks/session_08/ProteinLM_practice.ipynb) · [Published version after merge](https://colab.research.google.com/github/nbrg-ppcu/appliedmedtech/blob/main/notebooks/session_08/ProteinLM_practice.ipynb)

Explore ESM-2 tokenization and model outputs, extract protein embeddings, visualize localization labels with UMAP, and calculate pseudo-log-likelihood and amino-acid replacement scores. The notebook contains 12 exercises.

## Run in Google Colab

1. Open the Colab link and choose **File → Save a copy in Drive**.
2. Choose **Runtime → Change runtime type → Python 3 → GPU**. A T4 is sufficient.
3. Run the setup cell first, then work through the cells in order.

The notebook downloads [ESM-2 8M](https://huggingface.co/facebook/esm2_t6_8M_UR50D) and the [protein localization dataset](https://huggingface.co/datasets/neuralbioinfo/Course-ProteinLoc) directly from Hugging Face. No manual file uploads, Drive mount, or Hugging Face token are required.

CPU execution is supported but slower. For a short CPU run, set `N_PROTEINS = 100`; keep at least 51 proteins for the UMAP comparison exercise. The default uses 1,000 proteins and batches of 16. Inputs are limited to 512 tokens, including the boundary tokens.

The setup cell installs Transformers 4.57.1, Datasets 3.x–4.x, UMAP 0.5.9.post2 and Seaborn. Colab provides PyTorch, NumPy, pandas, Matplotlib and tqdm. If packages were already imported, restart the session after installation and rerun from the top. A missing `HF_TOKEN` warning is harmless for these public resources.

Exercise cells are intentionally incomplete. All supplied examples can run in order without completing the exercises. Mutation scores describe model preferences and are not measurements of biological function.
