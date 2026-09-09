# Applications of Machine Learning in Biotech and Medtech

**P-ITSZT-0061** · Course materials for MSc students in Bioinformatics and Info-bionics.

Exercises use Google Colab. Start with the Session 1 notebook below.

## Session 1 — Can this gene be switched on?

Build a promoter classifier with k-mer counts and logistic regression, inspect its predictions, and investigate how the choice of negative examples changes its performance.

[Open the notebook in Google Colab](https://colab.research.google.com/github/nbrg-ppcu/appliedmedtech/blob/session-01-promoters/notebooks/session_01/promoters.ipynb) · [View the notebook on GitHub](notebooks/session_01/promoters.ipynb)

1. Open the notebook in Colab.
2. Save a personal copy to Google Drive.
3. Run the cells from the top, then change one setting at a time and record the results.

The notebook installs its additional dependency and downloads its dataset from Hugging Face. No local setup is required for Colab.

## Repository structure

| Path | Purpose |
| --- | --- |
| `notebooks/session_01/` | Session 1 practical notebook |
| `materials/session_01/slides/` | Session 1 slide decks |
| `materials/session_01/handouts/` | Session 1 instructions and worksheets |
| `data/` | Dataset sources and retrieval notes |
| `assets/` | Shared figures, diagrams, and reusable course assets |

Add later sessions using the same two-digit convention: `session_02`, `session_03`, and so on.

The slide and handout folders are ready for materials to be added. Keep downloaded data, local environments, caches, and generated experiment outputs out of Git; the `.gitignore` covers these locations.

## Dataset

Session 1 uses [neuralbioinfo/bacterial_promoters](https://huggingface.co/datasets/neuralbioinfo/bacterial_promoters). See [data/README.md](data/README.md) for the revision used by the notebook.

## License

See [LICENSE](LICENSE) for the repository license. Refer to each external dataset's source for its own license and usage terms.
