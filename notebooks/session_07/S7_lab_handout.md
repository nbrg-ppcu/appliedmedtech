# A bone cell atlas

---

## Where we are

Every session so far came with the answer: a label the model had to learn. Today there is none. The
lecture was about what you can still do without labels: **make the data visible** (dimension
reduction) and **find its groups** (clustering). Then comes the hard part: deciding what those groups
mean.

## The question

Four hospitals pooled single-cell RNA-seq from human bone and bone marrow: 1,000 cells from 16 donors.
Nobody has labelled the cells. The consortium asks three things:

1. **Which cell types are in bone, and how are they related?**
2. Facility D reports a group of cells that no other facility has. **Is it a new cell type?**
3. **Does osteoporosis shift the balance between the cells that dissolve bone and the cells that build
   it?**

At the end, an expert panel's labels are revealed, so you can check how close your clusters came.

## The data

Each row is one cell, with 62 genes and some information about where the cell came from. The notebook
starts with a short introduction to the cells of bone; no biology background is needed.

| File | What it is |
|---|---|
| `bone_cells.csv` | 1,000 cells from 16 donors, 4 facilities: metadata and 62 genes, **no cell types** |
| `bone_labels.csv` | the expert panel's cell types, at 3 levels: used **only in the last part** |
| `bone_DICTIONARY.md` | what every column means, who the donors of each facility are, and a marker table for naming clusters: **read it first** |

## What you will do, and why

**Look at the table.** Before any method: who are the donors, and how did each facility handle its
tissue? Some of this can make cells look different for reasons that have nothing to do with their type.

**Make it visible.** 62 genes give 1,891 possible scatter plots. You reduce them to two axes with four
methods (PCA, PCoA, NMDS and UMAP) and see what each one keeps and what it loses.

**Build the tree.** Hierarchical clustering gives a whole tree of groupings. You choose where to cut it,
name the clusters from their marker genes, and find out which cell types are relatives.

**Look for density, and use the graph.** HDBSCAN finds dense groups by itself and may call a cell
"noise"; Leiden, the default of single-cell tools, splits a graph of neighbouring cells into
communities. You compare them with the tree, score them with the silhouette - a score that needs no
labels - and find out which cells HDBSCAN leaves out, and why.

**Question every cluster.** One cluster comes almost only from facility D. You decide whether it is a
cell type or something the handling of the tissue did to the cells.

**The reveal.** You score your clusterings against the expert labels at three levels with the
adjusted Rand index, and answer the osteoporosis question, counting the samples you really have and
asking what else could explain the difference.

## Questions you will be able to answer

- What does each dimension reduction method keep, and how should you (not) read a UMAP map?
- How do you choose the number of clusters, what does the linkage change, and what can a silhouette
  score tell you without labels?
- When is "noise" the right answer for a cell?
- How can you tell a real cell type from a cluster made by the experiment?
- When is the number of cells not the number of samples?

## What to hand in

- the **notebook**, with your conclusions, your answer to Task 10 and the three-line report, uploaded
  to Moodle

---

> **Note.** The data are synthetic. They were built to behave like real bone single-cell data, but no
> real person is in them and no number in them is a scientific result.
