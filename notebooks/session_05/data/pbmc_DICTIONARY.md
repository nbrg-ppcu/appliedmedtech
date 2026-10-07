# PBMC single-cell data - data dictionary

**Files:** `pbmc_train.csv` · `pbmc_external.csv`
**Technology:** 10x Genomics 3′ gene expression, chemistry v3, with CITE-seq surface-protein
measurement

---

## What one row is

**One row is one cell.** Each cell is a peripheral blood mononuclear cell (PBMC) from one donor:
a T cell, an NK cell, a B cell or a monocyte. Each donor contributed a few hundred cells.

| File | Cells | Donors | Labs |
|---|---|---|---|
| `pbmc_train.csv` | 4,007 | 16 (D01–D16) | A and B |
| `pbmc_external.csv` | 1,533 | 6 (D17–D22) | C |

No donor appears in both files.

## The three labs

| Lab | Who | How the blood was handled | What was loaded on the chip |
|---|---|---|---|
| **A** | university research lab | processed fresh, within a few hours of the blood draw | total PBMC |
| **B** | collaborating research lab | cryopreserved, then shipped overnight before processing | enriched: CD45RO+ T cells and CD16+ monocytes |
| **C** | hospital biobank | cryopreserved biobank sample, thawed for processing | total PBMC |

All three labs used the same chemistry (10x 3′ v3) and the same sequencing depth target.

## The label

`cell_type` - the cell's type, decided by **surface-protein gating** in the CITE-seq experiment
(antibodies against CD3, CD4, CD8, CD45RA, CD45RO, CD56, CD16, CD14 and CD19), not by its RNA.

| Value | Protein definition |
|---|---|
| CD4 naive T | CD3+ CD4+ CD45RA+ CD45RO− |
| CD4 memory T | CD3+ CD4+ CD45RO+ |
| CD8 T | CD3+ CD8+ |
| NK | CD3− CD56+ |
| B | CD19+ |
| CD14 monocyte | CD14++ CD16− |
| CD16 monocyte | CD16+ (non-classical and intermediate) |

---

## Columns

### Metadata

| Column | Type | Description |
|---|---|---|
| `cell_id` | text | 10x cell barcode (16 bases and `-1`) joined to the donor ID, e.g. `CCTGATTTAGTGGCCG-1_D01`; unique |
| `donor_id` | text | the person the cell came from, D01–D22 |
| `sex` | F, M | recorded sex of the donor |
| `age` | integer | age of the donor in years |
| `condition` | healthy, SLE | **SLE** - the donor has systemic lupus erythematosus, an autoimmune disease; **healthy** - no known autoimmune disease |
| `lab` | A, B, C | the lab that processed the sample - see above |
| `chemistry` | text | 10x chemistry; the same for every cell |
| `processing` | text | `fresh` · `cryopreserved, shipped overnight` · `cryopreserved biobank sample` |
| `sorting` | text | `total PBMC` · `enriched: CD45RO+ T cells and CD16+ monocytes` |
| `n_counts` | integer | $s_c$ - total UMIs of the cell, over the whole transcriptome (definition below) |
| `n_genes` | integer | $n_c$ - number of genes detected in the cell, over the whole transcriptome |
| `pct_mito` | number | $m_c$ - percentage of the cell's UMIs from mitochondrial genes |

**UMI** stands for *unique molecular identifier*. It is a short random barcode attached to each mRNA
molecule when it is captured, before amplification. Amplification (PCR) makes many copies of each
molecule; reads with the same UMI are counted once. So a UMI count is the number of distinct mRNA
molecules captured, not the number of sequencing reads.

For cell $c$, with $x_{gc}$ the UMI count of gene $g$, $G$ the genes of the whole transcriptome and
$\mathcal{M}$ the mitochondrial genes:

$$
s_c = \sum_{g=1}^{G} x_{gc}, \qquad
n_c = \left|\{\, g : x_{gc} > 0 \,\}\right|, \qquad
m_c = 100 \cdot \frac{\sum_{g \in \mathcal{M}} x_{gc}}{s_c}
$$

### Gene expression - 40 columns

Each value is the **shifted logarithm** of the count, scaled to a fixed total per cell:

$$
y_{gc} = \ln\!\left(1 + \frac{x_{gc}}{s_c}\, L\right), \qquad L = 10^{4}
$$

- $x_{gc}$ - UMI count of gene $g$ in cell $c$
- $s_c$ - total UMIs of cell $c$ over the whole transcriptome (`n_counts`)
- $L$ - a fixed target sum: $x_{gc}/s_c \cdot L$ is counts per 10,000 (CP10k)
- $\ln$ - the natural logarithm

This is `scanpy.pp.normalize_total(target_sum=1e4)` followed by `scanpy.pp.log1p`, and Seurat's
`LogNormalize` with `scale.factor = 10000`. The normalisation was computed on the whole
transcriptome; the file keeps 40 of the genes. Values range from 0 to about 7.

Because $s_c$ comes from cell $c$ alone and $L$ is a constant, $y_{gc}$ depends on no other cell.

$y_{gc} = 0$ exactly when $x_{gc} = 0$: no transcript of the gene was captured in that cell. About
half of all values here are 0. A zero does not always mean the gene is off - single-cell
sequencing captures only a fraction of each cell's RNA.

| Gene | Full name | Commonly used as a marker of |
|---|---|---|
| CD3D | CD3 delta subunit of T-cell receptor complex | T cells |
| CD3E | CD3 epsilon subunit of T-cell receptor complex | T cells |
| IL7R | interleukin 7 receptor | T cells |
| CCR7 | C-C motif chemokine receptor 7 | naive and central-memory T cells |
| SELL | selectin L (CD62L) | naive T cells, B cells, monocytes |
| LEF1 | lymphoid enhancer binding factor 1 | naive T cells |
| S100A4 | S100 calcium binding protein A4 | memory T cells, NK cells, monocytes |
| CD8A | CD8 subunit alpha | CD8 T cells |
| CD8B | CD8 subunit beta | CD8 T cells |
| GZMK | granzyme K | memory CD8 T cells |
| CCL5 | C-C motif chemokine ligand 5 | memory and effector T cells, NK cells |
| NKG7 | natural killer cell granule protein 7 | NK cells, cytotoxic T cells |
| GNLY | granulysin | NK cells, cytotoxic T cells |
| GZMB | granzyme B | NK cells, cytotoxic T cells |
| KLRD1 | killer cell lectin like receptor D1 (CD94) | NK cells |
| PRF1 | perforin 1 | NK cells, cytotoxic T cells |
| FCGR3A | Fc gamma receptor IIIa (CD16) | NK cells, CD16 monocytes |
| MS4A1 | membrane spanning 4-domains A1 (CD20) | B cells |
| CD79A | CD79a molecule | B cells |
| HLA-DRA | major histocompatibility complex, class II, DR alpha | B cells, monocytes |
| CD14 | CD14 molecule | CD14 monocytes |
| LYZ | lysozyme | monocytes |
| S100A8 | S100 calcium binding protein A8 | monocytes |
| S100A9 | S100 calcium binding protein A9 | monocytes |
| FCN1 | ficolin 1 | monocytes |
| MS4A7 | membrane spanning 4-domains A7 | monocytes |
| CST3 | cystatin C | monocytes |
| ACTB | actin beta | — |
| B2M | beta-2-microglobulin | — |
| MALAT1 | metastasis associated lung adenocarcinoma transcript 1 (non-coding RNA) | — |
| FOS | Fos proto-oncogene, AP-1 transcription factor subunit | — |
| JUN | Jun proto-oncogene, AP-1 transcription factor subunit | — |
| DUSP1 | dual specificity phosphatase 1 | — |
| HSPA1A | heat shock protein family A (Hsp70) member 1A | — |
| XIST | X inactive specific transcript (non-coding RNA) | — |
| RPS4Y1 | ribosomal protein S4 Y-linked 1 | — |
| ISG15 | ISG15 ubiquitin like modifier | — |
| IFI6 | interferon alpha inducible protein 6 | — |
| GSTM1 | glutathione S-transferase mu 1 | — |
| ERAP2 | endoplasmic reticulum aminopeptidase 2 | — |

### Label

| Column | Description |
|---|---|
| `cell_type` | the cell's type from surface-protein gating - see *The label* above |
