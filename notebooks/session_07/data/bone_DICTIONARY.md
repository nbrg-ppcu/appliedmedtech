# Bone single-cell data - data dictionary

**Files:** `bone_cells.csv` · `bone_labels.csv`
**Technology:** 10x Genomics 3′ gene expression, chemistry v3, on cells released from bone and marrow
by enzymatic digestion

---

## What one row is

**One row is one cell** (or, rarely, two cells captured together - see *Doublets* below). The cells
come from human bone and bone marrow: bone-forming cells, cartilage cells, the support cells of the
marrow, blood vessels and blood cells. Each donor contributed between about 30 and 90 cells.

| File | What it is |
|---|---|
| `bone_cells.csv` | 1,000 cells from 16 donors (D01–D16), 4 facilities. Metadata and 62 genes. **No cell types.** |
| `bone_labels.csv` | the cell type of every cell, from the expert panel - **only for the last part of the lab** |

The two files share the column `cell_id`.

## The four facilities

Each facility sent the cells of 4 donors. The facilities are different hospitals, and their donors are
different kinds of people.

| Facility | Who | Why the bone was available | Bone | How the tissue was handled |
|---|---|---|---|---|
| **A** | university orthopaedic clinic | hip replacement for **osteoarthritis** (worn joint cartilage); donors 60–80 years | femoral head | processed **fresh**, within 1–2 hours |
| **B** | trauma centre | surgery after a **hip fracture**; all four donors were treated for an **osteoporotic fracture** (osteoporosis: fragile, thinned bone); donors 70–90 years | femoral head | processed **fresh**, within 1–3 hours |
| **C** | young-adult orthopaedics | corrective bone surgery (osteotomy) on **healthy** bone; donors 16–29 years | iliac crest (hip bone) or tibia (shin bone) | processed **fresh**, within 1 hour |
| **D** | tissue bank | **deceased** donors, **healthy** bone; donors 30–65 years | femoral shaft | collected 12–24 hours after death, **cryopreserved** (frozen), thawed for processing |

All four used the same chemistry and the same sequencing depth target, and the same digestion protocol.

**Condition and facility go together.** Every osteoarthritis donor comes from facility A, every
osteoporotic-fracture donor from facility B, and every healthy donor from facility C or D. So a
difference between the conditions is also a difference between facilities: in age, handling, and
anything else that differs between the hospitals. The data cannot separate the two.

## The labels (`bone_labels.csv`)

An expert panel assigned a type to every cell. They used **the whole transcriptome** (about 20,000
genes, not just the 62 in your file) together with the tissue sections of each sample. Their labels
come at three levels, from coarse to fine:

| Column | Levels |
|---|---|
| `compartment` | 3 groups, plus `doublet` |
| `lineage` | 7 groups, plus `doublet` |
| `cell_type` | 16 types, plus `doublet` |

Use this file only when the notebook tells you to. Until then the task is the one every new atlas
faces: there are no labels.

**Doublets.** Sometimes two cells are caught in the same droplet and are read as one. Such a "cell"
carries the genes of two cell types at once and usually has more UMIs than a single cell. The panel
labelled these `doublet` at every level.

---

## Columns of `bone_cells.csv`

### Metadata

| Column | Type | Description |
|---|---|---|
| `cell_id` | text | 10x cell barcode (16 bases and `-1`) joined to the donor ID, e.g. `TGTCCTATGACAATTA-1_D01`; unique |
| `donor_id` | text | the person the cell came from, D01–D16 |
| `facility` | A, B, C, D | the facility that collected and processed the sample - see above |
| `procedure` | text | why the bone was available: `hip replacement` · `hip fracture surgery` · `corrective osteotomy` · `tissue bank donor (deceased)` |
| `bone_site` | text | which bone the sample came from |
| `age` | integer | age of the donor in years |
| `sex` | F, M | recorded sex of the donor |
| `condition` | text | `osteoarthritis` · `osteoporotic fracture` · `healthy` |
| `handling` | text | `fresh` or `cryopreserved` |
| `ischaemia_h` | number | hours between the loss of blood supply to the bone and the start of processing |
| `n_counts` | integer | $s_c$ - total UMIs of the cell, over the whole transcriptome (definition below) |
| `n_genes` | integer | $n_c$ - number of genes detected in the cell, over the whole transcriptome |
| `pct_mito` | number | $m_c$ - percentage of the cell's UMIs from mitochondrial genes; high values often mean a damaged or stressed cell |

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

### Gene expression - 62 columns

Each value is the **shifted logarithm** of the count, scaled to a fixed total per cell:

$$
y_{gc} = \ln\!\left(1 + \frac{x_{gc}}{s_c}\, L\right), \qquad L = 10^{4}
$$

- $x_{gc}$ - UMI count of gene $g$ in cell $c$
- $s_c$ - total UMIs of cell $c$ over the whole transcriptome (`n_counts`)
- $L$ - a fixed target sum: $x_{gc}/s_c \cdot L$ is counts per 10,000 (CP10k)
- $\ln$ - the natural logarithm

This is `scanpy.pp.normalize_total(target_sum=1e4)` followed by `scanpy.pp.log1p`, and Seurat's
`LogNormalize` with `scale.factor = 10000`. Values range from 0 to about 7. All genes are on the same
log scale, so the notebook uses them as they are, without rescaling.

$y_{gc} = 0$ exactly when $x_{gc} = 0$: no transcript of the gene was captured in that cell. Most
values here are 0. A zero does not always mean the gene is off - single-cell sequencing captures only
a fraction of each cell's RNA.

| Gene | Full name | Commonly used as a marker of |
|---|---|---|
| PTPRC | protein tyrosine phosphatase receptor type C (CD45) | all blood (haematopoietic) cells |
| PECAM1 | platelet and endothelial cell adhesion molecule 1 (CD31) | endothelial cells |
| CDH5 | cadherin 5 (VE-cadherin) | endothelial cells |
| PDGFRB | platelet derived growth factor receptor beta | stromal (mesenchymal) cells |
| COL1A1 | collagen type I alpha 1 chain | stromal cells; very high in bone-forming cells |
| RUNX2 | RUNX family transcription factor 2 | the bone-forming lineage |
| SP7 | Sp7 transcription factor (osterix) | osteoblasts |
| LEPR | leptin receptor | CAR cells (CXCL12-abundant reticular cells, the marrow's support cells) |
| CXCL12 | C-X-C motif chemokine ligand 12 | CAR cells; the signal that keeps blood stem cells in the marrow |
| SOX9 | SRY-box transcription factor 9 | chondrocytes |
| ALPL | alkaline phosphatase | osteoprogenitors, osteoblasts, Osteo-CAR cells |
| BGLAP | bone gamma-carboxyglutamate protein (osteocalcin) | osteoblasts, osteocytes |
| IBSP | integrin binding sialoprotein | osteoblasts |
| SPP1 | secreted phosphoprotein 1 (osteopontin) | osteoblasts |
| SOST | sclerostin | osteocytes |
| DMP1 | dentin matrix acidic phosphoprotein 1 | osteocytes |
| PHEX | phosphate regulating endopeptidase X-linked | osteocytes |
| ADIPOQ | adiponectin | Adipo-CAR cells, fat cells |
| LPL | lipoprotein lipase | Adipo-CAR cells, fat cells |
| COL2A1 | collagen type II alpha 1 chain | chondrocytes |
| ACAN | aggrecan | chondrocytes |
| COL10A1 | collagen type X alpha 1 chain | maturing (hypertrophic) chondrocytes |
| GJA5 | gap junction protein alpha 5 | arterial endothelial cells |
| SEMA3G | semaphorin 3G | arterial endothelial cells |
| STAB2 | stabilin 2 | sinusoidal endothelial cells |
| FLT4 | fms related receptor tyrosine kinase 4 (VEGFR3) | sinusoidal endothelial cells |
| TYROBP | transmembrane immune signaling adaptor TYROBP | myeloid cells |
| CSF1R | colony stimulating factor 1 receptor | monocytes, macrophages, osteoclasts |
| LYZ | lysozyme | monocytes, macrophages |
| CD14 | CD14 molecule | classical monocytes |
| FCN1 | ficolin 1 | classical monocytes |
| C1QA | complement C1q A chain | macrophages |
| CD163 | CD163 molecule | macrophages |
| S100A8 | S100 calcium binding protein A8 | neutrophils, monocytes |
| S100A9 | S100 calcium binding protein A9 | neutrophils, monocytes |
| CSF3R | colony stimulating factor 3 receptor | neutrophils |
| ACP5 | acid phosphatase 5, tartrate resistant (TRAP) | osteoclasts |
| CTSK | cathepsin K | osteoclasts |
| MMP9 | matrix metallopeptidase 9 | osteoclasts |
| DCSTAMP | dendrocyte expressed seven transmembrane protein | osteoclasts (needed for their fusion) |
| CD79A | CD79a molecule | B cells |
| MS4A1 | membrane spanning 4-domains A1 (CD20) | B cells |
| CD3E | CD3 epsilon subunit of T-cell receptor complex | T cells |
| CD3D | CD3 delta subunit of T-cell receptor complex | T cells |
| IL7R | interleukin 7 receptor | T cells, early lymphoid cells |
| CD34 | CD34 molecule | haematopoietic stem and progenitor cells |
| KIT | KIT proto-oncogene, receptor tyrosine kinase | haematopoietic stem and progenitor cells |
| HBB | haemoglobin subunit beta | erythroid (red blood cell) progenitors |
| GYPA | glycophorin A | erythroid progenitors |
| ACTB | actin beta | — (in every cell) |
| B2M | beta-2-microglobulin | — (in every cell) |
| MALAT1 | metastasis associated lung adenocarcinoma transcript 1 (non-coding RNA) | — (in every cell) |
| FOS | Fos proto-oncogene, AP-1 transcription factor subunit | — (rises within minutes when a cell is stressed) |
| JUN | Jun proto-oncogene, AP-1 transcription factor subunit | — (rises within minutes when a cell is stressed) |
| HSPA1A | heat shock protein family A (Hsp70) member 1A | — (rises when a cell is stressed) |
| DUSP1 | dual specificity phosphatase 1 | — (rises when a cell is stressed) |
| XIST | X inactive specific transcript (non-coding RNA) | — (cells of female donors) |
| RPS4Y1 | ribosomal protein S4 Y-linked 1 | — (cells of male donors) |
| MKI67 | marker of proliferation Ki-67 | dividing cells |
| TOP2A | DNA topoisomerase II alpha | dividing cells |
| CDKN2A | cyclin dependent kinase inhibitor 2A (p16) | ageing (senescent) cells |
| IL1B | interleukin 1 beta | inflammation |

---

## Marker table - which genes point to which cell

The cell types a bone atlas usually finds, with the genes the expert panel looks for first. Use it to
**name your clusters**. A cluster can match more than one row, or none. The names in the first column are
the ones used in `bone_labels.csv`.

**CAR cells** are *CXCL12-abundant reticular cells*: net-like support cells of the marrow that make the
signal CXCL12, which keeps blood stem cells in place. **EC** stands for *endothelial cell*, the cells
that line a vessel.

| Cell | What it does | Look for |
|---|---|---|
| osteoprogenitor | early bone-forming cell | RUNX2, ALPL, COL1A1; SP7 low |
| osteoblast | builds new bone | RUNX2, SP7, BGLAP, IBSP, SPP1, COL1A1 high |
| osteocyte | an osteoblast walled into the bone; senses load | SOST, DMP1, PHEX, BGLAP |
| Adipo-CAR | marrow support cell, leaning towards fat | LEPR, CXCL12, ADIPOQ, LPL |
| Osteo-CAR | marrow support cell, leaning towards bone | LEPR, CXCL12, ALPL, SPP1 |
| chondrocyte | cartilage cell | SOX9, COL2A1, ACAN, COL10A1 |
| arterial EC | lines arteries | PECAM1, CDH5, GJA5, SEMA3G |
| sinusoidal EC | lines the marrow's thin-walled vessels | PECAM1, CDH5, STAB2, FLT4 |
| classical monocyte | circulating myeloid cell | PTPRC, TYROBP, CD14, LYZ, FCN1 |
| macrophage | tissue myeloid cell that eats debris | PTPRC, TYROBP, C1QA, CD163 |
| neutrophil | most common white blood cell | PTPRC, S100A8, S100A9, CSF3R |
| osteoclast | resorbs (dissolves) bone | PTPRC, TYROBP, CSF1R, ACP5, CTSK, MMP9, DCSTAMP |
| B cell | makes antibodies | PTPRC, CD79A, MS4A1 |
| T cell | cellular immunity | PTPRC, CD3E, CD3D, IL7R |
| HSPC | haematopoietic stem and progenitor cell; makes all blood cells | PTPRC, CD34, KIT; MKI67 |
| erythroid progenitor | becomes a red blood cell | PTPRC (low), HBB, GYPA; MKI67 |

---

> **Note.** The data are synthetic. They were built to behave like real bone single-cell data, but no
> real person is in them and no number in them is a scientific result.
