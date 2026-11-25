# Protein localization data dictionary

Session 10 uses [neuralbioinfo/Course-ProteinLoc](https://huggingface.co/datasets/neuralbioinfo/Course-ProteinLoc).
The notebook downloads it with `load_dataset`. Inspect `dataset["train"].features`, example rows and
label counts before using it.

## One row, one protein

| Column | Meaning | Use in this notebook |
|---|---|---|
| `Protein ID` | Protein identifier | identify a row; not a model input |
| `Sequence` | Amino-acid sequence in one-letter codes | the only model input |
| `Localization` | Protein location in the cell | color the embedding map; not used to calculate UMAP coordinates |
| `Full Header` | Original sequence header, including label information | inspection only; using it as a feature would expose label information |
| `Is Test` | Whether the protein is marked for testing | exclude `True` rows from this session's embedding sample |

The Hugging Face split is named `train`, but it contains rows marked `Is Test=True`. The split name
alone is therefore not the evaluation boundary. This notebook filters with `not row["Is Test"]`,
shuffles with seed 42, and selects up to `N_PROTEINS = 1000` rows. It does not evaluate a classifier.

The notebook prints class counts for both the full split and the selected sample. Unequal counts
affect how prominent a class looks in the map. For later prediction tasks, investigate sequence
similarity between training and test proteins before claiming performance on unfamiliar proteins.

## From sequence to model input

| Array or setting | Meaning |
|---|---|
| `input_ids` | integer token IDs; each input also has boundary tokens |
| `attention_mask` | 1 for input tokens, including boundaries; 0 for padding |
| `special_tokens_mask` | marks special tokens to exclude from residue pooling |
| `MAX_LENGTH = 512` | the total token budget, including the two boundary tokens |
| `residue_mask` | positions with attention-mask value 1 that are not special tokens |

For dataset embeddings, proteins longer than 510 amino acids are truncated to their first 510
residues. Padding fills shorter inputs. The notebook reports the number of truncated proteins.
Neither padding nor boundary tokens contribute to the mean embedding.

## Arrays created during the session

Let `N` be the number of selected proteins, `T` the token count, `H` the model's hidden size, and `V`
its vocabulary size. Read `H` and `V` from `model.config` in Exercise 1.

| Variable | Shape or contents | Meaning |
|---|---|---|
| `outputs.logits` | batch × T × V | token scores before softmax |
| `last_hidden_state` | batch × T × H | contextual vectors from the final transformer layer |
| `embeddings` | N × H | one residue-mean vector per protein |
| `labels` | N | localization labels, in the same order as the vectors |
| `umap_embeddings` | N × 2 | visualization coordinates; labels are used only for coloring |
| `position_scores` | 108 rows: `position`, `original`, `log_probability` | recovery score for each original residue when that residue is masked |
| `pll` | scalar | sum of the 108 original-residue log probabilities |
| `mean_pll` | scalar | `pll / sequence_length` |
| `heatmap` | 20 × 108 | single-site replacement log-probability differences |

## Protein positions and replacement scores

For the unpadded scoring example, token 0 is the start token. Residue 1 is token 1, but its letter is
at Python string index 0. Row `position - 1` of the stored masked-site predictions corresponds to
that protein position.

Each masked copy hides exactly one residue. The heatmap compares each of the 20 standard amino acids
with the original residue using the same masked context. Columns are positions, starting at 1;
rows are candidate replacements in `ACDEFGHIKLMNPQRSTVWY` order. A self-replacement scores zero.
Positive values mean the model prefers the replacement. They do not measure functional improvement
and are not generally equal to rescoring the entire mutant sequence and subtracting the original PLL.
