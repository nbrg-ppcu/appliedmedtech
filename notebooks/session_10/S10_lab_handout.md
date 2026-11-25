# Protein language models

---

## Where we are

We have used numerical features to describe biological samples and explored dimension reduction.
Today, a pretrained protein language model supplies the representation. You will inspect what goes
into ESM-2, what comes out, and what its predictions can tell you about an amino-acid replacement.

## The question

A protein is given as a sequence of amino-acid letters. Can a model trained to recover missing
letters help us understand that sequence? We will ask three questions:

1. **How does ESM-2 turn a protein sequence into numerical vectors?**
2. **Do proteins with similar locations in the cell appear near each other in an embedding map?**
3. **Which amino-acid replacements does the model prefer, and what does that preference establish?**

The model is already trained. In this notebook you use its outputs without training a classifier or
fine-tuning its parameters.

## The data

Each row of the localization dataset is one protein, with its amino-acid sequence and a location
label. A short fragment makes the tokens easy to inspect; a separate 108-residue example is used for
sequence and mutation scoring. No prior knowledge of amino-acid codes is needed: the notebook starts
with a reference table.

| Resource | What it is |
|---|---|
| `neuralbioinfo/Course-ProteinLoc` | protein sequences, localization labels and a test flag; downloaded from Hugging Face |
| `facebook/esm2_t6_8M_UR50D` | the pretrained ESM-2 checkpoint and tokenizer used throughout |
| [protein_DICTIONARY.md](data/protein_DICTIONARY.md) | the dataset columns, selection rules and output arrays: **read it first** |

Only `Sequence` is passed to the model. Localization labels color the plots. The full sequence header
contains label information, so it must not become an input feature. Rows marked `Is Test=True` are
reserved for later evaluation.

## What you will do, and why

**Inspect the model.** Find the checkpoint specifications, tokenize a fragment, and identify the
boundary, padding and mask tokens. Token IDs are indices, not amino-acid similarity scores.

**Read the outputs.** Compare logits with probabilities and last-layer embeddings. Mask a residue
before asking whether the model can recover it; a visible residue is not a prediction test.

**Represent whole proteins.** Average the last-layer vectors over amino-acid positions, excluding
padding and boundary tokens. Process batches and check that each vector stays with its label.

**Make a map.** Project the protein vectors with UMAP and compare neighborhood settings. Ask which
observations persist and whether apparent separation is enough to demonstrate predictive accuracy.

**Score the sequence.** Hide one residue at a time in separate copies and recover its original-residue
log probability. Sum these values for pseudo-log-likelihood (PLL), then compute the per-residue mean.

**Compare replacements.** At each masked site, subtract the original residue's log probability from
each candidate replacement's log probability. Read the heatmap and distinguish this single-site score
from a full mutant-versus-original PLL difference.

## Questions you will be able to answer

- What is the difference between a token ID, a logit, a probability and an embedding?
- Which positions should contribute to a protein's mean embedding, and what does truncation discard?
- Why does a separated UMAP plot not establish localization prediction accuracy?
- Why mask one residue per copy, and how do protein positions relate to tensor indices?
- What do summed PLL, mean PLL and a single-site replacement score measure?
- Why does a positive replacement score not prove improved protein function?

## What to hand in

- the **completed notebook**, with your code and short written answers to the 12 numbered exercises,
  uploaded to Moodle

---

> **Note.** The localization map is exploratory and the mutation scores are model predictions.
> Biological claims require independent evidence; this notebook does not measure protein function.
