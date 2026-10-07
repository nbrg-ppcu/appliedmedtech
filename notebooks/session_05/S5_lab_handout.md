# Which cell is it?

---

## Where we are

The lecture was about one question: **is my test set new in the same way the real use will be
new?** In this lab you answer it for a real deployment, on single-cell data - and, unusually, you can
check your answer against reality.

## The question

A hospital biobank wants a classifier that labels the cell types in single-cell RNA-seq data from
its stored blood samples. Its samples are **cryopreserved**, **total PBMC**, from **patients the
model has never seen**.

We have data from two research labs. **How good will the classifier be in the biobank?**

Normally you could only estimate that from the training data. Here the biobank has labelled some of
its cells, so you can compare every estimate with the real performance - and find out which estimate
told the truth, and why the others did not.

## The data

Single cells from human blood - T cells, NK cells, B cells and monocytes. Each row is one cell, with
40 genes and some information about where the cell came from. The cell types were decided by surface
proteins, not by the RNA.

| File | What it is |
|---|---|
| `pbmc_train.csv` | 4,007 cells from 16 donors, processed by two research labs - the training data |
| `pbmc_external.csv` | 1,533 labelled cells from 6 new donors, processed by the biobank - the real performance |
| `pbmc_DICTIONARY.md` | what every column means, and how each lab handled its samples - **read it first** |

## What you will do, and why

**Look at the data.** Where did the training cells come from? Many problems in real data are visible
before any model is trained.

**Measure the model every way you can.** On the cells it was trained on, with random folds, and with
folds that hold out whole donors - because the biobank's patients are new. Then compare each estimate
with the biobank. Which one came closest?

**Find out why.** If an estimate was wrong, something in the training data made the problem look
easier than it is. You look for it: which genes differ between the labs, and what the model leans on.

**Repair it, and measure again.** Change the model, rebuild the same comparison, and see whether the
estimate now tells the truth - for the average, and for each cell type.

**Change the question.** The same cells, a new label: does the donor have **SLE** (systemic lupus
erythematosus, an autoimmune disease)? Every cell of a donor has the same answer. You see what that
does to the evaluation, and count how many samples you really have.

**Answer the question.** How good will the classifier be in the biobank, where can it be trusted,
and which estimate from the training data would have told you so?

## Questions you will be able to answer

- Why can a model look better in cross-validation than it will be in practice - even with the right
  split?
- Which estimate from the training data should you trust, and when?
- How can you tell, from the training data alone, whether a model learned what you think it learned?
- When is the number of cells not the number of samples?

These are the questions you will face in your own thesis, whatever the data.

## What to hand in

- the **Design Sheet**, filled in pen, in pairs, for the biobank question - before you open the
  notebook
- the **notebook**, with your conclusions and your answer to Task 10 - uploaded to Moodle

---

> **Note.** The data are synthetic. They were built to behave like real blood single-cell data, but
> no real person is in them and no number in them is a scientific result.
