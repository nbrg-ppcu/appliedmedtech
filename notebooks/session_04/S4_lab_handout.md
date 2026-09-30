# Would you approve it?

---

## Where we are

Last week you built a model. It looks at a person's gut bacteria and blood results,
and it predicts whether they will be diagnosed with type 2 diabetes within one year.
The AUC was about 0.79.

I have a question I cannot answer with that number: **could this be a test?**

Not a result in a paper. A real test. A doctor orders it, a laboratory runs it, a
result comes back, and somebody acts on that result. That is very different from "a
model with a good AUC". The difference between the two is what we do this afternoon.

## Four questions

Anyone who evaluates a diagnostic test asks four questions. You heard them this
morning. Keep them in mind all afternoon.

- In whom will it be used?
- Where is the line drawn?
- How often is a positive result actually true?
- Compared to what?

The same model can give good answers in one place and bad answers in another. That
is not a problem with the model. It is the reason we have to ask.

## The files

**`t2d_microbiome_clean.csv`** - last week's cohort, already cleaned for you. 700
people. This is where the model was built.

**`t2d_screening.csv`** - new. 2,000 adults with prediabetes, found at routine GP
screening and followed for one year. This is where a test would actually be used.
It has the same columns as the cleaned file, so your model runs on it directly.

**`ed_deterioration.csv`** - the emergency department table from two weeks ago, used
in one short section at the end.

The data dictionaries are in the same folder.

## How the afternoon runs

1. **Your card.** Each pair gets one intended-use card, below. It says how the test
   would be used and what a mistake would cost.
2. **The Design Sheet**, in pen, for your card - not for the study. Fill it in once,
   before you start the notebook. Where nothing has changed since last week, "same
   as last week" is a fine answer. Think hardest about the decision, the baseline
   and the metric.
3. **The notebook.** One model, fitted once. After that, every section changes
   something the model did not control, and you watch what happens to the numbers.
4. **The memo**, at the end of the notebook. Give a verdict for your card:
   **approve, approve with a condition, or reject.**
5. **Upload the notebook to Moodle.** Hand in the Design Sheet.

At report-back the four verdicts go up side by side. Same model, four cards. Let us
see if we get four different answers.

---

## Intended-use cards

**Card 1 · Population screening.**
Offered to every adult over 50 at their GP. About 1 in 200 will progress within a
year. A false alarm means an extra clinic visit and a worried month. A miss means a
person who never knew.

**Card 2 · Triage in a diabetes clinic.**
Used on people who were already referred for raised blood sugar. About 1 in 10 will
progress. A miss means a patient sent home who should have been monitored. A false
alarm means a monitoring appointment used for someone who did not need it.

**Card 3 · Monitoring known prediabetics.**
Used on the same kind of people the model was built on. About 1 in 4 progress. Both
kinds of mistake cost something, and neither is much worse than the other.

**Card 4 · Research biomarker discovery.**
Not used on anybody. The only question is whether the gut bacteria carry a signal
that is worth studying further.

---

> **Note.** All three datasets are synthetic. They were built to behave like a study,
> a screening programme and an emergency department of this kind, but no real person
> is in them and no number here is a scientific result.
