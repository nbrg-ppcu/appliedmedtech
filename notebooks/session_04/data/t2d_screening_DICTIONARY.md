# t2d_screening - data dictionary

**Programme:** regional GP prediabetes screening · harmonised export
**Dictionary version:** 1.0 · prepared with the export
**Contact:** screening programme data office

> **These data are synthetic**, generated for teaching. No real person is in them.

---

## What this file is

2,000 adults found to have prediabetes at routine GP screening, and followed for twelve
months. **One row is one person.** They were not recruited into a research study. They are
the kind of people a test like this would actually be used on.

The columns are exactly those of `t2d_microbiome_clean.csv`, in the same order, so a model
fitted on that file runs on this one without changes.

## How it compares with the study cohort

| | `t2d_microbiome_clean.csv` | `t2d_screening.csv` |
|---|---|---|
| How people were found | recruited into the PREDIGUT research study | routine GP screening |
| Rows | 700 | 2,000 |
| Diagnosed within twelve months | 28% | 6.9% |
| Eligibility | prediabetes, untreated | prediabetes, untreated |
| Export | cleaned for this course | harmonised by the programme's data office |

## The outcome

`progressed_12m` - diagnosis of type 2 diabetes within twelve months of the screening visit:
HbA1c of 6.5% or above, or fasting glucose of 7.0 mmol/L or above, on a confirmatory test.
Follow-up is complete.

## Columns

Same names, units and meanings as `t2d_microbiome_clean.csv` - see that dictionary for the
full list and the reference ranges. In short:

- identifiers and centre: `participant_id`, `site`, `extraction_kit`, `protocol_version`
- participant: `age`, `sex_c`, `smoking_c`, `family_history_t2d`, `physical_activity`,
  `antibiotics_last_3m`
- body and blood: `bmi`, `waist_cm`, `hba1c_pct`, `fasting_glucose_mmol_l`,
  `fasting_insulin_mu_l`, `homa_ir`, `hdl_mmol_l`, `diet_fibre_score`
- gut microbiome: the same ten `clr_` genera, CLR-transformed, no further scaling needed
- study record: `metformin_initiated`, `clinic_referral`
- outcome: `progressed_12m`

## Notes on use

- The programme's centres are the same five as the study's, with the same extraction kit
  at each.
- The two study-record fields are present for the same reason as in the study file.
