# t2d_microbiome_clean - data dictionary

**Study:** PREDIGUT · gut microbiome and progression to type 2 diabetes
**Design:** prospective cohort, one baseline sample, twelve-month follow-up
**Centres:** Aarhus · Leuven · Tartu · Uppsala · Porto
**This file:** the cleaned version of last week's `t2d_microbiome.csv`

> **These data are synthetic.** They were generated for teaching and contain no real
> participants. Do not cite any number in this file as a scientific result.

---

## What one row is

**One row is one participant.** 700 participants, 700 rows. Each gave one fasting blood
sample and one stool sample at enrolment, and was followed for twelve months.

## Who is in the study

Adults aged 35 to 79 with **prediabetes** at screening. Anyone already diagnosed with type 2
diabetes, or already taking a glucose-lowering medicine such as metformin, was excluded
before this file was produced.

## The outcome

`progressed_12m` - was this participant diagnosed with type 2 diabetes within twelve months
of enrolment? Diagnosis means HbA1c of 6.5% (48 mmol/mol) or above, or fasting plasma glucose
of 7.0 mmol/L or above, on a confirmatory test. Follow-up is complete for everyone.

## What was cleaned

These are last week's fixes, already applied. The column names now match
`t2d_screening.csv` exactly.

| Last week | This file | What changed |
|---|---|---|
| `bmi`, text, comma decimals at Leuven | `bmi` | numeric |
| `hba1c`, mmol/mol at Tartu, percent elsewhere | `hba1c_pct` | percent at every site |
| `sex`, six spellings | `sex_c` | `M` / `F` |
| `smoking_status`, seven spellings | `smoking_c` | `never` / `former` / `current` |

Nothing else was changed. No rows were removed.

---

## Columns

Ranges are clinical reference values for adults, not descriptions of this file.

### Identifiers and centre

| Column | Description | Values |
|---|---|---|
| `participant_id` | unique key | `P####` |
| `site` | recruiting centre | five centres |
| `extraction_kit` | DNA extraction protocol for the stool sample - one per centre | KA-100, KB-200 |
| `protocol_version` | study protocol in force | `MB-2.0` |

### Participant

| Column | Description | Values |
|---|---|---|
| `age` | age at enrolment, years | 35–79 |
| `sex_c` | recorded sex | M, F |
| `smoking_c` | smoking status | never, former, current |
| `family_history_t2d` | a parent or sibling with type 2 diabetes | 0, 1 |
| `physical_activity` | self-reported activity level | low, moderate, high |
| `antibiotics_last_3m` | an antibiotic course in the three months before the sample | 0, 1 |

### Body measurements and blood chemistry

All blood drawn fasting, at enrolment.

| Column | Units | Reference |
|---|---|---|
| `bmi` | kg/m² | 18.5–25 normal · over 30 obese |
| `waist_cm` | cm | central obesity above 102 (men) or 88 (women) |
| `hba1c_pct` | % | normal below 5.7 · prediabetes 5.7–6.4 · diabetes 6.5 and above |
| `fasting_glucose_mmol_l` | mmol/L | normal below 5.6 · prediabetes 5.6–6.9 · diabetes 7.0 and above |
| `fasting_insulin_mu_l` | mU/L | 2–25 |
| `homa_ir` | index | insulin resistance, calculated at export from fasting glucose and insulin |
| `hdl_mmol_l` | mmol/L | desirable above 1.0 (men) or 1.3 (women) |
| `diet_fibre_score` | g/day | from a food questionnaire; 30 g/day recommended |

Not every assay is run on every sample, and not every questionnaire is returned.

### Gut microbiome

Ten genera from 16S sequencing of the stool sample, **centred log-ratio (CLR)
transformed** by the sequencing core. **The transform is already applied. These columns
share a common scale and need no further standardisation.** Typical values lie between
−4 and 4.

| Column | Genus | Reported association with type 2 diabetes |
|---|---|---|
| `clr_Faecalibacterium` | *Faecalibacterium* | butyrate producer · reported lower |
| `clr_Roseburia` | *Roseburia* | butyrate producer · reported lower |
| `clr_Eubacterium` | *Eubacterium* | butyrate producer · reported lower |
| `clr_Bifidobacterium` | *Bifidobacterium* | reported lower |
| `clr_Akkermansia` | *Akkermansia* | mucin degrader · reported lower |
| `clr_Escherichia` | *Escherichia* | reported higher |
| `clr_Lactobacillus` | *Lactobacillus* | reported higher |
| `clr_Streptococcus` | *Streptococcus* | reported higher |
| `clr_Bacteroides` | *Bacteroides* | reported associations inconsistent |
| `clr_Prevotella` | *Prevotella* | strongly diet-related · reported associations inconsistent |

Published associations are modest and have not replicated well between cohorts. Recovery of
some genera - *Bifidobacterium* in particular - depends on the DNA extraction protocol.

### Study record

Kept in the participant's study record and entered when the event happens. One value per
participant, not measured at the enrolment visit.

| Column | Description | Values |
|---|---|---|
| `metformin_initiated` | metformin was started for this participant | 0, 1 |
| `clinic_referral` | a referral was made for this participant | none, dietitian, diabetes service |

### Outcome

| Column | Description | Values |
|---|---|---|
| `progressed_12m` | diagnosed within twelve months of enrolment | 0, 1 |
