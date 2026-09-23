# t2d_microbiome - data dictionary

**Study:** PREDIGUT · gut microbiome and progression to type 2 diabetes
**Design:** prospective cohort, single baseline sample, 12-month follow-up
**Centres:** Aarhus (177) · Leuven (163) · Tartu (131) · Uppsala (135) · Porto (94)
**Dictionary version:** 2.0 · prepared before first data lock · *not revised since*
**Contact:** study data management

> **These data are synthetic.** They were generated for teaching and contain no
> real participants. Do not cite any number in this file as a scientific finding.

---

## What one row is

**One row is one participant.** 700 participants, 700 rows. Each gave **one**
fasting blood draw and **one** stool sample at enrolment, and was then followed for
twelve months. There are no repeat visits in this extract.

## Who is in the study

Adults aged 35 to 79 with **prediabetes** at screening - fasting glucose or HbA1c
above the normal range and below the diagnostic threshold for type 2 diabetes.

**Excluded at screening:** anyone already meeting diagnostic criteria for type 2
diabetes, and anyone already taking a glucose-lowering medication, metformin
included. **No participant in this file was diabetic or on treatment when their
sample was taken.** The exclusions were applied before the extract was produced.

## The outcome

`progressed_12m` - did this participant meet diagnostic criteria for type 2
diabetes during the twelve months after enrolment? Criteria are HbA1c of 6.5%
(48 mmol/mol) or above, or fasting plasma glucose of 7.0 mmol/L or above, on a
confirmatory test. Follow-up is complete for every participant.

---

## Columns

Ranges below are **clinical reference values for adults**, not descriptions of this
file.

### Identifiers and centre

| Column | Type | Description | Values |
|---|---|---|---|
| `participant_id` | text | Unique key | `P####` |
| `site` | text | Recruiting centre | Aarhus, Leuven, Tartu, Uppsala, Porto |
| `extraction_kit` | text | DNA extraction protocol used for the stool sample | KA-100, KB-200 |
| `protocol_version` | text | Study protocol in force | `MB-2.0` |

### Participant

| Column | Type | Units | Description | Reference |
|---|---|---|---|---|
| `age` | integer | years | Age at enrolment | 35-79 |
| `sex` | text | — | Recorded sex | M, F |
| `smoking_status` | text | — | Smoking status at enrolment | never, former, current |
| `family_history_t2d` | integer | — | First-degree relative with type 2 diabetes | 0, 1 |
| `physical_activity` | text | — | Self-reported activity level | low, moderate, high |
| `antibiotics_last_3m` | integer | — | Systemic antibiotic course in the three months before the sample | 0, 1 |

### Anthropometry

| Column | Type | Units | Description | Reference |
|---|---|---|---|---|
| `bmi` | decimal | kg/m² | Body mass index | **18.5-25 normal** · overweight >25 · obese >30 |
| `waist_cm` | decimal | cm | Waist circumference | central obesity above 102 (men) or 88 (women) |

### Blood chemistry

All drawn fasting, at enrolment.

| Column | Type | Units | Description | Reference |
|---|---|---|---|---|
| `hba1c` | decimal | % | Glycated haemoglobin, average glycaemia over ~3 months | **normal <5.7** · prediabetes 5.7-6.4 · diabetes ≥6.5 |
| `fasting_glucose_mmol_l` | decimal | mmol/L | Fasting plasma glucose | **normal <5.6** · prediabetes 5.6-6.9 · diabetes ≥7.0 |
| `fasting_insulin_mu_l` | decimal | mU/L | Fasting insulin | **2-25** |
| `homa_ir` | decimal | index | Insulin resistance index, calculated at export from fasting glucose and fasting insulin | **<2 normal** · >2.5 suggests resistance |
| `hdl_mmol_l` | decimal | mmol/L | HDL cholesterol | desirable above 1.0 (men) or 1.3 (women) |

### Diet

| Column | Type | Units | Description | Reference |
|---|---|---|---|---|
| `diet_fibre_score` | decimal | g/day | Dietary fibre, from a self-completed questionnaire returned by post | 30 g/day recommended |

### Gut microbiome

Ten genera from 16S profiling of the stool sample. Values are **centred log-ratio
(CLR) transformed** by the sequencing core: each is the log abundance of that genus
relative to the geometric mean of the sample.

> **The transform has already been applied. These ten columns share a common scale,
> are centred near zero, and require no further standardisation.** Typical values
> fall between −4 and 4.

| Column | Genus | Reported association with type 2 diabetes |
|---|---|---|
| `clr_Faecalibacterium` | *Faecalibacterium* | butyrate producer · reported **lower** in T2D |
| `clr_Roseburia` | *Roseburia* | butyrate producer · reported **lower** |
| `clr_Eubacterium` | *Eubacterium* | butyrate producer · reported **lower** |
| `clr_Bifidobacterium` | *Bifidobacterium* | reported **lower** |
| `clr_Akkermansia` | *Akkermansia* | mucin degrader · reported **lower**, associated with better metabolic health |
| `clr_Escherichia` | *Escherichia* | reported **higher** |
| `clr_Lactobacillus` | *Lactobacillus* | reported **higher** |
| `clr_Streptococcus` | *Streptococcus* | reported **higher** |
| `clr_Bacteroides` | *Bacteroides* | abundant genus · reported associations inconsistent |
| `clr_Prevotella` | *Prevotella* | strongly diet-associated · reported associations inconsistent |

> Published associations are modest and have not replicated consistently between
> cohorts. A large reanalysis attributed a substantial part of the reported
> signature to **metformin treatment** rather than to diabetes itself. This study's
> eligibility rules were written with that finding in mind.

> **Extraction protocol.** Two centres use kit KB-200 and three use KA-100.
> Recovery of some genera - *Bifidobacterium* in particular, which has a tough cell
> wall - is known to vary between protocols.

### Study record

Fields maintained in the participant's **study record**, not measured at the
enrolment visit. A study record holds one value per participant, entered when the
event occurs.

| Column | Type | Description | Values |
|---|---|---|---|
| `metformin_initiated` | integer | Metformin was initiated for this participant | 0, 1 |
| `clinic_referral` | text | Referral raised for this participant | none, dietitian, diabetes service |

### Outcome

| Column | Type | Description | Values |
|---|---|---|---|
| `progressed_12m` | integer | Diagnosis within 12 months of enrolment | 0, 1 |

---

## Notes on use

- Centres operate independent laboratories. Reporting conventions are harmonised
  centrally before export.
- `homa_ir` is calculated at the point of export.
- Not every assay is run on every sample, and not every questionnaire is returned.
- One extraction protocol is used per centre.

## Synthetic data notes

Two properties differ deliberately from a real cohort of this kind:

1. **The progression rate is high.** About 28% of participants are diagnosed within
   twelve months. In a real prediabetes cohort the figure over that horizon is
   considerably lower. It was raised so the extract contains enough events to model.
2. **The cohort is small.** 700 participants is a plausible size for a single study
   and an implausible one for a deployable risk score. Expect any performance
   estimate from this file to move by around 0.1 depending on which participants
   land in the test set.
