# ed_deterioration — data dictionary

**Extract:** ED Deterioration Cohort · three sites · admissions 2023-01 to 2024-10
**Dictionary version:** 1.2 · prepared 2023-02 · *not revised since*
**Contact:** clinical data office

> **These data are synthetic.** They were generated for teaching

---

## What one row is

**One row is one emergency department admission.** Each admission has its own
`admission_id`. A patient may appear more than once.

## The outcome

`deterioration_24h` — the patient was transferred to intensive care, started on
vasopressors, or died, within 24 hours of arrival.

---

## Columns

Ranges below are **clinical reference ranges for adults**, not descriptions of
this file. Values outside them are possible and often clinically meaningful.

### Identifiers and administrative

| Column | Type | Description | Expected values |
|---|---|---|---|
| `admission_id` | text | Unique key for this admission | `ADM-#####` |
| `patient_id` | text | Patient identifier | `PT-####` |
| `site` | text | Hospital site | `A`, `B`, `C` |
| `admit_date` | date | Date of arrival | 2023-01-01 to 2024-10-31 |
| `consent_version` | text | Research consent form in force | `v2.1` |

### Demographics and anthropometry

| Column | Type | Units | Description | Reference range |
|---|---|---|---|---|
| `age` | integer | years | Age at admission | 18–98 (adults only) |
| `sex` | text | — | Recorded sex | `M`, `F` |
| `height_cm` | decimal | cm | Standing height | 145–201 |
| `weight_kg` | decimal | kg | Body weight | 40–155 |
| `bmi` | decimal | kg/m² | Body mass index | 16.5–47 · underweight <18.5 · obese >30 |

### Vital signs at triage

| Column | Type | Units | Description | Normal adult range |
|---|---|---|---|---|
| `heart_rate` | integer | bpm | Pulse | **60–100** · bradycardia <60 · tachycardia >100 |
| `resp_rate` | integer | breaths/min | Respiratory rate | **12–20** · tachypnoea >20 · a strong early warning sign |
| `sbp` | integer | mmHg | Systolic blood pressure | **90–140** · hypotension <90 |
| `temp_c` | decimal | °C | Body temperature | **36.1–37.8** · fever >38.0 · hypothermia <36.0 |
| `spo2` | integer | % | Peripheral oxygen saturation | **94–100** · hypoxaemia <92 |

### Laboratory results

| Column | Type | Units | Description | Normal adult range |
|---|---|---|---|---|
| `wbc` | decimal | ×10⁹/L | White cell count | **4.0–11.0** · neutropenia <1.0 · marked leucocytosis >20 |
| `creatinine` | decimal | mg/dL | Serum creatinine, renal function | **0.6–1.2** · acute kidney injury raises this sharply |
| `crp` | decimal | mg/L | C-reactive protein, inflammation | **<5** · bacterial infection 50–200 · severe >200 |
| `lactate` | decimal | mmol/L | Serum lactate, tissue perfusion | **0.5–2.2** · >2 concerning · >4 severe |

### Comorbidity

| Column | Type | Description | Expected values |
|---|---|---|---|
| `comorbidity_count` | integer | Number of recorded chronic conditions | 0–6 |
| `charlson_band` | text | Comorbidity burden band | `low`, `medium`, `high` |

### Admission course

| Column | Type | Description | Expected values |
|---|---|---|---|
| `abx_escalation_flag` | integer | Antibiotic therapy was escalated during the admission | 0, 1 |
| `discharge_disposition` | text | Patient status at discharge | `home`, `ward`, `ICU`, `died` |

### Outcome

| Column | Type | Description | Expected values |
|---|---|---|---|
| `deterioration_24h` | integer | Deterioration within 24 hours of arrival | 0, 1 |

---

## Notes on use

- Laboratory tests are ordered at clinical discretion. Not every admission has
  every result.
- Sites operate independent laboratory information systems. Reporting
  conventions are harmonised centrally before export.
- `bmi` is calculated at the point of export.
- `charlson_band` is derived from `comorbidity_count`.

## Synthetic data notes

Two properties of this file differ deliberately from a real ED cohort, and both
matter if you are tempted to read anything clinical into it:

1. **The outcome rate is inflated.** Roughly 30% of admissions deteriorate here.
   In a real emergency department the 24-hour figure is nearer 1–5%. The rate was
   raised so that 640 admissions contain enough events to model at all.
2. **The cohort is small by design.** 640 admissions is a realistic size for a
   single-site pilot study and an unrealistic one for a deployed model.
