# The coronary study: what is in the table



## The study

600 patients with diagnosed coronary artery disease and 600 people without it, from 2 cohorts.
**Blood was taken at enrolment, for the patients, after their diagnosis and after their treatment had
started.** A DNA sample was taken at the same time.

| Cohort | Patients | People without the disease (the controls) |
|---|---|---|
| **hospital** | patients of a cardiology clinic | **regular blood donors** of the same hospital - healthy adults under 60, non-smokers, BMI under 30 |
| **population** | patients drawn from the population | people drawn from the same population |

## The columns

| Column | Unit / values | What it is |
|---|---|---|
| `patient_id` | P0001 - P1200 | identifier |
| `cohort` | hospital, population | which cohort the person came from |
| `cad` | 0, 1 | **the label:** 1 = diagnosed coronary artery disease |
| `age` | years | age at enrolment |
| `female` | 0, 1 | sex |
| `bmi` | kg/m² | body-mass index - weight divided by height squared |
| `smoking` | 0, 1 | current smoker |
| `ldl_c` | mmol/L | LDL cholesterol: the cholesterol carried in LDL particles |
| `apob` | g/L | apolipoprotein B: 1 molecule sits on every LDL-type particle, so it counts the particles |
| `lpa` | mg/dL | lipoprotein(a): an LDL-like particle; its level is set mostly by the LPA gene |
| `hdl` | mmol/L | HDL cholesterol |
| `crp` | mg/L | C-reactive protein - made by the liver when there is inflammation |
| `troponin` | ng/L | cardiac troponin: a protein that leaks from damaged heart muscle |
| `ntprobnp` | pg/mL | NT-proBNP: a protein released by stretched or strained heart muscle |
| `hemolysis` | index | how many red blood cells broke in the sample tube; this can disturb some laboratory tests |
| `statin` | 0, 1 | takes a **statin**: a drug that lowers LDL cholesterol. Nearly every patient is prescribed one at diagnosis. |
| `g_pcsk9` | 0, 1 | carries a variant of the PCSK9 gene that lowers LDL cholesterol for life |
| `g_lpa` | 0, 1, 2 | number of LPA gene variants that raise lipoprotein(a) |
| `g_hdl` | 0, 1, 2 | number of copies of a variant that raises HDL cholesterol only |
| `protein_001`: `protein_200` | various | **200 plasma proteins measured on the same sample.** Their names are hidden for this session; units differ between proteins. |

All measurements come from the same blood sample, taken at enrolment. The 3 genetic variants were
read from the DNA sample.
