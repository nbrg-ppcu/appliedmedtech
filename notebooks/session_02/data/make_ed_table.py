"""
Generate ed_deterioration.csv — the S2 clinical table.

Scenario: adults presenting to an emergency department. Outcome is clinical
deterioration within 24 hours of arrival (ICU transfer, vasopressors, or death).

Planted pathologies, in the order the nine-check audit finds them:
  shape          n is small (640 rows) -> noise floor around +-0.05
  types          admit_date is a string; several numerics read as object
  missingness    lactate 68% missing, NOT at random (ordered when worried)
  cardinality    sex has six spellings for two categories
  duplicates     52 patients appear twice (readmission) -> rows not independent
  constants      consent_version is a single value
  suspicious     admission_id, patient_id, admit_date
  derived        bmi is a function of height_cm and weight_kg
  units          creatinine is umol/L at site C, mg/dL at sites A and B

  TRAP 1 (obvious) discharge_disposition  -- recorded at discharge
  TRAP 2 (subtle)  abx_escalation_flag    -- escalation happens BECAUSE of decline
"""

import numpy as np
import pandas as pd

RNG = np.random.default_rng(20260916)
N = 640


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


# ---------------------------------------------------------------- true state
age = np.clip(RNG.normal(66, 17, N), 18, 98).round().astype(int)
comorbidity_count = RNG.poisson(1.6, N).clip(0, 6)

# Latent severity drives both the physiology and the outcome.
severity = (
    0.030 * (age - 66)
    + 0.35 * (comorbidity_count - 1.6)
    + RNG.normal(0, 1.0, N)
)

# Clinician gestalt: the thing the nurse noticed that never became a column.
gestalt = RNG.normal(0, 1.0, N)

resp_rate = np.clip(RNG.normal(19 + 2.3 * severity, 2.9, N), 8, 46).round().astype(int)
heart_rate = np.clip(RNG.normal(88 + 7.5 * severity, 13, N), 40, 175).round().astype(int)
sbp = np.clip(RNG.normal(132 - 9.0 * severity, 16, N), 62, 205).round().astype(int)
temp_c = np.round(np.clip(RNG.normal(37.0 + 0.30 * severity, 0.75, N), 34.5, 41.0), 1)
spo2 = np.clip(RNG.normal(96.5 - 1.7 * severity, 2.1, N), 74, 100).round().astype(int)
wbc = np.round(np.clip(RNG.normal(9.5 + 2.2 * severity, 3.4, N), 0.8, 42), 1)
lactate_true = np.round(np.clip(RNG.normal(1.6 + 0.55 * severity, 0.75, N), 0.3, 12), 1)

# ------------------------------------------------------------------- outcome
logit = -1.55 + 1.45 * severity + 0.40 * (lactate_true - 1.6) + 1.45 * gestalt
y = RNG.binomial(1, sigmoid(logit))

# ------------------------------------------------- informative missingness
# A lactate is drawn when the clinician is already concerned.
p_measured = sigmoid(-0.95 + 0.70 * severity + 2.45 * gestalt)
lactate_measured = RNG.binomial(1, p_measured).astype(bool)
lactate = np.where(lactate_measured, lactate_true, np.nan)

# CRP: missing largely at random, mild signal.
crp_true = np.round(np.clip(RNG.normal(48 + 22 * severity, 34, N), 0.4, 400), 1)
crp = np.where(RNG.random(N) < 0.55, crp_true, np.nan)

# --------------------------------------------------------- site and units
site = RNG.choice(["A", "B", "C"], N, p=[0.45, 0.35, 0.20])
creat_mgdl = np.round(np.clip(RNG.normal(1.02 + 0.20 * severity, 0.34, N), 0.3, 6.5), 2)
creatinine = np.where(site == "C", np.round(creat_mgdl * 88.4, 1), creat_mgdl)
creatinine = np.where(RNG.random(N) < 0.12, np.nan, creatinine)

# ------------------------------------------------------- messy categoricals
sex_true = RNG.choice(["M", "F"], N)
sex = np.array([
    RNG.choice(["M", "Male", "m"], p=[0.62, 0.26, 0.12]) if s == "M"
    else RNG.choice(["F", "Female", "f"], p=[0.60, 0.28, 0.12])
    for s in sex_true
])

charlson_band = np.where(comorbidity_count <= 1, "low",
                  np.where(comorbidity_count <= 3, "medium", "high"))

# --------------------------------------------------------- derived columns
# Body habitus must be internally coherent: height depends on sex, and weight
# is generated FROM a plausible BMI rather than drawn independently. Drawing
# height and weight separately produces 188 cm adults weighing 43 kg.
height_cm = np.where(
    sex_true == "M",
    RNG.normal(178.0, 7.0, N),
    RNG.normal(165.0, 6.5, N),
)
height_cm = np.round(np.clip(height_cm, 145, 201), 1)
bmi_true = np.round(np.clip(RNG.normal(27.0, 4.8, N), 16.5, 47), 1)
weight_kg = np.round(bmi_true * (height_cm / 100.0) ** 2, 1)

# Height is missing for ~8% of admissions, which also removes the derived BMI.
height_cm = np.where(RNG.random(N) < 0.08, np.nan, height_cm)
bmi = np.round(weight_kg / (height_cm / 100.0) ** 2, 1)

# ------------------------------------------------------------------- TRAPS
# TRAP 2 (subtle): antibiotics escalated because the patient was declining.
abx_escalation_flag = np.where(
    y == 1, RNG.binomial(1, 0.88, N), RNG.binomial(1, 0.09, N)
)

# TRAP 1 (obvious): recorded when the patient leaves.
disposition = np.empty(N, dtype=object)
for i in range(N):
    if y[i] == 1:
        disposition[i] = RNG.choice(["ICU", "died", "ward"], p=[0.62, 0.17, 0.21])
    else:
        disposition[i] = RNG.choice(["home", "ward", "ICU"], p=[0.68, 0.30, 0.02])

# ------------------------------------------------------------- identifiers
admit_date = pd.to_datetime("2023-01-01") + pd.to_timedelta(
    RNG.integers(0, 640, N), unit="D"
)
admission_id = [f"ADM-{i:05d}" for i in range(1, N + 1)]

# 52 patients are readmissions: rows are not independent.
n_repeat = 52
patient_numbers = list(range(1, N - n_repeat + 1))
repeat_of = RNG.choice(patient_numbers, n_repeat, replace=False)
all_patients = patient_numbers + list(repeat_of)
RNG.shuffle(all_patients)
patient_id = [f"PT-{p:04d}" for p in all_patients]

# --------------------------------------------------------------- assemble
df = pd.DataFrame({
    "admission_id": admission_id,
    "patient_id": patient_id,
    "site": site,
    "admit_date": admit_date.strftime("%Y-%m-%d"),
    "age": age,
    "sex": sex,
    "height_cm": height_cm,
    "weight_kg": weight_kg,
    "bmi": bmi,
    "heart_rate": heart_rate,
    "resp_rate": resp_rate,
    "sbp": sbp,
    "temp_c": temp_c,
    "spo2": spo2,
    "wbc": wbc,
    "creatinine": creatinine,
    "crp": crp,
    "lactate": lactate,
    "comorbidity_count": comorbidity_count,
    "charlson_band": charlson_band,
    "consent_version": "v2.1",
    "abx_escalation_flag": abx_escalation_flag,
    "discharge_disposition": disposition,
    "deterioration_24h": y,
})

# Site B's export writes a comma decimal, so the whole column arrives as text.
df["temp_c"] = [
    f"{t:.1f}".replace(".", ",") if s == "B" else f"{t:.1f}"
    for t, s in zip(df["temp_c"].astype(float), df["site"])
]


# ---------------------------------------------------------------------------
# Readmissions must behave like the same person twice, not like two strangers
# sharing an ID. Without this, a repeated patient has a different sex and
# height at each visit, and a patient-grouped split shows no effect at all.
# ---------------------------------------------------------------------------
STABLE = ["sex", "height_cm", "weight_kg", "bmi", "comorbidity_count", "charlson_band"]
CARRY  = ["heart_rate", "resp_rate", "sbp", "spo2", "wbc", "crp", "lactate", "creatinine"]

order = df.groupby("patient_id").cumcount()
first_of = df[order == 0].set_index("patient_id")

for idx in df.index[order > 0]:
    pid = df.at[idx, "patient_id"]
    base = first_of.loc[pid]
    for col in STABLE:                       # identity does not change
        df.at[idx, col] = base[col]
    df.at[idx, "age"] = int(np.clip(base["age"] + RNG.integers(0, 2), 18, 98))
    for col in CARRY:                        # physiology partly persists
        if pd.notna(base[col]) and pd.notna(df.at[idx, col]):
            blended = 0.65 * base[col] + 0.35 * df.at[idx, col]
            df.at[idx, col] = (int(round(blended))
                               if pd.api.types.is_integer_dtype(df[col])
                               else round(blended, 2))
        elif pd.isna(base[col]):
            df.at[idx, col] = np.nan         # the same test goes unordered again
    if RNG.random() < 0.75:                  # outcome repeats more often than not
        df.at[idx, "deterioration_24h"] = base["deterioration_24h"]
        df.at[idx, "abx_escalation_flag"] = base["abx_escalation_flag"]
        df.at[idx, "discharge_disposition"] = base["discharge_disposition"]

df.to_csv("/home/claude/ed_deterioration.csv", index=False)

print("rows, cols:", df.shape)
print("prevalence:", round(df.deterioration_24h.mean(), 3))
print("lactate missing:", round(df.lactate.isna().mean(), 3))
print("crp missing:", round(df.crp.isna().mean(), 3))
print("unique patients:", df.patient_id.nunique())
print("sex levels:", sorted(df.sex.unique()))
print("dispositions:", sorted(df.discharge_disposition.unique()))
