# BIOSTAT 707 Checkpoint 1: Cohort Characterization and EDA

## Data and prediction time

The data contain 4,000 adult ICU records from PhysioNet Challenge 2012 set-a. Each record includes general information at ICU admission and clinical measurements from the first 48 hours. The admission variables include age, gender, height, weight, and ICU type. The time-series variables include vital signs, laboratory tests, urine output, mechanical ventilation, and repeated weight measurements.

A future model can only use information available during these 48 hours. Hospital length of stay, survival time, and in-hospital death cannot be predictors because they happen after the prediction window. I used SAPS-I and SOFA only to describe the cohort in this checkpoint.

## Data preparation and quality checks

The script combines all patient files into `output/set-a_long.csv`, which has 1,753,980 rows and keeps the measurement time. It also creates `output/set-a_wide.csv`, which has one row for each ICU record. The wide table summarizes each time-series variable using count, first, last, minimum, maximum, and mean over the full 48 hours. I chose one 48-hour summary because it is simple for this checkpoint, but separate 0–24 and 24–48 hour summaries may be useful later.

PhysioNet uses `-1` for missing or unknown values, so I changed these values to missing when appropriate. The long table keeps the original value in `ValueRaw` and the analysis value in `Value`. This makes the cleaning process easier to check. I also used conservative screening limits for variables with clear definitions, such as GCS from 3 to 15, FiO2 from 0 to 1, pH from 6.5 to 8.0, height from 100 to 250 cm, and temperature from 20 to 45 degrees Celsius.

The range check found 16 invalid height values, 235 invalid weight values, 129 invalid temperature values, and 11 invalid pH values. Examples included a height of 431.8 cm, a temperature of -17.8 degrees Celsius, and a pH of 735. I treated these values as missing in the summaries. I did not automatically remove other high laboratory values because they may be real values from very sick patients. The full check is in `output/measurement_range_check.csv`.

## Table 1 and outcomes

Table 1 is saved in `output/table1.csv`. The cohort included 3,446 survivors and 554 patients who died in the hospital. Patients who died were older, with a median age of 74 years compared with 66 years among survivors. They also had higher SAPS-I and SOFA scores. Medical ICU admission was more common among patients who died, while the gender percentages were similar between the groups. Admission height had a large amount of missing data, with 1,910 records missing a valid value.

The in-hospital mortality rate was 13.85%. Hospital length of stay had a median of 10 days and an interquartile range of 6 to 17 days. SAPS-I was missing for 190 records, SOFA was missing for 144 records, and length of stay was missing for 60 records. The complete results are in `output/outcome_summary.csv`.

## Missingness

Missingness was very different across variables. HR, temperature, GCS, BUN, and creatinine were measured in about 98% of records. Urine output was measured in 97.1%. In comparison, Troponin-I was measured in only 5.1%, cholesterol in 7.6%, Troponin-T in 21.6%, and respiration rate in 27.5%. These results are shown in `output/variable_coverage.csv` and `output/variable_coverage.png`.

The patient-by-variable map in `output/missingness_map_wide.png` shows that some variables were often measured together. Liver-related tests such as AST, ALT, bilirubin, and ALP had similar patterns. PaO2, PaCO2, and pH also often appeared together. Patients who died had a median of 29 variables measured, compared with 27 among survivors. This may mean that sicker patients received more tests.

The temporal map in `output/missingness_temporal.png` shows how measurement availability changed during the first 48 hours. There were 956,153 valid measurements during hours 0–24 and 781,126 during hours 24–48. Heart rate, urine output, temperature, and blood pressure were recorded repeatedly. Lactate was available for 1,957 patients in the first 24 hours but only 1,077 in the second 24 hours. This shows why time should be included when studying missingness.

## Is missingness informative?

The results in `output/missingness_vs_death.csv` suggest that missingness contains information. The death rate was 27.8% when Troponin-I was measured and 13.1% when it was not measured. For lactate, the rates were 17.6% and 9.3%. However, respiration rate showed the opposite pattern: mortality was 7.9% when it was measured and 16.1% when it was not measured.

These results do not mean that measuring a variable causes mortality. The decision to order a test may depend on disease severity, ICU type, equipment, or clinical workflow. In a future model, I may include whether a variable was measured and how many times it was measured instead of using only an imputed value.

## AI-use statement

I used OpenAI Codex substantially for this checkpoint. Codex read the assignment and PhysioNet documentation, helped set up the R project, and created or edited `.gitignore`, `README.md`, `pixi.toml`, `pixi.lock`, and `checkpoint1.R`. It helped write the long and wide table code, Table 1, outcome summaries, quality checks, missingness analyses, and figures. It downloaded the public data only for local testing; the submitted script does not download data.

Codex also debugged and tested the workflow. The checks confirmed 4,000 unique records, 4,000 rows in the wide table, no duplicate RecordID values, and times between 0 and 48 hours. The tables and figures were inspected, and the script was tested with `pixi run --locked checkpoint1` in a clean temporary directory.

The main advantages were saving time on repetitive code, debugging, and pixi setup. The disadvantage is that working AI-generated code can still make a wrong statistical or clinical decision. I am responsible for reading the final script, understanding the cleaning choices, and being able to explain the results.
