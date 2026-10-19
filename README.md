# BIOSTAT 707 Checkpoints

This project analyzes PhysioNet Challenge 2012 set-a with an R script. Put the 4,000 patient files in `data/set-a/` and put `Outcomes-a.txt` in `data/`. Run `pixi run --locked checkpoint1` from the repository root. The script creates the long and wide analysis tables, Table 1, the outcome summary, data-quality checks, missingness summaries, and missingness figures in `output/`.

For Checkpoint 2, run `pixi run --locked checkpoint2`. This single command rebuilds the Checkpoint 1 tables, creates 48-hour predictors, runs five-fold outer and three-fold inner cross-validation, compares elastic net, a GAM, and random forest with SAPS-I and SOFA, and writes `output/checkpoint2.pdf`. The English report contains all figures, formatted tables, and interpretation. Allow about 3 minutes for Checkpoint 2 and 10-20 seconds for Checkpoint 1 on this laptop, after the first environment installation. The code reads set-a only. No data download or other command is needed during analysis.

Checkpoint 2 source files are `checkpoint2.R`, `R/cp2_models.R`, and `R/cp2_report.R`. The saved final pipeline is `output/cp2_final_model.rds`, and held-out fold predictions are in `output/cp2_oof_predictions.csv`. Both are local derived outputs for Checkpoint 3. The PDF report is tracked in Git; input data, other outputs, and the local environment are ignored. `cp2_predict_final()` in `R/cp2_models.R` applies the saved preprocessing and model to a wide table built with the same rules. CV audits, tuning results, coefficients, and a session record are also written under `output/`.

On macOS, use a repository path without spaces. The conda-forge R launcher may fail before R starts when its environment path contains a space.
