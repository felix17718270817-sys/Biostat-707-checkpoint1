# BIOSTAT 707 Checkpoint 1

This project analyzes PhysioNet Challenge 2012 set-a with an R script. Put the 4,000 patient files in `data/set-a/` and put `Outcomes-a.txt` in `data/`. Run `pixi run --locked checkpoint1` from the repository root. The script creates the long and wide analysis tables, Table 1, the outcome summary, data-quality checks, missingness summaries, and missingness figures in `output/`.

On macOS, use a repository path without spaces. The conda-forge R launcher may fail before R starts when its environment path contains a space.
