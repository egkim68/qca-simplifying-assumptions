# How Reliable Are Simplifying Assumptions in Qualitative Comparative Analysis?

Code and data to reproduce the paper. Base R 4.3 or later; no packages needed.

```bash
Rscript simulate.R main        # 220,000 studies -> raw/main/       (several hours)
Rscript simulate.R extension   # 96,000 studies  -> raw/extension/
Rscript analyze.R              # all results     -> results/
Rscript figures.R              # Figures 1-4     -> figures/
```

- `simulate.R`: true functions, cases, truth tables, minimization, safeguards and scoring of each study.
- `analyze.R`: every number, table and figure input in the paper.
- `figures.R`: the figures.
- `data/`: the six published data sets.
- `results/`: the output of `analyze.R`, one file per table and figure (`table2.csv`, `tableS1.csv`, `figure2.csv`, ...) plus `text_numbers.csv` for numbers quoted in the text.

Each study has its own random-number stream (seed 60261110 or 70261110 plus the cell number), so any cell can be rerun alone with identical results.

MIT License.
