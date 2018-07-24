# Paper Repository

This repository contains the code to exactly reproduce all results and figures in our paper: [AdaPT: An interactive procedure for multiple testing with side information](https://arxiv.org/abs/1609.06035). 
After the acceptance of the paper, we developed an R package [`adaptMT`](https://github.com/lihualei71/adaptMT) that significantly improves the implementation and enriches the flexibility of AdaPT. Unless the reader just aims at reproducing the results in the paper, we recommend checking out our [new repository](https://github.com/lihualei71/adaptPaper) which re-implements all the examples in the paper, though producing slightly different yet improved results.

## R files
The folder `R/` contains all R files:

- `AdaPT.R` implements AdaPT as well as its variants (using GLM, GAM and glmnet),
- `AdaPT_mix_model/*.R` implements the EM pipeline to estimate local FDR. `AdaPT_mix_model.R` sources all files in `AdaPT_mix_model/`,
- `AdaPT_mix_summary.R` and `AdaPT_mix_plot.R` provides auxiliary functions associated with AdaPT,
- `useful_functions.R` contains helper functions in developing AdaPT,
- `*_expr.R` contains the code to implement all examples: `GEOquery_expr.R`, `GEOquery_random_expr.R`, `Bottomly_expr.R`, `airway_expr.R`, `pasilla_expr.R`, `proteomics_expr.R`, `simul1.R` and `simul2.R`,
- `AdaPT_expr_plot.R` produces all figures in the paper,
- Other files provide auxiliary functions.

## Replicating the experiments

- `GEOquery_expr.R`, `Bottomly_expr.R`, `airway_expr.R`, `pasilla_expr.R` and `proteomics_expr.R` can be executed in a laptop,
- `GEOquery_random_expr.R`, `simul1.R` and `simul2.R` must be executed on a cluster. To do so, submit `bash/gen_GEOquery_jobs.sh`, `bash/gen_jobs_simul1.sh` and `bash/gen_jobs_simul2.sh` to a cluster, respectively.
