# Bayesian hierarchical occupancy models with heterogeneous detection

Code for the UCT Statistics Honours thesis *Bayesian Hierarchical Models for Occupancy with Heterogeneous Detection* (Liam Erasmus, supervised by Dr Murray Christian). The case study is the Cape Rockjumper (*Chaetops frenatus*) using SABAP2 citizen science data.

The models add a latent visit-level access state `a_ij` between occupancy and detection. A visit only records the species if the pentad is occupied, the visit reached suitable habitat, and the bird is detected there. Co-occurring indicator species inform access.

## Models (`stan/`)

| File | Model |
|---|---|
| `standard_occupancy_ragged.stan` | Standard occupancy model with ragged visit counts |
| `composite_indicator_model_ragged.stan` | Model B: composite indicator `s_ij` drives access, `theta_ij = expit(beta0 + beta1 s_ij)` |
| `two_indicator_model.stan` | Model A: two indicator species |

Latent states are summed out analytically. Models are fitted with `cmdstanr` (4 chains, 1000 warmup + 1000 sampling).

## Scripts (`R/`)

All scripts assume the working directory is the project root. They use relative paths only.

| Script | Purpose |
|---|---|
| `fit_standard_ragged.R` | Standard model on SABAP2 data |
| `fit_SABAP_ragged_composite.R` | Model B on SABAP2 data (original priors; edit the beta priors for the wide-prior fit) |
| `run_recovery_modelB_ragged.R` | 100-replicate simulation study for Model B using the real visit structure |
| `run_recovery_two_indicators_sim100.R` | Recovery study for Model A |
| `check_recovery_two_indicators.R` | Summaries of the Model A recovery study |
| `figures/plot_*.R` | Thesis figures (study region, identifiability ridge, occupancy map and visits curve) |

## Data

No data are included. SABAP2 data are available from the Southern African Bird Atlas Project (https://sabap2.birdmap.africa) and can be downloaded with the [ABAP R package](https://github.com/AfricaBirdData/ABAP). The scripts expect the cleaned object at `data/clean/cleaned_multispp_occ_data.RDS`. Data preparation code is available from the author on request.

## Requirements

R, CmdStan, and the R packages `cmdstanr`, `posterior`, `ggplot2`, `patchwork`, `sf` and `rnaturalearth`. Record your own `sessionInfo()` when reproducing results.

## Licence

Code is released under the MIT licence (see `LICENSE`).
