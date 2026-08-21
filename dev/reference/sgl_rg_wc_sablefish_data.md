# West Coast Sablefish data for single region case study

A dataset containing the necessary elements for the 2025 West Coast
sablefish case study, including model inputs (dimensions, biologicals,
catch, four survey indices, a recruitment index, and sexed and unsexed
compositions), the assessment's selectivity parameters and time blocks,
its maximum likelihood estimates, and the Stock Synthesis derived
quantities and likelihood components used for bridge verification.

## Usage

``` r
sgl_rg_wc_sablefish_data
```

## Format

A list with model inputs at the top level,
`sel_fish`/`sel_srv`/`sel_male` (the assessment's double normal
selectivity parameters, which of them it estimates, and which blocks
share one), `mle` (Stock Synthesis maximum likelihood estimates, read
from the parameter file at twelve significant digits), and `ss3` (Stock
Synthesis derived quantities and likelihood components for bridge
comparison)

## Source

2025 Pacific Fishery Management Council West Coast Sablefish Assessment,
Stock Synthesis 3.30.23
