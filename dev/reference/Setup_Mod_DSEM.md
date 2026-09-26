# Set up a dynamic structural equation model on any deviation process

Deviations and covariate series become the columns of one year by series
grid, linked by arrow and lag lines. A linked series takes its penalty
from the dsem density instead of SPoRC's own, through its map mirror.
There is no separate objective:
[`SPoRC_rtmb`](https://chengmatt.github.io/SPoRC/dev/reference/SPoRC_rtmb.md)
evaluates the density whenever a dsem is set up. Fit with
`random = c(<linked arrays>, "dsem_x")`.

## Usage

``` r
Setup_Mod_DSEM(
  input_list,
  dsem_arrows,
  dsem_data,
  dsem_processes = NULL,
  dsem_family = NULL,
  dsem_link = NULL,
  dsem_fixed_sd = NULL,
  dsem_mu_spec = "est",
  covs = NULL,
  dsem_delta0_spec = "none",
  mod_var_logscale = FALSE,
  dsem_variance = "conditional"
)
```

## Arguments

- input_list:

  List with `data`, `par` and `map`.

- dsem_arrows:

  Arrow lines, one per element of a character vector or one per line of
  a string. Series are named for their process and every index dim, for
  example `rec`, `rec_Pop_1_Region_2` or
  `NAA_Pop_1_Region_1_Seas_1_Age_3_Sex_1`, with the bare label used when
  every index dim has one level. Covariates are named by the columns of
  `dsem_data`. A name of `NA` with a value in the fourth field fixes
  that arrow, an sd on the natural scale; to fix a parameter that
  several arrows share, write the value on each of them.

- dsem_data:

  Data frame with a `year` column and one column per covariate, `NA`
  where a covariate is not observed. `NULL` for a dsem among deviation
  series alone.

- dsem_processes:

  Processes whose series the arrows may name, from `dsem_process_table`.
  `NULL` (default) takes the processes whose module declared `"dsem"`:
  `RecDevs_model` in
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md);
  `NAA_re`, `growth_tv_model` and `growth_semipar` in
  [`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md);
  `cont_vary_movement` in
  [`Setup_Mod_Movement`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md);
  `fish_q_model` in
  [`Setup_Mod_Fishsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md);
  and `srv_q_model` in
  [`Setup_Mod_Srvsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md).
  It is `"rec"` when none did.

  A declared process has to have every series with an estimated cell in
  the arrows, since its own penalty is off. A process linked without a
  declaration keeps its penalty on the cells left out, but its sigma
  cannot stay estimated once every cell is linked.

  A linked series covers every year of its array, so a cell mapped off
  (`NA` in the map) inside one is not left out: the dsem reads it at its
  fixed starting value and still evaluates its innovation. A recruitment
  series with such a mix is refused; for the other processes the mix is
  reported.

  Under CTMC movement the deviations are the year to year part of
  habitat preference, so a `preference_formula` term that varies over
  years is refused once a movement series is linked.

- dsem_family:

  Named character vector, one entry per covariate, the distribution of
  its observations given the grid cells. `"fixed"` (default) takes an
  observed year as the cell itself, known, and leaves a missing year
  latent. `"normal"` leaves every year latent, each observation normal
  about the cell with sd `exp(ln_dsem_obs_sd)`. `"bernoulli"` (or
  `"binomial"`) takes 0/1 observations with the cell the logit of the
  probability; `"poisson"` takes counts with the cell the log mean;
  `"gamma"` takes positive values with the cell the log mean and
  `exp(ln_dsem_obs_sd)` the CV, shape \\1/CV^2\\; `"lognormal"` takes
  positive values with the cell the log median and `exp(ln_dsem_obs_sd)`
  the log-scale sd; `"tweedie"` takes non-negative values including
  zeros with the cell the log mean, `exp(ln_dsem_obs_sd)` the dispersion
  and `logit_dsem_tweedie_p` the power as \\1 +
  \mathrm{plogis}(\cdot)\\, in (1, 2) and starting at 1.5; and
  `"gaussian_fixed_sd"` is normal about the cell with a known sd per
  observed year from `dsem_fixed_sd`. The names `"gaussian"` and
  `"Gamma"` are accepted as well. Under a family whose link is not the
  identity the arrows, the mean under `dsem_mu_spec` and the grid are
  all on the link scale, and every latent cell starts at the series
  mean. The sd and power parameters are mapped off for the families with
  none.

- dsem_link:

  Optional named character vector, one entry per covariate, the link
  from the cell to the observation's mean: `"identity"`, `"log"`,
  `"logit"` or `"cloglog"`. Defaults to each family's own: identity for
  fixed, normal and the fixed-sd normal, logit for bernoulli, log for
  the rest. A link applies whatever the family, so a Poisson under the
  identity can be handed a negative mean.

- dsem_fixed_sd:

  Data frame with a `year` column and one column per
  `"gaussian_fixed_sd"` covariate, holding that year's known sd. Needed
  on every observed year of such a covariate.

- dsem_mu_spec:

  `"est"` (default) estimates every covariate's mean, `"fix"` holds them
  all at the observed mean, a character vector estimates only those
  named, and a named numeric vector fixes those covariates at the values
  given with the rest at the observed mean. A linked series' mean is
  always fixed at zero, since its process already sits under a level
  parameter.

- covs:

  Passed to `read_dsem_arrows`: series groups whose innovations may
  covary.

- dsem_delta0_spec:

  `"none"` (default), `"est"`, or the names of the series whose first
  year offset is estimated.

- mod_var_logscale:

  Passed to `read_dsem_arrows`: whether a moderated sd is the
  exponential of its series.

- dsem_variance:

  Passed to `read_dsem_arrows` as `variance`. `"conditional"` (default)
  reads each sd line as the innovation sd; `"diagonal"` and `"marginal"`
  read it as the series' marginal sd, solving the innovation sd so the
  series comes out at that spread, and differ only when covariance lines
  are present. A linked recruitment cell's lognormal correction is half
  its marginal variance under the arrows, which `dsem_margvar_grid`
  reports and the marginal forms hold at the sd line squared; the
  initial age deviations read its settled value. A random walk cannot be
  written under them, since its paths alone carry a cell past any fixed
  spread after year one, and setup checks the starting values and
  refuses.

## Value

`input_list` with the dsem data, parameters (`dsem_beta`, `ln_dsem_sd`,
`dsem_mu`, `dsem_x`, `ln_dsem_obs_sd`, `logit_dsem_tweedie_p`,
`dsem_delta0`) and their map.

## Details

A linked recruitment cell is a random effect and takes the full
lognormal correction whenever its own penalty would take one: its mean
drops by half its variance under the arrows, that variance given the
covariate values the model is handed
([`get_dsem_margvar`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_margvar.md)),
so \\R_0\\ scales mean recruitment with or without the arrows. A ramp at
zero means none, and a nonzero ramp on a linked year is refused.
