# Starting Values and Fixing (and Sharing) Parameters

In addition to the `SPoRC::Setup_xxx` functions, users can access
several advanced model features. These include:

- Specifying starting values (closely tied to parameter fixing)  
- Fixing parameters  
- Sharing parameters across model partitions

## Starting Values

Starting values can be specified in two ways:

1.  Directly within the `SPoRC::Setup_xxx` functions. Each setup
    function accepts starting values through the `...` argument. The
    inputs must match the model’s parameter names and dimensions. For
    more information on parameter dimensions, see the *Description of
    Model Parameters* vignette.

- [`SPoRC::Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Rec.md)  
  *ln_global_R0*, *Rec_prop*, *steepness_h*, *ln_InitDevs*,
  *ln_RecDevs*, *ln_sigmaR*

- [`SPoRC::Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Biologicals.md)  
  *ln_M*, *M_offset*

- [`SPoRC::Setup_Mod_Movement`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Movement.md)  
  *move_pars*, *move_devs*, *move_pe_pars*, *log_move_diffusion_pars*,
  *move_preference_pars*

- [`SPoRC::Setup_Mod_Tagging`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Tagging.md)  
  *ln_Init_Tag_Mort*, *ln_Tag_Shed*, *ln_tag_theta*,
  *Tag_Reporting_Pars*

- [`SPoRC::Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Catch_and_F.md)  
  *ln_sigmaC*, *ln_sigmaF*, *ln_sigmaF_agg*, *ln_F_mean*, *ln_F_devs*,
  *ln_F_mean_AggCatch*, *ln_F_devs_AggCatch*

- [`SPoRC::Setup_Mod_FishIdx_and_Comps`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_FishIdx_and_Comps.md)  
  *ln_FishAge_theta*, *FishAge_corr_pars*, *ln_FishAge_theta_agg*,
  *FishAge_corr_pars_agg*,  
  *ln_FishLen_theta*, *FishLen_corr_pars*, *ln_FishLen_theta_agg*,
  *FishLen_corr_pars_agg*

- [`SPoRC::Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_SrvIdx_and_Comps.md)  
  *ln_SrvAge_theta*, *SrvAge_corr_pars*, *ln_SrvAge_theta_agg*,
  *SrvAge_corr_pars_agg*,  
  *ln_SrvLen_theta*, *SrvLen_corr_pars*, *ln_SrvLen_theta_agg*,
  *SrvLen_corr_pars_agg*

- [`SPoRC::Setup_Mod_Fishsel_and_Q`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Fishsel_and_Q.md)  
  *ln_fish_fixed_sel_pars*, *ln_fish_q*, *fishsel_pe_pars*,
  *ln_fishsel_devs*

- [`SPoRC::Setup_Mod_Srvsel_and_Q`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Srvsel_and_Q.md)  
  *ln_srv_fixed_sel_pars*, *ln_srv_q*, *srvsel_pe_pars*,
  *ln_srvsel_devs*

2.  Post-hoc modification of starting values. Alternatively, you can
    first call the setup functions without specifying starting values,
    then access and modify the internally created parameter list
    (`input_list$par`) before running the model.

In the following, we will illustrate both methods using the recruitment
module (`Setup_Mod_Rec`) and specify starting values for `ln_global_R0`
and `ln_sigmaR`. First let us load in the package and define the model
dimensions.

``` r
# Load in packages
library(SPoRC) 
data("sgl_rg_sable_data") # load in data

input_list <- Setup_Mod_Dim(years = 1:length(sgl_rg_sable_data$years), # vector of years 
                            # (corresponds to year 1960 - 2024)
                            ages = 1:length(sgl_rg_sable_data$ages), # vector of ages
                            lens = seq(41,99,2), # number of lengths
                            n_regions = 1, # number of regions
                            n_sexes = sgl_rg_sable_data$n_sexes, # number of sexes == 1,
                            # female, == 2 male
                            n_fish_fleets = sgl_rg_sable_data$n_fish_fleets, # number of fishery
                            # fleet == 1, fixed gear, == 2 trawl gear
                            n_srv_fleets = sgl_rg_sable_data$n_srv_fleets, # number of survey fleets
                            verbose = FALSE
                            )
```

We can specify starting values directly using the `Setup_Mod_Rec`
function. Note that all inputs passed via the `...` argument must
exactly match the parameter names and their expected dimensions in the
model (see the *Description of Model Parameters* vignette for details).

``` r
input_list <- Setup_Mod_Rec(
  input_list = input_list,        # input data list from above
  # Model options
  do_rec_bias_ramp = 0,       # disable bias ramp
  sigmaR_switch = as.integer(length(1960:1975)),  # switch from early to late sigmaR
  dont_est_recdev_last = 1,       # do not estimate last recruitment deviate
  rec_model = "mean_rec",         # recruitment model type
  init_age_strc = 1,              # geometric series for initial age structure
  
  # Specify starting values
  ln_global_R0 = log(30),         # starting value for global R0
  ln_sigmaR = c(log(1.5), log(1.5))  # starting values for early and late sigmaR
)
```

In this example, the starting value for *ln_global_R0* is set to
`log(30)`, while *ln_sigmaR* is set to `log(1.5)` for both the early
(first element) and late (second element) periods.

Alternatively, starting values can be assigned after running the setup
functions without initial specifications. Users can extract the internal
parameter list and modify starting values as needed:

``` r
input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                            # Model options
                            do_rec_bias_ramp = 0, # don't do bias ramp 
                            sigmaR_switch = as.integer(length(1960:1975)), # when to switch from early to late sigmaR
                            dont_est_recdev_last = 1, # don't estimate last recruitment deviate
                            rec_model = "mean_rec", # recruitment model
                            init_age_strc = 1 # geometric series to derive age structure
                            )

# Specify starting values post-hoc
# R0
input_list$par$ln_global_R0 # default starting value
#> [1] 2.70805
input_list$par$ln_global_R0 <- log(30) # user specified starting value

# sigmaR
input_list$par$ln_sigmaR # default starting value
#> [1] 0 0
input_list$par$ln_sigmaR[] <- c(log(1.5), log(1.5)) # user specified starting value
```

## Mapping

Mapping is a core feature of TMB and RTMB models. It allows users to
either fix parameters at known values or to share parameters across
different parts of the model. In the sections below, we first
demonstrate how to use mapping to fix parameters. We then show how
mapping can be used to share parameters across model partitions.

### Fixing Parameters

The `SPoRC::Setup_xxx` functions include arguments that allow certain
parameters to be fixed, meaning they are not estimated during model
fitting. These arguments are not available for all parameters, but we
note that all parameters can be specified to be fixed. As an example, we
use `Setup_Mod_Rec` to show how the *ln_sigmaR* parameter can be fixed.
The argument `sigmaR_spec = "fix"`.

``` r
input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                            # Model options
                            do_rec_bias_ramp = 0, # don't do bias ramp 
                            sigmaR_switch = as.integer(length(1960:1975)), # when to switch from early to late sigmaR
                            dont_est_recdev_last = 1, # don't estimate last recruitment deviate
                            rec_model = "mean_rec", # recruitment model
                            init_age_strc = 1, # geometric series to derive age structure
                            
                            # Parameter Fixing
                            sigmaR_spec = 'fix'
                            )

input_list$map$ln_sigmaR # both values are fixed and not estimated (specified as factor(rep(NA, 2)))
#> [1] <NA> <NA>
#> Levels:
input_list$par$ln_sigmaR # ln_sigmaR is then fixed at the default starting value
#> [1] 0 0
```

To fix *ln_sigmaR* at a specific value, simply supply that value when
calling the function:

``` r
input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                            # Model options
                            do_rec_bias_ramp = , # don't do bias ramp 
                            sigmaR_switch = as.integer(length(1960:1975)), # when to switch from early to late sigmaR
                            dont_est_recdev_last = 1, # don't estimate last recruitment deviate
                            rec_model = "mean_rec", # recruitment model
                            init_age_strc = 1, # geometric series to derive age structure
                            
                            # Parameter Fixing
                            sigmaR_spec = 'fix',
                            ln_sigmaR = c(log(1.5), log(1.5)) # user specified starting value
                            )

input_list$map$ln_sigmaR # both values are fixed and not estimated (specified as factor(rep(NA, 2)))
#> [1] <NA> <NA>
#> Levels:
input_list$par$ln_sigmaR # ln_sigmaR is then fixed at the user specified starting value
#> [1] 0.4054651 0.4054651
```

However, not all parameters include a convenience argument like
`sigmaR_spec = "fix"`. For example, *ln_global_R0* does not. In such
cases, you can fix the parameter manually by modifying the map list
directly and specifying the desired starting value:

``` r
input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                            # Model options
                            do_rec_bias_ramp = 0, # don't do bias ramp 
                            sigmaR_switch = as.integer(length(1960:1975)), # when to switch from early to late sigmaR
                            dont_est_recdev_last = 1, # don't estimate last recruitment deviate
                            rec_model = "mean_rec", # recruitment model
                            init_age_strc = 1, # geometric series to derive age structure
                            )

input_list$map$ln_global_R0 <- factor(NA)
input_list$par$ln_global_R0 <- log(30)
```

In this example, *ln_global_R0* is fixed at `log(30)` by setting its map
entry to `NA` and providing the desired value in the parameter list.

### Sharing Parameters

The `SPoRC::Setup_xxx` functions also support sharing parameters across
model partitions. For example, selectivity parameters can be shared
between sexes using built-in convenience arguments. These options
simplify common sharing structures, but users can also implement more
customized behavior when needed. Below, we demonstrate two approaches.
First, we show how to share fishery selectivity parameters between sexes
using a convenience flag in `Setup_Mod_Fishsel_and_Q`. Then, we outline
how users can manually configure parameter sharing for finer control.

In the example below, we define two fishery fleets. Fleet 1 uses a
logistic selectivity model (2 parameters), while Fleet 2 uses a gamma
dome-shaped model (2 parameters). We first use the
`fish_fixed_sel_pars_spec` argument to specify that selectivity
parameters for both fleets should be estimated across all model
partitions (`est_all`). Thus, when we inspect the map, we should expect
a total of 8 unique numbers, which represent unique selectivity
parameters to be estimated for each sex and fleet (4 for each sex and
fleet combination).

``` r
input_list$data$Selex_Type <- 0 # specifying age-based selectivity for demonstration purposes
input_list <- SPoRC::Setup_Mod_Fishsel_and_Q(input_list = input_list,

                                      # Model options
                                      cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
                                      fish_sel_blocks = c("none_Fleet_1", "none_Fleet_2"),
                                      fish_sel_model = c("logist1_Fleet_1", "gamma_Fleet_2"),
                                      fish_q_blocks = c("none_Fleet_1", "none_Fleet_2"),
                                      fish_q_spec = c("fix", "fix"),

                                      # Share selectivity parameters across all partitions
                                      fish_fixed_sel_pars_spec = c("est_all", "est_all"))

input_list$map$ln_fish_fixed_sel_pars # 8 unique numbers, 4 for each sex and fleet combination
#> [1] <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA>
#> Levels:
```

Next, we use the `fish_fixed_sel_pars_spec` argument to indicate that
selectivity parameters should be shared across sex partitions for both
fleets. Specifically, we use the `"est_shared_s"` setting to link
parameters between sexes, even though the model is sex-structured. In
this setup, we still estimate separate parameters for each fleet, but
the same selectivity parameters are used for both sexes within each
fleet. As a result, when we inspect the mapping, we should expect four
unique values—two per fleet—despite having two sexes.

``` r
input_list$data$Selex_Type <- 0 # specifying age-based selectivity for demonstration purposes
input_list <- SPoRC::Setup_Mod_Fishsel_and_Q(input_list = input_list,

                                      # Model options
                                      cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
                                      fish_sel_blocks = c("none_Fleet_1", "none_Fleet_2"),
                                      fish_sel_model = c("logist1_Fleet_1", "gamma_Fleet_2"),
                                      fish_q_blocks = c("none_Fleet_1", "none_Fleet_2"),
                                      fish_q_spec = c("fix", "fix"),

                                      # Share selectivity parameters across all partitions
                                      fish_fixed_sel_pars_spec = c("est_shared_s", "est_shared_s"))

input_list$map$ln_fish_fixed_sel_pars # 4 unique numbers, 2 for each sex and fleet combination
#> [1] <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA>
#> Levels:
```

The parameter array ln_fish_fixed_sel_pars has the following structure:
`[n_regions, n_max_sel_pars, n_max_sel_blocks, n_sexes, n_fish_fleets]`
In this example, its dimensions are `1 x 2 x 1 x 2 x 2`.

This means that:

- For each fleet, there are two estimated selectivity parameters

- These are shared across the two sexes

- When flattened into a vector (as required by the mapping), this
  results in four distinct values (two per fleet), even though the
  underlying array spans two sexes

Below is a minimal example that replicates what
`Setup_Mod_Fishsel_and_Q` is doing internally:

``` r
# define dimensions
n_regions <- 1
n_max_sel_pars <- 2
n_max_sel_blocks <- 1
n_sexes <- 2
n_fish_fleets <- 2

# define empty parameter array for demonstration
ln_fish_fixed_sel_pars <- array(0, dim = c(n_regions,n_max_sel_pars,n_max_sel_blocks,n_sexes,n_fish_fleets))
ln_fish_fixed_sel_pars[1,,1,,1] <- c(1,2) # same parameters for each sex for fleet 1
ln_fish_fixed_sel_pars[1,,1,,2] <- c(3,4) # same parameters for each sex for fleet 2
ln_fish_fixed_sel_pars
#> , , 1, 1, 1
#> 
#>      [,1] [,2]
#> [1,]    1    2
#> 
#> , , 1, 2, 1
#> 
#>      [,1] [,2]
#> [1,]    1    2
#> 
#> , , 1, 1, 2
#> 
#>      [,1] [,2]
#> [1,]    3    4
#> 
#> , , 1, 2, 2
#> 
#>      [,1] [,2]
#> [1,]    3    4

# flatten array to a vector of factors for the map
custom_map <- factor(as.vector(ln_fish_fixed_sel_pars))
```

We can now compare the custom map object we constructed manually with
the one generated by `Setup_Mod_Fishsel_and_Q`:

``` r
custom_map
#> [1] 1 2 1 2 3 4 3 4
#> Levels: 1 2 3 4
input_list$map$ln_fish_fixed_sel_pars
#> [1] <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA>
#> Levels:
```

Both should match, confirming that parameters are being correctly shared
across sexes for each fleet.

Lastly, we demonstrate how to manually specify advanced mapping and
parameter-sharing options for fishery selectivity parameters. In this
example, we want to:

- Estimate logistic selectivity for Fleet 1, where:
  - The $a_{50}$ parameter (inflection point) is sex-specific  
  - The slope parameter is shared across sexes  
- Estimate gamma selectivity for Fleet 2, where:
  - All parameters are sex-specific

This level of control is not supported by the built-in
`fish_fixed_sel_pars_spec` convenience arguments, so we must construct
the mapping manually.

We begin by calling `Setup_Mod_Fishsel_and_Q` with the `"est_all"`
option for both fleets. This creates a parameter array with the correct
dimensions and assumes that all selectivity parameters are uniquely
estimated across model partitions.

``` r
input_list$data$Selex_Type <- 0 # specifying age-based selectivity for demonstration purposes
input_list <- SPoRC::Setup_Mod_Fishsel_and_Q(
  input_list = input_list,

  # Model options
  cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
  fish_sel_blocks = c("none_Fleet_1", "none_Fleet_2"),
  fish_sel_model = c("logist1_Fleet_1", "gamma_Fleet_2"),
  fish_q_blocks = c("none_Fleet_1", "none_Fleet_2"),
  fish_q_spec = c("fix", "fix"),

  # Start with all parameters estimated independently
  fish_fixed_sel_pars_spec = c("est_all", "est_all")
)
```

We then extract the parameter array created internally for
`ln_fish_fixed_sel_pars.` This array has dimensions:
`[n_regions, n_max_sel_pars, n_max_sel_blocks, n_sexes, n_fish_fleets]`.
In this example, its dimensions are `1 x 2 x 1 x 2 x 2`. We can then
assign integer values to this array to indicate which parameters are
shared (same value) or independently estimated (different values).

``` r
map_ln_fish_fixed_sel_pars <- input_list$par$ln_fish_fixed_sel_pars  # extract the array
```

We now assign specific values to define the desired sharing structure:

``` r
# Fleet 1 (logistic): 
# - sel_par 1 = a50 (sex-specific)
# - sel_par 2 = slope (shared across sex)
map_ln_fish_fixed_sel_pars[1,1,1,,1] <- c(1,2)   # a50: unique for each sex
map_ln_fish_fixed_sel_pars[1,2,1,,1] <- c(3)     # slope: same for both sexes

# Fleet 2 (gamma): 
# - sel_par 1 = amax (sex-specific)
# - sel_par 2 = slope (sex-specific)
map_ln_fish_fixed_sel_pars[1,1,1,,2] <- c(4,5)   # amax: unique per sex
map_ln_fish_fixed_sel_pars[1,2,1,,2] <- c(6,7)   # slope: unique per sex
```

We then flatten this array into a vector of factors and assign it to the
map list. This tells the model how parameters should be linked when
building the model.

``` r
# Flatten the mapping array and assign it to the map
input_list$map$ln_fish_fixed_sel_pars <- factor(as.vector(map_ln_fish_fixed_sel_pars))

# View the final map
input_list$map$ln_fish_fixed_sel_pars
#> [1] 1 3 2 3 4 6 5 7
#> Levels: 1 2 3 4 5 6 7
```

This custom mapping approach allows full control over parameter sharing
structures beyond what is available through the high-level setup
arguments.
