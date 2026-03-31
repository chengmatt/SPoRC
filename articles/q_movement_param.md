# Alternative Approaches for Movement Parameterization

The `SPoRC` package offers multiple approaches for modeling fish
movement between regions, each with different complexity and flexibility
trade-offs. This vignette demonstrates how to configure movement
parameterization using the `Setup_Mod_Movement` function with various
structural options. In general, there are two primary movement model
types, which include:

1.  **Unstructured Markov**: Discrete movement matrices with optional
    blocking structures across age, year, and sex dimensions
2.  **Continuous-Time Markov Chain**: CTMC-based movement with diffusion
    and preference parameters estimated using formula-based approaches

The following demonstrations use the three region sablefish dataset as a
basis (`three_rg_sable_data`). The initial setup establishes the genearl
model dimensions.

``` r
# Load in packages
library(SPoRC) 
data("three_rg_sable_data") # load in data

# setup model dimensions
input_list <- Setup_Mod_Dim(years = 1:length(three_rg_sable_data$years), # vector of years 
                            # (corresponds to year 1960 - 2024)
                            ages = 1:length(three_rg_sable_data$ages), # vector of ages
                            lens = seq(41,99,2), # number of lengths
                            n_regions = three_rg_sable_data$n_regions, # number of regions
                            n_sexes = three_rg_sable_data$n_sexes, # number of sexes == 1,
                            # female, == 2 male
                            n_fish_fleets = three_rg_sable_data$n_fish_fleets, # number of fishery
                            # fleet == 1, fixed gear, == 2 trawl gear
                            n_srv_fleets = three_rg_sable_data$n_srv_fleets, # number of survey fleets
                            verbose = TRUE
                            )
```

Like previous vignettes, the
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Dim.md)
function initializes the model structure with 62 years of data
(1960-2021), 30 age classes, 3 regions, 2 sexes, 2 fishery fleets (fixed
gear and trawl), and multiple survey fleets.

## Unstructured Markov

In the ‘simplest’ case, movement can be parameterized as an unstructured
Markov model (`move_type = 0`), where movement parameters are constant
across model partitions and are estimated as discrete transitions
between regions, with `nregions x nregions - 1` parameters estimated.

``` r
# Setup movement
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  do_recruits_move = 0, # recruits don't move
  move_type = 0 # unstructured markov
)
#> Movement is: Estimated
#> Movement priors are: Not Used
#> Recruits are: Not Moving
#> Continuous movement specification is: none
#> Continuous movement process error specification is: none
#> Movement type is: Unstructured Markov
#> Movement fixed effect blocks are sex-invariant
#> Movement fixed effect blocks are time-invariant
#> Movement fixed effect blocks are age-invariant
length(unique(input_list$map$move_pars)) # number of parameters estimated
#> [1] 6
```

### Age Blocks

The unstructured Markov model can be extended to incorporate blocking
structures, with parameter sharing within blocks to reduce parameter
load. In the following, we specify an unstructured Markov model with 2
age blocks, but constant movement across years and sexes. This results
in `nregions x nregions - 1 x 2` parameters estimated.

``` r
# define age block
age_blk <- list(c(1:15), c(16:30))
age_blk
#> [[1]]
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
#> 
#> [[2]]
#>  [1] 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30

# Setup movement
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  do_recruits_move = 0,  # recruits don't move
  move_type = 0, # unstructured markov
  Movement_ageblk_spec = age_blk, # age blocks
  Movement_yearblk_spec = "constant", # constant movement across years
  Movement_sexblk_spec = "constant" # constant movement across sexes
)
#> Movement is: Estimated
#> Movement priors are: Not Used
#> Recruits are: Not Moving
#> Continuous movement specification is: none
#> Continuous movement process error specification is: none
#> Movement type is: Unstructured Markov
#> Movement fixed effect blocks are sex-invariant
#> Movement fixed effect blocks are time-invariant
#> Movement fixed effect blocks are specified with 2 age blocks

length(unique(input_list$map$move_pars)) # number of parameters estimated
#> [1] 12
```

### Year Blocks

Year blocks can be specified in a similar fashion. Here, 5 year blocks
are specified, resulting in `nregions x nregions - 1 x 5` parameters
estimated.

``` r
# define year blocks
yr_blk <- list(c(1:15), c(16:30), c(31:45), c(46:60), c(61:62))
yr_blk
#> [[1]]
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
#> 
#> [[2]]
#>  [1] 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
#> 
#> [[3]]
#>  [1] 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45
#> 
#> [[4]]
#>  [1] 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60
#> 
#> [[5]]
#> [1] 61 62

# Setup movement
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  do_recruits_move = 0, # recruits dont move
  move_type = 0, # unstructured markov
  Movement_ageblk_spec = "constant", # constant movement across ages
  Movement_yearblk_spec = yr_blk, # time blocks for movement
  Movement_sexblk_spec = "constant" # constant movement across sexes
)
#> Movement is: Estimated
#> Movement priors are: Not Used
#> Recruits are: Not Moving
#> Continuous movement specification is: none
#> Continuous movement process error specification is: none
#> Movement type is: Unstructured Markov
#> Movement fixed effect blocks are sex-invariant
#> Movement fixed effect blocks are specified with 5 year blocks
#> Movement fixed effect blocks are age-invariant

length(unique(input_list$map$move_pars)) # number of parameters estimated
#> [1] 30
```

### Sex Blocks

Sex blocks are specified similarly. The following example shows
sex-specific movement, resulting in `nregions x nregions - 1 x nsexes`
parameters estimated.

``` r
# define sex blocks
sx_blk <- as.list(1:2)
sx_blk
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] 2

# Setup movement
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  do_recruits_move = 0, # recr4uits don't move
  move_type = 0, # unstructured markov
  Movement_ageblk_spec = "constant",  # constant movement across ages
  Movement_yearblk_spec = 'constant', # constant movement across years
  Movement_sexblk_spec = sx_blk # sex-specific movement
)
#> Movement is: Estimated
#> Movement priors are: Not Used
#> Recruits are: Not Moving
#> Continuous movement specification is: none
#> Continuous movement process error specification is: none
#> Movement type is: Unstructured Markov
#> Movement fixed effect blocks are specified with 2 sex blocks
#> Movement fixed effect blocks are time-invariant
#> Movement fixed effect blocks are age-invariant
```

### Blocks Across all Dimensions

Building on the principles described above, movement blocks can also be
defined simultaneously across ages, years, and sexes. The example below
specifies two age blocks, five year blocks, and sex-specific movement.
This configuration results in
`nregions × (nregions - 1) × 2 × 5 × nsexes` movement parameters being
estimated. However, this parameterization is likely excessive and may
lead to an unstable model solution.

``` r
# define age block
age_blk <- list(c(1:15), c(16:30))
age_blk
#> [[1]]
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
#> 
#> [[2]]
#>  [1] 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30

# define year blocks
yr_blk <- list(c(1:15), c(16:30), c(31:45), c(46:60), c(61:62))
yr_blk
#> [[1]]
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
#> 
#> [[2]]
#>  [1] 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
#> 
#> [[3]]
#>  [1] 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45
#> 
#> [[4]]
#>  [1] 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60
#> 
#> [[5]]
#> [1] 61 62

# define sex blocks
sx_blk <- as.list(1:2)
sx_blk
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] 2

# Setup movement
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  do_recruits_move = 0, # recruits dont move 
  move_type = 0, # unstructured markov
  
  # blocks across model partitions
  Movement_ageblk_spec = age_blk, 
  Movement_yearblk_spec = yr_blk,
  Movement_sexblk_spec = sx_blk
)
#> Movement is: Estimated
#> Movement priors are: Not Used
#> Recruits are: Not Moving
#> Continuous movement specification is: none
#> Continuous movement process error specification is: none
#> Movement type is: Unstructured Markov
#> Movement fixed effect blocks are specified with 2 sex blocks
#> Movement fixed effect blocks are specified with 5 year blocks
#> Movement fixed effect blocks are specified with 2 age blocks
```

## Continuous-Time Markov Chain (CTMC)

One potential approach to reduce model parameterization is to represent
movement as a continuous-time Markov chain (CTMC) process, in which
transitions are governed by a mechanistic framework composed of
diffusion (random dispersal) and taxis (directed preference) components.
Unlike the unstructured Markov model, CTMC movement processes are
defined in continuous time and are converted to annual movement
fractions using the matrix exponential, thereby allowing for sequential
transitions among regions. To implement a CTMC-based movement model, an
adjacency matrix must first be defined to specify which regions are
connected. In the example below, all regions are assumed to be
connected, permitting individuals to move among any spatial strata
within a given period. An accompanying data frame is then created to
define the `regions`, `years`, `ages`, `sexes`, and any additional
covariates associated with movement preference.

``` r
adjacency <- igraph::as_adjacency_matrix(
  igraph::make_graph(
    ~ 1 - 2,
    2 - 3,
    1 - 3
  )
)

# make ctmc data
ctmc_data <- expand.grid(
  regions = 1:three_rg_sable_data$n_regions,
  years = 1:length(three_rg_sable_data$years),
  ages = 1:length(three_rg_sable_data$ages),
  sexes = 1:three_rg_sable_data$n_sexes
)
```

### Constant Movement

In the code chunk below, movement is defined as arising from a purely
diffusive process, which effectively represents constant movement across
regions. This formulation results in the estimation of `nregions`
parameters, providing a more parsimonious alternative to the
unstructured Markov approach that estimates `nregions × (nregions - 1)`
parameters. Movement among adjacent areas is determined by the defined
adjacency matrix, and no directional preference is specified in this
case.

The argument `area_r = rep(1, 3)` specifies how diffusive processes
scale with area size. Here, all areas are assumed to be equal. However,
when areas differ in size, `area_r` should be defined as proportional to
area, such that smaller areas are associated with higher diffusion
rates.

``` r
# setup formulas for CTMC
diffusion_formula = ~0 + factor(regions) # constant diffusion
preference_formula = ~0

# Setup movement
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  do_recruits_move = 0, # recruits dont move
  move_type = 1, # ctmc movement
  ctmc_move_dat = ctmc_data, # ctmc data
  adjacency_mat = adjacency, # adjacency matrix
  area_r = rep(1, 3), # equal areas

  # formulas for CTMC
  diffusion_formula = diffusion_formula,
  preference_formula = preference_formula
)
#> Movement is: Estimated
#> Movement priors are: Not Used
#> Recruits are: Not Moving
#> Continuous movement specification is: none
#> Continuous movement process error specification is: none
#> Movement type is: Continuous Time Markov Chain

length(input_list$par$log_move_diffusion_pars)
#> [1] 3
length(input_list$par$move_preference_pars)
#> [1] 0
```

### Age-Varying

Movement can also be specified to vary across model partitions. The
examples below illustrate how age-varying movement can be represented
within the CTMC framework using both linear and spline-based
relationships.

#### Linear

In this example, diffusion is assumed constant across ages, while
preference varies linearly by age within each region. This specification
results in the estimation of a single diffusion parameter, along with
`nregions` additional parameters describing age-specific movement
preferences.

``` r
# setup formulas for CTMC
diffusion_formula = ~1 # constant diffusion
preference_formula = ~0 + factor(regions):ages # linear age effects

# Setup movement
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  do_recruits_move = 0,
  move_type = 1,
  ctmc_move_dat = ctmc_data,
  adjacency_mat = adjacency,
  area_r = rep(1, 3), # equal areas

  # formulas for CTMC
  diffusion_formula = diffusion_formula,
  preference_formula = preference_formula
)
#> Movement is: Estimated
#> Movement priors are: Not Used
#> Recruits are: Not Moving
#> Continuous movement specification is: none
#> Continuous movement process error specification is: none
#> Movement type is: Continuous Time Markov Chain

length(input_list$par$log_move_diffusion_pars)
#> [1] 1
length(input_list$par$move_preference_pars)
#> [1] 3
```

#### Spline

In some cases, age-specific movement patterns may be more complex than a
simple linear relationship can represent. To capture nonlinear
variation, spline-based age-specific movement can be modeled using a
spline basis, providing smoother and more flexible relationships between
age and movement preference. In the example below, one diffusion
parameter is estimated, along with `nregions × 4` parameters describing
age-specific movement preferences.

``` r
# setup formulas for CTMC
diffusion_formula = ~1 # constant diffusion
preference_formula = ~0 + factor(regions):splines2::bSpline(ages, df = 4, intercept = TRUE) # spline based age effect

# Setup movement
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  do_recruits_move = 0,
  move_type = 1,
  ctmc_move_dat = ctmc_data,
  adjacency_mat = adjacency,
  area_r = rep(1, 3), # equal areas

  # formulas for CTMC
  diffusion_formula = diffusion_formula,
  preference_formula = preference_formula
)
#> Movement is: Estimated
#> Movement priors are: Not Used
#> Recruits are: Not Moving
#> Continuous movement specification is: none
#> Continuous movement process error specification is: none
#> Movement type is: Continuous Time Markov Chain

length(input_list$par$log_move_diffusion_pars)
#> [1] 1
length(input_list$par$move_preference_pars)
#> [1] 12
```

#### Movement Across all Dimensions

Movement can also be specified to vary continuously across all model
dimensions. In this example, movement is modeled as a function of
region, age, year, and sex. This specification results in the estimation
of `nregions × 4 × 6 × nsexes` parameters.

``` r
# setup formulas for CTMC
diffusion_formula = ~1 # constant diffusion
preference_formula = ~0 + factor(regions):
  splines2::bSpline(ages, df = 4, intercept = TRUE):
  splines2::bSpline(years, df = 6, intercept = TRUE):
  factor(sexes)

# Setup movement
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  do_recruits_move = 0,
  move_type = 1,
  ctmc_move_dat = ctmc_data,
  adjacency_mat = adjacency,
  area_r = rep(1, 3), # equal areas

  # formulas for CTMC
  diffusion_formula = diffusion_formula,
  preference_formula = preference_formula
)
#> Movement is: Estimated
#> Movement priors are: Not Used
#> Recruits are: Not Moving
#> Continuous movement specification is: none
#> Continuous movement process error specification is: none
#> Movement type is: Continuous Time Markov Chain

length(input_list$par$log_move_diffusion_pars)
#> [1] 1
length(input_list$par$move_preference_pars)
#> [1] 144
```

## Process Error

Lastly, process error deviations can also be incorporated into movement
estimates. In the example below, movement is modeled using a CTMC
framework, although process error can similarly be applied to an
unstructured Markov model. Movement is specified to vary smoothly across
ages using a spline function, while allowing independent and identically
distributed (`iid`) deviations across years for each source region
(`cont_vary_movement = 'iid_y'`). Process error variance parameters are
specified to be shared across source regions, ages, and sexes
(`Movement_cont_pe_pars_spec = 'est_shared_r_a_s'`). Additional options
for `cont_vary_movement` and `Movement_cont_pe_pars_spec` is described
in the function documentation
([`?Setup_Mod_Movement`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Movement.md)).
Users may alternatively treat these variance parameters as random
effects (integrated out via the Laplace approximation) or as penalized
likelihood terms, depending on how the argument
`Movement_cont_pe_pars_spec` is defined. When
`Movement_cont_pe_pars_spec = 'fix'`, users can supply a fixed variance
value directly through `input_list$par$move_pe_pars`. Note that
deviations are only estimated for destination regions (i.e.,
`nregions - 1`), as no deviation term is defined for the source region.

``` r
# setup formulas for CTMC
diffusion_formula = ~1 # constant diffusion
preference_formula = ~0 + factor(regions):splines2::bSpline(ages, df = 4, intercept = TRUE) # spline based age effect

# Setup movement
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  do_recruits_move = 0,
  move_type = 1,
  ctmc_move_dat = ctmc_data,
  adjacency_mat = adjacency,
  area_r = rep(1, 3), # equal areas

  # formulas for CTMC
  diffusion_formula = diffusion_formula,
  preference_formula = preference_formula,
  cont_vary_movement = 'iid_y', 
  Movement_cont_pe_pars_spec = 'est_shared_r_a_s'
)
#> Movement is: Estimated
#> Movement priors are: Not Used
#> Recruits are: Not Moving
#> Continuous movement specification is: iid_y
#> Continuous movement process error specification is: est_shared_r_a_s
#> Movement type is: Continuous Time Markov Chain

length(input_list$par$log_move_diffusion_pars)
#> [1] 1
length(input_list$par$move_preference_pars)
#> [1] 12
length(unique(input_list$map$move_devs))
#> [1] 372
```
