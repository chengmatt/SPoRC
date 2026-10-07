# Arguments to Get_Selex_Array for one fleet type

The same call the objective makes for fishery, retention or survey
selectivity, read off the fit's data and a parameter list.

## Usage

``` r
sim_selex_args(type, data, pars)
```

## Arguments

- type:

  `"fish"`, `"ret"` or `"srv"`.

- data:

  Data list of the fit.

- pars:

  Parameter list.

## Value

Named list of `Get_Selex_Array` arguments.
