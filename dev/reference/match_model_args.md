# Arguments of a model function read from the data and parameter lists

Matches each argument of `fn` by name: a value given in `...` first,
then the parameter list, then the data list. Arguments found nowhere
keep the function's own default.

## Usage

``` r
match_model_args(fn, data, pars, ...)
```

## Arguments

- fn:

  A model function such as `Get_Growth`.

- data:

  Data list of the fit.

- pars:

  Parameter list at the fitted values.

- ...:

  Values for arguments the lists lack or hold under another name.

## Value

Named list ready for `do.call(fn, ...)`.
