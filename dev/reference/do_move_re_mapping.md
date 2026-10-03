# Map movement deviations and their process error parameters

Internal helper called by
[`Setup_Mod_Movement`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md).
The deviations are a field over origin-destination pair, population,
year, season, age and sex, every pair with its own, and every other dim
is cut into blocks by its switch
([`move_dim_blocks`](https://chengmatt.github.io/SPoRC/dev/reference/move_dim_blocks.md)):
cells whose blocks agree on every dim share one deviation, and a cell
with an inactive level on any dim, or a pair with no edge, holds none.
The process error parameters are one log sd and two AR1 correlations per
process error block of pairs, from `move_pe_spec`, with a correlation
mapped off unless its dim is `"ar1"`, all three mapped off under a dsem,
and all three mapped off under `"fix"` (held at their starting value;
any deviations are then estimated against that fixed sd); the
unstructured correlations are estimated only under `"us"`. Nothing is
built when movement is fixed, the model has one region, or every switch
is `"none"`.

## Usage

``` r
do_move_re_mapping(
  input_list,
  move_year_re,
  move_age_re,
  move_pop_re,
  move_seas_re,
  move_sex_re,
  move_pe_spec
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`, with the `move_re_*`
  active sets already in `$data`.

- move_year_re, move_age_re, move_pop_re, move_seas_re, move_sex_re,
  move_pe_spec:

  As in
  [`Setup_Mod_Movement`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md).

## Value

The input `input_list` with the factor maps of `move_devs`,
`move_pe_pars` and the three correlation vectors, the integer mirror
`$data$map_move_devs`, the pair table `$data$move_pairs`, the block of
every population, year, season, age and sex in `$data$move_*_block` and
the process error block of every pair in `$data$move_pe_block`.
