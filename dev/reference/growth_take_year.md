# Copy one year of the growth module's output into the model's own arrays

Under cohort growth a year's size-age transition, lengths, weights and
growth parameters are only known once the population loop reaches that
year, so they are written into the arrays the rest of the model reads
one year at a time.

## Usage

``` r
growth_take_year(dest, g, y, derive_waa)
```

## Arguments

- dest:

  Named list of the model's growth arrays, as `growth_take_year` returns
  them.

- g:

  The growth module's output, from
  [`Get_Growth_Year`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Growth_Year.md).

- y:

  Year index to copy.

- derive_waa:

  Integer (0/1); whether the weight arrays are present.

## Value

`dest` with year `y` filled.
