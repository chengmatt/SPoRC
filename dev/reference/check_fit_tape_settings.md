# Check the tape settings passed to fit_model

[`RTMB::TapeConfig`](https://rdrr.io/pkg/RTMB/man/Tape.html) stops on a
bad value only once the fit has begun, and
[`TMB::config`](https://rdrr.io/pkg/TMB/man/config.html) keeps a
misspelled name without a word, so both lists are checked before
anything is built and the error names the argument.

## Usage

``` r
check_fit_tape_settings(tape_config, tmb_config)
```

## Arguments

- tape_config:

  `NULL` or a named list of `comparison`, `atomic` and `vectorize`, as
  in
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).

- tmb_config:

  `NULL` or a named list of TMB settings, as in
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).

## Value

`NULL`, invisibly. Stops on any name or value not allowed.
