# Fill in defaults for input lists built by older versions of SPoRC

Assigns any data or parameter objects that a current `SPoRC_rtmb` call
expects but that older input lists predate, so previously built objects
keep evaluating unchanged. Values are written into `env` only when
absent and are never overwritten.

## Usage

``` r
maintain_backwards_compatibility(env = parent.frame())
```

## Arguments

- env:

  Environment holding the unpacked data and parameters, i.e. the
  `SPoRC_rtmb` frame after
  [`RTMB::getAll`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).
  Defaults to the caller.

## Value

`NULL`, invisibly. Called for its side effects on `env`.
