# Set up catchability deviations for the operating model

Gives each fleet a process error on annual catchability, matching
`srv_q_model` and `fish_q_model` in the estimation setup. Catchability
in replicate \\i\\ becomes \\q\_{r,y,f} \exp(\epsilon\_{r,y,f,i})\\,
where \\q\\ is the mean `Setup_Sim_Fishery` and `Setup_Sim_Survey`
supplied.

## Usage

``` r
Setup_Sim_q_devs(
  sim_list,
  fish_q_model = NULL,
  srv_q_model = NULL,
  sigma_fish_q = 0,
  sigma_srv_q = 0,
  fish_q_rho = 0,
  srv_q_rho = 0
)
```

## Arguments

- sim_list:

  List of operating model inputs, with `fish_q` and `srv_q` already set.

- fish_q_model, srv_q_model:

  Character vectors, one per fleet, of `"none"` (default), `"iid"`,
  `"rw"`, `"ar1"` or `"dsem"`. A `"dsem"` fleet draws nothing here and
  takes its series from
  [`Setup_Sim_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_DSEM.md).

- sigma_fish_q, sigma_srv_q:

  Numeric arrays `[n_regions, n_fleets]` of deviation standard
  deviations on the log scale, or a single value used for every region
  and fleet. Default zero.

- fish_q_rho, srv_q_rho:

  Numeric arrays `[n_regions, n_fleets]` of ar1 correlations on the
  natural scale, read only under `"ar1"`. Default zero.

## Value

`sim_list` with the deviation settings stored.

## Details

Conditioning years read their deviations from the fit and draw nothing,
so a random walk or ar1 in the projection continues from the fit's last
deviation rather than restarting at zero.
