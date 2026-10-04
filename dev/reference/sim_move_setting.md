# One replicate's movement setting

The replicate's own under a joint self test, where `Setup_Sim_Movement`
stored one per replicate, else the one every replicate shares, with
`move_args` falling back on the dsem's.

## Usage

``` r
sim_move_setting(sim_env, name, sim)
```

## Arguments

- sim_env:

  Simulation environment.

- name:

  `"move_args"`, `"move_devs_fit"`, `"move_pe_pars"`,
  `"move_pop_corr_pars"`, `"move_seas_corr_pars"` or
  `"move_sex_corr_pars"`.

- sim:

  Replicate.

## Value

The setting.
