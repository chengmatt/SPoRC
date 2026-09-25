# Processes a dsem can link to

Processes a dsem can link to

## Usage

``` r
dsem_process_table()
```

## Value

List, one entry per process: the `label` its series names start with,
the parameter array, the operating model array a drawn series is written
into (`sim_par`, the log state becomes the innovation `naa_eta_all`),
every dim of the parameter array in order (`"Yr"` marks the year dim),
the sigma its own penalty reads, and whether that penalty reads the
array's map mirror, which is what the dsem switches on.
