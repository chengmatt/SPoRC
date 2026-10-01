# Record which SPoRC build produced a fit

Reads the package version and, when the install came from GitHub through
remotes or pak, the commit it was built from. A package loaded with
[`devtools::load_all`](https://devtools.r-lib.org/reference/load_all.html)
from a git checkout reports that checkout's HEAD instead, and a plain
local install has neither, so the commit is `NA`. A saved fit then says
which code made it, so a later package version can tell which report
names to expect.

## Usage

``` r
fit_provenance()
```

## Value

Named list: `package_version`, `commit_sha` (40 character hash or `NA`),
`commit_source` (`"github"`, `"local_git"`, or `"unknown"`),
`r_version`, `rtmb_version`, and `fit_time`.
