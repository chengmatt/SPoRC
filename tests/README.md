# SPoRC tests

Read this before you change a test to make it pass.

Tests here fall into two kinds, and they call for opposite responses when they fail.
Most of the suite is the first kind, where a failure describes itself. A small set is
the second kind, where a failure is ambiguous unless you know where the numbers came
from. That set is listed below.

## Self-validating tests

These check a property rather than a stored number, so a failure tells you what broke.

- `test-sim_selftest_*.R` simulate data from known parameters, refit, and check recovery.
  A failure means the model no longer recovers the parameters it generated.
- `test-integration_*.R` compare two paths that must agree, such as an operating model
  against the estimation model. A failure means the two have diverged.
- `expect_jnLL_decomposes()` (see `helper-jnll_decomposition.R`) checks that the likelihood
  components sum to the reported total. A failure means a component is missing from the sum
  or double counted, which is the usual symptom of adding a likelihood without wiring it in.
- `test-setup_*.R` and `test-utils_*.R` check input validation, dimension handling, and
  mapping. A failure means a setup function accepts something it should reject, or rejects
  something it should accept.
- `test-model_population_dynamics.R` writes its expectations as arithmetic, for example
  `40 * exp(-0.3)`, so each assertion carries its own derivation.


## Pinned regression tests

These assert against stored numeric vectors:

- `test-regression_dusky.R`
- `test-regression_ebs_pollock_sgl.R`
- `test-regression_sabie_sgl.R`
- `test-regression_sabie_three_rg.R`
- `test-refpts_sgl_rg_spr.R`

The stored values are output from previous SPoRC fits of these assessments, taken from runs
we had checked and trusted at the time. They are not hand-derived, and they are not imported
from an external model. They exist to catch a change that silently moves a fitted result.

When one fails:

1. Assume you introduced a bug. This is by far the most common cause. Something in the
   likelihood, the population dynamics, or the setup path changed the fit without you
   intending it.
2. If the numerical change really was intended, re-baseline on purpose. Regenerate the
   vectors, confirm the new fit is one you would defend, and record the change in `NEWS.md`
   with the reason. A reviewer should be able to see from the diff that the move was
   deliberate.
3. Never regenerate because the diff looked small. A change of a few percent in SSB is
   exactly what a real bug looks like.

If a failure appears on one platform only and sits in the last digit or two, that is
optimizer or BLAS sensitivity rather than a code change. Give a pinned comparison enough
tolerance to absorb that, and do not tighten one to `tolerance = 0`.

## Adding a new model or assessment

Prefer a self-test or a parity check over a new pinned vector. Pin a fit only when the point
of the test is that this specific configuration keeps producing this specific answer. If you
do pin one, add a header comment saying where the values came from and add the file to the
list above.
