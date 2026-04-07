# version 1.2.0 (dev-popn-seasons)
- Follows the dev-movement branch with updates to incorporate population-specific and seasonal dynamics.

## Major changes
- Incorporated ability to simulate and estimate both population-specific (natal homing) and seasonal dynamics.
- Recoded tagging module to allow fleet-specific tag reporting rates, as well as missing attributes in tagged fish. 
- Reference points for local density-dependence in both meta-population and natal homing contexts.
- Added ability to simulate and fit to population-specific catches, indices, compositions, and tagging data. 

## Improvements
- Allowed both numeric and character codes for specifying dynamics in simulator, so as to be more consistent with how estimation models are specified, while maintaining backwards compatibility. 

# version 1.1.0
Release Date: 2026-3-31

## Major changes
- Added capability to estimate movement using continuous time Markov chains (CTMC) using preference functions. 

## Bug Fixes
- Fixed bug to allow for size-based selectivity in the estimation model. 

# version 1.0.0
Release Date: 2025-11-24

## Major changes
- First stable public release of SPoRC.
- Introduced a modular framework for estimating age-, sex-, season-, and region-structured population dynamics.
- Added closed-loop simulation and management strategy evaluation (MSE) capabilities.
