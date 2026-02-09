# SPoRC (dev-seasons)
- Follows the dev-movement branch with updates to incorporate seasonality.

## Major changes
- Incorporated ability to simulate seasonal dynamics.

## Improvements
- Allowed both numeric and character codes for specifying dynamics in simulator, so as to be more consistent with how estimation models are specified, while maintaining backwards compatibility. 


# SPoRC (dev-movement)

## Major changes
- Added capability to estimate movement using continuous time Markov chains (CTMC) using preference functions. 

## Bug Fixes
- Fixed bug to allow for size-based selectivity in the estimation model. 

# SPoRC 1.0.0
Release Date: 2025-11-24

## Major changes
- First stable public release of SPoRC.
- Introduced a modular framework for estimating age-, sex-, season-, and region-structured population dynamics.
- Added closed-loop simulation and management strategy evaluation (MSE) capabilities.
