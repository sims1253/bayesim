# Default vars_of_interest for a brmsfit: population-level effects, plus the residual scale `sigma` when present (E5: previously defaulted to population effects only, silently excluding sigma/auxiliary parameters from SBC). brms names effects "b\_"; strip the "b\_" prefix for true_params names.

Default vars_of_interest for a brmsfit: population-level effects, plus
the residual scale `sigma` when present (E5: previously defaulted to
population effects only, silently excluding sigma/auxiliary parameters
from SBC). brms names effects "b\_"; strip the "b\_" prefix for
true_params names.

## Usage

``` r
.default_prior_vars(fit)
```
