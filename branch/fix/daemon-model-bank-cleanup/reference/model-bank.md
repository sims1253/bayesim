# Model Bank for BrmsFitter

Internal infrastructure that compiles each distinct brms model spec once
(via `brms::brm(chains = 0)`) and reuses the resulting prefit across
simulation tasks via `stats::update(recompile = FALSE)`. This eliminates
catastrophic recompilation in studies that fit the same model to
thousands of datasets.

The bank is transported to mirai daemons via a session option
(`bayesim.model_bank`) rather than as an S7 fitter property, so it does
not corrupt the config fingerprint
([`capture_fitter_spec()`](https://sims1253.github.io/bayesim/reference/capture_fitter_spec.md)
hashes S7 properties only). See
[`set_model_bank()`](https://sims1253.github.io/bayesim/reference/set_model_bank.md)
/
[`get_model_bank()`](https://sims1253.github.io/bayesim/reference/get_model_bank.md).
