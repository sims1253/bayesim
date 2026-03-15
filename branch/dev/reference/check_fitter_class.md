# Validate Fitter Class Hierarchy

Validates that a fitter object is an S7 instance of the Fitter class.
This is a lightweight class-hierarchy check used internally by
[`validate_simulation_config()`](https://sims1253.github.io/bayesim/reference/validate_simulation_config.md).
For full interface validation including method existence and optional
smoke testing, use
[`validate_fitter()`](https://sims1253.github.io/bayesim/reference/validate_fitter.md).

## Usage

``` r
check_fitter_class(fitter)

validate_fitter_interface(fitter)
```

## Arguments

- fitter:

  An S7 object to validate as a Fitter.

## Value

The input `fitter`, invisibly, if validation passes.

## Details

The fitter must satisfy the following requirements:

- Must be an S7 object (checked via S7::S7_inherits())

- Must inherit from the "Fitter" class

Method existence is not checked here because S7 methods are dispatched
via generics, not stored as properties. The Fitter base class uses
S7::stop_method_not_implemented() for abstract methods, so subclasses
that don't override will raise errors when methods are called.

## Errors

Throws a `bayesim_contract_error` condition if validation fails.

## See also

[Fitter](https://sims1253.github.io/bayesim/reference/Fitter.md),
[`validate_fitter()`](https://sims1253.github.io/bayesim/reference/validate_fitter.md)

## Examples

``` r
# Validate the mock fitter
mock_fitter <- MockFitter()
check_fitter_class(mock_fitter)
```
