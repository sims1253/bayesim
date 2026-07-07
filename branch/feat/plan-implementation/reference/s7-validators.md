# S7 Property Validator Helpers

Factory functions that return validator closures suitable for use as the
`validator =` argument of
[`S7::new_property()`](https://rconsortium.github.io/S7/reference/new_property.html).
Each factory captures a small, commonly-repeated validation rule and
returns a `function(value)` that follows the S7 convention: it returns
`NULL` when `value` is valid, or a single character string describing
why `value` is invalid.

These helpers exist to deduplicate the small S7 property validators that
are otherwise copy-pasted across class definitions (B5).
