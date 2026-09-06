# Determine the response dependency sequence of a brms model.

S3 generic dispatched on the model object type. Returns a list of
character vectors: each element is the set of response names at a given
dependency depth, to be simulated together (topologically ordered so
that a response depending on another response is simulated after its
predecessor). This is the value iterated over by
[`brms_full_ppred()`](https://sims1253.github.io/bayesim/reference/brms_full_ppred.md)'s
`for (vars in resp)` loop, so each element must be a *response-name*
group, not a dependency list.

## Usage

``` r
brms_response_sequence(x)
```

## Details

Ported from the 0.x `brms_response_sequence.bform` method. Methods are
provided for both `brmsformula` (univariate) and `mvbrmsformula`
(multivariate); both feed into the adjacency-matrix +
[`nodes_by_depth()`](https://sims1253.github.io/bayesim/reference/nodes_by_depth.md)
construction (the 0.x `bform` logic), so dispatch works regardless of
whether the `"bform"` class alias is registered.
