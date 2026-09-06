# Extract a named parameter vector from a draws matrix at a given draw index. Errors if a requested variable cannot be found (neither as the cleaned name nor as "b\_"), since a silent NA would corrupt downstream SBC ranks.

Extract a named parameter vector from a draws matrix at a given draw
index. Errors if a requested variable cannot be found (neither as the
cleaned name nor as "b\_"), since a silent NA would corrupt downstream
SBC ranks.

## Usage

``` r
.extract_truth(draws_mat, draw_id, vars_of_interest)
```
