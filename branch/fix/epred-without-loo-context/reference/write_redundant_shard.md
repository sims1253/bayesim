# Write one immutable shard with an independently checksummed mirror.

The descriptor is committed last, so readers never observe a partially
published shard. Mirroring keeps later checkpoints readable after damage
to any single historical shard without rewriting accumulated history.

## Usage

``` r
write_redundant_shard(object, path)
```
