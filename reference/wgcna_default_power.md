# WGCNA's recommended default soft-thresholding power

Used when scale-free topology fit does not reach the target (e.g. small
sample sizes). Values follow the WGCNA FAQ recommendations.

## Usage

``` r
wgcna_default_power(n_samples, network_type = "signed")
```

## Arguments

- n_samples:

  number of samples

- network_type:

  network type

## Value

a numeric power
