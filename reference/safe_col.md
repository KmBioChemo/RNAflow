# Safely coerce a value to a hex color, with fallback

Used everywhere we read a user input that should be a color. Guarantees
the return value is always a valid hex code, never NA or invalid.

## Usage

``` r
safe_col(val, fb = "#888888")
```

## Arguments

- val:

  candidate value

- fb:

  fallback color if val is invalid

## Value

a valid hex color string
