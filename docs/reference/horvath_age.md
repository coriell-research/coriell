# Transform ages using Horvath's method

Transform age in years using the log transformation defined in Horvath
(2013). This function can also be used to convert transformed ages back
to the original scale.

## Usage

``` r
horvath_age(x, adult_age = 20, inverse = FALSE)
```

## Arguments

- x:

  numeric vector of ages

- adult_age:

  Age of the adult of the species. From the original manuscript, 20 for
  humans and 15 for chimpanzees.

- inverse:

  If TRUE, return transformed ages to the original scale

## Value

transformed ages

## Details

[Horvath, S. DNA methylation age of human tissues and cell types. Genome
Biol 14, 3156 (2013).
https://doi.org/10.1186/gb-2013-14-10-r115](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115)

## Examples

``` r
age <- 0:100
transformed <- horvath_age(age)
plot(age, transformed, xlab="Age", ylab="Transformed Age", type="l")


# Conversion back to original scale
all.equal(age, horvath_age(transformed, inverse = TRUE))
#> [1] TRUE
```
