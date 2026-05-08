# Extract variable from an environment and remove that environment

This function will extract all variables from the given environment and
assign them to the global environment and then optionally remove the
environment.

## Usage

``` r
env2global(x, remove = TRUE)
```

## Arguments

- x:

  Environment to extract variables from

- remove:

  Should the Environment be removed after extracting variables? Default
  TRUE

## Value

variables from x are assigned to the GlobalEnv after execution

## Examples

``` r
my_env <- new.env()
my_env$X <- 1000
my_env$Y <- 1:10

# Extract X and Y to global environment and remove my_env
env2global(my_env)
#> Warning: object 'my_env' not found

# Extract X and Y to global environment and keep my_env
env2global(my_env, remove = FALSE)
```
