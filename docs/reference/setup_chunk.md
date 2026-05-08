# Setup Chunk Generator

This function prints R code to the console that can be copied and pasted
into an analysis script or quarto file. It sets up standard libraries
and directory paths. If installed, `clipr` is used to copy the code
chuck to the clipboard.

## Usage

``` r
setup_chunk(extra_libs = NULL, dir_name = "01", quiet = FALSE)
```

## Arguments

- extra_libs:

  Character vector. Additional library names to load.

- dir_name:

  Character string. The specific subdirectory for results (e.g., "01",
  "02").

- quiet:

  boolean. If TRUE do not print the code chunk to the console (only
  attempt to copy to clipboard).

## Value

Prints formatted code to console. Returns the character vector
invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
 setup_chunk()
} # }
```
