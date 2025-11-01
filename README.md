
<!-- README.md is generated from README.Rmd. Please edit that file -->

# structures

<!-- badges: start -->

<!-- badges: end -->

An R package that contains a class used to store and represents
molecular structures (Molecule3D). Includes parser for the mol2 file
format.

## Installation

You can install the development version of structures like so:

``` r
if (!require("remotes"))
    install.packages("remotes")

remotes::install_github("selkamand/structures")
```

## Quick Start

``` r
library(structures)

# Read a mol2 file into a Molecule3D object
path <- system.file(package="structures", "benzene.mol2")
molecule <- read_mol2(path)

# Set anchor to the first element, and translate atom so the anchor atom is at the origin (0, 0, 0)
molecule |>
  set_anchor_by_atom(eleno = 1) |>
  translate_molecule_to_origin()
#> ===================
#> Chemical Molecule3D
#> ===================
#> Name: benzene
#> Atoms: 12
#> Bonds: 12
#> 
#> See @atoms paramater for atom positions and @bonds paramater for bonds
```
