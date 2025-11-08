
<!-- README.md is generated from README.Rmd. Please edit that file -->

# structures

<!-- badges: start -->

<!-- badges: end -->

An R package for working with 3D molecular structures in R.  
It provides the **`Molecule3D`** class for representing molecular
geometry, together with parsers for the `.mol2` file format and a
growing set of geometric manipulation tools.

------------------------------------------------------------------------

## Installation

You can install the development version of **structures** like so:

``` r
if (!require("remotes"))
  install.packages("remotes")

remotes::install_github("selkamand/structures")
```

------------------------------------------------------------------------

## Quick Start

``` r
library(structures)

# Read a mol2 file into a Molecule3D object
path <- system.file(package="structures", "benzene.mol2")
molecule <- read_mol2(path)

print(molecule)
#> ===================
#> Chemical Molecule3D
#> ===================
#> Name: benzene
#> Atoms: 12
#> Bonds: 12
#> Symmetry Axes: 0
#> 
#> -------------------
#> See @atoms paramater for atom positions
#> See @bonds paramater for bond positions
#> See @symmetry_axes for symmetry axes

# Set anchor to the first element and translate the molecule so that
# the anchor atom is at the origin (0, 0, 0)
molecule_centered <- molecule |>
  set_anchor_by_atom(eleno = 1) |>
  translate_molecule_to_origin()
```

------------------------------------------------------------------------

## Symmetry Axes

The `structures` package now supports defining and storing **rotational
symmetry axes** within a `Molecule3D` object.

### The `SymAxis` class

A **`SymAxis`** object represents a proper rotation axis (Cₙ) defined
by: - `Cn`: the fold/order of the axis (e.g., 2, 3, 4, …), - `posA` and
`posB`: two points in 3D space (`c(x, y, z)`) defining the axis line.

``` r
# Create a simple C3 symmetry axis through the Z axis
axis_C3 <- SymAxis(Cn = 3L, posA = c(0, 0, 0), posB = c(0, 0, 1))
print(axis_C3)
#> ===================
#> Symmetry Axis
#> ===================
#> Fold Symmetry (Cn): C3
#> PosA: 0, 0, 0
#> PosB: 0, 0, 1
```

### Adding symmetry axes to molecules

Symmetry axes can be attached to any `Molecule3D` object using
\[`add_symmetry_axis()`\]:

``` r
# Add a symmetry axis to the benzene molecule
molecule <- add_symmetry_axis(
  molecule,
  SymAxis(Cn = 6L, posA = c(0, 0, -1), posB = c(0, 0, 1))
)

# Inspect updated molecule summary
print(molecule)
#> ===================
#> Chemical Molecule3D
#> ===================
#> Name: benzene
#> Atoms: 12
#> Bonds: 12
#> Symmetry Axes: 1
#> Symmetry Orders (Cn): 6
#> 
#> -------------------
#> See @atoms paramater for atom positions
#> See @bonds paramater for bond positions
#> See @symmetry_axes for symmetry axes
```

The molecule keeps track of: - all attached symmetry axes
(`@symmetry_axes`), - unique orders present (`@symmetry_axes_orders`), -
and whether any symmetry axes exist (`@contains_symmetry_axes`).

You can also retrieve all axes of a given order:

``` r
fetch_all_symmetry_axes_with_order(molecule, Cn = 6L)
#> [[1]]
#> ===================
#> Symmetry Axis
#> ===================
#> Fold Symmetry (Cn): C6
#> PosA: 0, 0, -1
#> PosB: 0, 0, 1
```

------------------------------------------------------------------------

## Features Summary

- 🧬 `Molecule3D`: robust container for atomic and bonding data  
- 📂 `.mol2` parser (`read_mol2()`)  
- ⚙️ Coordinate manipulation (`translate_molecule_*`,
  `transform_molecule()`)  
- 🧭 Persistent anchors for positioning molecules  
- 🔁 **New:** Rotational symmetry support via `SymAxis` and
  `add_symmetry_axis()`  
- 🧩 Connectivity and geometric utilities (`fetch_eleno_*`,
  `compute_distance_between_atoms()`)

------------------------------------------------------------------------

## Citation

If you use **structures** in research or publications, please cite the
package using:

``` r
citation("structures")
#> To cite package 'structures' in publications use:
#> 
#>   El-Kamand S, Wallis M (????). _structures: Classes and Tools for 3D
#>   Molecular Structures and Symmetry_. R package version 0.0.0.9000,
#>   <https://github.com/selkamand/structures>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {structures: Classes and Tools for 3D Molecular Structures and Symmetry},
#>     author = {Sam El-Kamand and Matthew Wallis},
#>     note = {R package version 0.0.0.9000},
#>     url = {https://github.com/selkamand/structures},
#>   }
```
