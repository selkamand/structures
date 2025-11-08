#' SymAxis: a simple S7 class for a rotational symmetry axis (C_n)
#'
#' @description
#' Represents a proper rotation axis in molecular or solid-state symmetry.
#' The axis is defined by two points in 3D space (`posA`, `posB`) and an
#' integer fold/order `Cn` (e.g., 1, 2, 3, 4, ...).
#'
#' @details
#' - `Cn` is the **order (fold)** of the axis, i.e. rotation by `360/Cn` degrees.
#' - `posA` and `posB` are two distinct points (length-3 numeric vectors) that define
#'   the axis direction in Cartesian coordinates (x, y, z).
#' - The class stores both endpoints as supplied and does not normalize or store a
#'   direction vector explicitly.
#' - A human-readable `label` can be assigned to identify or describe the axis.
#'   This label is optional but useful when managing multiple symmetry axes within
#'   a molecule.
#'
#' @param Cn Numeric (length 1). The axis fold/order (must be \eqn{\ge} 1).
#' @param posA Numeric (length 3). A point on the axis: `c(x, y, z)`.
#' @param posB Numeric (length 3). A different point on the axis: `c(x, y, z)`.
#' @param label Character scalar. An optional user-defined name or description
#'   for the symmetry axis.
#'
#' @return An S7 object of class `SymAxis`.
#'
#' @examples
#' # Minimal valid axis (C3 about the line through (0,0,0) and (0,0,1))
#' ax <- SymAxis(Cn = 3L, posA = c(0, 0, 0), posB = c(0, 0, 1))
#' print(ax)
#'
#' # You can also add a label for clarity
#' SymAxis(Cn = 3L, posA = c(0,0,0), posB = c(0,0,1), label = "principal_z")
#'
#' # C1 is allowed (the trivial/identity axis)
#' SymAxis(Cn = 1L, posA = c(1, 0, 0), posB = c(1, 0, 2))
#'
#' # Coercion to integer for Cn if passed as numeric
#' SymAxis(Cn = 4, posA = c(0, 1, 0), posB = c(0, 2, 0))
#'
#' @export
SymAxis <- S7::new_class(
  name = "SymAxis",
  properties = list(
    label = S7::new_property(
      class = S7::class_character,
      validator = function(value){
        if (length(value) != 1) return(sprintf("Symmetry axis label must be a string, not a character vector of length [%s]", length(value)))
        if(is.na(value)) return("Symmetry axis label must NOT be a missing (NA) value")
        if(nchar(value) == 0) return("Symmetry axis label must NOT be an empty string")
      }

    ),
    Cn = S7::new_property(
      class = S7::class_numeric,
      setter = function(self, value){
        if(is_integerlike(value)) {
          value = as.integer(value)
        }
        S7::prop(self, "Cn") <- value
        return(self)
      },
      validator = function(value) {
        if (length(value) != 1L)
          return(
            sprintf(
              "SymAxis fold symmetry (Cn) must be a single number (length 1 integer), not a vector of length [%s]",
              length(value)
            )
          )
        if (is.na(value) || value < 1L)
          return(sprintf("SymAxis fold symmetry (Cn) must be >= 1. Not [%s]", value))

        if (is.infinite(value))
          return(sprintf("SymAxis fold symmetry (Cn) can not be an infinite value"))

        if (!is_integerlike(value))
          return(sprintf("SymAxis fold symmetry (Cn) must be a whole number, not [%g]", value))
      },
    ),
    posA = S7::new_property(
      class = S7::class_numeric,
      validator = function(value) {
        if (length(value) != 3L)
          return(
            sprintf(
              "posA must be a numeric vector of length 3: c(x, y, z). Supplied length [%s]",
              length(value)
            )
          )
        if (any(!is.finite(value)) || any(is.na(value)))
          return("posA must contain only finite, non-NA values.")
      },
      setter = function(self, value){
        names(value) <- c("x", "y", "z")[seq_along(value)]
        S7::prop(self, "posA") <- value
        return(self)
      }
    ),
    posB = S7::new_property(
      class = S7::class_numeric,
      validator = function(value) {
        if (length(value) != 3L)
          return(
            sprintf(
              "posB must be a numeric vector of length 3: c(x, y, z). Supplied length [%s]",
              length(value)
            )
          )
        if (any(!is.finite(value)) || any(is.na(value)))
          return("posB must contain only finite, non-NA values.")
      },
      setter = function(self, value){
        names(value) <- c("x", "y", "z")[seq_along(value)]
        S7::prop(self, "posB") <- value
        return(self)
      }
    )
  ),
  constructor = function(Cn, posA, posB, label = "unnamed") {

    # Return the S7 object of the correct class
    S7::new_object(
      S7::S7_object(),
      Cn   = Cn,
      posA = posA,
      posB = posB,
      label = label
    )
  },

  # Validate inter-property relationships
  validator = function(self){

    # Check posA is not equal to posB
    if (isTRUE(all(self@posA == self@posB))) {
        return("posA and posB must be distinct points to define an axis.")
    }
  }
)


# Generics ----------------------------------------------------------------
#' @export
S7::method(print, SymAxis) <- function(x, ...) {
  cat(
    "===================\n",
    "Symmetry Axis\n",
    "===================\n",
    sprintf("Fold Symmetry / Order (Cn): C%g\n", x@Cn),
    sprintf("PosA: %s\n", toString(x@posA)),
    sprintf("PosB: %s\n", toString(x@posB)),
    sprintf("Label: %s\n", toString(x@label)),
    sep = ""
  )
  return(invisible(x))
}


#' @export
S7::method(as.data.frame, SymAxis) <- function(x, ...) {
  posA = x@posA
  posB = x@posB
  label = x@label
  Cn = x@Cn

  data.frame(
    label = label,
    Cn = Cn,
    x=posA["x"],
    y=posA["y"],
    z=posA["z"],
    xend=posB["x"],
    yend=posB["y"],
    zend=posB["z"],
    row.names = NULL
  )
}


# Helpers -----------------------------------------------------------------
is_integerlike <- function(x, tol = .Machine$double.eps^0.5){
  if (length(x) == 0L) return(FALSE)
  if (anyNA(x) || any(!is.finite(x))) return(FALSE)
  all(abs(x - round(x)) <= tol)
}



# Non-Generic Methods -----------------------------------------------------

#' Check if an object is a SymAxis
#'
#' Tests whether an object inherits from the [`structures::SymAxis`] class.
#'
#' @param x An object to test.
#'
#' @return A logical scalar: `TRUE` if `x` is a `SymAxis`, otherwise `FALSE`.
#'
#' @examples
#' ax <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
#' is_symmetry_axis(ax)
#' is_symmetry_axis("not_an_axis")
#'
#' @export
is_symmetry_axis <- function(x){
  inherits(x, "structures::SymAxis")
}

#' Apply a 3D transformation to a SymAxis
#'
#' @description
#' Applies a user-supplied transformation function to the two endpoints (`posA`, `posB`)
#' that define a [`SymAxis`] line in 3D. This is useful for translating, rotating, or
#' otherwise mapping the axis to a new position/orientation in space. The symmetry order
#' `Cn` is left unchanged.
#'
#' @param x A [`structures::SymAxis`] object.
#' @param transformation A function that takes a length-3 coordinate and returns a
#'   transformed length-3 coordinate. The input will be passed as a named numeric
#'   vector `c(x=..., y=..., z=...)`. The function must return either:
#'   - a named numeric vector with names `x`, `y`, `z`, or
#'   - a list with elements `x`, `y`, `z` coercible to numeric.
#' @param ... Additional arguments forwarded to `transformation`.
#'
#' @return The same `SymAxis` object `x` with updated `posA` and `posB`.
#'
#' @details
#' The transformation is applied independently to `posA` and `posB`. No attempt is made
#' to renormalize or reorder points; the function simply replaces both endpoints with
#' their transformed coordinates. It is the caller's responsibility to ensure the
#' transformation preserves distinct endpoints and yields finite numeric results.
#'
#' @examples
#' # Define an axis along +Z through the origin
#' ax <- SymAxis(Cn = 3L, posA = c(0,0,0), posB = c(0,0,1))
#'
#' # 1) Translate by (+1, +2, +3)
#' translate <- function(p, dx=0, dy=0, dz=0) {
#'   c(x = p["x"] + dx, y = p["y"] + dy, z = p["z"] + dz)
#' }
#' ax_t <- transform_symmetry_axis(ax, translate, dx = 1, dy = 2, dz = 3)
#'
#' # 2) Rotate around Z by 90 degrees about the origin
#' rotate_z <- function(p, theta) {
#'   c(x =  cos(theta)*p["x"] - sin(theta)*p["y"],
#'     y =  sin(theta)*p["x"] + cos(theta)*p["y"],
#'     z =  p["z"])
#' }
#' ax_r <- transform_symmetry_axis(ax, rotate_z, theta = pi/2)
#'
#' # 3) A transform that returns a list is also accepted
#' as_list <- function(p) list(x = p["x"] + 1, y = p["y"], z = p["z"])
#' ax_l <- transform_symmetry_axis(ax, as_list)
#'
#'@export
transform_symmetry_axis <- function(x, transformation, ...){
  assertions::assert_class(x, "structures::SymAxis")
  assertions::assert_function(transformation)

  posA <- x@posA
  posB <- x@posB

  posAnew <- transformation(posA, ...)
  posBnew <- transformation(posB, ...)

  posAnew <- if(is.list(posAnew)) unlist(posAnew) else posAnew
  posBnew <- if(is.list(posBnew)) unlist(posBnew) else posBnew

  new_axis <- S7::set_props(x, posA=posAnew, posB = posBnew)

  return(new_axis)
}

