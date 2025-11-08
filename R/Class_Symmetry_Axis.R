#' SymAxis: a simple S7 class for a rotational symmetry axis (C_n)
#'
#' @description
#' Represents a proper rotation axis in molecular/solid-state symmetry.
#' The axis is specified by two points in 3D space (`posA`, `posB`) and an
#' integer fold/order `Cn` (e.g., 1, 2, 3, 4, ...).
#'
#' @details
#' - `Cn` is the **order (fold)** of the axis, i.e. rotation by `360/Cn` degrees.
#' - `posA` and `posB` are two distinct points (length-3 numeric vectors) that define
#'   the axis direction in the order x, y, z.
#' - This class does not normalize or store a direction vector; it stores the points as given.
#'
#' @param Cn Numeric (length 1). The axis fold/order (must be \eqn{\ge} 1).
#' @param posA Numeric (length 3). A point on the axis: c(x, y, z).
#' @param posB Numeric (length 3). A different point on the axis: c(x, y, z).
#'
#' @return An S7 object of class `SymAxis`.
#'
#' @examples
#' # Minimal valid axis (C3 about the line through (0,0,0) and (0,0,1))
#' ax <- SymAxis(Cn = 3L, posA = c(0, 0, 0), posB = c(0, 0, 1))
#' print(ax)
#'
#' # C1 is allowed (the trivial/identity axis)
#' SymAxis(Cn = 1L, posA = c(1, 0, 0), posB = c(1, 0, 2))
#'
#' # Coercion to integer for Cn if passed as numeric
#' SymAxis(Cn = 4, posA = c(0, 1, 0), posB = c(0, 2, 0))
#'
#'
#' @export
SymAxis <- S7::new_class(
  name = "SymAxis",
  properties = list(
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

      }
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
      }
    )
  ),
  constructor = function(Cn, posA, posB) {

    # Return the S7 object of the correct class
    S7::new_object(
      S7::S7_object(),
      Cn   = Cn,
      posA = posA,
      posB = posB
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
    sprintf("Fold Symmetry (Cn): C%g\n", x@Cn),
    sprintf("PosA: %s\n", toString(x@posA)),
    sprintf("PosB: %s\n", toString(x@posB)),
    sep = ""
  )
}


# Helpers -----------------------------------------------------------------
is_integerlike <- function(x, tol = .Machine$double.eps^0.5){
  if (length(x) == 0L) return(FALSE)
  if (anyNA(x) || any(!is.finite(x))) return(FALSE)
  all(abs(x - round(x)) <= tol)
}
