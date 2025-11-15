
# ProperRotationAxis Class Definition -------------------------------------


#' ProperRotationAxis: Rotational Symmetry Axis (C_n)
#'
#' The ProperRotationAxis object represents a proper rotation axis in molecular or solid-state symmetry.
#' The axis is defined by a point \code{posA}, an orientation \code{direction}
#' (unit vector). A secondary point on this defined axis (\code{posB}) can also be
#' derived as \code{posA + direction * L} and is read-only.
#'
#' @details
#' - \code{Cn} is the order (fold) of the axis, i.e. rotation by \code{360/Cn} degrees.
#' - Construction options:
#'   * Supply \code{posA} and \code{posB}: \code{direction} and \code{L} are inferred.
#'   * Supply \code{posA} and \code{direction}: \code{direction} is normalized; if
#'     \code{L} is missing it defaults to \code{1}.
#' - Setting \code{direction} normalizes it but does not change \code{L}; change \code{L}
#'   to adjust the span from \code{posA} to \code{posB}.
#'
#' @param Cn Numeric (length 1). Axis fold/order (\eqn{\ge} 1, whole number).
#' @param posA Numeric (length 3). A point on the axis: \code{c(x,y,z)}.
#' @param posB Numeric (length 3). A different point on the axis; mutually exclusive
#'   with \code{direction}. When supplied, \code{direction} and \code{L} are inferred.
#' @param direction Numeric (length 3). Orientation vector; normalized on set.
#'   Mutually exclusive with \code{posB}. If used without \code{L}, \code{L} defaults to 1.
#' @param L Numeric (length 1). Positive length from \code{posA} to \code{posB}.
#'   Defaults to \code{1} when \code{direction} is supplied and \code{L} is omitted.
#' @param label Character scalar. Optional user-defined name/description.
#'
#' @return An S7 object of class \code{ProperRotationAxis}.
#'
#' @note \code{posB} is read-only and computed as \code{posA + direction * L}.
#'
#' @examples
#' # Define via two points
#' ProperRotationAxis(Cn = 3L, posA = c(0,0,0), posB = c(0,0,1))
#'
#' # Define via direction; L defaults to 1
#' ProperRotationAxis(Cn = 2L, posA = c(0,0,0), direction = c(0,0,5))
#'
#' # Direction + explicit L
#' ProperRotationAxis(Cn = 2L, posA = c(1,0,0), direction = c(0,0,1), L = 2)
#'
#' @export
ProperRotationAxis <- S7::new_class(
  name = "ProperRotationAxis",
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
              "ProperRotationAxis fold symmetry (Cn) must be a single number (length 1 integer), not a vector of length [%s]",
              length(value)
            )
          )
        if (is.na(value) || value < 1L)
          return(sprintf("ProperRotationAxis fold symmetry (Cn) must be >= 1. Not [%s]", value))

        if (is.infinite(value))
          return(sprintf("ProperRotationAxis fold symmetry (Cn) can not be an infinite value"))

        if (!is_integerlike(value))
          return(sprintf("ProperRotationAxis fold symmetry (Cn) must be a whole number, not [%g]", value))
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

    # Never set directly
    posB = S7::new_property(
      class = S7::class_numeric,
      setter = function(self, value){
        stop("@posB is read-only; set @direction or @L instead")
      },
      getter = function(self){
        out <- self@posA + self@direction * self@L
        names(out) <- c("x","y","z")
        return(out)
      },
    ),



    # A unit direction vector representing the axis (computed as direction from posA to posB)
    direction = S7::new_property(
      class = S7::class_numeric,
      validator = function(value) {
        if (length(value) != 3L)
          return(
            sprintf(
              "direction must be a numeric vector of length 3: c(x, y, z). Supplied length [%s]",
              length(value)
            )
          )
        if (any(!is.finite(value)) || any(is.na(value)))
          return("direction must contain only finite, non-NA values.")
      },
      getter = function(self){ names(self@direction) <- c("x", "y", "z"); return(self@direction) },
      setter = function(self, value){
        value = normalise(value)
        self@direction <- value
        return(self)
      }
    ),

    # Length between posA and posB
    L = S7::new_property(
      class = S7::class_numeric,
      validator = function(value){
        if (length(value) != 1L) return("L must be a length-1 numeric")
        if (!is.finite(value) || is.na(value)) return("L must be finite")
        if (value <= 0) return("L must be > 0 to define an axis")
      }
    )
  ),

  constructor = function(Cn, posA, posB = NULL, direction = NULL, L = NULL, label = "unnamed"){
    if (is.null(posB) && is.null(direction)) stop("Supply either `posB` or `direction`")
    if (!is.null(posB) && !is.null(direction)) stop("Supply only one of `posB` or `direction`")
    assertions::assert_length(posA, length = 3, msg = "posA must be a numeric vector of length 3: c(x, y, z)")
    assertions::assert_no_missing(posA)
    assertions::assert(all(is.finite(posA)), msg = "posA must have no infinite values")

    names(posA) <- c("x","y","z")[seq_along(posA)]

    if (!is.null(posB)) {
      assertions::assert_length(posB, length = 3, msg = "posB must be a numeric vector of length 3: c(x, y, z)")
      assertions::assert_no_missing(posB)
      assertions::assert(all(is.finite(posB)), msg = "posB must have no infinite values")

      names(posB) <- c("x","y","z")[seq_along(posB)]
      vec <- posB - posA
      len <- sqrt(sum(vec^2))
      if (len <= 0) stop("posA and posB must be distinct points")
      direction <- vec / len
      L <- len
    } else {
      assertions::assert_length(direction, length = 3, msg = "direction must be a numeric vector of length 3: c(x, y, z)")
      assertions::assert_no_missing(direction)
      assertions::assert(all(is.finite(direction)), msg = "direction must have no infinite values")
      # direction provided
      direction <- normalise(direction)
      L <- if (is.null(L)) 1 else L
    }
    names(direction) <- c("x","y","z")

    S7::new_object(
      S7::S7_object(),
      Cn   = Cn,
      posA = posA,
      direction = direction,
      L = L,
      label = label
    )
  },

  # Validate inter-property relationships
  validator = function(self){
    # With derivation, this ensures L > 0 and endpoints differ
    if (identical(self@posA, self@posB))
      return("posA and posB must be distinct (L must be > 0)")
      # browser()
    NULL
  }
)


# Generics ----------------------------------------------------------------
#' @export
S7::method(print, ProperRotationAxis) <- function(x, ...) {
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
S7::method(as.data.frame, ProperRotationAxis) <- function(x, ...) {
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

#' Apply a 3D transformation to a ProperRotationAxis
#'
#' @description
#' Applies a user-supplied transformation function to the two endpoints (`posA`, `posB`)
#' that define a [`ProperRotationAxis`] line in 3D. This is useful for translating, rotating, or
#' otherwise mapping the axis to a new position/orientation in space. The symmetry order
#' `Cn` is left unchanged.
#'
#' @param x A [`structures::ProperRotationAxis`] object.
#' @param transformation A function that takes a length-3 coordinate and returns a
#'   transformed length-3 coordinate. The input will be passed as a named numeric
#'   vector `c(x=..., y=..., z=...)`. The function must return either:
#'   - a named numeric vector with names `x`, `y`, `z`, or
#'   - a list with elements `x`, `y`, `z` coercible to numeric.
#' @param ... Additional arguments forwarded to `transformation`.
#'
#' @return The same `ProperRotationAxis` object `x` with updated `posA` and `posB`.
#'
#' @details
#' The transformation is applied independently to `posA` and `posB`. No attempt is made
#' to renormalize or reorder points; the function simply replaces both endpoints with
#' their transformed coordinates. It is the caller's responsibility to ensure the
#' transformation preserves distinct endpoints and yields finite numeric results.
#'
#' @examples
#' # Define an axis along +Z through the origin
#' ax <- ProperRotationAxis(Cn = 3L, posA = c(0,0,0), posB = c(0,0,1))
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
  assertions::assert_class(x, c("structures::ProperRotationAxis", "structures::ImproperRotationAxis"))
  assertions::assert_function(transformation)

  posA <- x@posA
  posB <- x@posB

  posAnew <- transformation(posA, ...)
  posBnew <- transformation(posB, ...)

  posAnew <- if(is.list(posAnew)) unlist(posAnew) else posAnew
  posBnew <- if(is.list(posBnew)) unlist(posBnew) else posBnew
  direction = create_vector_from_start_end(start = posAnew, end = posBnew, unit = TRUE)
  L = compute_distance(posAnew, posBnew)

  new_axis <- S7::set_props(x, posA=posAnew, direction = direction, L=L)

  return(new_axis)
}

