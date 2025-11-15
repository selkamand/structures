# Symmetry Element (Abstract Base Class) -------------------------------------------------------
#' SymmetryElement
#'
#' Abstract parent class for all symmetry elements used in point groups and
#' molecular symmetry. Concrete subclasses include
#' [ProperRotationAxis()], [MirrorPlane()], [CentreOfInversion()],
#' and [ImproperRotationAxis()]. All symmetry elements share a user-facing
#' `label` and a derived `type`.
#'
#' @param label Character scalar. Optional user-defined name or description
#'   for the symmetry element.
#' @param type Character scalar (read-only). A derived string describing the
#'   element type. Must be one of
#'   `"Mirror Plane"`, `"Centre of Inversion"`,
#'   `"Improper Rotation Axis"`, `"Proper Rotation Axis"`.
#'
#' @return An S7 object of class `SymmetryElement` (or a subclass).
#'
#' @examples
#' # Typically you construct concrete subclasses, not SymmetryElement directly:
#' mp <- MirrorPlane(normal = c(0, 0, 1), position = c(0, 0, 0), label = "σ_xy")
#' mp@label
#' mp@type   # "Mirror Plane"
#'
#' @export
SymmetryElement <- S7::new_class(
  name = "SymmetryElement",
  properties = list(
    label = S7::new_property(
      class = S7::class_character,
      validator = function(value){
        if (length(value) != 1) return(sprintf("Symmetry element @label must be a string, not a character vector of length [%s]", length(value)))
        if(is.na(value)) return("Symmetry element @label must NOT be a missing (NA) value")
        if(nchar(value) == 0) return("Symmetry element @label must NOT be an empty string")
      }
    ),
    type = S7::new_property(
      class = S7::class_character,
      setter = function(self, value) { stop("@type is a read only property") },
      getter = function(self){ identify_symmetry_type(self) },
      validator =  function(value) {

        if(length(value) != 1)
          return(sprintf("Symmetry element @type must be a string, not a character vector of length [%s]. If you're seeing this message there is probably an error in the identify_symmetry_type function. Please alert package authors", length(value)))

        valid_symmetry_element_types <- c("Mirror Plane", "Centre of Inversion", "Improper Rotation Axis", "Proper Rotation Axis")
        if(!value %in% valid_symmetry_element_types){
          return(sprintf(
            "Symmetry element @type must be one of [%s], not [%s]. If you're seeing this message there is probably an error in the identify_symmetry_type function. Please alert package authors",
            toString(valid_symmetry_element_types), toString(value)
          ))
        }

        return(NULL)
      }
    )
  )

  # No properties. Just useful as parent to all other Symmetry elements
)

# Symmetry Element Collection -------------------------------------------------------
#
# adn allowing jus to dispatch generics like as.data.frame to the Symmetry Element List

#' SymmetryElementCollection
#'
#'A list of symmetry elements. Important for managing IDs of each element (ensuring they're unique)
#' Forces all symmetry elements to have label and type properties (see Arguments section).
#'
#' @param label Character scalar. Optional user-defined name/description.
#' @param type Character scalar. Read only. A string describing the element type. Must be one of  "Mirror Plane", "Centre of Inversion", "Improper Rotation Axis", "Proper Rotation Axis"
#'
SymmetryElementCollection <- S7::new_class(
  name = "SymmetryElementCollection",
  properties = list(
    symmetry_elements = S7::class_list,
    ids = S7::class_numeric,
    counts = S7::new_property(
      class = S7::class_data.frame,
      setter = function(self, value) { stop("@counts is a read only property") },
      getter = function(self) {
        warning("@counts paramater is not yet implemented.")
      }
    )
  ),
  validator = function(self){
    ids=self@ids
    elements=self@elements

    if(length(ids) != length(elements))
      return(sprintf("@ids must have one entry for each of the @symmetry_elements. Found [%s] IDs and [%s] symmetry_elements", length(ids),length(elements)))

    if(any(duplicated(ids)))
      return(sprintf("Duplicate @ids are not allowed. Found [%s] duplicates", sum(duplicated(ids))))

    if(anyNA(ids))
      return(sprintf("Missing (NA) @ids are not allowed. Found [%s].", sum(is.na(ids))))

    for (element in elements){
      if(!inherits(element, "SymmetryElement")) return(sprintf("All @symmetry_elements must inherit from: structures::SymmetryElement"))
    }

    return(NULL)
  }
)


## Generics ----------------------------------------------------------------
# S7::method(print, SymmetryElements) <- function(x, ...) {
#   symmetry_orders_string <- if(is.null(x@symmetry_axes_orders)) "" else sprintf("Symmetry Orders (Cn): %s\n", toString(x@symmetry_axes_orders))
#   cat(sep = "",
#       "===================\n",
#       "Chemical Molecule3D\n",
#       "===================\n",
#       sprintf("Name: %s\n", x@name),
#       sprintf("Atoms: %d\n", nrow(x@atoms)),
#       sprintf("Bonds: %d\n", nrow(x@bonds)),
#       sprintf("Symmetry Axes: %d\n", length(x@symmetry_axes)),
#       symmetry_orders_string,
#       "\n-------------------\n",
#       "See @atoms paramater for atom positions\n",
#       "See @bonds paramater for bond positions\n",
#       "See @symmetry_axes for symmetry axes"
#   )
#
#   return(invisible(x))
# }

# as.data.frame <- S7::new_generic("as.data.frame", "x")

#' @export
S7::method(as.data.frame, SymmetryElementCollection) <- function(x, ...) {
  ls_info <- lapply(x, function(x){
    data.frame(
      type = x@type,
      label = x@label,
    )
  })

  df_data <- as.data.frame(do.call(rbind, ls_info))
  df_data$ids <- x@ids
  df_data <- df_data[c("ids", setdiff(colnames(df_data), "ids"))]

  return(df_data)
}

# # as.matrix <- S7::new_generic("as.matrix", "x")
# S7::method(as.matrix, Molecule3D) <- function(x, ...) {
#   mx <- as.matrix(x@atoms[c("x", "y", "z")])
#   rownames(mx) <- x@atoms[["eleno"]]
#   return(mx)
# }



# Create an Class for each symmetry element that inherits from 'SymmetryElement'


# ProperRotationAxis  -------------------------------------

## Class  -------------------------------------
#' ProperRotationAxis: rotational symmetry axis (C_n)
#'
#' Represents a proper rotation axis (C\eqn{_n}) in 3D. The axis is defined by a
#' point `posA` and an orientation `direction` (unit vector). A second point
#' `posB` lies along the axis at distance `L` from `posA` and is derived as
#' `posA + direction * L`.
#'
#' @details
#' * `Cn` is the order (fold) of the axis, i.e. rotation by `360 / Cn` degrees.
#' * You can construct an axis either by:
#'   - supplying `posA` and `posB` (direction and `L` are inferred), or
#'   - supplying `posA` and `direction` (direction is normalised; `L` defaults to 1
#'     if not given).
#' * `posB` is read-only and always computed from `posA`, `direction`, and `L`.
#'
#' @param Cn Numeric scalar. Axis fold/order (whole number >= 1).
#' @param posA Numeric length-3 vector giving a point on the axis (`c(x, y, z)`).
#' @param posB Optional numeric length-3 vector giving a second point on the axis.
#'   Mutually exclusive with `direction`.
#' @param direction Optional numeric length-3 vector giving the axis direction.
#'   Normalised on set. Mutually exclusive with `posB`.
#' @param L Numeric scalar. Positive distance between `posA` and `posB`. If
#'   `direction` is supplied and `L` is omitted, `L` defaults to 1.
#' @inheritParams SymmetryElement
#'
#' @return An S7 object of class `ProperRotationAxis`.
#'
#' @examples
#' # Define via two points (C3 axis along +Z through the origin)
#' ax1 <- ProperRotationAxis(
#'   Cn   = 3L,
#'   posA = c(0, 0, 0),
#'   posB = c(0, 0, 1),
#'   label = "C3(z)"
#' )
#'
#' # Define via direction (C2 axis along +Z; L defaults to 1)
#' ax2 <- ProperRotationAxis(
#'   Cn        = 2L,
#'   posA      = c(0, 0, 0),
#'   direction = c(0, 0, 5),
#'   label     = "C2(z)"
#' )
#'
#' ax1@type   # "Proper Rotation Axis"
#'
#' @export
ProperRotationAxis <- S7::new_class(
  name = "ProperRotationAxis",
  parent = SymmetryElement,
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

  # Validate inter-property relationships
  validator = function(self){
    # With derivation, this ensures L > 0 and endpoints differ
    if (identical(self@posA, self@posB))
      return("posA and posB must be distinct (L must be > 0)")
    # browser()
    NULL
  },
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
  }
)


## Generics ----------------------------------------------------------------
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


## Helpers -----------------------------------------------------------------
is_integerlike <- function(x, tol = .Machine$double.eps^0.5){
  if (length(x) == 0L) return(FALSE)
  if (anyNA(x) || any(!is.finite(x))) return(FALSE)
  all(abs(x - round(x)) <= tol)
}



## Non-Generic Methods -----------------------------------------------------

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

# ImproperRotationAxis ----------------------------------------------------

## Class -------------------------------------------------------------

#' ImproperRotationAxis: improper rotational symmetry axis (S_n)
#'
#' Represents an improper rotation axis (S\eqn{_n}) in 3D, combining a proper
#' rotation (C\eqn{_n}) about the axis with reflection in a plane perpendicular
#' to that axis. Geometrically, the axis is described exactly as in
#' [`ProperRotationAxis()`] (via `Cn`, `posA`, `direction`, `L`), with an
#' additional point `plane_point` that fixes the location of the mirror plane.
#'
#' The mirror plane is defined as the set of points `x` satisfying
#' `(x - plane_point) · direction = 0`, i.e. a plane with normal equal to
#' the axis direction, passing through `plane_point`. In typical molecular
#' symmetry, `plane_point` will also lie on the axis (e.g. at `posA`).
#'
#' @param Cn Numeric scalar. Axis fold/order (whole number >= 1).
#' @param posA Numeric length-3 vector giving a point on the axis (`c(x, y, z)`).
#' @param posB Optional numeric length-3 vector giving a second point on the axis.
#'   Mutually exclusive with `direction`. When supplied, `direction` and `L` are
#'   inferred from `posA` and `posB`.
#' @param direction Optional numeric length-3 vector giving the axis direction.
#'   Normalised on set. Mutually exclusive with `posB`. If used without `L`,
#'   `L` defaults to 1.
#' @param L Numeric scalar. Positive distance between `posA` and `posB`. If
#'   `direction` is supplied and `L` is omitted, `L` defaults to 1.
#' @param plane_point Numeric length-3 vector giving a point in the mirror
#'   plane. The plane normal is taken to be the axis direction. By default this
#'   is `posA`, i.e. the plane cuts the axis at `posA`.
#' @inheritParams SymmetryElement
#'
#' @return An S7 object of class `ImproperRotationAxis`.
#'
#' @examples
#' # S4 axis along +Z through the origin; mirror plane perpendicular to the axis
#' # and cutting it at the origin (plane_point = posA).
#' iax1 <- ImproperRotationAxis(
#'   Cn          = 4L,
#'   posA        = c(0, 0, 0),
#'   posB        = c(0, 0, 1),
#'   plane_point = c(0, 0, 0),
#'   label       = "S4(z)"
#' )
#'
#' iax1@type   # "Improper Rotation Axis"
#'
#' # S6 axis along +Z, with the mirror plane shifted along the axis
#' iax2 <- ImproperRotationAxis(
#'   Cn          = 6L,
#'   posA        = c(0, 0, 0),
#'   direction   = c(0, 0, 1),
#'   L           = 2,
#'   plane_point = c(0, 0, 0.5),
#'   label       = "S6(z)"
#' )
#'
#' @export
ImproperRotationAxis <- S7::new_class(
  parent = ProperRotationAxis,
  name   = "ImproperRotationAxis",
  properties = list(
    plane_point = S7::new_property(
      class = S7::class_numeric,
      validator = function(value){
        if (length(value) != 3L)
          return(
            sprintf(
              "ImproperRotationAxis @plane_point must be a numeric vector of length 3: not [%s]",
              length(value)
            )
          )
        if (any(!is.finite(value)) || any(is.na(value)))
          return("ImproperRotationAxis @plane_point must contain only finite, non-NA values.")
        NULL
      },
      setter = function(self, value){
        names(value) <- c("x", "y", "z")[seq_along(value)]
        S7::prop(self, "plane_point") <- value
        self
      }
    )
  ),

  constructor = function(Cn,
                         posA,
                         plane_point = posA,
                         posB      = NULL,
                         direction = NULL,
                         L         = NULL,
                         label     = "unnamed") {

    if (is.null(posB) && is.null(direction))
      stop("Supply either `posB` or `direction`")
    if (!is.null(posB) && !is.null(direction))
      stop("Supply only one of `posB` or `direction`")

    # --- Axis definition (as in ProperRotationAxis) ---
    assertions::assert_length(posA, length = 3, msg = "posA must be a numeric vector of length 3: c(x, y, z)")
    assertions::assert_no_missing(posA)
    assertions::assert(all(is.finite(posA)), msg = "posA must have no infinite values")
    names(posA) <- c("x", "y", "z")[seq_along(posA)]

    if (!is.null(posB)) {
      assertions::assert_length(posB, length = 3, msg = "posB must be a numeric vector of length 3: c(x, y, z)")
      assertions::assert_no_missing(posB)
      assertions::assert(all(is.finite(posB)), msg = "posB must have no infinite values")

      names(posB) <- c("x", "y", "z")[seq_along(posB)]
      vec <- posB - posA
      len <- sqrt(sum(vec^2))
      if (len <= 0) stop("posA and posB must be distinct points")
      direction <- vec / len
      L <- len
    } else {
      assertions::assert_length(direction, length = 3, msg = "direction must be a numeric vector of length 3: c(x, y, z)")
      assertions::assert_no_missing(direction)
      assertions::assert(all(is.finite(direction)), msg = "direction must have no infinite values")
      direction <- normalise(direction)
      L <- if (is.null(L)) 1 else L
    }
    names(direction) <- c("x", "y", "z")

    # --- Mirror plane anchor ---
    assertions::assert_length(plane_point, length = 3,
                              msg = "plane_point must be a numeric vector of length 3: c(x, y, z)")
    assertions::assert_no_missing(plane_point)
    assertions::assert(all(is.finite(plane_point)),
                       msg = "plane_point must have no infinite values")
    names(plane_point) <- c("x", "y", "z")[seq_along(plane_point)]

    # Enforce that plane_point lies on the axis
    v <- plane_point - posA
    if (magnitude(v) > 0) {
      if (max(abs(normalise(v) - direction)) > 1e-6)
        stop("plane_point must lie on the rotation axis (collinear with direction from posA).")
    }

    S7::new_object(
      S7::S7_object(),
      Cn = Cn,
      posA = posA,
      direction = direction,
      L = L,
      plane_point = plane_point,
      label = label
    )
  }
)


# MirrorPlane -------------------------------------------------------------

## Class -------------------------------------------------------------

#' MirrorPlane: reflection symmetry element (σ)
#'
#' Represents a mirror plane in 3D space, defined by a unit normal vector and
#' a point lying on the plane. The plane contains all points `p` such that
#' `(p - position) · normal = 0`.
#'
#' @param normal Numeric length-3 vector giving the plane's normal. Must have
#'   non-zero magnitude.
#' @param position Numeric length-3 vector giving a point on the plane
#'   (`c(x, y, z)`).
#' @inheritParams SymmetryElement
#'
#' @return An S7 object of class `MirrorPlane`.
#'
#' @examples
#' # Mirror plane xy: normal along +Z, passing through the origin
#' mp <- MirrorPlane(
#'   normal   = c(0, 0, 1),
#'   position = c(0, 0, 0),
#'   label = "σ_xy"
#' )
#'
#' mp@type    # "Mirror Plane"
#'
#' @export
MirrorPlane  <-  S7::new_class(
  parent = SymmetryElement,
  name = "MirrorPlane",
  properties = list(
    normal = S7::new_property(
      class = S7::class_numeric,
      validator = function(value){
        tol=1e-18
        if(magnitude(value) <= tol) return(sprintf("Mirror plane @normal must be a vector with magnitude > 0, not [%d]", magnitude(value)))
        if(length(value) != 3) return(sprintf("Mirror plane @normal must be a numeric vector with length 3: Not [%s]", length(value)))
        if(length(value) != 3) return(sprintf("Mirror plane @normal must be a numeric vector with length 3: Not [%s]", length(value)))
        return(NULL)
      }
    ),
    position = S7::new_property(
      class = S7::class_numeric,
      validator = function(value){
        if(length(value) != 3) return(sprintf("Mirror plane @position must be a numeric vector with length 3: Not [%s]", length(value)))
        return(NULL)
      }
    )
  ),
  constructor = function(normal, position, label = "unnamed"){
    S7::new_object(
      S7::S7_object(),
      normal = normal,
      position = position,
      label = label
    )
  }
)


# CentreOfInversion -------------------------------------------------------

## Class -------------------------------------------------------------

#' CentreOfInversion: inversion symmetry element (i)
#'
#' Represents a centre of inversion in 3D space. Every point `p` is mapped to
#' `p' = 2 * position - p` under inversion through this centre.
#'
#' @param position Numeric length-3 vector giving the inversion centre
#'   (`c(x, y, z)`).
#' @inheritParams SymmetryElement
#'
#' @return An S7 object of class `CentreOfInversion`.
#'
#' @examples
#' # Inversion centre at the origin
#' ci <- CentreOfInversion(
#'   position = c(0, 0, 0),
#'   label = "i"
#' )
#'
#' ci@type    # "Centre of Inversion"
#'
#' @export
CentreOfInversion <- S7::new_class(
  parent = SymmetryElement,
  name = "CentreOfInversion",
  properties = list(
    position = S7::new_property(
      class = S7::class_numeric,
      validator = function(value){
        if(length(value) != 3) return(sprintf("Center of Inversion @position must be a numeric vector with length 3: Not [%s]", length(value)))
        return(NULL)
      }
    )
  ),
  constructor = function(position, label = "unnamed"){
    S7::new_object(
      S7::S7_object(),
      position = position,
      label = label
    )
  }
)






# Class Checkers ----------------------------------------------------------

#' Check if an object is a MirrorPlane
#'
#' @param x An object to test.
#'
#' @return Logical. `TRUE` if `x` inherits from [`structures::MirrorPlane`], otherwise `FALSE`.
#'
#' @examples
#' mp <- MirrorPlane(normal = c(0,0,1), position = c(0,0,0))
#' is_mirror_plane(mp)        # TRUE
#' is_mirror_plane("nope")    # FALSE
#'
#' @export
is_mirror_plane <- function(x){
  inherits(x, "structures::MirrorPlane")
}

#' Check if an object is a CentreOfInversion
#'
#' @param x An object to test.
#'
#' @return Logical. `TRUE` if `x` inherits from [`structures::CentreOfInversion`], otherwise `FALSE`.
#'
#' @examples
#' ci <- CentreOfInversion(position = c(0,0,0))
#' is_centre_of_inversion(ci)     # TRUE
#' is_centre_of_inversion(123)    # FALSE
#'
#' @export
is_centre_of_inversion <- function(x){
  inherits(x, "structures::CentreOfInversion")
}

#' Check if an object is an ImproperRotationAxis
#'
#' @param x An object to test.
#'
#' @return Logical. `TRUE` if `x` inherits from [`structures::ImproperRotationAxis`], otherwise `FALSE`.
#'
#' @examples
#' # Assuming ImproperRotationAxis() will have a constructor similar to other elements:
#' iax <- ImproperRotationAxis()
#' is_improper_rotation_axis(iax)   # TRUE
#' is_improper_rotation_axis(NULL)  # FALSE
#'
#' @export
is_improper_rotation_axis <- function(x){
  inherits(x, "structures::ImproperRotationAxis")
}

#' Check if an object is a ProperRotationAxis
#'
#' @param x An object to test.
#'
#' @return Logical. `TRUE` if `x` inherits from [`structures::ProperRotationAxis`]
#'   (and not from `ImproperRotationAxis`), otherwise `FALSE`.
#'
#' @examples
#' ax <- ProperRotationAxis(Cn = 3L, posA = c(0,0,0), posB = c(0,0,1))
#' is_proper_rotation_axis(ax)        # TRUE
#' is_proper_rotation_axis("nope")    # FALSE
#'
#' @export
is_proper_rotation_axis <- function(x){
  inherits(x, "structures::ProperRotationAxis") & !inherits(x, "structures::ImproperRotationAxis")
}


#' Check if an object is a symmetry axis
#'
#' Tests whether an object inherits from the [`structures::ProperRotationAxis`] or [`structures::ImproperRotationAxis`] class.
#'
#' @param x An object to test.
#'
#' @return A logical scalar: `TRUE` if `x` is a `ProperRotationAxis`, otherwise `FALSE`.
#'
#' @examples
#' ax <- ProperRotationAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
#' is_symmetry_axis(ax)
#' is_symmetry_axis("not_an_axis")
#'
#' @export
is_symmetry_axis <- function(x){
  inherits(x, "structures::ProperRotationAxis") | inherits(x, "structures::ImproperRotationAxis")
}

is_symmetry_element <- function(x){
  inherits(x, "structures::SymmetryElement")
}

identify_symmetry_type <- function(x, strict = TRUE){
  if(is_mirror_plane(x))
    return("Mirror Plane")
  else if(is_centre_of_inversion(x))
    return("Centre of Inversion")
  else if(is_improper_rotation_axis(x))
    return("Improper Rotation Axis")
  else if(is_proper_rotation_axis(x))
    return("Proper Rotation Axis")

  # else if(is_symmetry_element(x))
  #   return("Other Symmetry element")

  warn_message = sprintf("Object is not a known symmetry element. Class is: [%s]", toString(class(x)))
  if(strict) stop(warn_message)
  else warning(warn_message)

  return("Unknown")
}

# Helpers -----------------------------------------------------------------
symmetry_element_abbreviations <- function(){
c(
    "Cn" =  "Proper Rotation Axes",
    "Sn" =  "Improper Rotation Axes",
    "sigma" = "Mirror Plane",
    "i" =  "Centre of Inversion",
  )

}




