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
#'   element type. See See [valid_symmetry_element_types()] for valid values.
#'
#' @return An S7 object of class `SymmetryElement` (or a subclass).
#'
#' @examples
#' # Typically you construct concrete subclasses, not SymmetryElement directly:
#' mp <- MirrorPlane(normal = c(0, 0, 1), position = c(0, 0, 0), label = "<U+03C3>_xy")
#' mp@label
#' mp@type # "Mirror Plane"
#'
#' @export
SymmetryElement <- S7::new_class(
  name = "SymmetryElement",
  properties = list(
    label = S7::new_property(
      class = S7::class_character,
      validator = function(value) {
        if (length(value) != 1) {
          return(sprintf("Symmetry element @label must be a string, not a character vector of length [%s]", length(value)))
        }
        if (is.na(value)) {
          return("Symmetry element @label must NOT be a missing (NA) value")
        }
        if (nchar(value) == 0) {
          return("Symmetry element @label must NOT be an empty string")
        }
      }
    ),
    type = S7::new_property(
      class = S7::class_character,
      setter = function(self, value) {
        stop("@type is a read only property")
      },
      getter = function(self) {
        identify_symmetry_type(self)
      },
      validator = function(value) {
        if (length(value) != 1) {
          return(sprintf("Symmetry element @type must be a string, not a character vector of length [%s]. If you're seeing this message there is probably an error in the identify_symmetry_type function. Please alert package authors", length(value)))
        }

        valid_symmetry_element_types <- valid_symmetry_element_types()
        if (!value %in% valid_symmetry_element_types) {
          return(sprintf(
            "Symmetry element @type must be one of [%s], not [%s]. If you're seeing this message there is probably an error in the identify_symmetry_type function. Please alert package authors",
            toString(valid_symmetry_element_types), toString(value)
          ))
        }

        return(NULL)
      }
    )
  )
)


## Generics ----------------------------------------------------------------

#' Apply a 3D transformation to a symmetry element
#'
#' @description
#' Generic wrapper that applies a coordinate transformation to any
#' [`structures::SymmetryElement`] subclass. This keeps symmetry elements
#' consistent with transformed molecular coordinates (see [transform_molecule()]).
#'
#' Concrete methods are provided for:
#' \itemize{
#'   \item [`structures::ProperRotationAxis`]
#'   \item [`structures::ImproperRotationAxis`]
#'   \item [`structures::MirrorPlane`]
#'   \item [`structures::CentreOfInversion`]
#' }
#'
#' @param x A [`structures::SymmetryElement`] object.
#' @param transformation A function that takes a length-3 coordinate and returns a
#'   transformed length-3 coordinate. The input will be passed as a named numeric
#'   vector `c(x = ..., y = ..., z = ...)`. The function must return either:
#'   \itemize{
#'     \item a named numeric vector with names `x`, `y`, `z`, or
#'     \item a list with elements `x`, `y`, `z` coercible to numeric.
#'   }
#' @param ... Additional arguments forwarded to `transformation`.
#'
#' @return A transformed symmetry element of the same class as `x`.
#'
#' @export
transform_symmetry_element <- S7::new_generic(
  name = "transform_symmetry_element", "x"
)

# If we end up adding more SymmetryElement subclasses later and forget to define a transform_symmetry_element methods
# the code below will force a loud error
#' @export
S7::method(transform_symmetry_element, SymmetryElement) <- function(x, transformation, ...) {
  stop(
    "No transform_symmetry_element() method implemented for symmetry element of class: ",
    paste(class(x), collapse = " / ")
  )
}
# Symmetry Element Collection -------------------------------------------------------


#' SymmetryElementCollection: a container for symmetry elements
#'
#' A lightweight container that groups multiple [`SymmetryElement`] objects
#' together with numeric IDs. This is useful for storing the symmetry content
#' of a molecule or point group and for coercing the collection to a data frame.
#'
#' @param elements A list of objects inheriting from
#'   [`structures::SymmetryElement`].
#' @param ids Numeric vector of the same length as `elements`,
#'   giving a unique identifier for each element.
#' @param n_elements Integer: Read-only, How many elements are in this collection.
#'
#' @return An S7 object of class `SymmetryElementCollection`.
#'
#' @examples
#' mp <- MirrorPlane(normal = c(0, 0, 1), position = c(0, 0, 0), label = "<U+03C3>_xy")
#' ci <- CentreOfInversion(position = c(0, 0, 0), label = "i")
#'
#' coll <- SymmetryElementCollection(
#'   elements = list(mp, ci),
#'   ids = c(1, 2)
#' )
#'
#' # Summarise elements as a data frame
#' as.data.frame(coll)
#'
#' @export
SymmetryElementCollection <- S7::new_class(
  name = "SymmetryElementCollection",
  properties = list(

    ### Writable --------------------------------------------------------------
    elements = S7::class_list,
    ids = S7::class_numeric,

    ### Read Only --------------------------------------------------------------
    n_elements = S7::new_property(
      class = S7::class_integer,
      setter = function(self, value) {
        stop("@n_elements is a read only property")
      },
      getter = function(self) {
        length(self@elements)
      }
    ),
    summary = S7::new_property(
      class = S7::class_data.frame,
      setter = function(self, value) { stop("@summary is a read only property") },
      getter =  function(self){ as.data.frame(self) }
    ),
    mirror_planes = S7::new_property(
      class = S7::S7_object,
      setter = function(self, value) {
        stop("@mirror_planes is a read only property")
      },
      getter = function(self) {
          Filter(function(el) { el@type == "Mirror Plane" }, self@elements)
      }
    ),
    proper_rotation_axes = S7::new_property(
      class = S7::class_list,
      setter = function(self, value) {
        stop("@proper_rotation_axes is a read only property")
      },
      getter = function(self) {
        Filter(function(el) { el@type == "Proper Rotation Axis" }, self@elements)
      }
    ),
    improper_rotation_axes = S7::new_property(
      class = S7::class_list,
      setter = function(self, value) {
        stop("@improper_rotation_axes is a read only property")
      },
      getter = function(self) {
        Filter(function(el) { el@type == "Improper Rotation Axis" }, self@elements)
      }
    ),
    centres_of_inversion = S7::new_property(
      class = S7::class_list,
      setter = function(self, value) {
        stop("@centres_of_inversion is a read only property")
      },
      getter = function(self) {
        Filter(function(el) { el@type == "Centre of Inversion" }, self@elements)
      }
    ),
    collections_split = S7::new_property(
      class = S7::class_list,
      getter = function(self) { split_collections_by_type(self) },
      setter = function(self, value) { stop("@collection_split is a read only property") }
    ),
    unique_proper_axis_orders = S7::new_property(
      class = S7::class_integer,
      setter = function(self, value) {
        stop("@unique_proper_axis_orders is a read only property")
      },
      getter = function(self) {
        orders <- vapply(self@proper_rotation_axes, function(element) {
          element@n
        }, FUN.VALUE = numeric(1))
        orders <- unique(orders)
        orders <- sort(orders, decreasing = FALSE)
        return(orders)
      }
    )
  ),
  validator = function(self) {
    ids <- self@ids
    elements <- self@elements

    if (length(ids) != length(elements)) {
      return(sprintf("@ids must have one entry for each of the @elements. Found [%s] IDs and [%s] elements", length(ids), length(elements)))
    }

    if (any(duplicated(ids))) {
      return(sprintf("Duplicate @ids are not allowed. Found [%s] duplicates", sum(duplicated(ids))))
    }

    if (anyNA(ids)) {
      return(sprintf("Missing (NA) @ids are not allowed. Found [%s].", sum(is.na(ids))))
    }

    if(is.null(elements)){
      return(sprintf("@elements of SymmetryElementColleciton must not be NULL"))
    }

    for (element in elements) {
      if (!inherits(element, "structures::SymmetryElement")) {
        return(sprintf("All @elements must inherit from: structures::SymmetryElement"))
      }
    }

    return(NULL)
  },
  constructor = function(elements = list(), ids = NULL) {
    if (is.null(ids)) {
      ids <- seq_along(elements)
    }
    S7::new_object(
      S7::S7_object(),
      elements = elements,
      ids = ids
    )
  }
)


## Generics ----------------------------------------------------------------
S7::method(print, SymmetryElementCollection) <- function(x, ...) {

  collections <- x@collections_split
  collection_types <- names(collections)
  string_counts <- vapply(X = seq_along(collections), function(i) {
    coll = collections[[i]]
    name = collection_types[[i]]
     sprintf("%s: %s\n", name, coll@n_elements)
  }, character(1))


  unique_proper_axis_orders = x@unique_proper_axis_orders
  unique_proper_axis_orders <- if(length(unique_proper_axis_orders) > 0) unique_proper_axis_orders else "No proper rotation axes"

  cat(
    sep = "",
    "===================\n",
    "Symmetry Element Collection\n",
    "===================\n",
    sprintf("Number of Elements: %s\n", x@n_elements),
    string_counts,
    sprintf("\nUnique Proper Axis Orders: %s", toString(unique_proper_axis_orders)),
    # sprintf("Proper Rotation Axes: %s\n", x@proper_rotation_axes@n_elements),
    # sprintf("Improper Rotation Axes: %s\n", x@improper_rotation_axes@n_elements),
    # sprintf("Mirror Planes: %s\n", x@mirror_planes@n_elements),
    # sprintf("Centres of Inversion: %s\n", x@centres_of_inversion@n_elements),
    # sprintf("Symmetry Orders (Cn): %s\n", toString(x@unique_proper_axis_orders)),
    "\n-------------------\n",
    "See @elements paramater for a list of symmetry elements\n"
  )

  return(invisible(x))
}

#' @export
S7::method(as.data.frame, SymmetryElementCollection) <- function(x, ...) {
  ls_info <- lapply(x@elements, function(el) {
    data.frame(
      type = el@type %||% character(0),
      label = el@label %||% character(0)
    )
  })

  df_data <- as.data.frame(do.call(rbind, ls_info))
  df_data$ids <- x@ids
  df_data <- df_data[c("ids", setdiff(colnames(df_data), "ids"))]

  return(df_data)
}

## Generic ------------------------------------------------------------------







## Non-Generic Methods -----------------------------------------------------

#' Filter a SymmetryElementCollection by symmetry element type
#'
#' Returns a new [`SymmetryElementCollection`] containing only those
#' symmetry elements whose `@type` matches one of the requested types.
#'
#' @param collection A [`structures::SymmetryElementCollection`] object.
#' @param types Character vector of symmetry element types to keep. See [valid_symmetry_element_types()] for valid values.
#'
#' @return A [`structures::SymmetryElementCollection`] containing the filtered
#'   subset of elements (with `ids` preserved for the retained elements).
#'
#' @examples
#' mp <- MirrorPlane(normal = c(0, 0, 1), position = c(0, 0, 0), label = "<U+03C3>_xy")
#' ci <- CentreOfInversion(position = c(0, 0, 0), label = "i")
#'
#' coll <- SymmetryElementCollection(
#'   elements = list(mp, ci),
#'   ids = c(1, 2)
#' )
#'
#' # Keep only mirror planes
#' filter_symmetry_element_collection_by_type(coll, "Mirror Plane")
#'
#' @export
filter_symmetry_element_collection_by_type <- function(collection, types = valid_symmetry_element_types()) {
  if (!all(types %in% valid_symmetry_element_types())) {
    stop(sprintf("Tried to filter symmetry collection by an invalid type. Types must only include the values [%s]. Unexpected values: [%s]", toString(valid_symmetry_element_types()), toString(setdiff(types, valid_symmetry_element_types()))))
  }

  keep <- vapply(
    X = collection@elements,
    function(element) {
      element@type %in% types
    },
    FUN.VALUE = logical(1)
  )
  S7::set_props(collection, elements = collection@elements[keep], ids = collection@ids[keep])
}

# Split a collection to a list of collections named by type
split_collections_by_type <- function(collection){

  ls_collections <- lapply(valid_symmetry_element_types(), function(type) {
    filter_symmetry_element_collection_by_type(collection, types = type)
  })

  names(ls_collections) <- valid_symmetry_element_types()
  return(ls_collections)
}

#' Add a symmetry element to a collection
#'
#' Appends a single [`SymmetryElement`] to an existing
#' [`SymmetryElementCollection`], assigning it the next available numeric ID.
#'
#' @param collection A [`structures::SymmetryElementCollection`] object.
#' @param new A [`structures::SymmetryElement`] to append to `collection`.
#'
#' @return A [`structures::SymmetryElementCollection`] with `new` appended and
#'   `ids` updated to include the newly assigned ID.
#'
#' @examples
#' mp <- MirrorPlane(normal = c(0, 0, 1), position = c(0, 0, 0), label = "<U+03C3>_xy")
#' coll <- SymmetryElementCollection()
#' coll2 <- add_symmetry_element_to_collection(coll, mp)
#' coll2@ids # contains 1
#'
#' @export
add_symmetry_element_to_collection <- function(collection, new) {
  # Check input types are appropriate
  assertions::assert_class(collection, class = "structures::SymmetryElementCollection")
  assertions::assert_class(new, class = "structures::SymmetryElement")

  # Figure out what the new ID should be
  new_id <- max(collection@ids, 0) + 1

  S7::set_props(
    collection,
    elements = c(collection@elements, new),
    ids = c(collection@ids, new_id)
  )
}

#' Combine two SymmetryElementCollections
#'
#' Merges the symmetry elements and IDs from two
#' [`SymmetryElementCollection`] objects into a single collection. The IDs
#' from `collection2` are offset so that all IDs in the combined collection
#' remain unique.
#'
#' @param collection1 A [`structures::SymmetryElementCollection`] object
#'   providing the base collection.
#' @param collection2 A [`structures::SymmetryElementCollection`] object whose
#'   elements will be appended to `collection1`.
#'
#' @return A [`structures::SymmetryElementCollection`] containing all elements
#'   from both inputs, with updated `ids`.
#'
#' @examples
#' mp <- MirrorPlane(normal = c(0, 0, 1), position = c(0, 0, 0), label = "<U+03C3>_xy")
#' ci <- CentreOfInversion(position = c(0, 0, 0), label = "i")
#'
#' coll1 <- SymmetryElementCollection(elements = list(mp), ids = 1)
#' coll2 <- SymmetryElementCollection(elements = list(ci), ids = 1)
#'
#' coll12 <- combine_symmetry_element_collections(coll1, coll2)
#' coll12@ids # e.g. c(1, 2)
#'
#' @export
combine_symmetry_element_collections <- function(collection1, collection2) {
  # Assertions
  assertions::assert_class(collection1, class = "SymmetryElementCollection")
  assertions::assert_class(collection2, class = "SymmetryElementCollection")

  # Ids
  ids_1 <- collection1@ids
  ids_2 <- collection2@ids

  # Elements
  elements1 <- collection1@elements
  elements2 <- collection2@elements

  # New ids for collection 2
  max_id <- max(ids_1, 0)
  ids_2_new <- ids_2 + max_id

  # New ids for collection 2
  S7::set_props(
    collection1,
    elements = c(elements1, elements2),
    ids = c(ids_1, ids_2_new)
  )
}

fetch_symmetry_element_from_collection <- function(collection, id, error_if_missing = TRUE){
  assertions::assert_class(collection, "structures::SymmetryElementCollection")
  assertions::assert_number(id)
  if(error_if_missing & !id %in% collection@ids) { stop("Failed to fetch element with [id=",id, "]. ID not found in collection. If you would prefer to return NULL when id is missing, set error_if_missing=FALSE")}
  idx = match(id, collection@ids)
  collection@elements[[idx]]
}


#' Fetch all proper rotation axes of a given order (Cn) from a collection
#'
#' @description
#' Returns all [`structures::ProperRotationAxis`] objects in a
#' [`structures::SymmetryElementCollection`] whose fold \code{n} matches
#' the requested \code{Cn}.
#'
#' This is the low-level worker; see
#' [fetch_all_proper_rotation_axes_with_order()] for the
#' `Molecule3D`-convenience wrapper.
#'
#' @param collection A [`structures::SymmetryElementCollection`] object.
#' @param Cn Integer (or integer-like numeric) fold/order to match
#'   (e.g. \code{2}, \code{3}, \code{6}).
#'
#' @return
#' \itemize{
#'   \item \code{NULL} if the collection contains no symmetry elements at all.
#'   \item Otherwise, a \emph{list} of [`structures::ProperRotationAxis`]
#'         objects whose \code{n} equals \code{Cn}. The list may be empty if
#'         no axes have the requested order.
#' }
#'
#' @examples
#' mp <- MirrorPlane(normal = c(0, 0, 1), position = c(0, 0, 0), label = "σ_xy")
#' ax2 <- ProperRotationAxis(n = 2L, posA = c(0, 0, 0), posB = c(0, 0, 1), label = "C2(z)")
#' ax3 <- ProperRotationAxis(n = 3L, posA = c(1, 0, 0), posB = c(1, 0, 1), label = "C3(z)")
#'
#' coll <- SymmetryElementCollection(
#'   elements = list(mp, ax2, ax3),
#'   ids      = c(1, 2, 3)
#' )
#'
#' fetch_all_proper_rotation_axes_with_order_from_collection(coll, Cn = 2L) # list of C2 axes
#' fetch_all_proper_rotation_axes_with_order_from_collection(coll, Cn = 6L) # empty list
#'
#' @export
fetch_all_proper_rotation_axes_with_order_from_collection <- function(collection, Cn) {
  assertions::assert_class(collection, "structures::SymmetryElementCollection")

  # No symmetry elements at all → NULL, to match older behaviour
  if (collection@n_elements == 0L) {
    return(NULL)
  }

  # Normalise order to integer-like
  if (!is_integerlike(Cn)) {
    stop("Cn must be an integer-like value, not: ", toString(Cn))
  }
  Cn <- as.integer(Cn)

  # Restrict to proper rotation axes only
  axes <- collection@proper_rotation_axes
  if (length(axes) == 0L) {
    # Collection has elements but none of type "Proper Rotation Axis"
    return(list())
  }

  keep <- vapply(
    X = axes,
    FUN = function(axis) axis@n == Cn,
    FUN.VALUE = logical(1)
  )

  axes[keep]
}






# Create an Class for each symmetry element that inherits from 'SymmetryElement'


# [Symmetry Element] ProperRotationAxis  -------------------------------------

## Class  -------------------------------------
#' ProperRotationAxis: rotational symmetry axis (C_n)
#'
#' Represents a proper rotation axis (C\eqn{_n}) in 3D. The axis is defined by a
#' point `posA` and an orientation `direction` (unit vector). A second point
#' `posB` lies along the axis at distance `L` from `posA` and is derived as
#' `posA + direction * L`.
#'
#' @details
#' * `n` is the order (fold) of the axis, i.e. rotation by `360 / n` degrees.
#' * You can construct an axis either by:
#'   - supplying `posA` and `posB` (direction and `L` are inferred), or
#'   - supplying `posA` and `direction` (direction is normalised; `L` defaults to 1
#'     if not given).
#' * `posB` is read-only and always computed from `posA`, `direction`, and `L`.
#'
#' @param n Numeric scalar. Axis fold/order (whole number >= 1).
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
#'   n = 3L,
#'   posA = c(0, 0, 0),
#'   posB = c(0, 0, 1),
#'   label = "C3(z)"
#' )
#'
#' # Define via direction (C2 axis along +Z; L defaults to 1)
#' ax2 <- ProperRotationAxis(
#'   n = 2L,
#'   posA = c(0, 0, 0),
#'   direction = c(0, 0, 5),
#'   label = "C2(z)"
#' )
#'
#' ax1@type # "Proper Rotation Axis"
#'
#' @export
ProperRotationAxis <- S7::new_class(
  name = "ProperRotationAxis",
  parent = SymmetryElement,
  properties = list(
    n = S7::new_property(
      class = S7::class_numeric,
      setter = function(self, value) {
        if (is_integerlike(value)) {
          value <- as.integer(value)
        }
        S7::prop(self, "n") <- value
        return(self)
      },
      validator = function(value) {
        if (length(value) != 1L) {
          return(
            sprintf(
              "ProperRotationAxis fold symmetry (@n) must be a single number (length 1 integer), not a vector of length [%s]",
              length(value)
            )
          )
        }
        if (is.na(value) || value < 1L) {
          return(sprintf("ProperRotationAxis fold symmetry (@n) must be >= 1. Not [%s]", value))
        }

        if (is.infinite(value)) {
          return(sprintf("ProperRotationAxis fold symmetry (@n) can not be an infinite value"))
        }

        if (!is_integerlike(value)) {
          return(sprintf("ProperRotationAxis fold symmetry (@n) must be a whole number, not [%g]", value))
        }
      },
    ),
    posA = S7::new_property(
      class = S7::class_numeric,
      validator = function(value) {
        if (length(value) != 3L) {
          return(
            sprintf(
              "posA must be a numeric vector of length 3: c(x, y, z). Supplied length [%s]",
              length(value)
            )
          )
        }
        if (any(!is.finite(value)) || any(is.na(value))) {
          return("posA must contain only finite, non-NA values.")
        }
      },
      setter = function(self, value) {
        names(value) <- c("x", "y", "z")[seq_along(value)]
        S7::prop(self, "posA") <- value
        return(self)
      }
    ),

    # Never set directly
    posB = S7::new_property(
      class = S7::class_numeric,
      setter = function(self, value) {
        stop("@posB is read-only; set @direction or @L instead")
      },
      getter = function(self) {
        out <- self@posA + self@direction * self@L
        names(out) <- c("x", "y", "z")
        return(out)
      },
    ),



    # A unit direction vector representing the axis (computed as direction from posA to posB)
    direction = S7::new_property(
      class = S7::class_numeric,
      validator = function(value) {
        if (length(value) != 3L) {
          return(
            sprintf(
              "direction must be a numeric vector of length 3: c(x, y, z). Supplied length [%s]",
              length(value)
            )
          )
        }
        if (any(!is.finite(value)) || any(is.na(value))) {
          return("direction must contain only finite, non-NA values.")
        }
      },
      getter = function(self) {
        names(self@direction) <- c("x", "y", "z")
        return(self@direction)
      },
      setter = function(self, value) {
        value <- normalise(value)
        self@direction <- value
        return(self)
      }
    ),

    # Length between posA and posB
    L = S7::new_property(
      class = S7::class_numeric,
      validator = function(value) {
        if (length(value) != 1L) {
          return("L must be a length-1 numeric")
        }
        if (!is.finite(value) || is.na(value)) {
          return("L must be finite")
        }
        if (value <= 0) {
          return("L must be > 0 to define an axis")
        }
      }
    )
  ),

  # Validate inter-property relationships
  validator = function(self) {
    # With derivation, this ensures L > 0 and endpoints differ
    if (identical(self@posA, self@posB)) {
      return("posA and posB must be distinct (L must be > 0)")
    }
    # browser()
    NULL
  },
  constructor = function(n, posA, posB = NULL, direction = NULL, L = NULL, label = "unnamed") {
    if (is.null(posB) && is.null(direction)) stop("Supply either `posB` or `direction`")
    if (!is.null(posB) && !is.null(direction)) stop("Supply only one of `posB` or `direction`")
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
      # direction provided
      direction <- normalise(direction)
      L <- if (is.null(L)) 1 else L
    }
    names(direction) <- c("x", "y", "z")

    S7::new_object(
      S7::S7_object(),
      n = n,
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
    "Proper Rotation Axis\n",
    "===================\n",
    sprintf("Fold Symmetry / Order (Cn): C%g\n", x@n),
    sprintf("PosA: %s\n", toString(x@posA)),
    sprintf("Direction: %s\n", toString(x@direction)),
    sprintf("PosB: %s\n", toString(x@posB)),
    sprintf("Label: %s\n", toString(x@label)),
    sep = ""
  )
  return(invisible(x))
}


#' @export
S7::method(as.data.frame, ProperRotationAxis) <- function(x, ...) {
  posA <- x@posA
  posB <- x@posB
  label <- x@label
  Cn <- x@n

  data.frame(
    label = label,
    Cn = n,
    x = posA["x"],
    y = posA["y"],
    z = posA["z"],
    xend = posB["x"],
    yend = posB["y"],
    zend = posB["z"],
    row.names = NULL
  )
}


#' @export
S7::method(transform_symmetry_element, ProperRotationAxis) <- function(x, transformation, ...) {
  transform_symmetry_axis(x, transformation = transformation, ...)
}





## Helpers -----------------------------------------------------------------
is_integerlike <- function(x, tol = .Machine$double.eps^0.5) {
  if (length(x) == 0L) {
    return(FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    return(FALSE)
  }
  all(abs(x - round(x)) <= tol)
}



## Non-Generic Methods -----------------------------------------------------

#' Apply a 3D transformation to a rotation axis
#'
#' @description
#' Applies a user-supplied transformation function to the geometric definition
#' of a [`ProperRotationAxis`] or [`ImproperRotationAxis`]. For a proper axis,
#' the endpoints (`posA`, `posB`) are transformed and the direction/length
#' recomputed. For an improper axis, the same transformation is also applied
#' to `plane_point` so that the associated mirror plane remains consistent
#' with the transformed axis.
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
#' ax <- ProperRotationAxis(n = 3L, posA = c(0, 0, 0), posB = c(0, 0, 1))
#'
#' # 1) Translate by (+1, +2, +3)
#' translate <- function(p, dx = 0, dy = 0, dz = 0) {
#'   c(x = p["x"] + dx, y = p["y"] + dy, z = p["z"] + dz)
#' }
#' ax_t <- transform_symmetry_axis(ax, translate, dx = 1, dy = 2, dz = 3)
#'
#' # 2) Rotate around Z by 90 degrees about the origin
#' rotate_z <- function(p, theta) {
#'   c(
#'     x = cos(theta) * p["x"] - sin(theta) * p["y"],
#'     y = sin(theta) * p["x"] + cos(theta) * p["y"],
#'     z = p["z"]
#'   )
#' }
#' ax_r <- transform_symmetry_axis(ax, rotate_z, theta = pi / 2)
#'
#' # 3) A transform that returns a list is also accepted
#' as_list <- function(p) list(x = p["x"] + 1, y = p["y"], z = p["z"])
#' ax_l <- transform_symmetry_axis(ax, as_list)
#'
#' @export
transform_symmetry_axis <- function(x, transformation, ...) {
  assertions::assert_class(
    x,
    c("structures::ProperRotationAxis", "structures::ImproperRotationAxis")
  )
  assertions::assert_function(transformation)

  # --- Transform the axis line (posA / posB) ---
  posA <- x@posA
  posB <- x@posB

  posAnew <- transformation(posA, ...)
  posBnew <- transformation(posB, ...)

  if (is.list(posAnew)) posAnew <- unlist(posAnew)
  if (is.list(posBnew)) posBnew <- unlist(posBnew)

  direction <- create_vector_from_start_end(
    start = posAnew,
    end   = posBnew,
    unit  = TRUE
  )
  L <- compute_distance(posAnew, posBnew)

  # --- Base update (works for ProperRotationAxis) ---
  new_axis <- S7::set_props(
    x,
    posA     = posAnew,
    direction = direction,
    L        = L
  )

  # --- Extra step for ImproperRotationAxis: transform plane_point too ---
  if (is_improper_rotation_axis(new_axis)) {
    plane_point <- new_axis@plane_point
    plane_point_new <- transformation(plane_point, ...)

    if (is.list(plane_point_new)) plane_point_new <- unlist(plane_point_new)

    # Optional: you could assert that plane_point_new lies on the new axis,
    # but that may be too strict for arbitrary transformations.
    new_axis <- S7::set_props(new_axis, plane_point = plane_point_new)
  }

  return(new_axis)
}


# [Symmetry Element] ImproperRotationAxis ----------------------------------------------------

## Class -------------------------------------------------------------

#' ImproperRotationAxis: improper rotational symmetry axis (S_n)
#'
#' Represents an improper rotation axis (S\eqn{_n}) in 3D, combining a proper
#' rotation (C\eqn{_n}) about the axis with reflection in a plane perpendicular
#' to that axis. Geometrically, the axis is described exactly as in
#' [`ProperRotationAxis()`] (via `n`, `posA`, `direction`, `L`), with an
#' additional point `plane_point` that fixes the location of the mirror plane.
#'
#' The mirror plane is defined as the set of points `x` satisfying
#' `(x - plane_point) <U+00B7> direction = 0`, i.e. a plane with normal equal to
#' the axis direction, passing through `plane_point`. In typical molecular
#' symmetry, `plane_point` will also lie on the axis (e.g. at `posA`).
#'
#' @param n Numeric scalar. Axis fold/order (whole number >= 1).
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
#'   n = 4L,
#'   posA = c(0, 0, 0),
#'   posB = c(0, 0, 1),
#'   plane_point = c(0, 0, 0),
#'   label = "S4(z)"
#' )
#'
#' iax1@type # "Improper Rotation Axis"
#'
#' # S6 axis along +Z, with the mirror plane shifted along the axis
#' iax2 <- ImproperRotationAxis(
#'   n = 6L,
#'   posA = c(0, 0, 0),
#'   direction = c(0, 0, 1),
#'   L = 2,
#'   plane_point = c(0, 0, 0.5),
#'   label = "S6(z)"
#' )
#'
#' @export
ImproperRotationAxis <- S7::new_class(
  parent = ProperRotationAxis,
  name = "ImproperRotationAxis",
  properties = list(
    plane_point = S7::new_property(
      class = S7::class_numeric,
      validator = function(value) {
        if (length(value) != 3L) {
          return(
            sprintf(
              "ImproperRotationAxis @plane_point must be a numeric vector of length 3: not [%s]",
              length(value)
            )
          )
        }
        if (any(!is.finite(value)) || any(is.na(value))) {
          return("ImproperRotationAxis @plane_point must contain only finite, non-NA values.")
        }
        NULL
      },
      setter = function(self, value) {
        names(value) <- c("x", "y", "z")[seq_along(value)]
        S7::prop(self, "plane_point") <- value
        self
      }
    )
  ),
  constructor = function(n,
                         posA,
                         plane_point = posA,
                         posB = NULL,
                         direction = NULL,
                         L = NULL,
                         label = "unnamed") {
    if (is.null(posB) && is.null(direction)) {
      stop("Supply either `posB` or `direction`")
    }
    if (!is.null(posB) && !is.null(direction)) {
      stop("Supply only one of `posB` or `direction`")
    }

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
    assertions::assert_length(plane_point,
      length = 3,
      msg = "plane_point must be a numeric vector of length 3: c(x, y, z)"
    )
    assertions::assert_no_missing(plane_point)
    assertions::assert(all(is.finite(plane_point)),
      msg = "plane_point must have no infinite values"
    )
    names(plane_point) <- c("x", "y", "z")[seq_along(plane_point)]

    # Enforce that plane_point lies on the axis
    v <- plane_point - posA
    if (magnitude(v) > 0) {
      if (max(abs(normalise(v) - direction)) > 1e-6) {
        stop("plane_point must lie on the rotation axis (collinear with direction from posA).")
      }
    }

    S7::new_object(
      S7::S7_object(),
      n = n,
      posA = posA,
      direction = direction,
      L = L,
      plane_point = plane_point,
      label = label
    )
  }
)


## Generics ----------------------------------------------------------------

#' @export
S7::method(print, ImproperRotationAxis) <- function(x, ...) {
  cat(
    "===================\n",
    "Improper Rotation Axis\n",
    "===================\n",
    sprintf("Fold Symmetry / Order (Sn): S%g\n", x@n),
    sprintf("Normal (direction): %s\n", toString(x@direction)),
    sprintf("Normal (position): %s\n", toString(x@posA)),
    sprintf("Plane (position): %s\n", toString(x@plane_point)),
    sprintf("Label: %s\n", toString(x@label)),
    sep = ""
  )
  return(invisible(x))
}

#' @export
S7::method(transform_symmetry_element, ImproperRotationAxis) <- function(x, transformation, ...) {
  transform_symmetry_axis(x, transformation = transformation, ...)
}

# [Symmetry Element] MirrorPlane -------------------------------------------------------------

## Class -------------------------------------------------------------

#' MirrorPlane: reflection symmetry element (<U+03C3>)
#'
#' Represents a mirror plane in 3D space, defined by a unit normal vector and
#' a point lying on the plane. The plane contains all points `p` such that
#' `(p - position) <U+00B7> normal = 0`.
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
#'   label = "<U+03C3>_xy"
#' )
#'
#' mp@type # "Mirror Plane"
#'
#' @export
MirrorPlane <- S7::new_class(
  parent = SymmetryElement,
  name = "MirrorPlane",
  properties = list(
    normal = S7::new_property(
      class = S7::class_numeric,
      validator = function(value) {
        tol <- 1e-18
        if (magnitude(value) <= tol) {
          return(sprintf("Mirror plane @normal must be a vector with magnitude > 0, not [%d]", magnitude(value)))
        }
        if (length(value) != 3) {
          return(sprintf("Mirror plane @normal must be a numeric vector with length 3: Not [%s]", length(value)))
        }
        if (length(value) != 3) {
          return(sprintf("Mirror plane @normal must be a numeric vector with length 3: Not [%s]", length(value)))
        }
        return(NULL)
      }
    ),
    position = S7::new_property(
      class = S7::class_numeric,
      validator = function(value) {
        if (length(value) != 3) {
          return(sprintf("Mirror plane @position must be a numeric vector with length 3: Not [%s]", length(value)))
        }
        return(NULL)
      }
    )
  ),
  constructor = function(normal, position, label = "unnamed") {
    S7::new_object(
      S7::S7_object(),
      normal = normal,
      position = position,
      label = label
    )
  }
)


## Generics ----------------------------------------------------------------
#' @export
S7::method(print, MirrorPlane) <- function(x, ...) {
  cat(
    "===================\n",
    "Mirror Plane\n",
    "===================\n",
    sprintf("Normal: %s\n", toString(x@normal)),
    sprintf("Position: %s\n", toString(x@position)),
    sprintf("Label: %s\n", toString(x@label)),
    sep = ""
  )
  return(invisible(x))
}

#' @export
S7::method(transform_symmetry_element, MirrorPlane) <- function(x, transformation, ...) {
  assertions::assert_function(transformation)

  # A point on the plane and a second point offset along the normal
  pos        <- x@position
  pos_normal <- pos + x@normal

  pos_new        <- transformation(pos,        ...)
  pos_normal_new <- transformation(pos_normal, ...)

  # Allow transforms that return lists
  if (is.list(pos_new))        pos_new        <- unlist(pos_new)
  if (is.list(pos_normal_new)) pos_normal_new <- unlist(pos_normal_new)

  # Recompute unit normal from the transformed points
  new_normal <- create_vector_from_start_end(
    start = pos_new,
    end   = pos_normal_new,
    unit  = TRUE
  )

  x <- S7::set_props(
    x,
    position = pos_new,
    normal   = new_normal
  )
  x
}



# [Symmetry Element] CentreOfInversion -------------------------------------------------------

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
#' ci@type
#'
#' print(ci)
#'
#' @export
CentreOfInversion <- S7::new_class(
  parent = SymmetryElement,
  name = "CentreOfInversion",
  properties = list(
    position = S7::new_property(
      class = S7::class_numeric,
      validator = function(value) {
        if (length(value) != 3) {
          return(sprintf("Center of Inversion @position must be a numeric vector with length 3: Not [%s]", length(value)))
        }
        return(NULL)
      }
    )
  ),
  constructor = function(position, label = "unnamed") {
    S7::new_object(
      S7::S7_object(),
      position = position,
      label = label
    )
  }
)


# Generics ----------------------------------------------------------------
#' @export
S7::method(print, CentreOfInversion) <- function(x, ...) {
  cat(
    "===================\n",
    "Centre of Inversion\n",
    "===================\n",
    sprintf("Position: %s\n", toString(x@position)),
    sprintf("Label: %s\n", toString(x@label)),
    sep = ""
  )
  return(invisible(x))
}

#' @export
S7::method(transform_symmetry_element, CentreOfInversion) <- function(x, transformation, ...) {
  assertions::assert_function(transformation)

  pos_new <- transformation(x@position, ...)
  if (is.list(pos_new)) pos_new <- unlist(pos_new)

  x <- S7::set_props(x, position = pos_new)
  x
}


# Helpers -----------------------------------------------------------------
## Class Checkers ----------------------------------------------------------

#' Check if an object is a MirrorPlane
#'
#' @param x An object to test.
#'
#' @return Logical. `TRUE` if `x` inherits from [`structures::MirrorPlane`], otherwise `FALSE`.
#'
#' @examples
#' mp <- MirrorPlane(normal = c(0, 0, 1), position = c(0, 0, 0))
#' is_mirror_plane(mp) # TRUE
#' is_mirror_plane("nope") # FALSE
#'
#' @export
is_mirror_plane <- function(x) {
  inherits(x, "structures::MirrorPlane")
}

#' Check if an object is a CentreOfInversion
#'
#' @param x An object to test.
#'
#' @return Logical. `TRUE` if `x` inherits from [`structures::CentreOfInversion`], otherwise `FALSE`.
#'
#' @examples
#' ci <- CentreOfInversion(position = c(0, 0, 0))
#' is_centre_of_inversion(ci) # TRUE
#' is_centre_of_inversion(123) # FALSE
#'
#' @export
is_centre_of_inversion <- function(x) {
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
#' is_improper_rotation_axis(iax) # TRUE
#' is_improper_rotation_axis(NULL) # FALSE
#'
#' @export
is_improper_rotation_axis <- function(x) {
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
#' ax <- ProperRotationAxis(n = 3L, posA = c(0, 0, 0), posB = c(0, 0, 1))
#' is_proper_rotation_axis(ax) # TRUE
#' is_proper_rotation_axis("nope") # FALSE
#'
#' @export
is_proper_rotation_axis <- function(x) {
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
#' ax <- ProperRotationAxis(n = 2L, posA = c(0, 0, 0), posB = c(0, 0, 1))
#' is_symmetry_axis(ax)
#' is_symmetry_axis("not_an_axis")
#'
#' @export
is_symmetry_axis <- function(x) {
  inherits(x, "structures::ProperRotationAxis") | inherits(x, "structures::ImproperRotationAxis")
}

is_symmetry_element <- function(x) {
  inherits(x, "structures::SymmetryElement")
}

identify_symmetry_type <- function(x, strict = TRUE) {
  if (is_mirror_plane(x)) {
    return("Mirror Plane")
  } else if (is_centre_of_inversion(x)) {
    return("Centre of Inversion")
  } else if (is_improper_rotation_axis(x)) {
    return("Improper Rotation Axis")
  } else if (is_proper_rotation_axis(x)) {
    return("Proper Rotation Axis")
  }

  # else if(is_symmetry_element(x))
  #   return("Other Symmetry element")

  warn_message <- sprintf("Object is not a known symmetry element. Class is: [%s]", toString(class(x)))
  if (strict) {
    stop(warn_message)
  } else {
    warning(warn_message)
  }

  return("Unknown")
}

## Miscellaneous Helpers -----------------------------------------------------------------
symmetry_element_abbreviations <- function() {
  c(
    "Cn" = "Proper Rotation Axes",
    "Sn" = "Improper Rotation Axes",
    "sigma" = "Mirror Plane",
    "i" = "Centre of Inversion",
  )
}

#' Valid symmetry element types
#'
#' Returns the set of recognised symmetry element type strings used by
#' [`SymmetryElement`] and its subclasses.
#'
#' @return A character vector of valid symmetry element types.
#'
#' @examples
#' valid_symmetry_element_types()
#'
#' @export
valid_symmetry_element_types <- function() {
  c("Proper Rotation Axis",  "Improper Rotation Axis", "Mirror Plane", "Centre of Inversion")
}



