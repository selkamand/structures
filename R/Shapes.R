# Shape -------------------------------------------------------------------


#' Shape: a simple 3D polyhedral container
#'
#' @description
#' `Shape` is a lightweight S7 class representing a polyhedral shape via a set of
#' vertices, edges, and faces. It is used by helpers such as `shape_tetrahedron()`,
#' `shape_cube()`, `shape_octahedron()`, `shape_icosahedron()`, and
#' `shape_dodecahedron()` to construct common Platonic solids. The class enforces
#' basic validation for internal consistency (e.g. edges and faces must reference
#' known vertex names).
#'
#' @param vertices A data.frame with columns `name` (character, unique),
#'   `x`, `y`, and `z` (numeric) giving vertex coordinates.
#' @param edges A data.frame with columns `vertex1`, `vertex2` (character)
#'   listing undirected edges by vertex name. May be empty.
#' @param faces A named list of character vectors, where each vector contains
#'   vertex names that form a polygonal face. May be empty.
#' @param symmetry_axes An optional list of [`SymAxis`] objects describing known
#'   symmetry axes. May be empty.
#'
#' @return A `Shape` S7 object.
#'
#' @examples
#' # Create a custom shape with two vertices and one edge
#' v <- data.frame(
#'   name = c("A", "B"),
#'   x = c(0, 1),
#'   y = c(0, 0),
#'   z = c(0, 0),
#'   stringsAsFactors = FALSE
#' )
#' e <- data.frame(vertex1 = "A", vertex2 = "B", stringsAsFactors = FALSE)
#' s <- Shape(vertices = v, edges = e, faces = list())
#' s
#'
#' @export
Shape <- S7::new_class(
  name = "Shape",
  properties = list(
    vertices = S7::new_property(
      class = S7::class_data.frame,
      validator = function(value) {
        required_cols <- c("name", "x", "y", "z")
        missing <- setdiff(required_cols, names(value))
        if (length(missing) > 0) {
          return(sprintf("vertices is missing required columns: %s", toString(missing)))
        }
        if (!is.character(value$name)) {
          return("vertices$name must be character")
        }
        if (anyNA(value$name)) {
          return("vertices$name must not contain NA")
        }
        if (anyDuplicated(value$name)) {
          return("vertices$name must be unique")
        }
        if (!is.numeric(value$x) || !is.numeric(value$y) || !is.numeric(value$z)) {
          return("vertices x,y,z must be numeric")
        }
      }
    ),
    edges = S7::new_property(
      class = S7::class_data.frame,
      validator = function(value) {
        if (nrow(value) == 0) {
          return(NULL)
        }
        required_cols <- c("vertex1", "vertex2")
        missing <- setdiff(required_cols, names(value))
        if (length(missing) > 0) {
          return(sprintf("edges is missing required columns: %s", toString(missing)))
        }
        if (!is.character(value$vertex1) || !is.character(value$vertex2)) {
          return("edges$vertex1 and edges$vertex2 must be character")
        }
      }
    ),
    faces = S7::new_property(
      class = S7::class_list,
      validator = function(value) {
        if (length(value) == 0) {
          return(NULL)
        }
        are_char_vecs <- vapply(value, function(x) is.character(x), logical(1))
        if (!all(are_char_vecs)) {
          return("faces must be a list of character vectors (vertex names)")
        }
      }
    ),
    symmetry_axes = S7::new_property(
      class = S7::class_list,
      validator = function(value) {
        if (length(value) == 0) {
          return(NULL)
        }
        ok <- vapply(value, is_symmetry_axis, logical(1))
        if (!all(ok)) {
          return("symmetry_axes must be a list of SymAxis objects")
        }
      }
    ),

    # Computed Properties (edge centers, face centers, triangles, quads, etcs)
    edge_positions = S7::new_property(
      setter = function(self, value) {
        stop("@edge_positions is a read only property")
      },
      getter = function(self) {
        vertex1 <- self@edges[[1]]
        vertex2 <- self@edges[[2]]
        df_vertices <- self@vertices
        idx_v1 <- match(vertex1, df_vertices$name)
        idx_v2 <- match(vertex2, df_vertices$name)
        x <- df_vertices[["x"]][idx_v1]
        y <- df_vertices[["y"]][idx_v1]
        z <- df_vertices[["z"]][idx_v1]
        xend <- df_vertices[["x"]][idx_v2]
        yend <- df_vertices[["y"]][idx_v2]
        zend <- df_vertices[["z"]][idx_v2]

        data.frame(
          vertex1 = vertex1,
          vertex2 = vertex2,
          x = x,
          y = y,
          z = z,
          xend = xend,
          yend = yend,
          zend = zend,
          x_middle = pmean(x, xend, na.rm = TRUE),
          y_middle = pmean(y, yend, na.rm = TRUE),
          z_middle = pmean(z, zend, na.rm = TRUE)
        )
      }
    )
  ),
  constructor = function(vertices = default_vertices(), edges = default_edges(), faces = list(), symmetry_axes = list()) {
    S7::new_object(
      S7::S7_object(),
      vertices = vertices,
      edges = edges,
      faces = faces,
      symmetry_axes = symmetry_axes
    )
  },
  validator = function(self) {
    vertex_names <- self@vertices$name
    # edges reference check
    if (nrow(self@edges) > 0) {
      v1_ok <- self@edges$vertex1 %in% vertex_names
      v2_ok <- self@edges$vertex2 %in% vertex_names
      if (!all(v1_ok & v2_ok)) {
        return("edges reference vertex names not present in vertices$name")
      }
    }
    # faces reference check
    if (length(self@faces) > 0) {
      all_faces <- unique(unlist(self@faces, use.names = FALSE))
      if (length(all_faces) > 0 && !all(all_faces %in% vertex_names)) {
        return("faces contain vertex names not present in vertices$name")
      }
    }
    NULL
  }
)

## Defaults
default_vertices <- function() {
  data.frame(
    name = character(0),
    x = numeric(0),
    y = numeric(0),
    z = numeric(0),
    stringsAsFactors = FALSE
  )
}

default_edges <- function() {
  data.frame(vertex1 = character(0), vertex2 = character(0))
}


#' Tetrahedron (regular)
#'
#' @description
#' Constructs a regular tetrahedron centered at the origin with vertices at
#' permutations of `(±1, ±1, ±1)` of odd parity. Edges and faces are populated.
#'
#' @return A [`Shape`] representing a tetrahedron.
#'
#' @examples
#' shp <- shape_tetrahedron()
#' shp
#'
#' @export
shape_tetrahedron <- function() {
  vertices <- matrix(
    c(
      1, 1, 1, # A
      -1, -1, 1, # B
      -1, 1, -1, # C
      1, -1, -1 # D
    ),
    nrow = 4L, byrow = TRUE,
    dimnames = list(c("A", "B", "C", "D"), c("x", "y", "z"))
  )
  df_vertices <- as.data.frame(vertices)
  df_vertices$name <- rownames(vertices)

  df_edges <- data.frame(
    vertex1 = c("A", "A", "A", "B", "B", "C"),
    vertex2 = c("B", "C", "D", "C", "D", "D"),
    stringsAsFactors = FALSE
  )

  faces <- list(
    BCD = c("B", "C", "D"),
    ACD = c("A", "C", "D"),
    ABD = c("A", "B", "D"),
    ABC = c("A", "B", "C")
  )

  Shape(vertices = df_vertices, edges = df_edges, faces = faces)
}


# Helpers ----------------------------------------------------------------

compute_shape_edges_by_min_distance <- function(vertices_df, tol = 1e-8) {
  coords <- as.matrix(vertices_df[, c("x", "y", "z")])
  n <- nrow(coords)
  if (n <= 1) {
    return(default_edges())
  }
  dmat <- as.matrix(dist(coords))
  min_nonzero <- min(dmat[dmat > 0])
  edges <- which(abs(dmat - min_nonzero) < tol & upper.tri(dmat), arr.ind = TRUE)
  if (nrow(edges) == 0) {
    return(default_edges())
  }
  data.frame(
    vertex1 = vertices_df$name[edges[, 1]],
    vertex2 = vertices_df$name[edges[, 2]],
    stringsAsFactors = FALSE
  )
}

# Cube -------------------------------------------------------------------

#' Cube (regular)
#'
#' @description
#' Constructs a unit cube centered at the origin with vertices at all
#' combinations of `(±1, ±1, ±1)`. Edges are determined by unit Hamming-distance.
#'
#' @return A [`Shape`] representing a cube.
#'
#' @examples
#' shp <- shape_cube()
#' shp
#'
#' @export
shape_cube <- function() {
  signs <- expand.grid(x = c(-1, 1), y = c(-1, 1), z = c(-1, 1), KEEP.ALIVE = FALSE, stringsAsFactors = FALSE)
  signs$name <- paste0("V", seq_len(nrow(signs)))
  df_vertices <- signs[, c("name", "x", "y", "z")]
  coords <- as.matrix(df_vertices[, c("x", "y", "z")])
  n <- nrow(coords)
  pairs <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  is_edge <- function(i, j) {
    sum(coords[i, ] != coords[j, ]) == 1
  }
  keep <- apply(pairs, 1, function(ij) is_edge(ij[1], ij[2]))
  pairs <- pairs[keep, , drop = FALSE]
  df_edges <- if (nrow(pairs)) {
    data.frame(
      vertex1 = df_vertices$name[pairs[, 1]],
      vertex2 = df_vertices$name[pairs[, 2]],
      stringsAsFactors = FALSE
    )
  } else {
    default_edges()
  }
  Shape(vertices = df_vertices, edges = df_edges, faces = list())
}

# Octahedron --------------------------------------------------------------

#' Octahedron (regular)
#'
#' @description
#' Constructs a regular octahedron centered at the origin with vertices on
#' the coordinate axes `(±1, 0, 0)`, `(0, ±1, 0)`, `(0, 0, ±1)`.
#'
#' @return A [`Shape`] representing an octahedron.
#'
#' @examples
#' shp <- shape_octahedron()
#' shp
#'
#' @export
shape_octahedron <- function() {
  verts <- rbind(
    c(1, 0, 0),
    c(-1, 0, 0),
    c(0, 1, 0),
    c(0, -1, 0),
    c(0, 0, 1),
    c(0, 0, -1)
  )
  rownames(verts) <- paste0("V", seq_len(nrow(verts)))
  colnames(verts) <- c("x", "y", "z")
  df_vertices <- as.data.frame(verts)
  df_vertices$name <- rownames(verts)
  df_edges <- compute_shape_edges_by_min_distance(df_vertices)
  Shape(vertices = df_vertices, edges = df_edges, faces = list())
}

# Icosahedron -------------------------------------------------------------

#' Icosahedron (regular)
#'
#' @description
#' Constructs a regular icosahedron using the golden ratio `phi = (1 + sqrt(5))/2`.
#' Vertices are placed at standard coordinates; edges are derived from minimal
#' equal edge lengths.
#'
#' @return A [`Shape`] representing an icosahedron.
#'
#' @examples
#' shp <- shape_icosahedron()
#' shp
#'
#' @export
shape_icosahedron <- function() {
  phi <- (1 + sqrt(5)) / 2
  verts <- rbind(
    c(0, 1, phi), c(0, -1, phi), c(0, 1, -phi), c(0, -1, -phi),
    c(1, phi, 0), c(-1, phi, 0), c(1, -phi, 0), c(-1, -phi, 0),
    c(phi, 0, 1), c(phi, 0, -1), c(-phi, 0, 1), c(-phi, 0, -1)
  )
  rownames(verts) <- paste0("V", seq_len(nrow(verts)))
  colnames(verts) <- c("x", "y", "z")
  df_vertices <- as.data.frame(verts)
  df_vertices$name <- rownames(verts)
  df_edges <- compute_shape_edges_by_min_distance(df_vertices, tol = 1e-8)
  Shape(vertices = df_vertices, edges = df_edges, faces = list())
}

# Dodecahedron ------------------------------------------------------------

#' Dodecahedron (regular)
#'
#' @description
#' Constructs a regular dodecahedron using golden ratio relationships, combining
#' cube corners and additional vertices. Edges are derived from minimal equal
#' edge lengths.
#'
#' @return A [`Shape`] representing a dodecahedron.
#'
#' @examples
#' shp <- shape_dodecahedron()
#' shp
#'
#' @export
shape_dodecahedron <- function() {
  phi <- (1 + sqrt(5)) / 2
  invphi <- 1 / phi
  verts <- rbind(
    c(1, 1, 1), c(1, 1, -1), c(1, -1, 1), c(1, -1, -1),
    c(-1, 1, 1), c(-1, 1, -1), c(-1, -1, 1), c(-1, -1, -1),
    c(0, invphi, phi), c(0, invphi, -phi), c(0, -invphi, phi), c(0, -invphi, -phi),
    c(invphi, phi, 0), c(-invphi, phi, 0), c(invphi, -phi, 0), c(-invphi, -phi, 0),
    c(phi, 0, invphi), c(phi, 0, -invphi), c(-phi, 0, invphi), c(-phi, 0, -invphi)
  )
  rownames(verts) <- paste0("V", seq_len(nrow(verts)))
  colnames(verts) <- c("x", "y", "z")
  df_vertices <- as.data.frame(verts)
  df_vertices$name <- rownames(verts)
  df_edges <- compute_shape_edges_by_min_distance(df_vertices, tol = 1e-8)
  Shape(vertices = df_vertices, edges = df_edges, faces = list())
}

# Generics ----------------------------------------------------------------


#' @export
S7::method(as.matrix, Shape) <- function(x, ...) {
  as.matrix(x@vertices[c("x", "y", "z")])
}
