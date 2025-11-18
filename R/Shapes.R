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
#' @param symmetry_axes An optional list of [`ProperRotationAxis`] objects describing known
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

    ## Writeable ---------------------------------------------------------------
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
    symmetry_elements = SymmetryElementCollection,

    ## Read Only ---------------------------------------------------------------
    face_centroids = S7::new_property(
      class = S7::class_data.frame,
      setter = function(self, value) {
        stop("@face_centers is a read only property")
      },
      getter = function(self) {
        faces <- self@faces
        if (length(faces) == 0) {
          return(data.frame(
            face = character(0),
            x = numeric(0),
            y = numeric(0),
            z = numeric(0),
            vertex_x = character(0),
            vertex_y = character(0),
            vertex_z = character(0),

            stringsAsFactors = FALSE
          ))
        }

        verts <- self@vertices


        centers <- lapply(names(faces), function(fname) {
          vnames <- faces[[fname]]
          idx <- match(vnames, verts$name)
          coords <- verts[idx, c("x", "y", "z"), drop = FALSE]
          data.frame(
            face = fname,
            x = mean(coords$x),
            y = mean(coords$y),
            z = mean(coords$z),
            vertices = paste0(vnames, collapse = "-")
          )
        })
        df_centers <- do.call(rbind, centers)
        rownames(df_centers) <- NULL
        return(df_centers)
      }
    ),
    geometric_center = S7::new_property(
      class = S7::class_numeric,
      setter = function(self, value) { stop("@geometric_center is a read only property") },
      getter = function(self){
        vertices = self@vertices
        x = mean(vertices$x)
        y = mean(vertices$y)
        z = mean(vertices$z)
        c("x" = x, "y" = y, "z"=z)
      }
    ),

    # symmetry_axes = S7::new_property(
    #   class = S7::class_list,
    #   validator = function(value) {
    #     if (length(value) == 0) {
    #       return(NULL)
    #     }
    #     ok <- vapply(value, is_symmetry_axis, logical(1))
    #     if (!all(ok)) {
    #       return("symmetry_axes must be a list of ProperRotationAxis objects")
    #     }
    #   }
    # ),

    # Computed Properties (edge centers, face centers, triangles, quads, etcs)
    edge_positions = S7::new_property(
      class = S7::class_data.frame,
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
    ),
    # Edges interleaved
    segments_interleaved = S7::new_property(
      class = S7::class_data.frame,
      setter = function(self, value) { stop("@edge_positions is a read only property")},
      getter = function(self){
        interleaved = to_interleaved(self@edge_positions)
        return(interleaved)
      }
    )
  ),
  constructor = function(vertices = default_vertices(), edges = default_edges(), faces = list(), symmetry_elements = SymmetryElementCollection()) {
    S7::new_object(
      S7::S7_object(),
      vertices = vertices,
      edges = edges,
      faces = faces,
      symmetry_elements = SymmetryElementCollection()
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
  signs <- expand.grid(
    x = c(-1, 1),
    y = c(-1, 1),
    z = c(-1, 1),
    KEEP.ALIVE = FALSE,
    stringsAsFactors = FALSE
  )
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

  faces <- list(
    F1 = c("V1", "V2", "V4", "V3"), # x = -1
    F2 = c("V5", "V6", "V8", "V7"), # x =  1
    F3 = c("V1", "V2", "V6", "V5"), # y = -1
    F4 = c("V3", "V4", "V8", "V7"), # y =  1
    F5 = c("V1", "V3", "V7", "V5"), # z = -1
    F6 = c("V2", "V4", "V8", "V6")  # z =  1
  )

  Shape(vertices = df_vertices, edges = df_edges, faces = faces)
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

  faces <- list(
    F1 = c("V1", "V3", "V5"),
    F2 = c("V1", "V3", "V6"),
    F3 = c("V1", "V4", "V5"),
    F4 = c("V1", "V4", "V6"),
    F5 = c("V2", "V3", "V5"),
    F6 = c("V2", "V3", "V6"),
    F7 = c("V2", "V4", "V5"),
    F8 = c("V2", "V4", "V6")
  )

  Shape(vertices = df_vertices, edges = df_edges, faces = faces)
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

  faces <- list(
    F1  = c("V1", "V2", "V9"),
    F2  = c("V1", "V2", "V11"),
    F3  = c("V1", "V5", "V6"),
    F4  = c("V1", "V5", "V9"),
    F5  = c("V1", "V6", "V11"),
    F6  = c("V2", "V7", "V8"),
    F7  = c("V2", "V7", "V9"),
    F8  = c("V2", "V8", "V11"),
    F9  = c("V3", "V4", "V10"),
    F10 = c("V3", "V4", "V12"),
    F11 = c("V3", "V5", "V6"),
    F12 = c("V3", "V5", "V10"),
    F13 = c("V3", "V6", "V12"),
    F14 = c("V4", "V7", "V8"),
    F15 = c("V4", "V7", "V10"),
    F16 = c("V4", "V8", "V12"),
    F17 = c("V5", "V9", "V10"),
    F18 = c("V6", "V11", "V12"),
    F19 = c("V7", "V9", "V10"),
    F20 = c("V8", "V11", "V12")
  )

  Shape(vertices = df_vertices, edges = df_edges, faces = faces)
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

  faces <- list(
    F1  = c("V1", "V9",  "V5",  "V14", "V13"),
    F2  = c("V1", "V9",  "V11", "V3",  "V17"),
    F3  = c("V1", "V13", "V2",  "V18", "V17"),
    F4  = c("V2", "V10", "V6",  "V14", "V13"),
    F5  = c("V2", "V10", "V12", "V4",  "V18"),
    F6  = c("V3", "V11", "V7",  "V16", "V15"),
    F7  = c("V3", "V15", "V4",  "V18", "V17"),
    F8  = c("V4", "V12", "V8",  "V16", "V15"),
    F9  = c("V5", "V9",  "V11", "V7",  "V19"),
    F10 = c("V5", "V14", "V6",  "V20", "V19"),
    F11 = c("V6", "V10", "V12", "V8",  "V20"),
    F12 = c("V7", "V16", "V8",  "V20", "V19")
  )

  Shape(vertices = df_vertices, edges = df_edges, faces = faces)
}

# Generics ----------------------------------------------------------------


#' @export
S7::method(as.matrix, Shape) <- function(x, ...) {
  as.matrix(x@vertices[c("x", "y", "z")])
}
