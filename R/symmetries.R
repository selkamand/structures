
# Shape -------------------------------------------------------------------


Shape <- S7::new_class(
  name = "Shape",
  properties = list(
    vertices = S7::class_data.frame,
    edges = S7::class_data.frame,
    faces = S7::class_list,
    symmetry_axes = S7::class_list,

    # Computed Properties (edge centers, face centers, triangles, quads, etcs)
    edge_positions = S7::new_property(
      setter = function(self, value){ stop("@edge_positions is a read only property") },
      getter = function(self){
        vertex1 = self@edges[[1]]
        vertex2 = self@edges[[2]]
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
  constructor = function(vertices = default_vertices(), edges = default_edges(), faces = list(), symmetry_axes = list()){
    S7::new_object(
      S7::S7_object(),
      vertices  = vertices,
      edges = edges,
      faces = faces,
      symmetry_axes = symmetry_axes
    )
  }
)

## Defaults
default_vertices <- function(){
  data.frame(
    name = character(0),
    x = numeric(0),
    y = numeric(0),
    z = numeric(0),
    stringsAsFactors = FALSE
  )
}

default_edges <- function(){
  data.frame(vertex1 = character(0), vertex2 = character(0))
}


tetrahedron <- function(){
  vertices <- matrix(
    c(
      1,  1,  1,   # A
      -1, -1,  1,   # B
      -1,  1, -1,   # C
      1, -1, -1    # D
    ),
    nrow = 4L, byrow = TRUE,
    dimnames = list(c("A","B","C","D"), c("x","y","z"))
  )
  df_vertices <- as.data.frame(vertices)
  df_vertices$vertex <- rownames(vertices)

  df_edges <- data.frame(
    vertex1 = c("A","A","A","B","B","C"),
    vertex2 = c("B","C","D","C","D","D"),
    stringsAsFactors = FALSE
  )

  faces <- list(
    BCD = c("B","C","D"),
    ACD = c("A","C","D"),
    ABD = c("A","B","D"),
    ABC = c("A","B","C")
  )

  Shape(vertices = df_vertices, edges = df_edges, faces = faces)
}



#' Symmetry Tetrahedron
#'
#' @returns a list of symmetry relevant information
#' @export
#'
#' @examples
#' symmetry_tetrahedron()
symmetry_tetrahedron <- function(){
  vertices <- matrix(
    c(-1, 1, 1, -1, -1, 1, -1, 1, -1, -1, 1, 1),
    nrow = 4L,
    ncol = 3L,
    dimnames = list(c("A", "B", "C", "D"), c("x", "y", "z"))
  )

  df_vertices <- as.data.frame(vertices)
  df_vertices$vertex <- rownames(vertices)

  # rownames(vertices) <- LETTERS[seq_along(vertices[,1])]

  faces <- list(
    BCD = c(2,3,4),  # face opposite A (B,C,D)
    ACD = c(1,3,4),  # opposite B (A,C,D)
    ABD = c(1,2,4),  # opposite C (A,B,D)
    ABC = c(1,2,3)   # opposite D (A,B,C)
  )

  # Compute face centroids
  face_centroids <- t(sapply(faces, function(idx) colMeans(vertices[idx, ])))

  # Compute edge centers
  edges <- list(
    AB = c("A","B"),
    AC = c("A","C"),
    AD = c("A","D"),
    BC = c("B","C"),
    BD = c("B","D"),
    CD = c("C","D")
  )

  edge_centers <- t(sapply(edges, function(e)
    colMeans(vertices[e, , drop = FALSE])
  ))

  c3_possibilities <- expand.grid(vertex = df_vertices$vertex, vertices,face=names(faces))

  C3 <- Map(f = function(vertex, face){
    vertex_coord <- vertices[vertex, ]
    face_centroid <- face_centroids[face, ]
    create_vector_from_start_end(start = vertex_coord, end = face_centroid)
    }, c3_possibilities$vertex, c3_possibilities$face
    )

  names(C3) <- paste(c3_possibilities$vertex,c3_possibilities$face, sep = "|")

  # Compute C3 Vector & C2 vector
  c3_tetrahedron_A_BCD <- create_vector_from_start_end(start = vertices["A",], end = face_centroids["BCD",])
  c2_tetrahedron <- create_vector_from_start_end(start = edge_centers["AB",], end = edge_centers["CD",])
  c3_tetrahedron_ACD_B <- create_vector_from_start_end(start=face_centroids["ACD",], end = vertices["B",])


  return(list(vertices = vertices, edges = edges, face_centroids = face_centroids, edge_centers = edge_centers, c3_axes = C3, c2 = c2_tetrahedron))
}


# Generics ----------------------------------------------------------------


#' @export
S7::method(as.matrix, Shape) <- function(x, ...) {
  as.matrix(x@vertices[c("x", "y", "z")])
}
