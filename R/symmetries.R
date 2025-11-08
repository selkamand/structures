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

  # browser() # CHECK too many C3 vecs
  names(C3) <- paste(c3_possibilities$vertex,c3_possibilities$face, sep = "|")

  # Compute C3 Vector & C2 vector
  c3_tetrahedron_A_BCD <- create_vector_from_start_end(start = vertices["A",], end = face_centroids["BCD",])
  c2_tetrahedron <- create_vector_from_start_end(start = edge_centers["AB",], end = edge_centers["CD",])
  c3_tetrahedron_ACD_B <- create_vector_from_start_end(start=face_centroids["ACD",], end = vertices["B",])


  return(list(vertices = vertices, edges = edges, face_centroids = face_centroids, edge_centers = edge_centers, c3_axes = C3, c2 = c2_tetrahedron))
}

create_vector_from_start_end <- function(start, end){
  vec <- end - start
  return(vec)
}
