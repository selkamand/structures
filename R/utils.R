
enrich_bonds_with_xyz_position <- function(bonds, atoms,
                                           origin = "origin",
                                           target = "target",
                                           atom_id = "eleno") {
  # indices of the origin/target atoms in `atoms`
  io <- match(bonds[[origin]], atoms[[atom_id]])
  it <- match(bonds[[target]], atoms[[atom_id]])

  # build output
  out <- bonds
  out$x    <- atoms[["x"]][io]
  out$y    <- atoms[["y"]][io]
  out$z    <- atoms[["z"]][io]
  out$xend <- atoms[["x"]][it]
  out$yend <- atoms[["y"]][it]
  out$zend <- atoms[["z"]][it]
  out
}

to_interleaved <- function(df, coord = c("x", "y", "z"), end_suffix = "end") {
  # Columns assumed to be coordinates
  start_cols <- coord
  end_cols <- paste0(coord, end_suffix)

  # Safety check
  missing_cols <- setdiff(c(start_cols, end_cols), names(df))
  if (length(missing_cols)) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  n <- nrow(df)

  # Point labels
  point <- rep(c("start", "end"), times = n)

  # Build interleaved coordinates
  xyz_mat <- t(cbind(df[start_cols], df[end_cols]))
  xyz_mat <- matrix(xyz_mat, ncol = length(coord), byrow = TRUE)
  xyz_df <- as.data.frame(xyz_mat, stringsAsFactors = FALSE)
  names(xyz_df) <- coord

  # Keep and replicate non-coordinate (meta) columns
  meta_cols <- setdiff(names(df), c(start_cols, end_cols))
  meta_df <- if (length(meta_cols)) {
    df[rep(seq_len(n), each = 2), meta_cols, drop = FALSE]
  } else {
    NULL
  }

  # Assemble output (segment info + meta + coords)
  out <- data.frame(
    point = point,
    meta_df,
    xyz_df,
    row.names = NULL,
    check.names = FALSE
  )

  out
}

# Convert elena to element
elena_to_element <- function(elena){
  element <- gsub(x=elena, pattern = "[0-9]", replacement = "")
  element <- gsub(x=element, pattern = "(.)[A-Z]+", replacement = "\\1")
  return(element)
}


#' Convert a Molecule3D's bonds to an igraph graph
#'
#' Builds an undirected \pkg{igraph} graph from a \code{Molecule3D}'s bond table.
#' Vertices correspond to atoms (\code{atoms$eleno}); edges correspond to bonds
#' (\code{bonds$origin_atom_id} -- \code{bonds$target_atom_id}). When there are
#' no bonds, a vertex-only graph is returned (one vertex per atom).
#'
#' @param molecule A \code{structures::Molecule3D} object.
#' @param vertex_attrs Character vector of atom columns to keep as vertex attributes.
#'   The atom ID (\code{eleno}) is always included as \code{name}. Default:
#'   \code{c("elena", "element", "x", "y", "z")}. Columns not present are ignored.
#' @param edge_attrs Character vector of bond columns to keep as edge attributes
#'   (besides endpoints). Default: \code{c("bond_id", "bond_type")}. Columns not
#'   present are ignored.
#' @param include_isolated Logical; if \code{TRUE} (default) include all atoms as
#'   vertices even if they have no bonds.
#'
#' @return An undirected \pkg{igraph} graph. Vertex attribute \code{name} holds
#'   the atom ID (\code{eleno}). Additional attributes may be present as requested.
#'
#' @examples
#' # g <- bonds_as_igraph(m)
#' # igraph::components(g)
#'
#' @importFrom igraph graph_from_data_frame
#' @export
as_igraph <- function(molecule,
                            vertex_attrs = c("elena", "element", "x", "y", "z"),
                            edge_attrs   = c("bond_id", "bond_type"),
                            include_isolated = TRUE) {

  rlang::check_installed(pkg = "igraph", reason = "to convert molecular structure to igraph format")
  assertions::assert_class(molecule, "structures::Molecule3D")

  # --- Vertices (atoms) ---
  atoms <- molecule@atoms
  if (!is.data.frame(atoms) || nrow(atoms) == 0) {
    stop("molecule@atoms is empty; cannot build a graph without atoms.")
  }

  # Build vertex data.frame: igraph expects a 'name' column if vertices are supplied
  v_keep <- intersect(vertex_attrs, colnames(atoms))
  vdf <- atoms[, c("eleno", v_keep), drop = FALSE]
  # igraph usually treats 'name' as character; keep a numeric copy too if you like
  vdf$name <- as.character(vdf$eleno)

  # Reorder so 'name' is first column (igraph convention when a data.frame is used)
  vdf <- vdf[, c("name", setdiff(colnames(vdf), c("name"))), drop = FALSE]

  # --- Edges (bonds) ---
  bonds <- molecule@bonds
  if (is.null(bonds) || !is.data.frame(bonds) || nrow(bonds) == 0) {
    # Vertex-only graph (all atoms as isolated vertices)
    g <- igraph::graph_from_data_frame(d = NULL, directed = FALSE, vertices = vdf)
    return(g)
  }

  # Keep requested edge attributes if present
  e_keep <- intersect(edge_attrs, colnames(bonds))
  edf <- data.frame(
    from = bonds$origin_atom_id,
    to   = bonds$target_atom_id,
    stringsAsFactors = FALSE
  )
  if (length(e_keep)) {
    edf <- cbind(edf, bonds[, e_keep, drop = FALSE])
  }

  # Make sure endpoints are character if vertices use 'name' as character
  edf$from <- as.character(edf$from)
  edf$to   <- as.character(edf$to)

  # Build graph; include all atoms as vertices when requested
  if (include_isolated) {
    g <- igraph::graph_from_data_frame(d = edf, directed = FALSE, vertices = vdf)
  } else {
    g <- igraph::graph_from_data_frame(d = edf, directed = FALSE)
  }

  return(g)
}

# Once its an igraph we can use igraph::components(mode = "weak")


