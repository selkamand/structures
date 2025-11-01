
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
