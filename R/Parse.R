#' Read Mol2 file to a 'Structure' class object
#'
#' @param path path to the mol2 file
#'
#' @returns a Structure class object
#' @export
#'
#' @examples
#' mol2 <- system.file(package="structures", "benzene.mol2")
#' read_mol2(mol2)
read_mol2 <- function(path){
  lines <- readLines(path)
  lines <- sub(x=lines, pattern = "^[[:space:]]+", replacement = "")
  atom_marker <-"@<TRIPOS>ATOM"
  atom_marker_index <- head(grep(pattern = atom_marker, x = lines, fixed = TRUE), n = 1)
  if(length(atom_marker_index) != 1) stop("Failed to find Atom data record header (@<TRIPOS>ATOM). Are you sure your file is mol2?")

  atom_start <- atom_marker_index + 1

  lines_with_at <- grep(pattern = "^@", x = lines)
  atom_end <- lines_with_at[lines_with_at>atom_marker_index][1]-1
  if(length(atom_end) == 0) atom_end <- length(lines)

  atom_lines <- lines[atom_start:atom_end]

  ls_atoms  <- lapply(atom_lines, function(row){
    elements <- unlist(strsplit(x = row, split = "[[:space:]]+"))
    data.frame(
      eleno = pluck(elements, 1),
      elena = pluck(elements, 2),
      x = pluck(elements, 3),
      y = pluck(elements, 4),
      z = pluck(elements, 5),
      atom_type = pluck(elements, 6),
      subst_id = pluck(elements, 7),
      subst_name = pluck(elements, 8),
      charge = pluck(elements, 9),
      status_bit = pluck(elements, 10)
    )
  })

  df_atoms <- do.call(rbind, args = ls_atoms)


  return(df_atoms)

  # bond_marker <-"@<TRIPOS>BOND"
  # bond_marker_index <- head(grep(pattern = bond_marker, x = lines, fixed = TRUE), n = 1)
  # bond_start <- bond_marker_index + 1
}

pluck <- function(vec, i, outside = NA){
  if(i < 0 | i > length(vec)) {return(outside)}
  return(vec[i])
}
