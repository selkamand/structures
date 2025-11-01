#' Read Mol2 file to a 'Structure' class object
#'
#' @param path path to the mol2 file
#' @param name what is the name of this molecule. If NULL will be inferred from filename
#'
#' @returns a Structure class object
#' @export
#'
#' @examples
#' mol2 <- system.file(package="structures", "benzene.mol2")
#' read_mol2(mol2)
read_mol2 <- function(path, name = NULL){
  contents <- readr::read_file(path)
  if(is.null(name)) name <- sub(x=basename(path), pattern = "\\..*$", replacement = "")

  # Remove carriage return character (\r) in case mol2 file was made on windows with CLRF line-endings
  contents <- gsub(x=contents, pattern = "\r", replacement = "")
  ls_contents <- unlist(strsplit(contents, split = "@"))
  ls_contents <- Filter(function(v){ f = nchar(v) > 0}, x=ls_contents)
  ls_contents_substplit <- lapply(ls_contents, function(x){
    val <- unlist(strsplit(x, split =  "\n"));
    val <- sub(x=val, pattern = "^[[:space:]]+", replacement = "")
    return(val)
  })

  names(ls_contents_substplit) <- vapply(ls_contents_substplit, FUN = function(v){v[1]}, FUN.VALUE = character(1))
  ls_contents_substplit_named <- lapply(ls_contents_substplit, FUN = function(v){utils::tail(v, n=-1)})

  names(ls_contents_substplit_named) <- sub(x=names(ls_contents_substplit_named), pattern = "<TRIPOS>", replacement = "")

  mol2_fields_found = names(ls_contents_substplit_named)

  if("ATOM" %in% mol2_fields_found){
    df_atom <- character_to_dataframe(ls_contents_substplit_named$ATOM, cols = c("eleno", "elena", "x", "y", "z", "atom_type", "subst_id", "subst_name", "charge", "status_bit"), header=FALSE)
  }
  else{
    stop("No ATOM field found in mol2 file: ", path)
  }

  if("BOND" %in% mol2_fields_found){
    df_bonds <- character_to_dataframe(ls_contents_substplit_named$BOND, cols = c("bond_id", "origin_atom_id", "target_atom_id", "bond_type", "status_bits"), header=FALSE)
  }
  else{
    warning("No bonds found in mol2 file: ", path)
    df_bonds <- minimal_bonds()
  }
  chemical <- Molecule3D(
    name = name,
    atoms = df_atom,
    bonds = df_bonds,
    misc = ls_contents_substplit_named
  )

  return(chemical)
}

character_to_dataframe <- function(char, cols, header = FALSE){

  df <- utils::read.table(
    text = char,
    header = header,
    sep = "",            # any amount of whitespace
    fill = TRUE,         # pad short rows
    quote = "",          # no quoting in your sample
    comment.char = "",   # keep everything
    stringsAsFactors = FALSE,
    col.names = if (!header) cols else NULL
  )

  return(df)

}
