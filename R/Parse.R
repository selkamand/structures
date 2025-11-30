
# Reading Mol2 ------------------------------------------------------------


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

  # Get mol_type and charge_type from molecule section
  mol_type <- NULL
  charge_type <- NULL
  if ("MOLECULE" %in% names(ls_contents_substplit_named)){
    char_vec <- ls_contents_substplit_named[["MOLECULE"]]
    if(length(char_vec) >= 4){
      mol_type <- char_vec[3]
      charge_type <- char_vec[4]
      # Note char_vec[1] is also the name of the molecule - so might be worth exposing a setting of how we choose name
      # Currently if name=NULL we infer from filename, which works well.
      # But it might be worth adding the option to infer from TRIPOS<MOLECULE> section
    }
  }

  chemical <- Molecule3D(
    name = name,
    atoms = df_atom,
    bonds = df_bonds,
    mol_type = mol_type,
    charge_type = charge_type,
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



# Writing Mol2 ------------------------------------------------------------

#' Write Mol2 file
#'
#' Writes a [`Molecule3D`] object to disk in mol2 format.
#'
#' @param molecule A [Molecule3D()] object.
#' @param path path to mol2 file.
#'
#' @returns invisibly returns NULL. This function is run for its side effects. .
#'
#' @export
#'
#' @examples
#' \dontshow{
#' .old_wd <- setwd(tempdir())
#' }
#'
#' # Read mol2 file
#' mol2_path <- system.file("benzene.mol2", package = "structures")
#' mol <- read_mol2(mol2_path)
#'
#' # Write minimal mol2 file
#' write_mol2(mol, "benzene.written.mol2")
#'
#' \dontshow{
#' setwd(.old_wd)
#' }
write_mol2 <- function(molecule, path){
  name <- molecule@name
  n_atoms <- nrow(molecule@atoms)
  n_bonds <- nrow(molecule@bonds)
  molecule_type <- molecule@mol_type
  charge_type <- molecule@charge_type

  creation_time <- Sys.time()

  # Write molecule header and minimal "@<TRIPOS>MOLECULE" section
  cat(
    sep = "\n",file = path, append = FALSE,
    sprintf("#\tName: %s", name),
    sprintf("#\tCreation Time: %s\n", creation_time),
    "@<TRIPOS>MOLECULE",
    name,
    sprintf("\t%d\t%d", n_atoms, n_bonds), # TODO (OPTIONAL): There are additional, optional outputs we could write here (i.e. num_substructures, features, and sets)
    molecule_type,
    charge_type,
    "" # Just to add an extra newline after molecule section
  )

  # Write @<TRIPOS>ATOM section
  write("@<TRIPOS>ATOM", file = path, append = TRUE)
  df_atoms <- molecule_to_mol2_atoms(molecule, digits = 4)

  utils::write.table(
    x = df_atoms,
    file = path,
    append = TRUE,
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep="\t"
  )

  # Write @<TRIPOS>BOND section
  write("@<TRIPOS>BOND", file = path, append = TRUE)
  df_bonds <- molecule_to_mol2_bonds(molecule)

  utils::write.table(
    x = df_bonds,
    file = path,
    append = TRUE,
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep="\t"
  )
}


# Extractors --------------------------------------------------------------

molecule_to_mol2_atoms <- function(molecule, digits = 4){
  assertions::assert_class(molecule, class = "structures::Molecule3D")
  df_atoms <- molecule@atoms
  df_atoms$x <- round(df_atoms$x, digits = digits)
  df_atoms$y <- round(df_atoms$y, digits = digits)
  df_atoms$z <- round(df_atoms$z, digits = digits)
  df_atoms <- format_numeric_columns(df_atoms, digits = digits)
  df_atoms$empty <- rep("", times = nrow(df_atoms))

  # Start with essential columns
  cols <- c("empty", "eleno", "elena", "x", "y", "z", "atom_type")

  # Add optional columns that are present in data and non-NA
  optional_columns <- c("subst_id", "subst_name", "charge", "status_bit")
  observed_cols <- colnames(df_atoms)
  for (opt in optional_columns){
    if(!opt %in% observed_cols) break
    if(all(is.na(df_atoms[[opt]]))) break
    cols <- c(cols, opt)
  }

  df_mol2_atoms <- df_atoms[, cols, drop=FALSE]

  return(df_mol2_atoms)
}


molecule_to_mol2_bonds <- function(molecule){
  assertions::assert_class(molecule, class = "structures::Molecule3D")
  df_bonds <- molecule@bonds
  df_bonds$empty <- rep("", times = nrow(df_bonds))

  # Start with essential columns
  cols <- c("empty", "bond_id", "origin_atom_id", "target_atom_id", "bond_type")

  if("status_bits" %in% colnames(df_bonds) & !all(is.na(df_bonds[["status_bits"]]))){
    cols <- c(cols, "status_bits")
  }

  df_mol2_bonds <- df_bonds[,cols, drop=FALSE]


  return(df_mol2_bonds)
}


# Utils --------------------------------------------------------------------
format_numeric_columns <- function(df, digits = 6) {
  is_num <- vapply(df, is.numeric, logical(1))
  df[is_num] <- lapply(df[is_num], function(col) {
    format(col, scientific = FALSE, trim = TRUE, digits = digits)
  })
  return(df)
}
