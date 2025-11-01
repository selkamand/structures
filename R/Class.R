#' Create a Molecule3D object
#'
#' Constructs an S7 object representing a single molecule with 3D coordinates.
#'
#' @param name Character scalar. Molecule name.
#' @param atoms Data frame with required columns:
#'   \itemize{
#'     \item \code{eleno} (character) atom ID
#'     \item \code{elena} (character) element symbol (e.g., "C", "H")
#'     \item \code{x}, \code{y}, \code{z} (numeric) Cartesian coordinates
#'   }
#'   Columns are type-normalized on construction.
#' @param bonds Data frame with required columns:
#'   \itemize{
#'     \item \code{bond_id} (character)
#'     \item \code{origin_atom_id}, \code{target_atom_id} (character; match \code{atoms$eleno})
#'   }
#'   Optional \code{bond_type} (character). If missing, it is added with value \code{"un"}.
#' @param misc List of arbitrary metadata (for debugging/annotations).
#'
#' @return An S7 object of class \code{"Molecule3D"} with computed properties
#'   (e.g., \code{atom_ids}, \code{atom_positions}, \code{bond_positions},
#'   \code{bond_positions_interleaved}, \code{center}).
#'
#' @examples
#' # Minimal empty object
#' Molecule3D()
#'
#' # From small tables
#' atoms <- data.frame(
#'   eleno = c("1","2"),
#'   elena = c("C","O"),
#'   x = c(0, 1.2), y = c(0, 0.1), z = c(0, -0.2)
#' )
#' bonds <- data.frame(
#'   bond_id = c("b1"),
#'   origin_atom_id = "1",
#'   target_atom_id = "2"
#' )
#' m <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds)
#' m@center
#'
#' @export
Molecule3D <- S7::new_class(
  name = "Molecule3D",
  properties = list(
    name = S7::class_character,
    atoms = S7::class_data.frame,
    bonds = S7::class_data.frame,
    misc  = S7::class_list,

    ## COMPUTED PROPERTIES
    # Atom Ids described by atoms data.frame
    atom_ids = S7::new_property(class = S7::class_numeric, getter = function(self){ unique(self@atoms$eleno) }),

    # atom position matrix (row names are atom ids: eleno)
    atom_positions = S7::new_property(class = S7::class_numeric, getter = function(self){
      mx <- as.matrix(self@atoms[c("x", "y", "z")])
      rownames(mx) <- self@atoms[["eleno"]]
      return(mx)
    }),

    bond_positions = S7::new_property(class = S7::class_numeric, getter = function(self){
      bonds = self@bonds
      atoms = self@atoms

      # indexes of origin atoms in atom data.frame (io)
      io <- match(bonds[["origin_atom_id"]], atoms[["eleno"]])
      it <- match(bonds[["target_atom_id"]], atoms[["eleno"]])

      # Add positions to bond data.frame
      bonds$x    <- atoms[["x"]][io]
      bonds$y    <- atoms[["y"]][io]
      bonds$z    <- atoms[["z"]][io]
      bonds$xend <- atoms[["x"]][it]
      bonds$yend <- atoms[["y"]][it]
      bonds$zend <- atoms[["z"]][it]
      return(bonds)
    }),

    # Interleaved bond positions (pairs of rows = start and end point of bonds) useful for plotting in rgl
    bond_positions_interleaved = S7::new_property(class = S7::class_numeric, getter = function(self){
      to_interleaved(self@bond_positions)
    }),

    # Center position of all atoms
    center = S7::new_property(class = S7::class_numeric, getter = function(self){
      mx_positions <- self@atom_positions
      center <- c(
        x = mean(mx_positions[1]),
        y = mean(mx_positions[2]),
        z=  mean(mx_positions[3])
      )
      return(center)
    })
  ),

  # Add/normalize columns as the object is being created
  constructor = function(name = "MyChemical", atoms = minimal_atoms(), bonds = minimal_bonds(),  misc = list()) {

    # Add bond_type if not present (with all bond types set to 'unknown') and fix column types
    bonds <- format_bonds(bonds)

    # Fix column types
    atoms <- format_atoms(atoms)

    # Return the S7 object
    S7::new_object(
      S7::S7_object(),
      name=name,
      atoms = atoms,
      bonds = bonds,
      misc = misc
    )
  },

  validator = function(self) {

    ## ---- Validate Chemical Name ----
    if(!is.character(self@name)) return(sprintf("@name must be a string (length 1 chacter vector), not a %s", toString(class(self@name))))
    if(length(self@name) > 1) return(sprintf("@name must be a string, not a character vector of length %d", length(self@name)))
    if(nchar(self@name) == 0) return(sprintf("@name cannot be an empty string"))


    ## ---- Validate Bond Columns ----
    required_bond_cols <- c("bond_id", "origin_atom_id", "target_atom_id", "bond_type")
    observed_bond_cols <- colnames(self@bonds)

    if (!all(required_bond_cols %in% observed_bond_cols)) {
      missing <- setdiff(required_bond_cols, observed_bond_cols)
      return(sprintf(
        "@bonds data.frame is missing required column(s): %s",
        toString(missing)
      ))
    }

    if(!is.character(self@bonds$bond_id)) {return(sprintf("@bonds data.frame bond_id column must be a character, not %s", toString(class(self@bonds$bond_id))))}

    ## ---- Validate Atom Columns ----
    required_atom_cols <- c("eleno", "elena", "x", "y", "z")
    observed_atom_cols <- colnames(self@atoms)

    if (!all(required_atom_cols %in% observed_atom_cols)) {
      missing <- setdiff(required_atom_cols, observed_atom_cols)
      return(sprintf(
        "@atoms data.frame is missing required column(s): %s",
        toString(missing)
      ))
    }

    ## ---- Ensure all origin/target atom Ids in bonds dataframe are in atom dataframe  ----
    atom_ids <- self@atoms$eleno
    atom_ids_in_bonds_data <- unique(c(self@bonds$origin_atom_id, self@bonds$target_atom_id))
    if(!all(atom_ids_in_bonds_data %in% atom_ids)){
      return(
        sprintf("@bonds describes atoms not present in @atoms dataframe. Missing eleno's: %s", toString(setdiff(atom_ids_in_bonds_data, atom_ids)))
      )
    }


    ## If no problems:
    NULL
  }
)

#' Create a Minimal Atoms Data Frame
#'
#' Generates an empty `data.frame` containing the required columns for atoms
#' in a `Molecule3D` object. This serves as the default `atoms` input when
#' creating a new `Molecule3D`.
#'
#' @returns A `data.frame` with columns:
#' \describe{
#'   \item{eleno}{Character vector of atom element numbers (unique identifiers).}
#'   \item{elena}{Character vector of atom element names (e.g., "C", "H", "O").}
#'   \item{x}{Numeric vector for the atom's x-coordinate.}
#'   \item{y}{Numeric vector for the atom's y-coordinate.}
#'   \item{z}{Numeric vector for the atom's z-coordinate.}
#' }
#' Each column is initialized as a zero-length vector.
#'
#' @examples
#' minimal_atoms()
#'
#' @export
minimal_atoms <- function(){
  data.frame(
    "eleno" = character(0),
    "elena" = character(0),
    "x" = numeric(0),
    "y" = numeric(0),
    "z" = numeric(0)
  )
}

#' Create a Minimal Bonds Data Frame
#'
#' Generates an empty `data.frame` containing the required columns for bonds
#' in a `Molecule3D` object. This serves as the default `bonds` input when
#' creating a new `Molecule3D`.
#'
#' @returns A `data.frame` with columns:
#' \describe{
#'   \item{bond_id}{Numeric vector of bond identifiers.}
#'   \item{origin_atom_id}{Character vector of origin atom IDs (by eleno).}
#'   \item{target_atom_id}{Character vector of target atom IDs (by eleno).}
#' }
#' Each column is initialized as a zero-length vector.
#'
#' @examples
#' minimal_bonds()
#'
#' @export
minimal_bonds <- function(){
  data.frame(
    "bond_id" = numeric(0),
    "origin_atom_id" = character(0),
    "target_atom_id" = character(0)
  )
}


# Format atoms data.frame (cast required columns as the required types)
format_atoms <- function(atoms){
  cols <- colnames(atoms)
  if("eleno" %in% cols){ atoms[["eleno"]] <- as.character(atoms[["eleno"]])}
  if("elena" %in% cols){ atoms[["elena"]] <- as.character(atoms[["elena"]])}
  if("x" %in% cols){ atoms[["x"]] <- as.numeric(atoms[["x"]])}
  if("y" %in% cols){ atoms[["y"]] <- as.numeric(atoms[["y"]])}
  if("z" %in% cols){ atoms[["z"]] <- as.numeric(atoms[["z"]])}
  return(atoms)
}

# Add bond_type column and recast column
format_bonds <- function(bonds){
  # Add bond_type column if missing
  if (!"bond_type" %in% names(bonds)) {
    bonds$bond_type <- rep("un", times = nrow(bonds))
  }

  #  Recast columns
  cols <- colnames(bonds)

  if("bond_id" %in% cols){ bonds[["bond_id"]] <- as.character(bonds[["bond_id"]])}
  if("origin_atom_id" %in% cols){ bonds[["origin_atom_id"]] <- as.character(bonds[["origin_atom_id"]])}
  if("target_atom_id" %in% cols){ bonds[["target_atom_id"]] <- as.character(bonds[["target_atom_id"]])}
  if("bond_type" %in% cols){ bonds[["bond_type"]] <- as.character(bonds[["bond_type"]])}

  return(bonds)
}


#' SYBYL bond types
#'
#' @returns a character vector with all valid SYBYL bond types
#' @export
#'
#' @examples
#' valid_bond_types()
valid_bond_types <- function(){
  c("single" = "1", "double" = "2", "triple" = "3", "amide" = "am", "aromatic" = "ar", "dummy" = "du", "unknown" = "un", "not connected" = "nc")
}



# Generics ----------------------------------------------------------------
# print <- S7::new_generic("print", "x")
S7::method(print, Molecule3D) <- function(x, ...) { cat(
  sprintf(
    "===================\nChemical Molecule3D\n===================\nName: %s\nAtoms: %d\nBonds: %d\n\nSee @atoms paramater for atom positions and @bonds paramater for bonds\n",
    x@name, nrow(x@atoms), nrow(x@bonds))
)}

# as.data.frame <- S7::new_generic("as.data.frame", "x")
S7::method(as.data.frame, Molecule3D) <- function(x, ...) {
  x@atoms
}

# as.matrix <- S7::new_generic("as.matrix", "x")
S7::method(as.matrix, Molecule3D) <- function(x, ...) {
  mx <- as.matrix(x@atoms[c("x", "y", "z")])
  rownames(mx) <- x@atoms[["eleno"]]
  return(mx)
}


#' Delete atoms from molecule
#'
#' Removes atoms from molecule by atom_id (eleno).
#' Will also drop bonds that involve deleted atoms.
#'
#' @param x a Molecule3D object
#' @param eleno a vector of atom IDs to delete
#'
#' @returns a Molecule3D object
#' @export
#'
#' @examples
#' path <- system.file(package="structures", "benzene.mol2")
#' molecule <- read_mol2(path)
#' remove_atoms(molecule, c(1,2,3))
remove_atoms <- function(x, eleno){
  assertions::assert_class(x, class = "structures::Molecule3D")
  in_eleno <- x@atoms$eleno %in% as.character(eleno)
  new_atom_data <- x@atoms[!in_eleno, ,drop=FALSE]
  new_bond_data <- x@bonds[x@bonds$origin_atom_id %in% new_atom_data$eleno & x@bonds$target_atom_id %in% new_atom_data$eleno, , drop=FALSE]

  # Modify existing object by setting multiple properties simultaneously
  x <- S7::set_props(x, atoms = new_atom_data, bonds = new_bond_data)

  return(x)
}


#' Filter for atoms molecule
#'
#' Filters a molecule to include only specific atoms by atom_id (eleno).
#' Will also drop orphaned bonds. To delete atoms from a molecule, see [remove_atoms()]
#'
#' @param x a Molecule3D object
#' @param eleno a vector of atom IDs to filter. Will automatically cast to character vector.
#'
#' @returns a Molecule3D object
#' @export
#'
#' @examples
#' path <- system.file(package="structures", "benzene.mol2")
#' molecule <- read_mol2(path)
#' filter_atoms(molecule, c(1,2,3))
filter_atoms <- function(x, eleno){
  assertions::assert_class(x, class = "structures::Molecule3D")
  in_eleno <- x@atoms$eleno %in% as.character(eleno)
  new_atom_data <- x@atoms[in_eleno, ,drop=FALSE]
  new_bond_data <- x@bonds[x@bonds$origin_atom_id %in% new_atom_data$eleno & x@bonds$target_atom_id %in% new_atom_data$eleno, , drop=FALSE]

  # Update Properties
  x <- S7::set_props(x, atoms = new_atom_data, bonds = new_bond_data)

  return(x)
}


#' Apply arbitratry 3D transformations to molecule3D objects
#'
#' Applies a user-supplied transformation function to atoms in a Molecule3D object.
#' This lets you apply any function from the move package to a Molecule3D object.
#'
#' @param x a Molecule3D object
#' @param transformation A function that takes any length 3 numeric vector / named list of length 3
#'   (`x`, `y`, `z`) and returns a structure (e.g., named numeric or list)
#'   with elements `x`, `y`, and `z` in that order.
#' @param ... additional arguments to transformation function
#' @return A data frame with the same columns as `table`, but with the
#'   `x`, `y`, and `z` coordinates replaced by the transformed values.
#'
#' @examples
#' # Example: rotate points around Z-axis
#'
#' path <- system.file(package="structures", "benzene.mol2")
#' molecule <- read_mol2(path)
#'
#' rotate_z <- function(p) {
#'   angle <- pi / 4
#'   c(x = p["x"] * cos(angle) - p["y"] * sin(angle),
#'     y = p["x"] * sin(angle) + p["y"] * cos(angle),
#'     z = p["z"])
#' }
#'
#' transform_molecule(molecule, rotate_z)
#'
#' @export
transform_molecule <- function(x, transformation, ...){
  assertions::assert_class(x, "structures::Molecule3D")
  assertions::assert_function(transformation)

  mx_coords <- x@atom_positions
  res = apply(X = mx_coords, MARGIN = 1, FUN = transformation, ..., simplify = FALSE)

  x@atoms$x <- vapply(res, FUN = \(d){d[1]}, FUN.VALUE = numeric(1))
  x@atoms$y <- vapply(res, FUN = \(d){d[2]}, FUN.VALUE = numeric(1))
  x@atoms$z <- vapply(res, FUN = \(d){d[3]}, FUN.VALUE = numeric(1))
  return(x)
}

#' Locate geometrical center of molecule
#'
#' Locate geometrical center of all atoms in Molecule3D object. Can also just use center property.
#' Included as separate function to allow for tidy pipelines.
#'
#' @param x a Molecule3D object
#'
#' @returns a Molecule3D object
#' @export
#'
#' @examples
#' path <- system.file(package="structures", "benzene.mol2")
#' molecule <- read_mol2(path)
#' locate_molecule_center(molecule)
locate_molecule_center <- function(x){
  assertions::assert_class(x, "structures::Molecule3D")
  x@center
}
