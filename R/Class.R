#' Create a Molecule3D object
#'
#' Constructs an S7 object representing a single molecule with 3D coordinates.
#' In most workflows, you will not call this constructor directly — molecule
#' objects are usually created by parsers such as [structures::read_mol2()]
#'
#' @param name Character scalar. Molecule name.
#'
#' @param atoms Data frame describing atoms in the molecule.
#'   **Essential columns** (must be present on input):
#'   \itemize{
#'     \item \code{eleno} (numeric) atom ID (unique per atom)
#'     \item \code{elena} (character) element label or atom name
#'     \item \code{x}, \code{y}, \code{z} (numeric) Cartesian coordinates
#'   }
#'   **Optional columns** (automatically added if not supplied):
#'   \itemize{
#'     \item \code{element} (character) element symbol; if missing, it is derived from
#'           \code{elena} by removing digits (e.g., \code{"C3"} → \code{"C"}).
#'   }
#'   Data types are enforced so that \code{eleno}, \code{elena}, and \code{element}
#'   are character, and \code{x}, \code{y}, \code{z} are numeric.
#'   When accessed as \code{molecule@atoms}, the data frame always includes all
#'   essential columns plus \code{element}.
#'
#' @param bonds Data frame describing bonds between atoms.
#'   **Essential columns**:
#'   \itemize{
#'     \item \code{bond_id} (numeric) unique bond identifier
#'     \item \code{origin_atom_id}, \code{target_atom_id} (character) atom IDs
#'           that correspond to \code{atoms$eleno}
#'   }
#'   **Optional columns**:
#'   \itemize{
#'     \item \code{bond_type} (character) SYBYL/Tripos bond code (e.g. "1", "2", "ar");
#'           if missing, it is added with the value \code{"un"} (unknown).
#'   }
#'   All other columns are retained as-is.
#'
#' @param misc A list containing any additional metadata (e.g., provenance,
#'   notes, or debug information). Stored without modification.
#'
#' @param anchor Numeric length-3 vector \code{c(x, y, z)} specifying the molecule’s
#'   reference point (default \code{c(0,0,0)}). Used by the translation helpers
#'   (e.g., \code{\link{translate_molecule_to_origin}}, \code{\link{translate_molecule_to_position}},
#'   \code{\link{translate_molecule_by_vector}}) to reposition the molecule relative
#'   to this point. Typically set to an atom’s coordinates via
#'   \code{\link{set_anchor_by_atom}} or to an arbitrary position via
#'   \code{\link{set_anchor_by_position}}.
#'
#' @details
#' During creation, \code{atoms} and \code{bonds} are processed using internal
#' helper functions (\code{format_atoms()} and \code{format_bonds()}) to ensure
#' required columns are present, types are correct, and that all bond endpoints
#' refer to valid atoms.
#'
#' @section Anchor:
#' The \emph{anchor} is a persistent reference position stored as a length-3 numeric
#' vector (\code{c(x, y, z)}). It does not constrain geometry by itself; instead it
#' supports convenient re-positioning:
#' \itemize{
#'   \item \code{\link{set_anchor_by_position}} — set the anchor to an arbitrary position.
#'   \item \code{\link{set_anchor_by_atom}} — set the anchor to the coordinates of a given atom (\code{eleno}).
#'   \item \code{\link{translate_molecule_to_origin}} — translate the molecule so the anchor moves to \code{c(0,0,0)}.
#'   \item \code{\link{translate_molecule_to_position}} — translate the molecule so the anchor moves to a specified position.
#'   \item \code{\link{translate_molecule_by_vector}} — translate the molecule by a specified vector.
#' }
#' Functions that transform coordinates (e.g., \code{transform_molecule()}) also
#' update the anchor to keep it consistent with the new coordinate frame.
#'
#'
#' @return An S7 object of class \code{"Molecule3D"} with the following slots and properties:
#' \itemize{
#'   \item \strong{name} — character scalar, molecule name.
#'   \item \strong{atoms} — data frame with atom information
#'         (\code{eleno}, \code{elena}, \code{element}, \code{x}, \code{y}, \code{z}).
#'   \item \strong{bonds} — data frame with bond information
#'         (\code{bond_id}, \code{origin_atom_id}, \code{target_atom_id}, \code{bond_type}).
#'   \item \strong{misc} — list of arbitrary metadata.
#'   \item \strong{anchor} — numeric vector \code{c(x, y, z)}: the molecule’s reference point.
#'   \item \strong{atom_ids} — character vector of atom IDs.
#'   \item \strong{bond_ids} — numeric vector of bond IDs.
#'   \item \strong{maximum_atom_id} — numeric scalar giving the highest atom ID present.
#'   \item \strong{maximum_bond_id} — numeric scalar giving the highest bond ID present.
#'   \item \strong{atom_positions} — numeric matrix of atom coordinates
#'         (rows = atoms, columns = x/y/z).
#'   \item \strong{bond_positions} — data frame of bonds with start and end
#'         coordinates for each atom pair.
#'   \item \strong{bond_positions_interleaved} — data.frame of bond positions in interleaved format
#'         (neighbouring rows = start/end points), useful for 3D plotting.
#'   \item \strong{center} — named numeric vector (\code{c(x, y, z)}) giving the
#'         geometric center of all atoms.
#' }
#'
#' @examples
#' # Typical use via parser
#' # mol <- structures::read_mol2("benzene.mol2")
#'
#' # Direct creation + anchor operations
#' atoms <- data.frame(
#'   eleno = c("1","2"),
#'   elena = c("C","O"),
#'   x = c(0, 1.2), y = c(0, 0.1), z = c(0, -0.2)
#' )
#' bonds <- data.frame(
#'   bond_id = "b1",
#'   origin_atom_id = "1",
#'   target_atom_id = "2"
#' )
#' m <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds, anchor = c(0,0,0))
#' m <- set_anchor_by_atom(m, "1")
#' m <- translate_molecule_to_origin(m) # moves so atom "1" sits at the origin
#' m@anchor
#'
#' @seealso [structures::read_mol2()], [structures::valid_bond_types()],
#'   set_anchor_by_position, set_anchor_by_atom,
#'   translate_molecule_to_origin, translate_molecule_to_position, translate_molecule_by_vector
#' @export
Molecule3D <- S7::new_class(
  name = "Molecule3D",
  properties = list(
    name = S7::class_character,
    atoms = S7::class_data.frame,
    bonds = S7::class_data.frame,
    misc  = S7::class_list,
    anchor = S7::class_numeric,

    ## COMPUTED PROPERTIES
    # List all atom ids (eleno) described by atoms data.frame
    atom_ids = S7::new_property(
      class = S7::class_numeric,
      getter = function(self){ unique(self@atoms$eleno) },
      setter = function(self, value){stop("@atom_ids is a read only property")}
      ),

    # List all bond ids (bond_id) described by atoms data.frame
    bond_ids = S7::new_property(
      class = S7::class_numeric,
      getter = function(self){ unique(self@bonds$bond_id) },
      setter = function(self, value){stop("@bond_ids is a read only property")}
    ),

    # Maximum atom id (useful to know when combining two different molecules together or adding new atoms)
    maximum_atom_id = S7::new_property(
      class = S7::class_numeric,
      getter = function(self){max(c(0, self@atom_ids), na.rm = TRUE)},
      setter = function(self, value){stop("@maximum_atom_id is a read only property")}
    ),

    # Maximum bond id (useful to know when combining two different molecules together or adding new atoms)
    maximum_bond_id = S7::new_property(
      class = S7::class_numeric,
      getter = function(self){max(c(0, self@bond_ids), na.rm = TRUE)},
      setter = function(self, value){stop("@maximum_bond_id is a read only property")}
    ),

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
        x = mean(mx_positions[,1]),
        y = mean(mx_positions[,2]),
        z=  mean(mx_positions[,3])
      )
      return(center)
    })
  ),

  # Add/normalize columns as the object is being created
  constructor = function(name = "MyChemical", atoms = minimal_atoms(), bonds = minimal_bonds(),  misc = list(), anchor = c(0, 0, 0)) {

    # Add bond_type if not present (with all bond types set to 'unknown') and fix column types
    bonds <- format_bonds(bonds)

    # Fix column types (and add 'element' column if not present)
    atoms <- format_atoms(atoms)

    # Return the S7 object
    S7::new_object(
      S7::S7_object(),
      name=name,
      atoms = atoms,
      bonds = bonds,
      misc = misc,
      anchor = anchor
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

    if(!is.numeric(self@bonds$bond_id)) {return(sprintf("@bonds data.frame bond_id column must be a numeric vector, not a [%s]", toString(class(self@bonds$bond_id))))}

    ## ---- Validate Atom Columns ----
    required_atom_cols <- c("eleno", "elena", "element", "x", "y", "z")
    observed_atom_cols <- colnames(self@atoms)

    if (!all(required_atom_cols %in% observed_atom_cols)) {
      missing <- setdiff(required_atom_cols, observed_atom_cols)
      return(sprintf(
        "@atoms data.frame is missing required column(s): %s",
        toString(missing)
      ))
    }

    eleno <- self@atoms[["eleno"]]
    if(!is.numeric(eleno)) return(sprintf("@atoms column 'eleno' must be a numeric vector, not [%s]", toString(class(eleno))))
    if(anyNA(eleno)) return(sprintf("@atoms column 'eleno' must NOT contain any missing values. Found: [%d]", sum(is.na(eleno))))
    if(any(eleno < 1)) return(sprintf("@atoms column 'eleno' must NOT contain any values < 1. Problematic values found: [%s]", toString(unique(eleno[eleno < 1]))))
    if(any(duplicated(eleno))) return(sprintf("@atoms column 'eleno' can NOT contain duplicates. Duplicates found: [%s]", toString(eleno[duplicated(eleno)])))

    elena <- self@atoms[["elena"]]
    if(!is.character(elena)) return(sprintf("@atoms column 'elena' must be a character vector, not [%s]", toString(class(elena))))
    if(anyNA(elena)) return(sprintf("@atoms column 'elena' must NOT contain any missing values. Found: [%d]", sum(is.na(elena))))


    ## ---- Ensure all origin/target atom Ids in bonds dataframe are in atom dataframe  ----
    atom_ids <- self@atoms$eleno
    atom_ids_in_bonds_data <- unique(c(self@bonds$origin_atom_id, self@bonds$target_atom_id))
    if(!all(atom_ids_in_bonds_data %in% atom_ids)){
      return(
        sprintf("@bonds describes atoms not present in @atoms dataframe. Missing eleno's: %s", toString(setdiff(atom_ids_in_bonds_data, atom_ids)))
      )
    }

    ## ---- Validate Anchor  ----
    if(!is.numeric(self@anchor))  return(sprintf("@anchor must be a numeric vector, not a [%s]", class(self@anchor)))
    if(length(self@anchor) != 3)  return(sprintf("@anchor must be a length 3 vector with values corresponding to x, y, and z. Supplied anchor only has [%d] dimensions", length(self@anchor)))

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
    "eleno" = numeric(0),
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
# Also add an 'element' column representing elena with numbers stripped out
format_atoms <- function(atoms){
  cols <- colnames(atoms)
  if("eleno" %in% cols){ atoms[["eleno"]] <- as.numeric(atoms[["eleno"]])}
  if("elena" %in% cols){ atoms[["elena"]] <- as.character(atoms[["elena"]])}
  if("x" %in% cols){ atoms[["x"]] <- as.numeric(atoms[["x"]])}
  if("y" %in% cols){ atoms[["y"]] <- as.numeric(atoms[["y"]])}
  if("z" %in% cols){ atoms[["z"]] <- as.numeric(atoms[["z"]])}

  # Add element column
  if(!"element" %in% cols){ atoms[["element"]] <- gsub(x=atoms[["elena"]], pattern = "[0-9]", replacement = "")}

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

  if("bond_id" %in% cols){ bonds[["bond_id"]] <- as.numeric(bonds[["bond_id"]])}
  if("origin_atom_id" %in% cols){ bonds[["origin_atom_id"]] <- as.numeric(bonds[["origin_atom_id"]])}
  if("target_atom_id" %in% cols){ bonds[["target_atom_id"]] <- as.numeric(bonds[["target_atom_id"]])}
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



# Non Generic Methods ---------------------------------------------------


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
  in_eleno <- x@atoms$eleno %in% eleno
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
  in_eleno <- x@atoms$eleno %in% eleno
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

  # Also apply transformation to anchor position to
  # preserve its relative position to the rest of the molecule
  x@anchor <- transformation(x@anchor, ...)

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

#' Compute distance between two points
#'
#' @param x a Molecule3D object.
#' @param eleno1,eleno2 atom IDs (eleno's) of the two atoms to compute distance between
#'
#' @returns numeric euclidean distance between two points
#' @export
#'
#' @examples
#' path <- system.file(package="structures", "benzene.mol2")
#' molecule <- read_mol2(path)
#' compute_distance_between_atoms(molecule, 1, 2)
compute_distance_between_atoms <- function(x, eleno1, eleno2){
  assertions::assert_class(x, "structures::Molecule3D")
  valid_ids <- x@atom_ids
  if(!eleno1 %in% valid_ids) stop("atom ID (eleno) [", eleno1, "] could not be found in molecule")
  if(!eleno2 %in% valid_ids) stop("atom ID (eleno) [", eleno2, "] could not be found in molecule")

  pos1 = fetch_atom_position(x, eleno1, careful = FALSE)
  pos2 = fetch_atom_position(x, eleno2, careful = FALSE)
  sqrt(sum((pos2-pos1)^2))
}

#' Fetch atom Cartesian coordinates
#'
#' Returns the 3D position(s) of one or more atoms by \code{eleno}.
#'
#' @param x A \code{Molecule3D} object.
#' @param eleno Character or numeric vector of atom IDs (matching \code{x@atoms$eleno}).
#' @param careful Logical. If \code{TRUE} (default), validates \code{x} and that all
#'   \code{eleno} exist; otherwise missing IDs yield \code{NA} coordinates.
#'
#' @return If a single \code{eleno} is supplied, a named numeric vector
#'   \code{c(x, y, z)}. If multiple, a matrix with columns \code{x}, \code{y}, \code{z}
#'   and row names equal to the requested \code{eleno}.
#'
#' @examples
#' # Read Data
#' path <- system.file(package = "structures", "benzene.mol2")
#' molecule <- read_mol2(path)
#'
#' # Single atom position
#' # fetch_atom_position(molecule, "1")
#'
#' # Multiple atoms (rows named by eleno)
#' # fetch_atom_position(molecule, c(1, 2, 3))
#'
#' @seealso \code{\link{compute_distance_between_atoms}}
#' @export
fetch_atom_position <- function(x, eleno, careful=TRUE){
  if(careful){
    assertions::assert_class(x, "structures::Molecule3D")
    if(!all(eleno %in% x@atom_ids)) stop("atom ID (eleno) [", eleno, "] could not be found in molecule")
  }
  mx_positions <- x@atom_positions
  idx <- match(eleno, rownames(mx_positions))
  positions <- mx_positions[idx, ,drop=TRUE]

  if(length(eleno) > 1){
    rownames(positions) <- eleno
  }

  return(positions)
}

#' Set the molecule anchor by position
#'
#' Sets the \code{anchor} of a \code{Molecule3D} to an arbitrary position
#' \code{c(x, y, z)}. The anchor is a reference point used by translation
#' helpers to reposition the molecule.
#'
#' @param x A \code{Molecule3D} object.
#' @param position Numeric length-3 vector \code{c(x, y, z)}.
#'
#' @return A \code{Molecule3D} object with \code{x@anchor} set to \code{position}.
#' @examples
#' # m <- set_anchor_by_position(m, c(1, 2, 3))
#' @seealso \code{\link{set_anchor_by_atom}},
#'   \code{\link{translate_molecule_to_origin}},
#'   \code{\link{translate_molecule_to_position}},
#'   \code{\link{translate_molecule_by_vector}}
#' @export
set_anchor_by_position <- function(x, position){
  x@anchor <- position
  return(x)
}

#' Set the molecule anchor using an atom ID
#'
#' Sets the \code{anchor} of a \code{Molecule3D} to the current coordinates
#' of a specified atom (\code{eleno}).
#'
#' @param x A \code{Molecule3D} object.
#' @param eleno Atom ID (matching \code{x@atoms$eleno}).
#'
#' @return A \code{Molecule3D} object with \code{x@anchor} set to the atom position.
#' @examples
#' # m <- set_anchor_by_atom(m, "5")
#' @seealso \code{\link{set_anchor_by_position}},
#'   \code{\link{translate_molecule_to_origin}},
#'   \code{\link{translate_molecule_to_position}},
#'   \code{\link{translate_molecule_by_vector}}
#' @export
set_anchor_by_atom <- function(x, eleno){
  assertions::assert_class(x, "structures::Molecule3D")
  pos = fetch_atom_position(x, eleno)
  x@anchor <- pos
  return(x)
}


#' Translate the molecule so its anchor is at the origin
#'
#' Moves all atom coordinates by a translation that places \code{x@anchor} at
#' \code{c(0, 0, 0)}. The anchor is updated accordingly via
#' \code{\link{transform_molecule}}.
#'
#' @param x A \code{Molecule3D} object.
#'
#' @return A \code{Molecule3D} object translated so \code{x@anchor == c(0,0,0)}.
#' @examples
#' # m <- translate_molecule_to_origin(m)
#' @seealso \code{\link{translate_molecule_to_position}},
#'   \code{\link{translate_molecule_by_vector}},
#'   \code{\link{set_anchor_by_position}}, \code{\link{set_anchor_by_atom}}
#' @export
translate_molecule_to_origin <- function(x){
  assertions::assert_class(x, "structures::Molecule3D")
  translate_molecule_to_position(x, c(0, 0, 0))
}

#' Translate the molecule so its anchor matches a target position
#'
#' Moves all atom coordinates by a translation that places \code{x@anchor} at
#' \code{new_position}. The anchor is updated accordingly via
#' \code{\link{transform_molecule}}.
#'
#' @param x A \code{Molecule3D} object.
#' @param new_position Numeric length-3 vector \code{c(x, y, z)}.
#'
#' @return A \code{Molecule3D} object translated so \code{x@anchor == new_position}.
#' @examples
#' # m <- translate_molecule_to_position(m, c(5, 0, -2))
#' @seealso \code{\link{translate_molecule_to_origin}},
#'   \code{\link{translate_molecule_by_vector}},
#'   \code{\link{set_anchor_by_position}}, \code{\link{set_anchor_by_atom}}
#' @export
translate_molecule_to_position <- function(x, new_position){
  assertions::assert_class(x, "structures::Molecule3D")
  translation_vec <- new_position-x@anchor
  translate_molecule_by_vector(x, translation_vec)
}

#' Translate the molecule by a vector
#'
#' Adds a translation vector to all atom coordinates using
#' \code{\link{transform_molecule}}. The anchor is translated by the same vector,
#' keeping it consistent with the new coordinates.
#'
#' @param x A \code{Molecule3D} object.
#' @param vector Numeric length-3 vector \code{c(dx, dy, dz)}.
#'
#' @return A \code{Molecule3D} object translated by \code{vector}.
#' @examples
#' # m <- translate_molecule_by_vector(m, c(1, 0, 0))
#' @seealso \code{\link{translate_molecule_to_origin}},
#'   \code{\link{translate_molecule_to_position}},
#'   \code{\link{set_anchor_by_position}}, \code{\link{set_anchor_by_atom}}
#' @export
translate_molecule_by_vector <- function(x, vector){
  assertions::assert_class(x, "structures::Molecule3D")
  transform_molecule(x = x, transformation = function(original) { original + vector})
}


#' Combine two Molecule3D objects
#'
#' Merges the atoms and bonds of two \code{Molecule3D} objects into a single
#' object. By default, identifiers from \code{molecule2} are offset so they do
#' not clash with \code{molecule1}. A \code{source} column is added to both the
#' atoms and bonds tables indicating the originating molecule name. The anchor
#' and all other properties of \code{molecule1} are preserved.
#'
#' @param molecule1,molecule2 \code{Molecule3D} objects to combine.
#'   \code{molecule1} provides the base object (including its \code{@anchor}).
#' @param update_ids Logical (default \code{TRUE}). When \code{TRUE}, offsets the
#'   atom IDs and bond IDs of \code{molecule2} by \code{molecule1@maximum_atom_id}
#'   and \code{molecule1@maximum_bond_id}, respectively, to avoid identifier
#'   collisions. When \code{FALSE}, IDs are used as-is (and may clash).
#'
#' @details
#' The function:
#' \itemize{
#'   \item Appends the atom tables and bond tables (adding a \code{source} column to each).
#'   \item When \code{update_ids = TRUE}, increments \code{molecule2} atom IDs
#'         (\code{eleno}) and bond IDs (\code{bond_id}) before binding.
#'   \item Leaves \code{origin_atom_id} and \code{target_atom_id} internally consistent
#'         with the updated \code{eleno} values.
#'   \item Preserves \code{molecule1}'s \code{@anchor} and \code{@name}.
#' }
#'
#' All non-essential columns in the atom/bond tables are preserved. The returned
#' object inherits \code{molecule1}'s non-tabular properties.
#'
#' @return A \code{Molecule3D} object containing the combined atoms and bonds.
#'
#' @examples
#' m1 <- read_mol2(system.file("extdata", "benzene.mol2", package = "structures"))
#' m2 <- read_mol2(system.file("extdata", "benzene.mol2", package = "structures"))
#'
#' # Combine
#' m12 <- combine_molecules(m1, m2)
#'
#' m12
#'
#' @seealso \code{\link{Molecule3D}}, \code{\link{translate_molecule_to_origin}},
#'   \code{\link{set_anchor_by_atom}}
#' @export
combine_molecules <- function(molecule1, molecule2, update_ids = TRUE){
  assertions::assert_class(molecule1, "structures::Molecule3D")
  assertions::assert_class(molecule2, "structures::Molecule3D")

  atoms1 <- molecule1@atoms
  atoms2 <- molecule2@atoms
  bonds1 <- molecule1@bonds
  bonds2 <- molecule2@bonds

  if(update_ids){
    old_ids <- atoms2$eleno
    new_ids <- atoms2$eleno + molecule1@maximum_atom_id
    atoms2$eleno <- atoms2$eleno + molecule1@maximum_atom_id

    bonds2$bond_id <- if(nrow(bonds2) > 0) bonds2$bond_id + molecule1@maximum_bond_id
    bonds2$origin_atom_id <- new_ids[match(bonds2$origin_atom_id, old_ids)]
    bonds2$target_atom_id <- new_ids[match(bonds2$target_atom_id, old_ids)]
  }

  atoms1$source = molecule1@name
  atoms2$source = molecule2@name

  bonds1$source <- if(nrow(bonds1) > 0) molecule1@name else character(0)
  bonds2$source <- if(nrow(bonds2) > 0) molecule2@name else character(0)

  atoms <- dplyr::bind_rows(atoms1, atoms2)
  bonds <- dplyr::bind_rows(bonds1, bonds2)

  new <- molecule1
  new <- S7::set_props(new, atoms = atoms, bonds=bonds)

  return(new)
}
