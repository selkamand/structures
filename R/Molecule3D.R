# Molecule3D Class -------------------------------------------------------------------


#' Create a Molecule3D object
#'
#' Constructs an S7 object representing a single molecule with 3D coordinates.
#' In most workflows, you will not call this constructor directly molecule
#' objects are usually created by parsers such as [structures::read_mol2()].
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
#'           \code{elena} by removing digits (e.g., \code{"C3"} <U+2192> \code{"C"}).
#'      \item \code{atom_type} valid SYBYL atom types. See [valid_atom_types()] for valid values. If not supplied will automatically be set to 'Any'.
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
#'     \item \code{origin_atom_id}, \code{target_atom_id} (numeric) atom IDs
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
#' @param anchor Numeric length-3 vector c(x, y, z) specifying the molecule<U+2019>s
#'   reference point. If NULL (default), it is set to the geometric center of
#'   all atoms at construction time; if no atoms are present it falls back to
#'   c(0, 0, 0). Used by translation helpers
#'   (e.g., \code{\link{translate_molecule_to_origin}}, \code{\link{translate_molecule_to_position}},
#'   \code{\link{translate_molecule_by_vector}}) to reposition the molecule relative
#'   to this point. Typically set to an atom<U+2019>s coordinates via
#'   \code{\link{set_anchor_by_atom}} or to an arbitrary position via
#'   \code{\link{set_anchor_by_position}}.
#'
#' @param symmetry_elements A [`structures::SymmetryElementCollection`] object describing a set of symmetry axes
#' embedded in the molecule
#'
#' @details
#' During creation, \code{atoms} and \code{bonds} are processed using internal
#' helpers (\code{format_atoms()} and \code{format_bonds()}) to ensure
#' required columns are present, types are correct, and that all bond endpoints
#' refer to valid atoms. If \code{symmetry_elements} is supplied, it is validated
#' to be a [`structures::SymmetryElementCollection`]. The class also exposes derived,
#' read-only properties related to symmetry:
#' \itemize{
#'   \item \code{@contains_symmetry_axes} logical flag indicating whether any symmetry elements have been added to this molecule.
#' }
#' Coordinate transforms applied via \code{\link{transform_molecule}} will also
#' transform the endpoints of each stored symmetry axis so they remain consistent
#' with atom coordinates.
#'
#' @section Anchor:
#' The \emph{anchor} is a persistent reference position stored as a length-3 numeric
#' vector (\code{c(x, y, z)}). It does not constrain geometry by itself; instead it
#' supports convenient re-positioning:
#' \itemize{
#'   \item \code{\link{set_anchor_by_position}} <U+2014> set the anchor to an arbitrary position.
#'   \item \code{\link{set_anchor_by_atom}} <U+2014> set the anchor to the coordinates of a given atom (\code{eleno}).
#'   \item \code{\link{translate_molecule_to_origin}} <U+2014> translate so the anchor moves to \code{c(0,0,0)}.
#'   \item \code{\link{translate_molecule_to_position}} <U+2014> translate so the anchor moves to a specified position.
#'   \item \code{\link{translate_molecule_by_vector}} <U+2014> translate by a specified vector.
#' }
#' Functions that transform coordinates (e.g., \code{transform_molecule()}) also
#' update the anchor to keep it consistent with the new coordinate frame.
#'
#' @return An S7 object of class \code{"Molecule3D"} with (among others) the properties:
#' \itemize{
#'   \item \strong{name}, \strong{atoms}, \strong{bonds}, \strong{misc}, \strong{anchor}
#'   \item \strong{atom_ids}, \strong{bond_ids}, \strong{maximum_atom_id}, \strong{maximum_bond_id}
#'   \item \strong{atom_positions}, \strong{bond_positions}, \strong{bond_positions_interleaved}
#'   \item \strong{center}
#'   \item \strong{symmetry_elements} A [`structures::SymmetryElementCollection`] object describing the symmetry axes of this molecule.
#'   \item \strong{contains_symmetry_elements} logical (read-only)
#' }
#'
#' @examples
#' # Read mol2 file into a Molecule3D object
#' benzene <- read_mol2(system.file("benzene.mol2", package = "structures"))
#' print(benzene)
#'
#' # Directly create Molecule3D
#' # object from atoms and bonds data.frames
#' atoms <- data.frame(
#'   eleno = c(1, 2),
#'   elena = c("C", "O"),
#'   x = c(0, 1.2), y = c(0, 0.1), z = c(0, -0.2)
#' )
#' bonds <- data.frame(
#'   bond_id = 1,
#'   origin_atom_id = 1,
#'   target_atom_id = 2
#' )
#'
#' m <- Molecule3D(
#'   name = "CO", atoms = atoms, bonds = bonds,
#'   anchor = c(0, 0, 0)
#' )
#'
#' # Add a proper rotation axis
#' proper_rotation_axis <- ProperRotationAxis(
#'   n = 3L,
#'   posA = c(0, 0, 0),
#'   posB = c(0, 0, 1),
#'   label = "C3(z)"
#' )
#' m <- add_symmetry_element_to_molecule(m, proper_rotation_axis)
#' print(m)
#'
#' @seealso
#' [structures::read_mol2()],
#' [structures::ProperRotationAxis],
#' [add_symmetry_element_to_molecule()],
#' [transform_molecule()]
#' @export
Molecule3D <- S7::new_class(
  name = "Molecule3D",


  ## Properties --------------------------------------------------------------
  properties = list(

    ### Core Configurable Properties   ---------------------------------------------
    name = S7::class_character,
    atoms = S7::class_data.frame,
    bonds = S7::class_data.frame,
    misc = S7::class_list,
    anchor = S7::new_property(
      S7::class_numeric,
      setter = function(self, value) {
        # If anchor vector has x, y, and z names sort vector in xyz order
        current_names <- names(value)
        if (all(c("x", "y", "z") %in% current_names)) {
          value <- value[c("x", "y", "z")]
        }
        # Otherwise set the names to x,y,z
        else {
          names(value) <- c("x", "y", "z")[seq_along(value)]
        }

        # Set anchor and return self
        self@anchor <- value
        return(self)
      }
    ),

    ### Computed, Read Only  ---------------------------------------------
    # List all atom ids (eleno) described by atoms data.frame
    atom_ids = S7::new_property(
      class = S7::class_numeric,
      getter = function(self) {
        unique(self@atoms$eleno)
      },
      setter = function(self, value) {
        stop("@atom_ids is a read only property")
      }
    ),

    # List all bond ids (bond_id) described by atoms data.frame
    bond_ids = S7::new_property(
      class = S7::class_numeric,
      getter = function(self) {
        unique(self@bonds$bond_id)
      },
      setter = function(self, value) {
        stop("@bond_ids is a read only property")
      }
    ),


    # Maximum atom id (useful to know when combining two different molecules together or adding new atoms)
    maximum_atom_id = S7::new_property(
      class = S7::class_numeric,
      getter = function(self) {
        max(c(0, self@atom_ids), na.rm = TRUE)
      },
      setter = function(self, value) {
        stop("@maximum_atom_id is a read only property")
      }
    ),

    # Maximum bond id (useful to know when combining two different molecules together or adding new atoms)
    maximum_bond_id = S7::new_property(
      class = S7::class_numeric,
      getter = function(self) {
        max(c(0, self@bond_ids), na.rm = TRUE)
      },
      setter = function(self, value) {
        stop("@maximum_bond_id is a read only property")
      }
    ),

    # atom position matrix (row names are atom ids: eleno)
    atom_positions = S7::new_property(
      class = S7::class_numeric,
      setter = function(self, value) {
        stop("@atom_positions is a read only property")
      },
      getter = function(self) {
        mx <- as.matrix(self@atoms[c("x", "y", "z")])
        rownames(mx) <- self@atoms[["eleno"]]
        return(mx)
      }
    ),
    bond_positions = S7::new_property(
      class = S7::class_numeric, getter = function(self) {
        bonds <- self@bonds
        atoms <- self@atoms

        # indexes of origin atoms in atom data.frame (io)
        io <- match(bonds[["origin_atom_id"]], atoms[["eleno"]])
        it <- match(bonds[["target_atom_id"]], atoms[["eleno"]])

        # Add positions to bond data.frame
        bonds$x <- atoms[["x"]][io]
        bonds$y <- atoms[["y"]][io]
        bonds$z <- atoms[["z"]][io]
        bonds$xend <- atoms[["x"]][it]
        bonds$yend <- atoms[["y"]][it]
        bonds$zend <- atoms[["z"]][it]

        # Add middle positions (useful for labelling)
        bonds[["x_middle"]] <- pmean(bonds$x, bonds$xend)
        bonds[["y_middle"]] <- pmean(bonds$y, bonds$yend)
        bonds[["z_middle"]] <- pmean(bonds$z, bonds$zend)

        return(bonds)
      },
      setter = function(self, value) {
        stop("@bond_positions is a read only property")
      }
    ),

    # Interleaved bond positions (pairs of rows = start and end point of bonds) useful for plotting in rgl
    bond_positions_interleaved = S7::new_property(
      class = S7::class_numeric,
      getter = function(self) {
        to_interleaved(self@bond_positions)
      },
      setter = function(self, value) {
        stop("@bond_positions_interleaved is a read only property")
      }
    ),

    # Center position of all atoms
    center = S7::new_property(class = S7::class_numeric, getter = function(self) {
      mx_positions <- self@atom_positions
      center <- c(
        x = mean(mx_positions[, 1], na.rm = TRUE),
        y = mean(mx_positions[, 2], na.rm = TRUE),
        z = mean(mx_positions[, 3], na.rm = TRUE)
      )
      return(center)
    }),

    # Center position of all atoms
    n_dummy_atoms = S7::new_property(
      class = S7::class_numeric,
      setter = function(self, value) { stop("@n_dummy_atoms is a read only property") },
      getter = function(self) {
        sum(self@atoms[["atom_type"]] %in% c("Du", "Du.C"))
      }),

    #### Connectivity ---------------------------------------------
    # Returns a list of connected clusters, each containing a numeric vector of eleno representing members of each cluster.
    # Requires igraph
    connectivity = S7::new_property(
      class = S7::class_numeric,
      getter = function(self) {
        graph <- as_igraph(self)
        components <- igraph::components(graph, mode = "weak")
        lapply(split(names(components$membership), components$membership), as.numeric)
      },
      setter = function(self, value) {
        stop("@components is a read only property")
      }
    ),


    ### Symmetry Related Properties ---------------------------------------------
    # Symmetry elements
    symmetry_elements = S7::new_property(
      class = SymmetryElementCollection
    ),

    # Boolean: does the atom have any symmetry elements
    contains_symmetry_elements = S7::new_property(
      class = S7::class_logical,
      getter = function(self) {
        self@symmetry_elements@n_elements > 0
      },
      setter = function(self, value) {
        stop("@contains_symmetry_elements is a read only property")
      }
    )
  ),

  ## Validator ---------------------------------------------
  validator = function(self) {
    ## ---- Validate Chemical Name ----
    if (!is.character(self@name)) {
      return(sprintf("@name must be a string (length 1 chacter vector), not a %s", toString(class(self@name))))
    }
    if (length(self@name) > 1) {
      return(sprintf("@name must be a string, not a character vector of length %d", length(self@name)))
    }
    if (nchar(self@name) == 0) {
      return(sprintf("@name cannot be an empty string"))
    }


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

    bond_ids <- self@bonds$bond_id
    if (!is.numeric(bond_ids)) {
      return(sprintf("@bonds data.frame bond_id column must be a numeric vector, not a [%s]", toString(class(self@bonds$bond_id))))
    }
    if (anyNA(bond_ids)) {
      return(sprintf("@bonds column 'bond_id' must NOT contain any missing values. Found: [%d]", sum(is.na(bond_ids))))
    }
    if (any(bond_ids < 1)) {
      return(sprintf("@bonds column 'bond_id' must NOT contain any values < 1. Problematic values found: [%s]", toString(unique(bond_ids[bond_ids < 1]))))
    }
    if (any(duplicated(bond_ids))) {
      return(sprintf("@bonds column 'bond_id' can NOT contain duplicates. Duplicates found: [%s]", toString(bond_ids[duplicated(bond_ids)])))
    }
    if (!is.numeric(bond_ids)) {
      return(sprintf("@bonds data.frame 'bond_id' column must be a numeric vector, not a [%s]", toString(class(self@bonds$bond_id))))
    }

    ## ---- Validate Atom Columns ----
    required_atom_cols <- c("eleno", "elena", "element", "x", "y", "z", "atom_type")
    observed_atom_cols <- colnames(self@atoms)

    if (!all(required_atom_cols %in% observed_atom_cols)) {
      missing <- setdiff(required_atom_cols, observed_atom_cols)
      return(sprintf(
        "@atoms data.frame is missing required column(s): %s",
        toString(missing)
      ))
    }

    eleno <- self@atoms[["eleno"]]
    if (!is.numeric(eleno)) {
      return(sprintf("@atoms column 'eleno' must be a numeric vector, not [%s]", toString(class(eleno))))
    }
    if (anyNA(eleno)) {
      return(sprintf("@atoms column 'eleno' must NOT contain any missing values. Found: [%d]", sum(is.na(eleno))))
    }
    if (any(eleno < 1)) {
      return(sprintf("@atoms column 'eleno' must NOT contain any values < 1. Problematic values found: [%s]", toString(unique(eleno[eleno < 1]))))
    }
    if (any(duplicated(eleno))) {
      return(sprintf("@atoms column 'eleno' can NOT contain duplicates. Duplicates found: [%s]", toString(eleno[duplicated(eleno)])))
    }

    elena <- self@atoms[["elena"]]
    if (!is.character(elena)) {
      return(sprintf("@atoms column 'elena' must be a character vector, not [%s]", toString(class(elena))))
    }
    if (anyNA(elena)) {
      return(sprintf("@atoms column 'elena' must NOT contain any missing values. Found: [%d]", sum(is.na(elena))))
    }

    atom_type <- self@atoms[["atom_type"]]
    if(!all(atom_type %in% valid_atom_types())) {
      return(sprintf("@atoms column 'atom_type' must only contain valid SYBYL atom types. Invalid values: [%s]. See `valid_atom_types()` for a valid list)", toString(setdiff(atom_type, valid_atom_types()))))
    }

    ## ---- Ensure all origin/target atom Ids in bonds dataframe are in atom dataframe  ----
    atom_ids <- self@atoms$eleno
    atom_ids_in_bonds_data <- unique(c(self@bonds$origin_atom_id, self@bonds$target_atom_id))
    if (!all(atom_ids_in_bonds_data %in% atom_ids)) {
      return(
        sprintf("@bonds describes atoms not present in @atoms dataframe. Missing eleno's: %s", toString(setdiff(atom_ids_in_bonds_data, atom_ids)))
      )
    }

    ## ---- Validate Anchor  ----
    if (!is.numeric(self@anchor)) {
      return(sprintf("@anchor must be a numeric vector, not a [%s]", class(self@anchor)))
    }
    if (length(self@anchor) != 3) {
      return(sprintf("@anchor must be a length 3 vector with values corresponding to x, y, and z. Supplied anchor only has [%d] dimensions", length(self@anchor)))
    }

    ## If no problems:
    NULL
  },

  ## Constructor ---------------------------------------------
  # Add/normalize columns as the object is being created
  constructor = function(name = "MyChemical", atoms = minimal_atoms(), bonds = minimal_bonds(), symmetry_elements = SymmetryElementCollection(), misc = list(), anchor = NULL) {
    # Add bond_type if not present (with all bond types set to 'unknown') and fix column types
    bonds <- format_bonds(bonds)

    # Fix column types (and add 'element' column if not present)
    atoms <- format_atoms(atoms)

    # If anchor is NULL, default to the geometric center of atoms (or 0,0,0 if empty)
    if (is.null(anchor)) {
      if (nrow(atoms) > 0) {
        center_x <- mean(atoms$x, na.rm = TRUE)
        center_y <- mean(atoms$y, na.rm = TRUE)
        center_z <- mean(atoms$z, na.rm = TRUE)
        anchor <- c(center_x, center_y, center_z)
      } else {
        anchor <- c(0, 0, 0)
      }
    }

    # Return the S7 object
    S7::new_object(
      S7::S7_object(),
      name = name,
      atoms = atoms,
      bonds = bonds,
      misc = misc,
      symmetry_elements = symmetry_elements,
      anchor = anchor
    )
  }
)


# Constructor helpers -----------------------------------------------------


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
minimal_atoms <- function() {
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
minimal_bonds <- function() {
  data.frame(
    "bond_id" = numeric(0),
    "origin_atom_id" = character(0),
    "target_atom_id" = character(0)
  )
}


# Format atoms data.frame (cast required columns as the required types)
# Adds an atom_type column if not already present - and defaults values to 'Any'
# Also add an 'element' column representing elena with numbers stripped out
format_atoms <- function(atoms) {
  cols <- colnames(atoms)
  if ("eleno" %in% cols) {
    atoms[["eleno"]] <- as.numeric(atoms[["eleno"]])
  }
  if ("elena" %in% cols) {
    atoms[["elena"]] <- as.character(atoms[["elena"]])
  }
  if ("x" %in% cols) {
    atoms[["x"]] <- as.numeric(atoms[["x"]])
  }
  if ("y" %in% cols) {
    atoms[["y"]] <- as.numeric(atoms[["y"]])
  }
  if ("z" %in% cols) {
    atoms[["z"]] <- as.numeric(atoms[["z"]])
  }

  # If atom_type column is not present, make one (and set all values to 'Any' - a valid SYBYL atom type)
  if(!"atom_type" %in% cols){
    atoms["atom_type"] <- rep("Any", times = nrow(atoms))
  }
  # Otherise, cast as character
  else{
   atoms[["atom_type"]] <- as.character(atoms[["atom_type"]])

   # Fix some commonly variable / off-spec atom types seen in the wild
   # (e.g. capitalise oxygen in sulfur type
   atoms[["atom_type"]] <- ifelse(atoms[["atom_type"]] == "S.o", "S.O", atoms[["atom_type"]])
   atoms[["atom_type"]] <- ifelse(atoms[["atom_type"]] == "S.o2", "S.O2", atoms[["atom_type"]])

   # Unexpected Atom Types
   unexpected_types <- setdiff(atoms[["atom_type"]], valid_atom_types())
   if(length(unexpected_types) > 0){
     warning("Unexpected atom types found: [", toString(unexpected_types), "]. Replacing these with 'Any'")
     atoms[["atom_type"]][atoms[["atom_type"]] %in% unexpected_types] <- "Any"
   }
  }

  # Add element column
  if (!"element" %in% cols) {
    atoms[["element"]] <- elena_to_element(atoms[["elena"]])
  }

  return(atoms)
}

# Add bond_type column and recast column
format_bonds <- function(bonds) {
  # Add bond_type column if missing
  if (!"bond_type" %in% names(bonds)) {
    bonds$bond_type <- rep("un", times = nrow(bonds))
  }

  #  Recast columns
  cols <- colnames(bonds)

  if ("bond_id" %in% cols) {
    bonds[["bond_id"]] <- as.numeric(bonds[["bond_id"]])
  }
  if ("origin_atom_id" %in% cols) {
    bonds[["origin_atom_id"]] <- as.numeric(bonds[["origin_atom_id"]])
  }
  if ("target_atom_id" %in% cols) {
    bonds[["target_atom_id"]] <- as.numeric(bonds[["target_atom_id"]])
  }
  if ("bond_type" %in% cols) {
    bonds[["bond_type"]] <- as.character(bonds[["bond_type"]])
  }

  return(bonds)
}


#' SYBYL bond types
#'
#' @returns a character vector with all valid SYBYL bond types
#' @export
#'
#' @examples
#' valid_bond_types()
valid_bond_types <- function() {
  c("single" = "1", "double" = "2", "triple" = "3", "amide" = "am", "aromatic" = "ar", "dummy" = "du", "unknown" = "un", "not connected" = "nc")
}

#' SYBYL bond types
#'
#' @returns a character vector with all valid SYBYL atom types
#' @export
#'
#' @examples
#' valid_atom_types()
valid_atom_types <- function(){
  c(
    "carbon sp3" = "C.3",
    "carbon sp2" = "C.2",
    "carbon sp" = "C.1",
    "carbon aromatic" = "C.ar",
    "carbocation (C+) used only in a guadinium group" = "C.cat",
    "dummy atom" = "Du",
    "nitrogen sp3" = "N.3",
    "nitrogen sp2" = "N.2",
    "nitrogen sp" = "N.1",
    "nitrogen aromatic" = "N.ar",
    "nitrogen amide" = "N.am",
    "nitrogen trigonal planar" = "N.pl3",
    "nitrogen sp3 positively charged" = "N.4",
    "hydrogen" = "H",
    "hydrogen in Single Point Charge (SPC) water model" = "H.spc",
    "hydrogen in Transferable intermolecular Potential (TIP3P) water model" = "H.t3p",
    "lone pair" = "LP",
    "dummy carbon" = "Du.C",
    "any atom" = "Any",
    "halogen" = "Hal",
    "heteroatom = N, O, S, P" = "Het",
    "heavy atom (non hydrogen)" = "Hev",
    "lithium" = "Li",
    "sodium" = "Na",
    "oxygen sp3" = "O.3",
    "oxygen sp2" = "O.2",
    "oxygen in carboxylate and phosphate groups" = "O.co2",
    "oxygen in Single Point Charge (SPC) water model" = "O.spc",
    "oxygen in Transferable Intermolecular Potential (TIP3P) water model" = "O.t3p",
    "sulfur sp3" = "S.3",
    "sulfur sp2" = "S.2",
    "sulfoxide sulfur" = "S.O",
    "sulfone sulfur" = "S.O2",
    "phosphorous sp3" = "P.3",
    "fluorine" = "F",
    "chlorine" = "Cl",
    "bromine" = "Br",
    "iodine" = "I",
    "tin" = "Sn",
    "magnesium" = "Mg",
    "aluminum" = "Al",
    "silicon" = "Si",
    "potassium" = "K",
    "calcium" = "Ca",
    "chromium (tetrahedral)" = "Cr.th",
    "chromium (octahedral)" = "Cr.oh",
    "manganese" = "Mn",
    "iron" = "Fe",
    "cobalt (octahedral)" = "Co.oh",
    "copper" = "Cu",
    "zinc" = "Zn",
    "selenium" = "Se",
    "molybdenum" = "Mo"
  )
}


# Generics ----------------------------------------------------------------
# print <- S7::new_generic("print", "x")
S7::method(print, Molecule3D) <- function(x, ...) {
  symmetry_collection <- x@symmetry_elements
  # proper_rotation_axes@unique_proper_axis_orders

  # message("HELLO: ", toString(proper_rotation_axes@unique_proper_axis_orders))

  cat(
    sep = "",
    "===================\n",
    "Chemical Molecule3D\n",
    "===================\n",
    sprintf("Name: %s\n", x@name),
    sprintf("Atoms: %d (%d dummy atoms)\n", nrow(x@atoms), x@n_dummy_atoms),
    sprintf("Bonds: %d\n", nrow(x@bonds)),
    sprintf("Symmetry Elements: %d\n", length(x@symmetry_elements@elements)),
    sprintf("Symmetry Axes: %d\n", length(x@symmetry_elements@proper_rotation_axes)),
    sprintf("Symmetry Orders (Cn): %s\n", toString(symmetry_collection@unique_proper_axis_orders)),
    "\n-------------------\n",
    "See @atoms paramater for atom positions\n",
    "See @bonds paramater for bond positions\n",
    "See @symmetry_elements for symmetry elements\n"
  )


  return(invisible(x))
}

# as.data.frame <- S7::new_generic("as.data.frame", "x")
S7::method(as.data.frame, Molecule3D) <- function(x, ...) {
  x@atoms
}

# as.matrix <- S7::new_generic("as.matrix", "x")
S7::method(as.matrix, Molecule3D) <- function(x, ...) {
  mx <- as.matrix(x@atoms[c("x", "y", "z")])
  rownames(mx) <- as.character(x@atoms[["eleno"]])
  return(mx)
}




# Non Generic Methods ---------------------------------------------------

#' Check if an object is a Molecule3D
#'
#' @description
#' Tests whether an object inherits from the [`structures::Molecule3D`] class.
#'
#' @param x An object to test.
#'
#' @return A logical scalar: `TRUE` if `x` is a `Molecule3D`, otherwise `FALSE`.
#'
#' @examples
#' atoms <- data.frame(eleno = c(1, 2), elena = c("C", "O"), x = c(0, 1), y = c(0, 0), z = c(0, 0))
#' bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)
#' mol <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds)
#' is_molecule(mol)
#' is_molecule(123)
#'
#' @export
is_molecule <- function(x) {
  inherits(x, "structures::Molecule3D")
}



## Modifying Atoms ---------------------------------------------------------

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
#' path <- system.file(package = "structures", "benzene.mol2")
#' molecule <- read_mol2(path)
#' remove_atoms(molecule, c(1, 2, 3))
remove_atoms <- function(x, eleno) {
  assertions::assert_class(x, class = "structures::Molecule3D")
  in_eleno <- x@atoms$eleno %in% eleno
  new_atom_data <- x@atoms[!in_eleno, , drop = FALSE]
  new_bond_data <- x@bonds[x@bonds$origin_atom_id %in% new_atom_data$eleno & x@bonds$target_atom_id %in% new_atom_data$eleno, , drop = FALSE]

  # Modify existing object by setting multiple properties simultaneously
  x <- S7::set_props(x, atoms = new_atom_data, bonds = new_bond_data)

  return(x)
}

#' Remove Dummy atoms
#'
#' Remove dummy atoms from a molecule.
#' Dummy are identified based on 'atom_type' (according to SYBYL mol2 specification should be set to 'Du' / 'Du.C' for dummy atoms).
#'
#' @param x a [Molecule3D()] object
#'
#' @returns a [Molecule3D()] object with dummy atoms (and any bonds involving dummy atoms) removed
#' @export
#'
#' @examples
#' path <- system.file(package = "structures", "fe_dummies.mol2")
#' molecule_with_dummies <- read_mol2(path)
#' remove_dummy_atoms(molecule_with_dummies)
remove_dummy_atoms <- function(x, dummy_type = c("Du", "Du.C")){
  assertions::assert_class(x, class = "structures::Molecule3D")
  dummy_eleno <- x@atoms$eleno[x@atoms$atom_type %in% dummy_type]
  remove_atoms(x, dummy_eleno)
}

#' Filter molecule for specific atoms
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
#' path <- system.file(package = "structures", "benzene.mol2")
#' molecule <- read_mol2(path)
#' filter_atoms(molecule, c(1, 2, 3))
filter_atoms <- function(x, eleno) {
  assertions::assert_class(x, class = "structures::Molecule3D")
  in_eleno <- x@atoms$eleno %in% eleno
  new_atom_data <- x@atoms[in_eleno, , drop = FALSE]
  new_bond_data <- x@bonds[x@bonds$origin_atom_id %in% new_atom_data$eleno & x@bonds$target_atom_id %in% new_atom_data$eleno, , drop = FALSE]

  # Update Properties
  x <- S7::set_props(x, atoms = new_atom_data, bonds = new_bond_data)

  return(x)
}


#' Filter molecule for specific atoms by name
#'
#' Filters a molecule to include only atoms whose `elena` (atom name) matches
#' one or more requested values. Orphaned bonds (i.e. bonds whose endpoints are
#' not both retained) are automatically dropped, in the same way as
#' [filter_atoms()].
#'
#' @param x A [`structures::Molecule3D`] object.
#' @param elena Character vector of atom names to keep (matching
#'   `x@atoms$elena`).
#'
#' @return A [`structures::Molecule3D`] object containing only atoms whose
#'   `elena` is in `elena`, and the bonds between them.
#'
#' @examples
#' path <- system.file(package = "structures", "benzene.mol2")
#' molecule <- read_mol2(path)
#'
#' # Keep only atoms whose name is "C"
#' molecule_C <- filter_atoms_by_name(molecule, elena = "C")
#'
#' @export
filter_atoms_by_name <- function(x, elena) {
  assertions::assert_class(x, class = "structures::Molecule3D")
  # Reuse existing helper to find matching IDs
  eleno <- fetch_eleno_by_name(x, elena)
  filter_atoms(x, eleno = eleno)
}


#' Filter molecule for specific atoms by element
#'
#' Filters a molecule to include only atoms whose `element` (chemical symbol)
#' matches one or more requested values. Orphaned bonds (i.e. bonds whose
#' endpoints are not both retained) are automatically dropped, in the same way
#' as [filter_atoms()].
#'
#' @param x A [`structures::Molecule3D`] object.
#' @param element Character vector of element symbols to keep (matching
#'   `x@atoms$element`, e.g. `"C"`, `"H"`, `"O"`).
#'
#' @return A [`structures::Molecule3D`] object containing only atoms whose
#'   `element` is in `element`, and the bonds between them.
#'
#' @examples
#' path <- system.file(package = "structures", "benzene.mol2")
#' molecule <- read_mol2(path)
#'
#' # Keep only carbon atoms
#' molecule_C <- filter_atoms_by_element(molecule, element = "C")
#'
#' @export
filter_atoms_by_element <- function(x, element) {
  assertions::assert_class(x, class = "structures::Molecule3D")
  # Reuse existing helper to find matching IDs
  eleno <- fetch_eleno_by_element(x, element)
  filter_atoms(x, eleno = eleno)
}

## Adding Symmetry Elements ---------------------------------------------------------



#' Append a symmetry axis to a Molecule3D
#'
#' @description
#' Adds a single [`structures::SymmetryElement`] object to the \code{@symmetry_elements} paramater
#' of a [`structures::Molecule3D`]. The element is appended (order preserved) and given a unique ID.
#'
#' @param molecule A [`structures::Molecule3D`] object.
#' @param symmetry_element A [`structures::SymmetryElement`] object to append.
#'
#' @return A \code{Molecule3D} object with \code{symmetry_element} appended to
#'   \code{@symmetry_elements}.
#'
#' @examples
#' atoms <- data.frame(
#'   eleno = c(1, 2),
#'   elena = c("C", "O"),
#'   x = c(0, 1), y = c(0, 0), z = c(0, 0)
#' )
#' bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)
#' m <- Molecule3D("CO", atoms = atoms, bonds = bonds)
#' ax <- ProperRotationAxis(n = 2, posA = c(0, 0, 0), posB = c(0, 0, 1))
#' m <- add_symmetry_element_to_molecule(m, ax)
#'
#' length(m@symmetry_elements@proper_rotation_axes) # 1
#' m@symmetry_elements@unique_proper_axis_orders # 2
#'
#' @export
add_symmetry_element_to_molecule <- function(molecule, symmetry_element){

  assertions::assert_class(molecule, "structures::Molecule3D")
  assertions::assert_class(symmetry_element, "structures::SymmetryElement")

  molecule@symmetry_elements <- add_symmetry_element_to_collection(
    collection = molecule@symmetry_elements,
    new = symmetry_element
  )

  return(molecule)
}

## Fetch -------------------------------------------------------------------

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
fetch_atom_position <- function(x, eleno, careful = TRUE) {
  if (careful) {
    assertions::assert_class(x, "structures::Molecule3D")
    if (!all(eleno %in% x@atom_ids)) stop("atom ID (eleno) [", eleno, "] could not be found in molecule")
  }
  mx_positions <- x@atom_positions
  idx <- match(eleno, rownames(mx_positions))
  positions <- mx_positions[idx, , drop = TRUE]

  if (length(eleno) > 1) {
    rownames(positions) <- eleno
  }

  return(positions)
}



#' Fetch atom identifiers by name
#'
#' @param x a Molecule3D object
#' @param elena element name
#'
#' @returns a numeric vector of atom identifiers (eleno) corresponding to element names (elena).
#' @export
#'
#' @examples
#' path <- system.file(package = "structures", "benzene.mol2")
#' molecule <- read_mol2(path)
#' fetch_eleno_by_name(molecule, "C")
#'
fetch_eleno_by_name <- function(x, elena) {
  assertions::assert_class(x, class = "structures::Molecule3D")
  unique(unlist(x@atoms[x@atoms$elena %in% elena, "eleno", drop = FALSE]))
}

#' Fetch atom identifiers by element
#'
#' @param x a Molecule3D object
#' @param element element name
#'
#' @returns a numeric vector of atom identifiers (eleno) corresponding to element names (elena).
#' @export
#'
#' @examples
#' path <- system.file(package = "structures", "benzene.mol2")
#' molecule <- read_mol2(path)
#' fetch_eleno_by_element(molecule, "C")
#'
fetch_eleno_by_element <- function(x, element) {
  assertions::assert_class(x, class = "structures::Molecule3D")
  unique(unlist(x@atoms[x@atoms$element %in% element, "eleno", drop = FALSE]))
}

#' Fetch atom identifiers by SYBYL atom type
#'
#' @param x A [`structures::Molecule3D`] object.
#' @param atom_type Character vector of SYBYL atom types to match
#'   (values in `x@atoms$atom_type`; see [valid_atom_types()]).
#'
#' @returns A numeric vector of atom identifiers (`eleno`) whose
#'   `atom_type` matches any value in `atom_type`.
#'
#' @export
#'
#' @examples
#' path <- system.file(package = "structures", "fe_dummies.mol2")
#' molecule <- read_mol2(path)
#'
#' # Get Iron and dummy atoms
#' fetch_eleno_by_atom_type(molecule, c("Fe", "Du", "Du.C"))
fetch_eleno_by_atom_type <- function(x, atom_type){
  assertions::assert_class(x, class = "structures::Molecule3D")
  unique(unlist(x@atoms[x@atoms$atom_type %in% atom_type, "eleno", drop = FALSE]))
}


#' Fetch all atoms on one side of a bond
#'
#' Given a bond ID and a "direction" atom on that bond, returns the set of atom
#' IDs (\code{eleno}) that remain connected to the \emph{direction} atom after
#' virtually breaking the bond (i.e., removing the other endpoint atom). This is
#' useful for identifying the substructure that would rotate or translate when
#' treating the chosen bond as a rotatable hinge.
#'
#' @param molecule A \code{Molecule3D} object.
#' @param bond_id Numeric bond identifier (must exist in \code{molecule@bond_ids}).
#' @param direction_atom_id Numeric atom ID (\code{eleno}) that is one endpoint of
#'   \code{bond_id}. The returned cluster is the connected component containing
#'   this atom after the opposite endpoint is removed.
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Validates that \code{bond_id} exists and locates its two endpoint atom IDs.
#'   \item Checks that \code{direction_atom_id} is one of the endpoints; otherwise
#'         an error is thrown.
#'   \item Removes the \emph{other} endpoint atom from the molecule via
#'         \code{\link{remove_atoms}} (this simulates breaking the bond).
#'   \item Uses the molecule's \code{@connectivity} property (list of connected
#'         components, as numeric \code{eleno} vectors) to select the component
#'         that contains \code{direction_atom_id}.
#' }
#'
#' If no disconnection occurs (i.e., only one connected component remains),
#' the function errors with \code{"All points remain connected."}
#'
#' @return A numeric vector of atom IDs (\code{eleno}) forming the downstream
#'   connected component that includes \code{direction_atom_id}.
#'
#' @note The \code{@connectivity} property relies on \pkg{igraph} (via an
#'   internal conversion to an igraph and \code{igraph::components()}).
#'
#' @seealso \code{\link{remove_atoms}}, \code{\link{compute_distance_between_atoms}}
#'
#' @examples
#' m <- read_mol2(system.file("fe_tripod.mol2", package = "structures"))
#' # Suppose bond 4 connects atoms 10 and 12, and we want the side attached to atom 10:
#' ids <- fetch_eleno_downstream_of_bond(molecule = m, bond_id = 16, direction_atom_id = 20)
#' ids
#'
#' @export
fetch_eleno_downstream_of_bond <- function(molecule, bond_id, direction_atom_id) {
  assertions::assert_class(molecule, "structures::Molecule3D")
  assertions::assert_includes(molecule@bond_ids, required = bond_id)

  connected_atom_ids <- fetch_eleno_connected_by_bond(molecule, bond_id = bond_id)
  if (!direction_atom_id %in% connected_atom_ids) stop(sprintf("direction atom id: [%s] is not connected by bond id [%s]. Valid atom IDs are: [%s]", direction_atom_id, bond_id, toString(connected_atom_ids)))

  other_atom_id <- setdiff(connected_atom_ids, direction_atom_id)
  molecule_broken <- remove_atoms(molecule, other_atom_id)
  clusters <- molecule_broken@connectivity

  if (length(clusters) == 1) stop("All points remain connected.")

  chosen_cluster <- clusters[vapply(clusters, function(ids) {
    direction_atom_id %in% ids
  }, FUN.VALUE = logical(1))]
  cluster_atom_ids <- unname(unlist(utils::head(chosen_cluster, n = 1)))
  return(cluster_atom_ids)
}

#' Fetch the two atoms connected by a bond
#'
#' Returns the numeric atom IDs (\code{eleno}) of the two atoms joined by a given
#' bond in a \code{Molecule3D} object. This is a convenience helper to quickly
#' inspect bond connectivity or to use as input for other geometric operations.
#'
#' @param molecule A \code{Molecule3D} object.
#' @param bond_id Numeric bond identifier (must exist in \code{molecule@bond_ids}).
#'
#' @return A numeric vector of length 2 giving the \code{origin_atom_id} and
#'   \code{target_atom_id} for the specified bond.
#'
#' @seealso \code{\link{fetch_eleno_downstream_of_bond}}, \code{\link{add_bonds}}
#'
#' @examples
#' m <- read_mol2(system.file("fe_tripod.mol2", package = "structures"))
#' fetch_eleno_connected_by_bond(m, bond_id = 16)
#'
#' @export
fetch_eleno_connected_by_bond <- function(molecule, bond_id) {
  assertions::assert_class(molecule, "structures::Molecule3D")
  unname(unlist(molecule@bonds[molecule@bonds$bond_id %in% bond_id, c("origin_atom_id", "target_atom_id")]))
}

#' Fetch a single symmetry element by ID
#'
#' @param molecule A [Molecule3D()] object.
#' @param id the ID of a symmetry element.
#' @inheritParams fetch_symmetry_element_from_collection
#' @returns A [ProperRotationAxis()] object
#' @export
#'
#' @examples
#' # Minimal molecule
#' atoms <- data.frame(
#'   eleno = c(1, 2),
#'   elena = c("C", "O"),
#'   x = c(0, 1), y = c(0, 0), z = c(0, 0)
#' )
#' bonds <- data.frame(
#'   bond_id = 1,
#'   origin_atom_id = 1,
#'   target_atom_id = 2
#' )
#' m <- Molecule3D("CO", atoms = atoms, bonds = bonds)
#'
#' # Append a symmetry axis (ID is assigned automatically)
#' m <- add_symmetry_element_to_molecule(m, ProperRotationAxis(n = 2, posA = c(0, 0, 0), posB = c(0, 0, 1)))
#'
#' # Inspect available axis IDs
#' m@symmetry_elements@summary
#'
#' # Fetch the first axis by its ID
#' ax_id <- m@symmetry_elements@summary$ids[1]
#' ax <- fetch_symmetry_element_from_molecule(m, ax_id)
#' ax
fetch_symmetry_element_from_molecule <- function(molecule, id, error_if_missing = TRUE) {
  collection <- molecule@symmetry_elements
  fetch_symmetry_element_from_collection(collection, id = id, error_if_missing = error_if_missing)
}


#' Fetch the ID of the first proper rotation axis of order Cn
#'
#' @description
#' Returns the identifier (ID) of the **first** stored
#' [`structures::ProperRotationAxis`] within a [`structures::Molecule3D`]
#' whose order (`@n`) matches the user-supplied value `Cn`.
#'
#' This is a convenience helper for quickly selecting a representative
#' symmetry axis of a given order (e.g. the “principal” C\eqn{_2}, C\eqn{_3},
#' C\eqn{_6}, etc.).
#'
#' @param molecule A [`structures::Molecule3D`] object containing one or more
#'   symmetry elements.
#' @param Cn Integer (or integer-like numeric). The order of the proper rotation
#'   axis to search for (e.g. `2`, `3`, `6`).
#'
#' @return
#' A **character scalar** giving the ID of the first matching proper rotation
#' axis, or `NA_character_` if none exist.
#'
#' @details
#' The function operates on the molecule's internal
#' `@symmetry_elements@elements` list, which stores each symmetry element under
#' a unique ID (the list name). It:
#'
#' \enumerate{
#'   \item Iterates over all stored symmetry elements.
#'   \item Filters only those that are
#'         [`structures::ProperRotationAxis`] objects.
#'   \item Checks whether their order (`axis@n`) matches `Cn`.
#'   \item Returns the ID (list name) of the **first** match.
#' }
#'
#' If the molecule contains no proper rotation axes with order `Cn`, the function
#' returns `NA_character_`.
#'
#' @examples
#' atoms <- data.frame(
#'   eleno = c(1,2),
#'   elena = c("C","O"),
#'   x = c(0,1), y = c(0,0), z = c(0,0)
#' )
#' m <- Molecule3D("CO", atoms = atoms, bonds = minimal_bonds())
#'
#' m <- add_symmetry_element_to_molecule(
#'   m,
#'   ProperRotationAxis(n = 2L, posA = c(0,0,0), posB = c(0,0,1), label = "C2(z)")
#' )
#'
#' fetch_id_of_first_proper_rotation_axis_with_order(m, Cn = 2)
#'
#' @seealso
#'   [`fetch_all_proper_rotation_axes_with_order()`],
#'   [`fetch_symmetry_element_from_molecule()`],
#'   [`structures::ProperRotationAxis`]
#'
#' @export
fetch_id_of_first_proper_rotation_axis_with_order <- function(molecule, Cn){
  elements <- molecule@symmetry_elements
  is_proper_rotation_axis_with_order <- vapply(elements@elements, function(el){
    if(!is_proper_rotation_axis(el)) return(FALSE)
    return(el@n %in% Cn)
  }, FUN.VALUE = logical(1))

  elements@ids[is_proper_rotation_axis_with_order][1]
}

#' Fetch all proper rotation axes of a given order (Cn) from a Molecule3D
#'
#' @description
#' Convenience wrapper around
#' [fetch_all_proper_rotation_axes_with_order_from_collection()] that operates
#' directly on a [`structures::Molecule3D`] by using its
#' \code{@symmetry_elements} collection.
#'
#' @param molecule A [`structures::Molecule3D`] object.
#' @param Cn Integer (or integer-like numeric) fold/order to match
#'   (e.g. \code{2}, \code{3}, \code{6}).
#'
#' @return
#' \itemize{
#'   \item \code{NULL} if the molecule contains no symmetry elements at all.
#'   \item Otherwise, a \emph{list} of [`structures::ProperRotationAxis`]
#'         objects of order \code{Cn}. The list may be empty if no such axes
#'         exist.
#' }
#'
#' @examples
#' atoms <- data.frame(
#'   eleno = c(1, 2),
#'   elena = c("C","O"),
#'   x = c(0, 1), y = c(0, 0), z = c(0, 0)
#' )
#' bonds <- data.frame(
#'   bond_id = 1,
#'   origin_atom_id = 1,
#'   target_atom_id = 2
#' )
#' m <- Molecule3D("CO", atoms = atoms, bonds = bonds)
#'
#' m <- add_symmetry_element_to_molecule(
#'   m,
#'   ProperRotationAxis(n = 2L, posA = c(0, 0, 0), posB = c(0, 0, 1), label = "C2(z)")
#' )
#'
#' fetch_all_proper_rotation_axes_with_order(m, Cn = 2L)
#'
#' @export
fetch_all_proper_rotation_axes_with_order <- function(molecule, Cn) {
  assertions::assert_class(molecule, "structures::Molecule3D")

  fetch_all_proper_rotation_axes_with_order_from_collection(
    collection = molecule@symmetry_elements,
    Cn = Cn
  )
}


## Computations ------------------------------------------------------------

#' Compute a best-fit plane through all atoms of a Molecule3D
#'
#' Extracts atomic Cartesian coordinates from a [`structures::Molecule3D`]
#' object and computes the best-fit plane using
#' [`move::compute_plane_from_points()`]. If only 3 atoms are present, it computes an exact solution.
#' If there are more than 3 points it computes the best-fit plane by SVD.
#'
#' @param molecule A [`structures::Molecule3D`] object.
#'
#' @return A list describing the plane (normal, offset, and point centroid). See [move::compute_plane_from_points()] for details,
#'
#' @examples
#' m <- read_mol2(system.file("benzene.mol2", package = "structures"))
#' plane <- compute_plane_from_atoms(m)
#' plane$normal
#'
#' @export
compute_plane_from_atoms <- function(molecule){
  assertions::assert_class(molecule, "structures::Molecule3D")
  mx <- molecule@atom_positions
  move::compute_plane_from_points(mx)
}

## Transformations ---------------------------------------------------------
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
#' path <- system.file(package = "structures", "benzene.mol2")
#' molecule <- read_mol2(path)
#'
#' rotate_z <- function(p) {
#'   angle <- pi / 4
#'   c(
#'     x = p["x"] * cos(angle) - p["y"] * sin(angle),
#'     y = p["x"] * sin(angle) + p["y"] * cos(angle),
#'     z = p["z"]
#'   )
#' }
#'
#' transform_molecule(molecule, rotate_z)
#'
#' @export
transform_molecule <- function(x, transformation, ...) {
  assertions::assert_class(x, "structures::Molecule3D")
  assertions::assert_function(transformation)

  mx_coords <- x@atom_positions
  res <- apply(X = mx_coords, MARGIN = 1, FUN = transformation, ..., simplify = FALSE)

  x@atoms$x <- vapply(res, FUN = \(d){
    d[1]
  }, FUN.VALUE = numeric(1))
  x@atoms$y <- vapply(res, FUN = \(d){
    d[2]
  }, FUN.VALUE = numeric(1))
  x@atoms$z <- vapply(res, FUN = \(d){
    d[3]
  }, FUN.VALUE = numeric(1))

  # Also apply transformation to anchor position to
  # preserve its relative position to the rest of the molecule
  x@anchor <- transformation(x@anchor, ...)

  # Also apply transformation to all symmetry elements to
  # preserve their relative position to the rest of the molecule.
  elements <- x@symmetry_elements@elements
  elements_new <- lapply(elements, transform_symmetry_element, transformation = transformation, ...)

  x@symmetry_elements <- S7::set_props(
    x@symmetry_elements,
    elements = elements_new
  )

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
#' m <- read_mol2(system.file("benzene.mol2", package = "structures"))
#' locate_molecule_center(m)
locate_molecule_center <- function(x) {
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
#' m <- read_mol2(system.file("benzene.mol2", package = "structures"))
#' compute_distance_between_atoms(m, eleno1 = 1, eleno2 = 2)
compute_distance_between_atoms <- function(x, eleno1, eleno2) {
  assertions::assert_class(x, "structures::Molecule3D")
  valid_ids <- x@atom_ids
  if (!eleno1 %in% valid_ids) stop("atom ID (eleno) [", eleno1, "] could not be found in molecule")
  if (!eleno2 %in% valid_ids) stop("atom ID (eleno) [", eleno2, "] could not be found in molecule")

  pos1 <- fetch_atom_position(x, eleno1, careful = FALSE)
  pos2 <- fetch_atom_position(x, eleno2, careful = FALSE)
  sqrt(sum((pos2 - pos1)^2))
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
#' m <- read_mol2(system.file("benzene.mol2", package = "structures"))
#' m <- set_anchor_by_position(m, c(1, 2, 3))
#' @seealso \code{\link{set_anchor_by_atom}},
#'   \code{\link{translate_molecule_to_origin}},
#'   \code{\link{translate_molecule_to_position}},
#'   \code{\link{translate_molecule_by_vector}}
#' @export
set_anchor_by_position <- function(x, position) {
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
#' m <- read_mol2(system.file("benzene.mol2", package = "structures"))
#' m <- set_anchor_by_atom(m, eleno = 5)
#' @seealso \code{\link{set_anchor_by_position}},
#'   \code{\link{translate_molecule_to_origin}},
#'   \code{\link{translate_molecule_to_position}},
#'   \code{\link{translate_molecule_by_vector}}
#' @export
set_anchor_by_atom <- function(x, eleno) {
  assertions::assert_class(x, "structures::Molecule3D")
  pos <- fetch_atom_position(x, eleno)
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
translate_molecule_to_origin <- function(x) {
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
#' m <- read_mol2(system.file("benzene.mol2", package = "structures"))
#' m <- translate_molecule_to_position(m, c(5, 0, -2))
#' @seealso \code{\link{translate_molecule_to_origin}},
#'   \code{\link{translate_molecule_by_vector}},
#'   \code{\link{set_anchor_by_position}}, \code{\link{set_anchor_by_atom}}
#' @export
translate_molecule_to_position <- function(x, new_position) {
  assertions::assert_class(x, "structures::Molecule3D")
  translation_vec <- new_position - x@anchor
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
#' m <- read_mol2(system.file("benzene.mol2", package = "structures"))
#' m <- translate_molecule_by_vector(m, c(1, 0, 0))
#'
#' @seealso \code{\link{translate_molecule_to_origin}},
#'   \code{\link{translate_molecule_to_position}},
#'   \code{\link{set_anchor_by_position}}, \code{\link{set_anchor_by_atom}}
#' @export
translate_molecule_by_vector <- function(x, vector) {
  assertions::assert_class(x, "structures::Molecule3D")
  transform_molecule(x = x, transformation = function(original) {
    original + vector
  })
}

#' Translate molecule in a given direction
#'
#' Translates all atom coordinates in a [`structures::Molecule3D`] object
#' along a specified direction by a given distance.
#'
#' @param x A [Molecule3D()] object
#' @inheritParams move::translate_position_in_direction
#'
#' @returns A [Molecule3D()] object translated along \code{direction} by \code{magnitude} units
#' @export
#'
#' @examples
#' # m <- translate_molecule_in_direction(m, direction = c(1, 0, 0), magnitude = 2)
translate_molecule_in_direction <- function(x, direction, magnitude){
  transform_molecule(
    x = x,
    transformation = move::translate_position_in_direction,
    direction = direction,
    magnitude = magnitude
  )
}

#' Rotate a Molecule3D around an arbitrary axis through a point
#'
#' Rotates all atom coordinates in a [`structures::Molecule3D`] about a 3D axis
#' that passes through a specified point, by a given angle (radians). The
#' rotation is applied via [`transform_molecule()`], using
#' `move::rotate_vector_around_axis_through_point()` per-row. The molecule's
#' `@anchor` and any stored symmetry axes are rotated consistently.
#'
#' @param molecule A [`structures::Molecule3D`] object.
#' @param axis Numeric length-3 vector giving the axis direction (need not be unit).
#' @param position Numeric length-3 vector giving a point on the axis (the pivot line).
#'   Defaults to \code{c(0, 0, 0)}.
#' @param angle Numeric scalar in radians (positive = right-hand rule about \code{axis}).
#'
#' @return A \code{Molecule3D} object with rotated coordinates, anchor, and symmetry axes.
#'
#' @details
#' This is a convenience wrapper around \code{\link{transform_molecule}()} that
#' delegates the point-wise rotation to
#' \code{move::rotate_vector_around_axis_through_point(p, rotation_axis, point_on_axis, angle)}.
#' The axis direction is interpreted up to scale; only its orientation matters.
#'
#' @examples
#' # Rotate benzene 45<U+00B0> around the global z-axis through the origin
#' # m <- structures::read_mol2(system.file("benzene.mol2", package = "structures"))
#' # m_rot <- rotate_molecule_around_vector(m, axis = c(0,0,1), position = c(0,0,0), angle = pi/4)
#'
#' @seealso \code{\link{transform_molecule}}, \code{move::rotate_vector_around_axis_through_point}
#' @export
rotate_molecule_around_vector <- function(molecule, axis, position = c(0, 0, 0), angle) {
  transform_molecule(
    molecule,
    transformation = move::rotate_vector_around_axis_through_point,
    rotation_axis = axis,
    point_on_axis = position,
    angle = angle
  )
}



## Transformations Around ProperRotationAxis ------------------------------------------

#' Rotate a molecule so that a symmetry axis aligns with a target vector
#'
#' Rotates an entire [`structures::Molecule3D`] object so that the direction of
#' a stored symmetry element (typically a [`structures::ProperRotationAxis`])
#' aligns with a user-supplied 3D target vector.
#'
#' This is achieved by:
#' \enumerate{
#'   \item Extracting the symmetry element specified by `symmetry_element_id`.
#'   \item Computing the rotation (axis + angle) required to align the element’s
#'         direction vector to `target`, using
#'         [`move::rotate_vector_to_align_with_target()`].
#'   \item Applying the resulting rotation to the entire molecule via
#'         [`rotate_molecule_around_vector()`], including atoms, the anchor, and
#'         all symmetry elements.
#' }
#'
#' @details
#' This helper is useful when a molecule has annotated symmetry elements and
#' needs to be re-oriented so one of its symmetry axes matches a desired global
#' direction (e.g., align a C\eqn{_n} axis with the z-axis).
#'
#' The symmetry element’s direction is taken from its internal `@direction`
#' vector, which is automatically normalised.
#' The **target vector** does not need to be unit length.
#'
#' The rotation behaviour (parallel, anti-parallel, numerical tolerance, etc.)
#' follows the rules in [`move::rotate_vector_to_align_with_target()`].
#'
#' @param molecule A [`structures::Molecule3D`] object.
#' @param symmetry_element_id Character ID of the symmetry element to align.
#'   Must exist in `molecule@symmetry_elements`.
#' @param target Numeric length-3 vector giving the desired direction in space.
#' @param careful Assert that inputs are as expected? (flag. Default: TRUE)
#' @return A [`structures::Molecule3D`] object rotated so that the selected symmetry
#' element’s direction matches `target`.
#'
#' @examples
#' m <- read_mol2(system.file("benzene.mol2", package = "structures"))
#'
#' # Add a C6 axis (id assigned automatically)
#' m <- add_symmetry_element_to_molecule(
#'   m,
#'   ProperRotationAxis(n = 6, posA = c(0,0,-1), posB = c(0,0,1))
#' )
#'
#' # Get the axis ID
#' ax_id <- m@symmetry_elements@summary$ids[1]
#'
#' # Align it with the global x-axis
#' m2 <- rotate_molecule_so_symmetry_axis_aligns_with_vector(
#'   molecule = m,
#'   symmetry_element_id = ax_id,
#'   target = c(1, 0, 0)
#' )
#'
#' @seealso
#'   [`move::rotate_vector_to_align_with_target()`],
#'   [`rotate_molecule_around_vector()`],
#'   [`fetch_symmetry_element_from_molecule()`]
#'
#' @export
rotate_molecule_so_symmetry_axis_aligns_with_vector <- function(molecule, symmetry_element_id, target, careful=TRUE){
  if(careful){
    assertions::assert_class(molecule, "structures::Molecule3D")
    assertions::assert_length(symmetry_element_id, length = 1)
    assertions::assert_length(target, length = 3)
  }

  sym_element <- fetch_symmetry_element_from_molecule(molecule, id = symmetry_element_id, error_if_missing = TRUE)

  if(careful) assertions::assert_class(sym_element, class = c("structures::ProperRotationAxis", "structures::ImproperRotationAxis"))

  # Rotate molecule to align with target
  params <- move::rotate_vector_to_align_with_target(sym_element@direction, target = target, return = "axis_plus_angle")
  rotate_molecule_around_vector(molecule = molecule, axis = params$axis, angle = params$angle)
}


#' Rotate a molecule around a stored symmetry axis
#'
#' Rotates all atoms in a [`structures::Molecule3D`] object about one of its
#' stored symmetry elements (typically a [`structures::ProperRotationAxis`])
#' by a specified angle, using the axis direction and a point on the axis
#' as the rotation line.
#'
#' Internally, this is a convenience wrapper around
#' [`rotate_molecule_around_vector()`], which in turn delegates to
#' `move::rotate_vector_around_axis_through_point()` to apply a
#' Rodrigues-style rotation to each atomic coordinate.
#'
#' @param molecule A [`structures::Molecule3D`] object to rotate.
#' @param symmetry_element_id Character ID of the symmetry element to use as the
#'   rotation axis. This must correspond to an element in
#'   `molecule@symmetry_elements`.
#' @param angle Numeric scalar giving the rotation angle in radians.
#'   Positive angles follow the right-hand rule about the symmetry axis
#'   direction (`sym_element@direction`).
#' @param careful Logical; if `TRUE` (default), performs basic input checks
#'   (e.g. that `molecule` is a `Molecule3D`). Set to `FALSE` to skip checks
#'   in performance-critical code.
#'
#' @return A [`structures::Molecule3D`] object with all coordinates, the anchor,
#'   and any symmetry elements rotated about the chosen symmetry axis by
#'   `angle` radians.
#'
#' @details
#' The symmetry axis is defined by a direction vector (`@direction`) and a point
#' on the axis (`@posA`). This function:
#' \enumerate{
#'   \item Retrieves the symmetry element via
#'         [`fetch_symmetry_element_from_molecule()`].
#'   \item Uses its `@direction` and `@posA` as the rotation axis and point.
#'   \item Calls [`rotate_molecule_around_vector()`] with the supplied `angle`
#'         to rotate the entire molecule.
#' }
#'
#' This is particularly useful for applying rotations that respect the molecule's
#' annotated point-group symmetry (e.g. performing a C\eqn{_n} rotation).
#'
#' @examples
#' m <- read_mol2(system.file("benzene.mol2", package = "structures"))
#'
#' # Add a C6 axis along z (ID assigned automatically)
#' m <- add_symmetry_element_to_molecule(
#'   m,
#'   ProperRotationAxis(
#'     n    = 6L,
#'     posA = c(0, 0, -1),
#'     posB = c(0, 0,  1),
#'     label = "C6(z)"
#'   )
#' )
#'
#' # Fetch the ID of the first symmetry element
#' ax_id <- m@symmetry_elements@summary$ids[1]
#'
#' # Rotate the molecule by 60 degrees (2*pi/6) about that axis
#' m_rot <- rotate_molecule_around_symmetry_axis(
#'   molecule           = m,
#'   symmetry_element_id = ax_id,
#'   angle              = 2 * pi / 6
#' )
#'
#' @seealso
#'   [`rotate_molecule_around_vector()`],
#'   [`fetch_symmetry_element_from_molecule()`],
#'   `move::rotate_vector_around_axis()`,
#'   `move::rotate_vector_around_axis_through_point()`
#'
#' @export
rotate_molecule_around_symmetry_axis <- function(molecule, symmetry_element_id, angle, careful=TRUE){
  if(careful) assertions::assert_class(molecule, class = "structures::Molecule3D")
  sym_element <- fetch_symmetry_element_from_molecule(
    molecule,
    id = symmetry_element_id,
    error_if_missing = TRUE
  )

  if(careful) assertions::assert_class(sym_element, class = c("structures::ProperRotationAxis", "structures::ImproperRotationAxis"))

  # Rotate about axis
  rotate_molecule_around_vector(
    molecule,
    axis = sym_element@direction,
    position = sym_element@posA,
    angle = angle
  )
}



## Operators (Adding and subtracting molecules)---------------------------------------------------------------

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
#' m1 <- read_mol2(system.file("benzene.mol2", package = "structures"))
#' m2 <- read_mol2(system.file("benzene.mol2", package = "structures"))
#'
#' # Combine
#' m12 <- combine_molecules(m1, m2)
#'
#' m12
#'
#' @seealso \code{\link{Molecule3D}}, \code{\link{translate_molecule_to_origin}},
#'   \code{\link{set_anchor_by_atom}}
#' @export
combine_molecules <- function(molecule1, molecule2, update_ids = TRUE) {
  assertions::assert_class(molecule1, "structures::Molecule3D")
  assertions::assert_class(molecule2, "structures::Molecule3D")

  atoms1 <- molecule1@atoms
  atoms2 <- molecule2@atoms
  bonds1 <- molecule1@bonds
  bonds2 <- molecule2@bonds

  if (update_ids) {
    old_ids <- atoms2$eleno
    new_ids <- atoms2$eleno + molecule1@maximum_atom_id
    atoms2$eleno <- atoms2$eleno + molecule1@maximum_atom_id

    bonds2$bond_id <- if (nrow(bonds2) > 0) bonds2$bond_id + molecule1@maximum_bond_id
    bonds2$origin_atom_id <- new_ids[match(bonds2$origin_atom_id, old_ids)]
    bonds2$target_atom_id <- new_ids[match(bonds2$target_atom_id, old_ids)]
  }

  atoms1$source <- molecule1@name
  atoms2$source <- molecule2@name

  bonds1$source <- if (nrow(bonds1) > 0) molecule1@name else character(0)
  bonds2$source <- if (nrow(bonds2) > 0) molecule2@name else character(0)

  atoms <- dplyr::bind_rows(atoms1, atoms2)
  bonds <- dplyr::bind_rows(bonds1, bonds2)

  # Combine Symmetry Elements
  new_symmetry_elements <- combine_symmetry_element_collections(
    molecule1@symmetry_elements,
    molecule2@symmetry_elements
  )

  new <- molecule1
  new <- S7::set_props(new, atoms = atoms, bonds = bonds, symmetry_elements = new_symmetry_elements)

  return(new)
}


#' Add one or more bonds to a Molecule3D object
#'
#' Appends one or more new bond records to a \code{Molecule3D} object. Each bond
#' connects a specified origin atom (\code{origin_atom_id}) to one or more
#' target atoms (\code{target_atom_ids}). Bond IDs are automatically assigned by
#' incrementing from the molecule<U+2019>s current \code{@maximum_bond_id}.
#'
#' @param molecule A \code{Molecule3D} object.
#' @param origin_atom_id Numeric atom ID (matching \code{atoms$eleno}) that serves
#'   as the origin of the bond(s).
#' @param target_atom_ids One or more **numeric** atom IDs (each matching
#'   \code{atoms$eleno}) to connect to \code{origin_atom_id}.
#' @param bond_type Character scalar specifying the bond type code
#'   (e.g., `"1"`, `"2"`, `"ar"`, `"am"`, `"un"`). Defaults to `"un"`.
#'
#' @details
#' This function extends the molecule<U+2019>s \code{@bonds} table. It creates sequential
#' numeric \code{bond_id}s starting from \code{molecule@maximum_bond_id + 1}.
#' If multiple targets are supplied, one bond row is added per target.
#'
#' The function validates that all supplied atom IDs exist in \code{molecule@atoms$eleno}.
#' The anchor and other properties are preserved.
#'
#' @return A \code{Molecule3D} object with additional rows in \code{@bonds}.
#'
#' @examples
#' atoms <- data.frame(
#'   eleno = c(1, 2, 3),
#'   elena = c("C", "O", "H"),
#'   x = c(0, 1.2, -0.8),
#'   y = c(0, 0, 0.5),
#'   z = c(0, 0, 0)
#' )
#' bonds <- data.frame(
#'   bond_id = 1,
#'   origin_atom_id = 1,
#'   target_atom_id = 2,
#'   bond_type = "1"
#' )
#' m <- Molecule3D(name = "COH", atoms = atoms, bonds = bonds)
#' m <- add_bonds(m, origin_atom_id = 1, target_atom_ids = c(3), bond_type = "1")
#' m@bonds
#'
#' @seealso [structures::valid_bond_types()], [combine_molecules()]
#' @export
add_bonds <- function(molecule, origin_atom_id, target_atom_ids, bond_type = "un") {
  assertions::assert_class(molecule, "structures::Molecule3D")
  max_bond_id <- molecule@maximum_bond_id
  origin_atom_id <- as.numeric(origin_atom_id)
  target_atom_ids <- as.numeric(target_atom_ids)

  df_additional_bonds <- data.frame(
    bond_id = seq_along(target_atom_ids) + max_bond_id,
    origin_atom_id = rep(origin_atom_id, times = length(target_atom_ids)),
    target_atom_id = target_atom_ids,
    bond_type = bond_type
  )

  molecule@bonds <- dplyr::bind_rows(molecule@bonds, df_additional_bonds)
  return(molecule)
}


#' Add a dummy atom defined by internal coordinates
#'
#' Inserts a <U+201C>dummy<U+201D> atom (default element label \code{"Du"}) into a
#' \code{Molecule3D} using three reference atoms and internal coordinates:
#' bond length to atom C, bond angle at B<U+2013>C<U+2013>D, and torsion A<U+2013>B<U+2013>C<U+2013>D. The new
#' atom is appended to the \code{@atoms} table and connected by a single bond
#' to atom C (via \code{\link{add_bonds}}). Identifiers are assigned using
#' \code{molecule@maximum_atom_id + 1}.
#'
#' Geometry is computed with \pkg{compas} (\code{compas::calCo()}), with angles
#' interpreted in degrees.
#'
#' @param molecule A \code{Molecule3D} object.
#' @param atom_id_a,atom_id_b,atom_id_c Numeric atom IDs (matching \code{atoms$eleno})
#'   that define the reference frame for the new atom D. The new atom is placed at
#'   a distance \code{bond_length} from \code{atom_id_c}, with bond angle
#'   \code{B<U+2013>C<U+2013>D = bond_angle} and torsion \code{A<U+2013>B<U+2013>C<U+2013>D = torsion_angle}.
#' @param bond_length Numeric scalar; distance (<U+00C5>) from atom C to the new atom D.
#' @param torsion_angle Numeric scalar; dihedral angle \code{A<U+2013>B<U+2013>C<U+2013>D} in degrees.
#' @param bond_angle Numeric scalar; bond angle \code{B<U+2013>C<U+2013>D} in degrees.
#' @param bond_type Character scalar specifying the SYBYL/Tripos bond code for the
#'   new C<U+2013>D bond (e.g., \code{"1"}, \code{"2"}, \code{"ar"}, \code{"un"}). Defaults to \code{"du"}.
#' @param elena Character scalar atom label for the dummy atom (default \code{"Du"}).
#'
#' @details
#' The function:
#' \itemize{
#'   \item Retrieves the Cartesian coordinates of atoms A, B, C via
#'         \code{\link{fetch_atom_position}}.
#'   \item Calls \code{compas::calCo()} to compute the D coordinates from
#'         \code{bond_length}, \code{bond_angle}, and \code{torsion_angle}
#'         (angles in degrees).
#'   \item Creates a one-row atoms table for D (\code{eleno = maximum\_atom\_id + 1},
#'         \code{elena = elena}, and computed \code{x,y,z}), then merges it into
#'         the input molecule via \code{\link{combine_molecules}} (with
#'         \code{update_ids = FALSE}).
#'   \item Adds a single bond between C (\code{atom_id_c}) and D using
#'         \code{\link{add_bonds}} with the supplied \code{bond_type}.
#' }
#'
#' The molecule<U+2019>s \code{@anchor} and other properties are preserved.
#'
#' @return A \code{Molecule3D} object containing the appended dummy atom and its bond to C.
#'
#' @examples
#' m <- read_mol2(system.file("benzene.mol2", package = "structures"))
#'
#' # Suppose atoms 1-2-3 exist and define a sensible frame.
#' m2 <- add_dummy_atom(
#'   molecule = m,
#'   atom_id_a = 1, atom_id_b = 2, atom_id_c = 3,
#'   bond_length = 1.5,
#'   torsion_angle = 60,
#'   bond_angle = 109.5,
#'   bond_type = "du",
#'   elena = "Du"
#' )
#' tail(m2@atoms)
#' tail(m2@bonds)
#'
#' @seealso \code{\link{add_bonds}}, \code{\link{combine_molecules}},
#'   \code{\link{fetch_atom_position}}, \code{\link{Molecule3D}}
#' @export
add_dummy_atom <- function(molecule, atom_id_a, atom_id_b, atom_id_c, bond_length, torsion_angle, bond_angle, bond_type = "du", elena = "Du") {
  assertions::assert_class(molecule, "structures::Molecule3D")

  # positions is a 3x3 (or 3xN) matrix of A,B,C coordinates (rows = eleno, cols x/y/z)
  positions <- fetch_atom_position(molecule, eleno = c(atom_id_a, atom_id_b, atom_id_c))

  # compas expects angles in degrees
  dummy_position <- compas::calCo(
    prev_atoms = positions,
    length = bond_length,
    bAngle = bond_angle,
    tAngle = torsion_angle
  )

  dummy_id <- molecule@maximum_atom_id + 1

  additional_atoms <- data.frame(
    eleno = dummy_id,
    elena = elena,
    x = dummy_position[1],
    y = dummy_position[2],
    z = dummy_position[3]
  )

  dummy_molecule <- Molecule3D(name = "Dummy", atoms = additional_atoms)
  combined <- combine_molecules(molecule, dummy_molecule, update_ids = FALSE)
  combined <- add_bonds(
    molecule = combined,
    origin_atom_id = dummy_id,
    target_atom_ids = atom_id_c,
    bond_type = bond_type
  )

  return(combined)
}
