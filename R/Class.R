#' Create Structure Object
#'
#' @param atoms data.frame of atoms with columns:  eleno (element number) and elena (element name)
#' @param bonds data.frame of bonds (1 row per bond) with columns: from & to where values represent eleno of the atoms they are connecting. Direction is meaningless
#'
#' @returns an object of class 'Structure'
#' @export
#'
Structure <- S7::new_class(
  name = "Structure",
  properties = list(
    atoms = S7::class_data.frame,
    bonds = S7::class_data.frame
  ),
  validator = function(self) {
    # Validate Bond Columns
    required_bond_cols <- c("from", "to")
    observed_bond_cols <- colnames(self@bonds)
    if (!any(required_bond_cols %in% observed_bond_cols)) {
      missing <- setdiff(required_bond_cols, observed_bond_cols)
      toString("@bonds data.frame is missing required column/s: ", missing)
    }

    # Validate Atom Columns
    required_atom_cols <- c("eleno", "elena", "x", "y", "z")
    observed_atom_cols <- colnames(self@atoms)
    if (!any(required_atom_cols %in% observed_atom_cols)) {
      missing <- setdiff(required_atom_cols, observed_atom_cols)
      toString("@atoms data.frame is missing required column/s: ", missing)
    }
  }
)
