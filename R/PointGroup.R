# pointgroups_bondy <- function(){
#   utils::read.csv(system.file("bondy_pointgroups.csv", path = "structures"), header = TRUE, sep = ",", row.names = FALSE)
# }


#' PointGroup3D: A 3D Point Group
#'
#' In geometry, a point group is a list of Symmetry Elements (see [SymmetryElement()]) about
#' which symmetry operations can be performed without changing the position of atoms.
#'
#' @param name String (length 1). What is the name of this point group. We do not restrict this since with no lattice restriction there are infinitely many point groups.
#' However, if you are creating one of the 32 Crystallographic 3D point groups (the ones which allow for periodic crystals) we recommend setting this based on the Schönflies notation. See [crystallographic_pointgroup_names()].
#' @param symmetry_elements A list of [SymmetryElement()] objects.
#' @export
PointGroup3D <- S7::new_class(
  name = "PointGroup",
  properties = list(
    name = S7::class_character,
    symmetry_elements = S7::new_property(
      class = S7::class_list,
      validator = function(value){
        for (element in symmetry_elements){
          if(!inherits(element, "SymmetryElement")) return(sprintf("All @symmetry_elements must inherit from: structures::SymmetryElement"))
        }
        return(NULL)
      }
    )
    # Add computed properties to pull out all of a specific type of symmetry axis (
  )
)

# Symmetry Elements -------------------------------------------------------
SymmetryElement <- S7::new_class(
  name = "SymmetryElement"
  # No properties. Just useful as parent to al
)

# Create an Class for each symmetry element that inherits from 'SymmetryElement'


# MirrorPlane -------------------------------------------------------------
MirrorPlane  <-  S7::new_class(
  parent = SymmetryElement,
  name = "MirrorPlane",
  properties = list(
    normal = S7::new_property(
      class = S7::class_numeric,
      validator = function(value){
        tol=1e-18
        if(magnitude(value) <= tol) return(sprintf("Mirror plane @normal must be a vector with magnitude > 0, not [%d]", magnitude(value)))
        if(length(value) != 3) return(sprintf("Mirror plane @normal must be a numeric vector with length 3: Not [%s]", length(value)))
        if(length(value) != 3) return(sprintf("Mirror plane @normal must be a numeric vector with length 3: Not [%s]", length(value)))
        return(NULL)
      }
    ),
    position = S7::new_property(
      class = S7::class_numeric,
      validator = function(value){
        if(length(value) != 3) return(sprintf("Mirror plane @position must be a numeric vector with length 3: Not [%s]", length(value)))
        return(NULL)
      }
    )
  )
)


# CentreOfInversion -------------------------------------------------------
CentreOfInversion <- S7::new_class(
  parent = SymmetryElement,
  name = "CentreOfInversion",
  properties = list(
    position = S7::new_property(
      class = S7::class_numeric,
      validator = function(value){
        if(length(value) != 3) return(sprintf("Center of Inversion @position must be a numeric vector with length 3: Not [%s]", length(value)))
        return(NULL)
      }
    )
  )
)


# ImproperRotationAxis ----------------------------------------------------
ImproperRotationAxis <- S7::new_class(
  parent = SymmetryElement,
  name = "ImproperRotationAxis",
  ## Flesh out: should closely resemble the properties of the Proper Rotation Axes
)


#' Schönflies notation of the 32 crystallographic point groups
#'
#' Returns the 32 crystallographic 3D point groups in Schönflies notation,
#' either as a named character vector (default) or grouped into a named list
#' by crystal system.
#'
#' @param as_list Logical. If `FALSE` (default), returns a single named vector
#'   with repeated names indicating crystal system. If `TRUE`, returns a list
#'   where each element is a vector of point groups in that system.
#'
#' @return A named character vector (default) or a named list of character vectors.
#' @export
#'
#' @examples
#' crystallographic_pointgroup_names()
crystallographic_pointgroup_names <- function(as_list=FALSE){
  point_groups_3d <- c(
    # Triclinic (2)
    triclinic = "C1",
    triclinic = "Ci",

    # Monoclinic (3)
    monoclinic = "C2",
    monoclinic = "Cs",
    monoclinic = "C2h",

    # Orthorhombic (3)
    orthorhombic = "D2",
    orthorhombic = "C2v",
    orthorhombic = "D2h",

    # Tetragonal (7)
    tetragonal = "C4",
    tetragonal = "S4",
    tetragonal = "C4h",
    tetragonal = "D4",
    tetragonal = "C4v",
    tetragonal = "D2d",
    tetragonal = "D4h",

    # Trigonal (5)
    trigonal = "C3",
    trigonal = "C3i",
    trigonal = "D3",
    trigonal = "C3v",
    trigonal = "D3d",

    # Hexagonal (7)
    hexagonal = "C6",
    hexagonal = "C3h",
    hexagonal = "C6h",
    hexagonal = "D6",
    hexagonal = "C6v",
    hexagonal = "D3h",
    hexagonal = "D6h",

    # Cubic (5)
    cubic = "T",
    cubic = "Th",
    cubic = "O",
    cubic = "Td",
    cubic = "Oh"
  )

  if(as_list){
    point_groups_3d <- split(point_groups_3d, names(point_groups_3d))
    point_groups_3d <- lapply(point_groups_3d, unname)
  }
  return(point_groups_3d)
}


# valid_symmetry_elements <- function(){
#   c(
#     "Cn" =  "Proper Rotation Axes",
#     "Sn" =  "Improper Rotation Axes",
#     "sigma" = "Mirror Plane",
#     "i" =  "Centre of Inversion",
#   )
# }




