#' @export
S7::method(plot, Molecule3D) <- function(x, ...) {
  message("For 3D visualisation of Molecule3D objects we highly reccomend the `chemviewR::plot_molecule()` function")
  mx = as.matrix(x)
  plot(x=mx[,"x"], y=mx[,"y"], type = "n", bty="n", asp = 1, axes = FALSE, ann=FALSE)
  with(
    x@bond_positions,
    segments(x0 = x, y0 = y, x1 = xend, y1 = yend)
  )

  with(
    x@atoms,
    points(x = x, y = y, pch=19, cex = 4)
  )

  return(invisible(x))
}
