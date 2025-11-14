
# ProperRotationAxis Class Tests -----------------------------------------------------


test_that("ProperRotationAxis constructs with integer Cn", {
  ax <- ProperRotationAxis(Cn = 3L, posA = c(0,0,0), posB = c(0,0,1))
  expect_s3_class(ax, "structures::ProperRotationAxis")
  expect_type(ax@Cn, "integer")
  expect_equal(ax@Cn, 3L)
})

test_that("ProperRotationAxis coerces integer-like numerics to integer", {
  ax1 <- ProperRotationAxis(Cn = 3, posA = c(0,0,0), posB = c(0,0,1))     # double 3
  ax2 <- ProperRotationAxis(Cn = 3.0, posA = c(0,0,0), posB = c(0,0,1))   # double 3.0
  expect_type(ax1@Cn, "integer")
  expect_type(ax2@Cn, "integer")
  expect_identical(ax1@Cn, 3L)
  expect_identical(ax2@Cn, 3L)
})

test_that("ProperRotationAxis rejects true decimals", {
  # browser()
  expect_error(
    ProperRotationAxis(Cn = 1.1, posA = c(0,0,0), posB = c(0,0,1)),
    "must be a whole number, not",
    ignore.case = TRUE
  )
  expect_error(
    ProperRotationAxis(Cn = 2.7, posA = c(0,0,0), posB = c(0,0,1)),
    "must be a whole number, not",
    ignore.case = TRUE
  )
})

test_that("Cn must be length 1, finite, and >= 1", {
  expect_error(ProperRotationAxis(Cn = integer(0), posA = c(0,0,0), posB = c(0,0,1)),
               "single value|length", ignore.case = TRUE)
  expect_error(ProperRotationAxis(Cn = c(1L,2L), posA = c(0,0,0), posB = c(0,0,1)),
               "single", ignore.case = TRUE)
  expect_error(ProperRotationAxis(Cn = NA_integer_, posA = c(0,0,0), posB = c(0,0,1)),
               "non-NA|>= 1", ignore.case = TRUE)
  expect_error(ProperRotationAxis(Cn = 0L, posA = c(0,0,0), posB = c(0,0,1)),
               ">= 1", ignore.case = TRUE)
  expect_error(ProperRotationAxis(Cn = -2L, posA = c(0,0,0), posB = c(0,0,1)),
               ">= 1", ignore.case = TRUE)

  expect_error(ProperRotationAxis(Cn = Inf, posA = c(0,0,0), posB = c(0,0,1)),
               "finite", ignore.case = TRUE)
})

test_that("posA / posB length and finiteness validated", {
  expect_error(ProperRotationAxis(Cn = 2L, posA = c(0,0), posB = c(0,0,1)),
               "posA.*length 3", ignore.case = TRUE)

  expect_error(ProperRotationAxis <- ProperRotationAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0)),
               "posB.*length 3", ignore.case = TRUE)
  expect_error(ProperRotationAxis(Cn = 2L, posA = c(NA,0,0), posB = c(0,0,1)),
               "must have no missing values", ignore.case = TRUE)
  expect_error(ProperRotationAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,Inf,1)),
               "must have no infinite values", ignore.case = TRUE)
})

test_that("posA and posB must be distinct", {
  expect_error(
    ProperRotationAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,0)),
    "must be distinct", ignore.case = TRUE
  )
})


test_that("setter works when assigning after construction", {
  ax <- ProperRotationAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  # integer-like assignment
  ax@Cn <- 6.0
  expect_identical(ax@Cn, 6L)
  # reject decimal
  expect_error({ ax@Cn <- 2.5 }, "must be a whole number", ignore.case = TRUE)
})


# tests/test-ProperRotationAxis-label.R

test_that("ProperRotationAxis label: defaults to a single, non-empty character", {
  ax <- ProperRotationAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))  # no label supplied

  expect_type(ax@label, "character")
  expect_length(ax@label, 1L)
  expect_false(is.na(ax@label))
  expect_true(nzchar(ax@label))      # don't assert an exact default string
})

test_that("ProperRotationAxis label: accepts and stores a custom label", {
  ax <- ProperRotationAxis(Cn = 3L, posA = c(0,0,0), posB = c(0,0,1), label = "principal_z")
  expect_identical(ax@label, "principal_z")
})

test_that("ProperRotationAxis label: can be updated via set_props and revalidated", {
  ax  <- ProperRotationAxis(Cn = 3L, posA = c(0,0,0), posB = c(0,0,1), label = "z1")
  ax2 <- S7::set_props(ax, label = "z2")
  expect_identical(ax2@label, "z2")
})

test_that("ProperRotationAxis label: rejects invalid values", {
  # empty string
  expect_error(
    ProperRotationAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1), label = ""),
    regexp = "empty string|must NOT be an empty string",
    ignore.case = TRUE
  )

  # NA
  expect_error(
    ProperRotationAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1), label = NA_character_),
    regexp = "must NOT be a missing|NA",
    ignore.case = TRUE
  )

  # length > 1
  expect_error(
    ProperRotationAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1), label = c("a","b")),
    regexp = "length|vector of length",
    ignore.case = TRUE
  )

  # non-character (class check from S7::class_character)
  expect_error(
    ProperRotationAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1), label = 123),
    regexp = "character|class",
    ignore.case = TRUE
  )
})

test_that("ProperRotationAxis label: preserved by transform_symmetry_axis", {
  ax <- ProperRotationAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1), label = "axis_A")

  translate <- function(p, dx=0, dy=0, dz=0) c(x = p["x"]+dx, y = p["y"]+dy, z = p["z"]+dz)
  ax_t <- transform_symmetry_axis(ax, translate, dx = 1, dy = -2, dz = 0.5)

  expect_identical(ax_t@label, "axis_A")
  # sanity check endpoints changed as expected (not strictly about label, but cheap to assert)
  expect_equal(ax_t@posA, c(x=1, y=-2, z=0.5))
  expect_equal(ax_t@posB, c(x=1, y=-2, z=1.5))
})



# Test Transformations Work Sensibly -----------------------------------------------------
test_that("translation transform updates endpoints and preserves Cn", {
  ax <- ProperRotationAxis(Cn = 3L, posA = c(0,0,0), posB = c(0,0,1))

  translate <- function(p, dx=0, dy=0, dz=0) {
    unname(c(p["x"] + dx, p["y"] + dy, p["z"] + dz))
  }

  ax2 <- transform_symmetry_axis(ax, translate, dx = 1, dy = 2, dz = 3)

  expect_s3_class(ax2, "structures::ProperRotationAxis")
  expect_identical(ax2@Cn, 3L)                 # Cn unchanged
  expect_equal(unname(ax2@posA), c(1,2,3), tolerance = 1e-12)
  expect_equal(unname(ax2@posB), c(1,2,4), tolerance = 1e-12)
})

test_that("rotation around Z by 90 deg moves X-axis endpoints as expected", {
  # Start with an axis along +X: from (1,0,0) to (2,0,0)
  ax <- ProperRotationAxis(Cn = 2L, posA = c(1,0,0), posB = c(2,0,0))

  rotate_z <- function(p, theta) {
    c(x =  cos(theta)*p["x"] - sin(theta)*p["y"],
      y =  sin(theta)*p["x"] + cos(theta)*p["y"],
      z =  p["z"])
  }

  axr <- transform_symmetry_axis(ax, rotate_z, theta = pi/2)

  expect_identical(axr@Cn, 2L)
  expect_equal(unname(axr@posA), c(0,1,0), tolerance = 1e-12)
  expect_equal(unname(axr@posB), c(0,2,0), tolerance = 1e-12)
})

test_that("transformation returning a list is accepted", {
  ax <- ProperRotationAxis(Cn = 4L, posA = c(0,0,0), posB = c(1,1,1))

  as_list <- function(p) list(x = p["x"] + 1, y = p["y"], z = p["z"])
  ax2 <- transform_symmetry_axis(ax, as_list)

  expect_identical(ax2@Cn, 4L)
  expect_equal(unname(ax2@posA), c(1,0,0), tolerance = 1e-12)
  expect_equal(unname(ax2@posB), c(2,1,1), tolerance = 1e-12)
})

test_that("bad transformation outputs error: unnamed, missing fields, or non-finite", {
  ax <- ProperRotationAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))

  # Missing 'z'
  missing_z <- function(p) c(x = p["x"], y = p["y"])
  expect_error(transform_symmetry_axis(ax, missing_z),
               "must be a numeric vector of length 3")

  # Non-finite
  nonfinite <- function(p) c(x = NA_real_, y = p["y"], z = p["z"])
  expect_error(transform_symmetry_axis(ax, nonfinite),
               "finite")
})


# Test Print Generic ------------------------------------------------------

test_that("print.ProperRotationAxis runs without error and returns invisible ProperRotationAxis", {
  ax <- ProperRotationAxis(Cn = 3L, posA = c(0, 0, 0), posB = c(0, 0, 1))

  out <- NULL
  invisible(capture.output({
    out <- expect_no_error(print(ax))
  }))

  expect_identical(out, ax)

})


# Test as.data.frame generic ---------------------------------------------------

test_that("as.data.frame(ProperRotationAxis) returns a single-row data.frame with expected columns and types", {
  ax <- ProperRotationAxis(Cn = 3L, posA = c(0, 1, 2), posB = c(3, 4, 5), label = "axis-A")
  df <- as.data.frame(ax)

  # Structure
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 1L)
  expect_identical(
    names(df),
    c("label", "Cn", "x", "y", "z", "xend", "yend", "zend")
  )

  # Types
  expect_type(df$label, "character")
  # Cn should remain integer-like; the class may be "integer" within a data.frame
  expect_true(is.integer(df$Cn) || is.numeric(df$Cn))

  # Values
  expect_identical(df$label, "axis-A")
  expect_identical(as.integer(df$Cn), 3L)
  expect_equal(unname(c(df$x, df$y, df$z)), c(0, 1, 2))
  expect_equal(unname(c(df$xend, df$yend, df$zend)), c(3, 4, 5))
})

test_that("as.data.frame(ProperRotationAxis) respects numeric Cn input via coercion and preserves label", {
  # Pass Cn as numeric; class constructor coerces to integer
  ax <- ProperRotationAxis(Cn = 4, posA = c(1, 0, 0), posB = c(0, 1, 0), label = "axis-B")
  df <- as.data.frame(ax)

  expect_identical(df$label, "axis-B")
  expect_identical(as.integer(df$Cn), 4L)
  expect_equal(unname(c(df$x, df$y, df$z)), c(1, 0, 0))
  expect_equal(unname(c(df$xend, df$yend, df$zend)), c(0, 1, 0))
})

test_that("rbind of multiple as.data.frame(ProperRotationAxis) rows preserves labels and values", {
  ax1 <- ProperRotationAxis(Cn = 2L, posA = c(0, 0, 0), posB = c(0, 0, 1), label = "A")
  ax2 <- ProperRotationAxis(Cn = 3L, posA = c(1, 0, 0), posB = c(1, 0, 1), label = "B")

  rows <- lapply(list(ax1, ax2), as.data.frame)
  df <- do.call(rbind, rows)

  expect_equal(nrow(df), 2L)
  expect_identical(df$label, c("A", "B"))
  expect_identical(as.integer(df$Cn), c(2L, 3L))

  expect_equal(df$x,    c(0, 1))
  expect_equal(df$y,    c(0, 0))
  expect_equal(df$z,    c(0, 0))
  expect_equal(df$xend, c(0, 1))
  expect_equal(df$yend, c(0, 0))
  expect_equal(df$zend, c(1, 1))
})

