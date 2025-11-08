test_that("SymAxis constructs with integer Cn", {
  ax <- SymAxis(Cn = 3L, posA = c(0,0,0), posB = c(0,0,1))
  expect_s3_class(ax, "structures::SymAxis")
  expect_type(ax@Cn, "integer")
  expect_equal(ax@Cn, 3L)
})

test_that("SymAxis coerces integer-like numerics to integer", {
  ax1 <- SymAxis(Cn = 3, posA = c(0,0,0), posB = c(0,0,1))     # double 3
  ax2 <- SymAxis(Cn = 3.0, posA = c(0,0,0), posB = c(0,0,1))   # double 3.0
  expect_type(ax1@Cn, "integer")
  expect_type(ax2@Cn, "integer")
  expect_identical(ax1@Cn, 3L)
  expect_identical(ax2@Cn, 3L)
})

test_that("SymAxis rejects true decimals", {
  # browser()
  expect_error(
    SymAxis(Cn = 1.1, posA = c(0,0,0), posB = c(0,0,1)),
    "must be a whole number, not",
    ignore.case = TRUE
  )
  expect_error(
    SymAxis(Cn = 2.7, posA = c(0,0,0), posB = c(0,0,1)),
    "must be a whole number, not",
    ignore.case = TRUE
  )
})

test_that("Cn must be length 1, finite, and >= 1", {
  expect_error(SymAxis(Cn = integer(0), posA = c(0,0,0), posB = c(0,0,1)),
               "single value|length", ignore.case = TRUE)
  expect_error(SymAxis(Cn = c(1L,2L), posA = c(0,0,0), posB = c(0,0,1)),
               "single", ignore.case = TRUE)
  expect_error(SymAxis(Cn = NA_integer_, posA = c(0,0,0), posB = c(0,0,1)),
               "non-NA|>= 1", ignore.case = TRUE)
  expect_error(SymAxis(Cn = 0L, posA = c(0,0,0), posB = c(0,0,1)),
               ">= 1", ignore.case = TRUE)
  expect_error(SymAxis(Cn = -2L, posA = c(0,0,0), posB = c(0,0,1)),
               ">= 1", ignore.case = TRUE)

  expect_error(SymAxis(Cn = Inf, posA = c(0,0,0), posB = c(0,0,1)),
               "finite", ignore.case = TRUE)
})

test_that("posA / posB length and finiteness validated", {
  expect_error(SymAxis(Cn = 2L, posA = c(0,0), posB = c(0,0,1)),
               "posA.*length 3", ignore.case = TRUE)
  expect_error(Symaxis <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0)),
               "posB.*length 3", ignore.case = TRUE)
  expect_error(SymAxis(Cn = 2L, posA = c(NA,0,0), posB = c(0,0,1)),
               "posA.*non-NA|finite", ignore.case = TRUE)
  expect_error(SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,Inf,1)),
               "posB.*non-NA|finite", ignore.case = TRUE)
})

test_that("posA and posB must be distinct", {
  expect_error(
    SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,0)),
    "distinct points|define an axis", ignore.case = TRUE
  )
})

test_that("print method includes key fields", {
  ax <- SymAxis(Cn = 4L, posA = c(1,2,3), posB = c(4,5,6))
  out <- capture.output(print(ax))
  expect_true(any(grepl("Symmetry Axis", out)))
  expect_true(any(grepl("Fold Symmetry \\(Cn\\): C4", out)))
  expect_true(any(grepl("PosA: 1, 2, 3", out)))
  expect_true(any(grepl("PosB: 4, 5, 6", out)))
})

test_that("setter works when assigning after construction", {
  ax <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  # integer-like assignment
  ax@Cn <- 6.0
  expect_identical(ax@Cn, 6L)
  # reject decimal
  expect_error({ ax@Cn <- 2.5 }, "must be a whole number", ignore.case = TRUE)
})
