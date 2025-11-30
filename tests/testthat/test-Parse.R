# Test write_mol2 / read_mol2 round-trip ----------------------------------

test_that("write_mol2 writes a valid minimal mol2 file with correct header", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0.1234, 1.2345),
    y = c(0, 0),
    z = c(0, 0),
    atom_type = c("C.3", "O.2")
  )
  bonds <- data.frame(
    bond_id = 1,
    origin_atom_id = 1,
    target_atom_id = 2,
    bond_type = "1",
    stringsAsFactors = FALSE
  )

  m <- Molecule3D(
    name   = "WriteTest",
    atoms  = atoms,
    bonds  = bonds,
    mol_type    = "BIOPOLYMER",
    charge_type = "MMFF94_CHARGES",
    anchor = c(0,0,0)
  )

  path <- tempfile(fileext = ".mol2")
  on.exit(unlink(path), add = TRUE)

  write_mol2(m, path)

  expect_true(file.exists(path))

  lines <- readLines(path)
  # Basic structure
  expect_true(any(grepl("^@<TRIPOS>MOLECULE", lines)))
  expect_true(any(grepl("^@<TRIPOS>ATOM",     lines)))
  expect_true(any(grepl("^@<TRIPOS>BOND",     lines)))

  idx_mol <- which(grepl("^@<TRIPOS>MOLECULE", lines))
  # Immediately following lines: name, counts, mol_type, charge_type
  expect_equal(lines[idx_mol + 1], "WriteTest")
  expect_match(lines[idx_mol + 2], "\\s*2\\s+1\\s*$")  # 2 atoms, 1 bond
  expect_equal(lines[idx_mol + 3], "BIOPOLYMER")
  expect_equal(lines[idx_mol + 4], "MMFF94_CHARGES")
})

test_that("write_mol2 uses non-scientific numeric formatting for coordinates", {
  atoms <- data.frame(
    eleno = c(1),
    elena = c("C"),
    x = 1e-6,          # would normally be in scientific format
    y = 1234567.89,    # also often formatted scientific
    z = 3.14159265,
    atom_type = "C.3"
  )
  m <- Molecule3D("FmtTest", atoms = atoms, bonds = minimal_bonds())

  path <- tempfile(fileext = ".mol2")
  on.exit(unlink(path), add = TRUE)

  write_mol2(m, path)
  lines <- readLines(path)

  idx_atom <- which(grepl("^@<TRIPOS>ATOM", lines))
  atom_lines <- lines[(idx_atom + 1):(idx_atom + nrow(m@atoms))]

  # Traverse each atom line and ensure no 'e' or 'E' notation is present
  expect_false(any(grepl("[eE]", atom_lines)))
})

test_that("write_mol2 and read_mol2 round-trip atoms, bonds, and types", {
  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("C","O","H"),
    x = c(0.123456, 1.234567, -2.345678),
    y = c(0, 1, 0),
    z = c(0, 0, 1),
    atom_type = c("C.3","O.2","H")
  )
  bonds <- data.frame(
    bond_id = c(1, 2),
    origin_atom_id = c(1, 1),
    target_atom_id = c(2, 3),
    bond_type = c("1", "1"),
    stringsAsFactors = FALSE
  )

  m <- Molecule3D(
    name   = "RoundTrip",
    atoms  = atoms,
    bonds  = bonds,
    mol_type    = "PROTEIN",
    charge_type = "GAUSS80_CHARGES"
  )

  path <- tempfile(fileext = ".mol2")
  on.exit(unlink(path), add = TRUE)

  write_mol2(m, path)
  m2 <- read_mol2(path)

  # Same counts
  expect_equal(nrow(m2@atoms), nrow(m@atoms))
  expect_equal(nrow(m2@bonds), nrow(m@bonds))

  # Name from filename in read_mol2() (since name inferred from path), so we check counts and coordinates
  expect_equal(round(m2@atoms$x, 4), round(m@atoms$x, 4))
  expect_equal(round(m2@atoms$y, 4), round(m@atoms$y, 4))
  expect_equal(round(m2@atoms$z, 4), round(m@atoms$z, 4))

  # Bond connectivity preserved
  expect_setequal(m2@bonds$origin_atom_id, m@bonds$origin_atom_id)
  expect_setequal(m2@bonds$target_atom_id, m@bonds$target_atom_id)

  # mol_type and charge_type come back from the MOLECULE section
  expect_equal(m2@mol_type, "PROTEIN")
  expect_equal(m2@charge_type, "GAUSS80_CHARGES")
})

test_that("read_mol2 correctly parses mol_type and charge_type from TRIPOS<MOLECULE> section", {
  # Construct a tiny mol2 file as pure text
  tmp <- tempfile(fileext = ".mol2")
  on.exit(unlink(tmp), add = TRUE)

  lines <- c(
    "#\tName: Custom",
    "#\tCreation Time: 2024-01-01 00:00:00",
    "@<TRIPOS>MOLECULE",
    "CustomName",
    " 2 1 0 0 0",
    "NUCLEIC_ACID",
    "MULLIKEN_CHARGES",
    "",
    "@<TRIPOS>ATOM",
    "      1 C1        0.0000  0.0000  0.0000 C.3  1 RES1  0.0000",
    "      2 O1        1.0000  0.0000  0.0000 O.2  1 RES1  0.0000",
    "@<TRIPOS>BOND",
    "     1    1    2 1"
  )

  writeLines(lines, tmp)

  m <- read_mol2(tmp)

  expect_equal(m@mol_type, "NUCLEIC_ACID")
  expect_equal(m@charge_type, "MULLIKEN_CHARGES")
})
