
# Test Combine ------------------------------------------------------------

test_that("combine_molecules reindexes IDs and keeps references consistent", {

  # Molecule 1: two atoms, one bond (1-2)
  atoms1 <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 1), y = c(0, 0), z = c(0, 0)
  )
  bonds1 <- data.frame(
    bond_id = 1,
    origin_atom_id = "1",
    target_atom_id = "2"
  )
  m1 <- Molecule3D(name = "A", atoms = atoms1, bonds = bonds1, anchor = c(1,2,3))

  # Molecule 2: two atoms, one bond (1-2) — IDs overlap with m1 on purpose
  atoms2 <- data.frame(
    eleno = c(1, 2),
    elena = c("H", "H"),
    x = c(2, 3), y = c(0, 0), z = c(0, 0)
  )
  bonds2 <- data.frame(
    bond_id = 1,
    origin_atom_id = "1",
    target_atom_id = "2"
  )
  m2 <- Molecule3D(name = "B", atoms = atoms2, bonds = bonds2, anchor = c(9,9,9))

  # Combine with reindexing (update_ids = TRUE)
  cmb <- combine_molecules(m1, m2, update_ids = TRUE)
  # browser()
  expect_true(inherits(cmb, "structures::Molecule3D"))

  # Expect atoms count sums and bonds count sums
  expect_equal(nrow(cmb@atoms), nrow(m1@atoms) + nrow(m2@atoms))
  expect_equal(nrow(cmb@bonds), nrow(m1@bonds) + nrow(m2@bonds))

  # Expect 'source' columns present and correctly populated
  expect_true("source" %in% names(cmb@atoms))
  expect_true("source" %in% names(cmb@bonds))
  expect_setequal(unique(cmb@atoms$source), c("A", "B"))
  expect_setequal(unique(cmb@bonds$source), c("A", "B"))

  # m2 atoms should be reindexed by max eleno in m1 (which is 2)
  offset_atom <- m1@maximum_atom_id
  atoms2_expected <- atoms2
  atoms2_expected$eleno <- atoms2_expected$eleno + offset_atom

  # Check that all m2 (source == "B") eleno were actually shifted
  eleno_B <- cmb@atoms$eleno[cmb@atoms$source == "B"]
  expect_true(all(sort(eleno_B) == sort(atoms2_expected$eleno)))

  # m2 bond IDs should be reindexed by max bond_id in m1 (which is 1)
  offset_bond <- m1@maximum_bond_id
  bonds2_expected <- bonds2
  bonds2_expected$bond_id <- bonds2_expected$bond_id + offset_bond

  bond_id_B <- cmb@bonds$bond_id[cmb@bonds$source == "B"]
  expect_true(all(sort(bond_id_B) == sort(bonds2_expected$bond_id)))

  # --- CRITICAL: bond endpoints for the reindexed m2 must reference shifted eleno
  # (origin_atom_id/target_atom_id are character; compare against atoms$eleno as character)
  atoms_eleno_chr <- as.character(cmb@atoms$eleno)
  endpoints_B <- unique(c(cmb@bonds$origin_atom_id[cmb@bonds$source == "B"],
                          cmb@bonds$target_atom_id[cmb@bonds$source == "B"]))
  expect_true(all(endpoints_B %in% atoms_eleno_chr))

  # Anchor should be preserved from m1 (the base object)
  expect_equal(cmb@anchor, m1@anchor)
  # Name should also remain from m1
  expect_equal(cmb@name, m1@name)
})

test_that("combine_molecules without reindexing can leave duplicates", {
  # Build small overlapping molecules
  atoms1 <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 1), y = c(0, 0), z = c(0, 0)
  )
  bonds1 <- data.frame(
    bond_id = 1,
    origin_atom_id = "1",
    target_atom_id = "2"
  )
  m1 <- Molecule3D(name = "A", atoms = atoms1, bonds = bonds1)

  atoms2 <- data.frame(
    eleno = c(1, 2),
    elena = c("H", "H"),
    x = c(2, 3), y = c(0, 0), z = c(0, 0)
  )
  bonds2 <- data.frame(
    bond_id = 1,
    origin_atom_id = "1",
    target_atom_id = "2"
  )
  m2 <- Molecule3D(name = "B", atoms = atoms2, bonds = bonds2)

  # No reindexing leads to error because of duplicated eleno eleno present (since both had 1,2)
  expect_error(combine_molecules(m1, m2, update_ids = FALSE), regexp = "duplicates")

})

test_that("combine_molecules preserves non-tabular columns and adds sources", {
  atoms1 <- data.frame(
    eleno = c(10),
    elena = c("C"),
    x = 0, y = 0, z = 0
  )
  bonds1 <- data.frame(
    bond_id = 10,
    origin_atom_id = "10",
    target_atom_id = "10"
  )
  m1 <- Molecule3D(name = "A", atoms = atoms1, bonds = bonds1, misc = list(foo = 1), anchor = c(1,1,1))

  atoms2 <- data.frame(
    eleno = c(1),
    elena = c("H"),
    x = 1, y = 0, z = 0
  )
  bonds2 <- data.frame(
    bond_id = 1,
    origin_atom_id = "1",
    target_atom_id = "1"
  )
  m2 <- Molecule3D(name = "B", atoms = atoms2, bonds = bonds2, misc = list(bar = 2), anchor = c(9,9,9))

  cmb <- combine_molecules(m1, m2, update_ids = TRUE)

  # 'source' exists and matches names
  expect_true(all(c("source") %in% names(cmb@atoms)))
  expect_true(all(c("source") %in% names(cmb@bonds)))
  expect_setequal(unique(cmb@atoms$source), c("A", "B"))
  expect_setequal(unique(cmb@bonds$source), c("A", "B"))

  # anchor and name preserved from m1
  expect_equal(cmb@anchor, c(1,1,1))
  expect_equal(cmb@name, "A")

  # maximum_* derived from combined tables
  expect_equal(cmb@maximum_atom_id, max(cmb@atoms$eleno))
  expect_equal(cmb@maximum_bond_id, max(cmb@bonds$bond_id))
})

test_that("combine_molecules keeps bonds aligned with atoms after reindex", {
  # Construct molecules where m2 references its atoms (1-2),
  # after reindex they should reference (1-2 + offset)
  atoms1 <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("C", "O", "N"),
    x = c(0, 1, 2), y = 0, z = 0
  )
  bonds1 <- data.frame(
    bond_id = c(1, 2),
    origin_atom_id = c("1", "2"),
    target_atom_id = c("2", "3")
  )
  m1 <- Molecule3D(name = "A", atoms = atoms1, bonds = bonds1)

  atoms2 <- data.frame(
    eleno = c(1, 2),
    elena = c("H", "H"),
    x = c(3, 4), y = 0, z = 0
  )
  bonds2 <- data.frame(
    bond_id = c(1),
    origin_atom_id = c("1"),
    target_atom_id = c("2")
  )
  m2 <- Molecule3D(name = "B", atoms = atoms2, bonds = bonds2)

  cmb <- combine_molecules(m1, m2, update_ids = TRUE)

  # For the rows from B, each endpoint must exist in atoms$eleno (after shift)
  eleno_chr <- as.character(cmb@atoms$eleno)
  b_rows <- cmb@bonds$source == "B"
  expect_true(all(cmb@bonds$origin_atom_id[b_rows] %in% eleno_chr))
  expect_true(all(cmb@bonds$target_atom_id[b_rows] %in% eleno_chr))
})



# Test Center -------------------------------------------------------------
test_that("center returns column means for x/y/z", {
  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("C", "O", "H"),
    x = c(0, 2, 4),
    y = c(-1, 1, 3),
    z = c(10, 12, 14)
  )
  bonds <- minimal_bonds()
  mol <- Molecule3D(name = "toy", atoms = atoms, bonds = bonds)

  expect_type(mol@center, "double")
  expect_length(mol@center, 3)
  expect_named(mol@center, c("x", "y", "z"))

  exp <- c(
    x = mean(atoms$x),
    y = mean(atoms$y),
    z = mean(atoms$z)
  )
  expect_equal(mol@center, exp, tolerance = 1e-12)
})

test_that("center is invariant to atom row order", {
  atoms <- data.frame(
    eleno = c(10, 2, 7, 5),
    elena = c("C", "O", "H", "N"),
    x = c(1, 5, 9, -3),
    y = c(0, 2, 4,  6),
    z = c(3, 3, 3,  3)
  )
  bonds <- minimal_bonds()

  set.seed(1)
  idx <- sample(seq_len(nrow(atoms)))
  mol1 <- Molecule3D("m1", atoms = atoms,         bonds = bonds)
  mol2 <- Molecule3D("m2", atoms = atoms[idx, ],  bonds = bonds)

  expect_equal(mol1@center, mol2@center, tolerance = 1e-12)
})

test_that("center does not depend on bonds table", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 2),
    y = c(0, 2),
    z = c(0, 2)
  )
  bonds_empty <- minimal_bonds()
  bonds_one   <- data.frame(
    bond_id = 1,
    origin_atom_id = "1",
    target_atom_id = "2",
    bond_type = "1"
  )

  mol_a <- Molecule3D("no-bonds", atoms = atoms, bonds = bonds_empty)
  mol_b <- Molecule3D("one-bond", atoms = atoms, bonds = bonds_one)

  expect_equal(mol_a@center, mol_b@center, tolerance = 1e-12)
})

test_that("changing anchor alone does not change center", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(10, 12),
    y = c(-2,  2),
    z = c( 5,  7)
  )
  mol <- Molecule3D("anchor-test", atoms = atoms, bonds = minimal_bonds(), anchor = c(0,0,0))
  c0 <- mol@center

  mol <- set_anchor_by_position(mol, c(100, 200, 300))
  expect_equal(mol@center, c0, tolerance = 1e-12)

  mol <- set_anchor_by_atom(mol, 2)
  expect_equal(mol@center, c0, tolerance = 1e-12)
})

test_that("center translates by the same vector as the molecule", {
  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("C", "O", "H"),
    x = c(0, 1, 2),
    y = c(0, 1, 2),
    z = c(0, 1, 2)
  )
  mol <- Molecule3D("shift", atoms = atoms, bonds = minimal_bonds(), anchor = c(0,0,0))
  c0 <- mol@center

  v <- c(5, -3, 10)
  mol2 <- translate_molecule_by_vector(mol, v)
  expect_equal(mol2@center, c0 + v, tolerance = 1e-12)

  # Using anchor-targeted translation yields the same shift of the center
  target <- c(100, 200, -50)
  mol3 <- set_anchor_by_position(mol, c(1, 1, 1))
  delta <- target - mol3@anchor
  mol3 <- translate_molecule_to_position(mol3, target)
  expect_equal(mol3@center, c0 + delta, tolerance = 1e-12)
})
