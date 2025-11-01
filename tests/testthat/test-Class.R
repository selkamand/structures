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
