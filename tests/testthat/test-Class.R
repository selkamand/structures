

# Test Add Bonds ----------------------------------------------------------------

test_that("add_bonds() adds a single bond and increments bond_id", {
  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("C", "O", "H"),
    x = c(0, 1.2, -0.8),
    y = c(0, 0, 0.5),
    z = c(0, 0, 0)
  )
  bonds <- data.frame(
    bond_id = 1,
    origin_atom_id = 1,
    target_atom_id = 2,
    bond_type = "1"
  )
  m <- Molecule3D("COH", atoms = atoms, bonds = bonds, anchor = c(0,0,0))

  before_anchor <- m@anchor
  m2 <- add_bonds(m, origin_atom_id = 1, target_atom_ids = 3, bond_type = "2")

  # structure
  expect_equal(nrow(m2@bonds), nrow(m@bonds) + 1)
  # new row is correct
  new_row <- m2@bonds[nrow(m2@bonds), ]
  expect_equal(new_row$origin_atom_id, 1)
  expect_equal(new_row$target_atom_id, 3)
  expect_equal(new_row$bond_type, "2")
  # id advanced from maximum_bond_id
  expect_equal(m2@maximum_bond_id, max(m@bond_ids) + 1)
  # anchor unchanged
  expect_equal(m2@anchor, before_anchor)
})

test_that("add_bonds() with multiple targets creates multiple rows with sequential IDs", {
  atoms <- data.frame(
    eleno = 1:4,
    elena = c("C","H","H","H"),
    x = c(0,1,0,-1),
    y = c(0,0,1,0),
    z = c(0,0,0,1)
  )
  # start with empty bonds to ensure IDs start at 1
  bonds <- data.frame(
    bond_id = numeric(0),
    origin_atom_id = character(0),
    target_atom_id = character(0),
    bond_type = character(0)
  )
  m <- Molecule3D("Methane-ish", atoms = atoms, bonds = bonds)

  m2 <- add_bonds(m, origin_atom_id = 1, target_atom_ids = c(2,3,4), bond_type = "1")

  expect_equal(nrow(m2@bonds), 3)
  # IDs should be 1,2,3
  expect_equal(sort(unique(m2@bonds$bond_id)), 1:3)
  # All origin IDs should be 1
  expect_true(all(m2@bonds$origin_atom_id == 1))
  # Targets should match
  expect_setequal(m2@bonds$target_atom_id, c(2,3,4))
  # bond_type recycled properly
  expect_true(all(m2@bonds$bond_type == "1"))
})

test_that("add_bonds() errors for missing origin or target IDs", {
  atoms <- data.frame(
    eleno = c(10, 20),
    elena = c("C","O"),
    x = c(0, 1),
    y = c(0, 0),
    z = c(0, 0)
  )
  bonds <- minimal_bonds()
  m <- Molecule3D("TwoAtoms", atoms = atoms, bonds = bonds)

  # origin missing
  expect_error(
    add_bonds(m, origin_atom_id = 999, target_atom_ids = 20),
    "describes atoms not present in @atoms dataframe"
  )

  # target missing
  expect_error(
    add_bonds(m, origin_atom_id = 10, target_atom_ids = c(20, 999)),
    "describes atoms not present in @atoms dataframe"
  )
})

test_that("add_bonds() respects existing maximum_bond_id even with gaps", {
  atoms <- data.frame(
    eleno = c(1,2,3),
    elena = c("C","O","H"),
    x = c(0,1,0),
    y = c(0,0,1),
    z = c(0,0,0)
  )
  # bond IDs are non-contiguous; max is 10
  bonds <- data.frame(
    bond_id = c(5, 10),
    origin_atom_id = c(1, 1),
    target_atom_id = c(2, 3),
    bond_type = c("1","1")
  )
  m <- Molecule3D("GappyIDs", atoms = atoms, bonds = bonds)

  m2 <- add_bonds(m, origin_atom_id = 2, target_atom_ids = 3, bond_type = "2")

  # new id should be 11
  expect_equal(m2@maximum_bond_id, 11)
  expect_true(any(m2@bonds$bond_id == 11))
})

test_that("add_bonds() allows vector bond_type of same length as targets", {
  atoms <- data.frame(
    eleno = 1:4,
    elena = c("C","H","H","H"),
    x = 0:3, y = 0:3, z = 0:3
  )
  m <- Molecule3D("VecType", atoms = atoms, bonds = minimal_bonds())

  types <- c("1","2","am")
  m2 <- add_bonds(m, origin_atom_id = 1, target_atom_ids = c(2,3,4), bond_type = types)

  # Expect types to appear exactly once each
  expect_setequal(m2@bonds$bond_type, types)
  expect_equal(nrow(m2@bonds), 3)
})

test_that("add_bonds() coerces ids to numeric and preserves atoms/anchor", {
  atoms <- data.frame(
    eleno = c(1,2),
    elena = c("C","O"),
    x = c(0,1), y = c(0,0), z = c(0,0)
  )
  m <- Molecule3D("Coerce", atoms = atoms, bonds = minimal_bonds(), anchor = c(9,9,9))
  before_atoms  <- m@atoms
  before_anchor <- m@anchor

  # pass ids as character; function coerces to numeric internally
  m2 <- add_bonds(m, origin_atom_id = "1", target_atom_ids = c("2"))

  expect_equal(nrow(m2@bonds), 1)
  expect_equal(m2@bonds$origin_atom_id[1], 1)
  expect_equal(m2@bonds$target_atom_id[1], 2)
  expect_identical(m2@atoms, before_atoms)
  expect_identical(m2@anchor, before_anchor)
})

test_that("add_bonds() works when starting with completely empty bonds table", {
  atoms <- data.frame(
    eleno = c(1,2,3),
    elena = c("C","O","H"),
    x = c(0,1,0), y = c(0,0,1), z = c(0,0,0)
  )
  m <- Molecule3D("EmptyBonds", atoms = atoms, bonds = minimal_bonds())

  m2 <- add_bonds(m, origin_atom_id = 1, target_atom_ids = c(2,3), bond_type = "1")

  expect_equal(nrow(m2@bonds), 2)
  expect_true(all(m2@bonds$bond_id %in% c(1,2)))
  expect_setequal(m2@bonds$target_atom_id, c(2,3))
  expect_true(all(m2@bonds$origin_atom_id == 1))
})

test_that("add_bonds() with one origin and multiple targets adds multiple bonds", {
  atoms <- data.frame(
    eleno = 1:5,
    elena = c("C","H","H","H","O"),
    x = c(0, 1, -1, 0, 0.5),
    y = c(0, 0,  0, 1, 0.5),
    z = c(0, 0,  0, 0, -0.5)
  )
  # start with an empty bonds table
  m <- Molecule3D("MultiTargets", atoms = atoms, bonds = minimal_bonds())

  origin <- 1
  targets <- c(2, 3, 4, 5)
  m2 <- add_bonds(m, origin_atom_id = origin, target_atom_ids = targets, bond_type = "1")

  # number of new rows equals number of targets
  expect_equal(nrow(m2@bonds), length(targets))
  # all origin IDs match the single origin provided
  expect_true(all(m2@bonds$origin_atom_id == origin))
  # targets are exactly those supplied (order-insensitive)
  expect_setequal(m2@bonds$target_atom_id, targets)
  # bond ids are sequential starting at 1 when bonds were empty
  expect_equal(sort(unique(m2@bonds$bond_id)), seq_len(length(targets)))
})


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

test_that("combine_molecules works when m1 has no bonds", {
  # m1: one atom, no bonds
  atoms1 <- data.frame(
    eleno = 1,
    elena = "C",
    x = 0, y = 0, z = 0
  )
  bonds1 <- minimal_bonds()
  m1 <- Molecule3D(name = "A", atoms = atoms1, bonds = bonds1, anchor = c(1,2,3))

  # m2: two atoms, one bond (1-2)
  atoms2 <- data.frame(
    eleno = c(1, 2),
    elena = c("H", "H"),
    x = c(1, 2), y = 0, z = 0
  )
  bonds2 <- data.frame(
    bond_id = 1,
    origin_atom_id = "1",
    target_atom_id = "2"
  )
  m2 <- Molecule3D(name = "B", atoms = atoms2, bonds = bonds2)

  cmb <- combine_molecules(m1, m2, update_ids = TRUE)

  # Atoms = sum; Bonds = from m2 only
  expect_equal(nrow(cmb@atoms), nrow(m1@atoms) + nrow(m2@atoms))
  expect_equal(nrow(cmb@bonds), nrow(m2@bonds))

  # Bond IDs from m2 shouldn't shift because m1@maximum_bond_id == 0
  expect_equal(unique(cmb@bonds$bond_id), m2@bonds$bond_id)

  # Sources: bonds only from B
  expect_true("source" %in% names(cmb@bonds))
  expect_setequal(unique(cmb@bonds$source), "B")

  # Anchor/name preserved from m1
  expect_equal(cmb@anchor, m1@anchor)
  expect_equal(cmb@name, m1@name)

  # All bond endpoints exist in the combined atom table
  eleno_chr <- as.character(cmb@atoms$eleno)
  expect_true(all(cmb@bonds$origin_atom_id %in% eleno_chr))
  expect_true(all(cmb@bonds$target_atom_id %in% eleno_chr))
})

test_that("combine_molecules works when m2 has no bonds", {
  # m1: two atoms, one bond (1-2)
  atoms1 <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 1), y = 0, z = 0
  )
  bonds1 <- data.frame(
    bond_id = 1,
    origin_atom_id = "1",
    target_atom_id = "2"
  )
  m1 <- Molecule3D(name = "A", atoms = atoms1, bonds = bonds1, anchor = c(4,5,6))

  # m2: one atom, no bonds
  atoms2 <- data.frame(
    eleno = 1,
    elena = "H",
    x = 2, y = 0, z = 0
  )
  bonds2 <- minimal_bonds()
  m2 <- Molecule3D(name = "B", atoms = atoms2, bonds = bonds2)

  cmb <- combine_molecules(m1, m2, update_ids = TRUE)

  # Atoms = sum; Bonds = from m1 only
  expect_equal(nrow(cmb@atoms), nrow(m1@atoms) + nrow(m2@atoms))
  expect_equal(nrow(cmb@bonds), nrow(m1@bonds))

  # Sources: bonds only from A
  expect_true("source" %in% names(cmb@bonds))
  expect_setequal(unique(cmb@bonds$source), "A")

  # m2's atom ID should be offset by m1@maximum_atom_id (which is 2)
  expect_true(any(cmb@atoms$source == "B"))
  eleno_B <- cmb@atoms$eleno[cmb@atoms$source == "B"]
  expect_setequal(eleno_B, 1 + m1@maximum_atom_id)

  # Anchor/name preserved from m1
  expect_equal(cmb@anchor, m1@anchor)
  expect_equal(cmb@name, m1@name)
})

test_that("combine_molecules works when both molecules have no bonds", {
  # m1: one atom, no bonds
  atoms1 <- data.frame(
    eleno = 10,
    elena = "C",
    x = 0, y = 0, z = 0
  )
  m1 <- Molecule3D(name = "A", atoms = atoms1, bonds = minimal_bonds(), anchor = c(1,1,1))

  # m2: one atom, no bonds
  atoms2 <- data.frame(
    eleno = 1,
    elena = "H",
    x = 1, y = 0, z = 0
  )
  m2 <- Molecule3D(name = "B", atoms = atoms2, bonds = minimal_bonds(), anchor = c(9,9,9))

  cmb <- combine_molecules(m1, m2, update_ids = TRUE)

  # Bond table should exist and be empty
  expect_true(is.data.frame(cmb@bonds))
  expect_equal(nrow(cmb@bonds), 0)

  # Columns should include the canonical set (and 'source' added by combine)
  expect_true(all(c("bond_id","origin_atom_id","target_atom_id","bond_type","source") %in% names(cmb@bonds)))

  # Atoms are combined and 'source' present
  expect_equal(nrow(cmb@atoms), 2)
  expect_true("source" %in% names(cmb@atoms))
  expect_setequal(unique(cmb@atoms$source), c("A", "B"))

  # Anchor/name preserved from m1
  expect_equal(cmb@anchor, c(1,1,1))
  expect_equal(cmb@name, "A")
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


# Test Dummy Atom Addition ------------------------------------------------
# tests/testthat/test-add_dummy_atom.R

test_that("add_dummy_atom: adds 1 atom with next ID and a single bond to C", {
  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("A", "B", "C"),
    x = c(0,  1,  1),
    y = c(0,  0,  1),
    z = c(0,  0,  0)
  )
  bonds <- data.frame(
    bond_id = numeric(0),
    origin_atom_id = numeric(0),
    target_atom_id = numeric(0),
    bond_type = character(0)
  )
  m <- Molecule3D("triad", atoms = atoms, bonds = bonds, anchor = c(0,0,0))
  max_atom <- m@maximum_atom_id
  max_bond <- m@maximum_bond_id

  m2 <- add_dummy_atom(
    molecule = m,
    atom_id_a = 1, atom_id_b = 2, atom_id_c = 3,
    bond_length = 1.5, torsion_angle = 60, bond_angle = 109.5,
    bond_type = "1", elena = "Du"
  )

  # Atom count increases by 1
  expect_equal(nrow(m2@atoms), nrow(m@atoms) + 1)

  # New atom has next ID
  new_id <- max_atom + 1
  expect_true(new_id %in% m2@atoms$eleno)

  # New atom has expected element label
  expect_equal(m2@atoms$elena[m2@atoms$eleno == new_id], "Du")

  # One new bond added, points to C with correct type and incremented bond_id
  expect_equal(nrow(m2@bonds), 1)
  expect_equal(m2@bonds$origin_atom_id[1], new_id)
  expect_equal(m2@bonds$target_atom_id[1], 3)
  expect_equal(m2@bonds$bond_type[1], "1")
  expect_equal(m2@bonds$bond_id[1], max_bond + 1)

  # Anchor preserved
  expect_equal(m2@anchor, m@anchor)
})

test_that("add_dummy_atom: works when molecule already has bonds (bond_id increments)", {
  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("A", "B", "C"),
    x = c(0,  1,  1),
    y = c(0,  0,  1),
    z = c(0,  0,  0)
  )
  bonds <- data.frame(
    bond_id = c(1L, 2L),
    origin_atom_id = c(1, 2),
    target_atom_id = c(2, 3),
    bond_type = c("1", "1")
  )
  m <- Molecule3D("triad", atoms = atoms, bonds = bonds)

  m2 <- add_dummy_atom(
    molecule = m,
    atom_id_a = 1, atom_id_b = 2, atom_id_c = 3,
    bond_length = 1.0, torsion_angle = 0, bond_angle = 90,
    bond_type = "ar", elena = "Du"
  )

  # Bond table got one extra row
  expect_equal(nrow(m2@bonds), nrow(m@bonds) + 1)

  # New bond_id increments from previous maximum
  expect_equal(max(m2@bonds$bond_id), max(m@bonds$bond_id) + 1)

  # New bond connects (new_id -> 3) and has the specified type
  new_id <- max(m2@atoms$eleno)
  expect_true(any(with(m2@bonds, origin_atom_id == new_id & target_atom_id == 3 & bond_type == "ar")))
})

test_that("add_dummy_atom: default elena can be overridden", {
  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("A", "B", "C"),
    x = c(0,  1,  1),
    y = c(0,  0,  1),
    z = c(0,  0,  0)
  )
  m <- Molecule3D("triad", atoms = atoms, bonds = minimal_bonds())

  m2 <- add_dummy_atom(
    molecule = m,
    atom_id_a = 1, atom_id_b = 2, atom_id_c = 3,
    bond_length = 1.2, torsion_angle = 30, bond_angle = 100,
    bond_type = "2", elena = "Xx"
  )

  new_id <- max(m2@atoms$eleno)
  expect_equal(m2@atoms$elena[m2@atoms$eleno == new_id], "Xx")
})

test_that("add_dummy_atom: errors when reference atoms are missing", {
  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("A", "B", "C"),
    x = c(0,  1,  1),
    y = c(0,  0,  1),
    z = c(0,  0,  0)
  )
  m <- Molecule3D("triad", atoms = atoms, bonds = minimal_bonds())

  expect_error(
    add_dummy_atom(
      molecule = m,
      atom_id_a = 1, atom_id_b = 2, atom_id_c = 999,   # <-- missing
      bond_length = 1.0, torsion_angle = 0, bond_angle = 90
    ),
    regexp = "could not be found|not.*found",
    ignore.case = TRUE
  )
})

test_that("add_dummy_atom: creates a bond when no bonds existed before (bond_id starts at 1)", {
  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("A", "B", "C"),
    x = c(0,  1,  1),
    y = c(0,  0,  1),
    z = c(0,  0,  0)
  )
  m <- Molecule3D("triad", atoms = atoms, bonds = minimal_bonds())

  m2 <- add_dummy_atom(
    molecule = m,
    atom_id_a = 1, atom_id_b = 2, atom_id_c = 3,
    bond_length = 1.0, torsion_angle = 0, bond_angle = 90
  )

  expect_equal(nrow(m2@bonds), 1)
  expect_identical(m2@bonds$bond_id, 1)  # starts at 1 when previously empty
})

test_that("add_dummy_atom: geometry matches compas::calCo (if available)", {
  testthat::skip_if_not_installed("compas")

  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("A", "B", "C"),
    x = c(0,  1,  1),
    y = c(0,  0,  1),
    z = c(0,  0,  0)
  )
  m <- Molecule3D("triad", atoms = atoms, bonds = minimal_bonds())

  # Reference positions for A, B, C
  ref_pos <- fetch_atom_position(m, c(1, 2, 3))
  # Choose simple internal coordinates
  L <- 1.5; tau <- 60; ang <- 109.5

  expected <- compas::calCo(prev_atoms = ref_pos, length = L, bAngle = ang, tAngle = tau)

  m2 <- add_dummy_atom(
    molecule = m,
    atom_id_a = 1, atom_id_b = 2, atom_id_c = 3,
    bond_length = L, torsion_angle = tau, bond_angle = ang
  )

  new_id <- max(m2@atoms$eleno)
  got <- with(m2@atoms[m2@atoms$eleno == new_id, ], c(x, y, z))

  expect_equal(got, as.numeric(expected), tolerance = 1e-7)
})

test_that("add_dummy_atom: types are correct (numeric IDs, character bond_type)", {
  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("A", "B", "C"),
    x = c(0,  1,  1),
    y = c(0,  0,  1),
    z = c(0,  0,  0)
  )
  m <- Molecule3D("triad", atoms = atoms, bonds = minimal_bonds())

  m2 <- add_dummy_atom(
    molecule = m,
    atom_id_a = 1, atom_id_b = 2, atom_id_c = 3,
    bond_length = 1.3, torsion_angle = 10, bond_angle = 100,
    bond_type = "am"
  )

  expect_type(m2@atoms$eleno, "double") # numeric in R is double
  expect_type(m2@bonds$bond_id, "double")
  expect_type(m2@bonds$origin_atom_id, "double")
  expect_type(m2@bonds$target_atom_id, "double")
  expect_type(m2@bonds$bond_type, "character")
})


# Test eleno by bond network ----------------------------------------------


# Helper to make a simple Molecule3D from an undirected edge list
.make_mol_from_edges <- function(edges, n_atoms = NULL, name = "testmol") {
  if (is.null(n_atoms)) {
    n_atoms <- max(if (length(edges)) unlist(edges) else 0)
  }
  atoms <- data.frame(
    eleno = seq_len(n_atoms),
    elena = rep("C", n_atoms),
    x = seq_len(n_atoms), y = 0, z = 0
  )
  if (is.null(edges) || nrow(edges) == 0) {
    bonds <- minimal_bonds()
  } else {
    bonds <- data.frame(
      bond_id = seq_len(nrow(edges)),
      origin_atom_id = edges[, 1],
      target_atom_id = edges[, 2],
      bond_type = "1"
    )
  }
  Molecule3D(name = name, atoms = atoms, bonds = bonds)
}

testthat::test_that("Linear chain splits correctly and honors direction", {
  # Chain 1-2-3-4; break bond between 2-3
  edges <- rbind(
    c(1, 2),
    c(2, 3),
    c(3, 4)
  )
  m <- .make_mol_from_edges(edges, n_atoms = 4)

  # Direction toward atom 2 => cluster should be {1,2}
  ids_left <- fetch_eleno_downstream_of_bond(m, bond_id = 2, direction_atom_id = 2)
  testthat::expect_setequal(ids_left, c(1, 2))
  testthat::expect_true(2 %in% ids_left)
  testthat::expect_type(ids_left, "double")

  # Direction toward atom 3 => cluster should be {3,4}
  ids_right <- fetch_eleno_downstream_of_bond(m, bond_id = 2, direction_atom_id = 3)
  testthat::expect_setequal(ids_right, c(3, 4))
  testthat::expect_true(3 %in% ids_right)
})

testthat::test_that("Errors when direction atom not on the bond", {
  # Chain 1-2-3-4; ask for bond (2) but direction atom 4 (not an endpoint of bond 2)
  edges <- rbind(c(1,2), c(2,3), c(3,4))
  m <- .make_mol_from_edges(edges, n_atoms = 4)

  testthat::expect_error(
    fetch_eleno_downstream_of_bond(m, bond_id = 2, direction_atom_id = 4),
    "is not connected by bond id"
  )
})

testthat::test_that("Errors when bond_id is invalid", {
  edges <- rbind(c(1,2), c(2,3))
  m <- .make_mol_from_edges(edges, n_atoms = 3)

  testthat::expect_error(
    fetch_eleno_downstream_of_bond(m, bond_id = 999, direction_atom_id = 2)
  )
})

testthat::test_that("Non-splitting break (cycle) raises 'All points remain connected.'", {
  # Triangle cycle 1-2-3-1; remove one endpoint of bond 1-2
  edges <- rbind(
    c(1, 2),
    c(2, 3),
    c(3, 1)
  )
  m <- .make_mol_from_edges(edges, n_atoms = 3)

  # Bond 1 (1-2), pick direction 1 -> removing atom 2 leaves {1,3} still connected => single component
  testthat::expect_error(
    fetch_eleno_downstream_of_bond(m, bond_id = 1, direction_atom_id = 1),
    "All points remain connected\\."
  )
})

testthat::test_that("Return includes the direction atom and is numeric", {
  edges <- rbind(c(1,2), c(2,3), c(3,4))
  m <- .make_mol_from_edges(edges, n_atoms = 4)

  res <- fetch_eleno_downstream_of_bond(m, bond_id = 2, direction_atom_id = 3)
  testthat::expect_true(3 %in% res)
  testthat::expect_true(is.numeric(res))
})



# Test eleno connected by class -------------------------------------------

test_that("fetch_eleno_connected_by_bond returns the two endpoints in order", {
  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("C", "O", "H"),
    x = c(0, 1, 2),
    y = c(0, 0, 0),
    z = c(0, 0, 0)
  )
  bonds <- data.frame(
    bond_id = c(10, 11),
    origin_atom_id = c(1, 2),
    target_atom_id = c(2, 3),
    bond_type = c("1", "2")
  )
  m <- Molecule3D(name = "test", atoms = atoms, bonds = bonds)

  expect_equal(fetch_eleno_connected_by_bond(m, 10), c(1, 2))
  expect_equal(fetch_eleno_connected_by_bond(m, 11), c(2, 3))
})

test_that("fetch_eleno_connected_by_bond accepts bond_id as character too", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 1),
    y = c(0, 0),
    z = c(0, 0)
  )
  bonds <- data.frame(
    bond_id = 10,
    origin_atom_id = 1,
    target_atom_id = 2,
    bond_type = "1"
  )
  m <- Molecule3D(name = "test", atoms = atoms, bonds = bonds)

  # `%in%` will coerce types, so this should still work
  expect_equal(fetch_eleno_connected_by_bond(m, "10"), c(1, 2))
})

test_that("fetch_eleno_connected_by_bond returns numeric(0) for missing bond_id", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 1),
    y = c(0, 0),
    z = c(0, 0)
  )
  bonds <- data.frame(
    bond_id = 10,
    origin_atom_id = 1,
    target_atom_id = 2,
    bond_type = "1"
  )
  m <- Molecule3D(name = "test", atoms = atoms, bonds = bonds)

  res <- fetch_eleno_connected_by_bond(m, 99)
  expect_identical(res, numeric(0))
})

test_that("fetch_eleno_connected_by_bond enforces Molecule3D class", {
  expect_error(
    fetch_eleno_connected_by_bond(list(), 10),
    regexp = "structures::Molecule3D"
  )
})



