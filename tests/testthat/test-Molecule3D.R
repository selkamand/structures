

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
  expect_equal(cmb@anchor, c(x=1,y=1,z=1))
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
  expect_equal(cmb@anchor, c(x=1,y=1,z=1))
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


# Test Anchor Defaults-------------------------------------------------------------

test_that("constructor defaults anchor to atom center when anchor is NULL", {
  atoms <- data.frame(
    eleno = c(1, 2, 10),
    elena = c("C", "O", "H"),
    x = c(0, 2, 4),
    y = c(-1, 1, 3),
    z = c(10, 12, 14)
  )
  bonds <- data.frame(
    bond_id = 1,
    origin_atom_id = 1,
    target_atom_id = 2,
    bond_type = "1"
  )

  m <- Molecule3D(name = "test", atoms = atoms, bonds = bonds, anchor = NULL)

  expected <- c(
    x = mean(atoms$x, na.rm = TRUE),
    y = mean(atoms$y, na.rm = TRUE),
    z = mean(atoms$z, na.rm = TRUE)
  )

  expect_equal(m@anchor, expected)
})

test_that("constructor anchor falls back to 0,0,0 when atoms are empty", {
  m <- Molecule3D(name = "empty", atoms = minimal_atoms(), bonds = minimal_bonds(), anchor = NULL)
  expect_equal(m@anchor, c(x=0, y=0, z=0))
})

test_that("explicit anchor is respected (no override)", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(100, 200),
    y = c(0, 0),
    z = c(-5, 5)
  )

  a <- c(x=1, y=2, z=3)
  m <- Molecule3D(name = "explicit", atoms = atoms, bonds = minimal_bonds(), anchor = a)
  expect_equal(m@anchor, a)
})

test_that("transform_molecule preserves original atom and bond IDs (non-contiguous) and moves anchor", {
  atoms <- data.frame(
    eleno = c(1, 2, 5, 8),           # non-contiguous IDs
    elena = c("C", "O", "H", "H"),
    x = c(0, 1, 2, 3),
    y = c(0, 0, 0, 0),
    z = c(0, 0, 0, 0)
  )
  bonds <- data.frame(
    bond_id = c(10, 20),
    origin_atom_id = c(1, 5),
    target_atom_id = c(2, 8),
    bond_type = c("1", "1")
  )

  m <- Molecule3D(name = "ids", atoms = atoms, bonds = bonds, anchor = NULL)

  # Save IDs and current anchor
  atom_ids_before <- m@atom_ids
  bond_ids_before <- m@bond_ids
  anchor_before <- m@anchor

  # Simple translation
  shift <- c(10, -2, 7)
  tfn <- function(p) c(x = p["x"] + shift[1], y = p["y"] + shift[2], z = p["z"] + shift[3])

  m2 <- transform_molecule(m, tfn)

  # IDs unchanged
  expect_equal(sort(m2@atom_ids), sort(atom_ids_before))
  expect_equal(sort(m2@bond_ids), sort(bond_ids_before))

  # browser()
  # Anchor moved consistently
  expect_equal(m2@anchor, anchor_before + shift)
})

test_that("anchor defaulting uses z coordinate (regression guard)", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 0),
    y = c(0, 0),
    z = c(5, 7)
  )
  m <- Molecule3D(name = "reg", atoms = atoms, bonds = minimal_bonds(), anchor = NULL)
  # If z is mistakenly computed from y, this would be 0 instead of 6
  expect_equal(m@anchor, c(x=0, y=0, z=6))
})

test_that("anchor default handles NA coordinates by ignoring them (if any present)", {
  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("C", "O", "H"),
    x = c(0, NA, 2),
    y = c(NA, 4, 6),
    z = c(1, 3, NA)
  )
  # When NA present and na.rm=TRUE, means are computed from available values
  m <- Molecule3D(name = "na", atoms = atoms, bonds = minimal_bonds(), anchor = NULL)
  expect_equal(
    m@anchor,
    c(x=mean(c(0, 2), na.rm = TRUE), y=mean(c(4, 6), na.rm = TRUE), z=mean(c(1, 3), na.rm = TRUE))
  )
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


# Test is_molecule --------------------------------------------------------

test_that("is_molecule returns TRUE for Molecule3D objects", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 1), y = c(0, 0), z = c(0, 0)
  )
  bonds <- data.frame(
    bond_id = 1,
    origin_atom_id = 1,
    target_atom_id = 2
  )

  mol <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds)
  expect_true(is_molecule(mol))
})

test_that("is_molecule returns FALSE for non-molecule inputs", {
  expect_false(is_molecule(list()))
  expect_false(is_molecule("not a molecule"))
  expect_false(is_molecule(1:3))

  # Also FALSE for other structures classes
  ax <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  expect_false(is_molecule(ax))
})



# Test Symmetry Axes -----------------------------------------------------------
test_that("Molecule3D initializes with empty symmetry axes", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C","O"),
    x = c(0, 1), y = c(0, 0), z = c(0, 0)
  )
  bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)

  m <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds, symmetry_axes = list())

  expect_false(m@contains_symmetry_axes)
  expect_null(m@symmetry_axes_orders)
  expect_equal(length(m@symmetry_axes), 0L)
})

test_that("add_symmetry_axis appends a valid SymAxis and updates derived properties", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C","O"),
    x = c(0, 1), y = c(0, 0), z = c(0, 0)
  )
  bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)

  m <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds)

  ax2 <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  ax3 <- SymAxis(Cn = 3L, posA = c(1,0,0), posB = c(1,0,1))

  m <- add_symmetry_axis(m, ax2)
  expect_true(m@contains_symmetry_axes)
  expect_equal(length(m@symmetry_axes), 1L)
  expect_setequal(m@symmetry_axes_orders, 2)

  m <- add_symmetry_axis(m, ax3)
  expect_equal(length(m@symmetry_axes), 2L)
  expect_setequal(m@symmetry_axes_orders, c(2, 3))
})

test_that("fetch_all_symmetry_axes_with_order returns only matching axes", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C","O"),
    x = c(0, 1), y = c(0, 0), z = c(0, 0)
  )
  bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)

  m <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds)

  ax2a <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  ax2b <- SymAxis(Cn = 2L, posA = c(0,1,0), posB = c(0,1,1))
  ax3  <- SymAxis(Cn = 3L, posA = c(1,0,0), posB = c(1,0,1))

  m <- add_symmetry_axis(m, ax2a)
  m <- add_symmetry_axis(m, ax2b)
  m <- add_symmetry_axis(m, ax3)

  only_C2 <- fetch_all_symmetry_axes_with_order(m, 2L)
  only_C3 <- fetch_all_symmetry_axes_with_order(m, 3L)
  only_C5 <- fetch_all_symmetry_axes_with_order(m, 5L)

  expect_equal(length(only_C2), 2L)
  expect_true(all(vapply(only_C2, function(ax) ax@Cn, integer(1)) == 2L))

  expect_equal(length(only_C3), 1L)
  expect_equal(only_C3[[1]]@Cn, 3L)

  expect_true(is.null(only_C5) || length(only_C5) == 0L)
})

test_that("constructor rejects non-SymAxis entries in symmetry_axes", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C","O"),
    x = c(0, 1), y = c(0, 0), z = c(0, 0)
  )
  bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)

  # symmetry_axes contains an invalid element (numeric 1)
  expect_error(
    Molecule3D(name = "CO", atoms = atoms, bonds = bonds, symmetry_axes = list(1)),
    "SymAxis|symmetry_axes", ignore.case = TRUE
  )
})

test_that("returns NULL when no symmetry axes are present", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C","O"),
    x = c(0, 1), y = c(0, 0), z = c(0, 0)
  )
  bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)

  m <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds, symmetry_axes = list())
  out <- fetch_all_symmetry_axes_with_order(m, 2L)
  expect_null(out)
})

test_that("returns only axes with matching Cn and preserves multiplicity", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C","O"),
    x = c(0, 1), y = c(0, 0), z = c(0, 0)
  )
  bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)
  m <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds)

  ax2a <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  ax2b <- SymAxis(Cn = 2L, posA = c(0,1,0), posB = c(0,1,1))
  ax3  <- SymAxis(Cn = 3L, posA = c(1,0,0), posB = c(1,0,1))

  m <- add_symmetry_axis(m, ax2a)
  m <- add_symmetry_axis(m, ax2b)
  m <- add_symmetry_axis(m, ax3)

  only_C2 <- fetch_all_symmetry_axes_with_order(m, 2L)
  only_C3 <- fetch_all_symmetry_axes_with_order(m, 3L)

  expect_equal(length(only_C2), 2L)
  expect_true(all(vapply(only_C2, function(ax) ax@Cn, integer(1)) == 2L))

  expect_equal(length(only_C3), 1L)
  expect_identical(only_C3[[1]]@Cn, 3L)
})

test_that("numeric vs integer Cn inputs behave equivalently", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C","O"),
    x = c(0, 1), y = c(0, 0), z = c(0, 0)
  )
  bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)
  m <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds)

  m <- add_symmetry_axis(m, SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1)))

  out_int <- fetch_all_symmetry_axes_with_order(m, 2L)
  out_num <- fetch_all_symmetry_axes_with_order(m, 2)    # double
  out_num0 <- fetch_all_symmetry_axes_with_order(m, 2.0) # double

  expect_equal(length(out_int), 1L)
  expect_equal(length(out_num), 1L)
  expect_equal(length(out_num0), 1L)
})

test_that("returns empty list (not NULL) when axes exist but none match Cn", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C","O"),
    x = c(0, 1), y = c(0, 0), z = c(0, 0)
  )
  bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)
  m <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds)

  m <- add_symmetry_axis(m, SymAxis(Cn = 3L, posA = c(0,0,0), posB = c(0,0,1)))

  out <- fetch_all_symmetry_axes_with_order(m, 2L)
  expect_true(is.list(out))
  expect_identical(length(out), 0L)
})


test_that("add_symmetry_axis rejects wrong classes", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C","O"),
    x = c(0, 1), y = c(0, 0), z = c(0, 0)
  )
  bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)

  m <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds)

  # wrong: symmetry_axis not a SymAxis
  expect_error(add_symmetry_axis(m, 123), "SymAxis", ignore.case = TRUE)

  # wrong: molecule not a Molecule3D
  ax <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  expect_error(add_symmetry_axis(list(), ax), "Molecule3D", ignore.case = TRUE)
})



# Test Symmetry Axes IDs ---------------------------------------------------------------------
test_that("unnamed symmetry_axes get sequential numeric IDs via constructor", {
  ax1 <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  ax2 <- SymAxis(Cn = 3L, posA = c(1,0,0), posB = c(1,0,1))

  m <- Molecule3D(
    name = "X",
    atoms = minimal_atoms(),
    bonds = minimal_bonds(),
    symmetry_axes = list(ax1, ax2)  # both unnamed
  )

  expect_true(m@contains_symmetry_axes)
  expect_identical(names(m@symmetry_axes), c("1","2"))
})

test_that("partially named symmetry_axes: blanks are filled with fresh numeric IDs", {
  ax1 <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  ax2 <- SymAxis(Cn = 3L, posA = c(1,0,0), posB = c(1,0,1))
  ax3 <- SymAxis(Cn = 3L, posA = c(2,0,0), posB = c(2,0,1))

  lst <- list(ax1, ax2, ax3)
  names(lst) <- c("A", "", "")  # assign blanks post-creation

  m <- Molecule3D(name = "X", atoms = minimal_atoms(), bonds = minimal_bonds(), symmetry_axes = lst)

  ids <- names(m@symmetry_axes)
  expect_length(ids, 3L)
  expect_true(all(nzchar(ids)))
  expect_true("A" %in% ids)
  expect_true(all(grepl("^[0-9]+$", setdiff(ids, "A"))))
  expect_length(unique(ids), 3L)
})

test_that("pre-existing numeric IDs cause fresh IDs to start at max(existing)+1", {
  ax1 <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  ax2 <- SymAxis(Cn = 3L, posA = c(1,0,0), posB = c(1,0,1))
  ax3 <- SymAxis(Cn = 4L, posA = c(2,0,0), posB = c(2,0,1))

  lst <- list(ax1, ax2, ax3)
  names(lst) <- c("2", "", "")  # keep "2", leave two blanks

  m <- Molecule3D(name = "X", atoms = minimal_atoms(), bonds = minimal_bonds(), symmetry_axes = lst)

  ids <- names(m@symmetry_axes)
  expect_setequal(ids, c("2","3","4"))
})

test_that("setting symmetry_axes with duplicate IDs errors (validator)", {
  ax1 <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  ax2 <- SymAxis(Cn = 3L, posA = c(1,0,0), posB = c(1,0,1))

  m <- Molecule3D(name = "X", atoms = minimal_atoms(), bonds = minimal_bonds())

  m@symmetry_axes <- list("1" = ax1)
  expect_error(
    { m@symmetry_axes <- list("1" = ax1, "1" = ax2) },
    regexp = "unique|duplicate|ID|identifier",
    ignore.case = TRUE
  )
})

test_that("add_symmetry_axis appends with next numeric ID when none exist", {
  m <- Molecule3D(name = "X", atoms = minimal_atoms(), bonds = minimal_bonds())
  ax <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))

  m <- add_symmetry_axis(m, ax)
  expect_identical(names(m@symmetry_axes), "1")

  m <- add_symmetry_axis(m, ax)
  expect_identical(names(m@symmetry_axes), c("1","2"))
})

test_that("add_symmetry_axis uses max numeric among existing IDs, ignoring non-numeric names", {
  ax <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))

  m <- Molecule3D(
    name = "X",
    atoms = minimal_atoms(),
    bonds = minimal_bonds(),
    symmetry_axes = list("foo" = ax, "7" = ax)
  )

  m <- add_symmetry_axis(m, ax)
  ids <- names(m@symmetry_axes)

  expect_true("8" %in% ids)          # next after 7
  expect_true(all(nzchar(ids)))
  expect_length(unique(ids), length(ids))
})

test_that("add_symmetry_axis starts at '1' if all existing IDs are non-numeric", {
  ax <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))

  m <- Molecule3D(
    name = "X",
    atoms = minimal_atoms(),
    bonds = minimal_bonds(),
    symmetry_axes = list("alpha" = ax, "beta" = ax)
  )

  m <- add_symmetry_axis(m, ax)
  expect_true("1" %in% names(m@symmetry_axes))
})

test_that("replacing symmetry_axes preserves uniqueness and non-emptiness of IDs", {
  ax <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  m <- Molecule3D(name = "X", atoms = minimal_atoms(), bonds = minimal_bonds())

  m@symmetry_axes <- list(ax, ax, ax)
  expect_setequal(names(m@symmetry_axes), c("1","2","3"))

  m@symmetry_axes <- list(ax, ax)
  names(m@symmetry_axes) <- c("X", "")  # assign one blank post-creation
  ids <- names(m@symmetry_axes)
  expect_true(all(nzchar(ids)))
  expect_length(unique(ids), 2L)
})

test_that("transform_molecule keeps the symmetry axis IDs intact", {
  ax1 <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  ax2 <- SymAxis(Cn = 3L, posA = c(1,0,0), posB = c(1,0,1))

  m <- Molecule3D(
    name = "X",
    atoms = data.frame(eleno = 1, elena = "C", x = 0, y = 0, z = 0),
    bonds = minimal_bonds(),
    symmetry_axes = list("a" = ax1, "2" = ax2)
  )

  ids_before <- names(m@symmetry_axes)
  trans <- function(p) c(x = p["x"] + 1, y = p["y"] + 2, z = p["z"] + 3)
  m2 <- transform_molecule(m, trans)

  expect_identical(names(m2@symmetry_axes), ids_before)
})



# Test Symmetry Axes Dataframe --------------------------------------------

test_that("symmetry_axes_dataframe returns zero-row df with expected columns when no axes", {
  atoms <- data.frame(eleno = c(1,2), elena = c("C","O"), x = c(0,1), y = c(0,0), z = c(0,0))
  bonds <- data.frame(bond_id = numeric(0), origin_atom_id = numeric(0), target_atom_id = numeric(0))
  m <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds, symmetry_axes = list())

  df <- m@symmetry_axes_dataframe

  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 0L)
  expect_identical(
    names(df),
    c("id", "label", "Cn", "x", "y", "z", "xend", "yend", "zend")
  )
})

test_that("symmetry_axes_dataframe returns one row with correct values and id", {
  atoms <- data.frame(eleno = c(1,2), elena = c("C","O"), x = c(0,1), y = c(0,0), z = c(0,0))
  bonds <- data.frame(bond_id = numeric(0), origin_atom_id = numeric(0), target_atom_id = numeric(0))

  ax <- SymAxis(Cn = 3L, posA = c(0,0,0), posB = c(0,0,1), label = "axis-A")
  m  <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds, symmetry_axes = list("5" = ax))

  df <- m@symmetry_axes_dataframe

  expect_equal(nrow(df), 1L)
  expect_identical(df$id, "5")
  expect_identical(df$label, "axis-A")
  expect_identical(as.integer(df$Cn), 3L)
  expect_equal(unname(c(df$x, df$y, df$z)), c(0,0,0))
  expect_equal(unname(c(df$xend, df$yend, df$zend)), c(0,0,1))
})

test_that("symmetry_axes_dataframe preserves order and IDs for multiple axes", {
  atoms <- data.frame(eleno = c(1,2), elena = c("C","O"), x = c(0,1), y = c(0,0), z = c(0,0))
  bonds <- data.frame(bond_id = numeric(0), origin_atom_id = numeric(0), target_atom_id = numeric(0))

  ax1 <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1), label = "A")
  ax2 <- SymAxis(Cn = 3L, posA = c(1,0,0), posB = c(1,0,1), label = "B")
  ax3 <- SymAxis(Cn = 4L, posA = c(2,0,0), posB = c(2,0,1), label = "C")

  m <- Molecule3D(
    name = "CO",
    atoms = atoms,
    bonds = bonds,
    symmetry_axes = list("10" = ax1, "20" = ax2, "30" = ax3)
  )

  df <- m@symmetry_axes_dataframe

  expect_equal(nrow(df), 3L)
  expect_identical(df$id, c("10","20","30"))
  expect_identical(df$label, c("A","B","C"))
  expect_identical(as.integer(df$Cn), c(2L, 3L, 4L))
  expect_equal(df$x, c(0,1,2))
  expect_equal(df$xend, c(0,1,2))
})

test_that("symmetry_axes_dataframe updates after transform_molecule and preserves IDs", {
  atoms <- data.frame(eleno = c(1,2), elena = c("C","O"), x = c(0,1), y = c(0,0), z = c(0,0))
  bonds <- data.frame(bond_id = numeric(0), origin_atom_id = numeric(0), target_atom_id = numeric(0))
  ax <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(1,0,0), label = "A")
  m  <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds, symmetry_axes = list("7" = ax))

  # Translate by (+1, +2, +3)
  translate <- function(p, dx=0, dy=0, dz=0) {
    c(x = p["x"] + dx, y = p["y"] + dy, z = p["z"] + dz)
  }
  m2 <- transform_molecule(m, translate, dx = 1, dy = 2, dz = 3)

  # Dataframe should reflect transformed axis endpoints; id unchanged
  df <- m2@symmetry_axes_dataframe
  expect_equal(df$id, "7")
  expect_equal(unname(c(df$x, df$y, df$z)), c(1,2,3))
  expect_equal(unname(c(df$xend, df$yend, df$zend)), c(2,2,3))
})

test_that("symmetry_axes_dataframe cooperates with add_symmetry_axis auto-ID generation", {
  atoms <- data.frame(eleno = c(1,2), elena = c("C","O"), x = c(0,1), y = c(0,0), z = c(0,0))
  bonds <- data.frame(bond_id = numeric(0), origin_atom_id = numeric(0), target_atom_id = numeric(0))
  m  <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds)

  ax1 <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1), label = "A")
  ax2 <- SymAxis(Cn = 3L, posA = c(0,1,0), posB = c(0,1,1), label = "B")

  m <- add_symmetry_axis(m, ax1)
  m <- add_symmetry_axis(m, ax2)

  df <- m@symmetry_axes_dataframe

  # IDs should be distinct, non-empty, and match names(m@symmetry_axes)
  expect_false(any(duplicated(df$id)))
  expect_true(all(nzchar(df$id)))
  expect_identical(df$id, names(m@symmetry_axes))
})


# Test Print Generic ------------------------------------------------------


test_that("print.Molecule3D runs without error and returns invisible Molecule3D", {

  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 1),
    y = c(0, 0),
    z = c(0, 0)
  )
  bonds <- data.frame(
    bond_id = 1,
    origin_atom_id = 1,
    target_atom_id = 2
  )
  m <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds)

  out <- NULL
  invisible(capture.output({
    out <- expect_no_error(print(m))
  }))

  expect_identical(out, m)
})

test_that("print.Molecule3D works even with symmetry axes present", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C","O"),
    x = c(0, 1),
    y = c(0, 0),
    z = c(0, 0)
  )
  bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)
  ax <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  m <- Molecule3D(name = "CO", atoms = atoms, bonds = bonds, symmetry_axes = list(ax))

  out <- NULL
  invisible(capture.output({
    out <- expect_no_error(print(m))
  }))

  expect_identical(out, m)
})


# Test Transformations ----------------------------------------------------

# tests/testthat/test-transform_molecule.R

test_that("transform_molecule() translates atoms, anchor, and symmetry axes", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 1),
    y = c(0, 0),
    z = c(0, 0)
  )
  bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)
  axis  <- SymAxis(Cn = 2L, posA = c(0, 0, 0), posB = c(0, 0, 1))

  m <- Molecule3D(
    name = "CO",
    atoms = atoms,
    bonds = bonds,
    symmetry_axes = list(axis),
    anchor = c(0, 0, 0)
  )

  translate_by_one <- function(p) p + 1

  m2 <- transform_molecule(m, translate_by_one)

  # atoms
  expected_positions <- m@atoms[, c("x", "y", "z")] + 1
  expect_equal(m2@atoms[, c("x", "y", "z")], expected_positions)

  # anchor
  expect_equal(m2@anchor, m@anchor + 1)

  # symmetry axes
  expect_length(m2@symmetry_axes, 1)
  ax_new <- m2@symmetry_axes[[1]]
  expect_equal(ax_new@posA, axis@posA + 1)
  expect_equal(ax_new@posB, axis@posB + 1)
  expect_identical(ax_new@Cn, axis@Cn)
})

test_that("transform_molecule() works when no symmetry axes exist", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 1),
    y = c(0, 0),
    z = c(0, 0)
  )
  bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)
  m <- Molecule3D("CO", atoms = atoms, bonds = bonds)

  translate_by_one <- function(p) p + 1
  m2 <- transform_molecule(m, translate_by_one)

  expect_equal(m2@atoms[, c("x", "y", "z")], m@atoms[, c("x", "y", "z")] + 1)
  expect_false(m2@contains_symmetry_axes)
})

test_that("transform_molecule() preserves Cn and rotates correctly about Z", {
  atoms <- data.frame(
    eleno = 1,
    elena = "C",
    x = 1, y = 0, z = 0
  )
  bonds <- minimal_bonds()
  axis  <- SymAxis(Cn = 3L, posA = c(0, 0, 0), posB = c(0, 0, 1))

  m <- Molecule3D("Rot", atoms = atoms, bonds = bonds,
                  symmetry_axes = list(axis),
                  anchor = c(0, 0, 0))

  rotate90_z <- function(p) {
    a <- pi/2
    c(
      x = p[1] * cos(a) - p[2] * sin(a),
      y = p[1] * sin(a) + p[2] * cos(a),
      z = p[3]
    )
  }

  m2 <- transform_molecule(m, rotate90_z)

  # atom moved (1,0,0) -> (0,1,0)
  expect_equal(unname(unlist(m2@atoms[1, c("x", "y", "z")])),
               c(0, 1, 0), tolerance = 1e-8)

  # axis along Z unchanged by rotation about Z
  expect_equal(m2@symmetry_axes[[1]]@posA, axis@posA, tolerance = 1e-12)
  expect_equal(m2@symmetry_axes[[1]]@posB, axis@posB, tolerance = 1e-12)
  expect_equal(m2@symmetry_axes[[1]]@Cn, axis@Cn)
})

test_that("transform_molecule() handles multiple symmetry axes and preserves orders", {
  atoms <- data.frame(
    eleno = c(1, 2, 3),
    elena = c("C", "O", "H"),
    x = c(0, 1, -1),
    y = c(0, 0,  0),
    z = c(0, 0,  0)
  )
  bonds <- data.frame(
    bond_id = c(1, 2),
    origin_atom_id = c(1, 1),
    target_atom_id = c(2, 3)
  )
  ax1 <- SymAxis(Cn = 2L, posA = c(0, 0, 0), posB = c(0, 0, 1))
  ax2 <- SymAxis(Cn = 3L, posA = c(1, 0, 0), posB = c(1, 1, 0))

  m <- Molecule3D(
    name = "COH",
    atoms = atoms,
    bonds = bonds,
    symmetry_axes = list(ax1, ax2),
    anchor = c(0, 0, 0)
  )
  expect_setequal(m@symmetry_axes_orders, c(2, 3))

  # translation by arbitrary vector
  v <- c(2, -1, 5)
  translate_vec <- function(p) p + v

  m2 <- transform_molecule(m, translate_vec)

  # atoms translated
  expect_equal(m2@atoms[, c("x", "y", "z")], m@atoms[, c("x", "y", "z")] +
                 matrix(rep(v, each = nrow(m@atoms)), ncol = 3, byrow = FALSE))

  # axes translated
  expect_equal(m2@symmetry_axes[[1]]@posA, ax1@posA + v)
  expect_equal(m2@symmetry_axes[[1]]@posB, ax1@posB + v)
  expect_equal(m2@symmetry_axes[[2]]@posA, ax2@posA + v)
  expect_equal(m2@symmetry_axes[[2]]@posB, ax2@posB + v)

  # orders unchanged
  expect_setequal(m2@symmetry_axes_orders, c(2, 3))

  # fetch_all_symmetry_axes_with_order remains coherent after transform
  c2_axes <- fetch_all_symmetry_axes_with_order(m2, 2)
  c3_axes <- fetch_all_symmetry_axes_with_order(m2, 3)
  expect_equal(length(c2_axes), 1)
  expect_equal(length(c3_axes), 1)
  expect_true(is_symmetry_axis(c2_axes[[1]]))
  expect_true(is_symmetry_axis(c3_axes[[1]]))
})

test_that("transform_molecule() forwards ... arguments to transformation function", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(1, 2),
    y = c(1, 2),
    z = c(1, 2)
  )
  bonds <- minimal_bonds()
  ax <- SymAxis(Cn = 2L, posA = c(1, 1, 1), posB = c(2, 2, 2))

  m <- Molecule3D("scale", atoms = atoms, bonds = bonds,
                  symmetry_axes = list(ax), anchor = c(1, 1, 1))

  # scale by factor passed through ...
  scale_by <- function(p, k) p * k

  m2 <- transform_molecule(m, scale_by, k = 3)

  # atoms scaled
  expect_equal(m2@atoms[, c("x", "y", "z")], m@atoms[, c("x", "y", "z")] * 3)

  # anchor scaled
  expect_equal(m2@anchor, m@anchor * 3)

  # axes scaled
  expect_equal(m2@symmetry_axes[[1]]@posA, ax@posA * 3)
  expect_equal(m2@symmetry_axes[[1]]@posB, ax@posB * 3)
})

test_that("transform_molecule() keeps bond_positions coherent with atoms after transform", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 1),
    y = c(0, 0),
    z = c(0, 0)
  )
  bonds <- data.frame(bond_id = 1, origin_atom_id = 1, target_atom_id = 2)
  m <- Molecule3D("CO", atoms = atoms, bonds = bonds)

  v <- c(0.5, -2, 7)
  translate_vec <- function(p) p + v
  m2 <- transform_molecule(m, translate_vec)

  # recompute (property is read-only getter, so reading it is enough)
  bp <- m2@bond_positions

  # Expect bond start equals translated atom 1, end equals translated atom 2
  expect_equal(unname(unlist(bp[1, c("x", "y", "z")])),
               unname(unlist(m2@atoms[m2@atoms$eleno == 1, c("x", "y", "z")])))
  expect_equal(unname(unlist(bp[1, c("xend", "yend", "zend")])),
               unname(unlist(m2@atoms[m2@atoms$eleno == 2, c("x", "y", "z")])))
})

test_that("transform_molecule() returns a Molecule3D and mutates as expected", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 1),
    y = c(0, 0),
    z = c(0, 0)
  )
  bonds <- minimal_bonds()
  m <- Molecule3D("CO", atoms = atoms, bonds = bonds, anchor = c(0, 0, 0))
  saved_m_atoms <- m@atoms

  trans <- function(p) p + 2
  m2 <- transform_molecule(m, trans)

  expect_true(inherits(m2, "structures::Molecule3D"))

  # Ensure transform_molecule returns new Molecule3D object without modifying the original object
  expect_identical(
    saved_m_atoms,
    m@atoms
  )
})

test_that("transform_molecule() errors on wrong inputs", {
  atoms <- data.frame(
    eleno = 1, elena = "C", x = 0, y = 0, z = 0
  )
  m <- Molecule3D("bad", atoms = atoms, bonds = minimal_bonds())

  not_a_function <- 123
  expect_error(transform_molecule(m, not_a_function))

  # wrong class for x
  expect_error(transform_molecule(list(), function(p) p))
})

test_that("transform_molecule() surfaces errors from malformed transformation return values", {
  atoms <- data.frame(
    eleno = c(1, 2),
    elena = c("C", "O"),
    x = c(0, 1),
    y = c(0, 0),
    z = c(0, 0)
  )
  m <- Molecule3D("CO", atoms = atoms, bonds = minimal_bonds(), anchor = c(0, 0, 0))
  # Return wrong length -> vapply should error
  bad_transform <- function(p) c(99, 100) # length 2 instead of 3
  expect_error(transform_molecule(m, bad_transform))
})


test_that("transform_molecule preserves original symmetry axis IDs even when non-contiguous or non-numeric", {
  # Build three axes
  axA <- SymAxis(Cn = 2L, posA = c(0,0,0), posB = c(0,0,1))
  axB <- SymAxis(Cn = 3L, posA = c(1,0,0), posB = c(1,0,1))
  axC <- SymAxis(Cn = 4L, posA = c(2,0,0), posB = c(2,0,1))

  # Name them with sparse, mixed IDs
  axes <- list(axA, axB, axC)
  names(axes) <- c("axis-99", "5", "foo")

  # Minimal molecule (one atom, no bonds) with those axes
  atoms <- data.frame(eleno = 1, elena = "C", x = 0, y = 0, z = 0)
  m <- Molecule3D(name = "X", atoms = atoms, bonds = minimal_bonds(), symmetry_axes = axes)

  ids_before <- names(m@symmetry_axes)

  # Simple translation so we know a transform occurred
  trans <- function(p) c(x = p["x"] + 1, y = p["y"] + 2, z = p["z"] + 3)

  m2 <- transform_molecule(m, trans)

  # IDs should be identical (order and values)
  expect_identical(names(m2@symmetry_axes), ids_before)

  # Sanity check: positions actually changed for one axis, but ID stayed the same
  before_posA <- m@symmetry_axes[["axis-99"]]@posA
  after_posA  <- m2@symmetry_axes[["axis-99"]]@posA
  expect_equal(after_posA, before_posA + c(1,2,3))
})


