test_that("list_to_matrix matches %in% membership for character sets", {
  sets <- list(set1 = letters[1:5], set2 = letters[2:6], set3 = letters[1:7])
  m <- list_to_matrix(sets)
  expect_type(m, "integer")
  expect_equal(colSums(m), lengths(sets))
  oracle <- sapply(sets, function(s) as.integer(letters[1:7] %in% s))
  rownames(oracle) <- letters[1:7]
  expect_identical(m, oracle)
})

test_that("list_to_matrix handles duplicates, empty sets, order and NA", {
  sets <- list(A = c("z", "y", "z"), B = character(0), C = c("x", "y"))
  m <- list_to_matrix(sets)
  expect_identical(rownames(m), c("z", "y", "x"))
  expect_equal(unname(m[, "A"]), c(1L, 1L, 0L))
  expect_equal(unname(m[, "B"]), c(0L, 0L, 0L))
  expect_equal(unname(m[, "C"]), c(0L, 1L, 1L))
  expect_message(m2 <- list_to_matrix(list(A = c("a", NA), B = "b")), "NA values")
  expect_identical(rownames(m2), c("a", "b"))
  expect_equal(sum(m2), 2)
  expect_equal(dim(list_to_matrix(list(A = "a"))), c(1L, 1L))
})

test_that("list_to_matrix indexes numeric sets by value, not position", {
  m <- list_to_matrix(list(a = c(2, 1), b = 1))
  expect_equal(m["1", "b"], 1L)
  expect_equal(m["2", "b"], 0L)
  m2 <- list_to_matrix(list(a = c(10, 20), b = 20))
  expect_equal(unname(m2[, "a"]), c(1L, 1L))
  expect_equal(unname(m2[, "b"]), c(0L, 1L))
})

test_that("list_to_matrix indexes factor sets by label", {
  m <- list_to_matrix(list(a = factor(c("x", "y")), b = factor("y")))
  expect_equal(m["y", "b"], 1L)
  expect_equal(m["x", "b"], 0L)
})
