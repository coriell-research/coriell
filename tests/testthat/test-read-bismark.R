# Write bismark coverage strings to temp files and return the named paths
write_cov <- function(x) {
  paths <- vapply(seq_along(x), function(i) {
    f <- tempfile(fileext = ".cov")
    writeLines(x[[i]], f, sep = "")
    f
  }, character(1))
  names(paths) <- names(x)
  paths
}

# chrom, start, end, percent, methylated, unmethylated
s1 <- "chr1\t100\t100\t50\t5\t5\nchr1\t200\t200\t100\t2\t0\nchr2\t50\t50\t0\t0\t10\n"
s2 <- "chr1\t100\t100\t80\t8\t2\nchr1\t200\t200\t100\t30\t0\nchr1\t300\t300\t25\t1\t3\n"
s3 <- "chr1\t100\t100\t20\t2\t8\nchr1\t200\t200\t100\t20\t0\nchr2\t50\t50\t40\t4\t6\n"

test_that("Default keeps only sites present in all samples", {
  files <- write_cov(c(S1 = s1, S2 = s2, S3 = s3))
  on.exit(unlink(files))
  dt <- suppressMessages(read_bismark(files))
  expect_s3_class(dt, "data.table")
  expect_identical(
    names(dt),
    c("Sample", "Chrom", "Start", "End", "Percent", "Me", "Un", "Coverage")
  )
  expect_setequal(unique(paste(dt$Chrom, dt$Start)), c("chr1 100", "chr1 200"))
  expect_equal(nrow(dt), 6L)
  expect_equal(dt$Coverage, dt$Me + dt$Un)
  expect_setequal(dt$Sample, c("S1", "S2", "S3"))
  r <- dt[dt$Sample == "S2" & dt$Start == 100, ]
  expect_equal(c(r$Percent, r$Me, r$Un, r$Coverage), c(80, 8, 2, 10))
})

test_that("Coverage threshold is inclusive", {
  files <- write_cov(c(S1 = s1, S2 = s2, S3 = s3))
  on.exit(unlink(files))
  # chr1:100 has coverage 10 everywhere; chr1:200 has 2 in S1
  expect_setequal(unique(suppressMessages(read_bismark(files, coverage = 10))$Start), 100L)
  expect_equal(nrow(suppressMessages(read_bismark(files, coverage = 11))), 0L)
  expect_setequal(unique(suppressMessages(read_bismark(files, coverage = 2))$Start), c(100L, 200L))
})

test_that("prop_samples threshold matches a base-R count of covered samples", {
  files <- write_cov(c(S1 = s1, S2 = s2, S3 = s3))
  on.exit(unlink(files))
  count_sites <- function(p, cov = 1) {
    dt <- suppressMessages(read_bismark(files, coverage = cov, prop_samples = p))
    length(unique(paste(dt$Chrom, dt$Start)))
  }
  expect_equal(count_sites(1), 2L)
  expect_equal(count_sites(2 / 3), 3L)
  expect_equal(count_sites(0.5), 3L)
  expect_equal(count_sites(1 / 3), 4L)
  all <- rbind(read.table(text = s1), read.table(text = s2), read.table(text = s3))
  all$cov <- all$V5 + all$V6
  for (cv in c(1, 5, 10, 25)) {
    for (p in c(1 / 3, 0.5, 2 / 3, 1)) {
      tab <- table(paste(all$V1, all$V2)[all$cov >= cv])
      expect_equal(count_sites(p, cv), sum(tab / 3 >= p - 1e-12), info = paste(cv, p))
    }
  }
})

test_that("prop_samples is robust to floating point error", {
  # 25 * 0.28 is 7 + 8.9e-16 in floating point; chr1:100 is covered in 7/25
  one <- "chr1\t100\t100\t50\t5\t5\nchr1\t200\t200\t50\t5\t5\n"
  two <- "chr1\t200\t200\t50\t5\t5\n"
  files <- write_cov(setNames(c(rep(one, 7), rep(two, 18)), paste0("S", 1:25)))
  on.exit(unlink(files))
  dt <- suppressMessages(read_bismark(files, prop_samples = 0.28))
  expect_setequal(unique(dt$Start), c(100L, 200L))
  expect_equal(sum(dt$Start == 100), 7L)
  # one more sample required -> chr1:100 dropped
  dt <- suppressMessages(read_bismark(files, prop_samples = 0.29))
  expect_setequal(unique(dt$Start), 200L)
})

test_that("remove_zero_var drops constant and single-sample sites", {
  files <- write_cov(c(S1 = s1, S2 = s2, S3 = s3))
  on.exit(unlink(files))
  # chr1:200 has Percent 100 in every sample
  dt <- suppressMessages(read_bismark(files, remove_zero_var = TRUE))
  expect_setequal(unique(dt$Start), 100L)
  expect_equal(nrow(dt), 3L)
  # chr1:300 is only in S2 (variance NA) and is dropped as well
  dt <- suppressMessages(read_bismark(files, prop_samples = 1 / 3, remove_zero_var = TRUE))
  expect_setequal(unique(paste(dt$Chrom, dt$Start)), c("chr1 100", "chr2 50"))
})

test_that("return_mats gives matrices consistent with the long output", {
  files <- write_cov(c(S1 = s1, S2 = s2, S3 = s3))
  on.exit(unlink(files))
  l <- suppressMessages(read_bismark(files, prop_samples = 1 / 3, return_mats = TRUE))
  expect_named(l, c("Coverage", "Percent", "Methylated", "Unmethylated"))
  for (m in l) {
    expect_true(is.matrix(m))
    expect_identical(colnames(m), c("S1", "S2", "S3"))
    expect_identical(rownames(m), c("chr1.100", "chr1.200", "chr1.300", "chr2.50"))
  }
  expect_equal(unname(l$Coverage["chr1.100", ]), c(10, 10, 10))
  expect_equal(unname(l$Percent["chr1.100", ]), c(50, 80, 20))
  expect_equal(unname(l$Methylated["chr1.200", ]), c(2, 30, 20))
  # missing sites: Coverage filled with 0, others NA
  expect_equal(unname(l$Coverage["chr1.300", ]), c(0, 4, 0))
  expect_true(all(is.na(l$Percent["chr1.300", c("S1", "S3")])))
  expect_true(is.na(l$Methylated["chr2.50", "S2"]))
  obs <- !is.na(l$Methylated)
  expect_equal(l$Coverage[obs], (l$Methylated + l$Unmethylated)[obs])
  dt <- suppressMessages(read_bismark(files, prop_samples = 1 / 3))
  expect_equal(sum(l$Coverage), sum(dt$Coverage))
  expect_equal(sum(obs), nrow(dt))
})

test_that("Matrix columns follow the input file order", {
  files <- write_cov(c(S2 = s1, S10 = s2, S1 = s3))
  on.exit(unlink(files))
  l <- suppressMessages(read_bismark(files, return_mats = TRUE))
  for (m in l) {
    expect_identical(colnames(m), c("S2", "S10", "S1"))
  }
  # values stay with their sample: S2 was read from s1
  expect_equal(unname(l$Percent["chr1.100", ]), c(50, 80, 20))
})

test_that("Unnamed files warn and are labelled by index", {
  files <- write_cov(list(s1, s2, s3))
  on.exit(unlink(files))
  expect_warning(dt <- suppressMessages(read_bismark(unname(files))), "not named")
  expect_setequal(dt$Sample, 1:3)
  expect_warning(
    l <- suppressMessages(read_bismark(unname(files), return_mats = TRUE)),
    "not named"
  )
  expect_identical(colnames(l$Percent), c("1", "2", "3"))
})

test_that("Messages report the filter settings", {
  files <- write_cov(c(S1 = s1, S2 = s2, S3 = s3))
  on.exit(unlink(files))
  expect_message(read_bismark(files), "Filtering for Coverage >= 1 in 100% of samples")
  expect_error(read_bismark(files, prop_samples = 0), "prop_samples")
})
