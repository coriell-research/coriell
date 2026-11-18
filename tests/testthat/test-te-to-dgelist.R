skip_if_not_installed("edgeR")

# Write strings to temp files and return the (named) paths
write_files <- function(x, fileext = ".cntTable") {
  paths <- vapply(seq_along(x), function(i) {
    f <- tempfile(fileext = fileext)
    writeLines(x[[i]], f, sep = "")
    f
  }, character(1))
  names(paths) <- names(x)
  paths
}

# Build a count table from ids and counts as TEtranscripts writes it
cnt_table <- function(ids, counts, bam = "sample.bam") {
  paste0("gene/TE\t", bam, "\n", paste0(ids, "\t", counts, "\n", collapse = ""))
}

tet_ids <- c(
  "\"ENSG01.1\"", "\"ENSG02.1\"", "\"ENSG03.1\"",
  "AluY:Alu:SINE", "Eulor10:DNA:DNA?", "(CATTC)n:Satellite:Satellite"
)
tel_ids <- c(
  "\"ENSG01.1\"", "\"ENSG02.1\"", "\"ENSG03.1\"",
  "AluY_dup1:AluY:Alu:SINE", "AluY_dup2:AluY:Alu:SINE",
  "(CATTC)n:(CATTC)n:Satellite:Satellite"
)

tet_files <- function() {
  write_files(c(
    S1 = cnt_table(tet_ids, c(10, 0, 3, 5, 0, 1), "S1.bam"),
    S2 = cnt_table(tet_ids, c(12, 0, 0, 9, 0, 2), "S2.bam")
  ))
}
tel_files <- function() {
  write_files(c(
    S1 = cnt_table(tel_ids, c(10, 0, 3, 5, 0, 1)),
    S2 = cnt_table(tel_ids, c(12, 0, 0, 9, 0, 2))
  ))
}

gene_gtf_lines <- c(
  "##description: test",
  "##provider: GENCODE",
  # Shuffled relative to the count tables; ENSG03.1 is absent
  "chr2\tHAVANA\tgene\t900\t1500\t.\t-\t.\tgene_id \"ENSG02.1\"; gene_type \"protein_coding\"; gene_name \"GENE2\"; level 2;",
  "chr2\tHAVANA\ttranscript\t900\t1400\t.\t-\t.\tgene_id \"ENSG02.1\"; transcript_id \"ENST02.1\"; gene_type \"protein_coding\"; gene_name \"GENE2\";",
  "chr1\tHAVANA\tgene\t100\t500\t.\t+\t.\tgene_id \"ENSG01.1\"; gene_type \"lncRNA\"; gene_name \"GENE1\"; level 2;"
)
te_gtf_lines <- c(
  "chr1\tUCSC_rmsk\texon\t1001\t1300\t10\t+\t.\tgene_id \"AluY\"; transcript_id \"AluY\"; family_id \"Alu\"; class_id \"SINE\"; gene_name \"AluY:TE\";",
  "chr3\tUCSC_rmsk\texon\t51\t100\t10\t+\t.\tgene_id \"(CATTC)n\"; transcript_id \"(CATTC)n\"; family_id \"Satellite\"; class_id \"Satellite\"; gene_name \"(CATTC)n:TE\";",
  "chr2\tUCSC_rmsk\texon\t2001\t2100\t10\t-\t.\tgene_id \"AluY\"; transcript_id \"AluY\"; family_id \"Alu\"; class_id \"SINE\"; gene_name \"AluY:TE\";",
  # Same subfamily name under a different family must not be merged
  "chr2\tUCSC_rmsk\texon\t5001\t5010\t10\t-\t.\tgene_id \"AluY\"; transcript_id \"AluY\"; family_id \"Other\"; class_id \"SINE\"; gene_name \"AluY:TE\";"
)
loc_lines <- c(
  "TE\tchromsome:start-stop:strand",
  "(CATTC)n\tchr3:51-100:+",
  "AluY_dup2\tchr2:2001-2100:-",
  "AluY_dup1\tchr1:1001-1300:+",
  "AluY_dup9\tchr9:1-10:+"
)

write_lines <- function(lines, fileext) {
  f <- tempfile(fileext = fileext)
  writeLines(lines, f)
  f
}
write_gz <- function(lines) {
  f <- tempfile(fileext = ".gtf.gz")
  con <- gzfile(f, "w")
  writeLines(lines, con)
  close(con)
  f
}

test_that("TEtranscripts tables are imported into a DGEList", {
  files <- tet_files()
  on.exit(unlink(files))
  y <- suppressMessages(TEtranscripts2DGEList(files))

  expect_s4_class(y, "DGEList")
  expect_identical(dim(y), c(6L, 2L))
  expect_identical(colnames(y), c("S1", "S2"))
  expect_true(is.integer(y$counts))
  expect_equal(unname(y$counts[, "S2"]), c(12, 0, 0, 9, 0, 2))
  expect_false(any(grepl("\"", rownames(y))))
  expect_identical(rownames(y$genes), rownames(y$counts))
  expect_identical(y$samples$bam, c("S1.bam", "S2.bam"))
  expect_equal(y$samples$lib.size, c(19, 23))

  g <- y$genes
  expect_identical(
    names(g),
    c("feature_id", "feature_type", "te_subfamily", "te_family", "te_class")
  )
  expect_identical(g$feature_type, rep(c("gene", "TE"), each = 3))
  expect_identical(g$te_subfamily[4:6], c("AluY", "Eulor10", "(CATTC)n"))
  expect_identical(g$te_family[4:6], c("Alu", "DNA", "Satellite"))
  expect_identical(g$te_class[4:6], c("SINE", "DNA?", "Satellite"))
  expect_true(all(is.na(g$te_subfamily[1:3])))
})

test_that("TElocal tables are imported and loci parsed", {
  files <- tel_files()
  on.exit(unlink(files))
  y <- suppressMessages(TEtranscripts2DGEList(files))

  expect_true("te_locus" %in% names(y$genes))
  expect_identical(
    y$genes$te_locus[4:6],
    c("AluY_dup1", "AluY_dup2", "(CATTC)n")
  )
  expect_identical(y$genes$te_subfamily[4:6], c("AluY", "AluY", "(CATTC)n"))
  expect_identical(y$genes$te_class[4:6], c("SINE", "SINE", "Satellite"))
  expect_message(TEtranscripts2DGEList(files), "Detected TElocal")
})

test_that("type is detected and validated", {
  tet <- tet_files()
  tel <- tel_files()
  mixed <- write_files(list(cnt_table(c(tet_ids, "AluY_dup1:AluY:Alu:SINE"), 1:7)))
  on.exit(unlink(c(tet, tel, mixed)))

  expect_message(TEtranscripts2DGEList(tet), "Detected TEtranscripts")
  expect_no_error(suppressMessages(TEtranscripts2DGEList(tet, type = "TEtranscripts")))
  expect_error(suppressMessages(TEtranscripts2DGEList(tet, type = "TElocal")), "look like TEtranscripts")
  expect_error(suppressMessages(TEtranscripts2DGEList(tel, type = "TEtranscripts")), "look like TElocal")
  expect_error(suppressMessages(TEtranscripts2DGEList(mixed)), "Unrecognized TE feature IDs")
  expect_error(TEtranscripts2DGEList(tet, type = "salmon"))
})

test_that("unnamed files are named by file name", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- file.path(dir, c("A.cntTable", "B.cntTable"))
  writeLines(cnt_table(tet_ids, 1:6), files[1], sep = "")
  writeLines(cnt_table(tet_ids, 6:1), files[2], sep = "")

  y <- suppressMessages(TEtranscripts2DGEList(files))
  expect_identical(colnames(y), c("A", "B"))

  expect_error(
    suppressMessages(TEtranscripts2DGEList(c(files[1], files[1]))),
    "unique"
  )
})

test_that("mismatched features across files are an error", {
  files <- write_files(c(
    S1 = cnt_table(tet_ids, 1:6),
    S2 = cnt_table(rev(tet_ids), 1:6),
    S3 = cnt_table(tet_ids[-1], 1:5)
  ))
  on.exit(unlink(files))
  expect_error(suppressMessages(TEtranscripts2DGEList(files[1:2])), "differ")
  expect_error(suppressMessages(TEtranscripts2DGEList(files[c(1, 3)])), "differ")
})

test_that("drop_zero removes all-zero features and keeps genes aligned", {
  files <- tet_files()
  on.exit(unlink(files))
  y0 <- suppressMessages(TEtranscripts2DGEList(files))
  y <- suppressMessages(TEtranscripts2DGEList(files, drop_zero = TRUE))

  expect_identical(
    rownames(y),
    c("ENSG01.1", "ENSG03.1", "AluY:Alu:SINE", "(CATTC)n:Satellite:Satellite")
  )
  expect_identical(rownames(y$genes), rownames(y$counts))
  expect_identical(y$genes$feature_id, rownames(y))
  expect_equal(y$samples$lib.size, y0$samples$lib.size)
  expect_equal(y$counts, y0$counts[rownames(y), ])
})

test_that("arguments are passed to DGEList", {
  files <- tet_files()
  on.exit(unlink(files))
  y <- suppressMessages(TEtranscripts2DGEList(files, group = c("a", "b")))
  expect_identical(as.character(y$samples$group), c("a", "b"))
})

test_that("gene_gtf and te_gtf annotate TEtranscripts features", {
  files <- tet_files()
  ggtf <- write_lines(gene_gtf_lines, ".gtf")
  tgtf <- write_lines(te_gtf_lines, ".gtf")
  tgz <- write_gz(te_gtf_lines)
  on.exit(unlink(c(files, ggtf, tgtf, tgz)))

  expect_message(
    y <- TEtranscripts2DGEList(files, gene_gtf = ggtf, te_gtf = tgtf),
    "1 genes not found"
  )
  g <- y$genes
  expect_identical(rownames(g), rownames(y$counts))
  expect_identical(g$gene_name, c("GENE1", "GENE2", NA, NA, NA, NA))
  expect_identical(g$gene_type, c("lncRNA", "protein_coding", NA, NA, NA, NA))
  expect_identical(g$seqnames, c("chr1", "chr2", NA, NA, NA, NA))
  expect_equal(g$start, c(100, 900, NA, NA, NA, NA))
  expect_equal(g$end, c(500, 1500, NA, NA, NA, NA))
  expect_identical(g$strand, c("+", "-", NA, NA, NA, NA))
  # AluY:Alu:SINE has 2 loci (300 + 100 bp); AluY:Other:SINE is separate
  expect_identical(g$te_n_loci, c(NA, NA, NA, 2L, NA, 1L))
  expect_equal(g$te_total_width, c(NA, NA, NA, 400, NA, 50))

  # gzipped TE GTF gives the same result, also when read in small chunks
  y_gz <- suppressMessages(TEtranscripts2DGEList(files, gene_gtf = ggtf, te_gtf = tgz))
  expect_identical(y_gz$genes, g)
  a <- te_read_te_gtf(tgz)
  b <- te_read_te_gtf(tgtf)
  expect_equal(a[order(a$te_key), ], b[order(b$te_key), ])
  chunked <- te_read_gtf(tgz, FUN = function(dt) dt, chunk_size = 1L)
  expect_identical(nrow(chunked), 4L)

  # gene annotation alone adds no TE columns
  y_gene <- suppressMessages(TEtranscripts2DGEList(files, gene_gtf = ggtf))
  expect_false("te_n_loci" %in% names(y_gene$genes))
})

test_that("TElocal locations add per-locus coordinates", {
  files <- tel_files()
  ggtf <- write_lines(gene_gtf_lines, ".gtf")
  loc <- write_lines(loc_lines, ".locations")
  tgtf <- write_lines(te_gtf_lines, ".gtf")
  on.exit(unlink(c(files, ggtf, loc, tgtf)))

  y <- suppressMessages(TEtranscripts2DGEList(files, gene_gtf = ggtf, te_gtf = loc))
  g <- y$genes
  expect_identical(g$seqnames, c("chr1", "chr2", NA, "chr1", "chr2", "chr3"))
  expect_equal(g$start, c(100, 900, NA, 1001, 2001, 51))
  expect_equal(g$end, c(500, 1500, NA, 1300, 2100, 100))
  expect_identical(g$strand, c("+", "-", NA, "+", "-", "+"))
  expect_identical(g$gene_name[4:6], rep(NA_character_, 3))

  # locations alone, with drop_zero
  y2 <- suppressMessages(TEtranscripts2DGEList(files, te_gtf = loc, drop_zero = TRUE))
  expect_identical(y2$genes$te_locus, c(NA, NA, "AluY_dup1", "(CATTC)n"))
  expect_identical(y2$genes$seqnames, c(NA, NA, "chr1", "chr3"))

  # a TE GTF gives subfamily-level summaries for TElocal
  y3 <- suppressMessages(TEtranscripts2DGEList(files, te_gtf = tgtf))
  expect_identical(y3$genes$te_n_loci, c(NA, NA, NA, 2L, 2L, 1L))

  # locations cannot be used for TEtranscripts
  tet <- tet_files()
  on.exit(unlink(tet), add = TRUE)
  expect_error(
    suppressMessages(TEtranscripts2DGEList(tet, te_gtf = loc)),
    "locations file"
  )
})

test_that("invalid inputs give informative errors", {
  files <- tet_files()
  ggtf <- write_lines(gene_gtf_lines, ".gtf")
  tgtf <- write_lines(te_gtf_lines, ".gtf")
  csv <- write_lines(c("a,b", "1,2"), ".csv")
  bad <- write_files(list("gene/TE\tA.bam\tB.bam\nG1\t1\t2\n"))
  genes_only <- write_files(list(cnt_table(tet_ids[1:3], 1:3)))
  on.exit(unlink(c(files, ggtf, tgtf, csv, bad, genes_only)))

  expect_error(TEtranscripts2DGEList("does_not_exist.cntTable"), "not found")
  expect_error(TEtranscripts2DGEList(character(0)), "non-empty")
  expect_error(TEtranscripts2DGEList(files, drop_zero = NA), "drop_zero")
  expect_error(TEtranscripts2DGEList(files, drop_zero = "yes"), "drop_zero")
  expect_error(TEtranscripts2DGEList(files, gene_gtf = c(ggtf, ggtf)), "single file path")
  expect_error(TEtranscripts2DGEList(files, te_gtf = 1), "single file path")
  expect_error(TEtranscripts2DGEList(files, gene_gtf = "missing.gtf"), "not found")
  expect_error(TEtranscripts2DGEList(files, gene_gtf = csv), "gene GTF")
  expect_error(TEtranscripts2DGEList(files, gene_gtf = tgtf), "gene GTF")
  expect_error(TEtranscripts2DGEList(files, te_gtf = csv), "TE GTF")
  expect_error(TEtranscripts2DGEList(files, te_gtf = ggtf), "TE GTF")
  expect_error(suppressMessages(TEtranscripts2DGEList(bad)), "count table")
  expect_error(suppressMessages(TEtranscripts2DGEList(genes_only)), "No TE features")
  expect_no_error(suppressMessages(TEtranscripts2DGEList(genes_only, type = "TElocal")))
})

test_that("fractional counts are stored as doubles", {
  files <- write_files(c(
    S1 = cnt_table(tet_ids, c(1.5, 0, 3, 5, 0, 1)),
    S2 = cnt_table(tet_ids, 1:6)
  ))
  on.exit(unlink(files))
  y <- suppressMessages(TEtranscripts2DGEList(files))
  expect_true(is.double(y$counts))
  expect_equal(unname(y$counts[1, ]), c(1.5, 1))
})

test_that("gtf_attr extracts attributes and returns NA when absent", {
  x <- c(
    'gene_id "G1"; gene_name "A"; level 2;',
    'gene_id "G2"; havana_gene_id "X";',
    'gene_id "G3"; gene_name "";'
  )
  expect_identical(gtf_attr(x, "gene_id"), c("G1", "G2", "G3"))
  expect_identical(gtf_attr(x, "gene_name"), c("A", NA, ""))
})
