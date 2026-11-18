#' Import TEtranscripts or TElocal Count Tables into a DGEList
#'
#' Read the \code{*.cntTable} files produced by TEtranscripts, TEcount or
#' TElocal into an \code{edgeR::DGEList}. The tool that produced the files is
#' detected from the format of the TE feature IDs. Information encoded in the
#' feature IDs is always added to the \code{genes} slot, and gene and TE
#' annotations can optionally be added from the same annotation files used to
#' run TEtranscripts/TElocal (\code{--GTF} and \code{--TE}).
#'
#' @details
#' TEtranscripts/TEcount report TE counts per subfamily with IDs of the form
#' \code{subfamily:family:class} (e.g. \code{AluY:Alu:SINE}). TElocal reports
#' counts per TE locus with IDs of the form
#' \code{locus:subfamily:family:class} (e.g. \code{AluY_dup12:AluY:Alu:SINE}).
#' Gene IDs (quoted in the raw tables) are returned without quotes.
#'
#' The \code{genes} slot always contains the columns \code{feature_id},
#' \code{feature_type} ("gene" or "TE"), \code{te_locus} (TElocal only),
#' \code{te_subfamily}, \code{te_family} and \code{te_class}. Additional columns
#' are added when annotation files are supplied:
#' \itemize{
#'   \item \code{gene_gtf}: \code{gene_name}, \code{gene_type} (or
#'   \code{gene_biotype}), \code{seqnames}, \code{start}, \code{end} and
#'   \code{strand} for gene features, taken from the "gene" rows of the GTF.
#'   \item \code{te_gtf} as a TElocal locations table
#'   (\code{*.locInd.locations}): \code{seqnames}, \code{start}, \code{end}
#'   and \code{strand} for each TE locus, reported as given in the file.
#'   \item \code{te_gtf} as a TE GTF (e.g. \code{GRCh38_GENCODE_rmsk_TE.gtf}):
#'   \code{te_n_loci}, the number of genomic copies of the subfamily, and
#'   \code{te_total_width}, their summed width in bp. These are summarized
#'   per \code{subfamily:family:class}, so for TElocal data they describe the
#'   subfamily of each locus.
#' }
#' Features absent from an annotation file are left as \code{NA} and the number
#' of unmatched features is reported with a message.
#'
#' TElocal tables contain several million mostly zero-count features. Setting
#' \code{drop_zero = TRUE} removes features with zero counts in every sample
#' before the annotation is added, which greatly reduces memory use. Library
#' sizes are unaffected. Reading a full GENCODE GTF can take about a minute.
#' Annotation files may be gzip-compressed.
#'
#' @param files Character vector of paths to \code{.cntTable} files. Sample
#'   names are taken from \code{names(files)} if present, otherwise from the
#'   file names with the \code{.cntTable} extension removed.
#' @param gene_gtf \code{NULL} (default) or the path to the gene GTF used to
#'   generate the counts (the \code{--GTF} argument).
#' @param te_gtf \code{NULL} (default) or the path to TE annotation (the
#'   \code{--TE} argument). For TEtranscripts this is the TE GTF. For TElocal
#'   this is the locations table matching the \code{.locInd} index, which gives
#'   per-locus coordinates; a TE GTF is also accepted. The file type is
#'   detected from its content.
#' @param type One of "auto" (default), "TEtranscripts" or "TElocal". "auto"
#'   detects the type from the TE feature IDs.
#' @param drop_zero Logical. Remove features with zero counts in all samples?
#'   Default \code{FALSE}.
#' @param ... Additional arguments passed to \code{edgeR::DGEList()}, e.g.
#'   \code{group}.
#' @return A \code{DGEList} with one row per feature and one column per file.
#'   The \code{samples} slot gains a \code{bam} column with the BAM file named
#'   in each table's header.
#' @export
#' @examples
#' \dontrun{
#' count_files <- list.files(
#'   path = "/path/to/TEtranscripts_out",
#'   pattern = "*.cntTable",
#'   recursive  = TRUE,
#'   full.names = TRUE
#' )
#'
#' y <-  TEtranscripts2DGEList(
#'   count_files,
#'   gene_gtf = "gencode.v48.basic.annotation.gtf",
#'   te_gtf = "GRCh38_GENCODE_rmsk_TE.gtf.gz",
#'   type = "TEtranscripts"
#' )
#'
#' }
#'
TEtranscripts2DGEList <- function(
  files,
  gene_gtf = NULL,
  te_gtf = NULL,
  type = c("auto", "TEtranscripts", "TElocal"),
  drop_zero = TRUE,
  ...
) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop(
      "edgeR package is required. Install it with BiocManager::install('edgeR')."
    )
  }
  type <- match.arg(type)
  if (!is.character(files) || length(files) == 0L) {
    stop("files must be a non-empty character vector of file paths.")
  }
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0L) {
    stop("File(s) not found: ", paste(missing_files, collapse = ", "))
  }
  if (!is.logical(drop_zero) || length(drop_zero) != 1L || is.na(drop_zero)) {
    stop("drop_zero must be TRUE or FALSE.")
  }
  te_check_path(gene_gtf, "gene_gtf")
  te_check_path(te_gtf, "te_gtf")

  # Check annotation file formats before the (slow) count reading
  if (!is.null(gene_gtf) && te_annotation_kind(gene_gtf) != "gene_gtf") {
    stop("gene_gtf does not look like a gene GTF file: ", gene_gtf)
  }
  te_kind <- if (is.null(te_gtf)) NULL else te_annotation_kind(te_gtf)
  if (!is.null(te_kind) && !te_kind %in% c("te_gtf", "te_locations")) {
    stop(
      "te_gtf does not look like a TE GTF or TElocal locations file: ",
      te_gtf
    )
  }

  samples <- te_sample_names(files)

  message("Reading in count tables...")
  cnt <- te_read_counts(files)
  colnames(cnt$counts) <- samples

  type <- te_detect_type(cnt$ids[!cnt$is_gene], type)
  message("Detected ", type, " count tables.")
  if (identical(te_kind, "te_locations") && type != "TElocal") {
    stop(
      "te_gtf is a TElocal locations file but the counts are from ",
      type,
      ". Supply the TE GTF instead."
    )
  }

  if (drop_zero) {
    keep <- rowSums(cnt$counts) > 0
    message(
      "Removing ",
      sum(!keep),
      " features with zero counts in all samples..."
    )
    cnt$counts <- cnt$counts[keep, , drop = FALSE]
    cnt$ids <- cnt$ids[keep]
    cnt$is_gene <- cnt$is_gene[keep]
  }
  rownames(cnt$counts) <- cnt$ids

  genes <- te_parse_ids(cnt$ids, cnt$is_gene, type)
  if (!is.null(gene_gtf) || !is.null(te_gtf)) {
    genes <- te_annotate(genes, gene_gtf, te_gtf, te_kind)
  }
  rownames(genes) <- cnt$ids
  # DGEList() overwrites rownames(genes) without checking them
  stopifnot(identical(rownames(genes), rownames(cnt$counts)))

  y <- edgeR::DGEList(counts = cnt$counts, genes = genes, ...)
  y$samples$bam <- cnt$bam
  message("Done.")

  y
}

# Validate an optional single file path argument
te_check_path <- function(x, arg) {
  if (is.null(x)) {
    return(invisible(NULL))
  }
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(arg, " must be NULL or a single file path.")
  }
  if (!file.exists(x)) {
    stop(arg, " file not found: ", x)
  }
  invisible(NULL)
}

# Sample names from names(files) or file names without the .cntTable extension
te_sample_names <- function(files) {
  nms <- names(files)
  if (is.null(nms) || any(is.na(nms) | nms == "")) {
    nms <- sub("\\.cntTable$", "", basename(files))
  }
  if (anyDuplicated(nms)) {
    stop(
      "Sample names must be unique. Duplicated: ",
      paste(unique(nms[duplicated(nms)]), collapse = ", "),
      ". Supply a named vector of files."
    )
  }
  unname(nms)
}

# Read count tables into a features x samples matrix. Gene IDs are quoted in
# the raw tables, which distinguishes them from TE IDs.
te_read_counts <- function(files) {
  bam <- character(length(files))
  raw_ids <- NULL
  m <- NULL
  for (i in seq_along(files)) {
    f <- files[[i]]
    header <- strsplit(readLines(f, n = 1L, warn = FALSE), "\t", fixed = TRUE)
    if (length(header) == 0L || length(header[[1]]) != 2L) {
      stop(
        "'",
        f,
        "' is not a TEtranscripts/TElocal count table ",
        "(expected a 2 column header line, e.g. 'gene/TE<TAB>sample.bam')."
      )
    }
    bam[i] <- header[[1]][2]
    dt <- data.table::fread(
      f,
      skip = 1L,
      header = FALSE,
      sep = "\t",
      quote = "",
      col.names = c("id", "count"),
      colClasses = list(character = 1L),
      integer64 = "double",
      showProgress = FALSE
    )
    if (nrow(dt) == 0L) {
      stop("'", f, "' contains no features.")
    }
    if (!is.numeric(dt$count) || anyNA(dt$count)) {
      stop("'", f, "' contains missing or non-numeric counts.")
    }
    if (is.null(raw_ids)) {
      raw_ids <- dt$id
      if (anyDuplicated(raw_ids)) {
        stop("'", f, "' contains duplicated feature IDs.")
      }
      m <- matrix(0L, nrow = length(raw_ids), ncol = length(files))
    } else if (!identical(dt$id, raw_ids)) {
      stop(
        "Features in '",
        f,
        "' differ from those in '",
        files[[1]],
        "'. All count tables must be generated with the same annotation."
      )
    }
    # Assigning doubles promotes the whole matrix to double
    m[, i] <- dt$count
  }
  ids <- gsub("\"", "", raw_ids, fixed = TRUE)
  is_gene <- startsWith(raw_ids, "\"") | !grepl(":", ids, fixed = TRUE)
  list(counts = m, ids = ids, is_gene = is_gene, bam = bam)
}

# Detect TEtranscripts (subfamily:family:class) vs TElocal
# (locus:subfamily:family:class) from the number of ':' separated fields
te_detect_type <- function(te_ids, type) {
  if (length(te_ids) == 0L) {
    if (type == "auto") {
      stop("No TE features found in the count tables. Specify type explicitly.")
    }
    return(type)
  }
  n_fields <- nchar(te_ids) - nchar(gsub(":", "", te_ids, fixed = TRUE)) + 1L
  detected <- if (all(n_fields == 3L)) {
    "TEtranscripts"
  } else if (all(n_fields == 4L)) {
    "TElocal"
  } else {
    stop(
      "Unrecognized TE feature IDs (e.g. '",
      te_ids[!n_fields %in% c(3L, 4L)][1],
      "'). Expected 'subfamily:family:class' (TEtranscripts) or ",
      "'locus:subfamily:family:class' (TElocal) for every TE."
    )
  }
  if (type != "auto" && type != detected) {
    stop(
      "type = '",
      type,
      "' but the TE feature IDs look like ",
      detected,
      " output."
    )
  }
  detected
}

# Build the genes data.frame from the feature IDs
te_parse_ids <- function(ids, is_gene, type) {
  n <- length(ids)
  te <- !is_gene
  parts <- if (any(te)) {
    tstrsplit(ids[te], ":", fixed = TRUE)
  } else {
    rep(list(character(0)), if (type == "TElocal") 4L else 3L)
  }
  fill <- function(x) {
    out <- rep(NA_character_, n)
    out[te] <- x
    out
  }
  genes <- data.frame(
    feature_id = ids,
    feature_type = ifelse(is_gene, "gene", "TE"),
    stringsAsFactors = FALSE
  )
  if (type == "TElocal") {
    genes$te_locus <- fill(parts[[1]])
    parts <- parts[-1]
  }
  genes$te_subfamily <- fill(parts[[1]])
  genes$te_family <- fill(parts[[2]])
  genes$te_class <- fill(parts[[3]])
  genes
}

# Detect an annotation file's type from its first data lines
te_annotation_kind <- function(path) {
  con <- file(path, "r")
  on.exit(close(con))
  lines <- readLines(con, n = 100L, warn = FALSE)
  lines <- lines[!startsWith(lines, "#") & nzchar(lines)]
  if (length(lines) == 0L) {
    return("unknown")
  }
  fields <- strsplit(lines[seq_len(min(5L, length(lines)))], "\t", fixed = TRUE)
  loc_pattern <- "^.+:[0-9]+-[0-9]+:[-+.]$"
  is_loc <- vapply(
    fields,
    function(x) length(x) == 2L && grepl(loc_pattern, x[2]),
    logical(1)
  )
  if (any(is_loc)) {
    return("te_locations")
  }
  first <- fields[[1]]
  if (length(first) == 9L) {
    if (grepl("class_id \"", first[9], fixed = TRUE)) {
      return("te_gtf")
    }
    if (grepl("gene_id \"", first[9], fixed = TRUE)) {
      return("gene_gtf")
    }
  }
  "unknown"
}

# Extract a GTF attribute value; NA where the attribute is absent
gtf_attr <- function(x, key) {
  out <- rep(NA_character_, length(x))
  pattern <- paste0("(^|;)\\s*", key, " \"")
  hit <- grepl(pattern, x, perl = TRUE)
  out[hit] <- sub(
    paste0("^.*?(?:^|;)\\s*", key, " \"([^\"]*)\".*$"),
    "\\1",
    x[hit],
    perl = TRUE
  )
  out
}

# Read the columns seqnames, feature, start, end, strand and attr of a GTF,
# optionally keeping only rows of one feature type, and apply FUN to the
# result. Gzipped files are streamed in chunks (fread() needs R.utils for
# them) and FUN is applied per chunk.
te_read_gtf <- function(
  path,
  feature = NULL,
  FUN = identity,
  chunk_size = 5e5L
) {
  cols <- c("seqnames", "feature", "start", "end", "strand", "attr")
  read_text <- function(lines) {
    fread(
      text = lines,
      sep = "\t",
      quote = "",
      header = FALSE,
      select = c(1L, 3L, 4L, 5L, 7L, 9L),
      col.names = cols,
      showProgress = FALSE
    )
  }
  keep_feature <- function(dt) {
    if (!is.null(feature)) {
      dt <- dt[dt$feature == feature, ]
    }
    dt$seqnames <- as.character(dt$seqnames)
    dt$start <- as.numeric(dt$start)
    dt$end <- as.numeric(dt$end)
    dt
  }

  if (!grepl("\\.gz$", path, ignore.case = TRUE)) {
    head_lines <- readLines(path, n = 1000L, warn = FALSE)
    n_skip <- which(!startsWith(head_lines, "#"))[1] - 1L
    if (is.na(n_skip)) {
      n_skip <- length(head_lines)
    }
    dt <- fread(
      path,
      skip = n_skip,
      sep = "\t",
      quote = "",
      header = FALSE,
      select = c(1L, 3L, 4L, 5L, 7L, 9L),
      col.names = cols,
      showProgress = FALSE
    )
    return(FUN(keep_feature(dt)))
  }

  con <- gzfile(path, "r")
  on.exit(close(con))
  feature_pattern <- if (is.null(feature)) {
    NULL
  } else {
    paste0("^[^\t]*\t[^\t]*\t", feature, "\t")
  }
  out <- list()
  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (length(lines) == 0L) {
      break
    }
    lines <- lines[!startsWith(lines, "#") & nzchar(lines)]
    if (!is.null(feature_pattern)) {
      lines <- lines[grepl(feature_pattern, lines, perl = TRUE)]
    }
    if (length(lines) > 0L) {
      out[[length(out) + 1L]] <- FUN(keep_feature(read_text(lines)))
    }
  }
  rbindlist(out)
}

# Gene-level annotation from a gene GTF
te_read_gene_gtf <- function(path) {
  dt <- te_read_gtf(path, feature = "gene")
  gene_type <- gtf_attr(dt$attr, "gene_type")
  if (all(is.na(gene_type))) {
    gene_type <- gtf_attr(dt$attr, "gene_biotype")
  }
  out <- data.frame(
    gene_id = gtf_attr(dt$attr, "gene_id"),
    gene_name = gtf_attr(dt$attr, "gene_name"),
    gene_type = gene_type,
    seqnames = dt$seqnames,
    start = dt$start,
    end = dt$end,
    strand = dt$strand,
    stringsAsFactors = FALSE
  )
  out[!is.na(out$gene_id) & !duplicated(out$gene_id), ]
}

# Number of loci and total width per subfamily:family:class from a TE GTF
te_read_te_gtf <- function(path) {
  summarize_chunk <- function(dt) {
    agg <- data.table(
      te_key = paste(
        gtf_attr(dt$attr, "gene_id"),
        gtf_attr(dt$attr, "family_id"),
        gtf_attr(dt$attr, "class_id"),
        sep = ":"
      ),
      w = dt$end - dt$start + 1
    )
    agg[, list(te_n_loci = .N, te_total_width = sum(w)), by = "te_key"]
  }
  dt <- te_read_gtf(path, FUN = summarize_chunk)
  # Re-aggregate subfamilies split across chunks
  dt <- dt[,
    list(te_n_loci = sum(te_n_loci), te_total_width = sum(te_total_width)),
    by = "te_key"
  ]
  dt$te_n_loci <- as.integer(dt$te_n_loci)
  dt
}

# Per-locus coordinates from a TElocal locations table
te_read_locations <- function(path, loci) {
  dt <- fread(
    path,
    sep = "\t",
    header = FALSE,
    quote = "",
    col.names = c("locus", "loc"),
    colClasses = "character",
    showProgress = FALSE
  )
  dt <- dt[dt$locus %chin% loci, ]
  pattern <- "^(.*):([0-9]+)-([0-9]+):([-+.])$"
  ok <- grepl(pattern, dt$loc)
  dt <- dt[ok, ]
  data.frame(
    locus = dt$locus,
    seqnames = sub(pattern, "\\1", dt$loc),
    start = as.numeric(sub(pattern, "\\2", dt$loc)),
    end = as.numeric(sub(pattern, "\\3", dt$loc)),
    strand = sub(pattern, "\\4", dt$loc),
    stringsAsFactors = FALSE
  )
}

# Add annotation columns to genes, preserving row order
te_annotate <- function(genes, gene_gtf, te_gtf, te_kind) {
  is_gene <- genes$feature_type == "gene"
  is_te <- !is_gene
  coords <- c("seqnames", "start", "end", "strand")
  add_coords <- function(genes) {
    if (!"seqnames" %in% names(genes)) {
      genes$seqnames <- NA_character_
      genes$start <- NA_real_
      genes$end <- NA_real_
      genes$strand <- NA_character_
    }
    genes
  }
  report <- function(n, what, file) {
    if (n > 0L) {
      message(n, " ", what, " not found in ", basename(file), ".")
    }
  }

  if (!is.null(gene_gtf)) {
    message("Reading gene annotation from ", basename(gene_gtf), "...")
    g <- te_read_gene_gtf(gene_gtf)
    idx <- match(genes$feature_id, g$gene_id)
    idx[is_te] <- NA_integer_
    genes$gene_name <- g$gene_name[idx]
    genes$gene_type <- g$gene_type[idx]
    genes <- add_coords(genes)
    for (col in coords) {
      genes[[col]][is_gene] <- g[[col]][idx[is_gene]]
    }
    report(sum(is.na(idx[is_gene])), "genes", gene_gtf)
  }

  if (identical(te_kind, "te_locations")) {
    message("Reading TE locus coordinates from ", basename(te_gtf), "...")
    loc <- te_read_locations(te_gtf, genes$te_locus[is_te])
    idx <- match(genes$te_locus, loc$locus)
    idx[is_gene] <- NA_integer_
    genes <- add_coords(genes)
    for (col in coords) {
      genes[[col]][is_te] <- loc[[col]][idx[is_te]]
    }
    report(sum(is.na(idx[is_te])), "TE loci", te_gtf)
  } else if (identical(te_kind, "te_gtf")) {
    message("Reading TE annotation from ", basename(te_gtf), "...")
    te <- te_read_te_gtf(te_gtf)
    key <- paste(genes$te_subfamily, genes$te_family, genes$te_class, sep = ":")
    idx <- match(key, te$te_key)
    idx[is_gene] <- NA_integer_
    genes$te_n_loci <- te$te_n_loci[idx]
    genes$te_total_width <- te$te_total_width[idx]
    report(sum(is.na(idx[is_te])), "TE features", te_gtf)
  }

  genes
}
