json_resp <- function(body, status = 200) {
  httr2::response(
    status_code = status,
    headers = list(`Content-Type` = "application/json"),
    body = charToRaw(body)
  )
}

# Decoded form fields of a captured request
form_of <- function(req) {
  vapply(req$body$data, function(v) utils::URLdecode(unclass(v)), character(1))
}

ok <- '{"results":{"result":[{"number_in_list":1,"term":{"id":"X","label":"x"}}]}}'

test_that("Multiple results are parsed to one row per term", {
  skip_if_not_installed("httr2")
  multi <- '{"results":{"result":[
    {"number_in_list":5,"fold_enrichment":10.5,"fdr":0.001,"expected":0.47,
     "number_in_reference":120,"pValue":1e-5,"plus_minus":"+",
     "term":{"id":"GO:0016055","label":"Wnt signaling pathway"}},
    {"number_in_list":1,"fold_enrichment":0.2,"fdr":1,"expected":5,
     "number_in_reference":2000,"pValue":0.5,"plus_minus":"-",
     "term":{"label":"UNCLASSIFIED"}}]}}'
  httr2::local_mocked_responses(function(req) json_resp(multi))
  res <- panther_go(c("AXIN1", "AXIN2"), "9606", "biological_process")
  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), 2L)
  expect_identical(res$term_id, c("GO:0016055", NA))
  expect_identical(res$term_label, c("Wnt signaling pathway", "UNCLASSIFIED"))
  expect_false("term" %in% names(res))
  expect_equal(res$fold_enrichment, c(10.5, 0.2))
  expect_equal(res$number_in_list, c(5L, 1L))
  expect_equal(res$pValue, c(1e-5, 0.5))
})

test_that("A single result object is returned as one row", {
  skip_if_not_installed("httr2")
  single <- '{"results":{"result":{"number_in_list":5,"fold_enrichment":10.5,
    "fdr":0.001,"pValue":1e-5,"plus_minus":"+",
    "term":{"id":"GO:0016055","label":"Wnt signaling pathway"}}}}'
  httr2::local_mocked_responses(function(req) json_resp(single))
  res <- panther_go("AXIN1", "9606", "biological_process")
  expect_equal(nrow(res), 1L)
  expect_identical(res$term_id, "GO:0016055")
  expect_identical(res$term_label, "Wnt signaling pathway")
})

test_that("Request is encoded with the PANTHER codes for every option", {
  skip_if_not_installed("httr2")
  map <- list(
    biological_process = "GO:0008150",
    molecular_function = "GO:0003674",
    cellular_component = "GO:0005575",
    panther_go_slim_mf = "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF",
    panther_go_slim_bp = "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP",
    panther_go_slim_cc = "ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC",
    panther_pc = "ANNOT_TYPE_ID_PANTHER_PC",
    panther_pathway = "ANNOT_TYPE_ID_PANTHER_PATHWAY",
    panther_reactome_pathway = "ANNOT_TYPE_ID_REACTOME_PATHWAY"
  )
  for (ds in names(map)) {
    for (tt in c("fisher", "binomial")) {
      for (cr in c("fdr", "bonferroni", "none")) {
        captured <- NULL
        httr2::with_mocked_responses(
          function(req) {
            captured <<- req
            json_resp(ok)
          },
          panther_go(c("A", "B", "C"), "10090", ds,
            enrichment_test_type = tt, correction = cr
          )
        )
        f <- form_of(captured)
        expect_identical(
          names(f),
          c("geneInputList", "organism", "annotDataSet", "enrichmentTestType", "correction")
        )
        expect_identical(
          unname(f),
          c("A,B,C", "10090", map[[ds]], toupper(tt), toupper(cr))
        )
        expect_identical(captured$method, "POST")
      }
    }
  }
  expect_identical(captured$headers$Accept, "application/json")
})

test_that("A reference list adds refInputList and refOrganism", {
  skip_if_not_installed("httr2")
  captured <- NULL
  httr2::local_mocked_responses(function(req) {
    captured <<- req
    json_resp(ok)
  })
  panther_go(c("A", "B"), "9606", "panther_pathway", ref_input_list = c("A", "B", "C", "D"))
  f <- form_of(captured)
  expect_identical(f[["refInputList"]], "A,B,C,D")
  expect_identical(f[["refOrganism"]], "9606")
})

test_that("PANTHER errors are reported with genes missing from the reference", {
  skip_if_not_installed("httr2")
  err <- '{"search":{"error":"Some ids not in reference"}}'
  httr2::local_mocked_responses(function(req) json_resp(err))
  expect_error(
    panther_go(c("A", "B", "Z"), "9606", "biological_process", ref_input_list = c("A", "B")),
    "PANTHER returned an error: Some ids not in reference.*1 not found in ref_input_list: Z"
  )
  expect_error(panther_go("A", "9606", "biological_process"), "PANTHER returned an error: Some ids")
  e <- tryCatch(
    panther_go(paste0("G", 1:12), "9606", "biological_process", ref_input_list = "Q"),
    error = function(e) conditionMessage(e)
  )
  expect_match(e, "12 not found in ref_input_list: G1, G2, G3, G4, G5, G6, G7, G8, G9, G10, ...", fixed = TRUE)
  expect_false(grepl("G11", e))
})

test_that("HTTP errors and invalid arguments raise errors", {
  skip_if_not_installed("httr2")
  httr2::local_mocked_responses(function(req) json_resp("{}", status = 500))
  expect_error(panther_go("A", "9606", "biological_process"))
  expect_error(panther_go(character(), "9606", "biological_process"), "Empty gene_list")
  expect_error(panther_go("A", "9606", "go_bp"), "Incorrectly specified dataset")
})

test_that("Live PANTHER results satisfy statistical invariants", {
  skip_if_not_installed("httr2")
  skip_on_cran()
  skip_if_offline("pantherdb.org")
  genes <- c(
    "AXIN1", "AXIN2", "CTNNB1", "DKK1", "FZD1", "FZD8", "LEF1", "TCF7",
    "WNT1", "WNT5B", "WNT6", "DVL2", "NKD1", "FRAT1"
  )
  res <- panther_go(genes, "9606", "panther_pathway")
  expect_true(all(c(
    "term_id", "term_label", "number_in_list", "number_in_reference",
    "expected", "fold_enrichment", "pValue", "fdr", "plus_minus"
  ) %in% names(res)))
  # PANTHER can report p-values of 1 + ~1e-10
  expect_true(all(res$pValue >= 0 & res$pValue <= 1 + 1e-6))
  expect_true(all(res$fdr >= res$pValue - 1e-6))
  expect_true(all(res$number_in_list <= length(unique(genes))))
  fe <- res[res$expected > 0 & res$number_in_list > 0, ]
  expect_equal(fe$fold_enrichment, fe$number_in_list / fe$expected, tolerance = 0.01)
  up <- res[res$plus_minus == "+", ]
  expect_true(all(up$number_in_list > up$expected))
  # Wnt signaling pathway is the top term
  expect_identical(res$term_id[which.min(res$pValue)], "P00057")
})
