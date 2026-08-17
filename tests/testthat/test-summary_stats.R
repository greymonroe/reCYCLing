# Regression tests for the summary-statistics catalog (schema g6).

sim_small <- function(seed, ...) {
  set.seed(seed)
  sim_genome(K = 3L, target_sizes = rep(250L, 3L), n_generations = 300L,
             mu_total = 1e-4, ...)
}

test_that("genome catalog is fixed-length, named, and versioned", {
  st <- repeat_genome_stats(sim_small(1))
  expect_identical(st$version, "g6")
  p <- repeat_genome_predictors(st)
  expect_identical(length(p), 87L)
  expect_true(all(nzchar(names(p))))
  expect_identical(sum(!is.na(p[grep("^gdv", names(p))])), 11L)
})

test_that("reduced ABC set selects the three families plus size context", {
  p  <- repeat_genome_predictors(repeat_genome_stats(sim_small(2)))
  ab <- repeat_abc_predictors(p)
  expect_identical(length(ab), 39L)
  expect_equal(sum(ab[grep("^gdv", names(ab))]), 1, tolerance = 1e-9)
  expect_equal(sum(ab[grep("^wxd", names(ab))]), 1, tolerance = 1e-9)
  expect_true(all(c("frac_pairs_xchrom", "frac_genome_xchrom",
                    "log_array_size_m", "log_n_chrom") %in% names(ab)))
})

test_that("no-exchange genomes score genuine zero derived sharing (not NA)", {
  p <- repeat_genome_predictors(repeat_genome_stats(sim_small(3, p_translocation = 0)))
  expect_identical(unname(p["frac_pairs_xchrom"]), 0)
  expect_identical(unname(p["frac_genome_xchrom"]), 0)
  expect_identical(unname(p["log_n_xchrom_fam"]), -12)  # absence sentinel
})

test_that("exchange raises the cross-chromosome sharing family", {
  p0 <- repeat_genome_predictors(repeat_genome_stats(sim_small(4, p_translocation = 0)))
  p1 <- repeat_genome_predictors(repeat_genome_stats(
    sim_small(4, p_translocation = 0.3, transloc_chunk_mean = 100)))
  expect_gt(p1["frac_genome_xchrom"], p0["frac_genome_xchrom"])
  expect_gt(p1["frac_pairs_xchrom"], 0)
})

test_that("divergence filter and clock behave: older genome is more diverged", {
  young <- repeat_genome_stats(sim_small(5))
  set.seed(5)
  older <- repeat_genome_stats(sim_genome(K = 3L, target_sizes = rep(250L, 3L),
                                        n_generations = 900L, mu_total = 1e-4))
  expect_gt(older$mean_div_bp_m, young$mean_div_bp_m)
})
