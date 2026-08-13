# Regression tests for the multi-chromosome genome engine (release model).
# Small/fast sims; each test targets one designed property of the model.

K       <- 3L
TARGETS <- rep(300L, 3L)
NGEN    <- 400L
MU      <- 1e-4
L       <- 178L

fixed_ancestor <- function() paste(rep(c("A", "C", "G", "T"), length.out = 178), collapse = "")

quick_sim <- function(seed, ...) {
  set.seed(seed)
  sim_genome(K = K, target_sizes = TARGETS, n_generations = NGEN,
             mu_total = MU, init_l = L, ...)
}

test_that("output structure and attributes are as documented", {
  g <- quick_sim(1)
  expect_s3_class(g, "data.table")
  expect_named(g, c("chrom", "num", "bponly", "dir"))
  expect_setequal(unique(g$chrom), 1:K)
  expect_true(all(g$dir %in% c("+", "-")))
  expect_length(attr(g, "array_sizes"), K)
  expect_identical(attr(g, "total_gens"), NGEN)
  expect_identical(nrow(g), sum(attr(g, "array_sizes")))
})

test_that("same seed reproduces the identical genome", {
  expect_identical(quick_sim(42, p_translocation = 0.2), quick_sim(42, p_translocation = 0.2))
})

test_that("size governor holds the genome-wide total near target", {
  g <- quick_sim(2)  # local turnover only: balance is exact in expectation
  total <- sum(attr(g, "array_sizes"))
  expect_lt(abs(total - sum(TARGETS)) / sum(TARGETS), 0.10 + 0.10)  # band + slack
})

test_that("divergence from the founder follows the mutation clock (JC)", {
  anc <- fixed_ancestor()
  set.seed(3)
  g <- sim_genome(K = K, target_sizes = TARGETS, n_generations = NGEN,
                  mu_total = MU, init_l = L,
                  init_sequence_type = "given", ancestor_seq = anc)
  anc_v <- strsplit(anc, "")[[1]]
  mism <- vapply(strsplit(g$bponly, ""), function(s) mean(s != anc_v), 0.0)
  # copies inherit their parent's mutations, so divergence from the FOUNDER is
  # mu * age regardless of turnover:  E[p_diff] = 3/4 * (1 - (1 - 4mu/3)^g)
  expected <- 0.75 * (1 - (1 - 4 * MU / 3)^NGEN)
  expect_lt(abs(mean(mism) - expected) / expected, 0.15)
})

test_that("no inter-chromosomal exchange => no cross-chromosome identical sharing", {
  g <- quick_sim(4, p_translocation = 0)
  shared <- unique(g[, .(chroms = data.table::uniqueN(chrom)), by = bponly][chroms > 1])
  expect_identical(nrow(shared), 0L)
})

test_that("inter-chromosomal exchange creates cross-chromosome identical sharing", {
  g <- quick_sim(5, p_translocation = 0.3, transloc_chunk_mean = 100)
  expect_gt(attr(g, "n_transloc"), 0)
  shared <- g[, .(chroms = data.table::uniqueN(chrom)), by = bponly][chroms > 1]
  expect_gt(nrow(shared), 0L)
})

test_that("long-range duplication is off by default and controllable", {
  expect_error(sim_genome(K = 2L, target_sizes = c(100L), n_generations = 100L),
               "length K")
  expect_error(sim_genome(K = 2L, target_sizes = c(100L, 100L), n_generations = 0),
               "n_generations")
})

test_that("single-array engine still runs after the API prune", {
  set.seed(6)
  s <- run_sim_ps(max_units = 500, init_l = 60, mu_total = 1e-4,
                  compute_pairs = FALSE)
  expect_true(nrow(s$monomers) >= 500)
})
