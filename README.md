# reCYCLing

![reCYCLing logo](extra/logo.gif)

**Repeat Element Cyclical Evolutionary Simulations in Genomes**

Long-read sequencing is revealing repeat-rich genome features—centromeres, satellite arrays, and other tandem repeats—that were previously difficult to assemble and interpret. Understanding what we see in these regions often comes down to building intuition for how tandem arrays evolve under duplication, deletion, and mutation.

`reCYCLing` is an R toolkit for **simulating** the evolution of tandem repeat arrays (satellite-like monomers derived from a common ancestor) and for **analyzing** the resulting arrays using diagnostics that mirror how we interrogate real repeat architectures (e.g., arrays identified with TRASH-like workflows). A core analysis theme is leveraging **pairs of identical monomers** to study spatial structure and duplication signatures.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("greymonroe/reCYCLing")
```

---

## Quick start

Simulate a single replicate from a fixed ancestral monomer sequence (so results are reproducible), then generate core diagnostic plots.

```r
library(reCYCLing)

# Define a fixed ancestral monomer so results are reproducible
ancestor <- paste(sample(c("A","C","G","T"), 178, replace = TRUE), collapse = "")

sim <- run_sim_ps(
  n = 1,                    # number of replicates
  init_l = 178,             # monomer length (bp)
  init_k0 = 10,             # initial copy number
  init_sequence_type = "given",
  ancestor_seq = ancestor,
  max_units = 10000,        # stop when array reaches this many monomers
  max_t = 1e6,              # max generations (time cap)
  mu_total = 5e-5,          # per-base mutation rate
  p_local_dup = 4e-4,       # local duplication probability per unit per generation
  p_distal_dup = 1e-5,      # distal duplication probability per unit per generation
  p_del_chunk = 1e-5,       # chunk deletion probability per unit per generation
  verbose = FALSE
)

plot_repeat_fingerprint(sim, i = 1, ptsize = 0.5)
plot_pair_distances(sim, i = 1)
plot_mutation_load(sim, i = 1)
```

---

## What `run_sim_ps()` returns

`run_sim_ps()` returns a **named list with five elements**, each a list of length `n` (one entry per replicate):

- `monomers_list`: final array as a `data.table` (one row per monomer; key columns include `num` position and `bponly` sequence)
- `ps_list`: all **pairs of identical monomers** (`pairs_identical()` output), with positions `num1` and `num2`
- `L_vec_list`: array length trajectory over generations
- `H_vec_list`: Shannon entropy trajectory over generations
- `N_vec_list`: unique-monomer-count trajectory over generations

Example:

```r
mono <- sim$monomers_list[[1]]
ps   <- sim$ps_list[[1]]

head(mono)
head(ps[, c("bponly","num1","num2")])

# Quick summaries
length(unique(mono$bponly))
shannon_entropy(mono$bponly)
```

---

## Core ideas and controls

Each generation, the simulator applies four forces to an array of monomer units:

| Process | Controls |
|---|---|
| **Local duplication** | tandem copy–paste of a nearby chunk |
| **Distal duplication** | copy–paste to a more distant position (optional inversion) |
| **Chunk deletion** | removal of a consecutive block of monomers |
| **Per-base mutation** | substitutions by default (optional indels) |

Key knobs you’ll typically adjust:

- **Rates (per unit per generation):** `p_local_dup`, `p_distal_dup`, `p_del_chunk`
- **Distal inversion:** `p_invert_distal`
- **Mutation (per base per generation):** `mu_total` and optionally `p_sub`, `p_ins`, `p_del_bp`
- **Chunk-size distributions:** `local_dist`, `distal_dist`, `del_dist`

Chunk-size distributions are specified as lists, e.g.:

```r
local_dist  <- list(type="gamma", shape=2, scale=15)
distal_dist <- list(type="gamma", shape=2, scale=500)
del_dist    <- list(type="geom",  prob=0.1)
```

Supported `type` values include: `"fixed"`, `"poisson"`, `"normal"`, `"geom"`, `"unif"`, `"gamma"`.

---

## Diagnostic plots

Each diagnostic view is available as a standalone plotting function:

- `plot_repeat_fingerprint()` — dot-plot of identical pairs (self-similarity / duplication structure)
- `plot_pair_distances()` — distance distribution of identical pairs (local vs distal signature)
- `plot_copy_numbers()` — copies per unique monomer sequence
- `plot_mutation_load()` — mismatches per monomer relative to consensus
- `plot_consensus_support()` — per-position consensus support across the array
- `plot_ps_summary()` — combined multi-panel overview

Tip: `plot_mutation_load()` and `plot_consensus_support()` both use per-position counts; compute once and reuse:

```r
counts <- counts_long_nogap(sim$monomers_list[[1]]$bponly)
plot_mutation_load(sim, i = 1, counts = counts)
plot_consensus_support(sim, i = 1, counts = counts)
```

---

## Vignette

For a worked walkthrough and more examples:

```r
vignette("reCYCLing-intro")
```

---

## Notes

This is a research codebase developed mainly for internal use in the Monroe Lab; it is evolving and documentation may change. Please feel free to try - contribute to Issues if you would like
