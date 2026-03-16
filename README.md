# reCYCLing

![reCYCLing logo](extra/logo.gif)

**Repeat Element Cyclical Evolutionary Simulations in Genomes**

Long-read sequencing is revealing repeat-rich genome features—centromeres, satellite arrays, and other tandem repeats—that were previously difficult to assemble and interpret. Understanding what we see in these regions often comes down to building intuition for how tandem arrays evolve under duplication, deletion, and mutation.

`reCYCLing` is an R toolkit for **simulating** the evolution of tandem repeat arrays (satellite-like monomers derived from a common ancestor), **analyzing** the resulting arrays using diagnostics that mirror how we interrogate real repeat architectures, and **inferring evolutionary parameters** from observed data via Approximate Bayesian Computation (ABC).

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("greymonroe/reCYCLing")
```

---

## Quick start

Simulate a single replicate, then generate core diagnostic plots:

```r
library(reCYCLing)

sim <- run_sim_ps(
  init_l = 178,             # monomer length (bp)
  init_k0 = 10,             # initial copy number
  max_units = 10000,        # stop when array reaches this many monomers
  max_t = 1e6,              # max generations (time cap)
  mu_total = 5e-5,          # per-base mutation rate
  p_local_dup = 4e-4,       # local duplication rate per unit per generation
  p_distal_dup = 1e-5,      # distal duplication rate per unit per generation
  p_del_chunk = 1e-5,       # chunk deletion rate per unit per generation
  verbose = FALSE
)

# Multi-panel summary
plot_ps_summary(sim)

# Individual plots
plot_repeat_fingerprint(sim)
plot_pair_distances(sim)
plot_mutation_load(sim)
```

---

## What `run_sim_ps()` returns

A named list with five elements:

| Element | Description |
|---------|-------------|
| `monomers` | `data.table` of the final array (columns: `num`, `bponly`, `pos`, `chrom`, `hap`, `sample`, `dir`) |
| `ps` | All pairs of identical monomers (`pairs_identical()` output). Empty if `compute_pairs = FALSE`. |
| `L_vec` | Array length trajectory over generations |
| `H_vec` | Shannon entropy trajectory (final value only with C++ backend) |
| `N_vec` | Unique monomer count trajectory (final value only with C++ backend) |

For large arrays (>10K units), use `compute_pairs = FALSE` to skip the O(n^2) pair enumeration. Plot functions will compute pairs on the fly when needed.

---

## Simulation engine

Each generation, the simulator applies four evolutionary forces:

| Process | Description |
|---------|-------------|
| **Local duplication** | Tandem copy-paste of a nearby chunk |
| **Distal duplication** | Copy-paste to a distant position (optional inversion) |
| **Chunk deletion** | Removal of a consecutive block of monomers |
| **Per-base mutation** | Substitutions (default), optional insertions/deletions |

Two backends are available:

- **C++ backend** (default): ~100x faster. Uses deferred mutations—mutations accumulate as a counter per monomer and are only materialized on duplication using the Jukes-Cantor model. Substitutions only.
- **R backend**: Slower but supports indels (`p_ins`, `p_del_bp`). Selected automatically when indels are requested.

Key parameters:

- **Rates (per unit per generation):** `p_local_dup`, `p_distal_dup`, `p_del_chunk`
- **Mutation (per base per generation):** `mu_total`
- **Chunk-size distributions:** `local_dist`, `distal_dist`, `del_dist`

Chunk-size distributions are specified as lists:

```r
local_dist  <- list(type = "gamma", shape = 2, scale = 15)
distal_dist <- list(type = "gamma", shape = 2, scale = 500)
del_dist    <- list(type = "geom",  prob = 0.1)
```

Supported types: `"fixed"`, `"poisson"`, `"normal"`, `"geom"`, `"unif"`, `"gamma"`.

---

## Diagnostic plots

| Function | What it shows |
|----------|---------------|
| `plot_repeat_fingerprint()` | Dot-plot of identical pairs (self-similarity / duplication structure) |
| `plot_pair_distances()` | Distance distribution of identical pairs (local vs distal signature) |
| `plot_copy_numbers()` | Copies per unique monomer sequence |
| `plot_mutation_load()` | Mismatches per monomer relative to consensus |
| `plot_consensus_support()` | Per-position consensus support across the array |
| `plot_ps_summary()` | Combined multi-panel overview |

---

## ABC parameter inference

`reCYCLing` includes an SMC-ABC (Sequential Monte Carlo Approximate Bayesian Computation) framework for inferring simulation parameters that reproduce observed summary statistics from real tandem repeat arrays.

### What is SMC-ABC?

ABC is used when you can simulate from a model but can't write down the likelihood function. SMC-ABC maintains a **population of parameter vectors** (particles) that evolves over generations:

1. **Generation 0**: Sample particles broadly from the prior
2. **Score**: Run a simulation for each particle, compute summary statistics, measure distance to the target
3. **Select**: Keep the top fraction of particles
4. **Perturb**: Jitter the survivors to propose new particles for the next generation
5. **Repeat**: Distance threshold tightens automatically as the population improves

The distance function is auto-scaled: after Generation 0, the SD of each summary statistic is computed from the pilot batch and used to normalize all distances. This ensures each stat contributes proportionally without manual weight tuning.

### Summary statistics

The ABC fits four summary statistics to the target:

- **Distance mode 1 & 2**: The two peaks of the bimodal identical-pair distance distribution (fit via 2-component Gaussian mixture in log10 space)
- **Weight of mode 1**: Relative contribution of the near-distance peak
- **Mean mutation load**: Average mismatches per monomer relative to the consensus

### Interactive app

```r
launch_abc()
```

Opens a Shiny dashboard for real-time ABC inference with live convergence plots, parameter posteriors, failure rate tracking, and a Best Fit Explorer for visualizing top-ranked simulations.

### Command-line batch script

```bash
Rscript inst/scripts/smc_abc.R \
  --target_mode1 100 --target_mode2 2000 --target_weight1 0.3 \
  --target_load 10 --particles 200 --max_generations 20 \
  --max_units 20000 --outdir results/run1
```

The CLI script runs serially (single core) with the same algorithm. Use `--help` for all options.

### Performance notes

- The C++ simulation backend with deferred mutations handles large arrays efficiently
- For ABC, `compute_pairs = FALSE` is used automatically to skip the O(n^2) pair table
- Pairwise distances are estimated via sampling (caps at 5K pairs per sequence group)
- The 2-component Gaussian mixture fit subsamples to 10K distances
- Failed particles are retried to fill all slots each generation

---

## Interactive explorer

```r
launch_explorer()
```

Opens a Shiny app for interactive parameter exploration — tweak sliders and immediately see the resulting array diagnostics.

---

## Vignette

For a worked walkthrough:

```r
vignette("reCYCLing-intro")
```

---

## Notes

This is a research codebase developed mainly for internal use in the Monroe Lab; it is evolving and documentation may change. Please feel free to try — contribute to Issues if you would like.
