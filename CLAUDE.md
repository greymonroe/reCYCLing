# reCYCLing - Development Guide

## What this package does
R package for simulating evolution of tandem repeat arrays (satellite DNA monomers) and analyzing them. Four evolutionary forces per generation: local duplication, distal duplication, chunk deletion, per-base mutation. Includes ABC (Approximate Bayesian Computation) for inferring parameters from real data.

## Architecture

### Simulation engine
- **C++ backend** (`src/sim_engine.cpp`): ~100x faster than R, uses deferred mutation model (Jukes-Cantor). Only supports substitutions.
- **R backend** (`R/run_sim.R`): Slower but supports indels. `backend="auto"` selects C++ when possible.
- `run_sim_ps()` is the main entry point. Returns list: `monomers`, `ps` (identical pairs), `L_vec`, `H_vec`, `N_vec`.

### Shiny apps
- **Explorer** (`inst/shiny/explorer/app.R`): Interactive parameter tweaking + visualization. Works fine.
- **ABC** (`inst/shiny/abc/app.R`): SMC-ABC inference dashboard. **Currently problematic** - see below.

### CLI ABC script
- `inst/scripts/smc_abc.R`: Batch ABC with `parallel::mclapply` for forked parallelism.

## Current work: ABC Shiny app

### Known issue: ABC crashes the machine
The ABC Shiny app can exhaust system resources. Root causes:
1. **Slow simulations**: Many parameter combos in ABC's prior space produce arrays that grow too slowly to reach `max_units` within `max_t` generations, causing individual sims to run for a long time.
2. **No parallelization in Shiny** (it's serial via `invalidateLater(0)` reactive loop), but the tight loop + expensive sims = 100% CPU for extended periods.
3. **Memory pressure**: Large arrays (up to `hard_cap=40000`) with `pairs_identical()` can generate huge pair tables.

### Fixes applied
- Per-particle wall-clock timeout via `setTimeLimit()` (default 10s, configurable in UI). No forking.
- Reduced defaults: `max_units` 10000->5000, `max_t` 2000->500, particles 100->50.
- `hard_cap` set to `min(max_units * 2, 20000)` to prevent memory blowup.
- `invalidateLater(0)` changed to `invalidateLater(50)` so the Stop button and UI actually work.
- Status display shows per-particle timing, skip count, and generation elapsed time.
- Do NOT use parallelization (forking/mclapply) in the Shiny context.

## Key files
| File | Purpose |
|------|---------|
| `R/run_sim.R` | Main simulation function + helpers |
| `R/launch_app.R` | `launch_explorer()` and `launch_abc()` |
| `src/sim_engine.cpp` | C++ simulation engine |
| `inst/shiny/abc/app.R` | ABC Shiny dashboard |
| `inst/shiny/explorer/app.R` | Parameter explorer Shiny app |
| `inst/scripts/smc_abc.R` | CLI batch ABC script |
| `R/entropy.R` | Shannon entropy functions |
| `R/pairs_identical.R` | Identical monomer pair detection |
| `R/counts_analysis.R` | Base counts, mutation load, consensus |
| `R/plot_*.R` | Individual plotting functions |

## Build & test
```bash
# Build and install locally
R CMD build . && R CMD INSTALL reCYCLing_*.tar.gz
# Or from R:
devtools::install()
devtools::document()
```

## Conventions
- Uses `data.table` throughout (not tibble/dplyr)
- Plot functions use `ggplot2` + `cowplot` for multi-panel
- Chunk-size distributions specified as lists: `list(type="gamma", shape=2, scale=15)`
- Parameters sampled in log10 space for rates, natural space for shape/scale
- C++ code uses `Rcpp` (no OpenMP, single-threaded)
- `run_sim_ps(compute_pairs = FALSE)` skips `pairs_identical()` for large arrays (ABC uses sampled distances instead)
- Plot functions auto-compute pairs on the fly if `ps` is empty
