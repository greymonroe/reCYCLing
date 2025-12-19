![](extra/logo.gif)

# reCYCLing  
## Repeat Element Cyclical Evolutionary Simulations in Genomes

**reCYCLing** is a small R toolkit for simulating the evolution of tandem repeat arrays (e.g., satellite DNA–like monomers) under:
- **Local (tandem) duplication** of chunks
- **Distal duplication** (copy–paste elsewhere to the right; optional inversion of the duplicated chunk)
- **Chunk deletion**
- **Per-base mutation** (substitutions by default; optional indels)

It was built for internal use in the Monroe Lab and is currently a **beta / research** codebase: it works for our use-cases, but is still evolving and documentation is lightweight.

For worked examples, see: **`extra/vignettecode.R`**.

---

## Installation (development / internal)

This repository isn’t currently a polished package. Typical usage is just:

```r
# from the repo root
source("R/functions.R")
```

### Suggested R packages
Some functions assume these are installed/available:

- `data.table`
- `ggplot2`
- `patchwork`
- `cowplot`
- `scales`
- `pbapply` (used by some helpers)

---

## Quick start

Run a small simulation from a **given** ancestral monomer sequence, then plot a compact summary:

```r
source("R/functions.R")

# define a starting monomer (single repeat unit)
startseq <- paste0(sample(c("A","C","G","T"), size = 178, replace = TRUE), collapse = "")

# simulate n replicates
ps_results <- run_sim_ps(
  init_sequence_type = "given",
  ancestor_seq       = startseq,
  init_l             = 178,
  init_k0            = 10,
  n                  = 5,
  max_t              = 1e6,
  mu_total           = 3e-5,
  p_del_chunk        = 0
)

# summarize replicate i
plot_ps_summary(
  ps_results, i = 1,
  title = "mu = 3e-5",
  rel_widths_bottom = c(1.2, 1, 1, 1),
  rel_heights = c(3, 1)
)
```

You can explore how different parameters affect the outcome by running multiple simulations (e.g., changing `mu_total`, turning off distal duplication, etc.)—see the vignette for patterns we’ve found useful.

---

## Core user-facing functions

### `run_sim_ps()`: run repeat-array evolution simulations

`run_sim_ps()` is the main entry point. It runs `n` independent replicates, each evolving a tandem array of monomer “units” through generations until an array-size or time cap is reached.

**Returns** a list with:
- `monomers_list`: list of per-replicate tables describing the final array (one row per unit)
- `ps_list`: list of per-replicate “pairwise identical” tables (all pairs of identical monomers)
- `L_vec_list`: list of length trajectories across generations (array copy number over time)
- `H_vec_list`: list of Shannon entropy trajectories across generations
- `N_vec_list`: list of unique-monomer-count trajectories across generations

**Key arguments (most commonly used):**

#### Replicates, stopping rules, and initial state
- `n`: number of simulation replicates
- `max_t`: maximum generations per replicate
- `max_units`: stop if the array reaches this many units
- `hard_cap`: safety stop if the array exceeds this many units
- `init_l`: length (bp) of one monomer unit
- `init_k0`: initial number of monomers in the array
- `init_sequence_type`: `"random"` (generate an ancestral monomer) or `"given"`
- `ancestor_seq`: required if `init_sequence_type="given"`; must have length `init_l`

#### Duplication and deletion rates (per unit per generation)
- `p_local_dup`: probability a unit triggers a **local** (tandem) duplication event
- `p_distal_dup`: probability a unit triggers a **distal** duplication event
- `p_invert_distal`: probability a distal-duplicated chunk is reversed in unit order
- `p_del_chunk`: probability a unit triggers a **deletion** event

#### Chunk-size distributions
Chunk sizes are drawn from *distribution spec lists*:
- `local_dist`: chunk-size distribution for local duplication
- `distal_dist`: chunk-size distribution for distal duplication
- `del_dist`: chunk-size distribution for deletion

Each distribution spec is a list with a `type` and parameters. Supported `type` values include:
- `"fixed"`: `list(type="fixed", value=10)`
- `"poisson"`: `list(type="poisson", lambda=20)`
- `"normal"`: `list(type="normal", mean=20, sd=5)`
- `"geom"`: `list(type="geom", prob=0.1)`
- `"unif"`: `list(type="unif", min=1, max=50)`
- `"gamma"`: `list(type="gamma", shape=2, scale=15)` (or `rate=`)

Chunk sizes are always coerced to an integer `>= 1` and clipped to what’s feasible at the sampled start position.

#### Mutation model (per base per generation)
- `mu_total`: per-base mutation probability (applied independently per position)
- `p_sub`, `p_ins`, `p_del_bp`: relative probabilities of substitution / insertion / deletion (default is substitutions only)

#### Other
- `verbose`: print generation progress

---

### `plot_ps_summary()`: one-figure overview for a replicate

`plot_ps_summary(ps_results, i=...)` produces a compact multi-panel diagnostic figure for one replicate:

- **Repeat fingerprint**: dot-plot of all identical monomer pairs (a “self-similarity” view)
- **Distance between identical pairs**: histogram of `|num1 - num2|` (log-scaled)
- **Copy number distribution**: histogram of copies per unique monomer sequence (log-scaled)
- **Mutation load distribution**: mismatches per monomer vs a consensus sequence
- **Consensus support distribution**: per-position consensus proportion across the array

Useful arguments:
- `title`: displayed above the fingerprint panel
- `rel_widths_bottom`: adjust widths of the small panels
- `rel_heights`: adjust fingerprint-vs-bottom-row height ratio

---

## Additional helpers you may use

These functions are handy if you want to dig deeper than `plot_ps_summary()`:

- `plot.repeat.fingerprint(ps, kb=..., rotate=FALSE, zoom=...)`  
  Fingerprint (dot-plot) view of identical monomer pairs; supports zoom windows.

- `pairs_identical(repeats_dt)`  
  Compute the table of all identical monomer pairs. (This is what `run_sim_ps()` stores in `ps_list`.)

- `counts_long_nogap(rawseqs)`  
  Per-position symbol counts/proportions across aligned sequences; returns a long table and consensus.

- `plot_circular_seqcounts(seqcounts)`  
  Circular visualization of per-position base composition (alpha = proportion).

- `get_sim_entropies(ps_results)`  
  Quick summary table (per replicate) of Shannon entropy, total units, and unique monomer count.

---

## Tips for interpreting output

- **High duplication + low mutation** tends to produce many identical monomers and dense fingerprints.
- **Higher mutation** increases unique monomer diversity, often raising entropy and broadening the mutation-load distribution.
- **Distal duplication** can create “off-diagonal” structure in fingerprints (depending on chunk sizes and insertion positions).
- `L_vec_list`, `H_vec_list`, and `N_vec_list` are useful for plotting trajectories (growth, diversification, etc.) across generations.

---

## Vignette / examples

See **`extra/vignettecode.R`** for:
- Parameter sweeps (e.g., varying `mu_total`)
- Turning mechanisms on/off (e.g., set `p_distal_dup = 0`)
- Example figures and interpretation patterns

---

## Status, feedback, and contributions

This is a **beta** research tool. If you try it and have ideas (bug reports, feature requests, API suggestions), please:
- open an Issue on the repository **Issues** page, or
- email: **gmonroe@ucdavis.edu**

We’re happy to accept feedback andcontributions—especially if it helps make the simulator easier for collaborators to use.

