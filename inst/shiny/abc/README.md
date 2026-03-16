# reCYCLing ABC Parameter Inference — User Guide

## Purpose

This Shiny app uses **Approximate Bayesian Computation (ABC)** to infer the evolutionary parameters that shaped an observed tandem repeat array. Given empirical measurements from a real repeat array (distance distribution, mutation load, array size), the ABC searches for combinations of duplication, deletion, and mutation rates that reproduce those observations in simulation.

The scientific goal: demonstrate that the bimodal distance distribution observed between identical monomers in real tandem repeat arrays (a near-distance peak from local/tandem duplication, and a far-distance peak from distal/crossover duplication) is consistent with a model involving two distinct duplication mechanisms operating at different spatial scales.

---

## Two kinds of "generations"

This is a common source of confusion. There are two completely different uses of the word "generation" in this system:

### Simulation generations (biological time)
Each simulation evolves a tandem repeat array over many **simulation generations**. In each generation, four evolutionary forces act:

1. **Local (tandem) duplication** — a nearby chunk of monomers is copied and inserted adjacent to the source. Rate is **per-unit per-generation**: each monomer independently has a probability of triggering a tandem duplication event. More monomers = more events.

2. **Distal duplication** — a chunk is copied and inserted at a distant position in the array. Models unequal crossover or similar chromosome-level events. Rate is **per-array per-generation**: at most one event per generation, regardless of array size.

3. **Chunk deletion** — a consecutive block of monomers is removed. Rate is **per-unit per-generation**.

4. **Per-base mutation** — substitutions accumulate in each monomer. In the C++ backend, mutations are deferred (tracked as a counter) and only materialized when a monomer is duplicated, using the Jukes-Cantor model. This is mathematically equivalent but much faster.

The number of simulation generations is determined by the **molecular clock** (see below). It is NOT a free parameter — it is derived from the ancestor age and generation time.

### ABC iterations (algorithm steps)
The ABC algorithm itself runs through **iterations** (sometimes called "ABC generations" in the interface). Each iteration:

1. **Propose** N parameter sets ("particles")
2. **Simulate** — run a full evolutionary simulation for each particle
3. **Score** — compute summary statistics and measure distance to target
4. **Select** — keep the best fraction of particles
5. **Perturb** — add noise to survivors to create the next round of proposals
6. **Repeat**

Typically 10-20 ABC iterations are sufficient for convergence.

---

## What is a "particle"?

A **particle** is one candidate set of simulation parameters. For example, a particle might be:

```
p_local_dup = 0.0003, p_distal_dup = 0.05, p_del_chunk = 0.00001
```

Each ABC iteration proposes many particles (e.g., 500), runs a simulation for each one, and scores how well the result matches the target data. The best-scoring particles survive to the next iteration.

---

## The molecular clock framework

Rather than treating evolutionary time as a free parameter (which creates an identifiability problem — fast rates over short time looks the same as slow rates over long time), the simulation is anchored to real evolutionary time using a molecular clock:

### Known quantities (from empirical data)
- **Ancestor age**: How old is the common ancestor of these repeats? (from phylogenetic dating)
- **Generation time**: Years per generation for the organism
- **Monomer length**: Base pairs per repeat unit (from sequence data)
- **Observed mutation load**: Mismatches per monomer relative to consensus

### Derived quantities
- **Real generations** = ancestor_age / generation_time
- **Real mutation rate** = observed_load / (monomer_length × real_generations)
- **Compressed mutation rate** = real_rate × compression_factor
- **Simulated generations** = real_generations / compression_factor

### Time compression
We cannot simulate millions of generations directly. Instead, all per-generation rates are multiplied by a compression factor (e.g., 1000×) and the number of generations is divided by the same factor. The product `rate × time` is preserved, so the result is mathematically equivalent.

Example for pistachio knob repeats:
- Ancestor age: 10 million years
- Generation time: 10 years
- Real generations: 1,000,000
- Observed load: ~10 mismatches per 178bp monomer
- Derived mu: 10 / (178 × 1,000,000) ≈ 5.6 × 10⁻⁸ per base per gen
- With 1000× compression: mu_eff = 5.6 × 10⁻⁵, simulated gens = 1,000

---

## Parameters being fitted

### Fitted by default (3 rate parameters)

| Parameter | Type | Prior range (log10) | Description |
|-----------|------|-------------------|-------------|
| `p_local_dup` | Per-unit per-gen | [-4.2, -3.0] | Tandem duplication rate. At 20K units, this produces ~1-20 dup events/gen |
| `p_distal_dup` | Per-array per-gen | [-2, -0.3] | Crossover/distal duplication rate. 1-50% chance per generation |
| `p_del_chunk` | Per-unit per-gen | [-6, -4] | Chunk deletion rate |

### Fixed by default (6 chunk-size parameters)

These control the gamma distribution of duplication/deletion chunk sizes. They are fixed because they have weak effects on the summary statistics and are difficult to identify from the data. They can be enabled for fitting via checkboxes.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `local_shape` | 2 | Gamma shape for local dup chunks |
| `local_scale` | 15 | Gamma scale for local dup chunks (mean = shape × scale = 30 units) |
| `distal_shape` | 2 | Gamma shape for distal dup chunks |
| `distal_scale` | 2000 | Gamma scale for distal dup chunks (mean = 4000 units) |
| `del_shape` | 2 | Gamma shape for deletion chunks |
| `del_scale` | 15 | Gamma scale for deletion chunks (mean = 30 units) |

The gamma distribution with shape=2 produces a peaked, right-skewed distribution. Shape < 1 is exponential (mostly small events). Shape > 5 approaches a bell curve.

### Fixed always (from molecular clock)

| Parameter | Source | Description |
|-----------|--------|-------------|
| `mu_total` | Derived from clock | Per-base per-generation mutation rate (compressed) |
| `max_t` | Derived from clock | Number of simulation generations |
| `init_l` | Sequence data | Monomer length in bp |
| `init_k0` | Fixed at 10 | Initial copy number |

---

## Summary statistics (what the ABC fits to)

The ABC compares each simulation's output to target values using 5 summary statistics:

| Statistic | How it's computed | What it constrains |
|-----------|-------------------|-------------------|
| **Median near-distance** | Median of identical-pair distances ≤ threshold | Local duplication chunk size and rate |
| **Median far-distance** | Median of identical-pair distances > threshold | Distal duplication chunk size |
| **Near fraction** | Fraction of pairs below threshold | Balance of local vs distal duplication |
| **Mean mutation load** | Mean mismatches to consensus per monomer | Evolutionary time (molecular clock) |
| **Array size** | Number of monomers in final array | Duplication-deletion equilibrium |

The **near/far threshold** (default 500 units) splits the distance distribution into local and distal components. Set this based on the valley in your empirical distance histogram.

---

## Distance function

The ABC scores each particle by measuring how far its summary statistics are from the target. The distance function uses:

1. **Auto-scaling**: After the first ABC iteration (Gen 0), the standard deviation of each summary statistic across the pilot batch is computed. All subsequent distances are normalized by these SDs, so each statistic contributes proportionally regardless of its natural scale.

2. **Relative error for mutation load**: Load difference is divided by the target load, making it a dimensionless ratio comparable to the log-space mode differences.

3. **Log-space for array size**: Array size enters the distance as `log10(sim_size) - log10(target_size)`, so a 2× error counts the same whether the array is 10K or 40K.

4. **Blended distance**: `0.5 × euclidean + 0.5 × max_deviation`. This combines overall fit (euclidean) with a penalty for the worst-fitting stat (max deviation). This prevents the ABC from optimizing one statistic at the expense of others.

### Mutation load gate

Before computing full summary statistics, a quick subsample of 200 monomers is used to estimate mutation load. If the estimate is outside [target/3, target×3], the particle is immediately rejected. This is a hard biological constraint — if the load is wildly wrong, the parameter combination is biologically impossible given the molecular clock.

---

## The SMC-ABC algorithm in detail

### Generation 0 (prior sampling)
- Sample N particles uniformly from the prior ranges
- Run a simulation for each particle
- Compute summary statistics and distances
- Many particles will fail (wrong load, array too small/large, no identical pairs)
- Compute stat scales from the valid particles (pilot batch normalization)
- Recompute all distances with proper scaling

### Subsequent generations (perturbation)
- **Select**: Keep the top fraction (e.g., 30%) of valid particles by distance
- **Compute adaptive perturbation**: For each parameter, perturbation SD = min(2 × empirical SD of survivors, 15% of prior range). This adapts to the actual spread of the posterior — wide when uncertain, narrow when converging.
- **Propose**: For each new particle, randomly pick a survivor and add Gaussian noise with the adaptive SD
- **Simulate and score**: Run simulation, compute stats, compute distance
- **No retry**: Failed particles count toward the generation total. If 500 particles are proposed and 200 fail, the generation proceeds with 300 valid particles.

### Convergence detection
The algorithm stops when either:
- The median distance hasn't improved by >1% for 3 consecutive iterations
- The maximum number of iterations is reached
- Fewer than 3 valid particles remain

---

## Common issues and challenges

### High failure rates
**Symptom**: 50%+ of particles fail, especially in later iterations.
**Causes**:
- Array explodes past `hard_cap` (= 1.5× target size). The per-unit local duplication rate creates positive feedback: more units → more dup events → even more units.
- Array goes extinct (deletion overwhelms duplication)
- Mutation load outside the [target/3, target×3] gate
- Too few identical pairs to compute statistics

**Solutions**:
- Use the tighter prior for `p_local_dup` ([-4.2, -3.0])
- Fix chunk-size parameters (reduce to 3 fitted params)
- Increase the number of particles per iteration to compensate

### Non-convergence / divergence
**Symptom**: Distance bounces around instead of decreasing.
**Causes**:
- Too many parameters (9D space is hard to search with 4-5 summary stats)
- Perturbation too wide (the adaptive kernel hits the 15% cap every iteration)
- Parameter interactions: changing one rate affects multiple summary stats simultaneously

**Solutions**:
- Fix chunk-size parameters (3-parameter ABC converges much better)
- Increase particles per iteration (500+)
- Lower retention fraction (0.2-0.3) for stronger selection

### Explosive growth
**Symptom**: Simulations create millions of monomers and take seconds each.
**Cause**: `p_local_dup` is too high. At 20K units, `p_local = 10^-2` means 200 tandem dup events per generation — the array doubles in one step.
**Solution**: The prior upper bound of `10^-3` limits this. The `hard_cap` (1.5× target) kills explosive sims quickly in the C++ engine.

### Identifiability
**Issue**: Multiple parameter combinations can produce similar summary statistics. For example, high local dup + high deletion ≈ low local dup + low deletion (both produce similar-sized arrays). The molecular clock breaks the most severe degeneracy (rate-time tradeoff), but some residual non-identifiability remains for deletion vs duplication rates.

---

## File structure

```
inst/shiny/abc/
├── app.R          # The Shiny application (UI + server)
├── README.md      # This file

inst/scripts/
├── smc_abc.R      # CLI batch version (same algorithm, no UI)

R/
├── run_sim.R      # Core simulation function (run_sim_ps)
├── pairs_identical.R  # Identical pair detection

src/
├── sim_engine.cpp # C++ simulation backend (deferred mutations)
```

---

## Running from the command line

The CLI script `inst/scripts/smc_abc.R` implements the same ABC algorithm without the Shiny interface. Useful for batch runs on HPC clusters.

```bash
Rscript inst/scripts/smc_abc.R \
  --target_med_near 20 --target_med_far 5000 \
  --target_near_frac 0.5 --target_load 10 \
  --near_far_threshold 500 \
  --particles 500 --max_generations 20 \
  --outdir results/run1
```

---

## References

- Beaumont, M.A. (2009). Adaptive approximate Bayesian computation. *Biometrika*, 96(4), 983-990. (Adaptive perturbation kernel)
- Sisson, S.A., Fan, Y., & Tanaka, M.M. (2007). Sequential Monte Carlo without likelihoods. *PNAS*, 104(6), 1760-1765. (SMC-ABC framework)
- Toni, T., et al. (2009). Approximate Bayesian computation scheme for parameter inference and model selection in dynamical systems. *J R Soc Interface*, 6(31), 187-202. (ABC tutorial)
