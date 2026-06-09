#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <cstdint>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Monomer representation: integer-encoded bases + deferred mutation counter
// ---------------------------------------------------------------------------
// Bases: 0=A, 1=C, 2=G, 3=T
// pending_gens: number of mutation rounds not yet applied
// dir: 1 = "+", -1 = "-"

struct Monomer {
    std::vector<int> seq;
    int pending_gens;
    int dir;

    Monomer() : pending_gens(0), dir(1) {}
    Monomer(const std::vector<int>& s, int d) : seq(s), pending_gens(0), dir(d) {}
};

// ---------------------------------------------------------------------------
// Chunk size sampling (mirrors R's .sample_chunk_size)
// ---------------------------------------------------------------------------

// Expected (mean) chunk size for a chunk-size distribution spec. Used by the
// genome steady-state governor to balance duplication mass against deletion
// mass (net-zero turnover). Falls back to 1.0 if the mean is undefined/<1.
double dist_mean(const List& dist) {
    std::string type = as<std::string>(dist["type"]);
    double m;
    if (type == "fixed") {
        m = as<double>(dist["value"]);
    } else if (type == "poisson") {
        m = as<double>(dist["lambda"]);
    } else if (type == "normal") {
        m = as<double>(dist["mean"]);
    } else if (type == "geom") {
        double p = as<double>(dist["prob"]);
        m = (p > 0.0) ? (1.0 - p) / p : 1.0;   // R::rgeom mean = (1-p)/p
    } else if (type == "unif") {
        m = 0.5 * (as<double>(dist["min"]) + as<double>(dist["max"]));
    } else if (type == "gamma") {
        double shape = as<double>(dist["shape"]);
        double scale_val = dist.containsElementNamed("scale")
            ? as<double>(dist["scale"]) : 1.0 / as<double>(dist["rate"]);
        m = shape * scale_val;
    } else {
        m = 1.0;
    }
    if (!(m >= 1.0)) m = 1.0;   // also catches NaN
    return m;
}

int sample_chunk_size(const List& dist, int max_k) {
    std::string type = as<std::string>(dist["type"]);
    int s;

    if (type == "fixed") {
        s = (int)as<double>(dist["value"]);
    } else if (type == "poisson") {
        s = (int)R::rpois(as<double>(dist["lambda"]));
    } else if (type == "normal") {
        s = (int)std::round(R::rnorm(as<double>(dist["mean"]),
                                      as<double>(dist["sd"])));
    } else if (type == "geom") {
        s = (int)R::rgeom(as<double>(dist["prob"]));
    } else if (type == "unif") {
        s = (int)std::round(R::runif(as<double>(dist["min"]),
                                      as<double>(dist["max"])));
    } else if (type == "gamma") {
        double shape = as<double>(dist["shape"]);
        double scale_val;
        if (dist.containsElementNamed("scale")) {
            scale_val = as<double>(dist["scale"]);
        } else {
            scale_val = 1.0 / as<double>(dist["rate"]);
        }
        s = (int)std::round(R::rgamma(shape, scale_val));
    } else {
        stop("Unknown chunk size distribution type: %s", type.c_str());
    }

    if (s < 1) s = 1;
    if (s > max_k) s = max_k;
    return s;
}

// ---------------------------------------------------------------------------
// Materialize pending mutations using Jukes-Cantor batch formula
// ---------------------------------------------------------------------------
// After g generations of JC substitution with per-base rate mu:
//   P(base unchanged) = 1/4 + 3/4 * (1 - 4*mu/3)^g
//   P(base -> specific other) = (1 - (1 - 4*mu/3)^g) / 4
//
// This is exactly equivalent in distribution to applying mutations
// one generation at a time.

void materialize(Monomer& m, double mu_total) {
    int g = m.pending_gens;
    if (g == 0 || mu_total <= 0.0) return;

    double r = std::pow(1.0 - 4.0 * mu_total / 3.0, (double)g);
    double p_same = 0.25 + 0.75 * r;

    int L = (int)m.seq.size();
    for (int i = 0; i < L; i++) {
        if (R::runif(0.0, 1.0) >= p_same) {
            // Change to a uniformly random different base
            int offset = 1 + (int)(R::runif(0.0, 1.0) * 3);
            if (offset > 3) offset = 3; // guard against runif == 1.0
            m.seq[i] = (m.seq[i] + offset) % 4;
        }
    }
    m.pending_gens = 0;
}

// ---------------------------------------------------------------------------
// AUTOCATALYTIC LOCAL REDUNDANCY VECTOR (rich-get-richer local duplication)
// ---------------------------------------------------------------------------
// For each monomer position i, r[i] in [0,1] = the LOCAL REDUNDANCY around i:
// the fraction of monomers within +/- `window` positions of i whose EXACT
// (current, possibly-not-yet-materialized) sequence occurs >= 2 times WITHIN
// that same window (i.e. locally repeated). The effective per-monomer local-dup
// trigger probability is then p_local_dup * (1 + alpha * r[i]).
//
// Cost: O(N) per call. We slide a fixed-half-width window [i-window, i+window]
// over the array, maintaining a hash multiset of sequence ids inside the window
// and a running count `n_red` = how many positions currently in the window have
// an id whose multiplicity is >= 2. As the window advances by one i, exactly one
// position may enter on the right and one may leave on the left; each add/remove
// is O(1) amortized (a hash lookup) and updates n_red by a constant. Then
// r[i] = n_red / (#positions currently in window). No RNG is consumed, so with
// alpha == 0 the caller's behavior is bit-identical to the original engine
// (the multiplier collapses to 1 and the same runif draws are compared to the
// same p_local_dup) -- this function is only invoked when alpha > 0.
//
// Sequence identity uses a 64-bit FNV-1a hash of the integer-encoded seq plus
// the length, to avoid building std::string keys for every monomer every recompute.

static inline uint64_t seq_hash(const std::vector<int>& s) {
    uint64_t h = 1469598103934665603ULL;            // FNV-1a offset basis
    for (size_t j = 0; j < s.size(); j++) {
        h ^= (uint64_t)(s[j] + 1);
        h *= 1099511628211ULL;                       // FNV-1a prime
    }
    h ^= (uint64_t)(s.size() + 1) * 1099511628211ULL; // fold in length
    return h;
}

void compute_local_redundancy(const std::vector<Monomer>& arr,
                              int window,
                              std::vector<double>& r_out) {
    int n = (int)arr.size();
    r_out.assign(n, 0.0);
    if (n == 0) return;
    if (window < 0) window = 0;

    // Precompute per-position sequence hash once.
    std::vector<uint64_t> ids(n);
    for (int i = 0; i < n; i++) ids[i] = seq_hash(arr[i].seq);

    std::unordered_map<uint64_t, int> cnt;
    cnt.reserve((size_t)std::min(n, 2 * window + 1) * 2 + 16);

    int n_red = 0;            // # positions in current window whose id multiplicity >= 2
    int lo = 0, hi = -1;      // current window is [lo, hi] inclusive (none yet)

    // helper lambdas adjust n_red as ids enter/leave the window
    auto add_pos = [&](int p) {
        int c = ++cnt[ids[p]];
        if (c == 2) n_red += 2;        // this one + the previously-singleton one
        else if (c > 2) n_red += 1;    // just this one becomes redundant
    };
    auto rem_pos = [&](int p) {
        int c = --cnt[ids[p]];
        if (c == 1) n_red -= 2;        // the remaining one is no longer redundant
        else if (c >= 2) n_red -= 1;   // only the removed one was counted
        // c == 0: it was a singleton (not counted), nothing to subtract
    };

    for (int i = 0; i < n; i++) {
        int want_lo = i - window; if (want_lo < 0) want_lo = 0;
        int want_hi = i + window; if (want_hi > n - 1) want_hi = n - 1;
        // grow/shrink the window to [want_lo, want_hi]
        while (hi < want_hi) { hi++; add_pos(hi); }
        while (lo < want_lo) { rem_pos(lo); lo++; }
        // (window only ever moves forward as i increases, so hi never shrinks and
        //  lo never moves left -- each position added/removed at most once total => O(N))
        int wsize = hi - lo + 1;
        r_out[i] = (wsize > 0) ? (double)n_red / (double)wsize : 0.0;
    }
}

// ---------------------------------------------------------------------------
// Per-array within-chromosome moves (shared by single-array and multi-chrom).
// Operates IN PLACE on one array. Mirrors exactly the move logic inside
// sim_core_cpp so both engines behave identically per chromosome.
// (sim_core_cpp keeps its own inline copy unchanged for safety; this helper is
// used only by the new multi-chromosome engine sim_genome_cpp.)
// ---------------------------------------------------------------------------

void apply_within_array_moves(std::vector<Monomer>& arr,
                              double p_local_dup,
                              const List& local_dist,
                              double p_distal_dup,
                              const List& distal_dist,
                              double p_invert_distal,
                              double p_del_chunk,
                              const List& del_dist,
                              double p_distal_del,
                              double p_conversion,
                              const List& conv_dist,
                              const List& conv_tract,
                              double mu_total,
                              double p_autocorr_alpha,
                              int autocorr_window,
                              bool autocorr_active) {
    int k = (int)arr.size();
    if (k == 0) return;

    // === LOCAL DUPLICATION ===
    // AUTOCATALYTIC: when autocorr_active (alpha>0 and this is a recompute gen),
    // the effective trigger prob at position i is p_local_dup*(1+alpha*r[i]),
    // where r[i] is the local redundancy. Otherwise the original uniform rule.
    // At alpha==0 the caller passes autocorr_active=false => bit-identical RNG.
    if (p_local_dup > 0.0) {
        k = (int)arr.size();
        std::vector<double> r_loc;
        bool use_ac = autocorr_active && p_autocorr_alpha > 0.0;
        if (use_ac) compute_local_redundancy(arr, autocorr_window, r_loc);
        std::vector<int> triggered;
        for (int i = 0; i < k; i++) {
            double p_eff = use_ac ? p_local_dup * (1.0 + p_autocorr_alpha * r_loc[i])
                                  : p_local_dup;
            if (p_eff > 1.0) p_eff = 1.0;
            if (R::runif(0.0, 1.0) < p_eff) triggered.push_back(i);
        }
        for (int ti = (int)triggered.size() - 1; ti >= 0; ti--) {
            int idx = triggered[ti];
            int ck = (int)arr.size();
            int max_chunk = ck - idx;
            if (max_chunk <= 0) continue;
            int chunk_size = sample_chunk_size(local_dist, max_chunk);
            int start = idx;
            int end = std::min(start + chunk_size - 1, ck - 1);
            for (int j = start; j <= end; j++) materialize(arr[j], mu_total);
            std::vector<Monomer> chunk;
            chunk.reserve(end - start + 1);
            for (int j = start; j <= end; j++)
                chunk.emplace_back(arr[j].seq, arr[j].dir);
            arr.insert(arr.begin() + end + 1,
                       std::make_move_iterator(chunk.begin()),
                       std::make_move_iterator(chunk.end()));
        }
    }

    // === DISTAL DUPLICATION (within array) ===
    if (p_distal_dup > 0.0 && R::runif(0.0, 1.0) < p_distal_dup) {
        k = (int)arr.size();
        if (k > 0) {
            int idx = (int)(R::runif(0.0, 1.0) * k);
            if (idx >= k) idx = k - 1;
            int max_chunk = k - idx;
            int chunk_size = sample_chunk_size(distal_dist, max_chunk);
            int start = idx;
            int end = std::min(start + chunk_size - 1, k - 1);
            for (int j = start; j <= end; j++) materialize(arr[j], mu_total);
            bool invert = R::runif(0.0, 1.0) < p_invert_distal;
            std::vector<Monomer> chunk;
            chunk.reserve(end - start + 1);
            if (invert) {
                for (int j = end; j >= start; j--)
                    chunk.emplace_back(arr[j].seq, -arr[j].dir);
            } else {
                for (int j = start; j <= end; j++)
                    chunk.emplace_back(arr[j].seq, arr[j].dir);
            }
            int n_positions = k + 1;
            int insert_at = (int)(R::runif(0.0, 1.0) * n_positions);
            if (insert_at >= n_positions) insert_at = n_positions - 1;
            arr.insert(arr.begin() + insert_at,
                       std::make_move_iterator(chunk.begin()),
                       std::make_move_iterator(chunk.end()));
        }
    }

    // === CHUNK DELETION ===
    if (p_del_chunk > 0.0) {
        k = (int)arr.size();
        std::vector<int> triggered;
        for (int i = 0; i < k; i++) {
            if (R::runif(0.0, 1.0) < p_del_chunk) triggered.push_back(i);
        }
        for (int ti = (int)triggered.size() - 1; ti >= 0; ti--) {
            int idx = triggered[ti];
            int ck = (int)arr.size();
            if (idx >= ck) continue;
            int max_chunk = ck - idx;
            if (max_chunk <= 0) continue;
            int chunk_size = sample_chunk_size(del_dist, max_chunk);
            int start = idx;
            int end = std::min(start + chunk_size - 1, ck - 1);
            arr.erase(arr.begin() + start, arr.begin() + end + 1);
        }
    }

    // === GENE CONVERSION (homogenization) ===
    if (p_conversion > 0.0) {
        k = (int)arr.size();
        if (k >= 2) {
            std::vector<int> recip;
            for (int i = 0; i < k; i++) {
                if (R::runif(0.0, 1.0) < p_conversion) recip.push_back(i);
            }
            for (size_t ri = 0; ri < recip.size(); ri++) {
                int i = recip[ri];
                int off = sample_chunk_size(conv_dist, k - 1);
                int sign = (R::runif(0.0, 1.0) < 0.5) ? -1 : 1;
                int donor = i + sign * off;
                if (donor < 0)      donor = i + off;
                if (donor > k - 1)  donor = i - off;
                if (donor < 0 || donor > k - 1 || donor == i) continue;
                int room = 1 + std::min(k - 1 - i, k - 1 - donor);
                int tract = sample_chunk_size(conv_tract, room);
                for (int j = 0; j < tract; j++) {
                    int rpos = i + j, dpos = donor + j;
                    if (rpos > k - 1 || dpos > k - 1) break;
                    materialize(arr[dpos], mu_total);
                    arr[rpos].seq = arr[dpos].seq;
                    arr[rpos].pending_gens = 0;
                }
            }
        }
    }

    // === TERMINAL (LONG-RANGE) DELETION ===
    if (p_distal_del > 0.0 && R::runif(0.0, 1.0) < p_distal_del) {
        k = (int)arr.size();
        if (k > 0) {
            int chunk_size = sample_chunk_size(distal_dist, k);
            if (chunk_size > k) chunk_size = k;
            arr.erase(arr.end() - chunk_size, arr.end());
        }
    }
}

// ---------------------------------------------------------------------------
// Core simulation engine
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List sim_core_cpp(IntegerVector ancestor_seq_r,
                  int init_k0,
                  int max_units,
                  int max_t,
                  int hard_cap,
                  double p_local_dup,
                  List local_dist,
                  double p_distal_dup,
                  List distal_dist,
                  double p_invert_distal,
                  double p_del_chunk,
                  List del_dist,
                  double p_distal_del,
                  double p_conversion,
                  List conv_dist,
                  List conv_tract,
                  double mu_total,
                  double p_autocorr_alpha,
                  int autocorr_window,
                  int autocorr_every,
                  bool verbose) {

    // Convert ancestor sequence: R 1-indexed -> C++ 0-indexed
    int L = ancestor_seq_r.size();
    std::vector<int> ancestor_seq(L);
    for (int i = 0; i < L; i++) {
        ancestor_seq[i] = ancestor_seq_r[i] - 1;
    }

    // Initialize array with init_k0 identical monomers
    std::vector<Monomer> arr;
    arr.reserve(std::min(max_units, hard_cap) + 1000);
    for (int i = 0; i < init_k0; i++) {
        arr.emplace_back(ancestor_seq, 1);
    }

    // Trajectory: array length per generation
    std::vector<int> l_vec;
    l_vec.reserve(100000);
    l_vec.push_back(init_k0);

    int t = 0;
    bool hit_hard_cap = false;

    while ((int)arr.size() < max_units && t < max_t) {
        t++;
        int k = (int)arr.size();

        if (k > hard_cap) {
            hit_hard_cap = true;
            break;
        }
        // If array went extinct, re-initialize and reset clock
        if (k == 0) {
            arr.clear();
            for (int i = 0; i < init_k0; i++) {
                arr.emplace_back(ancestor_seq, 1);
            }
            t = 0;  // reset clock — new knob starts fresh
            continue;
        }

        if (verbose && t % 500 == 0) {
            Rprintf("  Gen %d - units: %d\n", t, k);
        }
        if (t % 100 == 0) Rcpp::checkUserInterrupt();

        // === LOCAL DUPLICATION ===
        // AUTOCATALYTIC (rich-get-richer): with autocorr enabled, the per-monomer
        // trigger prob at position i is p_local_dup*(1 + alpha*r[i]) where r[i] is
        // local redundancy (see compute_local_redundancy). To keep cost bounded for
        // big arrays over many generations, r is recomputed only every
        // `autocorr_every` generations (cadence); between recomputes the original
        // uniform rule is used. At alpha==0 we never enter the autocatalytic branch,
        // so RNG consumption is bit-identical to the original engine.
        if (p_local_dup > 0.0) {
            k = (int)arr.size();
            bool use_ac = (p_autocorr_alpha > 0.0) &&
                          (autocorr_every <= 1 || (t % autocorr_every) == 0);
            std::vector<double> r_loc;
            if (use_ac) compute_local_redundancy(arr, autocorr_window, r_loc);
            std::vector<int> triggered;
            for (int i = 0; i < k; i++) {
                double p_eff = use_ac
                    ? p_local_dup * (1.0 + p_autocorr_alpha * r_loc[i])
                    : p_local_dup;
                if (p_eff > 1.0) p_eff = 1.0;
                if (R::runif(0.0, 1.0) < p_eff) {
                    triggered.push_back(i);
                }
            }
            for (int ti = (int)triggered.size() - 1; ti >= 0; ti--) {
                int idx = triggered[ti];
                int ck = (int)arr.size();
                int max_chunk = ck - idx;
                if (max_chunk <= 0) continue;

                int chunk_size = sample_chunk_size(local_dist, max_chunk);
                int start = idx;
                int end = std::min(start + chunk_size - 1, ck - 1);

                // Materialize source monomers before copying
                for (int j = start; j <= end; j++) {
                    materialize(arr[j], mu_total);
                }

                // Create copies
                std::vector<Monomer> chunk;
                chunk.reserve(end - start + 1);
                for (int j = start; j <= end; j++) {
                    chunk.emplace_back(arr[j].seq, arr[j].dir);
                }

                // Insert right after source (tandem)
                arr.insert(arr.begin() + end + 1,
                           std::make_move_iterator(chunk.begin()),
                           std::make_move_iterator(chunk.end()));
            }
        }

        // === DISTAL DUPLICATION ===
        // Per-ARRAY probability (not per-unit): models chromosome-level events
        // like unequal crossover. At most one event per generation.
        if (p_distal_dup > 0.0 && R::runif(0.0, 1.0) < p_distal_dup) {
            k = (int)arr.size();
            if (k > 0) {
                int idx = (int)(R::runif(0.0, 1.0) * k);
                if (idx >= k) idx = k - 1;
                int max_chunk = k - idx;

                int chunk_size = sample_chunk_size(distal_dist, max_chunk);
                int start = idx;
                int end = std::min(start + chunk_size - 1, k - 1);

                // Materialize source
                for (int j = start; j <= end; j++) {
                    materialize(arr[j], mu_total);
                }

                // Create copies, possibly inverted
                bool invert = R::runif(0.0, 1.0) < p_invert_distal;
                std::vector<Monomer> chunk;
                chunk.reserve(end - start + 1);
                if (invert) {
                    for (int j = end; j >= start; j--) {
                        chunk.emplace_back(arr[j].seq, -arr[j].dir);
                    }
                } else {
                    for (int j = start; j <= end; j++) {
                        chunk.emplace_back(arr[j].seq, arr[j].dir);
                    }
                }

                // Insert at random position (anywhere in array)
                int n_positions = k + 1;
                int insert_at = (int)(R::runif(0.0, 1.0) * n_positions);
                if (insert_at >= n_positions) insert_at = n_positions - 1;

                arr.insert(arr.begin() + insert_at,
                           std::make_move_iterator(chunk.begin()),
                           std::make_move_iterator(chunk.end()));
            }
        }

        // === CHUNK DELETION ===
        if (p_del_chunk > 0.0) {
            k = (int)arr.size();
            std::vector<int> triggered;
            for (int i = 0; i < k; i++) {
                if (R::runif(0.0, 1.0) < p_del_chunk) {
                    triggered.push_back(i);
                }
            }
            for (int ti = (int)triggered.size() - 1; ti >= 0; ti--) {
                int idx = triggered[ti];
                int ck = (int)arr.size();
                if (idx >= ck) continue;
                int max_chunk = ck - idx;
                if (max_chunk <= 0) continue;

                int chunk_size = sample_chunk_size(del_dist, max_chunk);
                int start = idx;
                int end = std::min(start + chunk_size - 1, ck - 1);

                arr.erase(arr.begin() + start, arr.begin() + end + 1);
            }
        }

        // === GENE CONVERSION (homogenization) ===
        // Non-allelic gene conversion: a recipient monomer's SEQUENCE is over-
        // written by a NEARBY donor's current sequence (double-strand-break repair
        // templated on a nearby repeat unit). No length change -- this homogenizes
        // (reduces divergence, creates identical copies) WITHOUT adding units.
        // Donor separation ~ conv_dist (local kernel); tract ~ conv_tract (short,
        // typically ~1 monomer). Recipient orientation (dir) is unchanged.
        if (p_conversion > 0.0) {
            k = (int)arr.size();
            if (k >= 2) {
                std::vector<int> recip;
                for (int i = 0; i < k; i++) {
                    if (R::runif(0.0, 1.0) < p_conversion) recip.push_back(i);
                }
                for (size_t ri = 0; ri < recip.size(); ri++) {
                    int i = recip[ri];
                    int off = sample_chunk_size(conv_dist, k - 1);   // donor separation >=1
                    int sign = (R::runif(0.0, 1.0) < 0.5) ? -1 : 1;
                    int donor = i + sign * off;
                    if (donor < 0)      donor = i + off;             // reflect at ends
                    if (donor > k - 1)  donor = i - off;
                    if (donor < 0 || donor > k - 1 || donor == i) continue;
                    int room = 1 + std::min(k - 1 - i, k - 1 - donor);
                    int tract = sample_chunk_size(conv_tract, room);
                    for (int j = 0; j < tract; j++) {
                        int rpos = i + j, dpos = donor + j;
                        if (rpos > k - 1 || dpos > k - 1) break;
                        materialize(arr[dpos], mu_total);   // resolve donor to current seq
                        arr[rpos].seq = arr[dpos].seq;      // overwrite recipient sequence
                        arr[rpos].pending_gens = 0;         // freshly homogenized copy
                    }
                }
            }
        }

        // === TERMINAL (LONG-RANGE) DELETION ===
        // Per-ARRAY probability (mirror of distal duplication): with prob
        // p_distal_del, remove a chunk of size ~distal_dist from the END of the
        // array (telomere-proximal truncation -- arrays lose their distal end).
        if (p_distal_del > 0.0 && R::runif(0.0, 1.0) < p_distal_del) {
            k = (int)arr.size();
            if (k > 0) {
                int chunk_size = sample_chunk_size(distal_dist, k);
                if (chunk_size > k) chunk_size = k;
                arr.erase(arr.end() - chunk_size, arr.end());
            }
        }

        // === INCREMENT PENDING MUTATIONS FOR ALL MONOMERS ===
        for (auto& m : arr) {
            m.pending_gens++;
        }

        l_vec.push_back((int)arr.size());
    }

    // === MATERIALIZE ALL REMAINING MUTATIONS ===
    int n_final = (int)arr.size();
    for (auto& m : arr) {
        materialize(m, mu_total);
    }

    // === CONVERT TO R OUTPUT ===
    const char bases[] = {'A', 'C', 'G', 'T'};
    CharacterVector seq_out(n_final);
    CharacterVector dir_out(n_final);

    for (int i = 0; i < n_final; i++) {
        std::string s(arr[i].seq.size(), ' ');
        for (int j = 0; j < (int)arr[i].seq.size(); j++) {
            s[j] = bases[arr[i].seq[j]];
        }
        seq_out[i] = s;
        dir_out[i] = (arr[i].dir > 0) ? "+" : "-";
    }

    return List::create(
        Named("seqs") = seq_out,
        Named("dirs") = dir_out,
        Named("l_vec") = wrap(l_vec),
        Named("total_gens") = t,
        Named("hit_hard_cap") = hit_hard_cap
    );
}

// ---------------------------------------------------------------------------
// Multi-chromosome (single-haplotype) simulation engine
// ---------------------------------------------------------------------------
// Holds K chromosome arrays in one genome. Each generation:
//   1. for each chromosome: apply the within-array moves (shared helper),
//   2. with prob p_translocation (per genome per generation), copy a contiguous
//      chunk from a random DONOR chromosome into a random RECIPIENT chromosome
//      at a random position, inverted w.p. p_invert_transloc.
// Each array grows toward its OWN target size; the genome loop runs until ALL
// arrays have reached their targets or max_t / hard cap is hit.
//
// One ancestral monomer seeds all chromosomes (shared origin), so without
// translocation, identical sequences shared ACROSS different chromosomes only
// arise from the (vanishing) chance of independent convergence — i.e. ~0.

// [[Rcpp::export]]
List sim_genome_cpp(IntegerVector ancestor_seq_r,
                    int K,
                    IntegerVector target_sizes_r,
                    int init_k0,
                    int max_t,
                    IntegerVector hard_caps_r,
                    double p_local_dup,
                    List local_dist,
                    double p_distal_dup,
                    List distal_dist,
                    double p_invert_distal,
                    double p_del_chunk,
                    List del_dist,
                    double p_distal_del,
                    double p_conversion,
                    List conv_dist,
                    List conv_tract,
                    double p_translocation,
                    List transloc_dist,
                    double p_invert_transloc,
                    double mu_total,
                    double p_autocorr_alpha,
                    int autocorr_window,
                    int autocorr_every,
                    int n_generations,
                    double size_band,
                    bool verbose) {

    if (K < 1) stop("K must be >= 1");
    if ((int)target_sizes_r.size() != K) stop("target_sizes must have length K");
    if ((int)hard_caps_r.size() != K) stop("hard_caps must have length K");

    // Convert ancestor: R 1-indexed -> C++ 0-indexed
    int L = ancestor_seq_r.size();
    std::vector<int> ancestor_seq(L);
    for (int i = 0; i < L; i++) ancestor_seq[i] = ancestor_seq_r[i] - 1;

    std::vector<int> target_sizes(K), hard_caps(K);
    for (int c = 0; c < K; c++) {
        target_sizes[c] = target_sizes_r[c];
        hard_caps[c]    = hard_caps_r[c];
    }

    // Initialize K arrays from the SAME ancestral monomer.
    // AGE MODE: seed each chromosome AT its target size (all identical ancestral
    // monomers, divergence 0). The growth transient from init_k0 -> target under
    // a per-unit duplication rate is geometric and would consume the entire fixed
    // age budget (10 -> 8000 takes tens of thousands of gens), starving the
    // steady-state turnover that actually sets divergence/redundancy. Seeding at
    // target puts every chromosome straight into the "fully grown young array"
    // state; the fixed n_generations of governed turnover then builds divergence
    // (~mu*age) while holding redundancy, exactly like a single grown array.
    // LEGACY MODE: seed init_k0 and grow (old behavior, back-compat).
    std::vector<std::vector<Monomer> > genome(K);
    bool age_mode = (n_generations > 0);
    for (int c = 0; c < K; c++) {
        int seed_n = age_mode ? target_sizes[c] : init_k0;
        genome[c].reserve(std::min(seed_n, hard_caps[c]) + 1000);
        for (int i = 0; i < seed_n; i++)
            genome[c].emplace_back(ancestor_seq, 1);
    }

    int t = 0;
    bool hit_hard_cap = false;
    long n_transloc = 0;

    // STEADY-STATE GENOME DYNAMICS (fixes freeze-and-over-age).
    // ---------------------------------------------------------------------
    // Two modes, selected by n_generations:
    //   * n_generations > 0  -> AGE-DRIVEN. The loop runs for EXACTLY
    //     n_generations generations (age is the divergence dial: with mu fixed,
    //     divergence ~= mu * age). A chromosome below its target band GROWS
    //     (dup-biased); once it reaches the band it enters STEADY-STATE TURNOVER
    //     -- local duplication keeps firing at p_local_dup (fresh identical
    //     twins) while a PER-CHROM SIZE GOVERNOR raises chunk deletion to match
    //     the duplication mass (net-zero size), with a proportional correction
    //     that pulls size back toward target (absorbing translocation overshoot).
    //     So redundancy is held at steady size (like a single growing-then-
    //     stopping array) instead of decaying as the clock ages a frozen array.
    //   * n_generations <= 0 -> LEGACY: grow until all arrays reach target then
    //     stop (the old behavior; over-target arrays freeze). Kept for back-compat.
    // The governor uses dup/del chunk MEANS to balance per-unit mass; in steady
    // state local dup + WITHIN-ARRAY DISTAL DUP + conversion fire, and the
    // per-unit chunk deletion is raised to absorb BOTH the local AND the distal
    // duplication mass (so within-array distal dup leaves an inferrable
    // off-diagonal signature in turnover instead of being suppressed). Only
    // long-range distal *deletion* (p_distal_del, a terminal truncation that
    // would fight the governor) stays off in steady state.
    int  t_limit  = age_mode ? n_generations : max_t;
    double band   = (size_band > 0.0) ? size_band : 0.10;  // +/- fraction of target
    double mean_local  = dist_mean(local_dist);
    double mean_del    = dist_mean(del_dist);
    double mean_distal = dist_mean(distal_dist);
    // proportional governor gain (per-unit deletion-prob correction per unit of
    // fractional size error); chosen so a full-band overshoot adds ~the balance
    // deletion again. Clamped downstream so it can never exceed a sane rate.
    double gov_gain = 0.5;

    while (t < t_limit) {
        // LEGACY mode termination: all arrays reached target -> stop.
        if (!age_mode) {
            bool all_done = true;
            for (int c = 0; c < K; c++) {
                if ((int)genome[c].size() < target_sizes[c]) { all_done = false; break; }
            }
            if (all_done) break;
        }

        t++;

        // Extinction guard per chromosome (re-seed) and cap flag (informational
        // only -- does NOT terminate the genome).
        for (int c = 0; c < K; c++) {
            if ((int)genome[c].size() > hard_caps[c]) hit_hard_cap = true;
            if ((int)genome[c].size() == 0) {
                // re-seed an extinct chromosome from ancestor
                for (int i = 0; i < init_k0; i++)
                    genome[c].emplace_back(ancestor_seq, 1);
            }
        }

        if (verbose && t % 500 == 0) {
            Rprintf("  Gen %d:", t);
            for (int c = 0; c < K; c++) Rprintf(" %d", (int)genome[c].size());
            Rprintf("\n");
        }
        if (t % 100 == 0) Rcpp::checkUserInterrupt();

        // === 1. WITHIN-ARRAY MOVES, per chromosome (governed) ===
        // Autocatalytic local-dup cadence gate (same scheme as sim_core_cpp):
        // recompute the local-redundancy weighting only every autocorr_every gens.
        bool ac_active_gen = (p_autocorr_alpha > 0.0) &&
                             (autocorr_every <= 1 || (t % autocorr_every) == 0);
        for (int c = 0; c < K; c++) {
            int    sz  = (int)genome[c].size();
            double tgt = (double)target_sizes[c];
            double lo_band = tgt * (1.0 - band);

            if (!age_mode) {
                // LEGACY: only grow until target, then freeze (old behavior).
                if (sz >= target_sizes[c]) continue;
                if (sz >= hard_caps[c])    continue;
                apply_within_array_moves(genome[c],
                                         p_local_dup, local_dist,
                                         p_distal_dup, distal_dist, p_invert_distal,
                                         p_del_chunk, del_dist, p_distal_del,
                                         p_conversion, conv_dist, conv_tract,
                                         mu_total,
                                         p_autocorr_alpha, autocorr_window,
                                         ac_active_gen);
                continue;
            }

            // AGE mode: governed growth-then-turnover.
            if (sz < (int)lo_band) {
                // BELOW BAND -> GROWTH (dup-biased): base rates, full move set.
                // Cap-respect: don't let a runaway distal event blow past hard cap.
                if (sz >= hard_caps[c]) continue;
                apply_within_array_moves(genome[c],
                                         p_local_dup, local_dist,
                                         p_distal_dup, distal_dist, p_invert_distal,
                                         p_del_chunk, del_dist, p_distal_del,
                                         p_conversion, conv_dist, conv_tract,
                                         mu_total,
                                         p_autocorr_alpha, autocorr_window,
                                         ac_active_gen);
            } else {
                // AT/ABOVE BAND -> STEADY-STATE TURNOVER.
                // Keep local dup AND within-array distal dup firing (fresh twins /
                // off-diagonal copies); raise the per-unit chunk-deletion prob to
                // MATCH the TOTAL duplication mass (local + distal), plus a
                // proportional pull toward target (delete-biased when above,
                // dup-biased when below within the band). Long-range distal
                // DELETION (terminal truncation) stays off so the governor stays
                // in control of the size.
                double p_local_eff = p_local_dup;
                // net-zero balance per generation:
                //   local dup mass  = p_local * sz * mean_local   (per-unit)
                //   distal dup mass = p_distal_dup * mean_distal   (per-array event)
                //   deletion mass   = p_del * sz * mean_del        (per-unit)
                // => p_del_balance = p_local*mean_local/mean_del
                //                  + p_distal_dup*mean_distal/(sz*mean_del)
                double p_del_balance = 0.0;
                if (mean_del > 0.0) {
                    p_del_balance = p_local_eff * mean_local / mean_del;
                    if (sz > 0)
                        p_del_balance += p_distal_dup * mean_distal
                                         / ((double)sz * mean_del);
                }
                // fractional size error (positive when above target)
                double err = (tgt > 0.0) ? ((double)sz - tgt) / tgt : 0.0;
                double p_del_eff = p_del_balance + gov_gain * err / mean_del;
                if (p_del_eff < 0.0) p_del_eff = 0.0;
                if (p_del_eff > 0.5) p_del_eff = 0.5;
                apply_within_array_moves(genome[c],
                                         p_local_eff, local_dist,
                                         /*p_distal_dup=*/p_distal_dup, distal_dist, p_invert_distal,
                                         p_del_eff, del_dist, /*p_distal_del=*/0.0,
                                         p_conversion, conv_dist, conv_tract,
                                         mu_total,
                                         p_autocorr_alpha, autocorr_window,
                                         ac_active_gen);
            }
        }

        // === 2. INTER-CHROMOSOMAL TRANSLOCATION (per genome) ===
        if (K >= 2 && p_translocation > 0.0 &&
            R::runif(0.0, 1.0) < p_translocation) {

            // pick donor and a distinct recipient
            int donor = (int)(R::runif(0.0, 1.0) * K);
            if (donor >= K) donor = K - 1;
            int recip = (int)(R::runif(0.0, 1.0) * K);
            if (recip >= K) recip = K - 1;
            if (recip == donor) recip = (recip + 1) % K;

            int dk = (int)genome[donor].size();
            // A recipient already at/over its hard cap can still DONATE but must
            // not GROW unbounded -- skip insertion into an over-cap recipient.
            if (dk > 0 && (int)genome[recip].size() < hard_caps[recip]) {
                int idx = (int)(R::runif(0.0, 1.0) * dk);
                if (idx >= dk) idx = dk - 1;
                int max_chunk = dk - idx;
                int chunk_size = sample_chunk_size(transloc_dist, max_chunk);
                int start = idx;
                int end = std::min(start + chunk_size - 1, dk - 1);

                // resolve donor monomers to current sequence before copying
                for (int j = start; j <= end; j++)
                    materialize(genome[donor][j], mu_total);

                bool invert = R::runif(0.0, 1.0) < p_invert_transloc;
                std::vector<Monomer> chunk;
                chunk.reserve(end - start + 1);
                if (invert) {
                    for (int j = end; j >= start; j--)
                        chunk.emplace_back(genome[donor][j].seq,
                                           -genome[donor][j].dir);
                } else {
                    for (int j = start; j <= end; j++)
                        chunk.emplace_back(genome[donor][j].seq,
                                           genome[donor][j].dir);
                }

                // insert at random position in recipient
                int rk = (int)genome[recip].size();
                int n_positions = rk + 1;
                int insert_at = (int)(R::runif(0.0, 1.0) * n_positions);
                if (insert_at >= n_positions) insert_at = n_positions - 1;
                genome[recip].insert(genome[recip].begin() + insert_at,
                                     std::make_move_iterator(chunk.begin()),
                                     std::make_move_iterator(chunk.end()));
                n_transloc++;
            }
        }

        // === 3. INCREMENT PENDING MUTATIONS FOR ALL MONOMERS, ALL CHROMS ===
        for (int c = 0; c < K; c++)
            for (auto& m : genome[c]) m.pending_gens++;
    }

    // === MATERIALIZE ALL REMAINING MUTATIONS ===
    for (int c = 0; c < K; c++)
        for (auto& m : genome[c]) materialize(m, mu_total);

    // === BUILD FLAT OUTPUT: chrom, num, bponly, dir ===
    int n_total = 0;
    for (int c = 0; c < K; c++) n_total += (int)genome[c].size();

    IntegerVector chrom_out(n_total);
    IntegerVector num_out(n_total);
    CharacterVector seq_out(n_total);
    CharacterVector dir_out(n_total);

    const char bases[] = {'A', 'C', 'G', 'T'};
    int row = 0;
    for (int c = 0; c < K; c++) {
        int nc = (int)genome[c].size();
        for (int i = 0; i < nc; i++) {
            std::string s(genome[c][i].seq.size(), ' ');
            for (int j = 0; j < (int)genome[c][i].seq.size(); j++)
                s[j] = bases[genome[c][i].seq[j]];
            chrom_out[row] = c + 1;        // 1..K
            num_out[row]   = i + 1;        // 1..nc within chromosome
            seq_out[row]   = s;
            dir_out[row]   = (genome[c][i].dir > 0) ? "+" : "-";
            row++;
        }
    }

    IntegerVector sizes_out(K);
    for (int c = 0; c < K; c++) sizes_out[c] = (int)genome[c].size();

    return List::create(
        Named("chrom")        = chrom_out,
        Named("num")          = num_out,
        Named("bponly")       = seq_out,
        Named("dir")          = dir_out,
        Named("array_sizes")  = sizes_out,
        Named("total_gens")   = t,
        Named("n_transloc")   = (double)n_transloc,
        Named("hit_hard_cap") = hit_hard_cap
    );
}
