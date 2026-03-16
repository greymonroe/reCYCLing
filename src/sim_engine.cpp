#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

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
                  double mu_total,
                  bool verbose) {

    // Convert ancestor sequence: R 1-indexed -> C++ 0-indexed
    int L = ancestor_seq_r.size();
    std::vector<int> ancestor_seq(L);
    for (int i = 0; i < L; i++) {
        ancestor_seq[i] = ancestor_seq_r[i] - 1;
    }

    // Initialize array with init_k0 identical monomers
    std::vector<Monomer> arr;
    arr.reserve(max_units + 1000);
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
        if (k == 0) break;

        if (verbose && t % 500 == 0) {
            Rprintf("  Gen %d - units: %d\n", t, k);
        }
        if (t % 100 == 0) Rcpp::checkUserInterrupt();

        // === LOCAL DUPLICATION ===
        if (p_local_dup > 0.0) {
            k = (int)arr.size();
            std::vector<int> triggered;
            for (int i = 0; i < k; i++) {
                if (R::runif(0.0, 1.0) < p_local_dup) {
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
        if (p_distal_dup > 0.0) {
            k = (int)arr.size();
            std::vector<int> triggered;
            for (int i = 0; i < k; i++) {
                if (R::runif(0.0, 1.0) < p_distal_dup) {
                    triggered.push_back(i);
                }
            }
            for (int ti = (int)triggered.size() - 1; ti >= 0; ti--) {
                int idx = triggered[ti];
                int ck = (int)arr.size();
                if (idx >= ck) continue;
                int max_chunk = ck - idx;
                if (max_chunk <= 0) continue;

                int chunk_size = sample_chunk_size(distal_dist, max_chunk);
                int start = idx;
                int end = std::min(start + chunk_size - 1, ck - 1);

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
                int n_positions = ck + 1;
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
