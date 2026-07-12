// RNAmotifs - Bootstrap FDR estimation (C++ with OpenMP)
// Copyright (C) 2014-2026 Matteo Cereda
// SPDX-License-Identifier: GPL-2.0-or-later
//
// Replaces the R bootstrap-FDR.R script with a fast, parallelised
// implementation. Produces identical TSV output.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// ===========================================================================
// Data structures
// ===========================================================================

struct TetData {
    string name;
    vector<int> type;     // -1, 0, 1
    vector<int> hits_r1;  // 0 or 1
    vector<int> hits_r2;
    vector<int> hits_r3;
};

// ===========================================================================
// Fisher exact test (one-sided, alternative = "greater")
// ===========================================================================

static double fisher_greater(int a, int b, int c, int d) {
    //  Contingency table:
    //         hit   no-hit
    //  group   a      b     -> row1 = a+b
    //  ctrl    c      d     -> row2 = c+d
    //         ---    ---
    //  col1  a+c    b+d     -> total = n
    //
    //  P(X >= a) under Hypergeometric(n, row1, col1)
    int n    = a + b + c + d;
    int row1 = a + b;
    int col1 = a + c;
    int row2 = c + d;

    if (n == 0) return 1.0;

    int x_max = min(row1, col1);

    // log-denominator: log C(n, row1)
    double log_denom = lgamma(n + 1) - lgamma(row1 + 1) - lgamma(n - row1 + 1);

    double pval = 0.0;
    for (int x = a; x <= x_max; ++x) {
        double log_num = lgamma(col1 + 1) - lgamma(x + 1) - lgamma(col1 - x + 1)
                       + lgamma(n - col1 + 1) - lgamma(row1 - x + 1)
                       - lgamma(n - col1 - row1 + x + 1);
        pval += exp(log_num - log_denom);
    }
    return min(pval, 1.0);
}

// ===========================================================================
// Benjamini-Hochberg p-value adjustment (matches R's p.adjust(method="BH"))
// ===========================================================================

static void bh_adjust(vector<double> &pvals) {
    int n = (int)pvals.size();
    if (n <= 1) return;

    // Sort indices by p-value descending
    vector<int> idx(n);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(),
         [&](int a, int b) { return pvals[a] > pvals[b]; });

    double cummin = 1.0;
    for (int i = 0; i < n; ++i) {
        int rank = n - i;  // ascending rank
        double adj = min(pvals[idx[i]] * (double)n / (double)rank, 1.0);
        cummin = min(cummin, adj);
        pvals[idx[i]] = cummin;
    }
}

// ===========================================================================
// Compute Fisher p-values for all tetramers (6 values each)
// Applies BH adjustment per column.
// ===========================================================================

static vector<double> compute_fisher_bh(
        const vector<int> &indices,
        const vector<TetData> &tets) {

    int n_tet = (int)tets.size();
    int n_exons = (int)indices.size();
    // Output: flattened n_tet * 6 (columns: r1enh, r1sil, r2enh, r2sil, r3enh, r3sil)
    vector<double> pvals(n_tet * 6);

    for (int t = 0; t < n_tet; ++t) {
        const auto &td = tets[t];
        const vector<int> *hits[3] = {&td.hits_r1, &td.hits_r2, &td.hits_r3};

        for (int r = 0; r < 3; ++r) {
            int enh_p = 0, enh_a = 0, ctl_p = 0, ctl_a = 0;
            int sil_p = 0, sil_a = 0;

            // Match R behavior: type from ORIGINAL position i,
            // hits from BOOTSTRAPPED index indices[i]
            for (int i = 0; i < n_exons; ++i) {
                int ty = td.type[i];            // original position
                int h  = (*hits[r])[indices[i]]; // bootstrapped index
                if (ty == 1) {
                    if (h > 0) enh_p++; else enh_a++;
                } else if (ty == 0) {
                    if (h > 0) ctl_p++; else ctl_a++;
                } else if (ty == -1) {
                    if (h > 0) sil_p++; else sil_a++;
                }
            }

            pvals[t * 6 + r * 2]     = fisher_greater(enh_p, ctl_p, enh_a, ctl_a);
            pvals[t * 6 + r * 2 + 1] = fisher_greater(sil_p, ctl_p, sil_a, ctl_a);
        }
    }

    // BH adjustment per column (6 columns, each of length n_tet)
    for (int col = 0; col < 6; ++col) {
        vector<double> column(n_tet);
        for (int t = 0; t < n_tet; ++t)
            column[t] = pvals[t * 6 + col];
        bh_adjust(column);
        for (int t = 0; t < n_tet; ++t)
            pvals[t * 6 + col] = column[t];
    }

    return pvals;
}

// ===========================================================================
// Read tetramer count files
// ===========================================================================

static TetData read_count_file(const string &filepath, const string &tet_name) {
    TetData td;
    td.name = tet_name;
    ifstream f(filepath);
    if (!f.is_open()) {
        cerr << "Error: cannot open " << filepath << endl;
        return td;
    }
    string line;
    getline(f, line); // skip header
    while (getline(f, line)) {
        istringstream ss(line);
        string tok;
        vector<string> cols;
        while (getline(ss, tok, '\t')) cols.push_back(tok);
        if (cols.size() < 6) continue;
        td.type.push_back(stoi(cols[2]));
        td.hits_r1.push_back(stoi(cols[3]));
        td.hits_r2.push_back(stoi(cols[4]));
        td.hits_r3.push_back(stoi(cols[5]));
    }
    return td;
}

static vector<string> read_filelist(const string &path) {
    vector<string> files;
    ifstream f(path);
    if (!f.is_open()) return files;
    string line;
    while (getline(f, line)) {
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (!line.empty() && line.find("filelist") == string::npos)
            files.push_back(line);
    }
    return files;
}

// ===========================================================================
// Main
// ===========================================================================

int main(int argc, char *argv[]) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0]
             << " <results_dir> <run_name> <n_bootstraps> [n_cores] [--estimate]" << endl;
        return 1;
    }

    string results_dir = argv[1];
    string run_name    = argv[2];
    int    n_boot      = stoi(argv[3]);
    int    n_cores     = (argc >= 5 && string(argv[4]) != "--estimate") ? stoi(argv[4]) : 1;
    bool   estimate    = false;
    for (int i = 1; i < argc; ++i)
        if (string(argv[i]) == "--estimate") estimate = true;

#ifdef _OPENMP
    omp_set_num_threads(n_cores);
#endif

    string pp = results_dir + "/" + run_name + "/";

    // ── Load tetramer count data ────────────────────────────────────────────
    vector<TetData> tets;
    vector<string> tet_names;

    for (const string &subdir : {"r/", "nr/"}) {
        string flist_path = pp + subdir + "filelist_count.tsv";
        vector<string> files = read_filelist(flist_path);
        for (const auto &fname : files) {
            string tet_name = fname.substr(0, fname.find('_'));
            string fpath = pp + subdir + fname;
            TetData td = read_count_file(fpath, tet_name);
            if (!td.type.empty()) {
                tets.push_back(move(td));
                tet_names.push_back(tet_name);
            }
        }
    }

    int n_tet = (int)tets.size();
    if (n_tet == 0) {
        cerr << "Error: no tetramer counts loaded" << endl;
        return 1;
    }
    int n_exons = (int)tets[0].type.size();

    cerr << "Loaded " << n_tet << " tetramers, " << n_exons << " exons" << endl;

    // ── Original Fisher p-values ────────────────────────────────────────────
    vector<int> all_idx(n_exons);
    iota(all_idx.begin(), all_idx.end(), 0);

    cerr << "Computing Fisher exact tests..." << endl;
    vector<double> pFisher = compute_fisher_bh(all_idx, tets);

    // ── Estimate mode ───────────────────────────────────────────────────────
    if (estimate) {
        // Run 3 bootstrap iterations and print time per iteration
        auto t0 = chrono::steady_clock::now();
        mt19937 rng(30580);
        uniform_int_distribution<int> dist(0, n_exons - 1);
        for (int b = 0; b < 3; ++b) {
            vector<int> boot_idx(n_exons);
            for (int i = 0; i < n_exons; ++i)
                boot_idx[i] = dist(rng);
            compute_fisher_bh(boot_idx, tets);
        }
        auto t1 = chrono::steady_clock::now();
        double elapsed = chrono::duration<double>(t1 - t0).count();
        // Print: seconds per iteration (will be parsed by Python)
        cout << fixed << setprecision(3) << elapsed / 3.0 << "\n";
        return 0;
    }

    // ── Bootstrap ───────────────────────────────────────────────────────────
    cerr << "Bootstrapping (" << n_boot << " iterations, "
         << n_cores << " cores)..." << endl;

    int n_vals = n_tet * 6;
    vector<long> global_counts(n_vals, 0);

    auto t_start = chrono::steady_clock::now();
    int progress_done = 0;

    #pragma omp parallel
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        vector<long> local_counts(n_vals, 0);

        #pragma omp for schedule(static)
        for (int b = 0; b < n_boot; ++b) {
            // Option 1 (core-count-independent): seed RNG per iteration, not per thread.
            mt19937 rng(30580u + (unsigned long)b * 2654435761ul);
            uniform_int_distribution<int> dist(0, n_exons - 1);
            // Resample
            vector<int> boot_idx(n_exons);
            for (int i = 0; i < n_exons; ++i)
                boot_idx[i] = dist(rng);

            // Fisher + BH
            vector<double> boot_p = compute_fisher_bh(boot_idx, tets);

            // Compare
            for (int j = 0; j < n_vals; ++j)
                if (boot_p[j] <= pFisher[j])
                    local_counts[j]++;

            // Progress (only thread 0 reports)
            if (tid == 0) {
                #pragma omp atomic
                progress_done++;
                int done;
                #pragma omp atomic read
                done = progress_done;
                // Report every 1% or every iteration if < 100 total
                if (n_boot < 100 || done % max(1, n_boot / 100) == 0)
                    cerr << "PROGRESS " << done << " " << n_boot << endl;
            } else {
                #pragma omp atomic
                progress_done++;
            }
        }

        // Merge
        #pragma omp critical
        for (int j = 0; j < n_vals; ++j)
            global_counts[j] += local_counts[j];
    }

    auto t_end = chrono::steady_clock::now();
    double elapsed = chrono::duration<double>(t_end - t_start).count();
    cerr << "PROGRESS " << n_boot << " " << n_boot << endl;
    cerr << "Bootstrap completed in " << fixed << setprecision(1)
         << elapsed << " seconds" << endl;

    // ── Empirical p-values ──────────────────────────────────────────────────
    vector<double> pEmpirical(n_vals);
    for (int j = 0; j < n_vals; ++j)
        pEmpirical[j] = (1.0 + global_counts[j]) / (1.0 + n_boot);

    // ── Write output TSV ────────────────────────────────────────────────────
    string out_path = pp + "bootstrap_" + to_string(n_boot) + ".tsv";
    ofstream out(out_path);
    if (!out.is_open()) {
        cerr << "Error: cannot write " << out_path << endl;
        return 1;
    }

    const char *col_names[] = {
        "r1enh", "r1sil", "r2enh", "r2sil", "r3enh", "r3sil"
    };

    out << "tetramer";
    for (int c = 0; c < 6; ++c)
        out << "\t" << col_names[c] << "_pFis"
            << "\t" << col_names[c] << "_pEmp";
    out << "\n";

    for (int t = 0; t < n_tet; ++t) {
        out << tet_names[t];
        for (int c = 0; c < 6; ++c) {
            int idx = t * 6 + c;
            out << "\t" << fixed << setprecision(6) << pFisher[idx]
                << "\t" << fixed << setprecision(6) << pEmpirical[idx];
        }
        out << "\n";
    }
    out.close();

    cerr << "Results saved to " << out_path << endl;
    return 0;
}
