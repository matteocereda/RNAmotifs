// RNAmotifs-MaRs - Association score computation (C++ with OpenMP)
// Copyright (C) 2014-2026 Matteo Cereda
// SPDX-License-Identifier: GPL-2.0-or-later
//
// Replaces compute_association_scores.R with a fast C++ implementation.
// Computes SCORE1 (Signal Recovery Rate) and SCORE2 (Cosine Similarity)
// between RNAmotifs enrichment profiles and eCLIP binding profiles.
//
// Usage:
//   rnamotifs_mars_score -i <input_exons> -n <name> -d <sweep_dir>
//       -r <mars_ref_dir> -c <cell_line> -p <cores> -e <eclip_dir>
//       -o <output_dir> [--in-exon N] [--in-intron N] [-b N]
//       [--summary] [--compute-peak --mars-exons-dir <dir>]

#include <algorithm>
#include <random>
#include <cmath>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// ═══════════════════════════════════════════════════════════════════════════
// Data structures
// ═══════════════════════════════════════════════════════════════════════════

struct Exon {
    int id;
    int second_id;
    string chrom;
    string strand;
    int v5, v6, v7, v8; // upstream_ee, exon_start, exon_end, downstream_es
    double dirank;
    int type; // 1=enh, -1=sil, 0=ctrl
    string event_type; // "SE" or "RI"
};

struct Peak {
    string chrom;
    int start, end;
};

struct OptimalParams {
    string rbp;
    int hw, ew;
    string params_key; // "hw_5_ew_100"
};

// ═══════════════════════════════════════════════════════════════════════════
// Fisher exact test (one-sided, alternative = "greater")
// ═══════════════════════════════════════════════════════════════════════════

static double fisher_greater(int a, int b, int c, int d) {
    int n = a + b + c + d;
    int row1 = a + b;
    int col1 = a + c;
    if (n == 0) return 1.0;
    int x_max = min(row1, col1);
    double log_denom = lgamma(n+1) - lgamma(row1+1) - lgamma(n-row1+1);
    double pval = 0.0;
    for (int x = a; x <= x_max; ++x) {
        double log_num = lgamma(col1+1) - lgamma(x+1) - lgamma(col1-x+1)
                       + lgamma(n-col1+1) - lgamma(row1-x+1)
                       - lgamma(n-col1-row1+x+1);
        pval += exp(log_num - log_denom);
    }
    return min(pval, 1.0);
}

// ═══════════════════════════════════════════════════════════════════════════
// Cosine similarity
// ═══════════════════════════════════════════════════════════════════════════

static double cosine_similarity(const vector<double> &a, const vector<double> &b) {
    if (a.size() != b.size() || a.empty()) return 0.0;
    double dot = 0, na = 0, nb = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        dot += a[i] * b[i];
        na += a[i] * a[i];
        nb += b[i] * b[i];
    }
    double denom = sqrt(na) * sqrt(nb);
    return denom > 0 ? dot / denom : 0.0;
}

// ═══════════════════════════════════════════════════════════════════════════
// File readers
// ═══════════════════════════════════════════════════════════════════════════

static vector<Exon> read_input_exons(const string &path) {
    vector<Exon> exons;
    ifstream f(path);
    if (!f.is_open()) { cerr << "Error: cannot open " << path << endl; return exons; }
    string line;
    while (getline(f, line)) {
        if (line.empty()) continue;
        // Split by ; to handle optional 10th field
        vector<string> fields;
        istringstream iss(line);
        string tok;
        while (getline(iss, tok, ';')) fields.push_back(tok);
        if (fields.size() < 9) continue;

        Exon e;
        e.id = stoi(fields[0]);
        e.second_id = stoi(fields[1]);
        e.chrom = fields[2];
        e.strand = fields[3];
        e.v5 = stoi(fields[4]);
        e.v6 = stoi(fields[5]);
        e.v7 = stoi(fields[6]);
        e.v8 = stoi(fields[7]);
        e.dirank = stod(fields[8]);
        e.event_type = (fields.size() > 9 && !fields[9].empty()) ? fields[9] : "SE";
        // dIRank classification: align with core.cpp (`>= dIRO`, dIRO=1.0) so the
        // canonical categorical dIRank {-1,0,+1} is scored as regulated. The
        // previous strict `> 1.0` rejected ±1 -> treated regulated exons as
        // control (zeroed SRR/PEAK on {-1,0,1} input). Safe for whole-exon input
        // (|dirank| >= 2 there, nothing sits exactly at ±1).
        if (e.dirank >= 1.0) e.type = 1;
        else if (e.dirank <= -1.0) e.type = -1;
        else e.type = 0;
        exons.push_back(e);
    }
    return exons;
}

static vector<Peak> read_peaks(const string &path) {
    vector<Peak> peaks;
    ifstream f(path);
    if (!f.is_open()) return peaks;
    string line;
    while (getline(f, line)) {
        if (line.empty()) continue;
        istringstream ss(line);
        Peak p;
        ss >> p.chrom >> p.start >> p.end;
        peaks.push_back(p);
    }
    return peaks;
}

static map<string, double> read_auc_tsv(const string &path) {
    map<string, double> auc;
    ifstream f(path);
    if (!f.is_open()) { cerr << "Warning: cannot open AUC file " << path << endl; return auc; }
    string line;
    getline(f, line); // header
    while (getline(f, line)) {
        istringstream ss(line);
        string name; double val;
        ss >> name >> val;
        auc[name] = val;
    }
    return auc;
}

static vector<OptimalParams> read_optimal_params(const string &path,
                                                  const vector<string> &rbps) {
    vector<OptimalParams> params;
    set<string> rbp_set(rbps.begin(), rbps.end());
    ifstream f(path);
    if (!f.is_open()) { cerr << "Error: cannot open " << path << endl; return params; }
    string line;
    getline(f, line); // header
    while (getline(f, line)) {
        // CSV with quotes: "params","true_rbp","hw","ew"
        for (char &c : line) if (c == '"') c = ' ';
        istringstream ss(line);
        string params_key, rbp, hw_s, ew_s;
        char comma;
        ss >> params_key >> comma >> rbp >> comma >> hw_s >> comma >> ew_s;
        // Trim whitespace
        auto trim = [](string s) {
            s.erase(0, s.find_first_not_of(" \t"));
            s.erase(s.find_last_not_of(" \t") + 1);
            return s;
        };
        params_key = trim(params_key);
        rbp = trim(rbp);
        hw_s = trim(hw_s);
        ew_s = trim(ew_s);

        if (rbp_set.count(rbp)) {
            OptimalParams op;
            op.rbp = rbp;
            op.hw = stoi(hw_s);
            op.ew = stoi(ew_s);
            op.params_key = "hw_" + to_string(op.hw) + "_ew_" + to_string(op.ew);
            params.push_back(op);
        }
    }
    return params;
}

static vector<string> read_enriched_tetramers(const string &path) {
    vector<string> tets;
    ifstream f(path);
    if (!f.is_open()) return tets;
    string line;
    while (getline(f, line)) {
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (!line.empty()) tets.push_back(line);
    }
    return tets;
}

// Read regions.csv: "","R1","R2","R3"  -> map<tetramer, {R1,R2,R3}>
struct RegionInfo { string R1, R2, R3; };
static map<string, RegionInfo> read_regions_csv(const string &path) {
    map<string, RegionInfo> regions;
    ifstream f(path);
    if (!f.is_open()) return regions;
    string line;
    getline(f, line); // header
    while (getline(f, line)) {
        // Remove quotes
        for (char &c : line) if (c == '"') c = ' ';
        istringstream ss(line);
        string tet, r1, r2, r3;
        char comma;
        ss >> tet >> comma >> r1 >> comma >> r2 >> comma >> r3;
        auto trim = [](string s) {
            s.erase(0, s.find_first_not_of(" \t"));
            s.erase(s.find_last_not_of(" \t") + 1);
            return s;
        };
        tet = trim(tet); r1 = trim(r1); r2 = trim(r2); r3 = trim(r3);
        regions[tet] = {r1, r2, r3};
    }
    return regions;
}

// Read score CSV (e-*.csv or s-*.csv): header = tetramer names, rows = position scores
struct ScoreMatrix {
    vector<string> tetramers;
    vector<double> data; // [n_positions x n_tets]
    int n_pos, n_tets;
};

static ScoreMatrix read_score_csv(const string &path) {
    ScoreMatrix sm;
    ifstream f(path);
    if (!f.is_open()) { sm.n_pos = 0; sm.n_tets = 0; return sm; }
    string line;
    // Header
    getline(f, line);
    istringstream hss(line);
    string tok;
    while (getline(hss, tok, ',')) {
        tok.erase(tok.find_last_not_of(" \t\r\n\"") + 1);
        tok.erase(0, tok.find_first_not_of(" \t\""));
        if (!tok.empty()) sm.tetramers.push_back(tok);
    }
    sm.n_tets = (int)sm.tetramers.size();
    // Data
    while (getline(f, line)) {
        istringstream ss(line);
        for (int i = 0; i < sm.n_tets; ++i) {
            double v;
            ss >> v;
            sm.data.push_back(v);
            char c; ss >> c; // comma
        }
    }
    sm.n_pos = sm.n_tets > 0 ? (int)sm.data.size() / sm.n_tets : 0;
    return sm;
}

// Read region_count.tsv: myRID, rowID, type, hits_region1, hits_region2, hits_region3
struct RegionCount {
    int myRID, type;
    int hits[3];
};

static vector<RegionCount> read_region_count(const string &path) {
    vector<RegionCount> counts;
    ifstream f(path);
    if (!f.is_open()) return counts;
    string line;
    getline(f, line); // header
    while (getline(f, line)) {
        istringstream ss(line);
        RegionCount rc;
        int rowID;
        ss >> rc.myRID >> rowID >> rc.type >> rc.hits[0] >> rc.hits[1] >> rc.hits[2];
        counts.push_back(rc);
    }
    return counts;
}

// ═══════════════════════════════════════════════════════════════════════════
// SEeCLIPpeaks: compute eCLIP binding overlap at each position around exons
// Returns: overlap matrix [n_exons x n_positions_total]
// where n_positions_total = 4 * (in_exon + in_intron + 1)
// ═══════════════════════════════════════════════════════════════════════════

static vector<int> compute_eclip_overlap(
        const vector<Exon> &exons,
        const vector<Peak> &peaks,
        int in_exon, int in_intron, int n_cores) {

    int region_len = in_exon + in_intron + 1;
    int n_pos_total = 4 * region_len;
    int n_exons = (int)exons.size();

    // Build chromosome-indexed hash set of peak positions for O(1) lookup
    // Each peak spans start..end, so insert all covered positions
    unordered_map<string, unordered_set<int>> chr_peak_set;
    for (const auto &p : peaks) {
        auto &s = chr_peak_set[p.chrom];
        for (int pos = p.start; pos <= p.end; ++pos)
            s.insert(pos);
    }

    vector<int> overlap(n_exons * n_pos_total, 0);

    #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 100)
    for (int ei = 0; ei < n_exons; ++ei) {
        const Exon &e = exons[ei];
        auto it = chr_peak_set.find(e.chrom);
        if (it == chr_peak_set.end()) continue;
        const auto &peak_set = it->second;

        // For each position in the 4 regions, check if it's in the peak set
        for (int p = 0; p < region_len; ++p) {
            int offset_ei = -in_exon + p;
            int offset_ie = -in_intron + p;

            int positions[4];
            if (e.event_type == "RI") {
                // IR: R1=upstream exon (in_exon from 5'SS), R2=retained intron (in_intron from 5'SS),
                //     R3=retained intron (in_intron from 3'SS), R4=downstream exon (in_exon from 3'SS)
                positions[0] = e.v6 + offset_ei;   // R1: 5'SS into upstream exon
                positions[1] = e.v6 + offset_ie;   // R2: 5'SS into retained intron
                positions[2] = e.v7 + offset_ei;   // R3: 3'SS into retained intron
                positions[3] = e.v7 + offset_ie;   // R4: 3'SS into downstream exon
            } else {
                // SE: original cassette exon positions
                positions[0] = e.v5 + offset_ei;   // R1
                positions[1] = e.v6 + offset_ie;   // R2
                positions[2] = e.v7 + offset_ei;   // R3
                positions[3] = e.v8 + offset_ie;   // R4
            }

            for (int r = 0; r < 4; ++r) {
                if (peak_set.count(positions[r]))
                    overlap[ei * n_pos_total + r * region_len + p] = 1;
            }
        }

        // Handle minus strand: reverse the overlap row
        if (e.strand == "-") {
            int *row = &overlap[ei * n_pos_total];
            reverse(row, row + n_pos_total);
        }
    }
    return overlap;
}

// ═══════════════════════════════════════════════════════════════════════════
// Find which exons contain a tetramer hit in a specific region
// ═══════════════════════════════════════════════════════════════════════════

static vector<int> get_exons_with_tetramer(
        const string &tet, const string &sweep_folder,
        int region_num, const vector<Exon> &exons) {

    vector<int> matching_exon_ids;

    // Find region_count file for this tetramer
    string rc_path;
    for (const string &sub : {"nr/", "r/"}) {
        string path = sweep_folder + "/" + sub + tet + "_region_count.tsv";
        ifstream test(path);
        if (test.is_open()) { rc_path = path; break; }
    }
    if (rc_path.empty()) return matching_exon_ids;

    vector<RegionCount> counts = read_region_count(rc_path);
    for (const auto &rc : counts) {
        if (rc.type != 0 && rc.hits[region_num - 1] > 0) {
            matching_exon_ids.push_back(rc.myRID);
        }
    }
    return matching_exon_ids;
}

// ═══════════════════════════════════════════════════════════════════════════
// Find eCLIP peak file for an RBP
// ═══════════════════════════════════════════════════════════════════════════

static string find_eclip_file(const string &eclip_dir, const string &rbp) {
    DIR *dir = opendir(eclip_dir.c_str());
    if (!dir) return "";
    struct dirent *ent;
    string match;
    while ((ent = readdir(dir)) != nullptr) {
        string fname = ent->d_name;
        if (fname.find(rbp) != string::npos &&
            fname.find("ordered_merged.bed") != string::npos) {
            match = eclip_dir + "/" + fname;
            break;
        }
    }
    closedir(dir);
    return match;
}

// ═══════════════════════════════════════════════════════════════════════════
// Write score matrices as TSV (portable replacement for RDS)
// ═══════════════════════════════════════════════════════════════════════════

static void write_score_tsv(const string &path,
                            const map<string, map<string, double>> &scores,
                            const string &score_type) {
    ofstream f(path);
    if (!f.is_open()) { cerr << "Error: cannot write " << path << endl; return; }

    // Collect all tetramers
    set<string> all_tets;
    for (const auto &rbp_entry : scores)
        for (const auto &tet_entry : rbp_entry.second)
            all_tets.insert(tet_entry.first);

    // Header
    f << "rbp";
    for (const auto &t : all_tets) f << "\t" << t;
    f << "\n";

    // Data
    for (const auto &rbp_entry : scores) {
        f << rbp_entry.first;
        for (const auto &t : all_tets) {
            auto it = rbp_entry.second.find(t);
            f << "\t" << (it != rbp_entry.second.end() ? it->second : 0.0);
        }
        f << "\n";
    }
    cerr << "[mars_score] Wrote " << path << " (" << score_type << ")" << endl;
}

// ═══════════════════════════════════════════════════════════════════════════
// Find a file in a directory matching a prefix (e.g. "e-" or "s-")
// Returns full path or empty string
// ═══════════════════════════════════════════════════════════════════════════

static string find_file_by_prefix(const string &dir, const string &prefix, const string &suffix) {
    DIR *d = opendir(dir.c_str());
    if (!d) return "";
    struct dirent *ent;
    string match;
    while ((ent = readdir(d)) != nullptr) {
        string fname = ent->d_name;
        if (fname.size() >= prefix.size() + suffix.size() &&
            fname.substr(0, prefix.size()) == prefix &&
            fname.substr(fname.size() - suffix.size()) == suffix) {
            match = dir + "/" + fname;
            break;
        }
    }
    closedir(d);
    return match;
}

// Locate the per-combo sweep sub-folder, tolerating an optional run-name suffix
// appended to the folder (e.g. "_wb", "_wb_bayes", "_wbval" from --disc-suffix).
// The rnamotifs run names the folder "<base><suffix>" while this binary only
// knows <base>; prefer an exact match, otherwise accept "<base>_*" / "<base>-*".
static string find_combo_folder(const string &sweep_dir, const string &base) {
    // Exact match wins (canonical, no-suffix case).
    {
        ifstream t(sweep_dir + base + "/enriched_tetramers.txt");
        if (t.is_open()) return sweep_dir + base;
    }
    DIR *d = opendir(sweep_dir.c_str());
    if (!d) return sweep_dir + base;
    struct dirent *ent;
    string match;
    while ((ent = readdir(d)) != nullptr) {
        string fn = ent->d_name;
        // "<base>" followed by a separator ('_' or '-') then the suffix.
        if (fn.size() > base.size() &&
            fn.compare(0, base.size(), base) == 0 &&
            (fn[base.size()] == '_' || fn[base.size()] == '-')) {
            match = sweep_dir + fn;
            break;
        }
    }
    closedir(d);
    return match.empty() ? (sweep_dir + base) : match;
}

// ═══════════════════════════════════════════════════════════════════════════
// Read bootstrap TSV and write summary_results TSV
// ═══════════════════════════════════════════════════════════════════════════

struct BootstrapRow {
    string tetramer;
    double r1enh_pFis, r2enh_pFis, r3enh_pFis;
    double r1sil_pFis, r2sil_pFis, r3sil_pFis;
};

static vector<BootstrapRow> read_bootstrap_tsv(const string &path) {
    vector<BootstrapRow> rows;
    ifstream f(path);
    if (!f.is_open()) return rows;
    string line;
    getline(f, line); // header

    // Parse header to find column indices
    istringstream hss(line);
    vector<string> cols;
    string tok;
    while (hss >> tok) cols.push_back(tok);

    auto find_col = [&](const string &name) -> int {
        for (int i = 0; i < (int)cols.size(); ++i)
            if (cols[i] == name) return i;
        return -1;
    };

    int ci_tet = find_col("tetramer");
    int ci_r1e = find_col("r1enh_pFis");
    int ci_r2e = find_col("r2enh_pFis");
    int ci_r3e = find_col("r3enh_pFis");
    int ci_r1s = find_col("r1sil_pFis");
    int ci_r2s = find_col("r2sil_pFis");
    int ci_r3s = find_col("r3sil_pFis");

    while (getline(f, line)) {
        istringstream ss(line);
        vector<string> fields;
        string field;
        while (ss >> field) fields.push_back(field);
        if ((int)fields.size() <= max({ci_tet, ci_r1e, ci_r2e, ci_r3e, ci_r1s, ci_r2s, ci_r3s}))
            continue;

        BootstrapRow br;
        br.tetramer = (ci_tet >= 0) ? fields[ci_tet] : "";
        auto safe_d = [&](int ci) { return ci >= 0 ? stod(fields[ci]) : 1.0; };
        br.r1enh_pFis = safe_d(ci_r1e);
        br.r2enh_pFis = safe_d(ci_r2e);
        br.r3enh_pFis = safe_d(ci_r3e);
        br.r1sil_pFis = safe_d(ci_r1s);
        br.r2sil_pFis = safe_d(ci_r2s);
        br.r3sil_pFis = safe_d(ci_r3s);
        rows.push_back(br);
    }
    return rows;
}

static void write_summary_results(const string &output_dir, const string &name,
                                   int hw, int ew,
                                   const vector<BootstrapRow> &bootstrap_rows,
                                   const map<string, RegionInfo> &regions) {
    string params_key = "hw_" + to_string(hw) + "_ew_" + to_string(ew);
    string path = output_dir + "summary_results_" + name + "_" + params_key + ".tsv";
    ofstream f(path);
    if (!f.is_open()) { cerr << "Error: cannot write " << path << endl; return; }

    f << "pr\thw\tew\ttetramer\tR1\tR2\tR3\tr1enh_pFis\tr2enh_pFis\tr3enh_pFis\tr1sil_pFis\tr2sil_pFis\tr3sil_pFis\n";
    for (const auto &br : bootstrap_rows) {
        auto ri_it = regions.find(br.tetramer);
        string r1 = "0", r2 = "0", r3 = "0";
        if (ri_it != regions.end()) {
            r1 = ri_it->second.R1;
            r2 = ri_it->second.R2;
            r3 = ri_it->second.R3;
        }
        f << name << "\t" << hw << "\t" << ew << "\t" << br.tetramer
          << "\t" << r1 << "\t" << r2 << "\t" << r3
          << "\t" << br.r1enh_pFis << "\t" << br.r2enh_pFis << "\t" << br.r3enh_pFis
          << "\t" << br.r1sil_pFis << "\t" << br.r2sil_pFis << "\t" << br.r3sil_pFis
          << "\n";
    }
    cerr << "[mars_score] Wrote summary: " << path << endl;
}

// ═══════════════════════════════════════════════════════════════════════════
// Compute-peak mode: compute binding_profile_PEAK_normalized.tsv
// ═══════════════════════════════════════════════════════════════════════════

static int run_compute_peak(const string &mars_exons_dir, const string &eclip_dir,
                            const string &cell_line, const string &output_dir,
                            int in_exon, int in_intron, int cores,
                            const vector<string> &rbps) {
    cerr << "\n[mars_score] COMPUTE-PEAK mode" << endl;
    cerr << "  mars_exons_dir: " << mars_exons_dir << endl;
    cerr << "  Processing " << rbps.size() << " RBPs..." << endl;

    int wd = 15;
    int r = wd / 2;
    int in_intron_adj = in_intron + r;
    int in_exon_adj = in_exon + r;
    int adj_region_len = in_exon_adj + in_intron_adj + 1;
    int adj_n_pos_total = 4 * adj_region_len;

    // Region boundaries within adj_n_pos_total:
    // R1: [0, adj_region_len)
    // R2: [adj_region_len, 2*adj_region_len)
    // R3: [2*adj_region_len, 3*adj_region_len)

    // Output: rows=RBPs, cols=R1_enh,R2_enh,R3_enh,R1_sil,R2_sil,R3_sil
    vector<string> col_names = {"R1_enh","R2_enh","R3_enh","R1_sil","R2_sil","R3_sil"};
    map<string, vector<double>> peak_profiles; // rbp -> 6 values
    // Per-position CES profiles (full map V5|V6|V7|V8 in transcript order), Fig 2D/E
    map<string, vector<double>> ces_perpos_enh, ces_perpos_sil;   // smoothed obs CES (z)
    map<string, vector<double>> ces_qu_enh, ces_qu_sil;           // 95th percentile (qU)
    map<string, double> srr_enh, srr_sil;  // Signal Recovery Rate (downsampling AUC), per direction

    for (const auto &rbp : rbps) {
        cerr << "\n  " << rbp << ":" << endl;

        // Load exon file for this RBP
        string exon_file = mars_exons_dir + "/" + rbp + "_input_rnamotifs.txt";
        vector<Exon> exons = read_input_exons(exon_file);
        if (exons.empty()) {
            cerr << "    WARNING: no exons for " << rbp << ", skipping" << endl;
            peak_profiles[rbp] = {0,0,0,0,0,0};
            continue;
        }
        cerr << "    " << exons.size() << " exons" << endl;

        // Load eCLIP peaks
        string peak_file = find_eclip_file(eclip_dir, rbp);
        if (peak_file.empty()) {
            cerr << "    WARNING: no eCLIP file for " << rbp << endl;
            peak_profiles[rbp] = {0,0,0,0,0,0};
            continue;
        }
        vector<Peak> peaks = read_peaks(peak_file);
        cerr << "    " << peaks.size() << " peaks" << endl;

        // Compute overlap
        vector<int> overlap = compute_eclip_overlap(exons, peaks, in_exon_adj, in_intron_adj, cores);
        int n_exons = (int)exons.size();

        // Separate enh/sil/ctrl exon indices
        vector<int> enh_idx, sil_idx, ctrl_idx;
        for (int i = 0; i < n_exons; ++i) {
            if (exons[i].type == 1) enh_idx.push_back(i);
            else if (exons[i].type == -1) sil_idx.push_back(i);
            else ctrl_idx.push_back(i);
        }

        // ── Binding score (BS) per region, per Methods (¶78) ────────────────
        // CES = -2 log(Fisher p) per position (target enh/sil vs constitutive).
        // Empirical background: 1,000 label-shuffles -> per-position mean + 95th pct.
        // A position is "significantly bound" if observed CES > 95th percentile.
        // ΔCES = obs - mean(empirical); BS_region = max ΔCES over significant
        // positions in the region; BS scaled to per-RBP max afterwards.
        // Regions (¶77), strand handled by the row-reversal in compute_eclip_overlap
        // (positions are in transcript order): block 1 = 3'ss side (V6), block 2 =
        // 5'ss side (V7).  R1 = upstream intron (in_intron nt from 3'ss); R2 = exon
        // (in_exon nt each side); R3 = downstream intron (in_intron nt from 5'ss).
        // The 5 nt nearest each splice site are excluded from every region.
        // (Per-exon clamping to intron/exon midpoints is added in a follow-up step.)
        vector<double> profile(6, 0.0);
        const int n_bootstrap_peak = 1000;
        const int EXCL = 5;

        // 15-nt centred scrolling-window average of a CES profile (Methods ¶75,
        // R stats::filter with rep(1/15,15)); divisor fixed at 15, edges (within
        // h of the block border) left at 0 — region ranges start at q=r=h so are
        // unaffected.  Per-exon short-feature handling is region-extent only
        // (the fixed R1/R2/R3 position ranges), NOT applied to the Fisher test.
        const int SMW = 15, SMH = SMW / 2;
        auto smooth15 = [&](vector<double> &v) {
            vector<double> s(v.size(), 0.0);
            for (int p = SMH; p + SMH < (int)v.size(); ++p) {
                double sum = 0;
                for (int k = p - SMH; k <= p + SMH; ++k) sum += v[k];
                s[p] = sum / SMW;
            }
            v.swap(s);
        };

        // Per-position smoothed obs CES, empirical mean, and 95th percentile for
        // one block (CES exactly as in the paper: all exons, fixed offsets).
        auto block_stats = [&](const vector<int> &target_idx, int block_idx,
                               vector<double> &obs, vector<double> &emp_mean,
                               vector<double> &emp_p95) {
            int n_target = (int)target_idx.size();
            int n_ctrl   = (int)ctrl_idx.size();
            int base = block_idx * adj_region_len;
            obs.assign(adj_region_len, 0.0);
            emp_mean.assign(adj_region_len, 0.0);
            emp_p95.assign(adj_region_len, 0.0);
            if (n_target == 0 || n_ctrl == 0) return;

            for (int p = 0; p < adj_region_len; ++p) {
                int pp = base + p;
                int t = 0, c = 0;
                for (int ei : target_idx) t += overlap[ei * adj_n_pos_total + pp];
                for (int ei : ctrl_idx)   c += overlap[ei * adj_n_pos_total + pp];
                double pv = fisher_greater(t, c, n_target - t, n_ctrl - c);
                obs[p] = -2.0 * ((pv > 0) ? log(pv) : -700.0);
            }
            smooth15(obs);  // 15-nt scrolling-window average (Methods ¶75)

            vector<int> idx(target_idx);
            idx.insert(idx.end(), ctrl_idx.begin(), ctrl_idx.end());
            int N = (int)idx.size();
            vector<vector<double>> boot(n_bootstrap_peak,
                                        vector<double>(adj_region_len, 0.0));
            mt19937 rng(1234u);
            for (int b = 0; b < n_bootstrap_peak; ++b) {
                for (int i = N - 1; i > 0; --i) { int j = (int)(rng() % (unsigned)(i + 1)); swap(idx[i], idx[j]); }
                for (int p = 0; p < adj_region_len; ++p) {
                    int pp = base + p;
                    int t = 0, c = 0;
                    for (int k = 0; k < n_target; ++k) t += overlap[idx[k] * adj_n_pos_total + pp];
                    for (int k = n_target; k < N; ++k)  c += overlap[idx[k] * adj_n_pos_total + pp];
                    double pv = fisher_greater(t, c, n_target - t, N - n_target - c);
                    boot[b][p] = -2.0 * ((pv > 0) ? log(pv) : -700.0);
                }
                smooth15(boot[b]);
            }
            vector<double> col(n_bootstrap_peak);
            for (int p = 0; p < adj_region_len; ++p) {
                double s = 0;
                for (int b = 0; b < n_bootstrap_peak; ++b) { s += boot[b][p]; col[b] = boot[b][p]; }
                emp_mean[p] = s / n_bootstrap_peak;
                int k95 = (int)(0.95 * (n_bootstrap_peak - 1));
                nth_element(col.begin(), col.begin() + k95, col.end());
                emp_p95[p] = col[k95];
            }
        };

        // BS = max ΔCES over [q_lo,q_hi), where ΔCES = obs - 95th-percentile(qU).
        // (Validated 15/15 against the paper's binding_profile derived from its
        //  own Nobby z/qU.) Positions with obs <= qU contribute negative ΔCES and
        //  are naturally excluded by the max (floored at 0).
        auto region_bs = [&](const vector<double> &obs, const vector<double> &p95,
                             int q_lo, int q_hi) {
            double m = 0.0;
            for (int q = max(0, q_lo); q < min(adj_region_len, q_hi); ++q)
                m = max(m, obs[q] - p95[q]);
            return m;
        };

        // Block-local index ranges (q) for each region, with 5 nt ss exclusion.
        // Block 1 (V6): offset from 3'ss = q - in_intron_adj.
        int r1_lo = r,                          r1_hi = in_intron_adj - EXCL;        // [-in_intron, -(EXCL+1)]
        int r2a_lo = in_intron_adj + EXCL,      r2a_hi = in_intron_adj + in_exon + 1;// [+EXCL, +in_exon]
        // Block 2 (V7): offset from 5'ss = q - in_exon_adj.
        int r2b_lo = r,                         r2b_hi = in_exon_adj - EXCL;         // [-in_exon, -(EXCL+1)]
        int r3_lo = in_exon_adj + EXCL,         r3_hi = in_exon_adj + in_intron + 1; // [+EXCL, +in_intron]

        for (int d = 0; d < 2; ++d) {
            const vector<int> &tgt = (d == 0) ? enh_idx : sil_idx;
            // Full RNA splicing map = 4 blocks (V5 upstream constitutive exon,
            // V6 alt 3'ss, V7 alt 5'ss, V8 downstream constitutive exon). The BS
            // is measured ONLY on R1/R2/R3 around the cassette exon (blocks V6/V7);
            // blocks V5/V8 are computed solely for the full-map CES dump (Fig 2D).
            vector<double> o0, m0v, p0, o1, m1v, p1, o2, m2v, p2, o3, m3v, p3;
            block_stats(tgt, 0, o0, m0v, p0);   // V5 / upstream constitutive exon (map only)
            block_stats(tgt, 1, o1, m1v, p1);   // V6 / alt 3'ss  (BS R1,R2)
            block_stats(tgt, 2, o2, m2v, p2);   // V7 / alt 5'ss  (BS R2,R3)
            block_stats(tgt, 3, o3, m3v, p3);   // V8 / downstream constitutive exon (map only)
            double R1 = region_bs(o1, p1, r1_lo, r1_hi);
            double R2 = max(region_bs(o1, p1, r2a_lo, r2a_hi),
                            region_bs(o2, p2, r2b_lo, r2b_hi));
            double R3 = region_bs(o2, p2, r3_lo, r3_hi);
            profile[d * 3 + 0] = R1;
            profile[d * 3 + 1] = R2;
            profile[d * 3 + 2] = R3;
            // per-position dump (Fig 2D/E): smoothed obs CES (z) and 95th-pct null
            // qU, across the full map in transcript order V5|V6|V7|V8.
            vector<double> &dst  = (d == 0) ? ces_perpos_enh[rbp] : ces_perpos_sil[rbp];
            vector<double> &dstq = (d == 0) ? ces_qu_enh[rbp]     : ces_qu_sil[rbp];
            for (int q = 0; q < adj_region_len; ++q) { dst.push_back(o0[q]); dstq.push_back(p0[q]); }
            for (int q = 0; q < adj_region_len; ++q) { dst.push_back(o1[q]); dstq.push_back(p1[q]); }
            for (int q = 0; q < adj_region_len; ++q) { dst.push_back(o2[q]); dstq.push_back(p2[q]); }
            for (int q = 0; q < adj_region_len; ++q) { dst.push_back(o3[q]); dstq.push_back(p3[q]); }
        }

        peak_profiles[rbp] = profile;
        cerr << "    PEAK profile: ";
        for (double v : profile) cerr << fixed << setprecision(2) << v << " ";
        cerr << endl;

        // ── Signal Recovery Rate (SRR) per direction (Methods ¶75) ──────────
        // Downsample target exons to 10..100% (500 iters), compute the CES
        // profile over the full map, cosine to the full-set profile; SRR = AUC
        // of the (fraction, mean cosine) curve. ctrl hits are fixed across
        // subsets, so precompute them once.
        {
            int nc = (int)ctrl_idx.size();
            vector<int> ctrl_hits(adj_n_pos_total, 0);
            for (int pp = 0; pp < adj_n_pos_total; ++pp) {
                int c = 0;
                for (int ei : ctrl_idx) c += overlap[ei * adj_n_pos_total + pp];
                ctrl_hits[pp] = c;
            }
            auto profile_ces = [&](const vector<int> &tgt) {
                vector<double> ces(adj_n_pos_total, 0.0);
                int nt = (int)tgt.size();
                if (nt == 0 || nc == 0) return ces;
                for (int pp = 0; pp < adj_n_pos_total; ++pp) {
                    int t = 0;
                    for (int ei : tgt) t += overlap[ei * adj_n_pos_total + pp];
                    double pv = fisher_greater(t, ctrl_hits[pp], nt - t, nc - ctrl_hits[pp]);
                    ces[pp] = -2.0 * ((pv > 0) ? log(pv) : -700.0);
                }
                // SRR cosine uses the raw -2log(Fisher) profile (NOT smoothed) — the
                // paper's downsampling (04) operates on pFis directly.
                return ces;
            };
            for (int d = 0; d < 2; ++d) {
                const vector<int> &tgt = (d == 0) ? enh_idx : sil_idx;
                double &srr = (d == 0) ? srr_enh[rbp] : srr_sil[rbp];
                srr = 0.0;
                if ((int)tgt.size() < 2) continue;
                vector<double> full = profile_ces(tgt);
                mt19937 rng(7u);
                double x0 = 0.0, y0 = 0.0, auc = 0.0;
                for (int pi = 1; pi <= 10; ++pi) {
                    double frac = pi / 10.0;
                    int k = max(1, (int)(frac * (double)tgt.size()));
                    int nit = (pi == 10) ? 1 : 500;       // 100% is deterministic (cos=1)
                    double cossum = 0.0;
                    for (int it = 0; it < nit; ++it) {
                        vector<int> sub = tgt;
                        for (int i = 0; i < k; ++i) {       // partial Fisher-Yates, take first k
                            int j = i + (int)(rng() % (unsigned)((int)sub.size() - i));
                            swap(sub[i], sub[j]);
                        }
                        sub.resize(k);
                        cossum += cosine_similarity(profile_ces(sub), full);
                    }
                    double y = cossum / nit;
                    auc += (y0 + y) * (frac - x0) / 2.0;   // trapezoid from (x0,y0); x in [0,1]
                    x0 = frac; y0 = y;
                }
                srr = auc;
            }
            cerr << "    SRR: enh=" << fixed << setprecision(3) << srr_enh[rbp]
                 << " sil=" << srr_sil[rbp] << endl;
        }
    }

    // Write ABSOLUTE (un-normalised) BS profile — the per-region max ΔCES over
    // significant positions, BEFORE the per-RBP row scaling below. Used by the
    // supplementary BS~MS scatter (absolute binding score).
    {
        string abs_path = output_dir + cell_line + "_binding_profile_PEAK_absolute.tsv";
        ofstream af(abs_path);
        if (af.is_open()) {
            af << "rbp";
            for (const auto &c : col_names) af << "\t" << c;
            af << "\n";
            for (const auto &rbp : rbps) {
                auto it = peak_profiles.find(rbp);
                if (it == peak_profiles.end()) continue;
                af << rbp;
                for (double v : it->second) af << "\t" << fixed << setprecision(6) << v;
                af << "\n";
            }
            cerr << "[mars_score] Wrote ABSOLUTE PEAK profile: " << abs_path << endl;
        } else {
            cerr << "Warning: cannot write " << abs_path << endl;
        }
    }

    // Normalize per-RBP row to [0,1]
    for (auto &entry : peak_profiles) {
        double max_val = *max_element(entry.second.begin(), entry.second.end());
        if (max_val > 0) {
            for (double &v : entry.second) v /= max_val;
        }
    }

    // Write output
    string out_path = output_dir + cell_line + "_binding_profile_PEAK_normalized.tsv";
    ofstream f(out_path);
    if (!f.is_open()) { cerr << "Error: cannot write " << out_path << endl; return 1; }

    // Header
    f << "rbp";
    for (const auto &c : col_names) f << "\t" << c;
    f << "\n";

    // Data - write in sorted RBP order
    for (const auto &rbp : rbps) {
        auto it = peak_profiles.find(rbp);
        if (it == peak_profiles.end()) continue;
        f << rbp;
        for (double v : it->second) f << "\t" << fixed << setprecision(6) << v;
        f << "\n";
    }

    cerr << "[mars_score] Wrote PEAK profile: " << out_path << endl;

    // Signal Recovery Rate (SCORE1), per RBP and direction (Methods ¶75).
    {
        string srr_path = output_dir + cell_line + "_SRR.tsv";
        ofstream sf(srr_path);
        if (sf.is_open()) {
            sf << "rbp\tauc\n";
            for (const auto &rbp : rbps) {
                if (srr_enh.count(rbp)) sf << rbp << "_enh\t" << fixed << setprecision(6) << srr_enh[rbp] << "\n";
                if (srr_sil.count(rbp)) sf << rbp << "_sil\t" << fixed << setprecision(6) << srr_sil[rbp] << "\n";
            }
            cerr << "[mars_score] Wrote SRR (SCORE1): " << srr_path << endl;
        } else {
            cerr << "Warning: cannot write " << srr_path << endl;
        }
    }

    // Per-position CES dump (long format) for Figure 2 panels D/E.
    // pos runs 0..(4*adj_region_len-1), concatenating the full RNA splicing map
    // in transcript order: V5 (upstream exon) | V6 (alt 3'ss) | V7 (alt 5'ss) | V8
    // (downstream exon). adj_region_len positions per block.
    string ces_path = output_dir + cell_line + "_CES_perpos.tsv";
    ofstream cf(ces_path);
    if (cf.is_open()) {
        cf << "rbp\tdirection\tpos\tobs\tqU\n";
        for (const auto &rbp : rbps) {
            for (int d = 0; d < 2; ++d) {
                const char *dir = (d == 0) ? "enh" : "sil";
                auto &m  = (d == 0) ? ces_perpos_enh : ces_perpos_sil;
                auto &mq = (d == 0) ? ces_qu_enh     : ces_qu_sil;
                auto it = m.find(rbp);
                if (it == m.end()) continue;
                const auto &qv = mq[rbp];
                for (size_t p = 0; p < it->second.size(); ++p)
                    cf << rbp << "\t" << dir << "\t" << p << "\t"
                       << fixed << setprecision(6) << it->second[p] << "\t"
                       << (p < qv.size() ? qv[p] : 0.0) << "\n";
            }
        }
        cerr << "[mars_score] Wrote per-position CES: " << ces_path << endl;
    } else {
        cerr << "Warning: cannot write " << ces_path << endl;
    }
    return 0;
}

// ═══════════════════════════════════════════════════════════════════════════
// Main
// ═══════════════════════════════════════════════════════════════════════════

int main(int argc, char *argv[]) {
    // Parse arguments
    string input_file, name, sweep_dir, mars_ref_dir, cell_line;
    string eclip_dir, output_dir, mars_exons_dir;
    int cores = 1, bootstraps = 1000;
    int in_exon = 30, in_intron = 300;
    bool do_summary = false;
    bool do_compute_peak = false;
    string score_mode = "full";   // "full" = CS x SRR ; "cs-only" = CS (SRR disabled)

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--score-mode" && i+1 < argc) score_mode = argv[++i];
        else if (arg == "-i" && i+1 < argc) input_file = argv[++i];
        else if (arg == "-n" && i+1 < argc) name = argv[++i];
        else if (arg == "-d" && i+1 < argc) sweep_dir = argv[++i];
        else if (arg == "-r" && i+1 < argc) mars_ref_dir = argv[++i];
        else if (arg == "-c" && i+1 < argc) cell_line = argv[++i];
        else if (arg == "-p" && i+1 < argc) cores = stoi(argv[++i]);
        else if (arg == "-e" && i+1 < argc) eclip_dir = argv[++i];
        else if (arg == "-o" && i+1 < argc) output_dir = argv[++i];
        else if (arg == "-b" && i+1 < argc) bootstraps = stoi(argv[++i]);
        else if (arg == "--in-exon" && i+1 < argc) in_exon = stoi(argv[++i]);
        else if (arg == "--in-intron" && i+1 < argc) in_intron = stoi(argv[++i]);
        else if (arg == "--summary") do_summary = true;
        else if (arg == "--compute-peak") do_compute_peak = true;
        else if (arg == "--mars-exons-dir" && i+1 < argc) mars_exons_dir = argv[++i];
    }

#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif

    // ── Define RBPs ─────────────────────────────────────────────────────
    vector<string> rbps;
    if (cell_line == "HepG2")
        rbps = {"HNRNPC","HNRNPK","HNRNPU","NCBP2","PRPF8","PTBP1","QKI",
                "RBFOX2","RBM22","SF3A3","SF3B4","SRSF1","U2AF1","U2AF2","UCHL5"};
    else if (cell_line == "K562")
        rbps = {"AGGF1","EFTUD2","FXR1","HNRNPU","PRPF8","PTBP1","PUS1",
                "RBM15","SF3B4","SRSF1","TARDBP","U2AF1","U2AF2"};
    else {
        cerr << "Error: unsupported cell line " << cell_line << endl;
        return 1;
    }

    // ── Compute-peak mode ─────────────────────────────────────────────
    if (do_compute_peak) {
        if (mars_exons_dir.empty() || eclip_dir.empty() || output_dir.empty()) {
            cerr << "Usage: " << argv[0]
                 << " --compute-peak --mars-exons-dir <dir> -e <eclip_dir>"
                 << " -o <output_dir> -c <cell_line> [-p cores]"
                 << " [--in-exon N] [--in-intron N]" << endl;
            return 1;
        }
        if (eclip_dir.back() != '/') eclip_dir += '/';
        if (output_dir.back() != '/') output_dir += '/';
        if (mars_exons_dir.back() != '/') mars_exons_dir += '/';
        return run_compute_peak(mars_exons_dir, eclip_dir, cell_line,
                                output_dir, in_exon, in_intron, cores, rbps);
    }

    // ── Standard scoring mode ─────────────────────────────────────────
    if (input_file.empty() || name.empty() || sweep_dir.empty() ||
        cell_line.empty() || eclip_dir.empty() || output_dir.empty()) {
        cerr << "Usage: " << argv[0]
             << " -i <input> -n <name> -d <sweep_dir> -r <mars_ref_dir>"
             << " -c <cell_line> -p <cores> -e <eclip_dir> -o <output_dir>"
             << " [--in-exon N] [--in-intron N] [-b N] [--summary]"
             << " [--score-mode full|cs-only]" << endl;
        return 1;
    }
    if (score_mode != "full" && score_mode != "cs-only") {
        cerr << "Error: --score-mode must be 'full' (AS = CS x SRR, default) or "
             << "'cs-only' (AS = CS); got '" << score_mode << "'" << endl;
        return 1;
    }
    cerr << "[mars_score] score-mode: " << score_mode
         << (score_mode == "cs-only" ? " (SRR weight disabled, AS = CS)"
                                      : " (AS = CS x SRR)") << endl;

    // Ensure trailing slashes
    if (sweep_dir.back() != '/') sweep_dir += '/';
    if (eclip_dir.back() != '/') eclip_dir += '/';
    if (output_dir.back() != '/') output_dir += '/';

    cerr << "====================================================\n"
         << "  RNAmotifs-MaRs Score Computation (C++)\n"
         << "  Input:      " << input_file << "\n"
         << "  Name:       " << name << "\n"
         << "  Sweep dir:  " << sweep_dir << "\n"
         << "  Cell line:  " << cell_line << "\n"
         << "  Cores:      " << cores << "\n"
         << "  in_exon:    " << in_exon << "\n"
         << "  in_intron:  " << in_intron << "\n"
         << "  Summary:    " << (do_summary ? "yes" : "no") << "\n"
         << "====================================================" << endl;

    // ── Read optimal parameters ─────────────────────────────────────────
    string opt_path = mars_ref_dir + "/Tables/RNAmotifs_optimal_parameters.csv";
    vector<OptimalParams> opt_params = read_optimal_params(opt_path, rbps);
    if (opt_params.empty()) {
        cerr << "Error: no optimal parameters loaded from " << opt_path << endl;
        return 1;
    }

    // Get unique (hw, ew) combos
    set<pair<int,int>> unique_combos;
    for (const auto &op : opt_params)
        unique_combos.insert({op.hw, op.ew});
    cerr << "[mars_score] " << opt_params.size() << " RBP entries, "
         << unique_combos.size() << " unique (hw,ew) combos" << endl;

    // Build RBP -> list of param combos mapping
    map<string, vector<pair<int,int>>> rbp_to_combos;
    for (const auto &op : opt_params)
        rbp_to_combos[op.rbp].push_back({op.hw, op.ew});

    // ── Read AUC scores (SCORE1) ────────────────────────────────────────
    string auc_path = mars_ref_dir + "/Rdata/" + cell_line + "_AUC.tsv";
    map<string, double> auc_map = read_auc_tsv(auc_path);
    cerr << "[mars_score] Loaded " << auc_map.size() << " AUC entries" << endl;

    // ── Read input exons ────────────────────────────────────────────────
    vector<Exon> exons = read_input_exons(input_file);
    cerr << "[mars_score] " << exons.size() << " input exons" << endl;

    // ── Position adjustment parameters ──────────────────────────────────
    int wd = 15;
    int r = wd / 2;
    int in_intron_adj = in_intron + r;
    int in_exon_adj = in_exon + r;
    int region_len = in_exon + in_intron + 1;
    int adj_region_len = in_exon_adj + in_intron_adj + 1;
    int adj_n_pos_total = 4 * adj_region_len;

    // pos_adj: adjusted position indices for profile matching
    vector<int> pos_adj;
    for (int k = 0; k < 4; ++k) {
        int base = k * (region_len + 2 * r);
        for (int p = 0; p < region_len; ++p)
            pos_adj.push_back(base + r + p);
    }

    // ── Pre-load per-combo data (tetramers, regions, score CSVs) ────────
    // These are small and shared across RBPs, so load once.
    struct ComboData {
        int hw, ew;
        string params_key;
        string folder_path;
        vector<string> enh_tets, sil_tets;
        map<string, RegionInfo> regions;
        ScoreMatrix sm_enh, sm_sil;
        vector<BootstrapRow> bootstrap_rows; // for --summary
        bool valid;
    };

    map<pair<int,int>, ComboData> combo_data;
    for (const auto &combo : unique_combos) {
        int hw = combo.first, ew = combo.second;
        ComboData cd;
        cd.hw = hw;
        cd.ew = ew;
        cd.params_key = "hw_" + to_string(hw) + "_ew_" + to_string(ew);
        cd.folder_path = find_combo_folder(
            sweep_dir, name + "_hw_" + to_string(hw) + "_ew_" + to_string(ew));
        cd.valid = false;

        // Check if folder exists
        ifstream test(cd.folder_path + "/enriched_tetramers.txt");
        if (!test.is_open()) {
            cerr << "  [" << cd.params_key << "] No results, skipping" << endl;
            combo_data[combo] = cd;
            continue;
        }

        // Read enriched tetramers and regions
        vector<string> enr_tets = read_enriched_tetramers(
            cd.folder_path + "/enriched_tetramers.txt");
        cd.regions = read_regions_csv(cd.folder_path + "/regions.csv");

        if (enr_tets.empty()) {
            combo_data[combo] = cd;
            continue;
        }

        // Classify tetramers as enh/sil
        for (const auto &tet : enr_tets) {
            auto it = cd.regions.find(tet);
            if (it == cd.regions.end()) continue;
            const auto &ri = it->second;
            bool has_enh = (ri.R1 == "Enhanced" || ri.R1 == "Both" ||
                           ri.R2 == "Enhanced" || ri.R2 == "Both" ||
                           ri.R3 == "Enhanced" || ri.R3 == "Both");
            bool has_sil = (ri.R1 == "Silenced" || ri.R1 == "Both" ||
                           ri.R2 == "Silenced" || ri.R2 == "Both" ||
                           ri.R3 == "Silenced" || ri.R3 == "Both");
            if (has_enh) cd.enh_tets.push_back(tet);
            if (has_sil) cd.sil_tets.push_back(tet);
        }

        cerr << "  [" << cd.params_key << "] " << enr_tets.size() << " enriched, "
             << cd.enh_tets.size() << " enh, " << cd.sil_tets.size() << " sil" << endl;

        // Read splicing map score CSVs (filenames contain timestamps, find by prefix)
        string enh_csv = find_file_by_prefix(cd.folder_path, "e-", ".csv");
        string sil_csv = find_file_by_prefix(cd.folder_path, "s-", ".csv");
        cd.sm_enh = read_score_csv(enh_csv);
        cd.sm_sil = read_score_csv(sil_csv);

        // Read bootstrap data for --summary
        if (do_summary) {
            string bs_path = cd.folder_path + "/bootstrap_" + to_string(bootstraps) + ".tsv";
            // Try the Rdata-converted TSV first, fall back to raw
            ifstream bs_test(bs_path);
            if (!bs_test.is_open()) {
                // Try reading from the original rnamotifs result
                bs_path = cd.folder_path + "/bootstrap_" + to_string(bootstraps) + ".tsv";
            }
            cd.bootstrap_rows = read_bootstrap_tsv(bs_path);
        }

        cd.valid = true;
        combo_data[combo] = cd;
    }

    // ── Score matrices ──────────────────────────────────────────────────
    // SCORE1[params_key][rbp][tet] = AUC weight
    // SCORE2[params_key][rbp][tet] = cosine similarity
    map<string, map<string, map<string, double>>> score1_enh, score1_sil;
    map<string, map<string, map<string, double>>> score2_enh, score2_sil;

    // ── Process RBPs one at a time (memory-safe) ────────────────────────
    cerr << "\n[mars_score] Processing RBPs sequentially (memory-safe)..." << endl;

    int n_exons = (int)exons.size();

    for (const auto &rbp : rbps) {
        auto combo_it = rbp_to_combos.find(rbp);
        if (combo_it == rbp_to_combos.end()) continue;

        cerr << "\n  RBP: " << rbp << endl;

        // Load eCLIP peaks for this RBP
        string peak_file = find_eclip_file(eclip_dir, rbp);
        if (peak_file.empty()) {
            cerr << "    WARNING: no eCLIP file for " << rbp << endl;
            continue;
        }
        cerr << "    Loading peaks: " << peak_file << endl;
        vector<Peak> peaks = read_peaks(peak_file);
        cerr << "    " << peaks.size() << " peaks, computing overlap..." << endl;

        // Compute overlap for this RBP
        vector<int> overlap = compute_eclip_overlap(
            exons, peaks, in_exon_adj, in_intron_adj, cores);

        // Process all param combos for this RBP
        for (const auto &combo : combo_it->second) {
            auto cd_it = combo_data.find(combo);
            if (cd_it == combo_data.end() || !cd_it->second.valid) continue;
            const ComboData &cd = cd_it->second;

            // SCORE1: AUC-based weight
            for (const auto &tet : cd.enh_tets) {
                double auc = auc_map.count(rbp + "_enh") ? auc_map[rbp + "_enh"] : 0.0;
                score1_enh[cd.params_key][rbp][tet] = auc;
            }
            for (const auto &tet : cd.sil_tets) {
                double auc = auc_map.count(rbp + "_sil") ? auc_map[rbp + "_sil"] : 0.0;
                score1_sil[cd.params_key][rbp][tet] = auc;
            }

            // SCORE2: Cosine similarity
            auto compute_score2 = [&](const vector<string> &tets,
                                      const ScoreMatrix &sm,
                                      const string &direction,
                                      map<string, map<string, map<string, double>>> &s2) {
                int target_type = (direction == "enh") ? 1 : -1;

                for (const auto &tet : tets) {
                    auto reg_it = cd.regions.find(tet);
                    if (reg_it == cd.regions.end()) continue;
                    const auto &ri = reg_it->second;

                    set<int> exon_ids;
                    for (int rnum = 1; rnum <= 3; ++rnum) {
                        string rval;
                        if (rnum == 1) rval = ri.R1;
                        else if (rnum == 2) rval = ri.R2;
                        else rval = ri.R3;
                        if (rval != "0") {
                            auto eid = get_exons_with_tetramer(
                                tet, cd.folder_path, rnum, exons);
                            exon_ids.insert(eid.begin(), eid.end());
                        }
                    }

                    if (exon_ids.empty()) {
                        s2[cd.params_key][rbp][tet] = 0.0;
                        continue;
                    }

                    // Compute binding profile via Fisher test
                    vector<double> bind_profile(adj_n_pos_total, 0.0);

                    vector<int> target_indices, ctrl_indices;
                    for (int ei = 0; ei < n_exons; ++ei) {
                        if (exon_ids.count(exons[ei].id) && exons[ei].type == target_type)
                            target_indices.push_back(ei);
                        if (exons[ei].type == 0)
                            ctrl_indices.push_back(ei);
                    }
                    int n_target = (int)target_indices.size();
                    int n_ctrl = (int)ctrl_indices.size();

                    if (n_target == 0 || n_ctrl == 0) {
                        s2[cd.params_key][rbp][tet] = 0.0;
                        continue;
                    }

                    // Sum eCLIP hits per position
                    vector<int> target_hits(adj_n_pos_total, 0);
                    vector<int> ctrl_hits(adj_n_pos_total, 0);
                    for (int ei : target_indices)
                        for (int p = 0; p < adj_n_pos_total; ++p)
                            target_hits[p] += overlap[ei * adj_n_pos_total + p];
                    for (int ei : ctrl_indices)
                        for (int p = 0; p < adj_n_pos_total; ++p)
                            ctrl_hits[p] += overlap[ei * adj_n_pos_total + p];

                    // Fisher test at each position
                    for (int p = 0; p < adj_n_pos_total; ++p) {
                        double pv = fisher_greater(target_hits[p], ctrl_hits[p],
                                                   n_target - target_hits[p],
                                                   n_ctrl - ctrl_hits[p]);
                        bind_profile[p] = -2.0 * ((pv > 0) ? log(pv) : -700.0);
                    }

                    // Extract adjusted positions and compute cosine similarity
                    vector<double> x1, x2;
                    for (int p : pos_adj)
                        if (p < (int)bind_profile.size())
                            x1.push_back(bind_profile[p]);

                    // Find tetramer column in score matrix
                    int tet_col = -1;
                    for (int ti = 0; ti < sm.n_tets; ++ti)
                        if (sm.tetramers[ti] == tet) { tet_col = ti; break; }

                    if (tet_col >= 0) {
                        for (int pi = 0; pi < sm.n_pos; ++pi)
                            x2.push_back(sm.data[pi * sm.n_tets + tet_col]);
                    }

                    // Match dimensions
                    int min_len = min(x1.size(), x2.size());
                    if (min_len > 0) {
                        x1.resize(min_len);
                        x2.resize(min_len);
                        s2[cd.params_key][rbp][tet] = cosine_similarity(x1, x2);
                    } else {
                        s2[cd.params_key][rbp][tet] = 0.0;
                    }
                }
            };

            compute_score2(cd.enh_tets, cd.sm_enh, "enh", score2_enh);
            compute_score2(cd.sil_tets, cd.sm_sil, "sil", score2_sil);
        }

        cerr << "    " << rbp << ": scores computed, freeing overlap" << endl;
        // overlap vector goes out of scope here → memory freed
    }

    // ── Create diagnostics subdirectory ─────────────────────────────────
    string diag_dir = output_dir + "diagnostics/";
    {
        string cmd = "mkdir -p " + diag_dir;
        system(cmd.c_str());
    }

    // ── Write summary results if requested ─────────────────────────────
    if (do_summary) {
        cerr << "\n[mars_score] Writing summary results..." << endl;
        for (const auto &cd_pair : combo_data) {
            const ComboData &cd = cd_pair.second;
            if (!cd.valid || cd.bootstrap_rows.empty()) continue;
            write_summary_results(diag_dir, name, cd.hw, cd.ew,
                                  cd.bootstrap_rows, cd.regions);
        }
    }

    // ── Write output scores ─────────────────────────────────────────────
    cerr << "\n[mars_score] Writing output files..." << endl;

    // Flatten scores: aggregate across parameter combos per RBP
    auto flatten_scores = [&](
            const map<string, map<string, map<string, double>>> &s1,
            const map<string, map<string, map<string, double>>> &s2,
            const string &direction) {

        map<string, map<string, double>> final_scores; // rbp -> tet -> score

        for (const auto &op : opt_params) {
            string pk = op.params_key;
            auto s1_it = s1.find(pk);
            auto s2_it = s2.find(pk);
            if (s1_it == s1.end() || s2_it == s2.end()) continue;

            auto rbp_s1 = s1_it->second.find(op.rbp);
            auto rbp_s2 = s2_it->second.find(op.rbp);
            if (rbp_s1 == s1_it->second.end() || rbp_s2 == s2_it->second.end())
                continue;

            for (const auto &tet_s1 : rbp_s1->second) {
                double sc1 = tet_s1.second;          // SCORE1 = SRR
                auto tet_s2 = rbp_s2->second.find(tet_s1.first);
                double sc2 = (tet_s2 != rbp_s2->second.end()) ? tet_s2->second : 0.0;  // SCORE2 = CS
                // AS = CS x SRR (full) or CS alone (cs-only, SRR disabled / SRR=1)
                final_scores[op.rbp][tet_s1.first] =
                    (score_mode == "cs-only") ? sc2 : sc1 * sc2;
            }
        }
        return final_scores;
    };

    auto final_enh = flatten_scores(score1_enh, score2_enh, "enh");
    auto final_sil = flatten_scores(score1_sil, score2_sil, "sil");

    // Scale AS to [0,1] range per direction (paper: Figure S3)
    auto scale_to_01 = [](map<string, map<string, double>> &scores) {
        double max_val = 0;
        for (const auto &rbp_entry : scores)
            for (const auto &tet_entry : rbp_entry.second)
                if (!isnan(tet_entry.second) && tet_entry.second > max_val)
                    max_val = tet_entry.second;
        if (max_val > 0)
            for (auto &rbp_entry : scores)
                for (auto &tet_entry : rbp_entry.second)
                    if (!isnan(tet_entry.second))
                        tet_entry.second /= max_val;
    };
    scale_to_01(final_enh);
    scale_to_01(final_sil);

    write_score_tsv(output_dir + "association_scores_enh_" + name + ".tsv",
                    final_enh, "enhanced association scores");
    write_score_tsv(output_dir + "association_scores_sil_" + name + ".tsv",
                    final_sil, "silenced association scores");

    // Also write SCORE1 and SCORE2 separately into diagnostics/
    auto write_nested = [&](const string &prefix,
                            const map<string, map<string, map<string, double>>> &scores) {
        for (const auto &pk_entry : scores) {
            write_score_tsv(diag_dir + prefix + "_" + pk_entry.first + "_" + name + ".tsv",
                           pk_entry.second, prefix + " " + pk_entry.first);
        }
    };

    write_nested("SCORE1_enh", score1_enh);
    write_nested("SCORE1_sil", score1_sil);
    write_nested("SCORE2_enh", score2_enh);
    write_nested("SCORE2_sil", score2_sil);

    // ── RBP ranking ─────────────────────────────────────────────────────
    cerr << "\n[mars_score] RBP ranking:" << endl;
    vector<pair<double, string>> ranking;
    for (const auto &rbp_entry : final_sil) {
        double mean_score = 0;
        int count = 0;
        for (const auto &tet_entry : rbp_entry.second) {
            if (!isnan(tet_entry.second)) {
                mean_score += tet_entry.second;
                count++;
            }
        }
        if (count > 0) mean_score /= count;
        ranking.push_back({mean_score, rbp_entry.first});
    }
    // Add enhanced scores
    for (const auto &rbp_entry : final_enh) {
        bool found = false;
        for (auto &r : ranking) {
            if (r.second == rbp_entry.first) {
                double enh_mean = 0;
                int count = 0;
                for (const auto &te : rbp_entry.second) {
                    if (!isnan(te.second)) { enh_mean += te.second; count++; }
                }
                if (count > 0) r.first = max(r.first, enh_mean / count);
                found = true;
                break;
            }
        }
        if (!found) {
            double mean_score = 0;
            int count = 0;
            for (const auto &te : rbp_entry.second) {
                if (!isnan(te.second)) { mean_score += te.second; count++; }
            }
            if (count > 0) mean_score /= count;
            ranking.push_back({mean_score, rbp_entry.first});
        }
    }

    sort(ranking.rbegin(), ranking.rend());
    for (const auto &r : ranking)
        cerr << "  " << r.second << "\t" << fixed << setprecision(4) << r.first << endl;

    // Write ranking
    {
        string path = output_dir + "rbp_ranking_" + name + ".tsv";
        ofstream out(path);
        out << "rbp\tmean_association_score\n";
        for (const auto &r : ranking)
            out << r.second << "\t" << fixed << setprecision(6) << r.first << "\n";
        cerr << "[mars_score] Wrote " << path << endl;
    }

    cerr << "\n[mars_score] Done." << endl;
    return 0;
}
