// RNAmotifs - Tetramer selection and RNA splicing map data (C++)
// Copyright (C) 2014-2026 Matteo Cereda
// SPDX-License-Identifier: GPL-2.0-or-later
//
// Replaces selection.R with a fast C++ implementation.
// Reads bootstrap TSV, applies significance thresholds, computes
// positional Fisher tests, clusters tetramers, and writes all outputs.
//
// Usage:
//   rnamotifs_selection <results_dir> <name> <n_bootstraps>
//       <p_fisher> <p_empirical> [top_n] [counts_dir]
//       [in_exon] [in_intron] [--mars] [--mars-groups]

#include <algorithm>
#include <cmath>
#include <cstdlib>
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

using namespace std;

// ═══════════════════════════════════════════════════════════════════════════
// IUPAC code expansion
// ═══════════════════════════════════════════════════════════════════════════

static const map<char, string> IUPAC_EXPAND = {
    {'A', "A"}, {'C', "C"}, {'G', "G"}, {'T', "T"},
    {'R', "(A/G)"}, {'Y', "(C/T)"}, {'S', "(G/C)"}, {'W', "(A/T)"},
    {'K', "(G/T)"}, {'M', "(A/C)"},
    {'B', "(C/G/T)"}, {'D', "(A/G/T)"}, {'H', "(A/C/T)"}, {'V', "(A/C/G)"}
};

static const map<char, vector<char>> IUPAC_BASES = {
    {'A', {'A'}}, {'C', {'C'}}, {'G', {'G'}}, {'T', {'T'}},
    {'R', {'A','G'}}, {'Y', {'C','T'}}, {'S', {'G','C'}}, {'W', {'A','T'}},
    {'K', {'G','T'}}, {'M', {'A','C'}},
    {'B', {'C','G','T'}}, {'D', {'A','G','T'}}, {'H', {'A','C','T'}}, {'V', {'A','C','G'}}
};

static string iupac_full(const string &tet) {
    string out;
    for (char c : tet) {
        auto it = IUPAC_EXPAND.find(c);
        if (it != IUPAC_EXPAND.end()) out += it->second;
        else out += c;
    }
    return out;
}

// ═══════════════════════════════════════════════════════════════════════════
// Fisher exact test (one-sided, alternative = "greater")
// ═══════════════════════════════════════════════════════════════════════════

static double fisher_greater(int a, int b, int c, int d) {
    int n    = a + b + c + d;
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
// Region spacing computation (matches C++ core and R selection.R)
// ═══════════════════════════════════════════════════════════════════════════

static int compute_spacing(int in_exon, int in_intron) {
    return 2 * max(in_exon, in_intron) + 10;
}

static vector<int> compute_regions(int spacing) {
    vector<int> r(4);
    for (int i = 0; i < 4; ++i) r[i] = spacing * (i + 1);
    return r;
}

// Build rown: the position indices for all 4 regions
// R1: regions[0]-in_exon .. regions[0]+in_intron
// R2: regions[1]-in_intron .. regions[1]+in_exon
// R3: regions[2]-in_exon .. regions[2]+in_intron
// R4: regions[3]-in_intron .. regions[3]+in_exon
static vector<int> compute_rown(int in_exon, int in_intron, const vector<int> &regions) {
    vector<int> rown;
    auto add_range = [&](int start, int end) {
        for (int i = start; i <= end; ++i) rown.push_back(i);
    };
    add_range(regions[0] - in_exon,   regions[0] + in_intron);
    add_range(regions[1] - in_intron, regions[1] + in_exon);
    add_range(regions[2] - in_exon,   regions[2] + in_intron);
    add_range(regions[3] - in_intron, regions[3] + in_exon);
    return rown;
}

// ═══════════════════════════════════════════════════════════════════════════
// Bootstrap TSV reader
// ═══════════════════════════════════════════════════════════════════════════

struct BootstrapRow {
    string tetramer;
    // 12 values: r1enh_pFis, r1enh_pEmp, r1sil_pFis, r1sil_pEmp,
    //            r2enh_pFis, r2enh_pEmp, r2sil_pFis, r2sil_pEmp,
    //            r3enh_pFis, r3enh_pEmp, r3sil_pFis, r3sil_pEmp
    double r_enh_pFis[3], r_enh_pEmp[3];
    double r_sil_pFis[3], r_sil_pEmp[3];
};

static vector<BootstrapRow> read_bootstrap_tsv(const string &path) {
    vector<BootstrapRow> rows;
    ifstream f(path);
    if (!f.is_open()) {
        cerr << "Error: cannot open " << path << endl;
        return rows;
    }
    string line;
    getline(f, line); // skip header
    while (getline(f, line)) {
        istringstream ss(line);
        BootstrapRow br;
        ss >> br.tetramer;
        for (int r = 0; r < 3; ++r) {
            ss >> br.r_enh_pFis[r] >> br.r_enh_pEmp[r]
               >> br.r_sil_pFis[r] >> br.r_sil_pEmp[r];
        }
        rows.push_back(br);
    }
    return rows;
}

// ═══════════════════════════════════════════════════════════════════════════
// getTables: read positional tetramer counts from .bed files
// Returns 3 matrices (enh, sil, cont) as flat vectors [n_pos x n_tets]
// ═══════════════════════════════════════════════════════════════════════════

struct PositionalData {
    vector<double> enh;   // [n_pos * n_tets]
    vector<double> sil;
    vector<double> cont;
    int n_pos;
    int n_tets;
    vector<string> tet_names;
};

static vector<string> read_filelist(const string &dir) {
    vector<string> files;
    // Try filelist.txt first, then filelist_count.tsv
    for (const string &fname : {"filelist.txt", "filelist_count.tsv"}) {
        string path = dir + fname;
        ifstream f(path);
        if (!f.is_open()) continue;
        string line;
        while (getline(f, line)) {
            line.erase(line.find_last_not_of(" \t\r\n") + 1);
            if (!line.empty() && line.find("filelist") == string::npos)
                files.push_back(line);
        }
        break;
    }
    return files;
}

// Read a .bed file and return a map: position -> {cat -> count}
static void read_bed_file(const string &path,
                          unordered_map<int, unordered_map<int, int>> &data) {
    ifstream f(path);
    if (!f.is_open()) return;
    string line;
    // .bed files have no header, columns: pos cat count
    while (getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        istringstream ss(line);
        int pos, cat, count;
        if (ss >> pos >> cat >> count) {
            data[pos][cat] = count;
        }
    }
}

static PositionalData get_tables(const string &counts_dir,
                                 const vector<string> &tets,
                                 const vector<int> &rown) {
    PositionalData pd;
    pd.n_pos = (int)rown.size();

    // Collect .bed files from both r/ and nr/ subdirs
    vector<pair<string, string>> found; // (tet_name, filepath)
    for (const string &sub : {"nr/", "r/"}) {
        string subdir = counts_dir + sub;
        vector<string> flist = read_filelist(subdir);
        for (const string &fname : flist) {
            if (fname.size() < 4 || fname.substr(fname.size()-4) != ".bed")
                continue;
            string tet_name = fname.substr(0, fname.find('.'));
            // Only include if tetramer is in our requested list
            bool wanted = false;
            for (const auto &t : tets)
                if (tet_name.find(t) == 0 && tet_name.size() == t.size()) {
                    wanted = true; break;
                }
            if (wanted) found.push_back({tet_name, subdir + fname});
        }
    }

    // Deduplicate: prefer nr/ over r/
    map<string, string> tet_to_file;
    for (const auto &p : found)
        tet_to_file[p.first] = p.second; // last wins (r/ overwrites nr/)

    // Filter to requested order
    for (const auto &t : tets) {
        if (tet_to_file.count(t))
            pd.tet_names.push_back(t);
    }
    pd.n_tets = (int)pd.tet_names.size();
    if (pd.n_tets == 0) return pd;

    pd.enh.assign(pd.n_pos * pd.n_tets, 0.0);
    pd.sil.assign(pd.n_pos * pd.n_tets, 0.0);
    pd.cont.assign(pd.n_pos * pd.n_tets, 0.0);

    for (int ti = 0; ti < pd.n_tets; ++ti) {
        const string &tet = pd.tet_names[ti];
        string fpath = tet_to_file[tet];

        unordered_map<int, unordered_map<int, int>> bed_data;
        read_bed_file(fpath, bed_data);

        for (int pi = 0; pi < pd.n_pos; ++pi) {
            int pos = rown[pi];
            int idx = pi * pd.n_tets + ti;
            auto it = bed_data.find(pos);
            if (it != bed_data.end()) {
                auto &cats = it->second;
                if (cats.count(1))  pd.enh[idx]  = cats[1];
                if (cats.count(-1)) pd.sil[idx]  = cats[-1];
                if (cats.count(0))  pd.cont[idx] = cats[0];
            }
        }
    }
    return pd;
}

// ═══════════════════════════════════════════════════════════════════════════
// Positional Fisher test (lmb.cluster.fisher.test equivalent)
// Returns p-value matrix [n_pos x n_tets]
// ═══════════════════════════════════════════════════════════════════════════

static vector<double> positional_fisher(const vector<double> &d, int tot_d,
                                        const vector<double> &contr, int tot_contr,
                                        int n_pos, int n_tets) {
    vector<double> pvals(n_pos * n_tets);
    for (int ti = 0; ti < n_tets; ++ti) {
        for (int pi = 0; pi < n_pos; ++pi) {
            int idx = pi * n_tets + ti;
            int a = (int)d[idx];
            int c = (int)contr[idx];
            int b = tot_d - a;
            int dd = tot_contr - c;
            pvals[idx] = fisher_greater(a, c, b, dd);
        }
    }
    return pvals;
}

// ═══════════════════════════════════════════════════════════════════════════
// Trimeric overlap and clustering (SortingAndType)
// ═══════════════════════════════════════════════════════════════════════════

static vector<string> get_trimers(const string &mot) {
    vector<string> trimers;
    auto expand = [](char c) -> vector<char> {
        auto it = IUPAC_BASES.find(c);
        if (it != IUPAC_BASES.end()) return it->second;
        return {c};
    };

    bool has_degen = false;
    for (char c : mot)
        if (string("RYSWKMBDHV").find(c) != string::npos) { has_degen = true; break; }

    if (!has_degen) {
        trimers.push_back(mot.substr(0, 3));
        trimers.push_back(mot.substr(1, 3));
    } else {
        for (char x : expand(mot[0]))
            trimers.push_back(string(1, x) + mot[1] + mot[2]);
        for (char x : expand(mot[3]))
            trimers.push_back(string(1, mot[1]) + mot[2] + x);
    }
    return trimers;
}

static vector<int> id_align_mot(const string &m, const vector<string> &mall) {
    vector<string> mot_trimers = get_trimers(m);
    set<string> mot_set(mot_trimers.begin(), mot_trimers.end());

    vector<int> hits;
    for (int i = 0; i < (int)mall.size(); ++i) {
        vector<string> other_trimers = get_trimers(mall[i]);
        set<string> other_set(other_trimers.begin(), other_trimers.end());
        int overlap = 0;
        for (const auto &t : other_set)
            if (mot_set.count(t)) overlap++;
        if (other_set.size() > 0 &&
            (double)overlap / (double)other_set.size() >= 0.5)
            hits.push_back(i);
    }
    return hits;
}

// Pearson correlation between two vectors
static double pearson_cor(const vector<double> &a, const vector<double> &b) {
    int n = (int)a.size();
    if (n == 0) return 0.0;
    double ma = 0, mb = 0;
    for (int i = 0; i < n; ++i) { ma += a[i]; mb += b[i]; }
    ma /= n; mb /= n;
    double num = 0, da = 0, db = 0;
    for (int i = 0; i < n; ++i) {
        double x = a[i] - ma, y = b[i] - mb;
        num += x * y;
        da += x * x;
        db += y * y;
    }
    double denom = sqrt(da * db);
    return denom > 0 ? num / denom : 0.0;
}

// Extract column from flat matrix
static vector<double> get_column(const vector<double> &mat, int n_rows, int n_cols, int col) {
    vector<double> out(n_rows);
    for (int i = 0; i < n_rows; ++i) out[i] = mat[i * n_cols + col];
    return out;
}

struct ClusterResult {
    vector<string> order;
    vector<int> cluster_ids;
};

static ClusterResult sorting_and_type(const vector<double> &score_sort,
                                      int n_pos,
                                      vector<string> tets) {
    ClusterResult cr;
    if (tets.empty()) return cr;
    if (tets.size() == 1) {
        cr.order = tets;
        cr.cluster_ids = {1};
        return cr;
    }

    // Build score columns for each tetramer (for correlation)
    int n_tets = (int)tets.size();
    map<string, int> tet_idx;
    for (int i = 0; i < n_tets; ++i) tet_idx[tets[i]] = i;

    int cluster_id = 1;
    vector<string> remaining = tets;

    while (!remaining.empty()) {
        string m = remaining[0];
        remaining.erase(remaining.begin());

        if (remaining.empty()) {
            cr.order.push_back(m);
            cr.cluster_ids.push_back(cluster_id);
            break;
        }

        vector<int> ids = id_align_mot(m, remaining);

        vector<string> add;
        add.push_back(m);

        if (ids.empty()) {
            // No overlap: singleton
        } else if (ids.size() == 1) {
            add.push_back(remaining[ids[0]]);
            remaining.erase(remaining.begin() + ids[0]);
        } else {
            // Sort by Pearson correlation with m
            int mi = tet_idx[m];
            vector<double> m_score = get_column(score_sort, n_pos, n_tets, mi);

            vector<pair<double, int>> cor_vals;
            for (int idx : ids) {
                int oi = tet_idx[remaining[idx]];
                vector<double> o_score = get_column(score_sort, n_pos, n_tets, oi);
                double cor = pearson_cor(m_score, o_score);
                cor_vals.push_back({cor, idx});
            }
            sort(cor_vals.begin(), cor_vals.end(),
                 [](const pair<double,int> &a, const pair<double,int> &b) {
                     return a.first > b.first;
                 });

            // Add in correlation order, then remove from remaining (reverse order)
            for (const auto &cv : cor_vals)
                add.push_back(remaining[cv.second]);

            // Remove matched indices (sorted descending to preserve indices)
            vector<int> to_remove;
            for (const auto &cv : cor_vals) to_remove.push_back(cv.second);
            sort(to_remove.rbegin(), to_remove.rend());
            for (int idx : to_remove)
                remaining.erase(remaining.begin() + idx);
        }

        for (const auto &a : add) {
            cr.order.push_back(a);
            cr.cluster_ids.push_back(cluster_id);
        }
        cluster_id++;
    }
    return cr;
}

// ═══════════════════════════════════════════════════════════════════════════
// Region classification (sign_reg_plot.R logic for --mars mode)
// ═══════════════════════════════════════════════════════════════════════════

struct RegionClass {
    string R1, R2, R3;
};

static map<string, RegionClass> classify_regions(
        const vector<BootstrapRow> &boot,
        const vector<string> &enriched_tets,
        double p_emp_cutoff) {
    // Build adaptive p_cutoff per region (1st percentile, capped at 0.05)
    // Combine enh and sil Fisher p-values
    vector<double> r1_pf, r2_pf, r3_pf;
    for (const auto &b : boot) {
        r1_pf.push_back(b.r_enh_pFis[0]); r1_pf.push_back(b.r_sil_pFis[0]);
        r2_pf.push_back(b.r_enh_pFis[1]); r2_pf.push_back(b.r_sil_pFis[1]);
        r3_pf.push_back(b.r_enh_pFis[2]); r3_pf.push_back(b.r_sil_pFis[2]);
    }

    auto quantile_01 = [](vector<double> v) -> double {
        if (v.empty()) return 0.05;
        sort(v.begin(), v.end());
        int idx = max(0, (int)(v.size() * 0.01) - 1);
        double q = v[idx];
        return q < 0.05 ? q : 0.05;
    };

    double pc1 = quantile_01(r1_pf);
    double pc2 = quantile_01(r2_pf);
    double pc3 = quantile_01(r3_pf);

    // Build tetramer lookup
    map<string, const BootstrapRow*> boot_map;
    for (const auto &b : boot) boot_map[b.tetramer] = &b;

    map<string, RegionClass> result;
    for (const auto &tet : enriched_tets) {
        auto it = boot_map.find(tet);
        if (it == boot_map.end()) {
            result[tet] = {"0", "0", "0"};
            continue;
        }
        const auto &b = *it->second;
        RegionClass rc = {"0", "0", "0"};

        auto classify = [&](int r, double pc, string &out) {
            bool enh = b.r_enh_pFis[r] < pc && b.r_enh_pEmp[r] <= p_emp_cutoff;
            bool sil = b.r_sil_pFis[r] < pc && b.r_sil_pEmp[r] <= p_emp_cutoff;
            if (enh && sil) out = "Both";
            else if (enh) out = "Enhanced";
            else if (sil) out = "Silenced";
        };

        classify(0, pc1, rc.R1);
        classify(1, pc2, rc.R2);
        classify(2, pc3, rc.R3);
        result[tet] = rc;
    }
    return result;
}

// ═══════════════════════════════════════════════════════════════════════════
// Main
// ═══════════════════════════════════════════════════════════════════════════

int main(int argc, char *argv[]) {
    if (argc < 6) {
        cerr << "Usage: " << argv[0]
             << " <results_dir> <name> <n_bootstraps> <p_fisher> <p_empirical>"
             << " [top_n] [counts_dir] [in_exon] [in_intron] [--mars] [--mars-groups]"
             << endl;
        return 1;
    }

    string results_dir = argv[1];
    string name        = argv[2];
    int    n_boot      = stoi(argv[3]);
    double cFisher     = stod(argv[4]);
    double cEmp        = stod(argv[5]);
    int    top_n       = (argc >= 7 && argv[6][0] != '-') ? stoi(argv[6]) : 0;
    string counts_dir  = (argc >= 8 && argv[7][0] != '-') ? argv[7] : "";
    int    in_exon     = (argc >= 9 && argv[8][0] != '-') ? stoi(argv[8]) : 30;
    int    in_intron   = (argc >= 10 && argv[9][0] != '-') ? stoi(argv[9]) : 300;

    bool mars_mode = false, mars_groups = false;
    for (int i = 1; i < argc; ++i) {
        if (string(argv[i]) == "--mars") mars_mode = true;
        if (string(argv[i]) == "--mars-groups") { mars_mode = true; mars_groups = true; }
    }

    string pp = results_dir + "/";
    string cp = counts_dir.empty() ? pp : counts_dir + "/";

    // ── Region spacing ──────────────────────────────────────────────────
    int spacing = compute_spacing(in_exon, in_intron);
    vector<int> regions = compute_regions(spacing);
    vector<int> rown = compute_rown(in_exon, in_intron, regions);
    int n_pos = (int)rown.size();

    cerr << "[selection] in_exon=" << in_exon << " in_intron=" << in_intron
         << " spacing=" << spacing << " n_positions=" << n_pos << endl;

    // ── Read bootstrap results ──────────────────────────────────────────
    string boot_path = pp + "bootstrap_" + to_string(n_boot) + ".tsv";
    cerr << "[selection] Loading " << boot_path << endl;
    vector<BootstrapRow> boot = read_bootstrap_tsv(boot_path);
    if (boot.empty()) {
        cerr << "Error: no bootstrap data loaded" << endl;
        return 1;
    }
    cerr << "[selection] " << boot.size() << " tetramers loaded" << endl;

    // ── Select significantly enriched tetramers ─────────────────────────
    vector<string> sig_tets;
    for (const auto &b : boot) {
        bool sig = false;
        for (int r = 0; r < 3; ++r) {
            if ((b.r_enh_pFis[r] <= cFisher && b.r_enh_pEmp[r] <= cEmp) ||
                (b.r_sil_pFis[r] <= cFisher && b.r_sil_pEmp[r] <= cEmp)) {
                sig = true; break;
            }
        }
        if (sig) sig_tets.push_back(b.tetramer);
    }

    if (sig_tets.empty()) {
        cerr << "[selection] No significant tetramers found." << endl;
        return 0;
    }
    cerr << "[selection] " << sig_tets.size() << " significant tetramers" << endl;

    // ── Exon counts from first region_count file ────────────────────────
    int CE_enh = 0, CE_sil = 0, CE_ctrl = 0;
    {
        // Read filelist_count.tsv to get the first region_count file
        for (const string &sub : {"nr/", "r/"}) {
            string flist_path = cp + sub + "filelist_count.tsv";
            ifstream fl(flist_path);
            if (!fl.is_open()) continue;
            string first_file;
            if (getline(fl, first_file)) {
                first_file.erase(first_file.find_last_not_of(" \t\r\n") + 1);
                string fpath = cp + sub + first_file;
                ifstream cf(fpath);
                if (cf.is_open()) {
                    string hdr;
                    getline(cf, hdr); // skip header
                    string line;
                    while (getline(cf, line)) {
                        istringstream ss(line);
                        string tok;
                        vector<string> cols;
                        while (getline(ss, tok, '\t')) cols.push_back(tok);
                        if (cols.size() >= 3) {
                            int ty = stoi(cols[2]);
                            if (ty == 1) CE_enh++;
                            else if (ty == -1) CE_sil++;
                            else if (ty == 0) CE_ctrl++;
                        }
                    }
                }
            }
            break;
        }
    }
    cerr << "[selection] Exons: " << CE_enh << " enh, " << CE_sil
         << " sil, " << CE_ctrl << " ctrl" << endl;

    // ── Positional Fisher tests ─────────────────────────────────────────
    cerr << "[selection] Computing positional Fisher tests..." << endl;
    PositionalData pd = get_tables(cp, sig_tets, rown);
    if (pd.n_tets == 0) {
        cerr << "Error: no tetramer positional data loaded" << endl;
        return 1;
    }

    vector<double> p_enh = positional_fisher(pd.enh, CE_enh, pd.cont, CE_ctrl,
                                              pd.n_pos, pd.n_tets);
    vector<double> p_sil = positional_fisher(pd.sil, CE_sil, pd.cont, CE_ctrl,
                                              pd.n_pos, pd.n_tets);

    // ── Fisher's method combined score ──────────────────────────────────
    // ES[pos,tet] = -2 * (log(p_enh) + log(p_sil))
    // e[pos,tet]  = -2 * log(p_enh)
    // s[pos,tet]  = -2 * log(p_sil)
    int total = pd.n_pos * pd.n_tets;
    vector<double> ES(total), e_score(total), s_score(total);
    for (int i = 0; i < total; ++i) {
        double le = (p_enh[i] > 0) ? log(p_enh[i]) : -700.0;
        double ls = (p_sil[i] > 0) ? log(p_sil[i]) : -700.0;
        e_score[i] = -2.0 * le;
        s_score[i] = -2.0 * ls;
        ES[i] = e_score[i] + s_score[i];
    }

    // ── Sort by AUC (column sum of ES) ──────────────────────────────────
    vector<pair<double, int>> auc_order;
    for (int ti = 0; ti < pd.n_tets; ++ti) {
        double auc = 0;
        for (int pi = 0; pi < pd.n_pos; ++pi)
            auc += ES[pi * pd.n_tets + ti];
        auc_order.push_back({auc, ti});
    }
    sort(auc_order.begin(), auc_order.end(),
         [](const pair<double,int> &a, const pair<double,int> &b) {
             return a.first > b.first;
         });

    vector<string> sorted_tets;
    for (const auto &ao : auc_order)
        sorted_tets.push_back(pd.tet_names[ao.second]);

    // ── Score sort matrix (e - s) for clustering ────────────────────────
    vector<double> score_sort(total);
    for (int i = 0; i < total; ++i)
        score_sort[i] = e_score[i] - s_score[i];

    // ── Cluster tetramers ───────────────────────────────────────────────
    cerr << "[selection] Clustering " << sorted_tets.size() << " tetramers..." << endl;
    ClusterResult cr = sorting_and_type(score_sort, pd.n_pos, sorted_tets);

    // Build tetramer -> cluster_id map
    map<string, int> tet_cluster;
    for (int i = 0; i < (int)cr.order.size(); ++i)
        tet_cluster[cr.order[i]] = cr.cluster_ids[i];

    // ── Build summary table (MRMs CSV) ──────────────────────────────────
    // Build bootstrap lookup
    map<string, const BootstrapRow*> boot_map;
    for (const auto &b : boot) boot_map[b.tetramer] = &b;

    auto sanitize = [](double v) -> string {
        ostringstream ss;
        ss << fixed << setprecision(6) << v;
        string s = ss.str();
        // Replace dots with underscores for filename
        return s;
    };

    string san_emp = to_string(cEmp);
    string san_fish = to_string(cFisher);
    // Remove trailing zeros
    auto trim_zeros = [](string s) {
        size_t dot = s.find('.');
        if (dot != string::npos) {
            size_t last = s.find_last_not_of('0');
            if (last == dot) last++;
            s = s.substr(0, last + 1);
        }
        // Replace . with _
        for (auto &c : s) if (c == '.') c = '_';
        return s;
    };
    san_emp = trim_zeros(san_emp);
    san_fish = trim_zeros(san_fish);

    string base_name = "MRMs_" + name + "_emp-" + san_emp +
                       "_fisher-" + san_fish + "_nBoot-" + to_string(n_boot);

    // ── Write tetramer_order.txt ────────────────────────────────────────
    {
        string path = pp + "tetramer_order.txt";
        ofstream out(path);
        out << "tetramer\tcluster_id\n";
        for (int i = 0; i < (int)cr.order.size(); ++i)
            out << cr.order[i] << "\t" << cr.cluster_ids[i] << "\n";
        cerr << "[selection] Wrote " << path << endl;
    }

    // Apply top_n filter on display order
    vector<string> display_order = cr.order;
    if (top_n > 0 && top_n < (int)display_order.size())
        display_order.resize(top_n);

    // ── Write MRMs CSV ──────────────────────────────────────────────────
    {
        string csv_path = pp + base_name + ".csv";
        ofstream out(csv_path);
        out << "tetramer,full_motifs,exonType,is.sign,where.sign,"
            << "r1_pf,r1_pe,r2_pf,r2_pe,r3_pf,r3_pe,cluster_id\n";

        set<string> sig_set(sig_tets.begin(), sig_tets.end());

        for (const string &tet : cr.order) {
            auto it = boot_map.find(tet);
            if (it == boot_map.end()) continue;
            const auto &b = *it->second;

            for (const string &exon_type : {"enh", "sil"}) {
                string where;
                bool is_enh = (exon_type == "enh");
                for (int r = 0; r < 3; ++r) {
                    double pf = is_enh ? b.r_enh_pFis[r] : b.r_sil_pFis[r];
                    double pe = is_enh ? b.r_enh_pEmp[r] : b.r_sil_pEmp[r];
                    if (pf <= cFisher && pe <= cEmp) {
                        if (!where.empty()) where += ",";
                        where += to_string(r + 1);
                    }
                }
                bool is_sig = !where.empty();

                out << tet << "," << iupac_full(tet) << "," << exon_type
                    << "," << (is_sig ? "TRUE" : "FALSE") << "," << where;
                if (is_enh) {
                    out << "," << b.r_enh_pFis[0] << "," << b.r_enh_pEmp[0]
                        << "," << b.r_enh_pFis[1] << "," << b.r_enh_pEmp[1]
                        << "," << b.r_enh_pFis[2] << "," << b.r_enh_pEmp[2];
                } else {
                    out << "," << b.r_sil_pFis[0] << "," << b.r_sil_pEmp[0]
                        << "," << b.r_sil_pFis[1] << "," << b.r_sil_pEmp[1]
                        << "," << b.r_sil_pFis[2] << "," << b.r_sil_pEmp[2];
                }
                out << "," << tet_cluster[tet] << "\n";
            }
        }
        cerr << "[selection] Wrote " << csv_path << endl;
    }

    // ── Mars mode outputs ───────────────────────────────────────────────
    if (mars_mode) {
        // enriched_tetramers.txt
        {
            string path = pp + "enriched_tetramers.txt";
            ofstream out(path);
            for (const auto &t : cr.order) out << t << "\n";
            cerr << "[selection] Wrote " << path << endl;
        }

        // cluster_ids.txt
        {
            string path = pp + "cluster_ids.txt";
            ofstream out(path);
            for (int i = 0; i < (int)cr.order.size(); ++i)
                out << cr.order[i] << "\t" << cr.cluster_ids[i] << "\n";
            cerr << "[selection] Wrote " << path << endl;
        }

        // regions.csv
        {
            auto region_map = classify_regions(boot, cr.order, cEmp);
            string path = pp + "regions.csv";
            ofstream out(path);
            out << "\"\",\"R1\",\"R2\",\"R3\"\n";
            for (const auto &tet : cr.order) {
                const auto &rc = region_map[tet];
                out << "\"" << tet << "\",\"" << rc.R1 << "\",\""
                    << rc.R2 << "\",\"" << rc.R3 << "\"\n";
            }
            cerr << "[selection] Wrote " << path << endl;
        }

        // e-*.csv, s-*.csv, ES-*.csv (positional score matrices)
        // Columns = tetramers, Rows = positions
        auto write_score_csv = [&](const string &prefix, const vector<double> &scores) {
            string path = pp + prefix + "-" + name + ".csv";
            ofstream out(path);
            // Header: tetramer names
            for (int ti = 0; ti < pd.n_tets; ++ti) {
                if (ti > 0) out << ",";
                out << pd.tet_names[ti];
            }
            out << "\n";
            // Data rows
            for (int pi = 0; pi < pd.n_pos; ++pi) {
                for (int ti = 0; ti < pd.n_tets; ++ti) {
                    if (ti > 0) out << ",";
                    out << fixed << setprecision(6) << scores[pi * pd.n_tets + ti];
                }
                out << "\n";
            }
            cerr << "[selection] Wrote " << path << endl;
        };

        write_score_csv("e", e_score);
        write_score_csv("s", s_score);
        write_score_csv("ES", ES);
    }

    // ── Mars groups mode outputs ────────────────────────────────────────
    if (mars_groups) {
        auto region_map = classify_regions(boot, cr.order, cEmp);

        // Partition tetramers into groups
        vector<string> all_tets = cr.order;
        vector<string> sil_tets, enh_tets;
        for (const auto &tet : cr.order) {
            const auto &rc = region_map[tet];
            bool has_sil = (rc.R1 == "Silenced" || rc.R1 == "Both" ||
                           rc.R2 == "Silenced" || rc.R2 == "Both" ||
                           rc.R3 == "Silenced" || rc.R3 == "Both");
            bool has_enh = (rc.R1 == "Enhanced" || rc.R1 == "Both" ||
                           rc.R2 == "Enhanced" || rc.R2 == "Both" ||
                           rc.R3 == "Enhanced" || rc.R3 == "Both");
            if (has_sil) sil_tets.push_back(tet);
            if (has_enh) enh_tets.push_back(tet);
        }

        // For each group, compute pooled positional Fisher and write CSVs
        auto compute_group = [&](const vector<string> &group_tets,
                                 const string &group_prefix) {
            if (group_tets.empty()) return;

            // Pool counts across all tetramers in group
            vector<double> pool_enh(n_pos, 0), pool_sil(n_pos, 0), pool_cont(n_pos, 0);
            for (const auto &tet : group_tets) {
                // Find tetramer index
                int ti = -1;
                for (int i = 0; i < pd.n_tets; ++i)
                    if (pd.tet_names[i] == tet) { ti = i; break; }
                if (ti < 0) continue;
                for (int pi = 0; pi < n_pos; ++pi) {
                    int idx = pi * pd.n_tets + ti;
                    pool_enh[pi]  += pd.enh[idx];
                    pool_sil[pi]  += pd.sil[idx];
                    pool_cont[pi] += pd.cont[idx];
                }
            }

            // Compute Fisher at each position for the pooled group
            int max_enh = (int)*max_element(pool_enh.begin(), pool_enh.end());
            int max_sil = (int)*max_element(pool_sil.begin(), pool_sil.end());
            int max_cont = (int)*max_element(pool_cont.begin(), pool_cont.end());
            int tot_enh = max(CE_enh, max_enh);
            int tot_sil = max(CE_sil, max_sil);
            int tot_cont = max(CE_ctrl, max_cont);

            vector<double> group_e(n_pos), group_s(n_pos);
            for (int pi = 0; pi < n_pos; ++pi) {
                double pe = fisher_greater((int)pool_enh[pi], (int)pool_cont[pi],
                                           tot_enh - (int)pool_enh[pi],
                                           tot_cont - (int)pool_cont[pi]);
                double ps = fisher_greater((int)pool_sil[pi], (int)pool_cont[pi],
                                           tot_sil - (int)pool_sil[pi],
                                           tot_cont - (int)pool_cont[pi]);
                group_e[pi] = -2.0 * ((pe > 0) ? log(pe) : -700.0);
                group_s[pi] = -2.0 * ((ps > 0) ? log(ps) : -700.0);
            }

            // Write CSVs
            auto write_group_csv = [&](const string &suffix, const vector<double> &vals) {
                string col_name = "Group_" + group_prefix;
                string path = pp + suffix + "_group_" + group_prefix + "-" + name + ".csv";
                ofstream out(path);
                out << col_name << "\n";
                for (int pi = 0; pi < n_pos; ++pi)
                    out << fixed << setprecision(6) << vals[pi] << "\n";
                cerr << "[selection] Wrote " << path << endl;
            };

            write_group_csv("Enh" == group_prefix ? "Enh" :
                           ("Sil" == group_prefix ? "Sil" : "All"), group_e);

            // Also write sil version
            string sil_prefix = (group_prefix == "both") ? "All" :
                               ((group_prefix == "sil") ? "Sil" : "Enh");
            string path_sil = pp + sil_prefix + "_group_sil-" + name + ".csv";
            ofstream out_sil(path_sil);
            out_sil << "Group_" << group_prefix << "\n";
            for (int pi = 0; pi < n_pos; ++pi)
                out_sil << fixed << setprecision(6) << group_s[pi] << "\n";

            string path_both = pp + sil_prefix + "_group_both-" + name + ".csv";
            ofstream out_both(path_both);
            out_both << "Group_" << group_prefix << "\n";
            for (int pi = 0; pi < n_pos; ++pi)
                out_both << fixed << setprecision(6) << group_e[pi] + group_s[pi] << "\n";
        };

        compute_group(all_tets, "both");
        compute_group(sil_tets, "sil");
        compute_group(enh_tets, "enh");
    }

    cerr << "[selection] Done." << endl;
    return 0;
}
