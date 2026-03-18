// RNAmotifs - Tetramer genome search (C++ with OpenMP)
// Tetramer genome search: regions, motif finding, clustering, thresholding.
// Copyright (C) 2014-2026 Matteo Cereda
// SPDX-License-Identifier: GPL-2.0-or-later

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
namespace fs = std::filesystem;

// ===========================================================================
// IUPAC codes
// ===========================================================================

static const unordered_map<char, vector<char>> IUPAC = {
    {'R', {'A', 'G'}}, {'Y', {'C', 'T'}},
    {'S', {'G', 'C'}}, {'W', {'A', 'T'}}
};

static const unordered_map<char, char> COMPLEMENT = {
    {'A','T'}, {'T','A'}, {'U','A'}, {'G','C'}, {'C','G'},
    {'R','Y'}, {'Y','R'}, {'K','M'}, {'M','K'}, {'S','S'}, {'W','W'},
    {'B','V'}, {'D','H'}, {'H','D'}, {'V','B'}, {'N','N'}
};

// ===========================================================================
// Config (mirrors config.py)
// ===========================================================================

struct SearchConfig {
    string splicing_file;
    string genome_dir;
    string genome;
    string regions_name;
    string root;
    int    cluster_hw = 15;
    int    h_min      = 4;
    vector<double> pth = {0.5};

    string regions_path() const {
        return root + "/regions/" + regions_name + "/" + regions_name + ".tab";
    }
    string bed_folder() const {
        return root + "/results/" + regions_name + "/motifs.raw/";
    }
    string bed_path(const string &motif) const {
        return bed_folder() + motif + ".bed";
    }
    string stats_path() const {
        return root + "/results/" + regions_name + "/motifs.final/stats.txt";
    }
    string pth_folder(double p) const {
        // Match Python's f"pth_{pth_val}" exactly
        ostringstream ss;
        ss << root << "/results/" << regions_name << "/motifs.final/pth_" << p << "/";
        return ss.str();
    }
    string pth_path(double p, const string &motif) const {
        return pth_folder(p) + motif + ".bed";
    }
};

// ===========================================================================
// Region data
// ===========================================================================

struct Region {
    int    rid;
    int    start;
    int    stop;
    string rclass;
};

using ChromStrand = pair<string, string>;
using RegionsDB = map<ChromStrand, vector<Region>>;

// ===========================================================================
// Chromosome cache (read once, share across all motifs)
// ===========================================================================

static unordered_map<string, string> chrom_cache;
static mutex chrom_mutex;

static void preload_chromosome(const string &genome_dir, const string &genome,
                               const string &chrom) {
    string path = genome_dir + "/genomes/" + genome + "/" + chrom + ".string";
    ifstream fh(path, ios::binary);
    if (!fh.is_open()) {
        cerr << "  Warning: " << path << " not found\n";
        lock_guard<mutex> lock(chrom_mutex);
        chrom_cache[chrom] = "";
        return;
    }
    fh.seekg(0, ios::end);
    size_t sz = fh.tellg();
    fh.seekg(0, ios::beg);
    string content(sz, '\0');
    fh.read(&content[0], sz);
    lock_guard<mutex> lock(chrom_mutex);
    chrom_cache[chrom] = move(content);
}

static string get_subseq(const string &chrom, int start, int end) {
    auto it = chrom_cache.find(chrom);
    if (it == chrom_cache.end() || it->second.empty()) return "";
    const string &seq = it->second;
    if (start < 0) start = 0;
    if (end >= (int)seq.size()) end = (int)seq.size() - 1;
    if (start > end) return "";
    return seq.substr(start, end - start + 1);
}

// ===========================================================================
// Reverse complement
// ===========================================================================

static string reverse_complement(const string &seq) {
    string result(seq.size(), 'N');
    for (int i = (int)seq.size() - 1, j = 0; i >= 0; i--, j++) {
        char c = toupper(seq[i]);
        auto it = COMPLEMENT.find(c);
        result[j] = (it != COMPLEMENT.end()) ? it->second : c;
    }
    return result;
}

// ===========================================================================
// Tetramer enumeration (matches genome.py order exactly)
// ===========================================================================

static vector<string> enumerate_tetramers() {
    const char bases[] = {'A', 'C', 'T', 'G'};
    const char degen[] = {'R', 'Y', 'S', 'W'};
    vector<string> result;
    result.reserve(512);
    // 256 non-redundant
    for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++)
            for (int c = 0; c < 4; c++)
                for (int d = 0; d < 4; d++)
                    result.push_back({bases[a], bases[b], bases[c], bases[d]});
    // 256 redundant (DEGEN + base + base + DEGEN)
    for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++)
            for (int x = 0; x < 4; x++)
                for (int y = 0; y < 4; y++)
                    result.push_back({degen[x], bases[a], bases[b], degen[y]});
    return result;
}

// ===========================================================================
// Motif search (matches search.py exactly)
// ===========================================================================

static vector<string> expand_motif(const string &motif) {
    vector<char> first_opts, last_opts;
    auto f = IUPAC.find(motif[0]);
    first_opts = (f != IUPAC.end()) ? f->second : vector<char>{motif[0]};
    auto l = IUPAC.find(motif.back());
    last_opts = (l != IUPAC.end()) ? l->second : vector<char>{motif.back()};

    string mid = motif.substr(1, motif.size() - 2);
    vector<string> result;
    for (char fc : first_opts)
        for (char lc : last_opts)
            result.push_back(string(1, fc) + mid + string(1, lc));
    return result;
}

static vector<int> find_all_indices(const string &seq, const string &sub) {
    vector<int> indices;
    size_t pos = seq.find(sub, 0);
    while (pos != string::npos) {
        indices.push_back((int)pos);
        pos = seq.find(sub, pos + 1);
    }
    return indices;
}

static vector<pair<int,int>> find_motif(const string &seq, const string &motif) {
    vector<pair<int,int>> pool;
    bool has_degen = false;
    for (char c : motif)
        if (c == 'R' || c == 'Y' || c == 'S' || c == 'W') { has_degen = true; break; }

    if (!has_degen) {
        for (int i : find_all_indices(seq, motif))
            pool.emplace_back(i, i + (int)motif.size() - 1);
    } else {
        for (const string &exp : expand_motif(motif))
            for (int i : find_all_indices(seq, exp))
                pool.emplace_back(i, i + (int)exp.size() - 1);
    }

    sort(pool.begin(), pool.end());
    if (pool.size() <= 1) return pool;

    vector<pair<int,int>> merged;
    merged.push_back(pool[0]);
    for (size_t k = 1; k < pool.size(); k++) {
        if (pool[k].first <= merged.back().second + 1)
            merged.back().second = max(merged.back().second, pool[k].second);
        else
            merged.push_back(pool[k]);
    }
    return merged;
}

// ===========================================================================
// Interval merging (matches utils.py)
// ===========================================================================

static vector<pair<int,int>> merge_intervals(vector<pair<int,int>> &data) {
    if (data.empty()) return data;
    sort(data.begin(), data.end());
    vector<pair<int,int>> result;
    auto saved = data[0];
    for (size_t i = 1; i < data.size(); i++) {
        if (data[i].first <= saved.second)
            saved.second = max(saved.second, data[i].second);
        else { result.push_back(saved); saved = data[i]; }
    }
    result.push_back(saved);
    return result;
}

// ===========================================================================
// Bedgraph (matches utils.py)
// ===========================================================================

// Dense array per chrom+strand for O(1) prefix-sum queries
struct DenseArray {
    vector<double> data;
    vector<double> prefix;
    int offset;  // min position
    int sz;      // number of positions

    DenseArray() : offset(0), sz(0) {}
    DenseArray(int min_pos, int max_pos)
        : offset(min_pos), sz(max_pos - min_pos + 1) {
        data.assign(sz, 0.0);
    }

    void add(int pos, double val) {
        int idx = pos - offset;
        if (idx >= 0 && idx < sz) data[idx] += val;
    }

    void build_prefix() {
        prefix.resize(sz + 1, 0.0);
        for (int i = 0; i < sz; i++)
            prefix[i + 1] = prefix[i] + data[i];
    }

    // Sum of values in positions [lo, hi] (inclusive)
    double region_sum(int lo, int hi) const {
        int a = max(0, lo - offset);
        int b = min(sz, hi - offset + 1);
        if (a >= b) return 0.0;
        return prefix[b] - prefix[a];
    }

    double get_value(int pos) const {
        int idx = pos - offset;
        if (idx >= 0 && idx < sz) return data[idx];
        return 0.0;
    }
};

using ChromStrandKey = pair<string, string>;

struct Bedgraph {
    map<ChromStrandKey, DenseArray> arrays;

    // Initialize dense arrays with known position ranges from the regions DB
    void init_from_regions(const RegionsDB &rdb, int hw) {
        for (auto &[key, regions] : rdb) {
            if (regions.empty()) continue;
            int min_pos = regions[0].start;
            int max_pos = regions[0].stop;
            for (auto &r : regions) {
                min_pos = min(min_pos, r.start - hw * 2);
                max_pos = max(max_pos, r.stop + hw * 2);
            }
            arrays[key] = DenseArray(min_pos, max_pos);
        }
    }

    void load(const string &filename) {
        ifstream f(filename);
        if (!f.is_open()) return;
        string line;
        while (getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            istringstream ss(line);
            string chrom, s1, s2, s3;
            getline(ss, chrom, '\t');
            getline(ss, s1, '\t');
            getline(ss, s2, '\t');
            getline(ss, s3, '\t');
            int start = stoi(s1), stop = stoi(s2);
            double value = stod(s3);
            string strand = (value >= 0) ? "+" : "-";
            ChromStrandKey key{chrom, strand};
            auto it = arrays.find(key);
            if (it == arrays.end()) continue;
            for (int p = start; p < stop; p++)
                it->second.add(p, fabs(value));
        }
    }

    void build_all_prefix() {
        for (auto &[key, arr] : arrays)
            arr.build_prefix();
    }

    double get_value(const string &chrom, const string &strand, int pos) const {
        auto it = arrays.find({chrom, strand});
        if (it == arrays.end()) return 0.0;
        return it->second.get_value(pos);
    }

    double region_sum(const string &chrom, const string &strand, int lo, int hi) const {
        auto it = arrays.find({chrom, strand});
        if (it == arrays.end()) return 0.0;
        return it->second.region_sum(lo, hi);
    }

    // Cluster using prefix sums: O(1) per position instead of O(2*hw)
    void cluster(int hw) {
        // First build prefix sums on the raw data
        build_all_prefix();

        // Now replace each position's value with the window sum
        for (auto &[key, arr] : arrays) {
            vector<double> new_data(arr.sz, 0.0);
            for (int i = 0; i < arr.sz; i++) {
                if (arr.data[i] > 0.0) {
                    int pos = i + arr.offset;
                    new_data[i] = arr.region_sum(max(0, pos - hw), pos + hw);
                }
            }
            arr.data = move(new_data);
            // Rebuild prefix for any subsequent queries
            arr.build_prefix();
        }
    }
};

// ===========================================================================
// write_chr_positions (matches search.py)
// ===========================================================================

static void write_chr_positions(const map<int, int> &pos_dict,
                                const string &chrom, const string &strand,
                                ofstream &fw) {
    if (pos_dict.empty()) return;
    vector<int> positions;
    positions.reserve(pos_dict.size());
    for (auto &[p, _] : pos_dict) positions.push_back(p);
    sort(positions.begin(), positions.end());

    int start = positions[0], stop = positions[0];
    for (size_t i = 1; i < positions.size(); i++) {
        if (positions[i] == stop + 1) {
            stop = positions[i];
        } else {
            fw << chrom << "\t" << start << "\t" << (stop + 1)
               << "\t" << strand << (stop - start + 1) << "\n";
            start = stop = positions[i];
        }
    }
    fw << chrom << "\t" << start << "\t" << (stop + 1)
       << "\t" << strand << (stop - start + 1) << "\n";
}

// ===========================================================================
// make_bed (matches search.py)
// ===========================================================================

static void make_bed(const string &motif, const SearchConfig &cfg,
                     const RegionsDB &rdb) {
    ofstream f(cfg.bed_path(motif));
    for (auto &[key, regions] : rdb) {
        auto &[chrom, strand] = key;
        vector<pair<int,int>> results;
        for (auto &r : regions) {
            int ext_start = r.start - cfg.cluster_hw * 2;
            int ext_stop  = r.stop  + cfg.cluster_hw * 2;
            string seq = get_subseq(chrom, ext_start, ext_stop);
            if (seq.empty()) continue;
            string query = (strand == "+") ? motif : reverse_complement(motif);
            auto hits = find_motif(seq, query);
            for (auto &[ms, me] : hits)
                results.emplace_back(ext_start + ms, ext_start + me);
        }
        auto merged = merge_intervals(results);
        for (auto &[s, e] : merged)
            f << chrom << "\t" << s << "\t" << (e + 1) << "\t" << strand << "1\n";
    }
}

// ===========================================================================
// cluster_threshold (matches search.py)
// ===========================================================================

static vector<string> cluster_threshold(const string &motif,
                                        const SearchConfig &cfg,
                                        const RegionsDB &rdb) {
    long rlen = 0;
    for (auto &[key, regions] : rdb)
        for (auto &r : regions)
            rlen += r.stop - r.start + 1;

    Bedgraph bg;
    bg.init_from_regions(rdb, cfg.cluster_hw);
    bg.load(cfg.bed_path(motif));
    bg.cluster(cfg.cluster_hw);

    map<int, int> hc;
    for (auto &[key, regions] : rdb) {
        auto &[chrom, strand] = key;
        for (auto &r : regions)
            for (int i = r.start; i < r.stop; i++) {
                double v = bg.get_value(chrom, strand, i);
                if (v > 0) hc[(int)v]++;
            }
    }

    // Precompute cumulative sums from highest to lowest key (optimization #9)
    // hc is a map<int,int> so keys are already sorted ascending
    vector<pair<int, int>> hc_sorted(hc.begin(), hc.end());
    // Compute ge (>=key) by accumulating from the right
    vector<long> ge_vals(hc_sorted.size());
    {
        long cumsum = 0;
        for (int i = (int)hc_sorted.size() - 1; i >= 0; i--) {
            cumsum += hc_sorted[i].second;
            ge_vals[i] = cumsum;
        }
    }

    vector<string> stats;
    stats.push_back("regions_length=" + to_string(rlen) + "\n");
    stats.push_back("motif=" + motif + "\n");
    for (size_t i = 0; i < hc_sorted.size(); i++) {
        char buf[256];
        snprintf(buf, sizeof(buf), "|h>=%d|=%ld (%.3f %%)\n",
                 hc_sorted[i].first, ge_vals[i],
                 (double)ge_vals[i] / rlen * 100.0);
        stats.push_back(buf);
    }
    stats.push_back("\n");

    for (double pth : cfg.pth) {
        vector<pair<double, int>> distances;
        for (size_t i = 0; i < hc_sorted.size(); i++) {
            double ge_pct = (double)ge_vals[i] / rlen * 100.0;
            distances.emplace_back(fabs(pth - ge_pct), hc_sorted[i].first);
        }
        sort(distances.begin(), distances.end());
        int h_chosen = distances.empty() ? cfg.h_min
                     : max(cfg.h_min, distances[0].second);

        map<string, map<int,int>> plus_data, minus_data;
        for (auto &[key, regions] : rdb) {
            auto &[chrom, strand] = key;
            for (auto &r : regions)
                for (int i = r.start; i < r.stop; i++) {
                    double v = bg.get_value(chrom, strand, i);
                    if (v >= h_chosen) {
                        if (strand == "+") plus_data[chrom][i] = 1;
                        else               minus_data[chrom][i] = 1;
                    }
                }
        }

        ofstream fw(cfg.pth_path(pth, motif));
        set<string> all_chroms;
        for (auto &[c, _] : plus_data)  all_chroms.insert(c);
        for (auto &[c, _] : minus_data) all_chroms.insert(c);
        for (const string &chrom : all_chroms) {
            write_chr_positions(plus_data[chrom], chrom, "+", fw);
            write_chr_positions(minus_data[chrom], chrom, "-", fw);
        }
    }
    return stats;
}

// ===========================================================================
// Region preparation (matches regions.py)
// ===========================================================================

static void prepare_regions(const string &splicing_file, const string &output) {
    ifstream fin(splicing_file);
    ofstream fout(output);
    string line;
    while (getline(fin, line)) {
        istringstream iss(line);
        vector<string> cols;
        string tok;
        while (getline(iss, tok, ';')) cols.push_back(tok);
        if ((int)cols.size() < 9) continue;
        string chrom = cols[2];
        if (chrom.find("random") != string::npos) continue;
        string strand = cols[3];
        int sk_from = stoi(cols[4]);
        int in_from = stoi(cols[5]);
        int in_to   = stoi(cols[6]);
        int sk_to   = stoi(cols[7]);
        int coords[][2] = {
            {sk_from - 100, sk_from + 500},
            {in_from - 500, in_from + 100},
            {in_to - 100,   in_to + 500},
            {sk_to - 500,   sk_to + 100}
        };
        for (auto &c : coords)
            fout << chrom << "\t" << c[0] << "\t" << c[1] << "\t" << strand << "\n";
    }
}

static void merge_overlapping(const string &input, const string &output) {
    map<string, vector<pair<int,int>>> regions;
    ifstream fin(input);
    string line;
    while (getline(fin, line)) {
        istringstream ss(line);
        string chrom, s1, s2, strand;
        getline(ss, chrom, '\t'); getline(ss, s1, '\t');
        getline(ss, s2, '\t');    getline(ss, strand, '\t');
        int a = stoi(s1), b = stoi(s2);
        if (a > b) swap(a, b);
        regions[chrom + "_" + strand].emplace_back(a, b);
    }

    ofstream fout(output);
    fout << "id\tchrom\tstrand\tstart\tstop\tclass\n";
    int rid = 1;
    for (auto &[key, ivs] : regions) {
        sort(ivs.begin(), ivs.end());
        vector<pair<int,int>> stack = {ivs[0]};
        for (size_t i = 1; i < ivs.size(); i++) {
            if (ivs[i].first <= stack.back().second)
                stack.back().second = max(stack.back().second, ivs[i].second);
            else
                stack.push_back(ivs[i]);
        }
        size_t pos = key.rfind('_');
        string chrom = key.substr(0, pos), strand = key.substr(pos + 1);
        for (auto &[s, e] : stack) {
            fout << rid++ << "\t" << chrom << "\t" << strand
                 << "\t" << s << "\t" << e << "\tcase\n";
        }
    }
}

static RegionsDB read_regions(const string &path) {
    RegionsDB db;
    ifstream fin(path);
    string line;
    getline(fin, line); // header
    while (getline(fin, line)) {
        istringstream ss(line);
        string sr, ch, st, s1, s2, cl;
        getline(ss, sr, '\t'); getline(ss, ch, '\t');
        getline(ss, st, '\t'); getline(ss, s1, '\t');
        getline(ss, s2, '\t'); getline(ss, cl, '\t');
        db[{ch, st}].push_back({stoi(sr), stoi(s1), stoi(s2), cl});
    }
    return db;
}

// ===========================================================================
// Main
// ===========================================================================

int main(int argc, char *argv[]) {
    if (argc < 9) {
        cerr << "Usage: " << argv[0]
             << " <splicing_file> <genome_dir> <genome> <regions_name>"
             << " <data_root> <cluster_hw> <h_min> <pth> [n_cores]\n";
        return 1;
    }

    SearchConfig cfg;
    cfg.splicing_file = argv[1];
    cfg.genome_dir    = argv[2];
    cfg.genome        = argv[3];
    cfg.regions_name  = argv[4];
    cfg.root          = argv[5];
    cfg.cluster_hw    = stoi(argv[6]);
    cfg.h_min         = stoi(argv[7]);
    cfg.pth           = {stod(argv[8])};
    int n_cores       = (argc >= 10) ? stoi(argv[9]) : 1;

#ifdef _OPENMP
    omp_set_num_threads(n_cores);
#endif

    // Create output directories
    fs::create_directories(cfg.bed_folder());
    fs::create_directories(fs::path(cfg.stats_path()).parent_path());
    for (double p : cfg.pth)
        fs::create_directories(cfg.pth_folder(p));

    // Step 1-2: Regions
    string raw_path = cfg.root + "/regions/" + cfg.regions_name + "/overlapping.tab";
    fs::create_directories(fs::path(raw_path).parent_path());
    prepare_regions(cfg.splicing_file, raw_path);
    merge_overlapping(raw_path, cfg.regions_path());
    fs::remove(raw_path);

    // Step 3: Read regions
    RegionsDB rdb = read_regions(cfg.regions_path());

    // Pre-load chromosomes
    set<string> chrom_set;
    for (auto &[key, _] : rdb) chrom_set.insert(key.first);
    cerr << "Pre-loading " << chrom_set.size() << " chromosomes...\n";
    for (const string &ch : chrom_set)
        preload_chromosome(cfg.genome_dir, cfg.genome, ch);
    cerr << "Chromosomes loaded.\n";

    // Enumerate motifs
    vector<string> motifs = enumerate_tetramers();
    int total = (int)motifs.size();
    vector<vector<string>> all_stats(total);

    // Process motifs in parallel
    int done = 0;
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < total; i++) {
        make_bed(motifs[i], cfg, rdb);
        all_stats[i] = cluster_threshold(motifs[i], cfg, rdb);

        #pragma omp critical
        {
            done++;
            if (done % max(1, total / 100) == 0 || done == total)
                cerr << "PROGRESS " << done << " " << total << "\n";
        }
    }

    // Write stats (sequential)
    ofstream stats_out(cfg.stats_path());
    for (int i = 0; i < total; i++)
        for (const string &line : all_stats[i])
            stats_out << line;

    cerr << "PROGRESS " << total << " " << total << "\n";
    cerr << "Search complete.\n";
    return 0;
}
