// RNAmotifs - Sequence conservation profiling (PhyloP)
// Copyright (C) 2014-2026 Matteo Cereda
// SPDX-License-Identifier: GPL-2.0-or-later
//
// For each enriched tetramer, compute position-specific average PhyloP
// conservation scores across all exons with a tetramer occurrence.
// Reads pre-converted per-chromosome .phylop.bin files (float32 arrays).

#include <algorithm>
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
// PhyloP score cache (per-chromosome float arrays)
// ===========================================================================

static unordered_map<string, vector<float>> phylop_cache;
static unordered_map<string, size_t> phylop_len;
static mutex phylop_mutex;

static void preload_phylop(const string &genome_dir, const string &genome,
                           const string &chrom) {
    string path = genome_dir + "/genomes/" + genome + "/" + chrom + ".phylop.bin";
    ifstream fh(path, ios::binary);
    if (!fh.is_open()) {
        lock_guard<mutex> lock(phylop_mutex);
        phylop_cache[chrom] = {};
        phylop_len[chrom] = 0;
        return;
    }
    fh.seekg(0, ios::end);
    size_t sz = fh.tellg();
    size_t n = sz / sizeof(float);
    fh.seekg(0, ios::beg);
    vector<float> scores(n);
    fh.read(reinterpret_cast<char*>(scores.data()), sz);
    lock_guard<mutex> lock(phylop_mutex);
    phylop_len[chrom] = n;
    phylop_cache[chrom] = move(scores);
}

static float get_phylop(const string &chrom, int pos) {
    auto it = phylop_cache.find(chrom);
    if (it == phylop_cache.end() || it->second.empty()) return 0.0f;
    if (pos < 0 || pos >= (int)it->second.size()) return 0.0f;
    return it->second[pos];
}

// ===========================================================================
// Splicing file and BED reading (same as rnamotifs_structure.cpp)
// ===========================================================================

struct Exon {
    string chrom, strand;
    int skip_start, in_start, in_stop, skip_stop;
    int category, exon_id;
};

static vector<Exon> read_splicing_file(const string &path) {
    vector<Exon> exons;
    ifstream fin(path);
    string line;
    while (getline(fin, line)) {
        istringstream iss(line);
        vector<string> cols;
        string tok;
        while (getline(iss, tok, ';')) cols.push_back(tok);
        if ((int)cols.size() < 9) continue;
        Exon e;
        e.exon_id    = stoi(cols[0]);
        e.chrom      = cols[2];
        e.strand     = cols[3];
        e.skip_start = stoi(cols[4]);
        e.in_start   = stoi(cols[5]);
        e.in_stop    = stoi(cols[6]);
        e.skip_stop  = stoi(cols[7]);
        double dirank = stod(cols[8]);
        if (dirank >= 1.0)       e.category = 1;
        else if (dirank <= -1.0) e.category = -1;
        else                     e.category = 0;
        exons.push_back(e);
    }
    return exons;
}

struct BedHit { string chrom; int start, end; bool forward; };

static vector<BedHit> read_bed(const string &path) {
    vector<BedHit> hits;
    ifstream fin(path);
    string line;
    while (getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        istringstream ss(line);
        string chrom, s1, s2, s3;
        getline(ss, chrom, '\t'); getline(ss, s1, '\t');
        getline(ss, s2, '\t');    getline(ss, s3, '\t');
        hits.push_back({chrom, stoi(s1), stoi(s2), s3.find('-') == string::npos});
    }
    return hits;
}

static bool exon_has_hit(const Exon &e, const vector<BedHit> &hits) {
    int r1s = e.skip_start - 100, r4e = e.skip_stop + 100;
    bool fwd = (e.strand == "+");
    for (auto &h : hits)
        if (h.chrom == e.chrom && h.forward == fwd &&
            h.start <= r4e && h.end >= r1s) return true;
    return false;
}

// ===========================================================================
// Main
// ===========================================================================

int main(int argc, char *argv[]) {
    if (argc < 6) {
        cerr << "Usage: " << argv[0]
             << " <splicing_file> <genome_dir> <genome> <tetramers_dir>"
             << " <enriched_tets_file> <output_tsv> [n_cores]\n";
        return 1;
    }

    string splicing_file = argv[1];
    string genome_dir    = argv[2];
    string genome        = argv[3];
    string tet_dir       = argv[4];
    string enriched_file = argv[5];
    string output_tsv    = argv[6];
    int    n_cores       = (argc >= 8) ? stoi(argv[7]) : 1;
    int    in_exon       = (argc >= 9) ? stoi(argv[8]) : 30;
    int    in_intron     = (argc >= 10) ? stoi(argv[9]) : 300;
    int    spacing       = 2 * max(in_exon, in_intron) + 10;
    int    regions_c[4]  = {spacing, 2*spacing, 3*spacing, 4*spacing};

#ifdef _OPENMP
    omp_set_num_threads(n_cores);
#endif

    // Read enriched tetramers
    vector<string> enriched_tets;
    {
        ifstream fin(enriched_file);
        string line;
        while (getline(fin, line)) {
            line.erase(line.find_last_not_of(" \t\r\n") + 1);
            if (!line.empty()) enriched_tets.push_back(line);
        }
    }
    if (enriched_tets.empty()) {
        cerr << "No enriched tetramers\n"; return 1;
    }
    cerr << "Enriched tetramers: " << enriched_tets.size() << "\n";

    // Read exons
    vector<Exon> exons = read_splicing_file(splicing_file);
    cerr << "Exons: " << exons.size() << "\n";

    // Pre-load PhyloP scores
    set<string> chrom_set;
    for (auto &e : exons) chrom_set.insert(e.chrom);
    cerr << "Loading PhyloP scores for " << chrom_set.size() << " chromosomes...\n";
    for (auto &ch : chrom_set) {
        preload_phylop(genome_dir, genome, ch);
        if (phylop_len[ch] == 0)
            cerr << "  Warning: no PhyloP data for " << ch << "\n";
    }

    // RNA map positions: computed from in_exon/in_intron and region centers
    vector<int> map_positions;
    for (int i = regions_c[0] - in_exon; i <= regions_c[0] + in_intron; i++) map_positions.push_back(i);
    for (int i = regions_c[1] - in_intron; i <= regions_c[1] + in_exon; i++) map_positions.push_back(i);
    for (int i = regions_c[2] - in_exon; i <= regions_c[2] + in_intron; i++) map_positions.push_back(i);
    for (int i = regions_c[3] - in_intron; i <= regions_c[3] + in_exon; i++) map_positions.push_back(i);
    int n_pos = (int)map_positions.size();

    ofstream out(output_tsv);
    out << "tetramer\tposition\tcategory\tphylop_score\tn_exons\n";

    int done = 0, total = (int)enriched_tets.size();

    for (const string &tet : enriched_tets) {
        // Load BED hits
        vector<BedHit> hits;
        for (const string &sub : {"nr/", "r/"}) {
            string p = tet_dir + "/" + sub + tet + ".bed";
            if (fs::exists(p)) {
                auto h = read_bed(p);
                hits.insert(hits.end(), h.begin(), h.end());
            }
        }

        // Find exons with hits by category
        map<int, vector<int>> cat_exons;
        for (int ei = 0; ei < (int)exons.size(); ei++)
            if (exon_has_hit(exons[ei], hits))
                cat_exons[exons[ei].category].push_back(ei);

        for (int cat : {1, 0, -1}) {
            auto it = cat_exons.find(cat);
            if (it == cat_exons.end()) {
                for (int p : map_positions)
                    out << tet << "\t" << p << "\t" << cat << "\t0\t0\n";
                continue;
            }

            const vector<int> &eidxs = it->second;
            vector<double> avg_score(n_pos, 0.0);
            vector<int> counts(n_pos, 0);

            #pragma omp parallel for schedule(dynamic, 4)
            for (int pi = 0; pi < n_pos; pi++) {
                int map_pos = map_positions[pi];
                double score_sum = 0;
                int n_valid = 0;

                for (int ei : eidxs) {
                    const Exon &ex = exons[ei];
                    bool fwd = (ex.strand == "+");
                    int bounds[] = {ex.skip_start, ex.in_start,
                                    ex.in_stop, ex.skip_stop};
                    int reg_lo[] = {regions_c[0] - in_exon,  regions_c[1] - in_intron,
                                    regions_c[2] - in_exon,  regions_c[3] - in_intron};
                    int reg_hi[] = {regions_c[0] + in_intron, regions_c[1] + in_exon,
                                    regions_c[2] + in_intron, regions_c[3] + in_exon};

                    int gpos = -1;
                    for (int r = 0; r < 4; r++) {
                        int ref_reg = fwd ? r : (3 - r);
                        int delta = map_pos - regions_c[ref_reg];
                        int candidate = bounds[r] + (fwd ? delta : -delta);
                        int lo = reg_lo[ref_reg];
                        int hi = reg_hi[ref_reg];
                        if (map_pos >= lo && map_pos <= hi) {
                            gpos = candidate;
                            break;
                        }
                    }
                    if (gpos < 0) continue;

                    float s = get_phylop(ex.chrom, gpos);
                    score_sum += s;
                    n_valid++;
                }

                avg_score[pi] = (n_valid > 0) ? score_sum / n_valid : 0;
                counts[pi] = n_valid;
            }

            for (int pi = 0; pi < n_pos; pi++)
                out << tet << "\t" << map_positions[pi] << "\t" << cat
                    << "\t" << avg_score[pi] << "\t" << counts[pi] << "\n";
        }

        done++;
        cerr << "PROGRESS " << done << " " << total << "\n";
    }

    cerr << "Conservation profiling complete: " << output_tsv << "\n";
    return 0;
}
