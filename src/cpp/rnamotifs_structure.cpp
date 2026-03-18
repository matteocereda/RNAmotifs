// RNAmotifs - RNA secondary structure profiling
// Copyright (C) 2014-2026 Matteo Cereda
// SPDX-License-Identifier: GPL-2.0-or-later
//
// For each enriched tetramer, compute position-specific RNA single-stranded
// scores (from ViennaRNA partition function) averaged across all exons with
// a tetramer occurrence. Produces a TSV profile for RNA map visualisation.

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

// Single-strandedness is computed via the constrained partition function
// ratio method, calling ViennaRNA directly.

extern "C" {
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/constraints/hard.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/utils/basic.h>
}

using namespace std;
namespace fs = std::filesystem;

// ===========================================================================
// Chromosome cache (same as rnamotifs_search.cpp)
// ===========================================================================

static unordered_map<string, string> chrom_cache;
static mutex chrom_mutex;

static void preload_chromosome(const string &genome_dir, const string &genome,
                               const string &chrom) {
    string path = genome_dir + "/genomes/" + genome + "/" + chrom + ".string";
    ifstream fh(path, ios::binary);
    if (!fh.is_open()) {
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
// RNA single-stranded score via ViennaRNA constrained partition function
// ===========================================================================

// Compute unconstrained partition function free energy
static double compute_unconstrained_pf(const char *sequence,
                                        const vrna_md_t *md) {
    vrna_fold_compound_t *fc = vrna_fold_compound(sequence, md,
                                                   VRNA_OPTION_PF);
    double energy = vrna_pf(fc, NULL);
    vrna_fold_compound_free(fc);
    return energy;
}

// Compute constrained partition function (central nucleotide forced unpaired)
static double compute_constrained_pf(const char *sequence,
                                      int unpaired_pos,
                                      const vrna_md_t *md) {
    vrna_fold_compound_t *fc = vrna_fold_compound(sequence, md,
                                                   VRNA_OPTION_PF);
    // Force the nucleotide at unpaired_pos (0-based) to be unpaired
    // vrna_hc_add_up uses 1-based indexing
    vrna_hc_add_up(fc, unpaired_pos + 1, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);
    double energy = vrna_pf(fc, NULL);
    vrna_fold_compound_free(fc);
    return energy;
}

static double single_stranded_score(const string &seq) {
    // Compute the probability that the central nucleotide is unpaired
    // using the constrained partition function ratio method.
    // P(unpaired) = exp((G_unconstrained - G_constrained) / kT)
    int n = (int)seq.size();
    if (n == 0) return 0.5;

    // Convert to uppercase RNA, replace non-ACGU with N
    string clean(n, 'N');
    for (int i = 0; i < n; i++) {
        char c = toupper(seq[i]);
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'U')
            clean[i] = (c == 'T') ? 'U' : c;  // DNA->RNA
    }

    int center = n / 2;

    // Set up ViennaRNA model details; disable base-pair probability
    // computation since we only need partition function energies
    vrna_md_t md;
    vrna_md_set_default(&md);
    md.compute_bpp = 0;

    // Get kT from ViennaRNA parameters
    vrna_exp_param_t *params = vrna_exp_params(&md);
    double KT = params->kT / 1000.0;
    free(params);

    // Unconstrained partition function
    double g_unconstrained = compute_unconstrained_pf(clean.c_str(), &md);

    // Constrained: force center nucleotide to be unpaired
    double g_constrained = compute_constrained_pf(clean.c_str(), center, &md);

    // P(unpaired) = exp((G_unconstrained - G_constrained) / kT)
    double p_ss = exp((g_unconstrained - g_constrained) / KT);
    return max(0.0, min(1.0, p_ss));
}

// ===========================================================================
// Read splicing file -> exon coordinates
// ===========================================================================

struct Exon {
    string chrom;
    string strand;
    int    skip_start;
    int    in_start;
    int    in_stop;
    int    skip_stop;
    int    category;  // 1=enh, 0=ctrl, -1=sil
    int    exon_id;
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

// ===========================================================================
// Read BED file -> set of (chrom, start, end, strand) intervals
// ===========================================================================

struct BedHit {
    string chrom;
    int    start;
    int    end;
    bool   forward;
};

static vector<BedHit> read_bed(const string &path) {
    vector<BedHit> hits;
    ifstream fin(path);
    string line;
    while (getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        istringstream ss(line);
        string chrom, s1, s2, s3;
        getline(ss, chrom, '\t');
        getline(ss, s1, '\t');
        getline(ss, s2, '\t');
        getline(ss, s3, '\t');
        BedHit h;
        h.chrom   = chrom;
        h.start   = stoi(s1);
        h.end     = stoi(s2);
        h.forward = (s3.find('-') == string::npos);
        hits.push_back(h);
    }
    return hits;
}

// ===========================================================================
// Check if exon has a tetramer hit in any of its regions
// ===========================================================================

static int s_regions[4];  // set from main
static int s_in_exon, s_in_intron;

static bool exon_has_hit(const Exon &e, const vector<BedHit> &hits) {
    int r1s = e.skip_start - s_in_exon, r1e = e.skip_start + s_in_intron;
    int r2s = e.in_start - s_in_intron, r2e = e.in_start + s_in_exon;
    int r3s = e.in_stop - s_in_exon,    r3e = e.in_stop + s_in_intron;
    int r4s = e.skip_stop - s_in_intron, r4e = e.skip_stop + s_in_exon;

    bool fwd = (e.strand == "+");
    for (auto &h : hits) {
        if (h.chrom != e.chrom || h.forward != fwd) continue;
        if (h.start <= r4e && h.end >= r1s) return true;
    }
    return false;
}

// ===========================================================================
// Map genomic position to RNA splicing map coordinate
// ===========================================================================

static int map_to_rna(const Exon &e, int gpos) {
    int *regions = s_regions;
    bool fwd = (e.strand == "+");

    struct { int bound; int region; } anchors[] = {
        {e.skip_start, 0}, {e.in_start, 1},
        {e.in_stop, 2},    {e.skip_stop, 3}
    };

    // Find closest anchor
    int best_region = -1;
    int best_dist = INT_MAX;
    for (auto &a : anchors) {
        int d = abs(gpos - a.bound);
        if (d < best_dist) {
            best_dist = d;
            best_region = a.region;
        }
    }
    if (best_region < 0) return -1;

    int delta = gpos - anchors[best_region].bound;
    if (fwd)
        return regions[best_region] + delta;
    else
        return regions[3 - best_region] - delta;
}

// ===========================================================================
// Main
// ===========================================================================

int main(int argc, char *argv[]) {
    if (argc < 7) {
        cerr << "Usage: " << argv[0]
             << " <splicing_file> <genome_dir> <genome> <tetramers_dir>"
             << " <enriched_tets_file> <output_tsv> [n_cores] [window]\n"
             << "\n"
             << "  enriched_tets_file: one tetramer name per line\n"
             << "  window: folding window size (default: 31)\n";
        return 1;
    }

    string splicing_file = argv[1];
    string genome_dir    = argv[2];
    string genome        = argv[3];
    string tet_dir       = argv[4];
    string enriched_file = argv[5];
    string output_tsv    = argv[6];
    int    n_cores       = (argc >= 8) ? stoi(argv[7]) : 1;
    int    window        = (argc >= 9) ? stoi(argv[8]) : 31;
    int    in_exon       = (argc >= 10) ? stoi(argv[9]) : 30;
    int    in_intron     = (argc >= 11) ? stoi(argv[10]) : 300;
    int    spacing       = 2 * max(in_exon, in_intron) + 10;
    int    regions_c[4]  = {spacing, 2*spacing, 3*spacing, 4*spacing};
    memcpy(s_regions, regions_c, sizeof(s_regions));
    s_in_exon = in_exon; s_in_intron = in_intron;
    int    half_w        = window / 2;

#ifdef _OPENMP
    omp_set_num_threads(n_cores);
#endif

    // Read enriched tetramer names
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
        cerr << "No enriched tetramers in " << enriched_file << "\n";
        return 1;
    }
    cerr << "Enriched tetramers: " << enriched_tets.size() << "\n";

    // Read splicing file
    vector<Exon> exons = read_splicing_file(splicing_file);
    cerr << "Exons: " << exons.size() << "\n";

    // Pre-load chromosomes
    set<string> chrom_set;
    for (auto &e : exons) chrom_set.insert(e.chrom);
    for (auto &ch : chrom_set)
        preload_chromosome(genome_dir, genome, ch);
    cerr << "Chromosomes loaded: " << chrom_set.size() << "\n";

    // RNA map positions: computed from in_exon/in_intron and region centers
    vector<int> map_positions;
    for (int i = regions_c[0] - in_exon; i <= regions_c[0] + in_intron; i++) map_positions.push_back(i);
    for (int i = regions_c[1] - in_intron; i <= regions_c[1] + in_exon; i++) map_positions.push_back(i);
    for (int i = regions_c[2] - in_exon; i <= regions_c[2] + in_intron; i++) map_positions.push_back(i);
    for (int i = regions_c[3] - in_intron; i <= regions_c[3] + in_exon; i++) map_positions.push_back(i);
    int n_pos = (int)map_positions.size();

    // Output: per tetramer, per category, per position -> average SS score
    // categories: 1 (enh), -1 (sil), 0 (ctrl)
    ofstream out(output_tsv);
    out << "tetramer\tposition\tcategory\tss_score\tn_exons\n";

    int done = 0;
    int total = (int)enriched_tets.size();

    for (const string &tet : enriched_tets) {
        // Load BED hits for this tetramer (check both nr/ and r/)
        vector<BedHit> hits;
        for (const string &sub : {"nr/", "r/"}) {
            string bed_path = tet_dir + "/" + sub + tet + ".bed";
            if (fs::exists(bed_path)) {
                auto h = read_bed(bed_path);
                hits.insert(hits.end(), h.begin(), h.end());
            }
        }

        // Find exons with hits, grouped by category
        map<int, vector<int>> cat_exon_indices;  // category -> exon indices
        for (int ei = 0; ei < (int)exons.size(); ei++) {
            if (exon_has_hit(exons[ei], hits)) {
                cat_exon_indices[exons[ei].category].push_back(ei);
            }
        }

        // For each category, compute position-specific SS scores
        for (int cat : {1, 0, -1}) {
            auto it = cat_exon_indices.find(cat);
            if (it == cat_exon_indices.end()) {
                // Write zeros
                for (int p : map_positions)
                    out << tet << "\t" << p << "\t" << cat
                        << "\t0\t0\n";
                continue;
            }

            const vector<int> &eidxs = it->second;
            int n_exons = (int)eidxs.size();

            // For each map position, collect SS scores across exons
            vector<double> avg_ss(n_pos, 0.0);
            vector<int>    counts(n_pos, 0);

            #pragma omp parallel for schedule(dynamic, 1)
            for (int pi = 0; pi < n_pos; pi++) {
                int map_pos = map_positions[pi];
                double ss_sum = 0;
                int    n_valid = 0;

                for (int ei : eidxs) {
                    const Exon &ex = exons[ei];

                    // Reverse-map: RNA map position -> genomic position
                    // Find which region this map_pos belongs to and
                    // convert back to genomic coordinates
                    int gpos = -1;
                    bool fwd = (ex.strand == "+");
                    int bounds[] = {ex.skip_start, ex.in_start,
                                    ex.in_stop, ex.skip_stop};
                    // Region extents: odd regions (0,2) = [-in_exon,+in_intron], even (1,3) = [-in_intron,+in_exon]
                    int reg_lo[] = {regions_c[0] - in_exon,  regions_c[1] - in_intron,
                                    regions_c[2] - in_exon,  regions_c[3] - in_intron};
                    int reg_hi[] = {regions_c[0] + in_intron, regions_c[1] + in_exon,
                                    regions_c[2] + in_intron, regions_c[3] + in_exon};

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

                    // Extract window
                    string seq = get_subseq(ex.chrom,
                                            gpos - half_w, gpos + half_w);
                    if ((int)seq.size() != window) continue;

                    double ss = single_stranded_score(seq);
                    ss_sum += ss;
                    n_valid++;
                }

                avg_ss[pi] = (n_valid > 0) ? ss_sum / n_valid : 0;
                counts[pi] = n_valid;
            }

            // Write results
            for (int pi = 0; pi < n_pos; pi++) {
                out << tet << "\t" << map_positions[pi] << "\t" << cat
                    << "\t" << avg_ss[pi] << "\t" << counts[pi] << "\n";
            }
        }

        done++;
        cerr << "PROGRESS " << done << " " << total << "\n";
    }

    cerr << "Structure profiling complete: " << output_tsv << "\n";
    return 0;
}
