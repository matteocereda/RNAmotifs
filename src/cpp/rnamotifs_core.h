// RNAmotifs - Core definitions
// Copyright (C) 2014-2026 Matteo Cereda
// Based on GeCo++ v0.2 (Cereda & Pozzoli, Bioinformatics 2011)
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef RNAMOTIFS_CORE_H
#define RNAMOTIFS_CORE_H

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

// Split a string by delimiter into a vector (replaces gString::split)
inline vector<string> split_string(const string &s, char delim) {
    vector<string> tokens;
    istringstream iss(s);
    string tok;
    while (getline(iss, tok, delim))
        tokens.push_back(tok);
    return tokens;
}

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

struct Config {
    bool   mouse             = false;
    string search;                         // "N" or "R"
    string results_folder;
    string tetramer_folder;
    string splicing_file;
    double dIRZ              = 0.1;
    double dIRO              = 1.0;
    unsigned int regions[4]  = {};
    unsigned int in_exon     = 30;
    unsigned int in_intron   = 300;
    void compute_regions() {
        unsigned int mx = (in_exon > in_intron) ? in_exon : in_intron;
        unsigned int spacing = 2 * mx + 10;
        for (int i = 0; i < 4; ++i)
            regions[i] = spacing * (i + 1);
    }
    unsigned int enrichment_window = 30;

    Config() { compute_regions(); }
    Config(int argc, char *argv[]);
    void print(ostream &out) const;
};

// ---------------------------------------------------------------------------
// BED record mapped to the RNA splicing map
// ---------------------------------------------------------------------------

struct BedRecord {
    char         chr[10]      = {};
    unsigned int chrom_start  = 0;
    unsigned int chrom_end    = 0;
    unsigned int map_start    = 0;
    unsigned int map_end      = 0;
    int          strand_num   = 0;
    bool         forward      = true;
    unsigned int row_id       = 0;
    unsigned int exon_id      = 0;
    int          category     = 2;   // 1=enh, 0=ctrl, -1=sil, 2=skip
    int          count        = 1;

    BedRecord() = default;
    explicit BedRecord(const string &line);

    void map_to_splicing_map(const Config &cfg, bool fwd, int region,
                             unsigned int reg_start, unsigned int reg_stop,
                             unsigned int boundary,
                             unsigned int intron_span, unsigned int exon_span);

    void map_midpoint(const Config &cfg, int midpoint, bool fwd, int region,
                      unsigned int reg_start, unsigned int reg_stop,
                      unsigned int boundary,
                      unsigned int intron_span, unsigned int exon_span);
};

// Sorting helpers
bool sort_by_row_id(const BedRecord &a, const BedRecord &b);
bool sort_by_map_cat(const BedRecord &a, const BedRecord &b);

// ---------------------------------------------------------------------------
// Region-count window (per-exon tetramer hits in R1/R2/R3/R4)
// ---------------------------------------------------------------------------

struct RegionWindow : public BedRecord {
    unsigned int sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;

    RegionWindow() = default;
    explicit RegionWindow(const BedRecord &b);
    void accumulate(const BedRecord &b, const Config &cfg);
};

// ---------------------------------------------------------------------------
// Pipeline entry points
// ---------------------------------------------------------------------------

int run_tetramer(Config &cfg);
int run_counting(Config &cfg);

#endif // RNAMOTIFS_CORE_H
