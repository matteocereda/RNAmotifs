// RNAmotifs - Core implementation
// Copyright (C) 2014-2026 Matteo Cereda
// SPDX-License-Identifier: GPL-2.0-or-later

#include "rnamotifs_core.h"
#include <sys/stat.h>
#include <dirent.h>
#include <filesystem>

namespace fs = std::filesystem;

// ===========================================================================
// Config
// ===========================================================================

Config::Config(int argc, char *argv[]) {
    string config_path;
    for (int i = 1; i < argc; ++i) {
        string arg(argv[i]);
        if ((arg == "-c" || arg == "--config") && i + 1 < argc) {
            config_path = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            cout << "Usage: " << argv[0] << " -c <config_file>\n";
            exit(0);
        }
    }
    if (config_path.empty()) {
        cerr << "Error: configuration file required (-c <file>)\n";
        exit(1);
    }

    ifstream in(config_path);
    if (!in.is_open()) {
        cerr << "Error: cannot open config file: " << config_path << "\n";
        exit(1);
    }

    string line;
    int idx = 0;
    while (getline(in, line)) {
        size_t eq = line.find('=');
        if (eq == string::npos) {
            cerr << "Warning: malformed config line (no '='): " << line << "\n";
            continue;
        }
        string val = line.substr(eq + 1);
        switch (idx) {
            case 0: mouse = (val == "Mouse"); break;
            case 1: splicing_file = val;       break;
            case 2: tetramer_folder = val;     break;
            case 3: results_folder = val;      break;
            case 4: search = val;              break;
            case 5: enrichment_window = stoi(val); break;
            case 6: in_exon = stoi(val); break;
            case 7: in_intron = stoi(val); break;
            case 8: event_type = val; break;
            default:
                cerr << "Warning: extra config line ignored: " << line << "\n";
                break;
        }
        ++idx;
    }

    if (idx < 8) {
        cerr << "Error: config file incomplete (expected at least 8 lines, got "
             << idx << ")\n";
        exit(1);
    }
    // event_type defaults to "SE" if not in config (line 9 is optional)

    compute_regions();

    if (search == "R") {
        tetramer_folder += "/r/";
        results_folder  += "/r/";
    } else {
        tetramer_folder += "/nr/";
        results_folder  += "/nr/";
    }

    fs::create_directories(results_folder);
}

void Config::print(ostream &out) const {
    out << "  organism         : " << (mouse ? "Mouse" : "Human") << "\n"
        << "  search mode      : " << search << "\n"
        << "  splicing file    : " << splicing_file << "\n"
        << "  tetramer folder  : " << tetramer_folder << "\n"
        << "  results folder   : " << results_folder << "\n"
        << "  enrichment window: " << enrichment_window << "\n"
        << "  in-exon span     : " << in_exon << "\n"
        << "  in-intron span   : " << in_intron << "\n\n";
}

// ===========================================================================
// List .bed files in a directory (replaces system("wc -l ..."))
// ===========================================================================

static vector<string> list_bed_files(const string &dir) {
    vector<string> files;
    try {
        for (const auto &entry : fs::directory_iterator(dir)) {
            if (entry.is_regular_file() &&
                entry.path().extension() == ".bed") {
                files.push_back(entry.path().string());
            }
        }
    } catch (const fs::filesystem_error &e) {
        cerr << "Error listing directory: " << e.what() << "\n";
    }
    sort(files.begin(), files.end());
    return files;
}

// ===========================================================================
// BedRecord
// ===========================================================================

BedRecord::BedRecord(const string &line) {
    istringstream ss(line);
    string tok;
    vector<string> cols;
    while (getline(ss, tok, '\t')) cols.push_back(tok);

    if (cols.size() >= 4) {
        strncpy(chr, cols[0].c_str(), sizeof(chr) - 1);
        chrom_start = stoul(cols[1]);
        chrom_end   = stoul(cols[2]);
        strand_num  = stoi(cols[3]);
        forward     = (strand_num > 0);
        count       = 1;
    }
}

void BedRecord::map_to_splicing_map(
        const Config &cfg, bool fwd, int region,
        unsigned int reg_start, unsigned int reg_stop,
        unsigned int boundary,
        unsigned int intron_span, unsigned int exon_span) {

    int delta_start = abs((int)boundary - (int)chrom_start);
    int delta_stop  = abs((int)chrom_end - (int)boundary);

    auto map_pos = [&](unsigned int pos, int delta, bool fwd, int region) -> unsigned int {
        if (fwd) {
            return (pos <= boundary)
                ? cfg.regions[region] - delta
                : cfg.regions[region] + delta;
        } else {
            return (pos <= boundary)
                ? cfg.regions[3 - region] + delta
                : cfg.regions[3 - region] - delta;
        }
    };

    auto clamp_start = [&](bool fwd, int region) -> unsigned int {
        if (fwd)
            return (region == 0 || region == 2)
                ? cfg.regions[region] - exon_span
                : cfg.regions[region] - intron_span;
        else
            return (region == 0 || region == 2)
                ? cfg.regions[3 - region] + exon_span
                : cfg.regions[3 - region] + intron_span;
    };

    auto clamp_end = [&](bool fwd, int region) -> unsigned int {
        if (fwd)
            return (region == 0 || region == 2)
                ? cfg.regions[region] + intron_span
                : cfg.regions[region] + exon_span;
        else
            return (region == 0 || region == 2)
                ? cfg.regions[3 - region] - intron_span
                : cfg.regions[3 - region] - exon_span;
    };

    // Fully contained
    if (reg_start <= chrom_start && chrom_end <= reg_stop) {
        map_start = map_pos(chrom_start, delta_start, fwd, region);
        map_end   = map_pos(chrom_end,   delta_stop,  fwd, region);
    }
    // Left contained, right overflow
    else if (reg_start <= chrom_start && chrom_end > reg_stop) {
        map_start = map_pos(chrom_start, delta_start, fwd, region);
        map_end   = clamp_end(fwd, region);
    }
    // Left overflow, right contained
    else if (chrom_start < reg_start && chrom_end <= reg_stop) {
        map_start = clamp_start(fwd, region);
        map_end   = map_pos(chrom_end, delta_stop, fwd, region);
    }
    // Both overflow
    else if (chrom_start < reg_start && chrom_end > reg_stop) {
        map_start = clamp_start(fwd, region);
        map_end   = clamp_end(fwd, region);
    }
}

void BedRecord::map_midpoint(
        const Config &cfg, int midpoint, bool fwd, int region,
        unsigned int reg_start, unsigned int reg_stop,
        unsigned int boundary,
        unsigned int intron_span, unsigned int exon_span) {

    int delta = abs((int)boundary - midpoint);
    int pos;
    if (fwd) {
        pos = ((unsigned int)midpoint <= boundary)
            ? cfg.regions[region] - delta
            : cfg.regions[region] + delta;
    } else {
        pos = ((unsigned int)midpoint <= boundary)
            ? cfg.regions[3 - region] + delta
            : cfg.regions[3 - region] - delta;
    }
    map_start = pos;
    map_end   = pos;
}

bool sort_by_row_id(const BedRecord &a, const BedRecord &b) {
    return a.row_id < b.row_id;
}

bool sort_by_map_cat(const BedRecord &a, const BedRecord &b) {
    auto weight = [](const BedRecord &r) -> int {
        int w = 0;
        if (r.category == 1)       w = 10;
        else if (r.category == 0)  w = 1;
        else if (r.category == -1) w = -10;
        return w * (int)r.map_end;
    };
    int va = weight(a), vb = weight(b);
    if (va < 0 && vb < 0) return va > vb;
    return va < vb;
}

// ===========================================================================
// RegionWindow
// ===========================================================================

RegionWindow::RegionWindow(const BedRecord &b) : BedRecord(b) {}

void RegionWindow::accumulate(const BedRecord &b, const Config &cfg) {
    unsigned int r1s = cfg.regions[0] - cfg.in_exon,
                 r1e = cfg.regions[0] + cfg.in_intron,
                 r2s = cfg.regions[1] - cfg.in_intron,
                 r2e = cfg.regions[1] + cfg.in_exon,
                 r3s = cfg.regions[2] - cfg.in_exon,
                 r3e = cfg.regions[2] + cfg.in_intron,
                 r4s = cfg.regions[3] - cfg.in_intron,
                 r4e = cfg.regions[3] + cfg.in_exon;

    if (r1s <= b.map_end && b.map_end <= r1e)      sum1 += b.count;
    else if (r2s <= b.map_end && b.map_end <= r2e)  sum2 += b.count;
    else if (r3s <= b.map_end && b.map_end <= r3e)  sum3 += b.count;
    else if (r4s <= b.map_end && b.map_end <= r4e)  sum4 += b.count;
}

// ===========================================================================
// Tetramer: positional tetramer occurrences along the RNA splicing map
// ===========================================================================

static int process_tetramer_file(
        Config &cfg,
        const string &sp_fname,
        const string &bed_fname,
        const string &out_fname,
        const string &tet_name,
        ofstream &stats) {
    try {
        ifstream in(sp_fname), fbed(bed_fname);
        ofstream fout(out_fname);

        if (!in.is_open() || !fbed.is_open()) {
            cerr << "Cannot open: " << sp_fname << " or " << bed_fname << "\n";
            return 1;
        }

        // Read BED
        vector<BedRecord> bed_lines;
        string line;
        while (getline(fbed, line))
            bed_lines.emplace_back(line);
        fbed.close();

        // Index BED records by (chromosome, strand) for O(log N) lookup
        typedef pair<string, bool> ChrStrandKey;
        map<ChrStrandKey, vector<BedRecord*>> bed_index;
        for (auto &bl : bed_lines) {
            bed_index[ChrStrandKey(string(bl.chr), bl.forward)].push_back(&bl);
        }
        // Sort each group by chrom_start for binary search
        for (auto &kv : bed_index) {
            sort(kv.second.begin(), kv.second.end(),
                 [](const BedRecord *a, const BedRecord *b) {
                     return a->chrom_start < b->chrom_start;
                 });
        }

        vector<BedRecord> bed_results;
        unsigned int CE1 = 0, CE0 = 0, CEm1 = 0;
        unsigned int rCE1 = 0, rCE0 = 0, rCEm1 = 0;

        // Process splicing events
        string f;
        while (getline(in, f)) {
            vector<string> pars = split_string(f, ';');

            string       chr        = pars[2];
            bool         fwd        = (pars[3] == "+");
            unsigned int exon_id    = stoul(pars[0]);
            unsigned int row_id     = stoul(pars[1]);
            unsigned int skip_start = stoul(pars[4]);
            unsigned int in_start   = stoul(pars[5]);
            unsigned int in_stop    = stoul(pars[6]);
            unsigned int skip_stop  = stoul(pars[7]);
            double       dIRank     = stod(pars[8]);

            int cat = 2;
            if (-cfg.dIRZ <= dIRank && dIRank <= cfg.dIRZ) { CE0++;  cat = 0;  }
            else if (dIRank >= cfg.dIRO)                    { CE1++;  cat = 1;  }
            else if (dIRank <= -cfg.dIRO)                   { CEm1++; cat = -1; }

            // Auto-detect event type from optional 10th field
            string evt = cfg.event_type;
            if (pars.size() > 9 && !pars[9].empty()) evt = pars[9];

            unsigned int ie = cfg.in_exon;
            unsigned int ii = cfg.in_intron;

            unsigned int r1s, r1e, r2s, r2e, r3s, r3e, r4s, r4e;
            unsigned int intron_left, exon_span, intron_right;

            if (evt == "RI") {
                // ── Intron Retention regions ──────────────────────────
                // v5=upstreamES, v6=upstreamEE (5'SS), v7=downstreamES (3'SS), v8=downstreamEE
                unsigned int upstream_exon_len   = in_start - skip_start;
                unsigned int retained_intron_len = in_stop  - in_start;
                unsigned int downstream_exon_len = skip_stop - in_stop;
                unsigned int half_intron = retained_intron_len / 2;

                // in_exon → extent into flanking exons (R1, R4)
                unsigned int ie_left  = (upstream_exon_len < ie * 2)   ? upstream_exon_len / 2   : ie;
                unsigned int ie_right = (downstream_exon_len < ie * 2) ? downstream_exon_len / 2 : ie;
                // in_intron → extent into retained intron (R2, R3), capped at half-intron
                unsigned int ii_actual = (half_intron < ii) ? half_intron : ii;

                // R1: upstream exon near 5'SS
                r1s = in_start - ie_left;
                r1e = in_start;
                // R2: retained intron from 5'SS
                r2s = in_start;
                r2e = in_start + ii_actual;
                // R3: retained intron toward 3'SS
                r3s = in_stop - ii_actual;
                r3e = in_stop;
                // R4: downstream exon near 3'SS
                r4s = in_stop;
                r4e = in_stop + ie_right;

                // Clamp R2/R3 at intron midpoint if they overlap
                unsigned int intron_mid = in_start + half_intron;
                if (r2e > intron_mid) r2e = intron_mid;
                if (r3s < intron_mid) r3s = intron_mid;

                intron_left  = ie_left;       // "intron" side of R1 boundary = exon extent
                exon_span    = r2e - in_start; // "exon" side of R2 boundary = intron extent
                intron_right = ie_right;       // "intron" side of R3 boundary = exon extent
            } else {
                // ── Cassette Exon (SE) regions (original logic) ──────
                unsigned int left_intron_len  = in_start - skip_start;
                unsigned int exon_len         = in_stop  - in_start;
                unsigned int right_intron_len = skip_stop - in_stop;

                unsigned int ie_actual = (exon_len < ie * 2) ? exon_len / 2 : ie;
                unsigned int ii_left   = (left_intron_len < ii * 2)  ? left_intron_len / 2  : ii;
                unsigned int ii_right  = (right_intron_len < ii * 2) ? right_intron_len / 2 : ii;

                r1s = skip_start - ie_actual;
                r1e = skip_start + ii_left;
                r2s = in_start   - ii_left;
                r2e = in_start   + ie_actual;
                r3s = in_stop    - ie_actual;
                r3e = in_stop    + ii_right;
                r4s = skip_stop  - ii_right;
                r4e = skip_stop  + ie_actual;

                unsigned int left_mid  = skip_start + (in_start - skip_start) / 2;
                unsigned int exon_mid  = in_start   + (in_stop  - in_start)  / 2;
                unsigned int right_mid = in_stop    + (skip_stop - in_stop)  / 2;

                if (left_mid  < r1e) r1e = left_mid;
                if (left_mid  > r2s) r2s = left_mid;
                if (exon_mid  < r2e) r2e = exon_mid;
                if (exon_mid  > r3s) r3s = exon_mid;
                if (right_mid < r3e) r3e = right_mid;
                if (right_mid > r4s) r4s = right_mid;

                intron_left  = r1e - skip_start;
                exon_span    = r2e - in_start;
                intron_right = r3e - in_stop;
            }

            bool has_hit = false;
            // Indexed lookup: find BED records on same chr+strand
            auto it = bed_index.find(ChrStrandKey(chr, fwd));
            if (it != bed_index.end()) {
                const auto &group = it->second;
                // Binary search: find first record with chrom_start <= r4e
                // upper_bound finds first record with chrom_start > r4e;
                // everything before that has chrom_start <= r4e
                BedRecord sentinel;
                sentinel.chrom_start = r4e;
                auto ub = upper_bound(group.begin(), group.end(), &sentinel,
                    [](const BedRecord *a, const BedRecord *b) {
                        return a->chrom_start < b->chrom_start;
                    });
                // Iterate from beginning of group up to ub, skip if chrom_end < r1s
                for (auto gi = group.begin(); gi != ub; ++gi) {
                    const BedRecord &bl = **gi;
                    if (bl.chrom_end < r1s)
                        continue;

                    BedRecord bed(bl);
                    bed.row_id   = row_id;
                    bed.exon_id  = exon_id;
                    bed.category = cat;

                    if (evt == "RI") {
                        // IR: R1=upstream exon, R2=5' intron, R3=3' intron, R4=downstream exon
                        // Boundaries are at in_start (5'SS) and in_stop (3'SS)
                        if (bl.chrom_start <= r1e && bl.chrom_end >= r1s) {
                            bed.map_to_splicing_map(cfg, fwd, 0, r1s, r1e, in_start, intron_left, intron_left);
                            bed_results.push_back(bed);
                        }
                        else if (bl.chrom_start <= r2e && bl.chrom_end >= r2s) {
                            bed.map_to_splicing_map(cfg, fwd, 1, r2s, r2e, in_start, intron_left, exon_span);
                            bed_results.push_back(bed);
                            has_hit = true;
                        }
                        else if (bl.chrom_start <= r3e && bl.chrom_end >= r3s) {
                            bed.map_to_splicing_map(cfg, fwd, 2, r3s, r3e, in_stop, exon_span, intron_right);
                            bed_results.push_back(bed);
                            has_hit = true;
                        }
                        else if (bl.chrom_start <= r4e && bl.chrom_end >= r4s) {
                            bed.map_to_splicing_map(cfg, fwd, 3, r4s, r4e, in_stop, intron_right, intron_right);
                            bed_results.push_back(bed);
                        }
                    } else {
                        // SE: original cassette exon mapping
                        if (bl.chrom_start <= r1e && bl.chrom_end >= r1s) {
                            bed.map_to_splicing_map(cfg, fwd, 0, r1s, r1e, skip_start, intron_left, cfg.in_exon);
                            bed_results.push_back(bed);
                        }
                        else if (bl.chrom_start <= r2e && bl.chrom_end >= r2s) {
                            bed.map_to_splicing_map(cfg, fwd, 1, r2s, r2e, in_start, intron_left, exon_span);
                            bed_results.push_back(bed);
                            if (bl.chrom_end >= r2s + 100 && bl.chrom_start <= r2e - 15) has_hit = true;
                        }
                        else if (bl.chrom_start <= r3e && bl.chrom_end >= r3s) {
                            bed.map_to_splicing_map(cfg, fwd, 2, r3s, r3e, in_stop, intron_right, exon_span);
                            bed_results.push_back(bed);
                            if (bl.chrom_end >= r3s + 15 && bl.chrom_start <= r3e - 130) has_hit = true;
                        }
                        else if (bl.chrom_start <= r4e && bl.chrom_end >= r4s) {
                            bed.map_to_splicing_map(cfg, fwd, 3, r4s, r4e, skip_stop, intron_right, cfg.in_exon);
                            bed_results.push_back(bed);
                        }
                    }
                }
            }
            if (has_hit) {
                if (cat == 0)       rCE0++;
                else if (cat == 1)  rCE1++;
                else if (cat == -1) rCEm1++;
            }
        }

        stats << tet_name << "\t" << rCE1 << "\t" << rCEm1 << "\t" << rCE0 << "\n";

        if (bed_results.empty()) {
            cout << "  " << tet_name << ": no results\n";
            return 0;
        }

        // Expand: each BED span gets one record per covered base
        vector<BedRecord> expanded;
        expanded.reserve(bed_results.size() * 4);
        for (auto &br : bed_results) {
            unsigned int gen = br.chrom_end;
            unsigned int rna = br.map_end;
            int n = abs(br.strand_num);
            for (int i = 0; i < n; ++i) {
                BedRecord bb(br);
                bb.chrom_end = gen - i;
                bb.map_end   = (br.strand_num > 0) ? rna - i : rna + i;
                expanded.push_back(bb);
            }
        }

        sort(expanded.begin(), expanded.end(), sort_by_map_cat);

        // Collapse to unique (map_end, category) counts
        vector<BedRecord> ranked;
        ranked.push_back(expanded[0]);
        ranked.back().count = 1;
        for (size_t i = 1; i < expanded.size(); ++i) {
            if (expanded[i].map_end == expanded[i - 1].map_end &&
                expanded[i].category == expanded[i - 1].category) {
                ranked.back().count++;
            } else {
                expanded[i].count = 1;
                ranked.push_back(expanded[i]);
            }
        }

        for (auto &r : ranked)
            fout << r.map_end << "\t" << r.category << "\t" << r.count << "\n";

    } catch (exception &e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}

int run_tetramer(Config &cfg) {
    vector<string> bed_files = list_bed_files(cfg.tetramer_folder);
    if (bed_files.empty()) {
        cerr << "No tetramer .bed files in " << cfg.tetramer_folder << "\n";
        return 1;
    }

    string stats_path = cfg.results_folder + "STATS.txt";
    string flist_path = cfg.results_folder + "filelist.txt";
    ofstream stats(stats_path), flist(flist_path);

    for (const auto &bed_path : bed_files) {
        string fname = fs::path(bed_path).filename().string();
        string tet = fname.substr(0, fname.size() - 4); // strip .bed
        string out = cfg.results_folder + fname;
        cout << "  tetramer " << tet << "\n";
        process_tetramer_file(cfg, cfg.splicing_file, bed_path, out, tet, stats);
        flist << fname << "\n";
    }
    return 0;
}

// ===========================================================================
// Counting: per-region hit presence/absence for each exon
// ===========================================================================

static int count_per_regions(
        Config &cfg,
        const string &sp_fname,
        const string &bed_fname,
        const string &out_fname) {
    try {
        ifstream in(sp_fname), fbed(bed_fname);
        ofstream fout(out_fname);

        if (!in.is_open() || !fbed.is_open()) {
            cerr << "Cannot open: " << sp_fname << " or " << bed_fname << "\n";
            return 1;
        }

        vector<BedRecord> bed_lines;
        string line;
        while (getline(fbed, line))
            bed_lines.emplace_back(line);
        fbed.close();

        // Index BED records by (chromosome, strand) for O(log N) lookup
        typedef pair<string, bool> ChrStrandKey;
        map<ChrStrandKey, vector<BedRecord*>> bed_index;
        for (auto &bl : bed_lines) {
            bed_index[ChrStrandKey(string(bl.chr), bl.forward)].push_back(&bl);
        }
        for (auto &kv : bed_index) {
            sort(kv.second.begin(), kv.second.end(),
                 [](const BedRecord *a, const BedRecord *b) {
                     return a->chrom_start < b->chrom_start;
                 });
        }

        fout << "myRID\trowID\ttype\thits_region1\thits_region2\thits_region3\n";

        string f;
        while (getline(in, f)) {
            vector<string> pars = split_string(f, ';');

            string       chr      = pars[2];
            bool         fwd      = (pars[3] == "+");
            unsigned int exon_id  = stoul(pars[0]);
            unsigned int row_id   = stoul(pars[1]);
            unsigned int in_start = stoul(pars[5]);
            unsigned int in_stop  = stoul(pars[6]);
            double       dIRank   = stod(pars[8]);

            int cat = 2;
            if (-cfg.dIRZ <= dIRank && dIRank <= cfg.dIRZ)  cat = 0;
            else if (dIRank >= cfg.dIRO)                     cat = 1;
            else if (dIRank <= -cfg.dIRO)                    cat = -1;

            // Auto-detect event type from optional 10th field
            string evt_c = cfg.event_type;
            if (pars.size() > 9 && !pars[9].empty()) evt_c = pars[9];

            unsigned int exon_len_c = in_stop - in_start;

            bool h1 = false, h2 = false, h3 = false;
            unsigned int lr1s, lr1e, lr3s, lr3e;
            unsigned int exon_r1e, exon_r2s;

            if (evt_c == "RI") {
                // IR counting: R1=upstream intron near 5'SS, R2=first half intron, R3=second half intron
                unsigned int half_intron = exon_len_c / 2;
                unsigned int ii_c = (half_intron < cfg.in_intron) ? half_intron : cfg.in_intron;

                // R1: upstream of 5'SS (intronic flank)
                lr1s = in_start - cfg.enrichment_window - 5;
                lr1e = in_start - 5;
                // R2: first half of retained intron
                exon_r1e = in_start + ii_c;
                // R3: second half of retained intron
                exon_r2s = in_stop - ii_c;
                // R3 intronic flank: downstream of 3'SS
                lr3s = in_stop + 10;
                lr3e = in_stop + 10 + cfg.enrichment_window;

                unsigned int intron_mid = in_start + half_intron;
                if (intron_mid < exon_r1e) exon_r1e = intron_mid;
                if (intron_mid > exon_r2s) exon_r2s = intron_mid;

                if (!fwd) {
                    lr1s = in_stop + 5;
                    lr1e = in_stop + 5 + cfg.enrichment_window;
                    lr3s = in_start - 10 - cfg.enrichment_window;
                    lr3e = in_start - 10;
                }
            } else {
                // SE counting (original)
                unsigned int ie_c = (exon_len_c < cfg.in_exon * 2) ? exon_len_c / 2 : cfg.in_exon;

                lr1s = in_start - 5 - cfg.enrichment_window;
                lr1e = in_start - 5;
                exon_r1e = in_start + ie_c;
                exon_r2s = in_stop  - ie_c;
                unsigned int exon_mid = in_start + exon_len_c / 2;
                lr3s = in_stop + 10;
                lr3e = in_stop + 10 + cfg.enrichment_window;

                if (exon_mid < exon_r1e) exon_r1e = exon_mid;
                if (exon_mid > exon_r2s) exon_r2s = exon_mid;

                if (!fwd) {
                    lr1s = in_stop + 5;
                    lr1e = in_stop + 5 + cfg.enrichment_window;
                    lr3s = in_start - 10 - cfg.enrichment_window;
                    lr3e = in_start - 10;
                }
            }

            // Determine the overall genomic range that could produce any hit
            unsigned int min_start = lr1s;
            if (in_start < min_start)  min_start = in_start;
            if (exon_r2s < min_start)  min_start = exon_r2s;
            if (lr3s < min_start)      min_start = lr3s;

            unsigned int max_end = lr1e;
            if (exon_r1e > max_end)  max_end = exon_r1e;
            if (in_stop > max_end)   max_end = in_stop;
            if (lr3e > max_end)      max_end = lr3e;

            // Indexed lookup: find BED records on same chr+strand
            auto it = bed_index.find(ChrStrandKey(chr, fwd));
            if (it != bed_index.end()) {
                const auto &group = it->second;
                // Find first record with chrom_start > max_end
                BedRecord sentinel;
                sentinel.chrom_start = max_end;
                auto ub = upper_bound(group.begin(), group.end(), &sentinel,
                    [](const BedRecord *a, const BedRecord *b) {
                        return a->chrom_start < b->chrom_start;
                    });
                for (auto gi = group.begin(); gi != ub; ++gi) {
                    const BedRecord &bl = **gi;
                    if (bl.chrom_end < min_start)
                        continue;

                    if (bl.chrom_start <= lr1e && bl.chrom_end >= lr1s)               h1 = true;
                    else if (bl.chrom_start <= exon_r1e && bl.chrom_end >= in_start)   h2 = true;
                    else if (bl.chrom_start <= in_stop && bl.chrom_end >= exon_r2s)    h2 = true;
                    else if (bl.chrom_start <= lr3e && bl.chrom_end >= lr3s)            h3 = true;

                    if (h1 && h2 && h3) break;  // Early exit: all flags set
                }
            }

            fout << exon_id << "\t" << row_id << "\t" << cat << "\t"
                 << h1 << "\t" << h2 << "\t" << h3 << "\n";
        }
    } catch (exception &e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}

int run_counting(Config &cfg) {
    vector<string> bed_files = list_bed_files(cfg.tetramer_folder);
    if (bed_files.empty()) {
        cerr << "No tetramer .bed files in " << cfg.tetramer_folder << "\n";
        return 1;
    }

    string flist_path = cfg.results_folder + "filelist_count.tsv";
    ofstream flist(flist_path);

    for (const auto &bed_path : bed_files) {
        string fname = fs::path(bed_path).filename().string();
        string tet = fname.substr(0, fname.size() - 4); // strip .bed
        string out = cfg.results_folder + tet + "_region_count.tsv";
        cout << "  counting " << tet << "\n";
        count_per_regions(cfg, cfg.splicing_file, bed_path, out);
        flist << tet << "_region_count.tsv\n";
    }
    return 0;
}
