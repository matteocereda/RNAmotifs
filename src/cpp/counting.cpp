// RNAmotifs - Counting: per-region tetramer hit presence for each exon
// Copyright (C) 2014-2026 Matteo Cereda
// SPDX-License-Identifier: GPL-2.0-or-later

#include "rnamotifs_core.h"

int main(int argc, char *argv[]) {
    try {
        Config cfg(argc, argv);
        cout << "[counting] Configuration:\n";
        cfg.print(cout);
        return run_counting(cfg);
    } catch (exception &e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
