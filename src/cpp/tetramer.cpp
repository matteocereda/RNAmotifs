// RNAmotifs - Tetramer: positional occurrences along the RNA splicing map
// Copyright (C) 2014-2026 Matteo Cereda
// SPDX-License-Identifier: GPL-2.0-or-later

#include "rnamotifs_core.h"

int main(int argc, char *argv[]) {
    try {
        Config cfg(argc, argv);
        cout << "[tetramer] Configuration:\n";
        cfg.print(cout);
        return run_tetramer(cfg);
    } catch (exception &e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
