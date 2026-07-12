cmake_minimum_required(VERSION 3.10)
project(RNAmotifs VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# --- OpenMP ---
find_package(OpenMP REQUIRED)

# --- Tetramer mapping ---
add_executable(rnamotifs_tetramer src/cpp/tetramer.cpp src/cpp/rnamotifs_core.cpp)
target_include_directories(rnamotifs_tetramer PRIVATE src/cpp)
target_link_libraries(rnamotifs_tetramer stdc++fs)

# --- Region counting ---
add_executable(rnamotifs_counting src/cpp/counting.cpp src/cpp/rnamotifs_core.cpp)
target_include_directories(rnamotifs_counting PRIVATE src/cpp)
target_link_libraries(rnamotifs_counting stdc++fs)

# --- Bootstrap FDR (standalone, OpenMP) ---
add_executable(rnamotifs_bootstrap src/cpp/bootstrap_fdr.cpp)
target_link_libraries(rnamotifs_bootstrap OpenMP::OpenMP_CXX)

# --- Motif search (OpenMP) ---
add_executable(rnamotifs_search src/cpp/rnamotifs_search.cpp)
target_include_directories(rnamotifs_search PRIVATE src/cpp)
target_link_libraries(rnamotifs_search stdc++fs OpenMP::OpenMP_CXX)

# --- Tetramer selection (standalone) ---
add_executable(rnamotifs_selection src/cpp/rnamotifs_selection.cpp)
target_link_libraries(rnamotifs_selection stdc++fs)

# --- MaRs association scoring (OpenMP) ---
add_executable(rnamotifs_mars_score src/cpp/rnamotifs_mars_score.cpp)
target_link_libraries(rnamotifs_mars_score stdc++fs OpenMP::OpenMP_CXX)

# --- Conservation profiling (OpenMP) ---
add_executable(rnamotifs_conservation src/cpp/rnamotifs_conservation.cpp)
target_include_directories(rnamotifs_conservation PRIVATE src/cpp)
target_link_libraries(rnamotifs_conservation stdc++fs OpenMP::OpenMP_CXX)

# --- RNA structure profiling (ViennaRNA + OpenMP) ---
find_library(RNA_LIB RNA PATHS /usr/lib /usr/local/lib)
find_path(RNA_INCLUDE ViennaRNA/fold_compound.h PATHS /usr/include /usr/local/include)

if(RNA_LIB AND RNA_INCLUDE)
    message(STATUS "ViennaRNA found: ${RNA_LIB}")
    add_executable(rnamotifs_structure src/cpp/rnamotifs_structure.cpp)
    target_include_directories(rnamotifs_structure PRIVATE ${RNA_INCLUDE} src/cpp)
    target_link_libraries(rnamotifs_structure
        ${RNA_LIB} stdc++fs OpenMP::OpenMP_CXX m)
    install(TARGETS rnamotifs_tetramer rnamotifs_counting rnamotifs_bootstrap
            rnamotifs_search rnamotifs_selection rnamotifs_conservation rnamotifs_structure
            DESTINATION bin)
else()
    message(STATUS "ViennaRNA not found — rnamotifs_structure will not be built")
    message(STATUS "  Install with: sudo apt-get install vienna-rna libvienna-rna-dev")
    install(TARGETS rnamotifs_tetramer rnamotifs_counting rnamotifs_bootstrap
            rnamotifs_search rnamotifs_selection rnamotifs_conservation DESTINATION bin)
endif()
