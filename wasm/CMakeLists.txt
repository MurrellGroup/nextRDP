cmake_minimum_required(VERSION 3.20)
project(rdp_core LANGUAGES CXX)

if(NOT EMSCRIPTEN)
  message(FATAL_ERROR "rdp-core is a WebAssembly target; configure it through emcmake.")
endif()

option(RDP_THREADS "Emit a cross-origin-isolated pthread module" OFF)

add_executable(
  rdp-core
  src/alignment.cpp
  src/bootscan.cpp
  src/burt_confidence.cpp
  src/geneconv.cpp
  src/maxchi.cpp
  src/phylpro.cpp
  src/phylogeny.cpp
  src/rdp_api.cpp
  src/rdp_method.cpp
  src/siscan.cpp
  src/threeseq.cpp
)

target_include_directories(rdp-core PRIVATE include src)
target_compile_features(rdp-core PRIVATE cxx_std_20)
# Modern evergreen browsers all expose WebAssembly SIMD. LTO lets Clang see
# through the fused triplet-method translation units, while SIMD accelerates
# the integer-heavy alignment/profile loops without enabling unsafe floating-
# point reassociation.
target_compile_options(
  rdp-core
  PRIVATE
    -O3
    -flto
    -msimd128
    -Wall
    -Wextra
    -Wpedantic
)

set(RDP_EXPORTED_FUNCTIONS
  "['_malloc','_free','_rdp_create','_rdp_destroy','_rdp_version','_rdp_set_worker_threads','_rdp_load_alignment','_rdp_get_summary_json','_rdp_scan_begin','_rdp_scan_batch','_rdp_reconcile','_rdp_cancel','_rdp_get_progress_json','_rdp_get_results_json','_rdp_get_signal_plot_json','_rdp_get_event_alignment_json','_rdp_get_event_trees_json','_rdp_get_event_phylpro_json','_rdp_set_review_state','_rdp_set_event_review_state','_rdp_update_event','_rdp_update_event_group','_rdp_reconcile_after','_rdp_restore_alignment_begin','_rdp_restore_alignment_record','_rdp_restore_alignment_finish','_rdp_restore_scan_begin','_rdp_restore_signal','_rdp_restore_maxchi_discovery','_rdp_restore_chimaera_discovery','_rdp_restore_geneconv_discovery','_rdp_restore_threeseq_discovery','_rdp_restore_bootscan_discovery','_rdp_restore_siscan_discovery','_rdp_restore_scan_finish','_rdp_restore_event_state','_rdp_restore_reconciliation_required_after','_rdp_export_csv','_rdp_export_enabled_sequences_fasta','_rdp_export_masked_or_disabled_sequences_fasta','_rdp_export_recombinant_sequences_removed_fasta','_rdp_export_recombinant_columns_removed_fasta','_rdp_export_recombination_free_fasta','_rdp_export_fragmented_fasta','_rdp_export_project_json','_rdp_get_error']"
)

target_link_options(
  rdp-core
  PRIVATE
    -O3
    -flto
    -msimd128
    --no-entry
    -sMODULARIZE=1
    -sEXPORT_ES6=1
    -sEXPORT_NAME=createRdpModule
    -sENVIRONMENT=worker
    -sALLOW_MEMORY_GROWTH=1
    -sFILESYSTEM=0
    -sASSERTIONS=0
    -sSTACK_SIZE=1048576
    "SHELL:-sEXPORTED_FUNCTIONS=${RDP_EXPORTED_FUNCTIONS}"
    "SHELL:-sEXPORTED_RUNTIME_METHODS=['UTF8ToString','HEAPU8']"
)

if(RDP_THREADS)
  target_compile_definitions(rdp-core PRIVATE RDP_ENABLE_THREADS=1)
  target_compile_options(rdp-core PRIVATE -pthread)
  target_link_options(
    rdp-core
    PRIVATE
      -pthread
      -sPTHREAD_POOL_SIZE=5
      -sPTHREAD_POOL_SIZE_STRICT=2
  )
  set(RDP_OUTPUT_NAME "rdp-core-threads")
else()
  set(RDP_OUTPUT_NAME "rdp-core")
endif()

set_target_properties(
  rdp-core
  PROPERTIES
    OUTPUT_NAME "${RDP_OUTPUT_NAME}"
    SUFFIX ".mjs"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../public/wasm"
)
