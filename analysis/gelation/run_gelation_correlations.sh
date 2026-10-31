#!/bin/bash
# Correlate the 2021-2 acquisitions the gelation figures need and the October
# 2026 run (SAXPCS_8id/run_babnigg202107_2_correlations.sh) did not cover:
#
#   D029    the water blank, all 100,000 frames            -> 2021-2_D029
#   part1   the (VPAVG)10 heating ramp B039, frames 0-49,999   -> 2021-2_B039_part1
#   part2   the same acquisitions, frames 50,000-99,999        -> 2021-2_B039_part2
#
# BEAMLINE-ONLY: reads the raw frames on /gdata, so it runs on Amber
# (10.54.116.10), which has the A100 and /gdata; see CLAUDE.md.  The settings
# are those of the October run -- the same q-map with the two pixels cleared,
# -t Multitau, --crop-ratio-threshold 0.5, -w -- so these results and the ones
# in reprocess_results differ only in the frames correlated.  The frame ranges
# are those of the August 2025 half-exposure runs (-b 0 -e 50000 and
# -b 50000 -e 100000; end exclusive).
#
# Results are written beside the August 2025 analysis, not into /gdata.
#
#   bash run_gelation_correlations.sh [D029|part1|part2|all]   (default all)
#   GPU=0 NWORKERS=4 MAXMEM=16 setsid nohup bash run_gelation_correlations.sh all &
#
# On Amber, 31 Oct 2026: 15 min for the blank, 2 h for each half (GPU-bound).
#
# A rerun skips acquisitions that already have a result.  To stop a detached
# run, kill its process group (kill -- -PGID): killing the script alone leaves
# the xargs workers running.

set -u

NEXUS=/gdata/s8id-dmdtn/2021-2/babnigg202107_2_nexus
QMAP=$NEXUS/babnigg202107_2_nexus_Sq270_Dq27_Sphi8_Dphi1_q0p0032_Lin_2_BadPixRm.hdf
OUT=/home/8-id-i/2024-3/2024_1228_qz_llps_analysis/2026_10_gdata_reanalysis
NWORKERS=${NWORKERS:-4}
CHUNK=${CHUNK:-10}
# -i -2 schedules one GPU per process.  Amber has a single A100, so the workers
# then queue for it (about 13 acquisitions a minute for the blank); GPU=0 has
# them share it instead.  Each process may then take at most MAXMEM GB of it:
# at the correlator's default of 36 GB, four processes run out of the 80 GB.
GPU=${GPU:--2}
MAXMEM=${MAXMEM:-36}
export PATH=/home/beams/8IDIUSER/bin:$PATH

run_job() {          # name, run directory, frame arguments
    local name=$1 run=$2 frames=$3
    local out=$OUT/2021-2_$name
    mkdir -p "$out"
    local list=$out/filelist.txt
    find "$NEXUS/$run" -mindepth 2 -maxdepth 2 -name '*.bin' | sort > "$list"
    # skip acquisitions already done, so a rerun picks up where it stopped
    local todo=$out/todo.txt
    : > "$todo"
    while read -r b; do
        [ -e "$out/$(basename "$b" .bin)_results.hdf" ] || echo "$b" >> "$todo"
    done < "$list"
    echo "$(date '+%F %T') $name: $(wc -l < "$list") acquisitions, $(wc -l < "$todo") to do"
    rm -rf "$out/chunks"; mkdir -p "$out/chunks"
    split -l "$CHUNK" -d -a 4 "$todo" "$out/chunks/c_"
    export QMAP out frames GPU MAXMEM
    ls "$out"/chunks/c_* 2>/dev/null | xargs -r -P "$NWORKERS" -I{} bash -c \
        'boost_corr_bin -r $(tr "\n" " " < "$1") -q "$QMAP" -o "$out" $frames \
             -i "$GPU" --max-memory "$MAXMEM" -t Multitau --crop-ratio-threshold 0.5 -w \
             >> "$out/boost.log" 2>&1' _ {}
    local n=$(ls "$out"/*_results.hdf 2>/dev/null | wc -l)
    echo "$(date '+%F %T') $name: $n of $(wc -l < "$list") results"
}

case "${1:-all}" in
    D029)  run_job D029  D029_Blank_Water_250C10p_att00_Rq0 "-b 0 -e -1" ;;
    part1) run_job B039_part1 B039_VPAVG_H02_500C10p_att00_Rq0 "-b 0 -e 50000" ;;
    part2) run_job B039_part2 B039_VPAVG_H02_500C10p_att00_Rq0 "-b 50000 -e 100000" ;;
    all)   run_job D029  D029_Blank_Water_250C10p_att00_Rq0 "-b 0 -e -1"
           run_job B039_part1 B039_VPAVG_H02_500C10p_att00_Rq0 "-b 0 -e 50000"
           run_job B039_part2 B039_VPAVG_H02_500C10p_att00_Rq0 "-b 50000 -e 100000" ;;
    *)     echo "usage: $0 [D029|part1|part2|all]"; exit 1 ;;
esac
