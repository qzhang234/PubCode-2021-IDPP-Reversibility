#!/bin/bash
# Correlate the 2021-2 acquisitions the gelation figures need and the October
# 2026 run (SAXPCS_8id/run_babnigg202107_2_correlations.sh) did not cover:
#
#   D029    the water blank, all 100,000 frames      -> reprocess_results
#   part1   the (VPAVG)10 heating ramp B039, frames 0-49,999   -> reprocess_results_part1
#   part2   the same acquisitions, frames 50,000-99,999        -> reprocess_results_part2
#
# Everything is on /gdata: the raw frames and the q-map are read from the 2021-2
# NeXus tree, and the results are written back into it, the blank beside the
# October run's results and the two halves in directories of their own.  Run
# logs, file lists and chunks go to correlation_run_20261101/.  So it runs on a
# GPU host that mounts /gdata: Adamite (10.54.116.17, four A100s); see CLAUDE.md.
#
# The settings are those of the October run -- the same q-map with the two
# pixels cleared, -t Multitau, --crop-ratio-threshold 0.5, -w -- so these
# results and the ones in reprocess_results differ only in the frames
# correlated.  The frame ranges are those of the August 2025 half-exposure runs
# (-b 0 -e 50000 and -b 50000 -e 100000; end exclusive).
#
#   bash run_gelation_correlations.sh [D029|part1|part2|all]   (default all)
#   MAXMEM=16 setsid nohup bash run_gelation_correlations.sh all &
#
# On Amber's single A100, Oct-Nov 2026: 17 min for the blank, 2 h for each half
# (GPU-bound).  On Adamite's four, 1 Nov 2026: 70 acquisitions a minute, 40 min
# for a half.
#
# A rerun skips acquisitions that already have a result.  To stop a detached
# run, kill its process group (kill -- -PGID): killing the script alone leaves
# the xargs workers running.

set -u

NEXUS=/gdata/s8id-dmdtn/2021-2/babnigg202107_2_nexus
QMAP=$NEXUS/babnigg202107_2_nexus_Sq270_Dq27_Sphi8_Dphi1_q0p0032_Lin_2_BadPixRm.hdf
RUNDIR=$NEXUS/correlation_run_20261101
NWORKERS=${NWORKERS:-4}
CHUNK=${CHUNK:-10}
# -i -2 schedules one GPU per process, so on Adamite each of the four workers
# has an A100 to itself.  On a single-GPU host such as Amber the workers queue
# for it instead; GPU=0 has them share it, each taking at most MAXMEM GB: at the
# correlator's default of 36 GB, four processes run out of Amber's 80 GB.
# MAXMEM sets how the frames are batched, and so the summation order: one
# acquisition of the blank correlated at 16 GB and at 36 GB agrees in I(Q) to
# 1e-12 and in g2 to within 1e-3, the largest differences in the sparsest bins.
GPU=${GPU:--2}
MAXMEM=${MAXMEM:-36}
export PATH=/home/beams/8IDIUSER/bin:$PATH

run_job() {          # name, run directory, output directory, frame arguments
    local name=$1 run=$2 out=$3 frames=$4
    local work=$RUNDIR/$name
    mkdir -p "$out" "$work"
    local list=$work/filelist.txt
    find "$NEXUS/$run" -mindepth 2 -maxdepth 2 -name '*.bin' | sort > "$list"
    # skip acquisitions already done, so a rerun picks up where it stopped
    local todo=$work/todo.txt
    : > "$todo"
    while read -r b; do
        [ -e "$out/$(basename "$b" .bin)_results.hdf" ] || echo "$b" >> "$todo"
    done < "$list"
    echo "$(date '+%F %T') $name: $(wc -l < "$list") acquisitions, $(wc -l < "$todo") to do -> $out"
    rm -rf "$work/chunks"; mkdir -p "$work/chunks"
    split -l "$CHUNK" -d -a 4 "$todo" "$work/chunks/c_"
    export QMAP out work frames GPU MAXMEM
    ls "$work"/chunks/c_* 2>/dev/null | xargs -r -P "$NWORKERS" -I{} bash -c \
        'boost_corr_bin -r $(tr "\n" " " < "$1") -q "$QMAP" -o "$out" $frames \
             -i "$GPU" --max-memory "$MAXMEM" -t Multitau --crop-ratio-threshold 0.5 -w \
             >> "$work/boost.log" 2>&1' _ {}
    local n=0
    while read -r b; do
        [ -e "$out/$(basename "$b" .bin)_results.hdf" ] && n=$((n + 1))
    done < "$list"
    echo "$(date '+%F %T') $name: $n of $(wc -l < "$list") results"
}

D029_RUN=D029_Blank_Water_250C10p_att00_Rq0
B039_RUN=B039_VPAVG_H02_500C10p_att00_Rq0
case "${1:-all}" in
    D029)  run_job D029  $D029_RUN "$NEXUS/reprocess_results" "-b 0 -e -1" ;;
    part1) run_job part1 $B039_RUN "$NEXUS/reprocess_results_part1" "-b 0 -e 50000" ;;
    part2) run_job part2 $B039_RUN "$NEXUS/reprocess_results_part2" "-b 50000 -e 100000" ;;
    all)   run_job D029  $D029_RUN "$NEXUS/reprocess_results" "-b 0 -e -1"
           run_job part1 $B039_RUN "$NEXUS/reprocess_results_part1" "-b 0 -e 50000"
           run_job part2 $B039_RUN "$NEXUS/reprocess_results_part2" "-b 50000 -e 100000" ;;
    *)     echo "usage: $0 [D029|part1|part2|all]"; exit 1 ;;
esac
