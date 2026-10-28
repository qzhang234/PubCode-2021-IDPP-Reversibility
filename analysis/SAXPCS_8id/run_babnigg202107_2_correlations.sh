#!/bin/bash
# Correlate the six chain-length ramp runs of the 2021-2 beamtime.
#
# BEAMLINE-ONLY. Needs /gdata and a GPU host (adamite, amazonite). Written for
# the 28 Oct 2026 shutdown window.
#
# Adapted from run_babnigg202110_correlations.sh. Two things differ.
#
# LAYOUT. The 2021-3 tree is flat: <nexus>/<acq>/<acq>.bin. The 2021-2 tree
# nests each acquisition inside its scan, <nexus>/<scan>/<acq>/<acq>.bin, so
# the file list is built by searching for the .bin rather than by listing one
# level. The flat assumption is what made an earlier pass over this tree see
# 262 of 24,217 files.
#
# ORDER. Acquisitions are listed sample by sample, smallest first, so that a
# deadline cut loses whole samples from the end rather than the top of every
# ramp. The top of a ramp is the assembled state, which is the part worth
# having, and a run that is 90% complete is missing exactly that. Order:
# (VPAVG)30 (2,900), (VPAVG)20 (4,780), (VPAVG)10 (5,561).
#
# QMAP. ..._q0p0032_Lin_2_BadPixRm.hdf: the qmap the existing results were made
# against, with (30, 340) and (63, 166) cleared in mask, static_roi_map and
# dynamic_roi_map. Those two are the pixels masked for 2021-3. Whether either
# is stuck in 2021-2 has NOT been established -- (63, 166) failed between 10
# and 17 Oct 2021, after this July beamtime, so it probably is not -- but
# masking a healthy pixel costs one part in 524,288 and leaves the valid
# fraction at 0.8601, while leaving a stuck one unmasked corrupts its whole q
# bin. The asymmetry decides it. sweep_stuck_pixels_2021-2.py settles the
# question separately.
#
# The other design notes of the 2021-3 runner apply unchanged: -i -2 takes one
# GPU per process so parallelism comes from NWORKERS invocations;
# --crop-ratio-threshold 0.5 keeps cropping off, which otherwise returns
# all-NaN g2; -w overwrites rather than appending _1.

set -u

HOST=$(hostname -s)
NEXUS=/gdata/s8id-dmdtn/2021-2/babnigg202107_2_nexus
QMAP=$NEXUS/babnigg202107_2_nexus_Sq270_Dq27_Sphi8_Dphi1_q0p0032_Lin_2_BadPixRm.hdf
OUT=$NEXUS/reprocess_results
RUNDIR=$NEXUS/correlation_run_20261028
NWORKERS=${NWORKERS:-4}
CHUNK=${CHUNK:-12}
# 2026-10-28 11:00 local, set by the user: everything on adamite stops then.
# No new chunk starts after it; chunks already running finish, which takes
# under a minute. At the measured 0.70 files/s a 06:00 start covers about
# 12,000 of the 13,241, so (VPAVG)30 and (VPAVG)20 should complete and
# (VPAVG)10 should not.
DEADLINE=${DEADLINE:-1793203200}
EMAIL=${EMAIL:-qzhang234@anl.gov}

# Sample by sample, smallest first. H04x = (VPAVG)30, T21x = (VPAVG)20,
# B03x/B04x = (VPAVG)10; the mapping is set out in the manuscript README.
RUNS=(H041_VPAVG_H05_2_500C10p_att00_Rq0
      H042_VPAVG_H05_2_200C10p_att00_Rq0
      T214_H04_2nd_500C10p_att00_Rq0
      T215_H04_2nd_200C10p_att00_Rq0
      B039_VPAVG_H02_500C10p_att00_Rq0
      B040_VPAVG_H02_200C10p_att00_Rq0)

export PATH=/home/beams/8IDIUSER/bin:$PATH
mkdir -p "$RUNDIR" "$OUT"
LOG=$RUNDIR/progress_${HOST}.log
STATE=$RUNDIR/state_${HOST}

log() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }

notify() {
    [ -n "$EMAIL" ] || return 0
    printf 'To: %s\nFrom: 8idiuser@%s.xray.aps.anl.gov\nSubject: %s\n\n%s\n' \
        "$EMAIL" "$HOST" "$1" "$2" | /usr/sbin/sendmail -t -i 2>>"$LOG"
}

FINISHED=0
on_exit() {
    rc=$?
    [ "$FINISHED" = 1 ] && return 0
    notify "babnigg202107_2 correlation run ENDED EARLY on $HOST (rc=$rc)" \
"The run stopped before reporting a normal finish.

host          : $HOST
run directory : $RUNDIR

last log lines
--------------
$(tail -20 "$LOG" 2>/dev/null)"
}
trap on_exit EXIT

log "=== start on $HOST, workers=$NWORKERS chunk=$CHUNK deadline=$(date -d @$DEADLINE '+%F %T')"
[ -f "$QMAP" ] || { log "FATAL: qmap missing: $QMAP"; exit 1; }
log "qmap $(basename "$QMAP")"

nd=$(ls "$OUT"/*_results_1.hdf 2>/dev/null | wc -l)
[ "$nd" -gt 0 ] && { rm -f "$OUT"/*_results_1.hdf; log "removed $nd stale _1 files"; }

# Build the list sample by sample; the .bin sits one level below the scan dir.
LIST=$RUNDIR/filelist.txt
if [ ! -s "$LIST" ]; then
    : > "$LIST"
    for r in "${RUNS[@]}"; do
        n=$(find "$NEXUS/$r" -mindepth 2 -maxdepth 2 -name '*.bin' 2>/dev/null | sort | tee -a "$LIST" | wc -l)
        log "  $r : $n acquisitions"
    done
fi
TOTAL=$(wc -l < "$LIST")
log "acquisitions to process: $TOTAL"

CHUNKDIR=$RUNDIR/chunks_${HOST}
rm -rf "$CHUNKDIR"; mkdir -p "$CHUNKDIR"
split -l "$CHUNK" -d -a 4 "$LIST" "$CHUNKDIR/c_"
log "split into $(ls "$CHUNKDIR" | wc -l) chunks of $CHUNK"

run_chunk() {
    c="$1"
    if [ "$(date +%s)" -ge "$DEADLINE" ]; then return 0; fi
    boost_corr_bin -r $(tr '\n' ' ' < "$c") -q "$QMAP" -o "$OUT" \
        -i -2 -t Multitau --crop-ratio-threshold 0.5 -w \
        >> "$RUNDIR/boost_${HOST}.log" 2>&1
    echo "$(date '+%F %T') $(basename "$c") rc=$?" >> "$STATE"
}
export -f run_chunk
export QMAP OUT RUNDIR HOST STATE DEADLINE PATH

T0=$(date +%s)
# Chunks are consumed in name order, which is list order, which is sample order.
ls "$CHUNKDIR"/c_* | sort | xargs -r -P "$NWORKERS" -I{} bash -c 'run_chunk "$@"' _ {}
T1=$(date +%s)

# What did not run, per sample, so Monday picks up exactly the remainder.
MISSING=$RUNDIR/missing_${HOST}.txt
: > "$MISSING"
while read -r b; do
    a=$(basename "$b" .bin)
    [ -e "$OUT/${a}_results.hdf" ] || echo "$b" >> "$MISSING"
done < "$LIST"
NMISS=$(wc -l < "$MISSING")

PERSAMPLE=""
for r in "${RUNS[@]}"; do
    want=$(grep -c "/$r/" "$LIST" || true)
    miss=$(grep -c "/$r/" "$MISSING" || true)
    PERSAMPLE="$PERSAMPLE
  $(printf '%-38s %5d of %5d done, %5d missing' "$r" "$((want-miss))" "$want" "$miss")"
done

DONE=$((TOTAL-NMISS))
log "=== finished in $(( (T1-T0)/60 )) min: $DONE/$TOTAL, $NMISS missing"
log "$PERSAMPLE"

FINISHED=1
if [ "$NMISS" -eq 0 ]; then STATUS="complete, $DONE/$TOTAL"; else STATUS="$DONE/$TOTAL done, $NMISS left for Monday"; fi
notify "babnigg202107_2 correlation run finished on $HOST: $STATUS" \
"Chain-length ramp correlations, 2021-2.

elapsed    : $(( (T1-T0)/60 )) min
qmap       : $(basename "$QMAP")
results in : $OUT
deadline   : $(date -d @$DEADLINE '+%F %T') (set to 11:00, user instruction)

per sample
----------$PERSAMPLE

What did not run is listed acquisition by acquisition in
  $MISSING
Re-running this script picks up only those: boost_corr is given the whole list
but -w overwrites, so the simplest Monday action is to point NWORKERS at the
free GPUs and run it again with a later DEADLINE.

last log lines
--------------
$(tail -25 "$LOG" 2>/dev/null)"
