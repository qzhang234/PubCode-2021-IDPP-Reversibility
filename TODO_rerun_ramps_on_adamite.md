# TODO: re-reduce the chain-length ramps on Adamite with the current qmap

Raised 28 Oct 2026. Target: **Monday 2 Nov 2026**. Needs `/gdata` and a GPU
host, so it cannot run from kouga. Use **adamite** or **amazonite**.

## Why

The (VPAVG)10 / 20 / 30 ramp figures in this paper are built from the
*existing* correlation results, produced in 2025 with the qmap of that time:

| sample | runs | beamtime | results read from |
|---|---|---|---|
| (VPAVG)10 | `B039_VPAVG_H02_500C10p`, `B040_..._200C10p` | 2021-2 | `/home/8-id-i/2024-3/2024_1228_qz_llps_analysis/2025_08_analysis/2021-2/cluster_results_all/` |
| (VPAVG)20 | `T214_H04_2nd_500C10p`, `T215_..._200C10p` | 2021-2 | same |
| (VPAVG)30 | `H041_VPAVG_H05_2_500C10p`, `H042_..._200C10p` | 2021-2 | same |

Those are NOT the qmap the 2021-3 isothermal series was reprocessed against on
26 Oct 2026 (`..._Lin_CrossRm2.hdf`, which masks the two stuck pixels
`(30, 340)` and `(63, 166)`). The second pixel was stuck at ~1.5 photons/frame
from some time between 10 and 17 Oct 2021, so whether it affects the 2021-2
ramps at all has not been checked.

## Scoped on adamite, 28 Oct 2026 03:40-04:20

Measured, not estimated:

| | |
|---|---|
| ramp acquisitions to reprocess | **13,241** (B039 2781, B040 2780, T214 2390, T215 2390, H041 1450, H042 1450) |
| adamite throughput, 26 Oct benchmark | 0.70 files/s on 4x A100, NWORKERS=4 |
| correlation time alone | **~5.2 h** |
| NeXus coverage | **complete**: 24,217 metadata files against 24,231 archive acquisitions. No conversion needed. |

So Monday is correlations only, ~5.2 h. It was not started in the 28 Oct
window because that window was 5.3 h and ended when people came back.

### Already done, 28 Oct
`restore_dx_metadata_2021-2.py --apply` repaired the whole tree: 193,736
comparisons over 24,217 NeXus files, 48,230 + 143,729 fields written across two
passes, and a third dry run reports everything already correct. Carried from
the Data Exchange archive: the incident ion chamber (`source_end/I0Monitor`),
the ring current, the sample temperature and set point from the ACTIVE QNW
zone, the zone index itself, and the real acquisition timestamps, which had
read 2025-03-25. Manifest: `analysis/SAXPCS_8id/dx_restore_manifest_2021-2.csv`.

The active zone is not always zone 1: (VPAVG)10 and (VPAVG)20 ramp in zone 1,
(VPAVG)30 in zone 3, with the idle zones parked at 6 C. The ramp ceilings that
come out -- 45.2, 41.4 and 33.1 C -- match SI Figures S2, S3 and S4.

### Only the downstream chamber is unrecoverable
`source_end/TransmissionMonitor` is absent from every 2021-2 acquisition, so
the sample transmission cannot be measured. `abs_xsec.abs_xsec_coef_2021_2()`
supplies it from 2022-1 instead; see the comment at `TRANSMISSION_2021_2`.

## OUTCOME, run of 28 Oct 2026

Started 04:28, stopped by its own DEADLINE at 11:01:40 as instructed. 392 min,
**12,539 of 13,238 done, 699 missing**. No process survived the cutoff.

| run | sample | done | missing | missing indices |
|---|---|---|---|---|
| H041 | (VPAVG)30 up | 1425 / 1449 | 24 | 373-396, contiguous |
| H042 | (VPAVG)30 down | 1450 / 1450 | 0 | - |
| T214 | (VPAVG)20 up | 2305 / 2389 | 84 | 246-989, scattered |
| T215 | (VPAVG)20 down | 2390 / 2390 | 0 | - |
| B039 | (VPAVG)10 up | 2720 / 2780 | 60 | 687-1790, scattered |
| B040 | (VPAVG)10 down | 2249 / 2780 | 531 | 2250-2780, contiguous |

Two different causes, and only one of them is the deadline:

- **531 are the clean tail** of B040, the last run in the queue. Expected.
- **168 are failed chunks**, scattered through H041, T214 and B039. The boost
  log carries **26 MemoryError**. Four concurrent workers exhausted memory
  intermittently. This is the one to fix before Monday: drop to NWORKERS=3, or
  keep 4 and accept a second pass.

The 699 are listed acquisition by acquisition in
`analysis/SAXPCS_8id/missing_2021-2_20261028.txt`. Re-running the script picks
them up: boost_corr gets the whole list but `-w` overwrites, so the finished
12,539 are simply redone, or feed it the missing file directly to save the time.

## What to do

1. Run `analysis/SAXPCS_8id/sweep_stuck_pixels_2021-3.py` (copy it over from
   PubCode-2021-IDPP) against the 2021-2 and 2022-1 prefixes, to find out
   whether either beamtime carries a stuck pixel of its own.
2. Build the matching qmap if it does, the way `make_qmap_babnigg202110.py
   --crossrm2` does: copy the qmap in use and clear the offending pixels, so
   only `mask`, `static_roi_map` and `dynamic_roi_map` change.
3. Re-run the correlations on the ramp acquisitions.
4. Re-run `analysis/ramps/reduce_ramps.py` and `make figures`.
5. Diff the stdout against the baseline saved alongside this file, and update
   any number in the manuscript that moves.

## What is safe meanwhile

The five fitted Q bins of the isothermal work sit far from where either stuck
pixel lands (Q = 0.0214 A^-1 for the 2021-3 one). Nothing in THIS paper is
fitted at that Q either: the ramps are read over 0.002-0.03 A^-1 as whole
profiles and as g2 at a single bin. The figures are expected to be unchanged
to the precision they are drawn at, but that has not been demonstrated, which
is the point of the rerun.
