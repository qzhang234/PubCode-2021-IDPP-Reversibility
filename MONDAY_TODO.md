# Monday 2 Nov 2026 — where the split left things

Two papers now, split on 28 Oct 2026 from one combined Letter. The combined
version is recoverable from tag `pre-split` in `PubCode-2021-IDPP`.

| | repository | branch | state |
|---|---|---|---|
| **A. Arrest / T_s** | `PubCode-2021-IDPP` | `pnas-split` | builds clean, PNAS-ready, needs a reviewer list |
| **B. Reversibility** | `PubCode-2021-IDPP-Reversibility` | `main` | skeleton builds, needs three figures and three paragraphs |

---

## 1. Finish the 2021-2 correlations — 699 acquisitions

The 28 Oct run stopped at its own 11:00 deadline having done **12,539 of
13,238**. The 699 that did not run are listed one per line in

    analysis/SAXPCS_8id/missing_2021-2_20261028.txt

**Two causes, and only one is the deadline.**

| run | sample | done | missing | indices | cause |
|---|---|---|---|---|---|
| H041 | (VPAVG)30 up | 1425 / 1449 | 24 | 373-396 contiguous | failed chunk |
| H042 | (VPAVG)30 down | 1450 / 1450 | 0 | | |
| T214 | (VPAVG)20 up | 2305 / 2389 | 84 | 246-989 scattered | failed chunks |
| T215 | (VPAVG)20 down | 2390 / 2390 | 0 | | |
| B039 | (VPAVG)10 up | 2720 / 2780 | 60 | 687-1790 scattered | failed chunks |
| B040 | (VPAVG)10 down | 2249 / 2780 | 531 | 2250-2780 contiguous | deadline |

531 are the clean tail of the last run queued. The other **168 are failed
chunks**, and the boost log carries **26 `MemoryError`**: four concurrent
workers intermittently exhausted memory.

**Already scheduled.** `at` job 6 on adamite fires at **01:00 Thu 29 Oct**,
running `~/llps_2021-2_restore/at_job_pass2.sh` with `NWORKERS=3` and a 06:00
safety deadline. `filelist.txt` in the run directory has been pointed at the
699; the full 13,238 is preserved as `filelist_full_20261028.txt`. Expected
~30 min, so check for the completion email, then:

```bash
ssh adamite 'atq'                     # should be empty by morning
ssh adamite 'ls /gdata/s8id-dmdtn/2021-2/babnigg202107_2_nexus/reprocess_results/*_results.hdf | wc -l'
# expect 13238; if short, a fresh missing_adamite.txt lists what is left
```

If it has to be re-run by hand:

```bash
ssh adamite
cd ~/llps_2021-2_restore
NWORKERS=3 DEADLINE=$(date -d '2026-11-02 18:00' +%s) ./run_babnigg202107_2_correlations.sh
```

Note the trap: a plain re-run with the ORIGINAL filelist redoes all 13,238
(~6.9 h), because `-w` overwrites. Only the swapped filelist gives the 30-minute
version.

`NWORKERS=3` is the fix for the MemoryErrors; at 4 it loses about 1.3% of
chunks. The script rebuilds its own file list and `-w` overwrites, so a plain
re-run is safe and simply redoes the 12,539 as well. To do only the remainder,
point `LIST` at the missing file instead — it is already in the right format
(one absolute `.bin` path per line).

It emails qzhang234@anl.gov on any exit, including a crash or a deadline cut.
**Note: the 28 Oct completion email was sent but delivery was never confirmed.
If nothing arrived, check `/usr/sbin/sendmail` on adamite before relying on it.**

---

## 2. Then reduce and plot the ramps — nothing exists yet

`analysis/ramps/` is an empty directory. Two scripts to write, modelled on
`SAXPCS_8id/temperature_series_2021-3.py`:

- **`reduce_ramps.py`** — read `reprocess_results/{H041,H042,T214,T215,B039,B040}_*`,
  group consecutive acquisitions, apply the two outlier cuts of
  `average_ranges.py`, write averages into `analysis/ramps/data/`.
- **`plot_ramps.py`** — Figure 2 of paper B: the three chain lengths, up and
  down ramps, I(Q) and g2, in the format of the old `H02.pdf` / `H04.pdf` /
  `H06.pdf` in `Archive/Old/*_SI/`.

**Everything they need is already settled:**

- The q and tau axes are in the result files: `/xpcs/multitau/delay_list` (61,)
  and `/xpcs/qmap/dynamic_v_list_dim0` (27,). The old pickles lacked them,
  which is what blocked this before.
- Per-acquisition temperature is at `/entry/sample/temperature`, repaired on
  28 Oct from the Data Exchange archive.
- **Use `/entry/sample/qnw_zone`, do not assume zone 1.** (VPAVG)10 and
  (VPAVG)20 ramp in zone 1, (VPAVG)30 in zone 3, idle zones parked at 6 C.
  Assuming zone 1 puts a flat 6 C on the (VPAVG)30 temperature axis.
- Absolute scale: call `abs_xsec.abs_xsec_coef_2021_2()`, **not**
  `abs_xsec_coef()`. 2021-2 never recorded the downstream ion chamber, so the
  transmission is taken from 2022-1 (`TRANSMISSION_2021_2 = 0.369`). Read the
  comment there before quoting any absolute intensity: it is good to about 6%,
  it is a common scale factor, and 2021-2 intensities should be quoted to one
  significant figure.

### Sample mapping, confirmed three ways

| | runs | ramp ceiling | old SI figure |
|---|---|---|---|
| **(VPAVG)10** | B039 up, B040 down | 45.2 C | S2 |
| **(VPAVG)20** | T214 up, T215 down | 41.4 C | S3 |
| **(VPAVG)30** | H041 up, H042 down | 33.1 C | S4 |

Confirmed by the ramp ceilings, by the run names, and by the cell photographed
in old SI Figure S5 being labelled "H 02" and captioned (VPAVG)10.

---

## 3. Paper B — what is written and what is not

`manuscript/main.tex` builds (9 pages, `make papers`, `make check` clean).

**Already in it**, carried over from the combined Letter:

- The seven-cycle repeatability experiment and **Figure 2**, added 28 Oct.
  `thermal_cycle.py` reproduces it from the 2022-1 `data/`, now owned by this
  repository.

- The amyloid introduction, both paragraphs. This is its proper home: in paper
  A it set up a structure-dynamics link the text then disclaimed four times;
  here it sets up the actual claim.
- The 12-ID-B results paragraph and the hierarchical-assembly paragraph.
- Figure 1, the 12-ID-B SAXS/WAXS panel with the beta-sheet peak fits.
- A draft abstract.

**Still to write** — marked with a comment block in `main.tex`:

- [ ] the chain-length ramp paragraph and **Figure 2** (needs §1 and §2 above)
- [ ] irreversibility of (VPAVG)10, and what chain length does
- [ ] **the automated cell assembly**: Video S2 and Figure S1b, from the AdvMat
      version. Requested for this build; the source is the AdvMat SI PDF, not
      in any repo here.
- [ ] `si.tex` — does not exist yet. When it does, re-enable
      `\externaldocument[SI-]{build/si}` at the top of `main.tex` and restore
      the `\ref{SI-fig:...}` pointers that are currently flattened to
      "SI Appendix".
- [ ] Materials and Methods — assemble from tag `pre-split` of the other repo.

### Loose ends

- `Plot_12ID.py` still writes `Figure2_SAXS_WAXS.pdf`; the Makefile copies it
  to `Figure1_SAXS_WAXS.pdf`. Rename the output in the script and drop the copy.
- `dx_restore_manifest_2021-2.csv` is 46 MB. Gitignore it and keep a summary.
- `Archive/` is 644 MB of provenance. Decide whether it belongs in the
  published repo or only locally.

---

## 4. Paper A — PNAS, essentially ready

Branch `pnas-split` in `PubCode-2021-IDPP`. `make papers` and `make check` clean.
**3,690 words, 26 main-text references, 3 figures** against PNAS's ~4,000 / 50 /
4 for a six-page Research Report. Initial submission is format-neutral, so no
document-class conversion is needed.

Done: title, 118-word significance statement, classification, keywords,
Introduction/Results/Discussion/Materials and Methods order.

- [ ] **Reviewer lists — only you can write these.** PNAS requires a minimum of
      three each of Editorial Board members, NAS members, and qualified
      reviewers. The old ACS cover letter's four suggestions do not cover it.
- [ ] Rewrite the cover letter for PNAS; mention the companion paper.
- [ ] **Put the PRL 2017 connection in it.** `Zhang2017-ip` (PRL 119, 178006),
      which introduced the exp(dT/T_s) scaling this paper applies, is your own
      first-author paper, and you are second author on `Bahadur2019-ln`. As
      written the cover letter cites both as other people's prior work.
- [ ] SI figure labels still read `fig:S5`..`fig:S11` for figures now numbered
      S2..S8. LaTeX numbers them correctly; only the source keys are stale.
- [x] The absolute-scale calibration figure now lives in **neither**. It was
      the last thing in A reading the 2022-1 data, so it came out of A's SI and
      A's two pointers cite this paper instead. **B must therefore show it**:
      `saxpcs.py` still writes `FigureS5_Calibration.pdf` from `data/`, so add
      it to B's SI when `si.tex` is written.
- [ ] Still open from the audit: `SOURCE_DATE_EPOCH` in the Makefile for
      byte-reproducible figures, and the main text's claim that the 2021-3
      loadings were "held at 6 C for more than 30 min, warmed to 20 C", which
      the SI does not record.

---

## 5. What changed on `/gdata` on 28 Oct

All in `/gdata/s8id-dmdtn/2021-2/babnigg202107_2_nexus`, none of it touching
the archive:

- **24,217 `.bin` symlinks repaired.** Every one pointed at
  `/rsyncd-munged//gdata/...`, which does not exist, so the first correlation
  attempt produced nothing while reporting `rc=0` per chunk. The stray prefix
  was stripped and each target verified to exist first.
- **Metadata restored into 24,217 NeXus files** by
  `analysis/SAXPCS_8id/restore_dx_metadata_2021-2.py --apply`: incident ion
  chamber, ring current, sample temperature and set point from the active QNW
  zone, the zone index, and real acquisition timestamps, which had read
  2025-03-25. Idempotent; a re-run reports everything already correct.
- **A new qmap**, `..._q0p0032_Lin_2_BadPixRm.hdf`, a copy of `_Lin_2` with
  pixels (30, 340) and (63, 166) cleared in `mask`, `static_roi_map` and
  `dynamic_roi_map`. Valid fraction unchanged at 0.8601. Verified separately:
  Q < 0.003 is fully masked, the cross-shaped parasitic mask and beamstop disc
  are present and centred on the direct beam, and **2021-2 has no hot pixels**
  — zero unmasked pixels above 20x their q-ring median over 8 acquisitions.
