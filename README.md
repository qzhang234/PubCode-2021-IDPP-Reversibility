# Thermally Reversible and Microscopically Tunable Gelation in Recombinant Polypeptides

SA-XPCS of the thermally induced gelation of the elastin-like polypeptides
(VPAVG)10, (VPAVG)20 and (VPAVG)30 on heating and cooling: a hierarchical
network forms in all three and its dynamics arrest through a two-step decay,
while whether the gel dissolves again depends on chain length.

Since 29 Oct 2026 the manuscript reproduces the figure set of the Overleaf
version of the paper (the Advanced Materials and ACS submissions, PDFs in the
untracked `Archive/Overleaf/`), with every data-derived figure regenerated from
committed data by the scripts here, which are ports of the original analysis in
`Archive/Old/` and `Archive/New/`. Two figures of the previous version (28 Oct)
were added back to the SI on 30 Oct: the cross-beta-like WAXS peak fits of
12-ID-B, now Figure S5a, and the seven-cycle thermal cycling of 2022-1, now
Figure S6. The rest of that version, chiefly its chain-length ramps in 2 C
windows, is in the git history and builds with `make previous-figures`. The dynamical companion paper is
[PubCode-2021-IDPP](https://github.com/qzhang234/PubCode-2021-IDPP).

```
make figures    regenerate every data-derived figure of main.tex and si.tex
make papers     compile main.pdf and si.pdf into manuscript/build/
make check      report undefined citations / references / missing figures
```

## Reproducibility

Every figure regenerates from a clone of this repository alone. The reduced
data are committed, and none of the figure scripts reads beamline storage:

| Data | Made by | From |
|---|---|---|
| `analysis/gelation/data/` (4 MB) | `gelation/reduce_gelation.py` | August 2025 correlation results on `/home/8-id-i/2024-3/...` (raw frames on `/gdata`) |
| `analysis/SAXPCS_8id/data/` | `SAXPCS_8id/average_ranges.py` | October 2026 reprocessing on `/gdata/s8id-dmdtn/2022-1/` |
| `analysis/SAXS_12id/reduced_data/` | `SAXS_12id/Read_12ID_SAWAXS.py` | 12-ID-B processed files on `/gdata/s8id-dmdtn/2021-1/` |
| `analysis/Rad_Dam_Check/cluster_results/` | copied from `Archive/New/` | 2025 flux test, raw on `/gdata/dm/8ID/8IDI/2025-1/` |

`Archive/` is local provenance and is never read by `make`.

`analysis/gelation/reduce_gelation.py` produced `gelation/data/`. It reads the
per-acquisition correlation results of the August 2025 reanalysis under
`/home/8-id-i/2024-3/2024_1228_qz_llps_analysis/2025_08_analysis/` and the
2021-2 restore manifest (`analysis/SAXPCS_8id/dx_restore_manifest_2021-2.csv`,
46 MB, untracked), so it runs only on an 8-ID machine, and `make` never runs
it. Its section averages reproduce the original pickles in `Archive/Old/`
exactly in g2 and g2_err for all eleven datasets, and in I(Q) everywhere except
the corner Q bins that the original reader corrupted (below).

## Figures

| Figure | Script | Data |
|---|---|---|
| 1 | artwork, `manuscript/figures/Setup.pdf` (from the ACS Overleaf project) | |
| 2 | `gelation/saxs_scaling.py` | E0110 heating, D0138 buffer; 12-ID-B R_g |
| 3 | `gelation/g2_arrest.py` | B039 heating |
| 4 and its table | `gelation/two_step_fit.py` | B039, T214, E0110 in 20 sections |
| S1, S7, S8, Videos S1-S2 | photographs and drawings, extracted from the Overleaf SI | |
| S2-S4 | `gelation/ramp_si.py` | all six ramps, D029 and D0138 |
| S5 | `SAXS_12id/reversibility_12id.py` (fits in `saxs_waxs_fits.py`, `guinier.py`) | 12-ID-B, 2021-1 |
| S6 | `SAXPCS_8id/thermal_cycle.py` | 2022-1 B0075-B0090, buffer D0077/D0080, static reference F0145 |
| S9 | typeset in `si.tex` | |
| S10 | `SAXPCS_8id/ic_calibration.py` | constants in `abs_xsec.py` |
| S11, S12 | `gelation/stability_si.py` | B039 subsets and half-exposure correlations |
| S13 | `Rad_Dam_Check/flux_control.py` | E0171-E0174 (2025) |

## Where the reproduction departs from the Overleaf figures

Each change is explained in the docstring of the script that makes it.

* **Absolute scale.** The original used one fixed pair of coefficients
  (6.93e4, 7.62e4) for both beamtimes, computed for a 1 mm path through a 2 mm
  cell and labelled cm^-1. Here each section is scaled by its own restored
  ion-chamber readings through `abs_xsec.py`, in mm^-1; for 2022-1 that is
  within 5 % of half the old values.
* **Figure 2, T0 and gamma.** The original weighted the power-law fit by the
  chi^2 gradient at the grid minimum and pinned it to the coldest point with an
  error of 1e-6. Fed the original inputs, the port returns the original
  T0 = 29.96 C and gamma = -0.55 exactly, but that weighting gives anything from
  -0.55 to -0.92 for trivially different inputs. Unweighted, every input gives
  T0 = 30.5 C and gamma = -0.80 +/- 0.04, which is what the paper reports. The
  hand-drawn eta = 2.8 guide is now a fit (eta = 2.92 +/- 0.02); the T - T0 axis
  label is corrected to T0 - T; R_g is fitted from the committed 12-ID-B data.
* **Figure 4.** The original weighted all three fits with the (VPAVG)10 error
  bars. With each sample's own, (VPAVG)10 is unchanged, (VPAVG)30 shifts
  slightly and the (VPAVG)20 slow time goes from 84 to 164 +/- 44 s.
* **Figures S2-S4.** Time axes use the recorded start times (5 s per
  acquisition in 2021-2, 6 s in 2022-1; the original assumed 6 s). The rule
  that hides unresolvable g2 curves, which the original applied to S4 only, is
  applied to all three and stated in the caption.
* **Corner Q bins.** The original outlier routine wrote 1.0 into missing or
  non-positive I(Q) values in place; for one buffer acquisition this made the
  background-subtracted I(Q) above 0.030 A^-1 hugely negative. The corrected
  routine in `common/utils.py` is used, and I(Q) panels stop at 0.032 A^-1,
  where the original stopped.
* **Units.** Q is in A^-1 and g2 is labelled with the actual bin, 0.0038 A^-1
  (the original read 0.03 nm^-1).

## Samples

| | heating | cooling | ceiling | beamtime | buffer |
|---|---|---|---|---|---|
| (VPAVG)10 | B039 | B040 | 45.2 C | 2021-2 | D029 (water) |
| (VPAVG)20 | T214 | T215 | 41.4 C | 2021-2 | D029 (water) |
| (VPAVG)30 | E0110 | E0111 | 32.8 C | 2022-1 | D0138 |

2021-2 never recorded the downstream ion chamber, so for (VPAVG)10 and
(VPAVG)20 the transmission is borrowed from 2022-1 (`abs_xsec.TRANSMISSION_2021_2`).
It is a common scale factor on every 2021-2 curve, good to about 6 %, and moves
no peak, exponent or within-profile ratio.
