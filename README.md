# Thermally Reversible Beta-Sheet-like Ordering and Gelation in Elastin-like Polypeptides

The structural half of what was one Letter until 28 Oct 2026. The dynamical
half — isothermal arrest and the T_s scaling — is a separate paper in
[PubCode-2021-IDPP](https://github.com/qzhang234/PubCode-2021-IDPP); the
combined version is recoverable there from the tag `pre-split`.

```
make figures    regenerate every data-derived figure
make papers     compile main.pdf into manuscript/build/
make check      report undefined citations / references / missing figures
```

Every data-derived figure regenerates from a clone of this repository alone:
the averaged 8-ID-I files both figures read are committed, 151 MB of them, so
no beamline account and no `/gdata` are needed. The two scripts that DO read
beamline storage, `SAXPCS_8id/restore_dx_metadata_2021-2.py` and
`ramps/reduce_ramps.py`, produced those committed files and are kept for
provenance; `make figures` never runs them.

**This paper is mid-rebuild.** `MONDAY_TODO.md` is the authoritative list of
what exists, what does not, and what has to happen next. Read it first.

## Layout

| Path | What it is |
|---|---|
| `manuscript/` | `main.tex` (builds, 9 pp), `reference.bib`, vendored `texmf/`. No `si.tex` yet. |
| `analysis/SAXS_12id/` | 12-ID-B reduction and the beta-sheet peak fits -> Figure 1 |
| `analysis/SAXPCS_8id/` | the 2022-1 seven-cycle data and `thermal_cycle.py`; the 2021-2 metadata repair and correlation runner |
| `analysis/ramps/` | the chain-length ramps: `reduce_ramps.py` (beamline-only), `plot_ramps.py` -> Figure 3, and the committed `data/` |
| `analysis/common/` | `acs_style.py`, shared with the other repository |

## Samples

| | runs (2021-2) | ramp ceiling |
|---|---|---|
| (VPAVG)10 | `B039` up, `B040` down | 45.2 C |
| (VPAVG)20 | `T214` up, `T215` down | 41.4 C |
| (VPAVG)30 | `H041` up, `H042` down | 33.1 C |

## One thing to know before using 2021-2 intensities

That beamtime never recorded the downstream ion chamber, so the sample
transmission cannot be measured. `abs_xsec.abs_xsec_coef_2021_2()` supplies it
from 2022-1 instead. The absolute scale is then good to about 6%, as a common
scale factor that moves no peak, exponent or within-profile ratio. Quote
2021-2 intensities to one significant figure.
