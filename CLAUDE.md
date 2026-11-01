# Project Instructions

## APS IT infrastructure

- `/gdata` (the beamline data, archival and backed up by APS) is not mounted on
  Kouga. It is mounted on Adamite, 10.54.116.17 (four A100 GPUs, 40 GB each),
  and on Amber, 10.54.116.10 (one A100). Run anything that reads `/gdata` on
  Adamite: `ssh adamite '<command>'` from Kouga. On 1 Nov 2026 the user moved
  all jobs off Amber.
- Amber has no access to the external network. For anything that needs it, such
  as GitHub (`git push`, `gh`), use Kouga, 164.54.116.55.
- `/home/8-id-i` (this repository) and `/home/beams/8IDIUSER` (the conda
  environment `env-2021-LLDP`) are mounted on all three machines, so the same
  tree and the same Python are used on each.

## Data sources

- Every input of the analysis comes from `/gdata/s8id-dmdtn`, and every new
  correlation is written there. Nothing is read from `/home/8-id-i`
  outside this repository: not the experiment directories, not earlier
  reanalyses, and not `Archive/`. If something exists only there, ask for it
  to be copied to `/gdata/s8id-dmdtn` first.
- The figure scripts never read `/gdata`. They read the small committed files
  that the reduction scripts (`reduce_gelation.py`, `average_ranges.py`,
  `Read_12ID_SAWAXS.py`, `reduce_flux_control.py`) write when run on Adamite.
