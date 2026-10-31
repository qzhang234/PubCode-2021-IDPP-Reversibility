# Project Instructions

## APS IT infrastructure

- `/gdata` (the beamline data, archival and backed up by APS) is mounted only on
  Amber, 10.54.116.10. Run anything that reads `/gdata` there:
  `ssh 10.54.116.10 '<command>'` from Kouga.
- Amber has no access to the external network. For anything that needs it, such
  as GitHub (`git push`, `gh`), use Kouga, 164.54.116.55.
- `/home/8-id-i` (this repository) and `/home/beams/8IDIUSER` (the conda
  environment `env-2021-LLDP`) are mounted on both machines, so the same tree
  and the same Python are used on either.
