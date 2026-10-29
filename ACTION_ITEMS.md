# Paper B — what to measure next, and why

Written 29 Oct 2026, after the split. This replaces MONDAY_TODO.md, whose
correlation and reduction tasks are all done.

Paper B is **"Thermally Reversible Beta-Sheet-like Ordering Controls Gelation
in Elastin-like Polypeptides"**: three figures, 13 pages, builds clean from a
clone. Its claim is that ramp hysteresis is universal and kinetic while
reversibility is sequence-encoded.

As it stands it is a **Biomacromolecules paper** — a sequence-property study of
a recombinant polypeptide, where the SAXS/WAXS plus XPCS combination is
unusual and would be received well. *Macromolecules* and *Soft Matter* are
equivalent fallbacks. It is **not** a Nano Letters or PNAS paper yet, because
the central result is a structure-property correlation on three samples whose
strongest piece of evidence is a photograph of a cloudy cell.

---

## 1. The experiment: paired solutions at 8-ID, dual detector, one shift

Settled 29 Oct. **Everything structural and dynamical in one beamtime**, using
the post-Upgrade dual-detector setup: Rigaku 3M at 22 m for SAXPCS and Lambda
2M at 2.2 m for WAXPCS. The Lambda is on a diffractometer arm, so it reaches
the beta-sheet peaks at 1.38 and 0.72 A^-1 despite the long flight path --
checked, because a centred detector at 2.2 m would fall short by a factor of
five and the whole structural half depends on it.

That is a larger Q range and better Q resolution than the 12-ID-B pair the
existing Figure 1 came from, and it removes the need for a separate 12-ID trip.

### Per chain length, two solutions from one stock

| | before the beamtime | at the beamtime |
|---|---|---|
| **A** | photographed, heated until it gels, photographed, held ~24 h at 4 C, photographed | measured cold, then heated |
| **B** | photographed, never heated, held alongside A | measured cold, then heated |

Then at the beamline, both measured by simultaneous SAXPCS and WAXPCS while
cold, and both taken up through gelation with both detectors running.
Staggering the heating across the QNW zones fits all six in well under a day.

### What each comparison answers

- **A against B while cold** -- did the assembly come apart? SAXS answers it at
  the mesoscale, WAXS answers it for the beta-sheet contacts. Measuring the
  pair side by side in one session removes the beamtime-to-beamtime
  normalisation the present Figure 1 depends on.
- **A against B on heating** -- does a recycled sample re-gel into the same
  thing a virgin one does? Matching cold states shows the assembly dissolved;
  matching gelation runs show it dissolved into the same starting material.
  This is a stronger reversibility statement than the paper currently makes,
  and the design gives it away for free.
- **Across chain lengths** -- whether the beta-sheet peaks track the recovery.
  That is the one that turns the paper's correlation into a mechanism.

### Hold time: one point, deliberately

Twenty-four hours at 4 C, not a staggered series. A staggered series was
considered and dropped: the question here is binary, the existing photograph
already covers 48 h, and samples that have not recovered within a few hours
empirically stay that way. 

**Write the claim to match**: what this supports is *does not recover within
24 h at 4 C*, not *irreversible*. No finite hold establishes the stronger word.

Optional, one extra cell per construct, from the argument that untangling is
thermally activated and deep cooling starves it: hold a third solution at
10-15 C instead of 4 C. If (VPAVG)10 has still not recovered with that much
more thermal energy available, the claim is much harder to argue with.

### Photographs have to be quantitative-capable

Fixed lighting, fixed camera position, grey card in every frame. With that,
transmittance comes out of pixel intensity afterwards and the photographs
become data. Without it they stay illustration, which is exactly the position
the current (VPAVG)10 evidence is in and the thing this is meant to fix.

---

## 2. The same shift should serve Paper A

Paper A's second stated limitation is that its high-Q ordering and low-Q
dynamics come from separate experiments, so no time-resolved causal relation
between them is established. The dual-detector setup removes that outright: one
isothermal hold at 28.50 C with both detectors running gives structure and
dynamics on the same loading, at the same time.

That is the single experiment that would lift Paper A toward a
higher-impact venue, and it is one extra hold on a shift that is happening
anyway. **Build the beamtime request around both papers, not one** -- it is a
far stronger case than either alone.

---

## 3. Raman / FTIR: deferred, and the reason is recorded

An independent spectroscopic confirmation of the beta-sheet assignment, which
currently rests on two peak positions, would close the obvious referee
question.

**Deferred.** Gyorgy's objections are practical rather than scientific: it
needs a lot of sample, and there is no Raman or FTIR instrument to hand. Both
are surmountable if a referee makes it the one thing between the manuscript
and acceptance, but neither justifies the effort before that happens.

If it does come up, the technical note worth remembering: the usual difficulty
with FTIR here is the water bend at ~1640 cm^-1 sitting on top of amide I,
normally dodged with D2O, which shifts the LCST and so changes the thing being
measured. **Raman does not have that problem** -- water is a weak Raman
scatterer and amide I is clean in H2O -- so Raman is the one to ask for, on a
borrowed instrument, rather than FTIR.

---

## 4. More chain lengths: deferred

Twelve, fifteen or twenty-five repeats would show whether reversibility
switches on sharply at some length or crosses over gradually, and a sharp
threshold would be a much stronger result.

**Deferred**: it needs new constructs from Gyorgy. If it ever becomes
possible, the paired-capillary design of item 1 extends to a new construct at
the cost of two more capillaries.

A substitute that needs no cloning: vary **concentration or salt** on the three
constructs already in hand. Both move the LCST and the assembly, and a chain
length x concentration grid is a stronger sequence-property statement than
three points on one axis.

---

## 5. Manuscript tasks, independent of any new measurement

- [ ] **The (VPAVG)10 irreversibility photograph.** Confirmed sound, and the
      text already states what it shows, but the image is in the AdvMat SI and
      no repository here. Drop it in as
      `manuscript/figures/FigureS1_Irreversibility.png` and uncomment the
      float that is already in `main.tex`.
- [ ] **`si.tex` does not exist.** When it does, re-enable
      `\externaldocument[SI-]{build/si}` at the top of `main.tex` and restore
      the `\ref{SI-fig:...}` pointers currently flattened to "SI Appendix".
- [ ] **Materials and Methods** is a placeholder. Assemble from tag
      `pre-split` in `PubCode-2021-IDPP`: the 12-ID-B paragraph, the 8-ID-I
      cell and beam description, and the robotic cell assembly.
- [ ] **The automation material** (Video S2 and the robotic-assembly panel)
      was requested for this build and is not in. Like the photograph, it
      exists only in the AdvMat SI.
- [ ] `Plot_12ID.py` still writes `Figure2_SAXS_WAXS.pdf` and the Makefile
      copies it to `Figure1_...`; `thermal_cycle.py` likewise writes
      `FigureS4_...` and it is moved to `Figure2_...`. Rename the outputs in
      the scripts and drop the copies.

---

## If the experiment works

| outcome | venue |
|---|---|
| Recovery resolved per chain length, and the beta-sheet peaks track it | **PNAS or Nano Letters** is arguable: it is a mechanism, not a correlation |
| Recovery resolved, beta-sheet peaks do not track it | **ACS Nano / Nature Communications**: still a real result, and knowing the two are not linked is worth publishing |
| Cold-state comparison only, no re-gelation data | **Biomacromolecules**, as now, but with the photograph replaced by data |

The risk worth naming: if the beta-sheet peaks do **not** track recovery, both
papers lose the framing they currently lean on. That is a reason to do the
measurement, not to avoid it, but better to know before the framing is built
into a submitted abstract.

## What is no longer in the plan

- **A dedicated turbidimeter.** Superseded: SAXS answers the same question at
  the right length scales and WAXS sees contacts turbidity is blind to. The
  photographs, taken with the discipline above, cover the optical record.
- **A separate 12-ID-B trip.** The dual-detector setup at 8-ID covers the same
  Q range with better resolution, in the same shift as the dynamics.
- **A staggered recovery series.** One 24 h hold, for the reasons in item 1.
