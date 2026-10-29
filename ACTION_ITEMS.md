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

## 1. The one experiment that matters: paired capillaries at 12-ID-B

**This is the whole plan.** For each of the three chain lengths, two capillaries
from one stock:

| | |
|---|---|
| **Cycled** | taken through the full thermal ramp, then returned to 4 C |
| **Reference** | never heated, held at 4 C throughout |

Measure both at 10 C and 30 C by simultaneous SAXS/WAXS, exactly as the
existing (VPAVG)30 pair in Figure 1 was measured. That generalises Figure 1
from one chain length to three and answers the question the paper actually
asks, at both length scales at once:

- **Mesoscale**, from SAXS: does the low-Q assembly dissolve?
- **Molecular**, from WAXS: do the beta-sheet-like peaks at 1.38 and
  0.72 A^-1 go away?

The second is the one that turns a correlation into a mechanism. The paper
currently reports beta-sheet peaks and reversibility side by side and is
careful not to link them. If the peaks return to the reference profile for
thirty repeats and persist for ten, the junctions ARE the beta-sheet contacts,
and that is a different and much stronger paper.

### Make it a time series, not just an endpoint

One capillary per condition answers "did it recover by day N". Several
capillaries of the same construct, pulled from the fridge at **1 h, 6 h, 24 h,
48 h and 1 week**, give a coarse recovery curve per chain length as well. Same
beamtime, more capillaries, no extra instrument.

That matters because of a trap in the endpoint-only design: if (VPAVG)10 in
fact recovers at two weeks, a 48 h endpoint calls it irreversible and the
claim is wrong. A time series shows the trend; a single endpoint cannot
distinguish "very slow" from "never".

### Cheap insurance against that trap

Leave one capillary of each construct on the bench and photograph it weekly
against a fixed background, for as long as you can be bothered. It costs
nothing, it runs long after the beamtime ends, and it is the only thing that
can distinguish a two-week time constant from a genuinely arrested state.

---

## 2. Dedicated turbidity: probably not needed, with one caveat

The earlier recommendation was a benchtop turbidimeter. With the 12-ID
experiment above, **it is largely redundant**, because SAXS/WAXS answers the
same question better: turbidity integrates scattering from roughly micron-scale
heterogeneity and cannot distinguish a network that dissolved from one that
coarsened, while SAXS measures the mesoscale assembly directly and WAXS sees
the contacts turbidity is blind to.

**Caveat.** Turbidity's one real advantage is unlimited, cheap time resolution
over weeks. The weekly photographs of item 1 recover most of that for free.
Photographs are quantitative enough IF they are taken under fixed lighting with
a grey reference in frame, so transmittance can be read from pixel intensity
afterwards; without that discipline they are illustration, not data.

So: skip the dedicated turbidimeter, but **do not skip the long-time
observation** in some form.

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

## If the 12-ID experiment works

A realistic reading of where it lands the paper:

| outcome | venue |
|---|---|
| Recovery curves per chain length, beta-sheet peaks tracking them | **PNAS or Nano Letters** is arguable: it is a mechanism, not a correlation |
| Recovery curves, beta-sheet peaks not tracking | **ACS Nano / Nature Communications**: still a real result, and knowing the two are not linked is worth publishing |
| Endpoints only, no kinetics | **Biomacromolecules**, as now, but with the photograph replaced by data |

The risk worth naming: if the beta-sheet peaks do **not** track recovery, both
papers lose the framing they currently lean on. That is a reason to do the
measurement, not to avoid it, but it is better to know before the framing is
built into a submitted abstract.
