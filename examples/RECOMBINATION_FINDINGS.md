# Within-Athila recombination — findings

Athila families (ATHILA0–9) are near-identical at the PROTEIN level (best/2nd family
bitscore ~1.05×), so recombination is only detectable at the NUCLEOTIDE level
(per-domain family typing vs TAIR12 ATHILA internal exemplars; INT margin ~1.68× nt).

Recombination always sits at the ELEMENT ENDS (GAG or INT), never the RT/RH core.

## Full-5 screen (complete elements, all 5 domains present — the strict screen)
- 3,825 full-5 Athila elements; 2,143 with all 4 long domains confidently nt-typed.
- **2 recombinants**, both `GAG:ATHILA7a ; PROT/RT/RH/INT:ATHILA7` (Toufl-1, Elh-2):
  a clean 5′ GAG swap between sister families ATHILA7/ATHILA7a.
  Figure: `recomb_Toufl1_full5_allfamilies.png` (INT present; all families shown).

## Relaxed screen (>=2 confident long domains; includes INT-degraded elements)
- **128 recombinants.** Odd domain: GAG 104, INT 24 (never RT/RH).
- Dominant class: `GAG:ATHILA2 ; RT/RH:ATHILA4c`, recurrent in **80 copies, all Chr3
  16–20 Mb** (shared/ancestral) — but these lack an intact INT, so they are NOT full-5.
  Figure: `recomb_9994_allfamilies.png`.

Trade-off: requiring full-5 is clean but discards degraded-but-real recombinants
(old TEs commonly lose INT). The strict full-5 signal is sparse (2); the recurrent
Chr3 ATHILA2/4c class is the strongest population-level signal but sub-full-5.
