## Table S7. Precedence under four definitions of the background rate

The background rate the observed count is judged against can be drawn in more
than one way, and the choice is not neutral. The rate as published counts every
membrane position in the protein, including positions whose own residue is the
target; that inflates the denominator and biases each ratio toward depletion.
The self-exclusive rate counts only positions whose residue differs from the
target, which places the substitutable positions and the background on the same
footing. The two matched rates are computed separately within strata of observed
alphabet size, and within strata of alphabet size crossed with conservation
grade, so that a site is compared only against positions of comparable
alignment breadth. Alphabet size and conservation grade are derived from the
same alignment column and are strongly correlated (Spearman ρ = −0.78 for
EAA1, −0.80 for VPAT), so the doubly matched estimate is a conservative bound
rather than an independent control. *p* values are two-sided and adjusted across
the sixteen tests by the Benjamini–Hochberg procedure. The number of
sites and the number already carrying the target do not depend on the
background definition: L→Q 7 of 87, I→T 33 of 56, V→T 48 of 64 and
F→Y 21 of 37. Only the expectation, and so the ratio, moves between
definitions.

| Background | Substitution | Expected | Ratio | 95% CI | *p* | *p* adj. |
| --- | --- | --- | --- | --- | --- | --- |
| As published | L→Q | 9.97 | 0.70 | 0.13–1.28 | 0.31 | 0.45 |
| As published | I→T | 28.90 | 1.14 | 0.89–1.39 | 0.27 | 0.43 |
| As published | V→T | 32.45 | 1.48 | 1.24–1.72 | 8.4 × 10⁻⁵ | 2.2 × 10⁻⁴ |
| As published | F→Y | 7.07 | 2.97 | 2.31–3.63 | 5.6 × 10⁻⁹ | 2.2 × 10⁻⁸ |
| Self-exclusive | L→Q | 9.42 | 0.74 | 0.15–1.34 | 0.40 | 0.48 |
| Self-exclusive | I→T | 27.14 | 1.22 | 0.95–1.49 | 0.11 | 0.20 |
| Self-exclusive | V→T | 30.56 | 1.57 | 1.32–1.82 | 1.1 × 10⁻⁵ | 3.6 × 10⁻⁵ |
| Self-exclusive | F→Y | 6.47 | 3.25 | 2.55–3.95 | 3.1 × 10⁻¹⁰ | 1.7 × 10⁻⁹ |
| Alphabet matched | L→Q | 9.62 | 0.73 | 0.16–1.29 | 0.35 | 0.46 |
| Alphabet matched | I→T | 30.38 | 1.09 | 0.87–1.30 | 0.42 | 0.48 |
| Alphabet matched | V→T | 35.76 | 1.34 | 1.16–1.53 | 3.3 × 10⁻⁴ | 6.7 × 10⁻⁴ |
| Alphabet matched | F→Y | 5.93 | 3.54 | 2.82–4.26 | 5.0 × 10⁻¹² | 8.1 × 10⁻¹¹ |
| Alphabet and grade | L→Q | 8.96 | 0.78 | 0.21–1.35 | 0.45 | 0.48 |
| Alphabet and grade | I→T | 33.02 | 1.00 | 0.81–1.18 | 0.99 | 0.99 |
| Alphabet and grade | V→T | 35.81 | 1.34 | 1.16–1.52 | 1.7 × 10⁻⁴ | 4.0 × 10⁻⁴ |
| Alphabet and grade | F→Y | 6.35 | 3.31 | 2.62–4.00 | 4.8 × 10⁻¹¹ | 3.8 × 10⁻¹⁰ |

Background definitions: *As published*, every membrane position counted;
*Self-exclusive*, positions whose own residue is the target omitted;
*Alphabet matched*, strata of observed alphabet size; *Alphabet and grade*,
strata of alphabet size crossed with conservation grade.

Both enrichments survive every definition. F→Y strengthens as matching is
added, from 2.97 to 3.31; V→T attenuates from 1.48 to 1.34 while remaining
significant by three orders of magnitude, so part of its excess is attributable
to alignment breadth and part is not. Neither null result becomes significant
under any definition, and the apparent depletion of L→Q moves toward unity as
matching is added.

## Table S8. Precedence within strata of observed alphabet size

A position whose alignment column admits many residues is more likely to
admit any particular target, so precedence rises with alphabet breadth for
purely arithmetic reasons. Each stratum below is judged against the background
rate computed within that same stratum. Strata containing no substitutable
position, and strata in which the background rate is identically zero, are
omitted. *p* values are adjusted across the twenty tests.

| Substitution | Alphabet size | Sites | Precedent / expected | Ratio | 95% CI | *p* | *p* adj. |
| --- | --- | --- | --- | --- | --- | --- | --- |
| L→Q | 4-5 | 24 | 2 / 1.98 | 1.01 | 0.00–2.35 | 0.99 | 0.99 |
| L→Q | 6-7 | 32 | 0 / 2.20 | 0.00 | 0.00–1.27 | 0.12 | 0.29 |
| L→Q | 8-9 | 18 | 2 / 2.50 | 0.80 | 0.00–1.93 | 0.73 | 0.91 |
| L→Q | 10+ | 8 | 3 / 2.95 | 1.02 | 0.12–1.92 | 0.97 | 0.99 |
| I→T | 1-3 | 3 | 0 / 0.13 | 0.00 | 0.00–5.41 | 0.72 | 0.91 |
| I→T | 4-5 | 11 | 3 / 2.48 | 1.21 | 0.12–2.31 | 0.71 | 0.91 |
| I→T | 6-7 | 20 | 13 / 10.96 | 1.19 | 0.79–1.58 | 0.36 | 0.62 |
| I→T | 8-9 | 15 | 11 / 10.59 | 1.04 | 0.71–1.36 | 0.82 | 0.91 |
| I→T | 10+ | 7 | 6 / 6.22 | 0.96 | 0.70–1.23 | 0.79 | 0.91 |
| V→T | 1-3 | 2 | 0 / 0.10 | 0.00 | 0.00–6.04 | 0.75 | 0.91 |
| V→T | 4-5 | 15 | 6 / 3.20 | 1.87 | 0.90–2.84 | 0.078 | 0.25 |
| V→T | 6-7 | 20 | 16 / 10.89 | 1.47 | 1.07–1.87 | 0.022 | 0.10 |
| V→T | 8-9 | 14 | 13 / 9.88 | 1.32 | 0.98–1.65 | 0.067 | 0.25 |
| V→T | 10+ | 13 | 13 / 11.69 | 1.11 | 0.93–1.29 | 0.23 | 0.47 |
| F→Y | 1-3 | 9 | 5 / 0.97 | 5.14 | 3.27–7.01 | 1.5 × 10⁻⁵ | 1.7 × 10⁻⁴ |
| F→Y | 4-5 | 7 | 2 / 0.71 | 2.82 | 0.61–5.02 | 0.11 | 0.29 |
| F→Y | 6-7 | 14 | 8 / 2.18 | 3.67 | 2.45–4.88 | 1.8 × 10⁻⁵ | 1.7 × 10⁻⁴ |
| F→Y | 8-9 | 4 | 4 / 0.96 | 4.18 | 2.44–5.93 | 3.5 × 10⁻⁴ | 0.002 |
| F→Y | 10+ | 3 | 2 / 1.11 | 1.80 | 0.33–3.28 | 0.29 | 0.54 |

## Table S9. Precedence within strata of conservation grade

ConSurf grades run from 1 to 9, with 9 the most conserved. As above, each
stratum is judged against the background rate within that stratum, so a
conserved position is compared only against other conserved positions. *p*
values are adjusted across the sixteen tests.

| Substitution | Grade | Sites | Precedent / expected | Ratio | 95% CI | *p* | *p* adj. |
| --- | --- | --- | --- | --- | --- | --- | --- |
| L→Q | 1-3 variable | 15 | 0 / 0.71 | 0.00 | 0.00–2.04 | 0.34 | 0.58 |
| L→Q | 4-6 intermediate | 31 | 5 / 5.22 | 0.96 | 0.20–1.72 | 0.91 | 1.00 |
| L→Q | 7-8 conserved | 32 | 1 / 2.67 | 0.38 | 0.00–1.52 | 0.29 | 0.58 |
| L→Q | 9 most conserved | 9 | 1 / 1.27 | 0.79 | 0.00–2.36 | 0.79 | 0.99 |
| I→T | 1-3 variable | 5 | 4 / 4.00 | 1.00 | 0.60–1.40 | 1.00 | 1.00 |
| I→T | 4-6 intermediate | 29 | 21 / 20.86 | 1.01 | 0.79–1.23 | 0.95 | 1.00 |
| I→T | 7-8 conserved | 21 | 8 / 8.73 | 0.92 | 0.42–1.42 | 0.74 | 0.99 |
| V→T | 1-3 variable | 12 | 11 / 9.87 | 1.11 | 0.88–1.35 | 0.35 | 0.58 |
| V→T | 4-6 intermediate | 26 | 21 / 17.68 | 1.19 | 0.93–1.44 | 0.15 | 0.45 |
| V→T | 7-8 conserved | 23 | 15 / 9.06 | 1.66 | 1.16–2.15 | 0.01 | 0.036 |
| V→T | 9 most conserved | 3 | 1 / 0.55 | 1.83 | 0.00–4.23 | 0.50 | 0.74 |
| F→Y | 1-3 variable | 2 | 1 / 0.42 | 2.36 | 0.00–5.03 | 0.32 | 0.58 |
| F→Y | 4-6 intermediate | 14 | 10 / 3.54 | 2.82 | 1.93–3.72 | 6.7 × 10⁻⁵ | 5.1 × 10⁻⁴ |
| F→Y | 7-8 conserved | 14 | 5 / 1.68 | 2.98 | 1.56–4.40 | 0.006 | 0.031 |
| F→Y | 9 most conserved | 7 | 5 / 0.78 | 6.43 | 4.33–8.52 | 3.8 × 10⁻⁷ | 5.7 × 10⁻⁶ |

F→Y is enriched in three of the four conservation strata after
adjustment, and most strongly in the most conserved stratum, where five of
seven phenylalanine columns at grade 9 already contain tyrosine against an
expectation of 0.78. The enrichment is therefore not an artefact of permissive,
poorly conserved columns; it is largest exactly where the alignment is most
constrained.

## Table S10. Every possible target residue ranked at each native position

Each of the nineteen alternatives to the native residue was scored by the same
background-controlled procedure, over the same membrane positions, giving
seventy-six estimates. The five best-supported targets are listed for each
native residue, together with the QTY target where it falls outside the top
five. *p* values are adjusted across all seventy-six tests. The complete table
is Supplementary Data S2.

| Native | Substitution | Rank | Sites | With precedent | Expected | Ratio | *p* adj. |
| --- | --- | --- | --- | --- | --- | --- | --- |
| L | L→M | 1 of 19 | 87 | 63 | 42.89 | 1.47 | 2.2 × 10⁻⁴ |
| L | L→I | 2 of 19 | 87 | 74 | 53.47 | 1.38 | 1.1 × 10⁻⁴ |
| L | L→F | 3 of 19 | 87 | 56 | 42.48 | 1.32 | 0.017 |
| L | L→V | 4 of 19 | 87 | 78 | 60.52 | 1.29 | 4.3 × 10⁻⁴ |
| L | L→W | 5 of 19 | 87 | 8 | 8.26 | 0.97 | 0.97 |
| L | L→Q † | 12 of 19 | 87 | 7 | 9.97 | 0.70 | 0.46 |
| I | I→V | 1 of 19 | 56 | 56 | 39.44 | 1.42 | 3.0 × 10⁻⁵ |
| I | I→L | 2 of 19 | 56 | 52 | 37.05 | 1.40 | 2.6 × 10⁻⁴ |
| I | I→F | 3 of 19 | 56 | 37 | 26.74 | 1.38 | 0.024 |
| I | I→M | 4 of 19 | 56 | 38 | 28.66 | 1.33 | 0.046 |
| I | I→T † | 5 of 19 | 56 | 33 | 28.90 | 1.14 | 0.41 |
| V | V→I | 1 of 19 | 64 | 62 | 39.72 | 1.56 | 3.4 × 10⁻⁷ |
| V | V→T † | 2 of 19 | 64 | 48 | 32.45 | 1.48 | 7.1 × 10⁻⁴ |
| V | V→L | 3 of 19 | 64 | 53 | 42.26 | 1.25 | 0.021 |
| V | V→C | 4 of 19 | 64 | 27 | 21.67 | 1.25 | 0.32 |
| V | V→M | 5 of 19 | 64 | 39 | 32.33 | 1.21 | 0.22 |
| F | F→Y † | 1 of 19 | 37 | 21 | 7.07 | 2.97 | 3.4 × 10⁻⁷ |
| F | F→W | 2 of 19 | 37 | 10 | 3.53 | 2.83 | 0.002 |
| F | F→L | 3 of 19 | 37 | 25 | 24.31 | 1.03 | 0.88 |
| F | F→V | 4 of 19 | 37 | 24 | 25.63 | 0.94 | 0.66 |
| F | F→M | 5 of 19 | 37 | 16 | 18.00 | 0.89 | 0.60 |

† the substitution specified by the QTY code.

At phenylalanine positions the QTY target is the best-supported of all nineteen
alternatives, with tryptophan second; the two aromatic targets are the only ones
enriched, while glycine, alanine and serine are significantly depleted, so the
alignment at these positions discriminates by side-chain bulk and aromaticity
rather than by polarity alone. At valine positions threonine ranks second, behind
isoleucine alone,
and is the best-supported polar target. At isoleucine positions threonine ranks
fifth, behind four hydrophobic targets, and is again the best-supported polar
one. At leucine positions the four best-supported targets are all hydrophobic
and no polar target is enriched; glutamine ranks twelfth. Threonine ranks ninth
at leucine positions and is present in far more leucine columns than glutamine
is (32 against 7 of the same 87 columns), but that difference reflects the much
higher overall abundance of threonine in membrane spans, and after background
control the two are indistinguishable (0.75 and 0.70, overlapping intervals).

## Table S11. Concordance with general amino-acid exchangeability

The seventy-six estimates of Table S10 were correlated against five
general-purpose substitution matrices. None is specific to membrane proteins,
which is the point of the comparison: it asks how much of the pattern observed
inside the bilayer is simply general exchangeability.

| Matrix | Pairs | Spearman ρ | *p* |
| --- | --- | --- | --- |
| BLOSUM62 | 76 | 0.83 | 1.1 × 10⁻²⁰ |
| PAM250 | 76 | 0.80 | 2.8 × 10⁻¹⁸ |
| GONNET1992 | 76 | 0.88 | 1.4 × 10⁻²⁵ |
| JTT (Jones–Taylor–Thornton) | 76 | 0.81 | 5.7 × 10⁻¹⁹ |
| RISLER | 76 | 0.41 | 1.9 × 10⁻⁴ |

Precedence inside membrane spans tracks general exchangeability closely.
Positioning each QTY substitution as a percentile of the seventy-six estimates,
on the precedence scale and on each matrix scale, isolates the part that general
exchangeability does not account for.

| Substitution | Ratio | Precedence percentile | Mean matrix percentile | Difference |
| --- | --- | --- | --- | --- |
| F→Y | 2.97 | 100 | 98 | +2 |
| I→T | 1.14 | 80 | 77 | +3 |
| L→Q | 0.70 | 43 | 54 | −11 |
| V→T | 1.48 | 96 | 79 | +17 |

F→Y and I→T sit where general exchangeability predicts. V→T is
better precedented inside membrane spans than any of the five matrices
anticipates, and L→Q slightly worse, which is the expected direction if the
bilayer penalises the introduction of an amide more than a general-purpose
matrix does.

## Table S12. Leave-one-segment-out robustness

Each of the twenty-two annotated membrane segments was removed in turn and the
pooled estimate recomputed from the remaining segments, with the background rate
recomputed over the same reduced scope. No single segment carries any of the
four results.

| Substitution | Segments dropped | Full estimate | Minimum | Attained on dropping | Maximum | Attained on dropping | Drops with *p* < 0.05 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| L→Q | 22 | 0.70 | 0.48 | EAA1 TM7 | 0.78 | EAA1 TM5 | 0 of 22 |
| I→T | 22 | 1.14 | 1.08 | EAA1 TM3 | 1.20 | VPAT TM10 | 0 of 22 |
| V→T | 22 | 1.48 | 1.42 | VPAT TM6 | 1.54 | EAA1 TM9 | 22 of 22 |
| F→Y | 22 | 2.97 | 2.69 | VPAT TM3 | 3.26 | EAA1 TM5 | 22 of 22 |

## Table S13. Precedence by position within the membrane segment

Each membrane position was assigned a relative depth, zero at the midpoint of
its segment and one at either end, and split at the median into the inner half
of the segment and the outer half abutting the headgroup region. Background
rates are computed separately within each zone. *p* values are adjusted across
the eight tests.

| Substitution | Zone | Sites | Precedent / expected | Ratio | 95% CI | *p* | *p* adj. |
| --- | --- | --- | --- | --- | --- | --- | --- |
| L→Q | Core | 40 | 2 / 2.96 | 0.68 | 0.00–1.75 | 0.55 | 0.55 |
| L→Q | Interfacial | 47 | 5 / 6.69 | 0.75 | 0.06–1.44 | 0.47 | 0.55 |
| I→T | Core | 24 | 16 / 11.55 | 1.39 | 0.97–1.80 | 0.069 | 0.11 |
| I→T | Interfacial | 32 | 17 / 15.30 | 1.11 | 0.76–1.47 | 0.54 | 0.55 |
| V→T | Core | 28 | 23 / 13.32 | 1.73 | 1.34–2.12 | 2.5 × 10⁻⁴ | 9.8 × 10⁻⁴ |
| V→T | Interfacial | 36 | 25 / 17.17 | 1.46 | 1.12–1.79 | 0.008 | 0.016 |
| F→Y | Core | 10 | 4 / 1.07 | 3.75 | 1.96–5.54 | 0.003 | 0.007 |
| F→Y | Interfacial | 27 | 17 / 6.45 | 2.64 | 1.96–3.31 | 1.9 × 10⁻⁶ | 1.5 × 10⁻⁵ |

Both enrichments are present in both zones and are somewhat larger in the
interior of the segment, where a polar substitution would be least expected to
be tolerated.

## Table S14. Precedence and predicted burial

ConSurf classifies positions as buried or exposed from the predicted relative
solvent accessibility of the query structure. The classification is reported
only at the substitutable positions, so the comparison below is between buried
and exposed positions carrying the same native residue, by Fisher's exact test,
rather than against a burial-matched background. *n* is the number of positions and *k* the number whose alignment column already contains the target. Twenty-six of the 351 pooled
positions carry no call and are excluded.

| Substitution | Buried, *k* / *n* (rate) | Exposed, *k* / *n* (rate) | Odds ratio | *p* | *p* adj. |
| --- | --- | --- | --- | --- | --- |
| L→Q | 5 / 71 (0.07) | 2 / 16 (0.12) | 0.53 | 0.61 | 0.81 |
| I→T | 19 / 37 (0.51) | 14 / 19 (0.74) | 0.38 | 0.15 | 0.31 |
| V→T | 31 / 45 (0.69) | 17 / 19 (0.90) | 0.26 | 0.12 | 0.31 |
| F→Y | 16 / 28 (0.57) | 5 / 9 (0.56) | 1.07 | 1.00 | 1.00 |

No substitution shows an association between precedence and predicted
burial after adjustment.

## Table S15. The same four pairs read in the opposite direction

Precedence for the reverse substitution was computed at positions carrying the
QTY target residue. The two directional estimates need not agree even for a
freely exchangeable pair, because each conditions on the query residue rather
than on membership of the alignment column; the symmetric lift in the final
column removes that asymmetry by using column membership on both sides. The
minimum detectable ratio is the smallest enrichment the available sites could
have resolved at 95% confidence.

| Pair | Sites | Precedent / expected | Reverse ratio (95% CI) | *p* adj. | Forward ratio | Lift | Min. detect. |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Q→L | 4 | 1 / 2.66 | 0.38 (0.00–1.07) | 0.16 | 0.70 | 0.93 | 1.70 |
| T→I | 29 | 18 / 18.45 | 0.98 (0.70–1.25) | 0.86 | 1.14 | 1.22 | 1.27 |
| T→V | 29 | 20 / 20.75 | 0.96 (0.73–1.19) | 0.86 | 1.48 | 1.19 | 1.23 |
| Y→F | 10 | 9 / 4.84 | 1.86 (1.22–2.50) | 0.034 | 2.97 | 1.45 | 1.64 |

Phenylalanine is enriched at tyrosine positions as tyrosine is at
phenylalanine positions, so the aromatic pair is exchangeable in both
directions. The threonine reversals are unremarkable, and the glutamine
reversal rests on four positions and could not have resolved an enrichment
below 1.70; it is reported for completeness and carries no weight.
