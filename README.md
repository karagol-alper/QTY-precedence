# Evolutionary precedent for the QTY code

Data and analysis code for *Unequal evolutionary standing of the four
QTY-code substitutions in two human membrane transporters* (Karagöl and
Karagöl).

The QTY code water-solubilises membrane proteins by replacing leucine with
glutamine, isoleucine and valine with threonine, and phenylalanine with
tyrosine throughout transmembrane segments, and is applied as a single uniform
rule. This analysis asks, position by position, whether evolution has already
sampled the residue the code installs there.

Every reported number is produced by the code in this repository from the
ConSurf outputs in `data/`. Nothing is entered by hand.

## Headline result

Within UniProt-annotated membrane segments of EAA1/SLC1A3 and VPAT/SLC18B1,
pooled over 244 substitutable positions:

| Substitution | Sites | With precedent | Ratio | 95% CI | *p* |
| --- | --- | --- | --- | --- | --- |
| F→Y | 37 | 21 | 2.97 | 2.31–3.63 | 5.6 × 10⁻⁹ |
| V→T | 64 | 48 | 1.48 | 1.24–1.72 | 8.4 × 10⁻⁵ |
| I→T | 56 | 33 | 1.14 | 0.89–1.39 | 0.27 |
| L→Q | 87 | 7 | 0.70 | 0.13–1.28 | 0.31 |

The ratio is the observed rate divided by the same protein's background rate
for the same target residue over the same sequence scope. A ratio above 1
means evolution already favours the substitution at the positions the code
rewrites.

## Layout

```
data/        ConSurf residue-variety tables, one row per position
analysis/    the analysis code
notebooks/   qty_precedence.ipynb, the primary notebook (Colab + Drive)
results/     the analysis run reproduced in the paper
supplementary/  the supplementary data files and figures
```

### `data/`

| File | Positions | MSA depth | Role |
| --- | --- | --- | --- |
| `eaa1_residue_variety.csv` | 542 | 150 | EAA1/SLC1A3, pooled |
| `eaa1_residue_variety_std.csv` | 457 | 150 | EAA1 SEQRES 43–499, with conservation grades |
| `vpat_residue_variety.csv` | 456 | 50 | VPAT/SLC18B1, pooled |
| `eaa2_residue_variety.csv` | 459 | 150 | EAA2/SLC1A2, reference, QTY-designed |
| `tlr5_residue_variety.csv` | 644 | 50 | TLR5S, reference, naturally soluble |
| `eaa1_consurf_43_499.py` | 457 | 150 | independent transcription of the same ConSurf job, used to cross-check the grades |

The run's EAA1 input carries no conservation grade, so grades are taken from
`eaa1_residue_variety_std.csv`. That file is identical to the run input
residue for residue over SEQRES 43–499, at the same MSA depth, and contains
the whole EAA1 membrane span (48–494); `core.load_tables()` asserts the
agreement rather than assuming it. Taking grades from the per-site table
instead would restrict them to L/I/V/F positions and make every
conservation-matched background hydrophobic-only.

Columns are `seqres_pos`, `aa`, `n_variety`, `variety` (the observed alphabet,
pipe-separated), and where present `grade` (1–9, 9 most conserved), `score`
(lower is more conserved) and `msa_total`.

### Provenance and the manifest

`results/run_2026-11-04_161147/run_manifest.json` records every input with its
SHA-256 digest, the alignment depths, the membrane boundaries and where they
came from, the full parameter set and the library versions, so any reported
number is traceable to the ConSurf file that produced it.

The absolute paths of those inputs have been reduced to bare filenames by
`analysis/redact_manifest.py`: they described the directory layout of the
machine that produced the run, including unrelated work, and the digest is
what anchors provenance. Every digest, count and parameter is untouched. To
confirm an input is the file the run used, hash it and compare against
`inputs[].sha256`.

### `analysis/`

| Script | What it does |
| --- | --- |
| `core.py` | loading, membrane masking, background rates, the pooled Poisson-binomial estimator, BH adjustment |
| `check_repro.py` | rebuilds the four primary ratios from the raw tables and compares them against the published run |
| `analyses.py` | the extended analyses, writing `results/extended/S*.csv` |
| `make_fig1.py` | Figure 1, the forest plot |
| `make_fig2.py` | Figure 2, the residue ranking at EAA1 leucine positions |
| `make_figS1.py` | Figure S1, the confound-control panels |
| `gen_supp.py` | emits the markdown for supplementary Tables S7–S15 from those CSVs |
| `verify.py` | checks every numeric claim in the paper against the CSVs |
| `diag.py` | covariate collinearity, stratum sizes, symmetric co-occurrence lift |
| `redact_manifest.py` | reduces absolute input paths in a manifest to filenames |
| `tmfetch.py` | retrieves UniProt membrane features and reconciles them against the ConSurf query |
| `build_nb.py` | regenerates the notebook with the data embedded |

## Reproducing

```bash
pip install -r requirements.txt
./run_all.sh
```

`run_all.sh` runs the reproduction check first and stops if it fails, since
nothing downstream is trustworthy if the primary estimates cannot be rebuilt.
It then recomputes the extended analyses into `results/extended/`, redraws
all three figures, regenerates the supplementary tables, and runs the
verification script. Every figure and every number in the paper is
produced by this one command.

To analyse a different run or a different set of ConSurf outputs without
editing the source:

```bash
QTY_RUN=/path/to/run QTY_DATA=/path/to/consurf python3 analysis/analyses.py
```

`verify.py` runs its numeric checks against the CSVs always. Its text checks
compare the manuscript against those numbers and are skipped unless the
manuscript sources are present at `manuscript/` or at `$QTY_MANUSCRIPT`.

### The notebook

`notebooks/qty_precedence.ipynb` is the primary analysis as run for the paper.
It mounts Google Drive, reads ConSurf `consurf_grades.txt` files from a
configured directory, and writes a timestamped run folder containing every
table, the figures, and a manifest recording each input with its SHA-256 hash,
the alignment depths, the membrane boundaries with their provenance, the
parameter set and the library versions. The four ConSurf tables are also
embedded in the notebook, so it runs and produces the published numbers
without Drive access; when it falls back to the embedded copies it says so in
the output, in the manifest, and in the figure title.

## Method in brief

For each native residue the code targets, the observed rate is the fraction of
positions carrying that residue whose alignment column already contains the
QTY target. The background is the fraction of all scored positions in the same
protein, over the same scope, whose column contains that target. Dividing by a
within-protein, within-scope background makes the estimate independent of
alignment depth, of composition, and of which regions are included.

Scope matters more than it might appear. Glutamine availability falls two- to
three-fold between whole sequence and membrane span in both proteins, while
threonine and tyrosine availability does not, so a whole-sequence ratio for
L→Q compares membrane-resident leucines against a background drawn largely
from soluble domains and reports depletion for a reason unrelated to the
positions being substituted.

Across proteins, counts are summed and compared against the sum of each
protein's own expectation — a Poisson-binomial summarised by its mean and
variance and evaluated as a normal deviate. With two proteins, Cochran's *Q*
carries one degree of freedom and *I*² is uninterpretable, so no meta-analytic
estimate is reported.

## Citation

Karagöl A, Karagöl T. Unequal evolutionary standing of the four QTY-code
substitutions in two human membrane transporters. *Manuscript in preparation.*

## License

No license has been applied yet. Add one before making the repository public:
without a license, others have no right to use or redistribute the contents.
