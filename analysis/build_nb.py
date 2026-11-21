#!/usr/bin/env python3
"""Generate qty_precedence.ipynb — Colab + Drive, self-contained, extensible."""
import json, gzip, base64, uuid, pathlib

DATA = pathlib.Path(__file__).resolve().parent.parent / "data"

blobs = {}
for key, name in [("EAA1 / SLC1A3", "eaa1_residue_variety_std.csv"),
                  ("VPAT / SLC18B1", "vpat_residue_variety.csv"),
                  ("EAA2 / SLC1A2 (QTY)", "eaa2_residue_variety.csv"),
                  ("TLR5S (soluble)", "tlr5_residue_variety.csv")]:
    blobs[key] = base64.b64encode(
        gzip.compress((DATA / name).read_bytes(), 9)).decode()

def _id(): return uuid.uuid4().hex[:8]
def md(src):  return {"cell_type": "markdown", "id": _id(), "metadata": {}, "source": src}
def code(src): return {"cell_type": "code", "id": _id(), "execution_count": None,
                       "metadata": {}, "outputs": [], "source": src}

cells = []

cells.append(md(r"""# Evolutionary precedent for the QTY code

For every site the QTY code rewrites, does the homolog alignment already sample the
residue the code puts there?

The code maps **L&rarr;Q, I&rarr;T, V&rarr;T, F&rarr;Y** inside transmembrane segments.
ConSurf reports, per position, the amino acids actually observed across the alignment
(its `RESIDUE VARIETY` column), which is enough to ask &mdash; per site &mdash; whether
the substitution is one evolution has already made there. Every rate is divided by that
protein's **own** background rate for the same target residue, so alignments of
different depth stay comparable.

The question can only be asked of a **native membrane protein**: a sequence that already
carries the QTY substitutions has no native residue left at the sites of interest, and a
soluble protein has no transmembrane segment for the code to act on. Section 4 classifies
every input from its composition and hydropathy and pools only the native membrane ones.

**Run order:** mount Drive &rarr; set paths &rarr; run all. Results are written to a new
timestamped folder in Drive, so nothing is overwritten between runs and nothing is lost
when the Colab VM is recycled."""))

cells.append(md("""## 0 · Mount Drive

Colab only. Running anywhere else is fine — the cell detects it, skips the mount, and
the notebook falls back to bundled copies of the four inputs.

A first mount that raises `ValueError` is almost always a stale or half-finished mount
left by an earlier session, so the cell retries once with `force_remount=True` before
giving up. **If it still fails the notebook analyses bundled demo data, not yours** —
which produces a complete, plausible set of numbers that are not your result. Every later
cell is marked when that happens; do not read a run that carries the mark."""))

cells.append(code(r'''import sys, os, io, re, gzip, json, base64, hashlib, shutil, datetime
from pathlib import Path

IN_COLAB = "google.colab" in sys.modules or os.path.exists("/content")
DRIVE_OK = False

if IN_COLAB:
    try:
        from google.colab import drive
        try:
            drive.mount("/content/drive", force_remount=False)
        except ValueError as exc:
            # usually a stale mountpoint from a previous session
            print(f"first mount attempt failed ({exc}); retrying with force_remount=True")
            drive.mount("/content/drive", force_remount=True)
        DRIVE_OK = Path("/content/drive/MyDrive").is_dir()
    except Exception as exc:                       # declined auth, or not really Colab
        print(f"Drive mount unavailable ({type(exc).__name__}: {exc})")
else:
    print("Not in Colab — skipping Drive mount.")

DRIVE = Path("/content/drive/MyDrive") if DRIVE_OK else None
print(f"Drive mounted: {DRIVE_OK}")
if IN_COLAB and not DRIVE_OK:
    print("\n" + "!" * 74)
    print("Drive is NOT mounted. The notebook will fall back to BUNDLED DEMO DATA and")
    print("results will be written to the VM, which is discarded when it recycles.")
    print("Fix the mount and re-run from the top before using any number below.")
    print("!" * 74)'''))

cells.append(md(r"""## 1 · Paths

`DATA_ROOTS` are searched **recursively** for `*consurf_grades*.txt`. The defaults are
the two branches that hold the ConSurf runs; add or replace freely. Point them at a
tight folder rather than at the whole of `MyDrive` &mdash; the Drive FUSE mount is slow
to walk, and a broad search picks up every stale copy.

Results go to a fresh `run_<timestamp>` folder under `RESULTS_ROOT`."""))

cells.append(code(r'''# ---- where the ConSurf runs live -------------------------------------------
if DRIVE_OK:
    # Folders to search for ConSurf output. Point these at your own run
    # directories; each is searched recursively for consurf_grades.txt and
    # no_model_consurf_grades.txt. Leave the list empty to use the ConSurf
    # tables bundled in this notebook, which reproduce the published numbers.
    DATA_ROOTS = [
        DRIVE / "QTY_precedence/consurf",
    ]
    RESULTS_ROOT = DRIVE / "QTY_precedence"
else:
    DATA_ROOTS = []                 # no Drive -> bundled data
    RESULTS_ROOT = Path("qty_precedence")

# ---- a new, timestamped results folder per run ------------------------------
RUN_ID = datetime.datetime.now().strftime("run_%Y-%m-%d_%H%M%S")
RESULTS_DIR = RESULTS_ROOT / RUN_ID
TABLES_DIR, FIGURES_DIR = RESULTS_DIR / "tables", RESULTS_DIR / "figures"
for d in (TABLES_DIR, FIGURES_DIR):
    d.mkdir(parents=True, exist_ok=True)

# ---- the substitution rule under test ---------------------------------------
QTY_MAP = {"L": "Q", "I": "T", "V": "T", "F": "Y"}

# Folder names are what identify a protein (every file is called consurf_grades.txt),
# and some are opaque timestamps. Map them to readable labels here; anything not
# listed keeps its folder name. Names in EXCLUDE are dropped.
NAME_MAP = {
    "EAA1plos":   "EAA1 / SLC1A3",
    "VPATplos":   "VPAT / SLC18B1",
    "EAA2QTY":    "EAA2 / SLC1A2 (QTY)",
    "TLR5sarxiv": "TLR5S (soluble)",
}
# Neither of the last two is a native membrane protein, and the classifier at the end of
# section 4 will say so from composition alone:
#   EAA2QTY    is the QTY-DESIGNED sequence — the substitutions are already in it.
#   TLR5sarxiv is TLR5S, a separate teleost gene encoding the LRR ectodomain alone:
#              no transmembrane segment, no TIR domain, secreted. A membrane receptor
#              evolution has already solubilised — by deleting the anchor, not by
#              polarising it.
# Both are measured and reported; neither is pooled.
EXCLUDE = set()

for r in DATA_ROOTS:
    print(f"{'ok ' if r.is_dir() else 'MISSING'}  {r}")
print(f"\nresults -> {RESULTS_DIR}")'''))

cells.append(md(r"""## 2 · ConSurf parser

Handles both layouts ConSurf emits:

* **with-model** &mdash; has a `3LATOM` column (`SER:1:A`), variety is bare letters
* **`no_model`** &mdash; no `3LATOM`, variety entries may carry frequencies (`D 80%, E 20%`)

It also absorbs markdown-escaped asterisks (`8\*`), which appear when a file is pulled
through a text-rendering API rather than downloaded raw."""))

cells.append(code(r'''ATOM_RE = re.compile(r"^[A-Z]{3}:-?\d+:\w+$")
MSA_RE  = re.compile(r"(\d+)\s*/\s*(\d+)")
NUM_RE  = re.compile(r"^-?\d+\.?\d*$")


def parse_consurf(text: str) -> "pd.DataFrame":
    """ConSurf grades file -> tidy DataFrame.

    Columns: seqres_pos, aa, atom_resid, score, grade, low_conf,
             msa_n, msa_total, buried, functional, n_variety, variety
    `variety` is pipe-delimited so it survives CSV round-trips.
    """
    rows = []
    for raw in text.splitlines():
        line = raw.replace("\\*", "*").replace("\t", " ").strip()
        if not line:
            continue
        tok = line.split()
        if len(tok) < 5 or not tok[0].isdigit() or len(tok[1]) != 1 or not tok[1].isalpha():
            continue
        m = MSA_RE.search(line)
        if not m:
            continue

        pos, aa = int(tok[0]), tok[1]
        head, tail = line[:m.start()].split(), line[m.end():].strip()

        i, atom_resid = 2, None
        if i < len(head) and (ATOM_RE.match(head[i]) or head[i] == "-"):
            if head[i] != "-":
                atom_resid = int(head[i].split(":")[1])
            i += 1
        if i >= len(head) or not NUM_RE.match(head[i]):
            continue
        score = float(head[i]); i += 1
        grade_tok = head[i] if i < len(head) else ""
        low_conf = grade_tok.endswith("*")
        grade = int(grade_tok.rstrip("*"))

        flags = set(head[i + 1:])
        buried = "b" if "b" in flags else ("e" if "e" in flags else "")
        funct  = "f" if "f" in flags else ("s" if "s" in flags else "")

        var = []
        for chunk in tail.split(","):
            c = chunk.strip()
            if c and c[0].isalpha() and c[0].isupper():
                var.append(c[0])
        var = list(dict.fromkeys(var))

        rows.append(dict(seqres_pos=pos, aa=aa, atom_resid=atom_resid,
                         score=score, grade=grade, low_conf=low_conf,
                         msa_n=int(m.group(1)), msa_total=int(m.group(2)),
                         buried=buried, functional=funct,
                         n_variety=len(var), variety="|".join(var)))

    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError("no data rows parsed — is this a ConSurf grades file?")
    df["atom_resid"] = df["atom_resid"].astype("Int64")   # nullable; no_model is all-NA

    # integrity: a native residue must appear in its own alignment column
    bad = (~df.apply(lambda r: r.aa in r.variety.split("|"), axis=1)).sum()
    if bad > 0.02 * len(df):
        raise ValueError(f"parse looks wrong: native residue missing from its own "
                         f"variety list in {bad}/{len(df)} rows")
    df.attrs["n_self_missing"] = int(bad)
    return df'''))

cells.append(md(r"""## 3 · Scientific stack"""))

cells.append(code(r'''import numpy as np
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt

C_ENR, C_DEP, C_NEU = "#2a78d6", "#e34948", "#c6c5be"
mpl.rcParams.update({"figure.dpi": 130, "savefig.dpi": 220, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": False, "savefig.bbox": "tight",
                     "figure.facecolor": "white"})
print(f"pandas {pd.__version__} · numpy {np.__version__}")'''))

cells.append(md(r"""## 4 · Discover and load

Two things this guards against, both of which are live hazards in a real Drive:

1. **Every file is called `consurf_grades.txt`**, so the *parent folder* names the
   protein, not the filename.
2. **The same run is often copied into several folders.** Pooling two copies of one
   protein would silently double-count it, so files are hashed and exact duplicates
   dropped, keeping the shallowest path."""))

cells.append(code(r'''def discover(roots, patterns=("*consurf_grades*.txt",)):
    """Recursively find ConSurf files under roots; drop byte-identical duplicates."""
    hits = []
    for root in roots:
        root = Path(root)
        if not root.is_dir():
            print(f"skip (missing): {root}")
            continue
        for pat in patterns:
            hits.extend(root.rglob(pat))

    seen, keep, dupes = {}, [], []
    for p in sorted(set(hits), key=lambda q: (len(q.parts), str(q))):
        try:
            h = hashlib.sha256(p.read_bytes()).hexdigest()
        except OSError as exc:
            print(f"unreadable, skipped: {p} ({exc})"); continue
        if h in seen:
            dupes.append((p, seen[h]))
        else:
            seen[h] = p; keep.append((p, h))

    if dupes:
        print(f"{len(dupes)} duplicate file(s) ignored (byte-identical to a kept copy):")
        for d, orig in dupes[:10]:
            print(f"   {d}\n     = {orig}")
        if len(dupes) > 10:
            print(f"   ... and {len(dupes) - 10} more")
    return keep


GENERIC = {"RMSD", "ORG", "DATA", "OUTPUT", "OUTPUTS", "RESULTS", "NO_MODEL", "COPY"}
OPAQUE  = re.compile(r"^consurf[_-]?\d+$", re.I)


def label_for(path: Path) -> str:
    """Protein label from the folder tree, since every file has the same name.

    A ConSurf output folder is often an opaque job id (`consurf_1761475466`) or a
    generic bucket (`RMSD`, `org`), so climb until something informative appears.
    """
    if path.parent.name in NAME_MAP:
        return NAME_MAP[path.parent.name]
    for part in path.parent.parts[::-1]:
        if part.upper() in GENERIC or OPAQUE.match(part):
            continue
        return NAME_MAP.get(part, part)
    return path.parent.name


def load_bundled() -> dict:
    out = {}
    for name, b64 in _BUNDLED.items():
        d = pd.read_csv(io.BytesIO(gzip.decompress(base64.b64decode(b64))))
        d["variety"] = d["variety"].astype(str)
        d.attrs["source"] = "bundled"; d.attrs["sha256"] = ""
        out[name] = d
    return out


_BUNDLED = ''' + json.dumps(blobs) + r'''

proteins, manifest, SUFFIXED = {}, [], set()
USING_BUNDLED = False
found = discover(DATA_ROOTS) if DATA_ROOTS else []

for path, sha in found:
    label = label_for(path)
    if label in EXCLUDE:
        print(f"excluded: {label}"); continue
    try:
        d = parse_consurf(path.read_text(errors="replace"))
    except Exception as exc:
        print(f"parse failed, skipped: {path}\n   {type(exc).__name__}: {exc}"); continue
    if label in proteins:
        label = f"{label} ({path.parent.parent.name})"      # disambiguate
        SUFFIXED.add(label)
    d.attrs["source"] = str(path); d.attrs["sha256"] = sha
    proteins[label] = d
    manifest.append(dict(label=label, path=str(path), sha256=sha,
                         n_positions=len(d),
                         msa_total=int(d.msa_total.mode().iloc[0]),
                         n_self_missing=d.attrs["n_self_missing"]))

if not proteins:
    USING_BUNDLED = True
    print("\n" + "!" * 74)
    print("!!  NO ConSurf FILES FOUND. Falling back to BUNDLED DEMO DATA.")
    print("!!  Nothing below describes your files. The EAA1 copy bundled here is the")
    print("!!  43-499 subrange, not the full sequence, so even its numbers differ from")
    print("!!  a real run. Fix DATA_ROOTS (or the Drive mount) and re-run.")
    print("!" * 74)
    proteins = load_bundled()
    manifest = [dict(label=k, path="bundled", sha256="", n_positions=len(v),
                     msa_total=int(v.msa_total.mode().iloc[0]) if "msa_total" in v else -1,
                     n_self_missing=0) for k, v in proteins.items()]

for name, d in proteins.items():
    d["vset"] = d["variety"].str.split("|").apply(set)

# ---- same protein, two ConSurf runs? byte-hashing cannot see this -------------
AUTO_DEDUP_SEQUENCES = True

def _seq(d):
    return "".join(d.sort_values("seqres_pos").aa)

def _depth(d):
    """Modal MSA depth, or -1 when the input does not record one."""
    return int(d.msa_total.mode().iloc[0]) if "msa_total" in d.columns else -1


def _depth_str(d):
    n = _depth(d)
    return "unknown" if n < 0 else str(n)

def _overlap(a, b):
    """Is one sequence contained in the other? Catches trimmed constructs."""
    s, l = (a, b) if len(a) <= len(b) else (b, a)
    return s in l


seqs = {n: _seq(d) for n, d in proteins.items()}
groups, assigned = [], set()
for n in list(proteins):
    if n in assigned:
        continue
    grp = [n]; assigned.add(n)
    for m in list(proteins):
        if m in assigned:
            continue
        if _overlap(seqs[n], seqs[m]):
            grp.append(m); assigned.add(m)
    groups.append(grp)

for names in groups:
    if len(names) < 2:
        continue
    exact = len({seqs[n] for n in names}) == 1
    kind = "IDENTICAL sequence" if exact else "NESTED sequence (one is a subrange)"
    ranked = sorted(names, key=lambda n: (-_depth(proteins[n]), -len(proteins[n])))
    keep, drop = ranked[0], ranked[1:]
    print(f"\n!! {kind} in {len(names)} inputs: {', '.join(names)}")
    for n in names:
        print(f"     {n}: {len(proteins[n])} positions, MSA depth {_depth_str(proteins[n])}")
    if AUTO_DEDUP_SEQUENCES:
        print(f"   -> keeping '{keep}' (deepest alignment); dropping {', '.join(drop)}."
              f"\n      These are one protein; pooling them would count it more than once.")
        for n in drop:
            proteins.pop(n)
            manifest[:] = [m for m in manifest if m["label"] != n]
    else:
        print("   -> AUTO_DEDUP_SEQUENCES is off; these are being pooled as independent.")

# a disambiguating suffix is only needed while the clash exists. Only suffixes THIS
# notebook added are removed — a parenthetical that is part of a NAME_MAP label
# ("EAA2 / SLC1A2 (QTY)") is meaning, not bookkeeping, and must survive.
for n in list(proteins):
    if n not in SUFFIXED:
        continue
    base = re.sub(r"\s*\([^)]*\)$", "", n)
    if base != n and base not in proteins:
        proteins[base] = proteins.pop(n)
        for m in manifest:
            if m["label"] == n:
                m["label"] = base

print(f"\n{len(proteins)} protein(s) loaded"
      + ("   [BUNDLED DEMO DATA — NOT YOUR FILES]" if USING_BUNDLED else "") + "\n")
print(pd.DataFrame(manifest)[["label", "n_positions", "msa_total", "n_self_missing"]]
      .to_string(index=False))'''))

cells.append(md(r"""**Check the list above before going on.** Two alignments of the *same*
protein (a 150-sequence and a 500-sequence run, say) are not independent, and pooling
them inflates `n`. Drop one with `EXCLUDE`, or rename via `NAME_MAP`, and re-run.
`n_self_missing` should be 0 everywhere."""))

cells.append(md(r"""### Optional: restrict to a residue range

Use this if only part of the sequence was modelled, so the sequence analysis covers the
same span as anything you later join to it."""))

cells.append(code(r'''RANGES = {}          # e.g. {"EAA1 / SLC1A3": (43, 499)}

for name, (lo, hi) in RANGES.items():
    if name in proteins:
        before = len(proteins[name])
        proteins[name] = proteins[name].query("@lo <= seqres_pos <= @hi").copy()
        print(f"{name}: {before} -> {len(proteins[name])} positions")'''))

cells.append(md(r"""### Is each input the native protein?

The QTY code is applied **to** a native membrane protein, so a precedence question can
only be asked of one. Two kinds of input pass every structural check and still cannot
answer it:

* a **QTY-designed** sequence &mdash; the substitutions are already present, so the native
  L/I/V/F the code would act on is missing from the denominator, and the target residue
  inflates its own background (ConSurf lists the query residue in its own alignment
  column, so the self-consistency check above cannot see this).
* a sequence with **no transmembrane segment** &mdash; the code has nowhere to apply, and
  every site measured is a soluble-domain site.

Composition and hydropathy separate the three cases with no external annotation. A native
polytopic transporter runs 33&ndash;37% L+I+V+F against 7&ndash;12% Q+T+Y; its QTY variant
inverts that. Set `FORCE_KIND` to overrule the call."""))

cells.append(code(r'''KD = dict(zip("AVLIPFMWGSTCYNQDEKRH",
              [1.8, 4.2, 3.8, 4.5, -1.6, 2.8, 1.9, -0.9, -0.4, -0.8,
               -0.7, 2.5, -1.3, -3.5, -3.5, -3.5, -3.5, -3.9, -4.5, -3.2]))
SIGNAL_PEPTIDE_LEN = 30        # an N-terminal hydrophobic run inside this is not a TM
TM_WINDOW, TM_CUT, TM_MIN = 19, 1.6, 17
DESIGNED_QTY_OVER_LIVF = 0.8   # native 0.19-0.36, QTY-designed ~1.6
DESIGNED_SELF_BG = 0.20        # native 0.06-0.13, QTY-designed ~0.28
FORCE_KIND = {}                # e.g. {"EAA2 / SLC1A2": "native membrane"}


def tm_segments(seq, skip=SIGNAL_PEPTIDE_LEN, window=None, cut=None, pad=0):
    """Kyte-Doolittle window scan. Crude — but none-vs-several is all it is asked for.
    window/cut/pad are exposed so the boundary definition can be swept, not assumed."""
    window = TM_WINDOW if window is None else window
    cut = TM_CUT if cut is None else cut
    if len(seq) < window:
        return []
    sm = np.convolve(np.array([KD.get(a, 0.0) for a in seq]),
                     np.ones(window) / window, mode="valid")
    runs, i = [], 0
    while i < len(sm):
        if sm[i] >= cut:
            j = i
            while j + 1 < len(sm) and sm[j + 1] >= cut:
                j += 1
            runs.append((i + 1, j + window)); i = j + 1
        else:
            i += 1
    merged = []
    for lo, hi in runs:                       # rejoin a run split by one dip
        if merged and lo - merged[-1][1] <= 5:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    keep = [(a, b) for a, b in merged if b - a + 1 >= TM_MIN and b > skip]
    return [(max(1, a - pad), min(len(seq), b + pad)) for a, b in keep]


def self_background(d):
    """Share of each target's background that the query residue itself supplies."""
    out = {}
    for t in set(QTY_MAP.values()):
        has = d[[t in s for s in d.vset]]
        out[t] = len(has[has.aa == t]) / len(has) if len(has) else 0.0
    return out


def classify(name, d):
    seq = "".join(d.sort_values("seqres_pos").aa)
    livf = sum(seq.count(a) for a in QTY_MAP) / len(seq)
    qty = sum(seq.count(t) for t in set(QTY_MAP.values())) / len(seq)
    sb = max(self_background(d).values())
    segs = tm_segments(seq)
    ratio = qty / livf if livf else np.inf
    kind = ("QTY-designed" if ratio > DESIGNED_QTY_OVER_LIVF or sb > DESIGNED_SELF_BG
            else "soluble (no TM)" if not segs
            else "native membrane")
    return dict(protein=name, n=len(seq), livf=100 * livf, qty=100 * qty,
                qty_over_livf=ratio, tm_segs=len(segs), self_bg=100 * sb,
                kind=FORCE_KIND.get(name, kind),
                forced=name in FORCE_KIND)


klass = {n: classify(n, d) for n, d in proteins.items()}
kt = pd.DataFrame(list(klass.values()))

print("input classification\n")
print(kt[["protein", "n", "livf", "qty", "qty_over_livf", "tm_segs", "self_bg", "kind"]]
      .to_string(index=False, formatters={
          "livf": "{:.1f}%".format, "qty": "{:.1f}%".format,
          "qty_over_livf": "{:.2f}".format, "self_bg": "{:.0f}%".format}))

POOL = [n for n in proteins if klass[n]["kind"] == "native membrane"]
REFERENCE = [n for n in proteins if n not in POOL]

for n in REFERENCE:
    k = klass[n]
    why = (f"the QTY substitutions are already present (Q+T+Y {k['qty']:.0f}% against "
           f"L+I+V+F {k['livf']:.0f}%), so the residues the code would act on are gone "
           f"from the denominator and each target supplies up to {k['self_bg']:.0f}% of "
           f"its own background"
           if k["kind"] == "QTY-designed" else
           "no transmembrane segment outside the signal peptide, so the code has nowhere "
           "to apply and every site scored is a soluble-domain site")
    print(f"\n!! {n} — {k['kind']}.\n   {why}.")
    print("   Measured and reported, but NOT pooled.")

print(f"\npooled   : {', '.join(POOL) if POOL else '(none)'}")
if REFERENCE:
    print(f"reference: {', '.join(REFERENCE)}")
if len(POOL) < 2:
    print("\n!! Fewer than two poolable proteins. Every cross-protein statement below "
          "rests on one sequence.")'''))

cells.append(md(r"""### Transmembrane regions from UniProt

Annotated topology beats a hydropathy scan, and for this family it is not a close call:
UniProt marks the re-entrant hairpins of the SLC1 fold as `Intramembrane`, and a
Kyte&ndash;Doolittle window cannot see a hairpin at all. Those are membrane-embedded
residues the QTY code would act on.

The fetch is guarded. Coordinates are only accepted after the UniProt sequence is
reconciled **residue for residue** against the sequence ConSurf actually scored, and the
offset between the two numbering schemes is derived from that match rather than assumed.
A mismatch, an ambiguous match or a network failure leaves the protein on the hydropathy
fallback and says so. Annotation pulled from the wrong accession is worse than a guess,
because it looks authoritative.

Give an accession or a gene symbol; anything not matching the accession pattern is
looked up as a reviewed human gene."""))

cells.append(code(r'''import urllib.request, urllib.parse

UNIPROT_IDS = {                      # accession, or gene symbol to search
    "EAA1 / SLC1A3":  "SLC1A3",
    "VPAT / SLC18B1": "SLC18B1",
}
FETCH_UNIPROT_TM = True
INCLUDE_INTRAMEM = True              # re-entrant hairpins count as membrane-embedded
ORGANISM_ID = 9606
UNIPROT = "https://rest.uniprot.org/uniprotkb"
ACC_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")
TM_REGIONS = {}


def _get_json(url, timeout=30):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def uniprot_lookup(token):
    """Accession or gene symbol -> (accession, entry)."""
    if ACC_RE.match(token):
        return token, _get_json(f"{UNIPROT}/{token}.json")
    q = f"gene_exact:{token} AND organism_id:{ORGANISM_ID} AND reviewed:true"
    hits = _get_json(f"{UNIPROT}/search?query={urllib.parse.quote(q)}"
                     f"&format=json&size=5").get("results", [])
    if not hits:
        raise LookupError(f"no reviewed entry for gene {token}")
    return hits[0]["primaryAccession"], hits[0]


def membrane_features(entry):
    out = []
    for f in entry.get("features", []):
        if f.get("type") not in {"Transmembrane", "Intramembrane"}:
            continue
        loc = f.get("location", {})
        s = loc.get("start", {}).get("value"); e = loc.get("end", {}).get("value")
        if s is None or e is None:
            continue                               # fuzzy position: unusable
        out.append((int(s), int(e), f["type"], f.get("description", "")))
    return sorted(out)


def reconcile(consurf_seq, uniprot_seq, first_pos):
    """UniProt position p maps to ConSurf position p + offset. None = do not use."""
    if not consurf_seq or not uniprot_seq:
        return None, "empty sequence"
    if consurf_seq == uniprot_seq:
        return first_pos - 1, "identical sequences"
    if consurf_seq in uniprot_seq:
        if uniprot_seq.count(consurf_seq) > 1:
            return None, "query matches UniProt in more than one place"
        i = uniprot_seq.index(consurf_seq)
        return first_pos - 1 - i, (f"query is UniProt {i + 1}-{i + len(consurf_seq)}; "
                                   f"{len(uniprot_seq) - len(consurf_seq)} aa not scored")
    if uniprot_seq in consurf_seq:
        return first_pos - 1 + consurf_seq.index(uniprot_seq), "UniProt seq is a subrange"
    if len(consurf_seq) == len(uniprot_seq):
        ident = sum(a == b for a, b in zip(consurf_seq, uniprot_seq)) / len(consurf_seq)
        if ident >= 0.98:
            return first_pos - 1, f"same length, {100 * ident:.1f}% identical"
        return None, f"same length but only {100 * ident:.1f}% identical"
    return None, (f"sequences differ (ConSurf {len(consurf_seq)}, "
                  f"UniProt {len(uniprot_seq)}) and neither contains the other")


if FETCH_UNIPROT_TM and UNIPROT_IDS:
    for name, token in UNIPROT_IDS.items():
        if name not in proteins:
            print(f"  {name}: not loaded, skipped"); continue
        d = proteins[name]
        seq = "".join(d.sort_values("seqres_pos").aa)
        first = int(d.seqres_pos.min())
        try:
            acc, entry = uniprot_lookup(token)
            useq = entry.get("sequence", {}).get("value", "")
            off, note = reconcile(seq, useq, first)
            feats = membrane_features(entry)
            if off is None:
                print(f"  {name}: {acc} REJECTED — {note}")
                print(f"      falling back to the hydropathy scan for this protein")
                continue
            lo, hi = first, first + len(seq) - 1
            regs, kinds = [], []
            for s, e, typ, _ in feats:
                if typ == "Intramembrane" and not INCLUDE_INTRAMEM:
                    continue
                a, b = max(s + off, lo), min(e + off, hi)
                if a <= b:
                    regs.append((a, b)); kinds.append(typ)
            if not regs:
                print(f"  {name}: {acc} has no usable membrane features"); continue
            TM_REGIONS[name] = sorted(regs)
            n_tm = kinds.count("Transmembrane"); n_im = kinds.count("Intramembrane")
            print(f"  {name}: {acc} ok — {note}")
            print(f"      {n_tm} transmembrane + {n_im} intramembrane -> "
                  f"{len(regs)} regions, offset {off:+d}")
        except Exception as exc:
            print(f"  {name}: lookup failed ({type(exc).__name__}: {exc})")
            print(f"      falling back to the hydropathy scan for this protein")
else:
    print("FETCH_UNIPROT_TM off — using the hydropathy scan for every protein.")

print(f"\nTM_REGIONS declared for {len(TM_REGIONS)}/{len(proteins)} input(s)"
      + (f": {', '.join(TM_REGIONS)}" if TM_REGIONS else ""))'''))

cells.append(md(r"""### Where does the code actually act?

The QTY code rewrites **transmembrane** segments. A ratio computed over a whole
polytopic transporter therefore mixes the sites the code would rewrite with soluble-loop
sites it would never touch, and the loops dilute whatever is happening in the bilayer.

`TM_REGIONS` holds real topology as SEQRES ranges per protein, normally filled by the
UniProt cell above; set it by hand to use boundaries read off a structure instead. The
fallback, for any protein it does not cover, is the same Kyte&ndash;Doolittle scan used for classification above, which is
good enough to tell a transporter from a secreted ectodomain and *not* good enough to
place a helix boundary: it merges adjacent helices, misses re-entrant hairpins, and will
put a residue a few positions either side of where a structure puts it.

So the primary scope is chosen for you: **TM-restricted only when every pooled protein
has declared regions**, whole-protein otherwise. Both are always computed and compared,
because the difference between them is itself a result."""))

cells.append(code(r'''# TM_REGIONS comes from the UniProt cell above; edit it there, or set it by hand
# here with e.g. TM_REGIONS["EAA1 / SLC1A3"] = [(48, 72), (92, 112), ...]
PRIMARY_SCOPE = "auto"               # "auto" | "whole protein" | "TM only"

tm_source = {}
for name, d in proteins.items():
    if name in TM_REGIONS:
        segs, src = list(TM_REGIONS[name]), "declared"
    else:
        seq = "".join(d.sort_values("seqres_pos").aa)
        off = int(d.seqres_pos.min()) - 1
        segs = [(a + off, b + off) for a, b in tm_segments(seq)]
        src = "predicted"
    tm_source[name] = src
    proteins[name]["in_tm"] = d.seqres_pos.apply(
        lambda p: any(lo <= p <= hi for lo, hi in segs))
    d.attrs["tm_regions"] = segs

print(f"  {'protein':<22}{'source':>11}{'segments':>10}{'% of residues':>15}")
for name, d in proteins.items():
    print(f"  {name:<22}{tm_source[name]:>11}{len(d.attrs['tm_regions']):>10}"
          f"{100 * d.in_tm.mean():>14.0f}%")

all_declared = bool(POOL) and all(tm_source[n] == "declared" for n in POOL)
SCOPE = ("TM only" if all_declared else "whole protein") \
        if PRIMARY_SCOPE == "auto" else PRIMARY_SCOPE
print(f"\nprimary scope: {SCOPE}" +
      ("" if all_declared else
       "   (TM_REGIONS not declared for every pooled protein — the TM-restricted\n"
       "    numbers below rest on a hydropathy guess and are a sensitivity check,\n"
       "    not the headline)"))


def scoped(name):
    d = proteins[name]
    return d[d.in_tm] if SCOPE == "TM only" else d'''))

cells.append(md(r"""## 5 · Precedence and the background control

Per native residue the code touches:

* **observed** &mdash; fraction of those sites whose alignment column contains the QTY target
* **background** &mdash; fraction of *all* positions in that protein containing that same residue

Their ratio is the quantity of interest. Below 1 the target is *depleted* exactly where
the code applies it; above 1 evolution already favours the swap there. The background
term is what stops a residue that is simply rare from looking meaningful."""))

cells.append(code(r'''PREC_COLS = ["native", "target", "n_sites", "n_prec", "obs", "bg", "ratio", "p"]


def precedence(df, qty_map=QTY_MAP):
    """Always returns the full column set, so an empty result stays a usable frame.

    Empty is a real case, not an error: restrict a soluble protein to its
    transmembrane segments and there is nothing left to score."""
    if len(df) == 0:
        return pd.DataFrame(columns=PREC_COLS)
    bg = {t: float(np.mean([t in s for s in df.vset])) for t in set(qty_map.values())}
    rows = []
    for native, target in qty_map.items():
        sites = df[df.aa == native]
        if len(sites) == 0:
            continue
        k, n, p0 = int(sum(target in s for s in sites.vset)), len(sites), bg[target]
        rows.append(dict(native=native, target=target, n_sites=n, n_prec=k,
                         obs=100 * k / n, bg=100 * p0,
                         ratio=(k / n) / p0 if p0 else np.nan,
                         p=stats.binomtest(k, n, p0).pvalue if p0 else np.nan))
    return pd.DataFrame(rows, columns=PREC_COLS)


per_protein = {name: precedence(scoped(name)) for name in proteins}
EMPTY = [n for n, t in per_protein.items() if not len(t)]

print(f"scope: {SCOPE}")
for name, t in per_protein.items():
    print(f"\n{name}")
    if not len(t):
        n_in = int(proteins[name].in_tm.sum())
        print(f"  no scorable sites under '{SCOPE}' "
              f"({n_in} residues inside the membrane regions) — dropped from here on")
        continue
    print(t.assign(swap=t.native + "->" + t.target)
           [["swap", "n_sites", "n_prec", "obs", "bg", "ratio", "p"]]
           .to_string(index=False,
                      formatters={"obs": "{:.0f}%".format, "bg": "{:.0f}%".format,
                                  "ratio": "{:.2f}".format, "p": "{:.1e}".format}))

if EMPTY:
    POOL = [n for n in POOL if n not in EMPTY]
    REFERENCE = [n for n in REFERENCE if n not in EMPTY]
    print(f"\n{len(EMPTY)} input(s) have nothing to score under '{SCOPE}' "
          f"and are excluded: {', '.join(EMPTY)}")'''))

cells.append(md(r"""## 6 · Pool across proteins

Only the native membrane proteins in `POOL` are pooled. Each carries its own
expectation, so a shallow alignment cannot drag the pooled estimate. The test is a
normal deviate on a sum of independent binomials (a Poisson-binomial, summarised by its
mean and variance)."""))

cells.append(code(r'''POOL_COLS = ["native", "target", "n_sites", "n_prec", "expected", "ratio",
             "ci_lo", "ci_hi", "det_lo", "det_hi", "z", "p",
             "n_proteins", "all_below_1", "all_above_1"]


def pool(per_protein):
    frames = [t.assign(protein=n) for n, t in per_protein.items() if len(t)]
    if not frames:
        return pd.DataFrame(columns=POOL_COLS)
    long = pd.concat(frames, ignore_index=True)
    rows = []
    for (native, target), g in long.groupby(["native", "target"], sort=False):
        n, k = int(g.n_sites.sum()), int(g.n_prec.sum())
        exp = float((g.n_sites * g.bg / 100).sum())
        var = float((g.n_sites * (g.bg / 100) * (1 - g.bg / 100)).sum())
        sd = np.sqrt(var) if var > 0 else np.nan
        z = (k - exp) / sd if var > 0 else np.nan
        # CI on the ratio, and the effect this design could have resolved at all:
        # a non-significant result means nothing without the second number.
        ci_lo = max(0.0, k - 1.96 * sd) / exp if exp else np.nan
        ci_hi = (k + 1.96 * sd) / exp if exp else np.nan
        det_lo = max(0.0, exp - 1.96 * sd) / exp if exp else np.nan
        det_hi = (exp + 1.96 * sd) / exp if exp else np.nan
        rows.append(dict(native=native, target=target, n_sites=n, n_prec=k,
                         expected=exp, ratio=k / exp if exp else np.nan,
                         ci_lo=ci_lo, ci_hi=ci_hi, det_lo=det_lo, det_hi=det_hi, z=z,
                         p=2 * stats.norm.sf(abs(z)), n_proteins=len(g),
                         all_below_1=bool((g.ratio < 1).all()),
                         all_above_1=bool((g.ratio > 1).all())))
    order = {k: i for i, k in enumerate(QTY_MAP)}
    return (pd.DataFrame(rows, columns=POOL_COLS)
              .sort_values("native", key=lambda s: s.map(order))
              .reset_index(drop=True))


pooled = pool({n: per_protein[n] for n in POOL})
pooled["swap"] = pooled.native + "->" + pooled.target
# "no effect" and "the proteins disagree" are different findings and get different
# words. Lumping them under one label makes a clean null read as a failure.
pooled["verdict"] = np.select(
    [pooled.p >= .05,
     pooled.all_below_1, pooled.all_above_1],
    ["no effect", "depleted in all", "enriched in all"],
    default="significant, but proteins differ in sign")

print(f"pooled over {len(POOL)} native membrane protein(s): {', '.join(POOL)}\n")
_show = pooled.assign(**{"95% CI": [f"[{a:.2f}, {b:.2f}]" for a, b
                                    in zip(pooled.ci_lo, pooled.ci_hi)]})
print(_show[["swap", "n_sites", "n_prec", "expected", "ratio", "95% CI", "p", "verdict"]]
      .to_string(index=False, formatters={"expected": "{:.1f}".format,
                                          "ratio": "{:.2f}".format,
                                          "p": "{:.1e}".format}))

_null = pooled[pooled.p >= .05]
if len(_null):
    print("\n  A null is only as strong as the effect the design could resolve.")
    print("  At these site counts the test reaches p < 0.05 only outside these bounds:")
    for _, r in _null.iterrows():
        print(f"    {r.swap:<8} n={int(r.n_sites):<4} would need ratio <= "
              f"{r.det_lo:.2f} or >= {r.det_hi:.2f}   "
              f"(observed {r.ratio:.2f}, CI [{r.ci_lo:.2f}, {r.ci_hi:.2f}])")
    print("  Where the CI still spans 1 by a wide margin, read 'underpowered',")
    print("  not 'no effect' — especially if every protein points the same way.")'''))

cells.append(md(r"""### Whole protein against transmembrane only

Both scopes, pooled, side by side. `shift` is the TM ratio divided by the whole-protein
one: above 1 the effect is *stronger* inside the membrane than across the protein,
below 1 it is weaker. A swap whose effect is real membrane chemistry should sharpen when
the soluble loops are dropped; one that is an artefact of overall composition should
not move."""))

cells.append(code(r'''def pool_scope(names, tm_only):
    per = {}
    for n in names:
        d = proteins[n]
        per[n] = precedence(d[d.in_tm] if tm_only else d)
    return pool(per)


if POOL:
    w, t = pool_scope(POOL, False), pool_scope(POOL, True)
    w["swap"] = w.native + "->" + w.target
    t["swap"] = t.native + "->" + t.target
    cmp_ = (w[["swap", "n_sites", "n_prec", "ratio", "p"]]
            .merge(t[["swap", "n_sites", "n_prec", "ratio", "p"]],
                   on="swap", suffixes=("_whole", "_tm")))
    cmp_["tm_over_whole"] = cmp_.ratio_tm / cmp_.ratio_whole
    print(f"pooled over {', '.join(POOL)}\n")
    print(f"  {'swap':<8}{'n':>5}{'k':>5}{'whole':>8}{'p':>10}   "
          f"{'n':>5}{'k':>5}{'TM':>8}{'p':>10}{'shift':>8}")
    for _, r in cmp_.iterrows():
        print(f"  {r.swap:<8}{int(r.n_sites_whole):>5}{int(r.n_prec_whole):>5}"
              f"{r.ratio_whole:>8.2f}{r.p_whole:>10.1e}   "
              f"{int(r.n_sites_tm):>5}{int(r.n_prec_tm):>5}{r.ratio_tm:>8.2f}"
              f"{r.p_tm:>10.1e}{r.tm_over_whole:>8.2f}")
    print("\n  k = sites whose alignment column already contains the QTY target.")
    if not all_declared:
        print("  TM columns rest on predicted boundaries — treat as a sensitivity check.")

    # Why the two scopes can disagree: the DENOMINATOR moves between them.
    # A target residue that is common in soluble domains and rare in the bilayer
    # makes the whole-protein ratio look depleted for a reason that has nothing
    # to do with the sites being substituted.
    print("\n\nbackground availability of each target, whole protein vs membrane")
    print(f"  {'protein':<22}{'target':>7}{'whole':>8}{'TM':>7}{'TM/whole':>11}")
    flags = []
    for n in POOL:
        d = proteins[n]
        for t in sorted(set(QTY_MAP.values())):
            bw = float(np.mean([t in s for s in d.vset]))
            bt = float(np.mean([t in s for s in d[d.in_tm].vset])) if d.in_tm.any() else np.nan
            r = bt / bw if bw else np.nan
            mark = ""
            if np.isfinite(r) and (r < 0.7 or r > 1.4):
                mark = "   <-- scope-sensitive"; flags.append((n, t, r))
            print(f"  {n:<22}{t:>7}{100 * bw:>7.0f}%{100 * bt:>6.0f}%{r:>11.2f}{mark}")
    if flags:
        tg = sorted({t for _, t, _ in flags})
        verb = "differs" if len(tg) == 1 else "differ"
        print(f"\n  {', '.join(tg)} {verb} markedly in availability between the two")
        print("  scopes, so for any swap targeting them the whole-protein ratio compares")
        print("  membrane-resident native residues against a background drawn largely")
        print("  from soluble domains. Prefer the TM-restricted estimate for those.")
    else:
        print("\n  No target's availability differs much between scopes, so the two")
        print("  ratios are comparable and any difference is about the sites, not")
        print("  the denominator.")
else:
    cmp_ = None
    print("Nothing pooled.")'''))

cells.append(md(r"""### Does the boundary definition change the answer?

The TM numbers above rest on one arbitrary boundary rule. Rather than defend it, sweep
it: every combination of smoothing window, hydropathy threshold and edge padding, from
a definition that keeps only helix cores to one that spills well into the loops. If a
swap keeps the same sign across the whole sweep, the conclusion does not depend on the
boundary and can be stated without exact topology. If it changes sign, nothing can be
said until real `TM_REGIONS` are supplied.

Magnitude is a separate question from sign, and the sweep reports both: a swap can be
robustly depleted and still have a ratio that no single number represents."""))

cells.append(code(r'''TM_SENS_GRID = dict(window=[15, 19, 21], cut=[1.2, 1.4, 1.6, 1.8, 2.0],
                    pad=[-5, 0, 5])
SENS_ALPHA = 0.05
import itertools


def _segs_under(name, window, cut, pad):
    """Declared regions ignore window/cut — only the padding applies to them."""
    d = proteins[name]
    if name in TM_REGIONS:
        lo_, hi_ = int(d.seqres_pos.min()), int(d.seqres_pos.max())
        return [(max(lo_, a - pad), min(hi_, b + pad)) for a, b in TM_REGIONS[name]]
    seq = "".join(d.sort_values("seqres_pos").aa)
    off = int(d.seqres_pos.min()) - 1
    return [(a + off, b + off)
            for a, b in tm_segments(seq, window=window, cut=cut, pad=pad)]


def tm_sensitivity(names, grid=TM_SENS_GRID, min_sites=10):
    """Distinct boundary definitions, not distinct grid points.

    Declared TM_REGIONS ignore window and threshold, so the full grid would count
    the same three paddings fifteen times each and overstate the evidence."""
    rows, seen = [], set()
    for window, cut, pad in itertools.product(grid["window"], grid["cut"], grid["pad"]):
        segs_by = {n: tuple(_segs_under(n, window, cut, pad)) for n in names}
        key = tuple(sorted(segs_by.items()))
        if key in seen:
            continue
        seen.add(key)
        per, cov, ok = {}, [], True
        for n in names:
            d = proteins[n]
            segs = segs_by[n]
            mask = d.seqres_pos.apply(lambda p: any(lo <= p <= hi for lo, hi in segs))
            if mask.sum() < 20:
                ok = False; break
            per[n] = precedence(d[mask]); cov.append(100 * mask.mean())
        if not ok:
            continue
        p = pool(per)
        for _, r in p.iterrows():
            if r.n_sites < min_sites or not np.isfinite(r.ratio):
                continue
            rows.append(dict(window=window, cut=cut, pad=pad,
                             coverage=float(np.mean(cov)),
                             swap=f"{r.native}->{r.target}", ratio=float(r.ratio),
                             n_sites=int(r.n_sites), p=float(r.p)))
    return pd.DataFrame(rows)


if POOL:
    sens = tm_sensitivity(POOL)
    n_defs = sens[["window", "cut", "pad"]].drop_duplicates().shape[0]
    print(f"{n_defs} boundary definitions, membrane coverage "
          f"{sens.coverage.min():.0f}%–{sens.coverage.max():.0f}% of residues\n")
    print(f"  {'swap':<9}{'min':>8}{'median':>9}{'max':>8}{'whole':>9}"
          f"{'defs':>6}{'sig':>5}   verdict")
    rows_sum = []
    for native, target in QTY_MAP.items():
        sw = f"{native}->{target}"
        g = sens[sens.swap == sw]
        if not len(g):
            continue
        w = pooled[pooled.swap == sw]
        wr = float(w.ratio.iloc[0]) if len(w) else np.nan
        # A ratio straddling 1 means one of two different things, and only one of
        # them is a problem: a null straddles 1 because there is no effect, while a
        # boundary-dependent result is significant in BOTH directions. Counting
        # significant definitions separates them; counting signs does not.
        sig_lo = int(((g.ratio < 1) & (g.p < SENS_ALPHA)).sum())
        sig_hi = int(((g.ratio > 1) & (g.p < SENS_ALPHA)).sum())
        if sig_lo and sig_hi:
            verdict = "CONTRADICTORY — significant both ways; needs real TM_REGIONS"
        elif not sig_lo and not sig_hi:
            verdict = f"null under every definition (min p = {g.p.min():.2f})"
        else:
            which = "depleted" if sig_lo else "enriched"
            verdict = f"{which}; significant in {max(sig_lo, sig_hi)}/{len(g)}"
        rows_sum.append(dict(swap=sw, lo=g.ratio.min(), med=g.ratio.median(),
                             hi=g.ratio.max(), whole=wr, n_defs=len(g),
                             n_sig_depleted=sig_lo, n_sig_enriched=sig_hi,
                             min_p=float(g.p.min()),
                             contradictory=bool(sig_lo and sig_hi), verdict=verdict))
        print(f"  {sw:<9}{g.ratio.min():>8.2f}{g.ratio.median():>9.2f}"
              f"{g.ratio.max():>8.2f}{wr:>9.2f}{len(g):>6}"
              f"{max(sig_lo, sig_hi):>5}   {verdict}")
    sens_summary = pd.DataFrame(rows_sum)
    print(f"\n  sig = definitions reaching p < {SENS_ALPHA} in the stated direction.")
    print("  A swap significant in BOTH directions is the only one that truly needs")
    print("  real topology. One that is never significant is a null, not a failure.")
    print("  The min-max spread is the price of not having exact boundaries: quote it,")
    print("  or supply TM_REGIONS and quote one number.")
else:
    sens = sens_summary = None'''))

cells.append(md(r"""## 7 · Do the pooled proteins agree?

A fixed-effect pool assumes every protein estimates the same quantity. That is an
assumption, not a given: a swap can be enriched in one fold and depleted in another.

With **three or more** proteins Cochran's *Q* and *I*&sup2; test it, and a
DerSimonian&ndash;Laird random-effects estimate gives an interval that admits
between-protein variation. Read *I*&sup2; as the share of variation that is real
disagreement rather than sampling noise: under ~40% pool freely, 40&ndash;75% report the
random-effects estimate, above 75% do not pool at all.

With **two**, *Q* has one degree of freedom and *I*&sup2; is not interpretable &mdash; it
is as likely to read 0% as 90% on the same underlying agreement. The notebook refuses to
compute it and prints the two estimates side by side instead, with the ratio between them
as the plain measure of agreement."""))

cells.append(code(r'''def _logratio_table(per_protein):
    """Per-protein log ratio with a delta-method SE, for meta-analysis."""
    rows = []
    for name, t in per_protein.items():
        for _, r in t.iterrows():
            if r.n_prec == 0:
                continue                      # log undefined; reported separately
            p = r.n_prec / r.n_sites
            rows.append(dict(protein=name, swap=f"{r.native}->{r.target}",
                             native=r.native, ratio=r.ratio, n=int(r.n_sites),
                             logr=np.log(r.ratio),
                             se=np.sqrt((1 - p) / (r.n_sites * p))))
    return pd.DataFrame(rows)


def meta(g):
    """Fixed- and random-effects pool of log ratios, with Q and I-squared."""
    w = 1 / g.se ** 2
    fe = float((w * g.logr).sum() / w.sum())
    Q = float((w * (g.logr - fe) ** 2).sum())
    df = len(g) - 1
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    denom = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - df) / denom) if denom > 0 else 0.0
    wr = 1 / (g.se ** 2 + tau2)
    re = float((wr * g.logr).sum() / wr.sum())
    se_re = float(np.sqrt(1 / wr.sum()))
    return dict(fe=np.exp(fe), re=np.exp(re),
                lo=np.exp(re - 1.96 * se_re), hi=np.exp(re + 1.96 * se_re),
                p=2 * stats.norm.sf(abs(re / se_re)) if se_re > 0 else np.nan,
                Q=Q, df=df, pQ=1 - stats.chi2.cdf(Q, df) if df > 0 else np.nan,
                I2=I2, tau2=tau2, k=len(g))


def heterogeneity(per_protein, subset=None, label="all proteins"):
    lr = _logratio_table(per_protein)
    if subset:
        lr = lr[lr.protein.isin(subset)]
    out = []
    for native, target in QTY_MAP.items():
        g = lr[lr.swap == f"{native}->{target}"]
        if len(g) < 2:
            continue
        m = meta(g); m["swap"] = f"{native}->{target}"
        m["verdict"] = ("consistent — pool" if m["I2"] < 40 else
                        "moderate — use RE" if m["I2"] < 75 else "HIGH — do not pool")
        out.append(m)
    out = pd.DataFrame(out)[["swap", "k", "fe", "re", "lo", "hi", "p",
                             "Q", "pQ", "I2", "verdict"]]
    print(f"\n{label}  (k = proteins per swap)")
    print(out.to_string(index=False, formatters={
        "fe": "{:.2f}".format, "re": "{:.2f}".format, "lo": "{:.2f}".format,
        "hi": "{:.2f}".format, "p": "{:.1e}".format, "Q": "{:.1f}".format,
        "pQ": "{:.3f}".format, "I2": "{:.0f}%".format}))
    return out


MIN_K_FOR_META = 3


def agreement(per_protein, subset, label):
    """Side-by-side ratios. What to report when k is too small for Q and I-squared.

    Built from per_protein rather than the log-ratio table, which drops k = 0.
    A dropped protein leaves a NaN, and max/min over one remaining value returns
    1.00 — a missing protein would otherwise be printed as perfect agreement."""
    rows = []
    for n in subset:
        t = per_protein.get(n)
        if t is None or not len(t):
            continue
        for _, r in t.iterrows():
            rows.append(dict(protein=n, swap=f"{r.native}->{r.target}",
                             ratio=float(r.ratio)))
    if not rows:
        print(f"\n{label} — nothing to compare"); return None
    wide = (pd.DataFrame(rows).pivot_table(index="swap", columns="protein",
                                           values="ratio")
            .reindex([f"{a}->{t}" for a, t in QTY_MAP.items()]).dropna(how="all"))
    out = wide.copy()
    n_have = wide.notna().sum(axis=1)
    mn, mx = wide.min(axis=1), wide.max(axis=1)
    out["n"] = n_have
    out["spread"] = np.where(n_have < 2, np.nan,
                             np.where(mn > 0, mx / mn, np.inf))
    out["same side"] = np.where(n_have < 2, False,
                                ((wide > 1).all(axis=1) | (wide < 1).all(axis=1)))
    print(f"\n{label} — per-protein ratios side by side")
    print(out.to_string(float_format=lambda v: f"{v:.2f}"))
    print(f"\n  n = proteins contributing that swap (of {len(subset)}).")
    print("  spread = larger ratio / smaller ratio; 1.00 is exact agreement, and")
    print("  blank where only one protein has the swap. inf means one ratio is 0.")
    print("  'same side' = every contributing protein falls on the same side of 1.")
    return out


if len(POOL) >= MIN_K_FOR_META:
    het_pool = heterogeneity(per_protein, subset=POOL,
                             label=f"POOL ({len(POOL)} proteins)")
    agree = None
else:
    het_pool = None
    print(f"{len(POOL)} protein(s) in the pool. Cochran's Q would carry "
          f"df = {max(0, len(POOL) - 1)}; I-squared is not interpretable at that size, "
          f"so no meta-analysis is computed.")
    agree = agreement(per_protein, POOL, "POOL") if len(POOL) >= 2 else None'''))

cells.append(md(r"""### The reference inputs

The inputs excluded from the pool still carry information, and the most useful kind is a
negative. Each is compared against the pooled estimate per swap. A reference protein that
tracks the pool says the effect is not peculiar to the pooled folds; one that departs
sharply localises where the effect comes from.

**TLR5S is the strongest case of the second kind.** It is a distinct teleost gene whose
product is the leucine-rich-repeat ectodomain alone &mdash; no transmembrane segment, no
TIR domain, secreted. It is a membrane receptor that evolution has already solubilised,
and it did so by *removing the anchor*, not by polarising it. Its leucines are the
structural L of the LRR consensus, buried in the solenoid rather than in a bilayer, and
they carry the most extreme Q depletion in the whole set. Read that way it is not an
outlier to be excused but a control on what the depletion means: wherever leucine is doing
structural work, glutamine is not an alternative evolution has sampled &mdash; bilayer or
no bilayer."""))

cells.append(code(r'''if REFERENCE and len(POOL):
    # Read ratios straight from per_protein, not from the log-ratio table: that
    # table drops k = 0 because log(0) is undefined, which would print the single
    # most extreme observation — the target never sampled — as missing data.
    _w = max(12, *(len(n) + 3 for n in REFERENCE))
    print(f"{'swap':<9}{'pool':>7}" + "".join(f"{n:>{_w}}" for n in REFERENCE))
    for native, target in QTY_MAP.items():
        sw = f"{native}->{target}"
        pr = pooled[pooled.swap == sw]
        if not len(pr):
            continue
        line = f"{sw:<9}{pr.ratio.iloc[0]:>7.2f}"
        for n in REFERENCE:
            t = per_protein.get(n)
            g = t[(t.native == native) & (t.target == target)] if t is not None \
                and len(t) else None
            if g is None or not len(g):
                line += f"{chr(8212):>{_w}}"           # swap absent in that input
            else:
                v, k = float(g.ratio.iloc[0]), int(g.n_prec.iloc[0])
                line += (f"{'0.00*':>{_w}}" if k == 0 else f"{v:>{_w}.2f}")
        print(line)
    print("\n  * 0.00 means the target is never sampled at any of those sites —")
    print("    an observed value, not a missing one. An em dash means the swap")
    print("    has no sites at all in that input under the current scope.")
    print("\n  A reference column is NOT evidence for the pooled number — it is a")
    print("  different kind of sequence. It shows where the effect does and does not")
    print("  travel.")
else:
    print("No reference inputs, or nothing pooled.")'''))

cells.append(md(r"""## 8 · Figure — enrichment and depletion per swap"""))

cells.append(code(r'''REF_MARKERS = ["o", "s", "^", "v", "P"]
DY_POOL, DY_REF = -0.26, 0.26        # keep the diamond's own line clear


def _ratio_at(t, row):
    m = t[(t.native == row.native) & (t.target == row.target)]
    return float(m.ratio.iloc[0]) if len(m) else None


def _short(n):
    """'EAA2 / SLC1A2 (QTY)' -> 'EAA2 (QTY)'."""
    head = n.split(" / ")[0]
    m = re.search(r"\(([^)]*)\)$", n)
    if m and m.group(0) not in head:
        head = f"{head} ({m.group(1)})"
    return head


def plot_ratios(per_protein, pooled, het=None, reference=None, alpha=0.05,
                fname="fig_precedence.png"):
    reference = reference or {}
    n_sw = len(pooled)
    n_leg = 3 + (1 if het is not None and len(het) else 0) + len(reference)
    leg_rows = int(np.ceil(n_leg / (2 if n_leg > 3 else max(1, n_leg))))
    fig, ax = plt.subplots(figsize=(7.1, 0.92 * n_sw + 1.05 + 0.26 * leg_rows))
    ax.axvline(1, color="#8a8984", ls="--", lw=1.2, zorder=1)
    zeros = []                 # ratios of exactly 0: no position on a log axis

    for i, row in pooled.iterrows():
        y = n_sw - 1 - i
        # pooled proteins sit BELOW the diamond so close agreement stays visible.
        # The bar spans them: its WIDTH is the disagreement, which is what Q and
        # I-squared would have measured had there been enough proteins to try.
        vs = [v for t in per_protein.values() if (v := _ratio_at(t, row)) is not None]
        pos = [v for v in vs if v > 0]
        zeros.extend((y + DY_POOL, "pool", None) for v in vs if v <= 0)
        if len(pos) > 1:
            ax.plot([min(pos), max(pos)], [y + DY_POOL] * 2, "-", color=C_NEU,
                    lw=1.4, alpha=.65, solid_capstyle="round", zorder=2.5)
        for v in pos:
            ax.plot(v, y + DY_POOL, "o", ms=5, color=C_NEU,
                    mec="white", mew=1.0, zorder=3)
        # excluded inputs sit ABOVE it, open, one marker shape each
        for j, (name, t) in enumerate(reference.items()):
            v = _ratio_at(t, row)
            if v is None:
                continue
            if v > 0:
                ax.plot(v, y + DY_REF, REF_MARKERS[j % len(REF_MARKERS)], ms=5.5,
                        mfc="none", mec="#8a8984", mew=1.3, zorder=3)
            else:
                zeros.append((y + DY_REF, "ref", REF_MARKERS[j % len(REF_MARKERS)]))
        col = C_NEU if row.p >= alpha else (C_DEP if row.ratio < 1 else C_ENR)
        if het is not None:
            h = het[het.swap == row.swap]
            if len(h):
                h = h.iloc[0]
                ax.plot([h.lo, h.hi], [y, y], "-", color=col, lw=2.0,
                        solid_capstyle="round", alpha=.55, zorder=1.5)
        if row.ratio > 0:
            ax.plot(row.ratio, y, "D", ms=11, color=col, zorder=2)
        else:
            zeros.append((y, "pooled", col))
        # the label shares the diamond's line, so only the diamond and its CI
        # can collide with it — the dots are on their own sub-rows
        ci = []
        if het is not None and len(het[het.swap == row.swap]):
            hh = het[het.swap == row.swap].iloc[0]; ci = [hh.lo, hh.hi]
        span = [row.ratio] + list(ci)
        if row.ratio >= 1:
            ax.text(max(span) * 1.08, y, f"{row.ratio:.2f}×", va="center",
                    ha="left", fontsize=9.5, fontweight="bold")
        else:
            ax.text(min(span) * 0.93, y, f"{row.ratio:.2f}×", va="center",
                    ha="right", fontsize=9.5, fontweight="bold")

    ax.set_yticks(range(n_sw))
    ax.set_yticklabels([f"{r.native} → {r.target}\n$n$={r.n_sites}, $p$={r.p:.0e}"
                        for _, r in pooled[::-1].iterrows()], fontsize=9)
    ax.set_ylim(-0.6, n_sw - 0.4)
    ax.set_xscale("log")
    # A ratio of exactly 0 (the target never appears at any of those sites) is a
    # real and extreme observation, but it has no position on a log axis. Keep it
    # out of the limits, then draw it clamped at the left edge as a bounded value.
    vals = [v for v in pooled.ratio if np.isfinite(v) and v > 0]
    for t in list(per_protein.values()) + list(reference.values()):
        vals += [v for v in t.ratio if np.isfinite(v) and v > 0]
    if het is not None and len(het):
        vals += [v for v in list(het.lo) + list(het.hi) if np.isfinite(v) and v > 0]
    if not vals:
        vals = [1.0]
    lo = min(0.45, min(vals) * 0.62)      # headroom for the value labels
    hi = max(2.2, max(vals) * 1.55)
    ax.set_xlim(lo, hi)

    # now that the axis exists, pin the zeros to its left edge as bounded values
    for yy, kind, extra in zeros:
        if kind == "pooled":
            ax.plot(lo * 1.03, yy, "<", ms=9, color=extra or C_DEP, zorder=4,
                    clip_on=False)
        elif kind == "ref":
            ax.plot(lo * 1.03, yy, "<", ms=5.5, mfc="none", mec="#8a8984",
                    mew=1.3, zorder=4, clip_on=False)
        else:
            ax.plot(lo * 1.03, yy, "<", ms=5, color=C_NEU, mec="white", mew=1.0,
                    zorder=4, clip_on=False)

    ticks = [t for t in (0.25, 0.5, 1, 2, 3, 4) if lo <= t <= hi]
    ax.xaxis.set_major_locator(mpl.ticker.FixedLocator(ticks))
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.xaxis.set_major_formatter(
        mpl.ticker.FixedFormatter([f"{t:g}×" for t in ticks]))
    ax.set_xlabel("observed rate ÷ within-protein background rate")

    ax.text(0.95, 1.02, "depleted", transform=ax.get_xaxis_transform(), ha="right",
            va="bottom", fontsize=8.5, style="italic", color="#8a8984")
    ax.text(1.06, 1.02, "already favoured", transform=ax.get_xaxis_transform(),
            ha="left", va="bottom", fontsize=8.5, style="italic", color="#8a8984")

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="D", ls="", ms=8, color=C_ENR,
                      label=f"pooled (grey if $p\\geq${alpha})"),
               Line2D([], [], marker="o", ls="-", ms=5, color=C_NEU, lw=1.4,
                      alpha=.8, label=f"each pooled protein, and their range")]
    if het is not None and len(het):
        handles.append(Line2D([], [], ls="-", lw=2, color=C_ENR, alpha=.55,
                              label="95% CI, random effects"))
    for j, name in enumerate(reference):
        handles.append(Line2D([], [], marker=REF_MARKERS[j % len(REF_MARKERS)], ls="",
                              ms=5.5, mfc="none", mec="#8a8984", mew=1.3,
                              label=f"{_short(name)} — not pooled"))
    if zeros:
        handles.append(Line2D([], [], marker="<", ls="", ms=7, color="#8a8984",
                              label="ratio 0 — target never sampled"))
    ncol = 2 if len(handles) > 3 else len(handles)
    ax.legend(handles=handles, loc="upper center",
              bbox_to_anchor=(0.5, -0.095 - 0.045 * leg_rows),
              ncol=ncol, frameon=False, fontsize=8.5)
    ax.set_title("Is the QTY target already sampled at the sites it replaces?"
                 + ("\nBUNDLED DEMO DATA — not a result" if USING_BUNDLED else ""),
                 fontsize=10.5, loc="left", pad=26)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    fig.savefig(FIGURES_DIR / fname)
    print(f"saved {FIGURES_DIR / fname}")
    return fig


plot_ratios({n: per_protein[n] for n in POOL}, pooled, het=het_pool,
            reference={n: per_protein[n] for n in REFERENCE});'''))

cells.append(md(r"""## 9 · What evolution actually puts there

The pooled ratio says a swap is unusual, not what the alternative is. This ranks every
amino acid by how often it appears at the sites of one native residue; anything above
the QTY target is a substitution the family tolerates more readily."""))

cells.append(code(r'''AAS = list("ACDEFGHIKLMNPQRSTVWY")


def alphabet_at(df, native):
    sites = df[df.aa == native]
    freq = {aa: 100 * float(np.mean([aa in s for s in sites.vset])) for aa in AAS}
    out = (pd.DataFrame({"aa": AAS, "pct": [freq[a] for a in AAS]})
             .sort_values("pct", ascending=False).reset_index(drop=True))
    out["rank"] = out.index + 1
    return out, len(sites)


def plot_alphabet(df, native, label, top=14):
    tbl, n = alphabet_at(df, native)
    target = QTY_MAP[native]
    show = tbl.head(top)
    tgt_rank = int(tbl.loc[tbl.aa == target, "rank"].iloc[0])
    pct = float(tbl.loc[tbl.aa == target, "pct"].iloc[0])
    cols = [C_DEP if a == target else (C_ENR if (a != native and rk < tgt_rank) else C_NEU)
            for a, rk in zip(show.aa, show["rank"])]
    fig, ax = plt.subplots(figsize=(5.6, 2.6))
    ax.bar(range(len(show)), show.pct, color=cols, width=.66)
    ax.set_xticks(range(len(show))); ax.set_xticklabels(show.aa)
    ax.set_ylabel(f"% of {native} sites"); ax.set_ylim(0, 104)
    ax.set_title(f"{label} — residues seen at the {n} {native} positions\n"
                 f"QTY target {target} ranks {tgt_rank}/20 at {pct:.0f}% of sites; "
                 f"blue = tolerated more readily", fontsize=9.5, loc="left")
    for i, (a, p) in enumerate(zip(show.aa, show.pct)):
        if a == target:
            ax.text(i, p + 3, f"{p:.0f}%", ha="center", fontsize=8.5,
                    color=C_DEP, fontweight="bold")
    safe = re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_")
    fig.savefig(FIGURES_DIR / f"fig_alphabet_{native}_{safe}.png")
    return tbl


# default to the deepest alignment — shallow ones give thin per-site alphabets
FOCUS = max(POOL or list(proteins),
            key=lambda n: (_depth(proteins[n]), len(proteins[n])))
print(f"FOCUS = {FOCUS}  (MSA depth {_depth_str(proteins[FOCUS])})\n")
tbl = plot_alphabet(proteins[FOCUS], "L", FOCUS)
print(tbl.head(10).to_string(index=False, formatters={"pct": "{:.0f}%".format}))'''))

cells.append(md(r"""## 10 · Triage

Splits the sites of one native residue by what its alignment column contains, in
decreasing order of evidence that the position tolerates a polar side chain. Only the
last class &mdash; a **strictly hydrophobic** alphabet &mdash; is an exemption candidate.

Two things the obvious version of this gets wrong:

* **A charged residue is evidence, and a stronger kind than polar.** A column containing
  D, E, K or R has already accepted a full formal charge at that position; glutamine is a
  smaller jump than any of them. Testing only S/T/N/Q/Y/H/C/W would send such a site to
  the exempt list, which inverts its meaning.
* **A small alphabet at a variable site is sampling, not constraint.** A position with a
  middling conservation grade can show four residues simply because the alignment is
  finite. The shortlist is therefore filtered on grade as well as on alphabet."""))

cells.append(code(r'''POLAR = set("STNQYHCW")
CHARGED = set("DEKR")
EXEMPT_MIN_GRADE = 8          # conserved sites only; 9 = most conserved


def triage(df, native, qty_map=QTY_MAP):
    target = qty_map[native]
    alt = {"L": "T", "I": "S", "V": "S", "F": "W"}.get(native)

    def call(v):
        if target in v:        return f"1 · {target} sampled — code as written"
        if alt and alt in v:   return f"2 · {alt} sampled, {target} not"
        if v & CHARGED:        return "3 · charged sampled — a larger jump than " + target
        if v & POLAR:          return "4 · other polar sampled"
        return "5 · strictly hydrophobic — exempt"

    d = df[df.aa == native].copy()
    d["call"] = d.vset.apply(call)
    d["swap"] = d.aa + d.seqres_pos.astype(str) + target
    d["charged_seen"] = d.vset.apply(lambda v: "".join(sorted(v & CHARGED)))
    summary = (d.groupby("call")
                 .agg(n=("seqres_pos", "size"), median_grade=("grade", "median"))
                 .assign(pct=lambda x: 100 * x.n / len(d)).reset_index()
                 .sort_values("call"))
    return d, summary


sites, summary = triage(proteins[FOCUS], "L")
_d = _depth(proteins[FOCUS])
print(f"{FOCUS} — all {len(sites)} leucine sites   (MSA depth {_depth_str(proteins[FOCUS])}, "
      f"median alphabet {proteins[FOCUS].n_variety.median():.0f} residues/site)")
if 0 <= _d < 100:
    print(f"   CAUTION: at depth {_d} the alphabets are thin, so the exempt class is\n"
          f"   inflated — a residue unobserved in {_d} sequences is not a residue\n"
          f"   evolution rejects.\n")
else:
    print()
print(summary[["call", "n", "pct", "median_grade"]]
      .to_string(index=False, formatters={"pct": "{:.0f}%".format,
                                          "median_grade": "{:.0f}".format}))

hydrophobic = sites[sites.call.str.startswith("5")]
exempt = hydrophobic[hydrophobic.grade >= EXEMPT_MIN_GRADE].sort_values(
    ["grade", "score"], ascending=[False, True])
dropped = len(hydrophobic) - len(exempt)

print(f"\nExemption candidates ({len(exempt)}): strictly hydrophobic alphabet AND "
      f"grade >= {EXEMPT_MIN_GRADE}\n")
print(exempt[["swap", "grade", "score", "n_variety", "variety", "in_tm"]]
      .to_string(index=False, formatters={"score": "{:+.2f}".format,
                                          "variety": lambda s: s.replace("|", ""),
                                          "in_tm": lambda b: "TM" if b else "— outside"}))
_out = exempt[~exempt.in_tm]
if len(_out):
    print(f"\n  {len(_out)} of these ({', '.join(_out.swap)}) fall OUTSIDE the "
          f"{tm_source[FOCUS]} transmembrane\n  regions. The code would never rewrite "
          f"them, so they cannot be exemptions from it.")
    print(f"  Acting list: {', '.join(exempt[exempt.in_tm].swap)} "
          f"({int(exempt.in_tm.sum())} sites).")
if dropped:
    loose = hydrophobic[hydrophobic.grade < EXEMPT_MIN_GRADE].sort_values("grade",
                                                                         ascending=False)
    print(f"\n{dropped} further site(s) have a hydrophobic alphabet but fall below "
          f"grade {EXEMPT_MIN_GRADE}.")
    print("A small alphabet at a variable site is a finite alignment, not a constraint:")
    print(loose[["swap", "grade", "score", "n_variety", "variety"]]
          .to_string(index=False, formatters={"score": "{:+.2f}".format,
                                              "variety": lambda s: s.replace("|", "")}))

charged = sites[sites.call.str.startswith("3")]
if len(charged):
    print(f"\n{len(charged)} site(s) have a CHARGED residue in the observed alphabet. "
          f"These are the\nopposite of exempt — the position has already accepted a "
          f"formal charge:")
    print(charged.sort_values("grade", ascending=False)
          [["swap", "grade", "score", "charged_seen", "variety"]]
          .to_string(index=False, formatters={"score": "{:+.2f}".format,
                                              "variety": lambda s: s.replace("|", "")}))'''))

cells.append(md(r"""## 11 · Write the results folder

Everything lands under this run's folder, with a manifest recording the exact input
files and their SHA-256 hashes &mdash; so a number in a manuscript can be traced back to
the ConSurf run that produced it."""))

cells.append(code(r'''per_site = []
for name, d in proteins.items():
    s = d[d.aa.isin(QTY_MAP)].copy()
    s["protein"] = name
    s["input_kind"] = klass[name]["kind"]
    s["pooled"] = name in POOL
    s["target"] = s.aa.map(QTY_MAP)
    s["precedented"] = [t in v for t, v in zip(s.target, s.vset)]
    s["swap"] = s.aa + s.seqres_pos.astype(str) + s.target
    per_site.append(s.drop(columns="vset"))
cols = list(dict.fromkeys(c for f in per_site for c in f.columns))
per_site = pd.concat([f.reindex(columns=cols) for f in per_site], ignore_index=True)

pooled.drop(columns=["all_below_1", "all_above_1"]).to_csv(
    TABLES_DIR / "pooled_by_swap.csv", index=False)
pd.concat([t.assign(protein=n) for n, t in per_protein.items() if len(t)],
          ignore_index=True).to_csv(TABLES_DIR / "by_protein.csv", index=False)
per_site.to_csv(TABLES_DIR / "per_site.csv", index=False)
per_site[~per_site.precedented].to_csv(TABLES_DIR / "unprecedented_sites.csv", index=False)
sites.drop(columns="vset").to_csv(TABLES_DIR / f"triage_{FOCUS[:4]}_L.csv", index=False)
kt.to_csv(TABLES_DIR / "input_classification.csv", index=False)
if cmp_ is not None:
    cmp_.to_csv(TABLES_DIR / "scope_whole_vs_tm.csv", index=False)
if sens is not None:
    sens.to_csv(TABLES_DIR / "tm_boundary_sensitivity.csv", index=False)
    sens_summary.to_csv(TABLES_DIR / "tm_boundary_summary.csv", index=False)
if het_pool is not None:
    het_pool.to_csv(TABLES_DIR / "heterogeneity_pool.csv", index=False)
if agree is not None:
    agree.to_csv(TABLES_DIR / "pool_agreement.csv")

run_manifest = dict(
    run_id=RUN_ID,
    created=datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
    qty_map=QTY_MAP,
    ranges={k: list(v) for k, v in RANGES.items()},
    excluded=sorted(EXCLUDE),
    data_roots=[str(r) for r in DATA_ROOTS],
    drive_mounted=bool(DRIVE_OK),
    data_source="bundled demo" if USING_BUNDLED else "discovered files",
    inputs=manifest,
    classification={n: dict(kind=k["kind"], livf=round(k["livf"], 1),
                            qty=round(k["qty"], 1), tm_segs=int(k["tm_segs"]),
                            self_bg=round(k["self_bg"], 1), forced=bool(k["forced"]))
                    for n, k in klass.items()},
    scope=SCOPE,
    tm_source=tm_source,
    tm_regions={n: [list(x) for x in proteins[n].attrs.get("tm_regions", [])]
                for n in proteins},
    scope_comparison=(None if cmp_ is None else
                      {r.swap: dict(whole=round(float(r.ratio_whole), 3),
                                    tm=round(float(r.ratio_tm), 3),
                                    shift=round(float(r.tm_over_whole), 3))
                       for _, r in cmp_.iterrows()}),
    tm_boundary_sensitivity=(None if sens_summary is None else
        {r.swap: dict(lo=round(float(r.lo), 3), median=round(float(r.med), 3),
                      hi=round(float(r.hi), 3), n_defs=int(r.n_defs),
                      n_sig_depleted=int(r.n_sig_depleted),
                      n_sig_enriched=int(r.n_sig_enriched),
                      min_p=float(r.min_p), contradictory=bool(r.contradictory))
         for _, r in sens_summary.iterrows()}),
    pooled_proteins=sorted(POOL),
    reference_proteins=sorted(REFERENCE),
    heterogeneity=({r.swap: dict(I2=round(float(r.I2), 1), re=round(float(r.re), 3),
                                 ci=[round(float(r.lo), 3), round(float(r.hi), 3)])
                    for _, r in het_pool.iterrows()} if het_pool is not None else
                   f"not computed: k = {len(POOL)} < {MIN_K_FOR_META}"),
    n_substitutable_sites=int(len(per_site)),
    n_unprecedented=int((~per_site.precedented).sum()),
    n_substitutable_sites_pooled=int(per_site.pooled.sum()),
    n_unprecedented_pooled=int((per_site.pooled & ~per_site.precedented).sum()),
    pooled={r.swap: dict(n=int(r.n_sites), ratio=round(float(r.ratio), 4),
                         ci=[round(float(r.ci_lo), 4), round(float(r.ci_hi), 4)],
                         detectable_outside=[round(float(r.det_lo), 4),
                                             round(float(r.det_hi), 4)],
                         p=float(r.p), verdict=r.verdict)
            for _, r in pooled.iterrows()},
    versions=dict(pandas=pd.__version__, numpy=np.__version__,
                  matplotlib=mpl.__version__, python=sys.version.split()[0]),
)
(RESULTS_DIR / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2))

_pooled_sites = per_site[per_site.pooled]
print(f"{len(per_site)} substitutable sites across {len(proteins)} input(s); "
      f"{len(_pooled_sites)} of them in the {len(POOL)} pooled protein(s); "
      f"{(~_pooled_sites.precedented).sum()} of those unprecedented\n")
print(f"results folder: {RESULTS_DIR}\n")
if USING_BUNDLED:
    print("!! these outputs are from BUNDLED DEMO DATA, not your files !!\n")
for f in sorted(RESULTS_DIR.rglob("*")):
    if f.is_file():
        print(f"  {f.relative_to(RESULTS_DIR)}  ({f.stat().st_size:,} bytes)")'''))

cells.append(md(r"""## 12 · Adding proteins

1. Run ConSurf, keep `consurf_grades.txt` (either layout).
2. Put it under one of the `DATA_ROOTS`. The label comes from the folder tree: opaque
   job folders (`consurf_1761475466`) and generic buckets (`RMSD`, `org`) are skipped
   over, so `EAA1plos/consurf_1761475466/` and `EAA1plos/org/` both label as
   `EAA1plos`. Override anything in `NAME_MAP`.
3. Re-run from the top.

### Getting ConSurf runs without the web form

There is **no public API for submitting a ConSurf job.** Three routes exist, and they are
not interchangeable:

1. **ConSurf standalone** &mdash; the authors' own Python pipeline, run locally or in this
   notebook. This is the only scriptable way to *compute* new runs, and the only one that
   lets you fix the alignment depth, the homolog database and the search parameters across
   every protein. For a pooled comparison that control is the whole game: it is what would
   let VPAT be re-run at the same depth as EAA1, and what would produce a native EAA2 run.
2. **ConSurf-DB** &mdash; precomputed profiles for PDB chains, fetchable (the helper below).
   No UniProt lookup and, as of the 2023 release, no AlphaFold models, so a protein without
   an experimental structure is simply absent.
3. The web form, one protein at a time.

The trap in 2 is that ConSurf-DB builds its alignments with its own fixed pipeline, which
is **not** the pipeline behind a run you submitted yourself. Mixing the two inside one pool
reintroduces exactly the between-protein heterogeneity the rest of this notebook works to
exclude. Use ConSurf-DB for every protein or for none."""))

cells.append(code(r'''CONSURFDB = "https://consurfdb.tau.ac.il"


def consurfdb_entry_ids(pdb_id, timeout=30):
    """{chain: entry_id} for a PDB entry. Scrapes the chain-selection page — this is
    not a documented API and the markup can change without notice."""
    url = f"{CONSURFDB}/scripts/chain_selection.php?pdb_ID={pdb_id.upper()}"
    req = urllib.request.Request(url, headers={"Accept": "text/html"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        html = r.read().decode(errors="replace")
    out = {}
    for value, label in re.findall(r'<option[^>]*value="([^"]+)"[^>]*>([^<]*)</option>',
                                   html, re.I):
        m = re.search(r"\b([A-Za-z0-9])\b", label) or re.search(r"\b([A-Za-z0-9])\b", value)
        if m and value.strip():
            out.setdefault(m.group(1).upper(), value.strip())
    return out


def fetch_consurfdb(pdb_id, chain, dest_dir, timeout=60):
    """Download one ConSurf-DB summary into dest_dir so discovery picks it up."""
    ids = consurfdb_entry_ids(pdb_id, timeout)
    if not ids:
        raise LookupError(f"{pdb_id}: no ConSurf-DB entries (no experimental structure?)")
    key = chain.upper()
    if key not in ids:
        raise LookupError(f"{pdb_id} chain {chain} not in ConSurf-DB; have {sorted(ids)}")
    url = f"{CONSURFDB}/DB/{ids[key]}/consurf_summary.txt"
    dest = Path(dest_dir) / f"{pdb_id.upper()}_{key}" / "consurf_grades.txt"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=timeout) as r:
        text = r.read().decode(errors="replace")
    dest.write_text(text)
    print(f"  {pdb_id} {key} -> {dest}  ({len(text):,} bytes, entry {ids[key]})")
    return dest


# Uncomment and point at a DATA_ROOT, then re-run the notebook from section 4.
# Remember: do not mix ConSurf-DB entries with your own submitted runs in one pool.
#
# for pdb, ch in [("5LLU", "A")]:
#     try:
#         fetch_consurfdb(pdb, ch, DATA_ROOTS[0])
#     except Exception as exc:
#         print(f"  {pdb} {ch}: {type(exc).__name__}: {exc}")
print("ConSurf-DB helper defined. See the commented block to use it.")'''))

cells.append(md(r"""### Two kinds of duplicate, and only one is caught by hashing

**Byte-identical copies** of one file in several folders are detected at load and
dropped automatically.

**The same protein run twice** &mdash; a 50-sequence and a 150-sequence alignment of the
same sequence &mdash; is not byte-identical, so hashing misses it. The load step now
hashes the *amino-acid sequence* as well; identical sequences are reported and the
deepest alignment kept. Set `AUTO_DEDUP_SEQUENCES = False` to pool them anyway, but
understand that doing so counts one protein twice and inflates every pooled `n` and
*p*-value.

### Alignment depth

Depth drives two different things. The **pooled ratio** is a ratio to that protein's own
background, so it survives a shallow alignment. The **exemption list** does not: at
depth 50 a site's observed alphabet is about 5 residues, at depth 150 about 8, so the
strictly-hydrophobic class grows from roughly a fifth to roughly two fifths of leucine
sites on the same protein. A residue absent from 50 sequences is not a residue
evolution rejects. Prefer depth &ge;150 for anything that becomes a candidate list, and
treat the exempt class from a shallow run as a shortlist to re-run deeper, not a result.
The grade filter (`EXEMPT_MIN_GRADE`) is the second guard on the same problem: it keeps
a site off the list unless conservation independently says the position is constrained.

### Scope

`QTY_MAP` is editable &mdash; set `"L": "T"` to score that variant of the code instead.

What the mapping *means* depends on the input, which is why section 4 classifies every
file before anything is pooled. The QTY code is applied to the transmembrane segments of
a native membrane protein, so only such a protein can answer a precedence question. Two
failure modes are silent without that check:

* **A QTY-designed sequence.** Its substituted sites no longer carry the native residue,
  so they drop out of the denominator and only the untouched soluble-domain sites remain;
  meanwhile the target residue appears in its own alignment column at every designed
  position and inflates the background it is being compared against. Both errors push the
  ratio toward 1 from opposite directions, and nothing in the file looks wrong.
* **A soluble protein.** TLR5S is the clear case: a leucine-rich-repeat ectodomain with
  no transmembrane segment at all, whose leucines are the structural L of the LRR
  consensus. It shows the most extreme L&rarr;Q depletion in the set for reasons that have
  nothing to do with bilayers. That makes it a genuinely useful **reference** &mdash; it
  shows the depletion tracks buried structural leucine rather than membrane leucine
  &mdash; but it is not a transporter, and averaging it into a transporter number merges
  two different claims.

Add a protein and the classifier will place it; `FORCE_KIND` overrules it if you
disagree."""))

nb = {"cells": cells,
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python",
                                  "name": "python3"},
                   "language_info": {"name": "python", "version": "3.11"},
                   "colab": {"provenance": [], "toc_visible": True}},
      "nbformat": 4, "nbformat_minor": 5}

import pathlib
_OUT = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "qty_precedence.ipynb"
_OUT.parent.mkdir(parents=True, exist_ok=True)
with open(_OUT, "w") as fh:
    json.dump(nb, fh, indent=1)
print("wrote notebook with", len(cells), "cells")
