"""Fetch transmembrane regions from UniProt and verify them against the ConSurf query.

The whole point is the verification. Pulling coordinates from an accession that does
not correspond, residue for residue, to the sequence ConSurf actually scored silently
shifts every boundary — which is worse than the hydropathy guess, because it looks
authoritative.
"""
import json
import urllib.request
import urllib.parse

UNIPROT = "https://rest.uniprot.org/uniprotkb"
MEMBRANE_TYPES = {"Transmembrane", "Intramembrane"}
TIMEOUT = 30


def _get_json(url, timeout=TIMEOUT):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def uniprot_entry(accession, timeout=TIMEOUT):
    """One entry by accession."""
    return _get_json(f"{UNIPROT}/{accession}.json", timeout)


def uniprot_search(gene, organism_id=9606, reviewed=True, timeout=TIMEOUT):
    """Best reviewed hit for a gene symbol. Returns (accession, entry) or (None, None)."""
    q = f"gene_exact:{gene} AND organism_id:{organism_id}"
    if reviewed:
        q += " AND reviewed:true"
    url = (f"{UNIPROT}/search?query={urllib.parse.quote(q)}"
           f"&format=json&size=5")
    hits = _get_json(url, timeout).get("results", [])
    if not hits:
        return None, None
    return hits[0]["primaryAccession"], hits[0]


def membrane_features(entry):
    """[(start, end, type, description)] for Transmembrane and Intramembrane."""
    out = []
    for f in entry.get("features", []):
        if f.get("type") not in MEMBRANE_TYPES:
            continue
        loc = f.get("location", {})
        s, e = loc.get("start", {}).get("value"), loc.get("end", {}).get("value")
        if s is None or e is None:
            continue                       # fuzzy/unknown position: unusable
        out.append((int(s), int(e), f["type"], f.get("description", "")))
    return sorted(out)


def entry_sequence(entry):
    return entry.get("sequence", {}).get("value", "")


def reconcile(consurf_seq, uniprot_seq, consurf_first_pos=1):
    """How do ConSurf's SEQRES numbers map onto UniProt's?

    Returns (offset, note). A UniProt position p corresponds to ConSurf position
    p + offset. None means the sequences cannot be reconciled and the annotation
    must not be used.
    """
    if not consurf_seq or not uniprot_seq:
        return None, "empty sequence"
    if consurf_seq == uniprot_seq:
        return consurf_first_pos - 1, "identical sequences"
    if consurf_seq in uniprot_seq:
        i = uniprot_seq.index(consurf_seq)          # 0-based start in UniProt
        if uniprot_seq.count(consurf_seq) > 1:
            return None, "query matches UniProt at more than one position"
        return consurf_first_pos - 1 - i, (
            f"ConSurf query is UniProt {i + 1}-{i + len(consurf_seq)} "
            f"({len(uniprot_seq) - len(consurf_seq)} residues not scored)")
    if uniprot_seq in consurf_seq:
        i = consurf_seq.index(uniprot_seq)
        return consurf_first_pos - 1 + i, (
            f"UniProt sequence is a subrange of the ConSurf query at {i + 1}")
    # last resort: same length, report identity so the caller can judge
    if len(consurf_seq) == len(uniprot_seq):
        ident = sum(a == b for a, b in zip(consurf_seq, uniprot_seq)) / len(consurf_seq)
        return (consurf_first_pos - 1, f"same length, {100 * ident:.1f}% identical") \
            if ident >= 0.98 else (None, f"same length but only {100 * ident:.1f}% identical")
    return None, (f"sequences differ (ConSurf {len(consurf_seq)} aa, "
                  f"UniProt {len(uniprot_seq)} aa) and neither contains the other")


def regions_for(consurf_seq, entry, consurf_first_pos=1, include_intramem=True):
    """(regions, note, features). regions is [] when it cannot be trusted."""
    useq = entry_sequence(entry)
    offset, note = reconcile(consurf_seq, useq, consurf_first_pos)
    feats = membrane_features(entry)
    if offset is None:
        return [], note, feats
    lo, hi = consurf_first_pos, consurf_first_pos + len(consurf_seq) - 1
    regions = []
    for s, e, typ, _ in feats:
        if typ == "Intramembrane" and not include_intramem:
            continue
        a, b = s + offset, e + offset
        a, b = max(a, lo), min(b, hi)
        if a <= b:
            regions.append((a, b))
    return sorted(regions), note, feats
