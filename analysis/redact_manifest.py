"""Reduce absolute input paths in a run manifest to bare filenames.

The manifest records where each ConSurf file sat on the machine that produced
the run. Those paths say nothing a reader needs — the SHA-256 digest is what
anchors provenance — while exposing the directory structure of unrelated work.
This rewrites `inputs[].path` to `inputs[].file` and drops `data_roots`,
leaving every digest, count and parameter untouched.

Usage: python3 redact_manifest.py <run_dir> [...]
"""
import json
import os
import sys


def redact(run_dir):
    path = os.path.join(run_dir, "run_manifest.json")
    with open(path) as fh:
        m = json.load(fh)
    changed = False
    for entry in m.get("inputs", []):
        if "path" in entry:
            entry["file"] = os.path.basename(entry.pop("path"))
            changed = True
    if "data_roots" in m:
        m["data_roots"] = "redacted; see inputs[].file and inputs[].sha256"
        changed = True
    if not changed:
        print(f"{path}: already redacted")
        return
    with open(path, "w") as fh:
        json.dump(m, fh, indent=2)
        fh.write("\n")
    print(f"{path}: redacted {len(m.get('inputs', []))} input paths")


if __name__ == "__main__":
    targets = sys.argv[1:] or [
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                     "results", "run_2026-11-04_161147")]
    for t in targets:
        redact(t)
