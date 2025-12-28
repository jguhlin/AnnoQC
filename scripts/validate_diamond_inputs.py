#!/usr/bin/env python3
"""Validate FASTA/query inputs and DIAMOND databases before running blastp.

Checks performed:
  * File existence & non-zero size
  * Optional gzip integrity
  * FASTA structural sanity (headers, sequence lines)
  * Non-IUPAC residues for protein FASTA (or nucleotide when requested)
  * Warn if protein FASTA appears nucleotide-like (or vice versa)
  * Optional DIAMOND dbinfo probe with hash reporting

Usage:
  pixi run python scripts/validate_diamond_inputs.py --fasta a.fa --db foo.dmnd
"""
import argparse
import gzip
import os
import re
import subprocess
import sys
from collections import Counter
from typing import Dict, List, Tuple

PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBZXJUO*-.?")
NUCLEOTIDE_ALPHABET = set("ACGTURYKMSWBDHVNX-.")
FASTA_HEADER = re.compile(r"^>.+")


def open_text(path: str):
    if path.endswith(".gz"):
        handle = gzip.open(path, "rt", encoding="utf-8", errors="replace")
        try:
            handle.read(1)
            handle.seek(0)
        except (OSError, EOFError):
            handle.close()
            raise
        return handle
    return open(path, "rt", encoding="utf-8", errors="replace")


def detect_mode(expect: str) -> Tuple[set, set]:
    if expect == "protein":
        return PROTEIN_ALPHABET, NUCLEOTIDE_ALPHABET
    return NUCLEOTIDE_ALPHABET, PROTEIN_ALPHABET


def analyze_fasta(path: str, expect: str, warn_threshold: float) -> Dict:
    """Scan a FASTA file for structure, alphabet validity, and mode warnings."""
    stats = {
        "path": path,
        "headers": 0,
        "sequences": 0,
        "total_residues": 0,
        "invalid_lines": [],
        "warnings": [],
        "errors": [],
    }
    allowed, alternate = detect_mode(expect)
    invalid_re = re.compile(f"[^{re.escape(''.join(sorted(allowed)))}]")
    nucleotide_like = 0
    seq_chars = Counter()
    current_seq_len = 0

    try:
        with open_text(path) as handle:
            for line_no, raw in enumerate(handle, 1):
                line = raw.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith(">"):
                    stats["headers"] += 1
                    if current_seq_len:
                        if current_seq_len == sum(seq_chars.values()):
                            share = sum(seq_chars[c] for c in alternate) / max(current_seq_len, 1)
                            if share >= warn_threshold:
                                nucleotide_like += 1
                        seq_chars.clear()
                        current_seq_len = 0
                    continue
                if stats["headers"] == 0:
                    stats["errors"].append(f"line {line_no}: sequence data encountered before first header")
                    break
                stats["sequences"] += 1
                line_upper = line.upper()
                seq_chars.update(line_upper)
                current_seq_len += len(line_upper)
                stats["total_residues"] += len(line_upper)
                invalid_chars = invalid_re.findall(line_upper)
                if invalid_chars:
                    stats["invalid_lines"].extend((line_no, ch) for ch in invalid_chars)
    except (OSError, EOFError) as exc:
        stats["errors"].append(f"failed to read FASTA: {exc}")
        return stats

    if stats["headers"] == 0:
        stats["errors"].append("no FASTA headers detected")
    if stats["total_residues"] == 0:
        stats["errors"].append("no sequence data detected")
    if stats["invalid_lines"]:
        bad_preview = ", ".join(f"line {ln}:{ch}" for ln, ch in stats["invalid_lines"][:10])
        stats["errors"].append(
            f"detected {len(stats['invalid_lines'])} residue(s) outside expected {expect} alphabet (e.g., {bad_preview})"
        )
    if nucleotide_like and expect == "protein":
        stats["warnings"].append(
            f"{nucleotide_like} sequence(s) look nucleotide-like (>={warn_threshold:.0%} ATGCN content)"
        )
    if stats["headers"] and stats["sequences"] == 0:
        stats["warnings"].append("headers present but no sequence lines found")
    return stats


def check_db(db_path: str, diamond_bin: str) -> Tuple[bool, str]:
    """Run `diamond dbinfo` and capture output for the given database."""
    try:
        proc = subprocess.run(
            [diamond_bin, "dbinfo", "--db", db_path],
            capture_output=True,
            text=True,
            shell=False,
            check=True,
        )
    except FileNotFoundError:
        return False, f"diamond binary '{diamond_bin}' not found"
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        stdout = (exc.stdout or "").strip()
        return False, stderr or stdout or f"diamond dbinfo failed with exit code {exc.returncode}"
    return True, proc.stdout.strip()


def main():
    parser = argparse.ArgumentParser(description="Validate DIAMOND inputs")
    parser.add_argument("--fasta", action="append", help="query FASTA to validate", required=False)
    parser.add_argument("--db", help="DIAMOND database (.dmnd) to probe", required=False)
    parser.add_argument(
        "--expect",
        choices=["protein", "nucleotide"],
        default="protein",
        help="expected alphabet for FASTA queries (default: protein)",
    )
    parser.add_argument("--diamond-bin", default="diamond", help="diamond executable (default: diamond)")
    parser.add_argument("--warn-threshold", type=float, default=0.9, help="fraction for nucleotide-like warning")
    args = parser.parse_args()

    had_error = False
    if args.fasta:
        for fasta in args.fasta:
            print(f"→ Checking FASTA {fasta}")
            if not os.path.exists(fasta):
                print(f"  ERROR: file not found")
                had_error = True
                continue
            if os.path.getsize(fasta) == 0:
                print(f"  ERROR: file is empty")
                had_error = True
                continue
            stats = analyze_fasta(fasta, args.expect, args.warn_threshold)
            for err in stats["errors"]:
                print(f"  ERROR: {err}")
            for warn in stats["warnings"]:
                print(f"  WARN: {warn}")
            print(
                f"  Summary: {stats['headers']} headers, {stats['total_residues']} residues, "
                f"{stats['sequences']} sequence lines"
            )
            if stats["errors"]:
                had_error = True
    else:
        print("No FASTA paths provided; skipping FASTA validation")

    if args.db:
        print(f"→ Checking DIAMOND DB {args.db}")
        if not os.path.exists(args.db):
            print("  ERROR: database file not found")
            had_error = True
        else:
            ok, msg = check_db(args.db, args.diamond_bin)
            if ok:
                print("  dbinfo output:\n" + "\n".join(f"    {line}" for line in msg.splitlines()))
            else:
                print(f"  ERROR: diamond dbinfo failed: {msg}")
                had_error = True
    else:
        print("No DIAMOND DB path provided; skipping dbinfo probe")

    sys.exit(1 if had_error else 0)


if __name__ == "__main__":
    main()
