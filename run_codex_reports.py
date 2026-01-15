#!/usr/bin/env python3
"""
Walks reports/ recursively, runs codex on each report file, and moves it to
completed_reports/ preserving relative directory structure.

Usage:
  python run_codex_reports.py
  python run_codex_reports.py --dry-run

Notes:
  - Requires `codex` CLI on PATH.
  - Files are only moved on successful codex exit (code 0).
  - Designed to be run from the repo root.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

#CODEX_MODEL = "mistralai/Devstral-Small-2-24B-Instruct-2512"
CODEX_ARGS = (
    "exec",
#    "--oss",
    "--dangerously-bypass-approvals-and-sandbox",
    "-s",
    "danger-full-access",
#    "-m",
#    CODEX_MODEL,
)


def iter_files(root: Path) -> list[Path]:
    return sorted([p for p in root.rglob("*") if p.is_file()])


def repo_preview(repo_root: Path, limit: int = 30) -> str:
    entries = sorted(repo_root.iterdir(), key=lambda p: p.name)
    preview = []
    for entry in entries[:limit]:
        name = f"{entry.name}/" if entry.is_dir() else entry.name
        preview.append(name)
    if len(entries) > limit:
        preview.append("...")
    return " ".join(preview)


def read_agents_instructions(repo_root: Path) -> str:
    agents_path = repo_root / "AGENTS.md"
    if not agents_path.exists():
        return "AGENTS.md not found."
    try:
        return agents_path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        return f"AGENTS.md unreadable: {exc}"


def build_prompt(
    file_path: Path,
    completed_dir: Path,
    repo_root: Path,
    target_path: Path,
    target_exists: bool,
    repo_listing: str,
) -> str:
    report_text = file_path.read_text(encoding="utf-8", errors="replace")
    agents_text = read_agents_instructions(repo_root)
    return (
        "You are working in the AnnoQC repo. "
        "Apply the fixes described by the report if they are still valid, safe, and "
        "do not regress other behavior. "
        f"Report file: {file_path}. "
        f"Target file (report path without 'reports/' and without the .md suffix): {target_path}. "
        f"Target exists: {'yes' if target_exists else 'no'}. "
        f"Repo top-level preview: {repo_listing}. "
        "Each report mirrors the repo path without the leading 'reports/' prefix. "
        "If the target file no longer exists, mark the report as done immediately. "
        f"After finishing, move the report to {completed_dir} preserving the relative path. "
        "Be mindful that files may have changed since the report was generated. "
        "Run only the smallest relevant tests (if any). Do not run fmt for this pass. "
        "Make a concise git commit per report. Only add the files you have touched this session, most likely no more than the 1 we adjusted. "
        "Do not reformat unrelated files or add generated data. "
        "Use pixi for builds/tests when relevant (e.g., `pixi run cargo test`). "
        "Keep ECS systems deterministic and streaming; follow existing Bevy 0.17 patterns. "
        "If you hit tool friction, log it in MCP_PAIN_POINTS.md. "
        "Speed: do not scan the whole repo; only open files referenced by the report and their direct dependencies. "
        "Use ripgrep for targeted searches. Avoid heavy or full test suites. "
        "Follow AGENTS.md instructions; they are authoritative. "
        "AGENTS.md contents begin below:\n\n"
        f"{agents_text}\n"
        "AGENTS.md contents end above.\n\n"
        "Report contents begin below:\n\n"
        f"{report_text}\n"
        "Report contents end above."
    )


def run_codex_on_file(
    file_path: Path,
    completed_dir: Path,
    repo_root: Path,
    target_path: Path,
    target_exists: bool,
    repo_listing: str,
) -> int:
    prompt = build_prompt(
        file_path,
        completed_dir,
        repo_root,
        target_path,
        target_exists,
        repo_listing,
    )
    cmd = ["codex", *CODEX_ARGS, prompt]
    proc = subprocess.run(cmd)
    return proc.returncode


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reports-dir",
        type=Path,
        default=Path("reports"),
        help="Root directory to scan (default: reports).",
    )
    parser.add_argument(
        "--completed-dir",
        type=Path,
        default=Path("completed_reports"),
        help="Root directory to move completed files into (default: completed_reports).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print actions without running codex or moving files.",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=None,
        help="Optional cap on number of files processed.",
    )
    args = parser.parse_args(argv)

    repo_root = Path(__file__).resolve().parent
    reports_dir = (repo_root / args.reports_dir).resolve()
    completed_dir = (repo_root / args.completed_dir).resolve()

    if not reports_dir.exists():
        print(f"reports dir not found: {reports_dir}", file=sys.stderr)
        return 2

    files = iter_files(reports_dir)
    if args.max_files is not None:
        files = files[: args.max_files]

    print(f"Found {len(files)} files under {reports_dir}")

    repo_listing = repo_preview(repo_root)

    for i, file_path in enumerate(files, start=1):
        print(f"[{i}/{len(files)}] Processing {file_path}")

        rel = file_path.relative_to(reports_dir)
        completed_path = completed_dir / rel
        target_rel = rel.with_suffix("") if rel.suffix == ".md" else rel
        target_path = (repo_root / target_rel).resolve()
        target_exists = target_path.exists()

        # Skip work when the report already exists in the completed directory.
        if completed_path.exists():
            print(
                f"  already completed at {completed_path}; skipping",
                file=sys.stderr,
            )
            continue

        if not target_exists:
            print(
                f"  target file missing at {target_path}; marking report completed",
                file=sys.stderr,
            )
            completed_path.parent.mkdir(parents=True, exist_ok=True)
            file_path.replace(completed_path)
            print(f"  moved to {completed_path}")
            continue

        if args.dry_run:
            prompt = build_prompt(
                file_path,
                completed_dir,
                repo_root,
                target_path,
                target_exists,
                repo_listing,
            )
            cmd_preview = f"codex {' '.join(CODEX_ARGS)} \"{prompt}\""
            print(f"  would run: {cmd_preview}")
            print(f"  would move to: {completed_dir / rel}")
            continue

        rc = run_codex_on_file(
            file_path,
            completed_dir,
            repo_root,
            target_path,
            target_exists,
            repo_listing,
        )
        if rc == 0:
            expected_dest = completed_dir / rel
            if file_path.exists():
                print(
                    f"  codex succeeded but file still present at {file_path}; "
                    f"expected move to {expected_dest}",
                    file=sys.stderr,
                )
            else:
                print(f"  codex succeeded; file moved to {expected_dest}")
        else:
            print(f"  codex failed (exit {rc}); leaving file in place", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
