# Repository Guidelines

## Project Structure & Module Organization
- `src/`: Rust CLI entry (`src/main.rs`) with subcommands `prepare` and `analyze`.
- `book/`: mdBook docs. Sources in `book/src`, config in `book/book.toml`.
- Generated data in repo root: `uniprot_sprot.*`, `clusters*` (ignored by `.gitignore`).
- Uses `needletail` for FASTA handling; DIAMOND runs via external binary.
- `config.example.toml`: sample `analyze` configuration with `scoring.*` overrides.
- ECS orchestration lives in `src/main.rs` via Bevy `Startup`/`Update` schedules; jobs, tasks, and results are tracked as resources.

## Build, Test, and Development Commands
- `cargo build` – compile debug; `cargo build --release` – optimized binary.
- `cargo run -- prepare` – download SwissProt, build DIAMOND DB, cluster/realign/recluster.
- `cargo run -- analyze <clusters.recluster>` – analyze clusters (expects artifacts from `prepare`).
- `cargo fmt --all` – format; `cargo clippy -- -D warnings` – lint strictly.
- Docs: `mdbook build book` or `mdbook serve book -o`.
- `cargo test` – run unit tests for scoring/config logic.
- `cargo run -- analyze --config config.example.toml` – example run with MAFFT + scoring overrides; override with flags like `--mafft-bin` or `--reference-fasta`.

### Pixi environment
- Prefer running tools via Pixi to ensure DIAMOND/HMMER are on PATH:
  - `pixi run cargo build`
  - `pixi run cargo test`
  - `pixi run cargo run -- analyze --config config.example.toml`
  - `pixi run diamond --version` / `pixi run diamond blastp ...`
  - `pixi run hmmscan --help` (for future Pfam usage)

## Coding Style & Naming Conventions
- Rust 2021; 4‑space indent; `snake_case` for modules/functions/vars; `PascalCase` for types.
- Keep CLI parsing in `main.rs`; factor analysis/helpers into small pure functions.
- Logging via `env_logger`: `RUST_LOG=info cargo run -- …` (use `debug` when investigating).
- Bevy crates (`bevy_app`, `bevy_ecs`, `bevy_tasks`) pinned to `0.17.x` for the ECS scheduler; keep systems lightweight and deterministic.

## Testing Guidelines
- Use Rust’s built‑in tests (`#[cfg(test)]` modules). Run with `cargo test`.
- Add tests for parsing, file I/O helpers, and any new analysis logic.
- Prefer deterministic fixtures; avoid relying on network or DIAMOND in unit tests.
- Keep tests fast by stubbing external tools; current coverage includes scoring + config overrides.
- When adding MAFFT/HMMER integrations, gate them behind config and add unit tests around parsing/metrics only (no external calls).

## Commit & Pull Request Guidelines
- Commits: short, imperative subject (e.g., “Add analysis metrics”); body explains rationale when needed.
- Commit cadence: make a commit immediately after implementing a new feature or fix (small, focused commits). Group only closely related changes; avoid mixing refactors with feature commits.
- Conventional commits are encouraged for clarity:
  - `feat: add DIAMOND pre-run caching`
  - `fix: handle empty DIAMOND output without panic`
  - `docs: expand metrics page`
  - `refactor: extract ecs scheduler into module`
  - `test: add smoke test for analyze`
  - `chore: bump dependencies`
- PRs include: clear description, linked issues, CLI usage examples/output, and docs updates under `book/` when behavior changes.
- If Bevy behavior changes or you hit recurring ECS quirks, add notes to `BEVY_GUIDE.md` in the PR.

## Security & Configuration Tips
- Prereqs: Rust stable, DIAMOND available in `PATH`. Needletail is provided via Cargo (no system install).
- `prepare` downloads large files and writes multi‑GB artifacts—do not commit outputs.
- Long runs currently use 16 threads in code; propose flags if you need configurability.
- Optional: MAFFT in `PATH` (or set `--mafft-bin`) to enable conserved-region evidence; ensure `reference_fasta` matches the DIAMOND database used in config.

## Agent‑Specific Instructions
- Scope: applies to the entire repo. Keep patches minimal and focused; don’t reformat unrelated files or add generated data.
- Prioritize modularity and maintainability: prefer small, testable Rust modules (e.g., `src/diamond.rs`, `src/ecs.rs`, `src/taxonomy.rs`) over monolithic files. Extract sub‑libraries into modules when it improves separation of concerns, reuse, and testability. Follow best practices: clear ownership of responsibilities, narrow interfaces, and zero shared mutable state unless behind resources or channels.
- Log MCP pain points in `MCP_PAIN_POINTS.md` (tool name + error/symptom + workaround).
- Prefer `graphrag` for file discovery, edits, and context gathering whenever possible (use it as the default editing workflow before falling back to ad hoc shell tools).
- When reporting long-running DIAMOND status, always verify the process is active with `pgrep` or `ps` before saying it is still running. Prefer using Codex background process mode for long runs and only report "still running" after confirmation.

### Graphrag Editing Template
Use this short spiel when starting edits:

"I'll use graphrag to locate and edit the relevant files for this change. If anything is unclear, I'll confirm the target before applying edits."

### ECS + Streaming Emit Notes
- ECS is the backbone: DIAMOND intake, MAFFT, HMMER, structvar, and now JSON/CSV rendering all run through Bevy systems. Any new heavy stage should follow the same pattern (seed jobs → async handles → collect → flush).
- Rendering is streaming: `run_render_pipeline` writes JSONL/CSV per gene as soon as its components are ready. No giant `HashMap` buffering; always prefer components/resources over ad‑hoc maps.
- Keep systems deterministic: cap async pools with `max_jobs`, maintain ordered flush via `next_index`, and despawn entities after flushing to keep memory flat.
- Diagnostics belong in ECS resources too (e.g., domains_arch_debug, structvar summaries). When adding new per-gene diagnostics, aggregate via systems instead of post-processing the whole dataset.
- When touching Bevy APIs, consult upstream docs for 0.17 patterns—avoid direct entity handles after despawn; use resources/queues instead.
