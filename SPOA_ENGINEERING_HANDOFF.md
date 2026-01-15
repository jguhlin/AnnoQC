# SPOA Alignment Integration — Engineering Handoff

This document describes how to integrate the `spoa` crate as an in‑process
multiple sequence alignment backend in AnnoQC, alongside (but not replacing)
the existing MAFFT-based alignment pipeline.

The goal is to let us experiment with SPOA for speed and reduced I/O while
keeping the current MAFFT behaviour as a stable, documented baseline. The
design below assumes we are still in “development mode” and can tolerate minor
score shifts, as long as behaviour is well understood and controlled.

---

## High‑Level Goals

- **Optional in‑process aligner**
  - Introduce a new aligner backend using `spoa` (Rust crate) that runs wholly
    in process (no temp FASTA, no subprocess), fed directly from ECS.
  - Keep MAFFT as the default aligner; SPOA is opt‑in via CLI/config.

- **Same logical output contract**
  - `mafft.rs` currently returns an alignment representation that downstream
    systems consume for:
    - `start_concordance` / `end_concordance` (termini concordance pillar),
    - conserved fraction, gap runs, missing exon / retained intron heuristics.
  - The SPOA path must produce *the same* alignment struct (or a minimal superset)
    so everything downstream is unchanged.

- **Minimal impact on ECS orchestration**
  - Reuse the existing MAFFT ECS pipeline (`run_alignment_pipeline` etc.) and
    simply swap out “how we build the MSA for this panel” depending on an
    `Aligner` enum.
  - Respect `CpuBudget` and existing concurrency knobs
    (`mafft_threads_per_job`, `mafft_max_jobs`, etc.) as much as possible.

- **Safe experimentation**
  - MAFFT remains the default; SPOA can be enabled via CLI/config for benchmark
    runs (e.g., head200/head1000/5k/full) without changing published defaults.
  - SPOA failures should **fall back** to MAFFT (with a warning) so we don’t
    hard‑fail runs while we are still stabilising behaviour.

---

## Current Alignment Flow (MAFFT)

This is a quick map of the current MAFFT path, so you know where to hook SPOA
in. Line numbers are approximate and may drift; use `rg` to locate symbols.

- `src/mafft.rs`
  - Contains helper(s) to:
    - Write panel sequences to a temporary FASTA.
    - Spawn `mafft` via `std::process::Command` with `--quiet`, optional
      `--thread` flags, and “fast” options (FFT‑NS‑1) when configured.
    - Parse the resulting MSA into an internal alignment struct that records:
      - aligned sequences (per gene),
      - per‑column conservation / gap info, or at least enough to derive
        conserved fraction, gap runs and termini concordance.
  - Public functions (names may differ slightly, check the file):
    - `run_mafft_alignment(...)` or similar: a synchronous helper invoked by
      ECS jobs.
    - A public `AlignmentStats` struct (or equivalent) used by scoring.

- `src/ecs.rs`
  - The MAFFT pipeline is wired here. Look for something like:
    - `MafftPipelineConfig { ... }` struct,
    - `run_alignment_pipeline(app, ecs_config, mafft_config, ...)`,
    - Systems that:
      1. Seed alignment jobs for genes with an eligible consensus panel.
      2. Launch MAFFT subprocesses in bounded parallelism (respecting
         `mafft_max_jobs` and `CpuBudget`).
      3. Collect MSA results and attach alignment components to gene entities.
  - This pipeline is already integrated with `CpuBudget` and the global ECS
    scheduler, so we want to **reuse this structure**, not build a parallel
    SPOA pipeline.

- `src/main.rs`
  - CLI / config definition for MAFFT:
    - Flags like `--mafft-bin`, `--mafft-fast`, `--mafft-threads-per-job`,
      `--mafft-max-jobs` are defined in the top‑level CLI struct.
    - Config overrides live in `ConfigFile` structs (e.g., `[mafft]` section)
      and merged with CLI in `main()`.
    - The final `EcsConfig` and `MafftPipelineConfig` are constructed and passed
      into `run_alignment_pipeline`.

- `src/scoring.rs` + `src/length.rs` + `src/structvar.rs`
  - These modules consume the alignment struct (via `GeneMetrics` / alignment
    components) to compute:
    - `start_concordance`, `end_concordance`,
    - `missing_exon_run`, `retained_intron_run`,
    - gap statistics, conserved fraction, etc.
  - As long as SPOA produces the same alignment representation, these modules
    do not need to change in any fundamental way.

---

## Target Design: Aligner Enum + SPOA Backend

### 1. Add an `Aligner` enum

Location: `src/mafft.rs` or a new `src/align.rs` if you prefer a neutral name.

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignerBackend {
    Mafft,
    Spoa,
}
```

Expose this enum via a small configuration struct:

```rust
#[derive(Debug, Clone)]
pub struct AlignerConfig {
    pub backend: AlignerBackend,
    pub mafft_bin: String,
    pub mafft_fast: bool,
    pub mafft_threads_per_job: usize,
    pub mafft_max_jobs: usize,
    // Optional: SPOA‑specific knobs later (gap penalties, matrix).
}
```

### 2. Wire CLI + config

Location: `src/main.rs` in the CLI struct and config merge logic.

**CLI changes**

Add an enum to the CLI args:

```rust
#[derive(Debug, Clone, Copy, clap::ValueEnum)]
pub enum AlignerCliBackend {
    #[clap(name = "mafft")]
    Mafft,
    #[clap(name = "spoa")]
    Spoa,
}

impl From<AlignerCliBackend> for mafft::AlignerBackend {
    fn from(v: AlignerCliBackend) -> Self {
        match v {
            AlignerCliBackend::Mafft => mafft::AlignerBackend::Mafft,
            AlignerCliBackend::Spoa => mafft::AlignerBackend::Spoa,
        }
    }
}
```

Then extend the top‑level `Args` struct:

```rust
#[derive(Parser, Debug)]
pub struct Args {
    // existing fields …

    /// Alignment backend (mafft or spoa)
    #[clap(long, value_enum, default_value = "mafft")]
    pub aligner: AlignerCliBackend,
}
```

**Config changes**

In the TOML config structs (still in `src/main.rs`), add an optional field:

```rust
#[derive(Debug, Deserialize)]
pub struct MafftConfig {
    pub bin: Option<String>,
    pub fast: Option<bool>,
    pub threads_per_job: Option<usize>,
    pub max_jobs: Option<usize>,
    pub backend: Option<String>, // "mafft" or "spoa"
}
```

When merging CLI + file config into `AlignerConfig`:

```rust
let backend = if let Some(b) = file_cfg.mafft.as_ref().and_then(|m| m.backend.as_deref()) {
    match b.to_ascii_lowercase().as_str() {
        "spoa" => mafft::AlignerBackend::Spoa,
        _ => mafft::AlignerBackend::Mafft,
    }
} else {
    args.aligner.into()
};

let aligner_cfg = mafft::AlignerConfig {
    backend,
    mafft_bin: args.mafft_bin.clone(),
    mafft_fast: args.mafft_fast
        || file_cfg.mafft.as_ref().and_then(|m| m.fast).unwrap_or(false),
    mafft_threads_per_job: args.mafft_threads_per_job
        .or_else(|| file_cfg.mafft.as_ref().and_then(|m| m.threads_per_job))
        .unwrap_or(2),
    mafft_max_jobs: args.mafft_max_jobs
        .or_else(|| file_cfg.mafft.as_ref().and_then(|m| m.max_jobs))
        .unwrap_or(6),
};
```

Finally, pass `aligner_cfg` into `EcsConfig` / `MafftPipelineConfig` instead of
passing individual MAFFT fields.

### 3. Implement SPOA backend

Location: `src/mafft.rs` (or new `src/align.rs`).

#### 3.1 Define a shared alignment struct

You likely already have something like this; if not, define it and migrate
MAFFT to use it explicitly:

```rust
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    pub gene_id: String,
    pub seq_ids: Vec<String>,          // panel sequence IDs, in order
    pub aligned_seqs: Vec<String>,     // same length; gapped
    // Optionally: precomputed stats if the file already does this.
}
```

MAFFT path should construct this struct from the FASTA + output MSA (no change
in observed behaviour).

#### 3.2 SPOA helper

Add a function using the `spoa` crate:

```rust
use spoa::{AlignmentEngine, AlignmentMode, Graph, Scoring};

pub fn run_spoa_alignment(
    gene_id: &str,
    ids_and_seqs: &[(String, String)],
) -> Result<AlignmentResult, String> {
    if ids_and_seqs.len() < 2 {
        // nothing to align; mimic MAFFT behaviour
        return Ok(AlignmentResult {
            gene_id: gene_id.to_string(),
            seq_ids: ids_and_seqs.iter().map(|(id, _)| id.clone()).collect(),
            aligned_seqs: ids_and_seqs.iter().map(|(_, s)| s.clone()).collect(),
        });
    }

    // Simple protein scoring scheme; can be tuned later.
    let scoring = Scoring::new(2, -1, -2, -1); // match, mismatch, gap_open, gap_extend
    let mut engine = AlignmentEngine::new(AlignmentMode::Global, scoring)
        .map_err(|e| format!("spoa engine init failed: {e}"))?;
    let mut graph = Graph::new();

    for (_, seq) in ids_and_seqs {
        let s = seq.as_bytes();
        let alignment = engine
            .align(s, &graph)
            .map_err(|e| format!("spoa align failed: {e}"))?;
        graph
            .add_sequence(&alignment, s)
            .map_err(|e| format!("spoa add_sequence failed: {e}"))?;
    }

    let msa = graph
        .generate_multiple_sequence_alignment()
        .map_err(|e| format!("spoa MSA generation failed: {e}"))?;

    let aligned_seqs: Vec<String> = msa
        .iter()
        .map(|row| String::from_utf8_lossy(row).into_owned())
        .collect();

    if aligned_seqs.len() != ids_and_seqs.len() {
        return Err(format!(
            "spoa MSA row count mismatch: expected {} got {}",
            ids_and_seqs.len(),
            aligned_seqs.len()
        ));
    }

    Ok(AlignmentResult {
        gene_id: gene_id.to_string(),
        seq_ids: ids_and_seqs.iter().map(|(id, _)| id.clone()).collect(),
        aligned_seqs,
    })
}
```

**Notes / pitfalls**

- This uses a global POA alignment; for our panel sizes that is fine.
- Scoring parameters are placeholder; we can tune later.
- We don’t try to emulate MAFFT’s exact gap pattern; minor shifts are acceptable
  per your guidance, but they will slightly perturb termini metrics.

#### 3.3 Unified align function

Add a single entry point used by ECS:

```rust
pub fn run_alignment_for_panel(
    cfg: &AlignerConfig,
    gene_id: &str,
    ids_and_seqs: &[(String, String)],
) -> Result<AlignmentResult, String> {
    match cfg.backend {
        AlignerBackend::Mafft => run_mafft_alignment(cfg, gene_id, ids_and_seqs),
        AlignerBackend::Spoa => {
            match run_spoa_alignment(gene_id, ids_and_seqs) {
                Ok(a) => Ok(a),
                Err(e) => {
                    log::warn!("SPOA failed for {}: {}; falling back to MAFFT", gene_id, e);
                    run_mafft_alignment(cfg, gene_id, ids_and_seqs)
                }
            }
        }
    }
}
```

Here `run_mafft_alignment` is your existing helper rewritten to accept the
same `(gene_id, ids_and_seqs)` signature.

### 4. ECS integration

Location: `src/ecs.rs`.

Look for the MAFFT pipeline, something like:

- `MafftPipelineConfig` struct;
- A system that, given `MafftPipelineConfig`, pulls sequences for each gene and
  launches asynchronous MAFFT jobs.

Changes:

1. Replace `MafftPipelineConfig` with (or extend it to hold) `AlignerConfig`.
2. Where the job system currently calls `mafft::run_mafft_alignment(...)`, call
   `mafft::run_alignment_for_panel(&aligner_cfg, gene_id, &ids_and_seqs)`.
3. Keep job scheduling identical (same `CpuBudget` plumbing, job throttling,
   and error propagation). The difference is only inside the job body.

No changes are needed to `run_heavy_pipelines` or `CpuBudget` so long as SPOA
jobs are treated as the same “heavy” category as MAFFT.

### 5. Downstream scoring

Because `AlignmentResult` is identical for both backends, scoring modules
should not need structural changes. That said, we should be aware of:

- **Termini concordance:** Slight alignment shifts may change the exact start
  position distribution; we should run a small comparison:
  - For head200/head1000, compute differences in `start_concordance`,
    `start_class`, `end_concordance`, `end_class` between MAFFT and SPOA.
  - If the majority of changes are small (e.g., <0.05 concordance shift) and
    classifications rarely flip between “LikelyComplete” vs truncation, we can
    accept that as a trade‑off.

- **Gap‑based warnings:** `missing_exon_run`, `retained_intron_run`, and
  overall gap_run statistics will change slightly. Worth checking the counts of
  each warning class before/after.

---

## Testing & Benchmarking Plan

### 1. Unit tests

Add targeted unit tests in `src/mafft.rs` (or a new `tests/align_spoa.rs`):

- **Round‑trip:** For a small hand‑crafted panel (3–4 sequences), verify SPOA
  returns:
  - the same number of rows as input sequences,
  - all rows have identical length,
  - no unexpected characters (only amino acids plus gap `-`).

- **Fallback path:** Force SPOA to error (e.g., by passing an empty panel) and
  assert that we fall back to `run_mafft_alignment` (can mock via a small helper
  that returns a known alignment).

### 2. CLI smoke test

Extend `tests/analyze_smoke.rs` to add a new test case that uses SPOA:

```rust
#[test]
fn analyze_smoke_spoa() {
    let tmp = TempDir::new().unwrap();
    let out_dir = tmp.path().join("out_spoa");
    let status = Command::cargo_bin("AnnoQC")
        .unwrap()
        .args([
            "--config", "config.weights_demo.toml",
            "analyze",
            "--fasta", "a9_head200.faa",
            "--db", "uniprot_sprot.dmnd",
            "--out", out_dir.to_str().unwrap(),
            "--threads", "4",
            "--aligner", "spoa",
            "--mafft-bin", "mafft", // still required for fallback
        ])
        .status()
        .unwrap();
    assert!(status.success());
}
```

### 3. Runtime comparison

For benchmarking, run (outside the test suite):

- MAFFT baseline:

  ```bash
  RUST_LOG=info pixi run cargo run -- \
    --config config.weights_demo.toml \
    analyze --fasta a9_head200.faa --db uniprot_sprot.dmnd \
    --out results_a9_head200_mafft \
    --threads 16 --mafft-bin mafft --aligner mafft
  ```

- SPOA variant:

  ```bash
  RUST_LOG=info pixi run cargo run -- \
    --config config.weights_demo.toml \
    analyze --fasta a9_head200.faa --db uniprot_sprot.dmnd \
    --out results_a9_head200_spoa \
    --threads 16 --mafft-bin mafft --aligner spoa
  ```

Then compare:

- `run_metrics.json` → `mafft_alignments` time vs SPOA, total ECS time.
- `qc_summary.csv` → classify differences in `final_score`, termini classes,
  and warning columns.

If the SPOA path is faster and classification changes are acceptable, we can
consider making it the recommended backend in docs (still not the default
until we are comfortable).

---

## Pitfalls & Open Questions

- **Scoring differences:** SPOA’s scoring model is not MAFFT’s; we currently
  choose a simple match/mismatch/gap scheme. This may slightly bias alignments
  toward different gap placements. If this becomes a problem, we can:
  - Tune SPOA scoring (or use BLOSUM‑like matrix support if available), and/or
  - Add a “strict” mode that keeps MAFFT for panels where termini metrics are
    especially sensitive.

- **Very small panels:** For panels of size 1–2, alignment is effectively
  trivial; SPOA is overkill. We already short‑circuit in
  `run_spoa_alignment` to return identity alignments in that case.

- **CPU budget accounting:** SPOA runs fully in process, so we should be
  careful not to spawn too many jobs per frame. Reusing the existing MAFFT
  job throttling is the right approach; avoid spinning up a separate SPOA
  pipeline that ignores `CpuBudget`.

- **Dependency footprint:** `spoa` pulls in `cxx` + `spoa-sys` and system
  toolchain requirements. We have already fetched and built it successfully
  using `CARGO_HOME=.cargo`. On other environments, ensure a C++ toolchain is
  present.

- **Future nucleotide mode:** Today we align proteins only. If we later support
  nucleotide‑level alignments (nt/ORF mode from `GV_PARITY.md`), we may need
  separate scoring/parameter sets for SPOA in nucleotide space.

---

## Summary

To integrate SPOA safely:

1. **Add `AlignerBackend` + `AlignerConfig` in `mafft.rs` (or `align.rs`).**
2. **Introduce a CLI enum `AlignerCliBackend` and config field** to choose
   between `mafft` and `spoa`, defaulting to `mafft`.
3. **Implement `run_spoa_alignment`** that consumes `(gene_id, ids_and_seqs)`
   and returns the same `AlignmentResult` struct as MAFFT.
4. **Provide `run_alignment_for_panel`** that dispatches to MAFFT or SPOA and
   falls back to MAFFT on SPOA errors.
5. **Swap ECS alignment jobs** to call `run_alignment_for_panel` instead of the
   MAFFT‑only helper.
6. **Add unit tests + a SPOA CLI smoke test** and run head200/head1000
   benchmarks to quantify runtime and scoring differences.

This keeps the system coherent with the existing ECS design (one alignment
pipeline, multiple backends) and gives us a clear path to evaluate whether
SPOA genuinely improves throughput in the presence of heavy DIAMOND/HMMER work.

