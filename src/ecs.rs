use std::collections::{BTreeMap, HashMap, VecDeque};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;

use bevy_app::{App, First, Last, PostUpdate, PreUpdate, Startup, Update};
use bevy_ecs::prelude::*;
use bevy_tasks::{AsyncComputeTaskPool, Task, TaskPoolBuilder};
use futures_lite::future::{block_on, poll_once};
use needletail::parse_fastx_file;
use polars::error::PolarsError;
use polars::export::arrow::datatypes::{ArrowDataType, ArrowSchema, PhysicalType};
use polars::export::arrow::types::PrimitiveType;
use polars::prelude::*;
use polars_parquet::arrow::write::{
    transverse, CompressionOptions, Encoding, FileWriter, RowGroupIterator, Version, WriteOptions,
};
use polars_parquet::parquet::metadata::KeyValue;
use std::time::{Duration, Instant};

use crate::hmmer::{run_hmmscan, HmmscanSummary};
use crate::mafft::{run_alignment_for_panel, AlignerConfig, AlignmentMetrics};
use crate::preflight;
use crate::render_gene_record;
use crate::{
    taxonomy, Checksums, RenderContext, RenderedRecord, ScoreCard, OUTPUT_SCHEMA_VERSION,
    TOOL_NAME, TOOL_VERSION,
};

const CPU_LABEL_ALIGN: &str = "align";
// CPU budget label for HMMER jobs.
const CPU_LABEL_HMMER: &str = "hmmer";
const ALIGNMENT_LAUNCH_CAP: usize = 5;
const HMMER_LAUNCH_CAP: usize = 5;
const PARQUET_ROW_GROUP_SIZE: usize = 5000;

#[derive(Resource)]
struct CpuBudget {
    capacity: usize,
    available: usize,
    reserved: usize,
    usage_threads: HashMap<&'static str, usize>,
    inflight_jobs: HashMap<&'static str, usize>,
    last_log: Instant,
    last_reported_used: usize,
    log_interval: Duration,
    log_json: bool,
    dirty: bool,
}

impl CpuBudget {
    fn new(total_threads: usize, reserve: usize, log_json: bool) -> Self {
        let capacity = total_threads.saturating_sub(reserve).max(1);
        Self {
            capacity,
            available: capacity,
            reserved: reserve,
            usage_threads: HashMap::new(),
            inflight_jobs: HashMap::new(),
            last_log: Instant::now(),
            last_reported_used: 0,
            log_interval: Duration::from_secs(60),
            log_json,
            dirty: true,
        }
    }

    fn try_acquire(&mut self, label: &'static str, threads: usize) -> bool {
        if threads == 0 {
            *self.inflight_jobs.entry(label).or_insert(0) += 1;
            self.dirty = true;
            return true;
        }
        if threads > self.available {
            return false;
        }
        self.available -= threads;
        *self.usage_threads.entry(label).or_insert(0) += threads;
        *self.inflight_jobs.entry(label).or_insert(0) += 1;
        self.dirty = true;
        true
    }

    fn release(&mut self, label: &'static str, threads: usize) {
        if threads > 0 {
            self.available = (self.available + threads).min(self.capacity);
            if let Some(entry) = self.usage_threads.get_mut(label) {
                *entry = entry.saturating_sub(threads);
            }
        }
        if let Some(entry) = self.inflight_jobs.get_mut(label) {
            *entry = entry.saturating_sub(1);
        }
        self.dirty = true;
    }

    fn snapshot_due(&self) -> bool {
        let used = self.capacity.saturating_sub(self.available);
        self.dirty
            || used != self.last_reported_used
            || Instant::now().duration_since(self.last_log) >= self.log_interval
    }

    fn record_snapshot(&mut self) {
        self.last_log = Instant::now();
        self.last_reported_used = self.capacity.saturating_sub(self.available);
        self.dirty = false;
    }

    fn log_snapshot(&self) {
        let used = self.capacity.saturating_sub(self.available);
        let mut labels: Vec<_> = self.usage_threads.keys().copied().collect();
        labels.sort_unstable();
        let mut per_label = Vec::new();
        for label in labels {
            let threads = *self.usage_threads.get(label).unwrap_or(&0);
            let jobs = *self.inflight_jobs.get(label).unwrap_or(&0);
            per_label.push(format!("{}:threads={} jobs={}", label, threads, jobs));
        }
        if self.log_json {
            log::info!(
                "{}",
                serde_json::json!({
                    "event": "cpu_budget",
                    "capacity": self.capacity,
                    "reserved": self.reserved,
                    "used": used,
                    "available": self.available,
                    "labels": per_label,
                })
            );
        } else {
            log::info!(
                "cpu budget: used={}/{} reserved={} available={} [{}]",
                used,
                self.capacity,
                self.reserved,
                self.available,
                per_label.join(" ")
            );
        }
    }
}

fn cpu_budget_log_system(mut budget: ResMut<CpuBudget>) {
    if budget.snapshot_due() {
        budget.record_snapshot();
        budget.log_snapshot();
    }
}

#[derive(Resource, Clone)]
pub struct EcsConfig {
    pub fasta_path: String,
    pub diamond_tsv: String,
    pub threads: usize,
    pub log_json: bool,
}

#[derive(Copy, Clone)]
pub struct RenderOutputConfig {
    pub jsonl: bool,
    pub csv: bool,
    pub parquet: bool,
    pub resume: bool,
}

#[derive(Debug, Clone)]
pub struct GeneMetrics {
    pub gene_id: String,
    // Sequence length from the input FASTA.
    #[allow(dead_code)]
    pub length: usize,
    pub hits: usize,
}

#[derive(Resource, Default)]
pub struct CompletionState {
    pub finished: bool,
}

#[derive(Clone)]
pub struct AlignmentPipelineConfig {
    pub aligner: AlignerConfig,
    pub jobs: Vec<(String, Vec<String>)>,
    pub query_map: HashMap<String, Vec<u8>>,
    pub ref_map: HashMap<String, Vec<u8>>,
}

#[derive(Clone)]
pub struct HmmerPipelineConfig {
    pub hmmscan_bin: String,
    pub db_path: String,
    pub items: Vec<(String, Vec<u8>)>,
    pub threads_per_job: usize,
    pub max_jobs: usize,
    pub top_n: usize,
    pub max_ievalue: Option<f64>,
}

pub struct HeavyPipelineConfig {
    pub cpu_threads: usize,
    pub reserve_threads: usize,
    pub log_json: bool,
    pub alignment: Option<AlignmentPipelineConfig>,
    pub hmmer: Option<HmmerPipelineConfig>,
}

pub struct HeavyPipelineResults {
    pub alignment_map: HashMap<String, AlignmentMetrics>,
    pub hmmer_map: HashMap<String, HmmscanSummary>,
    pub alignment_secs: HashMap<String, f64>,
    pub hmmer_secs: HashMap<String, f64>,
}

#[derive(Resource, Default)]
struct WorkQueue {
    queue: Vec<(String, usize)>,
    hit_counts: HashMap<String, usize>,
}

#[derive(Resource, Default)]
struct Results {
    records: Vec<GeneMetrics>,
}

#[derive(Resource, Default)]
struct InFlight {
    tasks: Vec<Task<GeneMetrics>>,
    max_in_flight: usize,
}

#[derive(Resource, Default)]
struct LengthMap(pub HashMap<String, usize>);

#[derive(Resource)]
struct Progress {
    start: Instant,
    last_log: Instant,
    processed: usize,
    total: usize,
    log_json: bool,
}

fn intake_system(mut commands: Commands, cfg: Res<EcsConfig>) {
    // Parse FASTA ids and lengths
    let mut ids: Vec<(String, usize)> = Vec::new();
    let mut lengths = HashMap::<String, usize>::new();
    let mut reader = parse_fastx_file(&cfg.fasta_path).expect("open fasta");
    while let Some(record) = reader.next() {
        let rec = record.expect("valid record");
        let id = String::from_utf8_lossy(rec.id()).to_string();
        let gene = id.split_whitespace().next().unwrap_or("").to_string();
        if gene.is_empty() {
            continue;
        }
        // predictive cost: length + 100*hit_count (filled later)
        ids.push((gene.clone(), 0));
        lengths.insert(gene, rec.seq().len());
    }

    // Parse DIAMOND tsv -> qseqid -> count
    let mut hit_counts: HashMap<String, usize> = HashMap::new();
    if !cfg.diamond_tsv.is_empty() && std::path::Path::new(&cfg.diamond_tsv).exists() {
        if let Ok(text) = std::fs::read_to_string(&cfg.diamond_tsv) {
            for line in text.lines() {
                if line.is_empty() {
                    continue;
                }
                if let Some((qseqid, _rest)) = line.split_once('\t') {
                    *hit_counts.entry(qseqid.to_string()).or_insert(0) += 1;
                }
            }
        }
    }

    // Build predictive costs and sort heavy-first
    for tup in ids.iter_mut() {
        let (ref gid, ref mut cost) = tup;
        let len = *lengths.get(gid).unwrap_or(&0);
        let hits = *hit_counts.get(gid).unwrap_or(&0);
        *cost = len + hits * 100;
    }
    ids.sort_by(|a, b| a.1.cmp(&b.1)); // ascending
    let total = ids.len();
    commands.insert_resource(WorkQueue {
        queue: ids,
        hit_counts,
    });
    commands.insert_resource(Results::default());
    commands.insert_resource(InFlight {
        tasks: Vec::with_capacity(cfg.threads.max(1)),
        max_in_flight: cfg.threads.max(1),
    });
    commands.insert_resource(LengthMap(lengths));
    let now = Instant::now();
    commands.insert_resource(Progress {
        start: now,
        last_log: now,
        processed: 0,
        total,
        log_json: cfg.log_json,
    });
}

fn schedule_system(
    mut queue: ResMut<WorkQueue>,
    lengths: Res<LengthMap>,
    mut inflight: ResMut<InFlight>,
) {
    let pool = AsyncComputeTaskPool::get();
    while inflight.tasks.len() < inflight.max_in_flight {
        let Some((gene, _cost)) = queue.queue.pop() else {
            break;
        };
        let length = *lengths.0.get(&gene).unwrap_or(&0);
        let hits = *queue.hit_counts.get(&gene).unwrap_or(&0);
        let gene_cloned = gene.clone();
        let task = pool.spawn(async move {
            // Placeholder for heavier per-gene work
            GeneMetrics {
                gene_id: gene_cloned,
                length,
                hits,
            }
        });
        inflight.tasks.push(task);
    }
}

fn collect_system(
    mut inflight: ResMut<InFlight>,
    mut results: ResMut<Results>,
    mut prog: ResMut<Progress>,
    queue: Res<WorkQueue>,
) {
    let mut i = 0;
    while i < inflight.tasks.len() {
        let done = block_on(poll_once(&mut inflight.tasks[i]));
        if let Some(metrics) = done {
            results.records.push(metrics);
            drop(inflight.tasks.swap_remove(i));
            prog.processed += 1;
        } else {
            i += 1;
        }
    }
    let now = Instant::now();
    if now.duration_since(prog.last_log) > Duration::from_secs(1)
        || prog.processed.is_multiple_of(100)
    {
        let dt = now.duration_since(prog.start).as_secs_f64().max(1e-6);
        let rate = (prog.processed as f64) / dt;
        if prog.log_json {
            log::info!(
                "{}",
                serde_json::json!({
                    "event":"progress",
                    "processed":prog.processed,
                    "total":prog.total,
                    "rate":format!("{:.2}",rate),
                    "inflight": inflight.tasks.len(),
                    "queue_remaining": queue.queue.len()
                })
            );
        } else {
            log::info!(
                "progress: {}/{} ({:.2} genes/s) inflight={} queue={}",
                prog.processed,
                prog.total,
                rate,
                inflight.tasks.len(),
                queue.queue.len()
            );
        }
        prog.last_log = now;
    }
}

fn finish_system(
    queue: Res<WorkQueue>,
    inflight: Res<InFlight>,
    mut state: ResMut<CompletionState>,
) {
    if queue.queue.is_empty() && inflight.tasks.is_empty() {
        state.finished = true;
    }
}

pub fn run_scheduler(cfg: EcsConfig) -> Vec<GeneMetrics> {
    // Initialize async task pool before scheduling systems use it
    AsyncComputeTaskPool::get_or_init(|| TaskPoolBuilder::new().num_threads(cfg.threads).build());
    let mut app = App::new();
    app.insert_resource(cfg)
        .insert_resource(CompletionState::default())
        .add_systems(Startup, intake_system)
        .add_systems(Update, (schedule_system, collect_system, finish_system));

    loop {
        app.update();
        let state = app.world().get_resource::<CompletionState>().unwrap();
        if state.finished {
            break;
        }
    }
    let results = app.world().get_resource::<Results>().unwrap();
    results.records.clone()
}

// --------- Alignment pipeline ---------

#[derive(Component)]
struct AlignmentJob {
    gene_id: String,
    panel_ids: Vec<String>,
}

#[derive(Resource)]
struct AlignmentSeedJobs(Vec<(String, Vec<String>)>);

#[derive(Component)]
struct AlignmentHandle {
    task: Task<Result<AlignmentMetrics, String>>,
    threads_used: usize,
    started_at: Instant,
}

#[derive(Resource)]
struct AlignmentConfig {
    aligner: AlignerConfig,
}

#[derive(Resource)]
struct AlignmentStores {
    query_map: Arc<HashMap<String, Vec<u8>>>,
    ref_map: Arc<HashMap<String, Vec<u8>>>,
}

#[derive(Resource, Default)]
struct AlignmentResultsRes {
    map: HashMap<String, AlignmentMetrics>,
}

#[derive(Resource, Default)]
struct AlignmentTimingRes {
    map: HashMap<String, f64>,
}

#[derive(Resource, Default)]
struct AlignmentInflight {
    current: usize,
}

#[derive(Resource)]
struct AlignmentProgress {
    total: usize,
    completed: usize,
}

#[derive(Resource)]
#[allow(dead_code)]
struct AlignmentDiag {
    total: usize,
    skipped_no_panel: usize,
    skipped_no_query: usize,
    skipped_no_refs: usize,
    scheduled: usize,
    last_log: std::time::Instant,
    next_pct: usize,
}

impl AlignmentDiag {
    fn new(total: usize) -> Self {
        Self {
            total,
            skipped_no_panel: 0,
            skipped_no_query: 0,
            skipped_no_refs: 0,
            scheduled: 0,
            last_log: std::time::Instant::now(),
            next_pct: 5,
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn alignment_schedule_system(
    mut commands: Commands,
    config: Res<AlignmentConfig>,
    stores: Res<AlignmentStores>,
    mut inflight: ResMut<AlignmentInflight>,
    mut progress: ResMut<AlignmentProgress>,
    mut diag: ResMut<AlignmentDiag>,
    mut budget: ResMut<CpuBudget>,
    hmmer_presence: Res<HmmerPresence>,
    mut query: Query<(Entity, &AlignmentJob), Without<AlignmentHandle>>,
) {
    if progress.total == 0 {
        return;
    }
    let max_jobs = config.aligner.mafft_max_jobs.max(1);
    if inflight.current >= max_jobs {
        return;
    }
    let pool = AsyncComputeTaskPool::get();
    let threads_needed = config.aligner.mafft_threads_per_job.max(1);
    let reserve_for_hmmer = if hmmer_presence.active && hmmer_presence.reserve_threads > 0 {
        hmmer_presence
            .reserve_threads
            .min(budget.capacity.saturating_sub(1))
    } else {
        0
    };
    let mut launched_this_frame = 0usize;
    for (entity, job) in query.iter_mut() {
        if inflight.current >= max_jobs {
            break;
        }
        if launched_this_frame >= ALIGNMENT_LAUNCH_CAP {
            break;
        }
        if job.panel_ids.is_empty() {
            progress.completed += 1;
            if let Ok(mut ecmd) = commands.get_entity(entity) {
                ecmd.insert(AlignmentDone);
                ecmd.despawn();
            }
            diag.skipped_no_panel += 1;
            continue;
        }
        let Some(qseq) = stores.query_map.get(&job.gene_id).cloned() else {
            log::warn!("alignment skipped {}, query sequence missing", job.gene_id);
            progress.completed += 1;
            if let Ok(mut ecmd) = commands.get_entity(entity) {
                ecmd.insert(AlignmentDone);
                ecmd.despawn();
            }
            diag.skipped_no_query += 1;
            continue;
        };
        let mut ids_and_seqs: Vec<(String, String)> = Vec::with_capacity(job.panel_ids.len() + 1);
        ids_and_seqs.push((
            job.gene_id.clone(),
            String::from_utf8_lossy(&qseq).into_owned(),
        ));
        for id in job.panel_ids.iter() {
            let key = taxonomy::canonical_accession(id);
            if let Some(seq) = stores.ref_map.get(&key) {
                ids_and_seqs.push((key.clone(), String::from_utf8_lossy(seq).into_owned()));
            }
        }
        if ids_and_seqs.len() <= 1 {
            progress.completed += 1;
            if let Ok(mut ecmd) = commands.get_entity(entity) {
                ecmd.insert(AlignmentDone);
                ecmd.despawn();
            }
            diag.skipped_no_refs += 1;
            continue;
        }
        if reserve_for_hmmer > 0 && budget.available <= reserve_for_hmmer {
            break;
        }
        if !budget.try_acquire(CPU_LABEL_ALIGN, threads_needed) {
            break;
        }
        let gene_id = job.gene_id.clone();
        let cfg = config.aligner.clone();
        let task =
            pool.spawn(async move { run_alignment_for_panel(&cfg, &gene_id, &ids_and_seqs) });
        if let Ok(mut ecmd) = commands.get_entity(entity) {
            ecmd.insert(AlignmentHandle {
                task,
                threads_used: threads_needed,
                started_at: Instant::now(),
            });
        }
        inflight.current += 1;
        diag.scheduled += 1;
        launched_this_frame += 1;
    }
}

#[allow(clippy::too_many_arguments)]
fn alignment_collect_system(
    mut commands: Commands,
    mut inflight: ResMut<AlignmentInflight>,
    mut results: ResMut<AlignmentResultsRes>,
    mut timings: ResMut<AlignmentTimingRes>,
    mut progress: ResMut<AlignmentProgress>,
    mut diag: ResMut<AlignmentDiag>,
    mut budget: ResMut<CpuBudget>,
    mut query: Query<(Entity, &AlignmentJob, &mut AlignmentHandle)>,
) {
    for (entity, job, mut handle) in query.iter_mut() {
        if let Some(res) = block_on(poll_once(&mut handle.task)) {
            let elapsed = handle.started_at.elapsed().as_secs_f64();
            inflight.current = inflight.current.saturating_sub(1);
            budget.release(CPU_LABEL_ALIGN, handle.threads_used);
            match res {
                Ok(metrics) => {
                    results.map.insert(job.gene_id.clone(), metrics);
                }
                Err(err) => {
                    log::warn!("alignment failed for {}: {}", job.gene_id, err);
                }
            }
            timings.map.insert(job.gene_id.clone(), elapsed);
            progress.completed += 1;
            if let Ok(mut ecmd) = commands.get_entity(entity) {
                ecmd.remove::<AlignmentHandle>();
                ecmd.insert(AlignmentDone);
                ecmd.despawn();
            }
        }
    }
    // Throttled progress logging: on each 5% completion or every 60s
    let now = std::time::Instant::now();
    let pct = if progress.total > 0 {
        (progress.completed * 100) / progress.total
    } else {
        100
    };
    let should_log_pct = pct >= diag.next_pct && diag.next_pct <= 100;
    let should_log_time = now.duration_since(diag.last_log).as_secs() >= 60;
    if should_log_pct || should_log_time {
        log::info!(
            "alignment progress: completed={}/{} ({}%), inflight={}, scheduled={}, skipped_no_panel={}, skipped_no_query={}, skipped_no_refs={}",
            progress.completed,
            progress.total,
            pct,
            inflight.current,
            diag.scheduled,
            diag.skipped_no_panel,
            diag.skipped_no_query,
            diag.skipped_no_refs
        );
        diag.last_log = now;
        if should_log_pct {
            diag.next_pct = (pct / 5 + 1) * 5; // next 5% bucket
        }
    }
}

#[derive(Component)]
struct AlignmentDone;
fn alignment_seed_system(mut commands: Commands, mut seeds: ResMut<AlignmentSeedJobs>) {
    for (gene_id, panel_ids) in seeds.0.drain(..) {
        commands.spawn(AlignmentJob { gene_id, panel_ids });
    }
}

// --------- HMMER pipeline ---------

#[derive(Component)]
struct HmmerJob {
    gene_id: String,
    seq: Vec<u8>,
}

#[derive(Component)]
struct HmmerHandle {
    task: Task<Result<HmmscanSummary, String>>,
    threads_used: usize,
    started_at: Instant,
}

#[derive(Component)]
struct HmmerDone;

#[derive(Resource)]
struct HmmerConfig {
    bin: String,
    db_path: String,
    max_jobs: usize,
    top_n: usize,
    max_ievalue: Option<f64>,
    threads_per_job: usize,
}

#[derive(Resource, Default)]
struct HmmerResultsRes {
    map: HashMap<String, HmmscanSummary>,
}

#[derive(Resource, Default)]
struct HmmerTimingRes {
    map: HashMap<String, f64>,
}

#[derive(Resource, Default)]
struct HmmerInflight {
    current: usize,
}

#[derive(Resource)]
struct HmmerProgress {
    total: usize,
    completed: usize,
}

#[derive(Resource)]
struct HmmerSeedJobs(Vec<(String, Vec<u8>)>);

#[derive(Resource, Default, Clone, Copy)]
struct HmmerPresence {
    active: bool,
    reserve_threads: usize,
}

pub fn run_heavy_pipelines(cfg: HeavyPipelineConfig) -> HeavyPipelineResults {
    let HeavyPipelineConfig {
        cpu_threads,
        reserve_threads,
        log_json,
        mut alignment,
        mut hmmer,
    } = cfg;
    AsyncComputeTaskPool::get_or_init(|| {
        TaskPoolBuilder::new()
            .num_threads(cpu_threads.max(1))
            .build()
    });

    let alignment_total = alignment.as_ref().map(|c| c.jobs.len()).unwrap_or(0);
    let hmmer_total = hmmer.as_ref().map(|c| c.items.len()).unwrap_or(0);
    let hmmer_reserve_threads = hmmer
        .as_ref()
        .map(|c| c.threads_per_job.max(1))
        .unwrap_or(0);

    let mut app = App::new();
    app.insert_resource(CompletionState::default())
        .insert_resource(CpuBudget::new(
            cpu_threads.max(1),
            reserve_threads,
            log_json,
        ))
        .insert_resource(AlignmentResultsRes::default())
        .insert_resource(AlignmentTimingRes::default())
        .insert_resource(AlignmentInflight::default())
        .insert_resource(AlignmentProgress {
            total: alignment_total,
            completed: 0,
        })
        .insert_resource(AlignmentDiag::new(alignment_total))
        .insert_resource(HmmerResultsRes::default())
        .insert_resource(HmmerTimingRes::default())
        .insert_resource(HmmerInflight::default())
        .insert_resource(HmmerProgress {
            total: hmmer_total,
            completed: 0,
        })
        .insert_resource(HmmerSeedJobs(Vec::new()))
        .insert_resource(HmmerPresence {
            active: hmmer_total > 0,
            reserve_threads: hmmer_reserve_threads,
        });

    if let Some(mcfg) = alignment.take() {
        app.insert_resource(AlignmentConfig {
            aligner: mcfg.aligner.clone(),
        })
        .insert_resource(AlignmentStores {
            query_map: Arc::new(mcfg.query_map),
            ref_map: Arc::new(mcfg.ref_map),
        })
        .insert_resource(AlignmentSeedJobs(mcfg.jobs))
        .add_systems(Startup, alignment_seed_system)
        .add_systems(First, alignment_collect_system)
        .add_systems(PreUpdate, alignment_schedule_system);
    } else {
        app.insert_resource(AlignmentConfig {
            aligner: AlignerConfig {
                backend: crate::mafft::AlignerBackend::Mafft,
                mafft_bin: String::new(),
                mafft_fast: false,
                mafft_threads_per_job: 1,
                mafft_max_jobs: 1,
            },
        })
        .insert_resource(AlignmentStores {
            query_map: Arc::new(HashMap::new()),
            ref_map: Arc::new(HashMap::new()),
        })
        .insert_resource(AlignmentSeedJobs(Vec::new()));
    }

    if let Some(hcfg) = hmmer.take() {
        app.insert_resource(HmmerConfig {
            bin: hcfg.hmmscan_bin.clone(),
            db_path: hcfg.db_path.clone(),
            max_jobs: hcfg.max_jobs.max(1),
            top_n: hcfg.top_n,
            max_ievalue: hcfg.max_ievalue,
            threads_per_job: hcfg.threads_per_job.max(1),
        })
        .insert_resource(HmmerSeedJobs(hcfg.items))
        .add_systems(Startup, hmmer_seed_system)
        .add_systems(First, hmmer_collect_system)
        .add_systems(Update, hmmer_schedule_system);
    } else {
        app.insert_resource(HmmerConfig {
            bin: String::new(),
            db_path: String::new(),
            max_jobs: 1,
            top_n: 1,
            max_ievalue: None,
            threads_per_job: 1,
        });
    }

    app.add_systems(PostUpdate, heavy_finish_system)
        .add_systems(Last, cpu_budget_log_system);

    loop {
        app.update();
        let state = app.world().get_resource::<CompletionState>().unwrap();
        if state.finished {
            break;
        }
    }

    let world = app.world_mut();
    let align_map = world
        .remove_resource::<AlignmentResultsRes>()
        .unwrap_or_default()
        .map;
    let align_secs = world
        .remove_resource::<AlignmentTimingRes>()
        .unwrap_or_default()
        .map;
    let hmmer_map = world
        .remove_resource::<HmmerResultsRes>()
        .unwrap_or_default()
        .map;
    let hmmer_secs = world
        .remove_resource::<HmmerTimingRes>()
        .unwrap_or_default()
        .map;

    HeavyPipelineResults {
        alignment_map: align_map,
        hmmer_map,
        alignment_secs: align_secs,
        hmmer_secs,
    }
}

fn hmmer_seed_system(mut commands: Commands, mut seeds: ResMut<HmmerSeedJobs>) {
    for (gene_id, seq) in seeds.0.drain(..) {
        commands.spawn(HmmerJob { gene_id, seq });
    }
}

fn heavy_finish_system(
    mut state: ResMut<CompletionState>,
    align_progress: Res<AlignmentProgress>,
    hmmer_progress: Res<HmmerProgress>,
) {
    let align_done = align_progress.total == 0 || align_progress.completed >= align_progress.total;
    let hmmer_done = hmmer_progress.total == 0 || hmmer_progress.completed >= hmmer_progress.total;
    if align_done && hmmer_done {
        state.finished = true;
    }
}

fn hmmer_schedule_system(
    mut commands: Commands,
    config: Res<HmmerConfig>,
    mut inflight: ResMut<HmmerInflight>,
    mut progress: ResMut<HmmerProgress>,
    mut budget: ResMut<CpuBudget>,
    mut query: Query<(Entity, &mut HmmerJob), Without<HmmerHandle>>,
) {
    if inflight.current >= config.max_jobs {
        return;
    }
    let pool = AsyncComputeTaskPool::get();
    let threads_needed = config.threads_per_job.max(1);
    let dynamic_cap = HMMER_LAUNCH_CAP.max((budget.capacity / threads_needed.max(1)).max(1));
    let mut launched_this_frame = 0usize;
    for (entity, mut job) in query.iter_mut() {
        if inflight.current >= config.max_jobs {
            break;
        }
        if launched_this_frame >= dynamic_cap {
            break;
        }
        if job.seq.is_empty() {
            progress.completed += 1;
            if let Ok(mut ecmd) = commands.get_entity(entity) {
                ecmd.insert(HmmerDone);
                ecmd.despawn();
            }
            continue;
        }
        if !budget.try_acquire(CPU_LABEL_HMMER, threads_needed) {
            break;
        }
        let gene_id = job.gene_id.clone();
        let seq = std::mem::take(&mut job.seq);
        let bin = config.bin.clone();
        let db = config.db_path.clone();
        let top_n = config.top_n;
        let max_ev = config.max_ievalue;
        let task = pool.spawn(async move {
            match run_hmmscan(&bin, &db, &gene_id, &seq) {
                Ok(mut summary) => {
                    if let Some(ev) = max_ev {
                        summary.hits.retain(|h| h.evalue <= ev);
                    }
                    if summary.hits.len() > top_n {
                        summary.hits.truncate(top_n);
                    }
                    summary.hits_count = summary.hits.len();
                    summary.top_accession = summary.hits.first().map(|h| h.accession.clone());
                    summary.top_evalue = summary.hits.first().map(|h| h.evalue);
                    Ok(summary)
                }
                Err(e) => Err(e),
            }
        });
        if let Ok(mut ecmd) = commands.get_entity(entity) {
            ecmd.insert(HmmerHandle {
                task,
                threads_used: threads_needed,
                started_at: Instant::now(),
            });
        }
        inflight.current += 1;
        launched_this_frame += 1;
    }
}

fn hmmer_collect_system(
    mut commands: Commands,
    mut inflight: ResMut<HmmerInflight>,
    mut results: ResMut<HmmerResultsRes>,
    mut timings: ResMut<HmmerTimingRes>,
    mut progress: ResMut<HmmerProgress>,
    mut budget: ResMut<CpuBudget>,
    mut query: Query<(Entity, &HmmerJob, &mut HmmerHandle)>,
) {
    for (entity, job, mut handle) in query.iter_mut() {
        if let Some(res) = block_on(poll_once(&mut handle.task)) {
            let elapsed = handle.started_at.elapsed().as_secs_f64();
            inflight.current = inflight.current.saturating_sub(1);
            budget.release(CPU_LABEL_HMMER, handle.threads_used);
            match res {
                Ok(summary) => {
                    results.map.insert(job.gene_id.clone(), summary);
                }
                Err(err) => {
                    log::warn!("hmmscan failed for {}: {}", job.gene_id, err);
                }
            }
            timings.map.insert(job.gene_id.clone(), elapsed);
            progress.completed += 1;
            if let Ok(mut ecmd) = commands.get_entity(entity) {
                ecmd.remove::<HmmerHandle>();
                ecmd.insert(HmmerDone);
                ecmd.despawn();
            }
        }
    }
}

// --------- Render pipeline ---------

#[derive(Resource)]
struct RenderConfig {
    ctx: Arc<RenderContext>,
}

#[derive(Resource, Default)]
struct RenderInflight {
    tasks: Vec<Task<Result<RenderedRecord, String>>>,
    max_in_flight: usize,
}

#[derive(Resource)]
struct RenderQueue(VecDeque<(usize, GeneMetrics)>);

#[derive(Resource, Default)]
struct RenderPending {
    map: BTreeMap<usize, RenderedRecord>,
}

#[derive(Resource)]
struct RenderBuffers {
    json: Option<BufWriter<File>>,
    csv: Option<BufWriter<File>>,
    next_index: usize,
}

#[derive(Resource)]
struct RenderProgress {
    total: usize,
    flushed: usize,
}

#[derive(Resource, Default)]
struct RenderErrors {
    first: Option<String>,
}

#[derive(Resource)]
struct RenderParquetBuffer {
    pending: Vec<ScoreCard>,
    writer: Option<ParquetMetaWriter>,
    file: Option<File>,
    row_group_size: usize,
    enabled: bool,
    metadata: Option<Vec<KeyValue>>,
}

#[derive(Resource, Clone)]
struct NotationStats {
    h_count: usize,
    h_complete: usize,
    h_fragmented: usize,
    m_count: usize,
    l_count: usize,
    l_novel: usize,
    l_artifact: usize,
    x_count: usize,
    total: usize,
    taxonomy_available: bool,
    taxonomy_status_counts: HashMap<String, usize>,
    taxonomy_considered: usize,
    taxonomy_domain_counts: HashMap<String, usize>,
    taxonomy_genus_counts: HashMap<String, usize>,
    taxonomy_expected_domain: Option<String>,
    taxonomy_warn_non_target_min_frac: f64,
    taxonomy_warn_non_target_min_hits: usize,
    taxonomy_warn_non_target_strong_frac: f64,
    taxonomy_warn_non_target_strong_hits: usize,
    taxonomy_warn_genus_min_frac: f64,
    taxonomy_warn_genus_min_hits: usize,
    taxonomy_low_coverage_frac: f64,
}

impl Default for NotationStats {
    fn default() -> Self {
        Self {
            h_count: 0,
            h_complete: 0,
            h_fragmented: 0,
            m_count: 0,
            l_count: 0,
            l_novel: 0,
            l_artifact: 0,
            x_count: 0,
            total: 0,
            taxonomy_available: true,
            taxonomy_status_counts: HashMap::new(),
            taxonomy_considered: 0,
            taxonomy_domain_counts: HashMap::new(),
            taxonomy_genus_counts: HashMap::new(),
            taxonomy_expected_domain: None,
            taxonomy_warn_non_target_min_frac: 0.05,
            taxonomy_warn_non_target_min_hits: 200,
            taxonomy_warn_non_target_strong_frac: 0.10,
            taxonomy_warn_non_target_strong_hits: 500,
            taxonomy_warn_genus_min_frac: 0.15,
            taxonomy_warn_genus_min_hits: 300,
            taxonomy_low_coverage_frac: 0.05,
        }
    }
}

impl NotationStats {
    fn from_ctx(ctx: &RenderContext) -> Self {
        let mut stats = Self::default();
        stats.taxonomy_expected_domain = ctx.taxonomy_expected_domain.clone();
        stats.taxonomy_warn_non_target_min_frac = ctx.taxonomy_warn_non_target_min_frac;
        stats.taxonomy_warn_non_target_min_hits = ctx.taxonomy_warn_non_target_min_hits;
        stats.taxonomy_warn_non_target_strong_frac = ctx.taxonomy_warn_non_target_strong_frac;
        stats.taxonomy_warn_non_target_strong_hits = ctx.taxonomy_warn_non_target_strong_hits;
        stats.taxonomy_warn_genus_min_frac = ctx.taxonomy_warn_genus_min_frac;
        stats.taxonomy_warn_genus_min_hits = ctx.taxonomy_warn_genus_min_hits;
        stats.taxonomy_low_coverage_frac = ctx.taxonomy_low_coverage_frac;
        stats
    }
}
#[derive(Debug, Clone)]
pub struct RenderSummary {
    pub total: usize,
    pub high: usize,
    #[allow(dead_code)]
    pub medium: usize,
    pub low: usize,
    pub x: usize,
    pub high_complete: usize,
    pub high_fragmented: usize,
    pub low_novel: usize,
    pub low_artifact: usize,
}

impl RenderSummary {
    fn from_stats(stats: &NotationStats) -> Self {
        Self {
            total: stats.total,
            high: stats.h_count,
            medium: stats.m_count,
            low: stats.l_count,
            x: stats.x_count,
            high_complete: stats.h_complete,
            high_fragmented: stats.h_fragmented,
            low_novel: stats.l_novel,
            low_artifact: stats.l_artifact,
        }
    }
}

const CSV_HEADER_VERBOSE: &str = "gene_id,hits_count,panel_swissprot,panel_refprot,panel_cluster,top_hit,top_bitscore,top_evalue,top_qcov,top_scov,bitscore_density,coverage_delta,coverage_ratio,subject_cov_score,subject_cov_penalty,fusion_split,structvar_multiplier,final_score,classification,homology_score,intrinsic_score,genomic_score,taxonomy_score,domains_score,domains_arch_score,orphan_domain_score,length_score,length_ratio,length_class,length_expected_min,length_expected_max,length_in_expected_range,length_panel_n,conserved_regions_score,termini_score,divergence_score,mafft_enabled,conserved_fraction,pairwise_identity,panel_pairwise_identity,divergence_ratio,sequences_aligned,query_gap_fraction,gap_run_count,max_gap_run,missing_exon_run,retained_intron_run,start_concordance,start_class,end_concordance,end_class,structvar_class,structvar_gap,structvar_left_len,structvar_right_len,structvar_cov_left,structvar_cov_right,orphan_status,taxonomy_contamination,taxonomy_support,taxonomy_considered,taxonomy_support_frac,consensus_taxon,taxonomy_consensus_rank,taxonomy_status,plugin_penalty,plugin_names,plugin_scores,plugin_penalties,plugin_metadata,warnings";

const CSV_HEADER_STANDARD: &str = "gene_id,hits_count,panel_swissprot,panel_refprot,panel_cluster,top_hit,top_bitscore,top_evalue,top_qcov,top_scov,bitscore_density,coverage_delta,coverage_ratio,subject_cov_score,subject_cov_penalty,fusion_split,structvar_multiplier,final_score,classification,mafft_enabled,conserved_fraction,pairwise_identity,panel_pairwise_identity,divergence_ratio,sequences_aligned,query_gap_fraction,gap_run_count,max_gap_run,domains_score,domains_arch_score,orphan_domain_score,structvar_class,structvar_gap,structvar_left_len,structvar_right_len,structvar_cov_left,structvar_cov_right,orphan_status,genomic_score,taxonomy_score,taxonomy_contamination,taxonomy_support,taxonomy_considered,taxonomy_support_frac,consensus_taxon,taxonomy_consensus_rank,taxonomy_status,plugin_penalty,plugin_names,plugin_scores,plugin_penalties,plugin_metadata,warnings";

fn classification_base(label: &str) -> &str {
    label.split([' ', '(', '[']).next().unwrap_or(label)
}

fn top_counts(counts: &HashMap<String, usize>, top: usize) -> Vec<(String, usize)> {
    let mut entries: Vec<(String, usize)> = counts.iter().map(|(k, v)| (k.clone(), *v)).collect();
    entries.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    entries.truncate(top);
    entries
}

fn format_top_counts(counts: &HashMap<String, usize>, denom: usize, top: usize) -> String {
    if denom == 0 || counts.is_empty() {
        return "NA".to_string();
    }
    top_counts(counts, top)
        .into_iter()
        .map(|(label, count)| {
            let pct = count as f64 / denom as f64 * 100.0;
            format!("{} {:.1}%", label, pct)
        })
        .collect::<Vec<_>>()
        .join(", ")
}

fn sanitize_meta_value(value: &str) -> String {
    value.trim().replace([',', '\n', '\r'], " ")
}

fn write_csv_metadata<W: Write>(
    writer: &mut W,
    tools: &preflight::ToolVersions,
    checksums: &Checksums,
) -> std::io::Result<()> {
    let mut parts = vec![
        format!("schema_version={}", OUTPUT_SCHEMA_VERSION),
        format!("tool={}", TOOL_NAME),
        format!("tool_version={}", TOOL_VERSION),
    ];
    if let Some(v) = &tools.diamond {
        parts.push(format!("diamond_version={}", sanitize_meta_value(v)));
    }
    if let Some(v) = &tools.mafft {
        parts.push(format!("mafft_version={}", sanitize_meta_value(v)));
    }
    if let Some(v) = &tools.hmmscan {
        parts.push(format!("hmmscan_version={}", sanitize_meta_value(v)));
    }
    if let Some(v) = &checksums.fasta_xx64 {
        parts.push(format!("fasta_xx64={}", sanitize_meta_value(v)));
    }
    if let Some(v) = &checksums.db_xx64 {
        parts.push(format!("db_xx64={}", sanitize_meta_value(v)));
    }
    writeln!(writer, "# {}", parts.join("; "))
}

fn write_json_metadata<W: Write>(
    writer: &mut W,
    tools: &preflight::ToolVersions,
    checksums: &Checksums,
) -> std::io::Result<()> {
    let meta = serde_json::json!({
        "type": "metadata",
        "schema_version": OUTPUT_SCHEMA_VERSION,
        "tool": TOOL_NAME,
        "tool_version": TOOL_VERSION,
        "tools": {
            "diamond": tools.diamond,
            "mafft": tools.mafft,
            "hmmscan": tools.hmmscan,
        },
        "checksums": {
            "fasta_xx64": checksums.fasta_xx64.clone(),
            "db_xx64": checksums.db_xx64.clone(),
        },
    });
    writeln!(writer, "{}", meta)
}

fn build_parquet_metadata(tools: &preflight::ToolVersions, checksums: &Checksums) -> Vec<KeyValue> {
    let mut meta = Vec::with_capacity(8);
    meta.push(KeyValue {
        key: "schema_version".to_string(),
        value: Some(OUTPUT_SCHEMA_VERSION.to_string()),
    });
    meta.push(KeyValue {
        key: "tool".to_string(),
        value: Some(TOOL_NAME.to_string()),
    });
    meta.push(KeyValue {
        key: "tool_version".to_string(),
        value: Some(TOOL_VERSION.to_string()),
    });
    if let Some(v) = &tools.diamond {
        meta.push(KeyValue {
            key: "diamond_version".to_string(),
            value: Some(sanitize_meta_value(v)),
        });
    }
    if let Some(v) = &tools.mafft {
        meta.push(KeyValue {
            key: "mafft_version".to_string(),
            value: Some(sanitize_meta_value(v)),
        });
    }
    if let Some(v) = &tools.hmmscan {
        meta.push(KeyValue {
            key: "hmmscan_version".to_string(),
            value: Some(sanitize_meta_value(v)),
        });
    }
    if let Some(v) = &checksums.fasta_xx64 {
        meta.push(KeyValue {
            key: "fasta_xx64".to_string(),
            value: Some(sanitize_meta_value(v)),
        });
    }
    if let Some(v) = &checksums.db_xx64 {
        meta.push(KeyValue {
            key: "db_xx64".to_string(),
            value: Some(sanitize_meta_value(v)),
        });
    }
    meta
}

fn parquet_encodings(schema: &ArrowSchema) -> Vec<Vec<Encoding>> {
    schema
        .fields
        .iter()
        .map(|f| {
            transverse(&f.data_type, |dt: &ArrowDataType| {
                match dt.to_physical_type() {
                    PhysicalType::Dictionary(_)
                    | PhysicalType::LargeBinary
                    | PhysicalType::LargeUtf8 => Encoding::RleDictionary,
                    PhysicalType::Primitive(dt) => match dt {
                        PrimitiveType::Float32
                        | PrimitiveType::Float64
                        | PrimitiveType::Float16 => Encoding::Plain,
                        _ => Encoding::RleDictionary,
                    },
                    _ => Encoding::Plain,
                }
            })
        })
        .collect()
}

struct ParquetMetaWriter {
    writer: FileWriter<File>,
    schema: ArrowSchema,
    encodings: Vec<Vec<Encoding>>,
    options: WriteOptions,
    metadata: Option<Vec<KeyValue>>,
}

impl ParquetMetaWriter {
    fn new(file: File, df: &DataFrame, metadata: Vec<KeyValue>) -> Result<Self, String> {
        let schema = df.schema().to_arrow();
        let encodings = parquet_encodings(&schema);
        let options = WriteOptions {
            write_statistics: false,
            version: Version::V2,
            compression: CompressionOptions::Zstd(None),
            data_pagesize_limit: None,
        };
        let writer = FileWriter::try_new(file, schema.clone(), options)
            .map_err(|e: PolarsError| e.to_string())?;
        Ok(Self {
            writer,
            schema,
            encodings,
            options,
            metadata: Some(metadata),
        })
    }

    fn write_batch(&mut self, df: &mut DataFrame) -> Result<(), String> {
        df.align_chunks();
        let iter = df.iter_chunks().map(Ok);
        let row_groups =
            RowGroupIterator::try_new(iter, &self.schema, self.options, self.encodings.clone())
                .map_err(|e: PolarsError| e.to_string())?;
        for group in row_groups {
            let row_group = group.map_err(|e: PolarsError| e.to_string())?;
            self.writer
                .write(row_group)
                .map_err(|e: PolarsError| e.to_string())?;
        }
        Ok(())
    }

    fn finish(&mut self) -> Result<(), String> {
        let meta = self.metadata.take();
        self.writer
            .end(meta)
            .map_err(|e: PolarsError| e.to_string())?;
        Ok(())
    }
}

pub fn run_render_pipeline(
    out_dir: &str,
    metrics: &[GeneMetrics],
    ctx: Arc<RenderContext>,
    max_jobs: usize,
    tools: &preflight::ToolVersions,
    checksums: &Checksums,
    output: RenderOutputConfig,
) -> Result<RenderSummary, Box<dyn std::error::Error>> {
    let json_path = Path::new(out_dir).join("qc_report.jsonl");
    let csv_path = Path::new(out_dir).join("qc_summary.csv");
    let parquet_path = Path::new(out_dir).join("qc_summary.parquet");
    let mut json_writer = None;
    let mut csv_writer = None;
    let mut parquet_file = None;
    let mut parquet_meta = None;

    let json_has_content = json_path.metadata().map(|m| m.len() > 0).unwrap_or(false);
    if output.jsonl {
        let json_file = if output.resume && json_path.exists() && json_has_content {
            std::fs::OpenOptions::new().append(true).open(&json_path)?
        } else {
            File::create(&json_path)?
        };
        let mut writer = BufWriter::new(json_file);
        if !(output.resume && json_has_content) {
            write_json_metadata(&mut writer, tools, checksums)?;
        }
        json_writer = Some(writer);
    }

    let csv_has_content = csv_path.metadata().map(|m| m.len() > 0).unwrap_or(false);
    if output.csv {
        let csv_file = if output.resume && csv_path.exists() && csv_has_content {
            std::fs::OpenOptions::new().append(true).open(&csv_path)?
        } else {
            File::create(&csv_path)?
        };
        let mut writer = BufWriter::new(csv_file);
        if !(output.resume && csv_has_content) {
            write_csv_metadata(&mut writer, tools, checksums)?;
            if ctx.csv_verbose {
                writeln!(writer, "{}", CSV_HEADER_VERBOSE)?;
            } else {
                writeln!(writer, "{}", CSV_HEADER_STANDARD)?;
            }
        }
        csv_writer = Some(writer);
    }

    let mut parquet_enabled = output.parquet;
    if output.parquet {
        let parquet_has_content = parquet_path
            .metadata()
            .map(|m| m.len() > 0)
            .unwrap_or(false);
        if output.resume && parquet_path.exists() && parquet_has_content {
            parquet_enabled = false;
            log::warn!("resume: parquet append not supported; skipping parquet output");
        } else {
            parquet_file = Some(File::create(&parquet_path)?);
            parquet_meta = Some(build_parquet_metadata(tools, checksums));
        }
    }
    let seeds: VecDeque<(usize, GeneMetrics)> = metrics.iter().cloned().enumerate().collect();
    AsyncComputeTaskPool::get_or_init(|| {
        TaskPoolBuilder::new().num_threads(max_jobs.max(1)).build()
    });
    let mut app = App::new();
    app.insert_resource(CompletionState::default())
        .insert_resource(RenderConfig { ctx: ctx.clone() })
        .insert_resource(RenderQueue(seeds))
        .insert_resource(RenderInflight {
            tasks: Vec::with_capacity(max_jobs.max(1)),
            max_in_flight: max_jobs.max(1),
        })
        .insert_resource(RenderPending::default())
        .insert_resource(RenderProgress {
            total: metrics.len(),
            flushed: 0,
        })
        .insert_resource(RenderErrors::default())
        .insert_resource(RenderBuffers {
            json: json_writer,
            csv: csv_writer,
            next_index: 0,
        })
        .insert_resource(RenderParquetBuffer {
            pending: Vec::new(),
            writer: None,
            file: parquet_file,
            row_group_size: PARQUET_ROW_GROUP_SIZE,
            enabled: parquet_enabled,
            metadata: parquet_meta,
        })
        .insert_resource(NotationStats::from_ctx(&ctx))
        .add_systems(
            Update,
            (
                render_schedule_system,
                render_collect_system,
                render_flush_system,
                render_finish_system,
            ),
        );

    loop {
        app.update();
        let state = app.world().get_resource::<CompletionState>().unwrap();
        if state.finished {
            break;
        }
    }

    {
        let world = app.world_mut();
        let mut buffers = world.remove_resource::<RenderBuffers>().unwrap();
        if let Some(ref mut json) = buffers.json {
            json.flush()?;
        }
        if let Some(ref mut csv) = buffers.csv {
            csv.flush()?;
        }
    }
    let world = app.world_mut();
    if let Some(errs) = world.remove_resource::<RenderErrors>() {
        if let Some(err) = errs.first {
            return Err(err.into());
        }
    }

    let mut render_summary = RenderSummary {
        total: 0,
        high: 0,
        medium: 0,
        low: 0,
        x: 0,
        high_complete: 0,
        high_fragmented: 0,
        low_novel: 0,
        low_artifact: 0,
    };

    // Print Notation
    if let Some(stats) = world.remove_resource::<NotationStats>() {
        render_summary = RenderSummary::from_stats(&stats);
        if stats.total > 0 {
            let t = stats.total as f64;
            let h_pct = stats.h_count as f64 / t * 100.0;
            let m_pct = stats.m_count as f64 / t * 100.0;
            let c_pct = stats.h_complete as f64 / t * 100.0;
            let f_pct = stats.h_fragmented as f64 / t * 100.0;
            let l_pct = stats.l_count as f64 / t * 100.0;
            let n_pct = stats.l_novel as f64 / t * 100.0;
            let a_pct = stats.l_artifact as f64 / t * 100.0;
            let x_pct = stats.x_count as f64 / t * 100.0;

            println!("----------------------------------------------------------------");
            if stats.taxonomy_available {
                println!(
                    "AnnoQC Summary (n={}): H:{:.2}% [C:{:.2}%, F:{:.2}%], M:{:.2}%, L:{:.2}% [N:{:.2}%, A:{:.2}%], X:{:.2}%",
                    stats.total, h_pct, c_pct, f_pct, m_pct, l_pct, n_pct, a_pct, x_pct
                );
            } else {
                println!(
                    "AnnoQC Summary (n={}): H:{:.2}% [C:{:.2}%, F:{:.2}%], M:{:.2}%, L:{:.2}% [N:{:.2}%, A:{:.2}%], X:NA (taxonomy unavailable)",
                    stats.total, h_pct, c_pct, f_pct, m_pct, l_pct, n_pct, a_pct
                );
            }
            println!(
                "Configuration ({} v{}): {}",
                TOOL_NAME, TOOL_VERSION, ctx.features_string
            );
            if stats.taxonomy_available {
                let considered = stats.taxonomy_considered;
                let coverage_pct = if stats.total > 0 {
                    considered as f64 / stats.total as f64 * 100.0
                } else {
                    0.0
                };
                let consensus = *stats.taxonomy_status_counts.get("Consensus").unwrap_or(&0);
                let borrowed = *stats.taxonomy_status_counts.get("Borrowed").unwrap_or(&0);
                let coarse = *stats
                    .taxonomy_status_counts
                    .get("CoarseConsensus")
                    .unwrap_or(&0);
                let nohits = *stats.taxonomy_status_counts.get("NoHits").unwrap_or(&0);
                let insufficient = *stats
                    .taxonomy_status_counts
                    .get("InsufficientHits")
                    .unwrap_or(&0);
                let mut status_parts = vec![
                    format!("Consensus {}", consensus),
                    format!("Borrowed {}", borrowed),
                    format!("NoHits {}", nohits),
                ];
                if insufficient > 0 {
                    status_parts.push(format!("Insufficient {}", insufficient));
                }
                if coarse > 0 {
                    status_parts.push(format!("Coarse {}", coarse));
                }
                let domain_summary =
                    format_top_counts(&stats.taxonomy_domain_counts, considered, 3);
                let genus_summary = format_top_counts(&stats.taxonomy_genus_counts, considered, 3);
                println!(
                    "Taxonomy (context): hits {}/{} ({:.2}%) | Status: {} | Domain: {} | Top genera: {}",
                    considered,
                    stats.total,
                    coverage_pct,
                    status_parts.join(", "),
                    domain_summary,
                    genus_summary
                );
                let mut notes = Vec::new();
                if considered > 0
                    && (considered as f64 / stats.total.max(1) as f64)
                        < stats.taxonomy_low_coverage_frac
                {
                    notes.push(format!("coverage low ({:.1}% hits)", coverage_pct));
                }
                if let Some(expected) = stats.taxonomy_expected_domain.as_ref() {
                    let mut non_target_hits = 0usize;
                    for (domain, count) in &stats.taxonomy_domain_counts {
                        if domain != expected && domain != "Unknown" {
                            non_target_hits += *count;
                        }
                    }
                    if considered > 0 {
                        let non_target_frac = non_target_hits as f64 / considered as f64;
                        if non_target_hits >= stats.taxonomy_warn_non_target_strong_hits
                            && non_target_frac >= stats.taxonomy_warn_non_target_strong_frac
                        {
                            notes.push(format!(
                                "non-target domain high ({:.1}% / {} hits)",
                                non_target_frac * 100.0,
                                non_target_hits
                            ));
                        } else if non_target_hits >= stats.taxonomy_warn_non_target_min_hits
                            && non_target_frac >= stats.taxonomy_warn_non_target_min_frac
                        {
                            notes.push(format!(
                                "non-target domain elevated ({:.1}% / {} hits)",
                                non_target_frac * 100.0,
                                non_target_hits
                            ));
                        }
                    }
                }
                if considered > 0 && !stats.taxonomy_genus_counts.is_empty() {
                    if let Some((genus, count)) =
                        top_counts(&stats.taxonomy_genus_counts, 1).first()
                    {
                        let frac = *count as f64 / considered as f64;
                        if *count >= stats.taxonomy_warn_genus_min_hits
                            && frac >= stats.taxonomy_warn_genus_min_frac
                        {
                            notes.push(format!(
                                "dominant genus {} ({:.1}% / {} hits)",
                                genus,
                                frac * 100.0,
                                count
                            ));
                        }
                    }
                }
                if !notes.is_empty() {
                    println!("Taxonomy note: {}", notes.join("; "));
                }
            } else {
                println!("Taxonomy (context): unavailable");
            }
            println!("----------------------------------------------------------------");
        }
    }

    // Write parquet
    if let Some(mut buffer) = world.remove_resource::<RenderParquetBuffer>() {
        if buffer.enabled {
            if let Err(e) = flush_parquet_buffer(&mut buffer) {
                return Err(e.into());
            }
            let wrote_rows = buffer.writer.is_some();
            if let Some(mut writer) = buffer.writer {
                if let Err(e) = writer.finish() {
                    return Err(e.into());
                }
            }
            if !wrote_rows {
                if let Some(file) = buffer.file.take() {
                    drop(file);
                    let path = Path::new(out_dir).join("qc_summary.parquet");
                    std::fs::remove_file(path).ok();
                }
            }
        }
    }

    Ok(render_summary)
}

fn cards_to_dataframe(cards: &[ScoreCard]) -> PolarsResult<DataFrame> {
    let gene_id = Series::new(
        "gene_id",
        cards.iter().map(|c| c.gene_id.as_str()).collect::<Vec<_>>(),
    );
    let hits_count = Series::new(
        "hits_count",
        cards
            .iter()
            .map(|c| c.hits_count as u32)
            .collect::<Vec<_>>(),
    );
    let final_score = Series::new(
        "final_score",
        cards.iter().map(|c| c.final_score).collect::<Vec<_>>(),
    );
    let classification = Series::new(
        "classification",
        cards
            .iter()
            .map(|c| c.classification.as_str())
            .collect::<Vec<_>>(),
    );

    DataFrame::new(vec![
        gene_id,
        hits_count,
        final_score,
        classification,
        Series::new(
            "top_hit",
            cards.iter().map(|c| c.top_hit.as_str()).collect::<Vec<_>>(),
        ),
        Series::new(
            "top_bitscore",
            cards.iter().map(|c| c.top_bitscore).collect::<Vec<_>>(),
        ),
        Series::new(
            "top_evalue",
            cards.iter().map(|c| c.top_evalue).collect::<Vec<_>>(),
        ),
        Series::new(
            "top_qcov",
            cards.iter().map(|c| c.top_qcov).collect::<Vec<_>>(),
        ),
        Series::new(
            "top_scov",
            cards.iter().map(|c| c.top_scov).collect::<Vec<_>>(),
        ),
        Series::new(
            "homology_score",
            cards.iter().map(|c| c.homology_score).collect::<Vec<_>>(),
        ),
        Series::new(
            "intrinsic_score",
            cards.iter().map(|c| c.intrinsic_score).collect::<Vec<_>>(),
        ),
        Series::new(
            "genomic_score",
            cards.iter().map(|c| c.genomic_score).collect::<Vec<_>>(),
        ),
        Series::new(
            "taxonomy_score",
            cards
                .iter()
                .map(|c| c.taxonomy_score.unwrap_or(0.0))
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "domains_score",
            cards
                .iter()
                .map(|c| c.domains_score.unwrap_or(0.0))
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "divergence_score",
            cards.iter().map(|c| c.divergence_score).collect::<Vec<_>>(),
        ),
        Series::new(
            "conserved_regions_score",
            cards
                .iter()
                .map(|c| c.conserved_regions_score)
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "length_score",
            cards.iter().map(|c| c.length_score).collect::<Vec<_>>(),
        ),
        Series::new(
            "mafft_enabled",
            cards.iter().map(|c| c.mafft_enabled).collect::<Vec<_>>(),
        ),
        Series::new(
            "conserved_fraction",
            cards
                .iter()
                .map(|c| c.conserved_fraction)
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "pairwise_identity",
            cards
                .iter()
                .map(|c| c.pairwise_identity)
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "panel_pairwise_identity",
            cards
                .iter()
                .map(|c| c.panel_pairwise_identity)
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "divergence_ratio",
            cards.iter().map(|c| c.divergence_ratio).collect::<Vec<_>>(),
        ),
        Series::new(
            "genomic_introns",
            cards
                .iter()
                .map(|c| c.genomic_introns.map(|v| v as u32))
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "genomic_splice_canonical",
            cards
                .iter()
                .map(|c| c.genomic_splice_canonical.map(|v| v as u32))
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "genomic_splice_noncanonical",
            cards
                .iter()
                .map(|c| c.genomic_splice_noncanonical.map(|v| v as u32))
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "genomic_splice_weird",
            cards
                .iter()
                .map(|c| c.genomic_splice_weird.map(|v| v as u32))
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "plugin_penalty",
            cards.iter().map(|c| c.plugin_penalty).collect::<Vec<_>>(),
        ),
        Series::new(
            "plugin_count",
            cards
                .iter()
                .map(|c| c.plugin_count as u32)
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "plugin_names",
            cards
                .iter()
                .map(|c| c.plugin_names.as_str())
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "plugin_scores",
            cards
                .iter()
                .map(|c| c.plugin_scores.as_str())
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "plugin_penalties",
            cards
                .iter()
                .map(|c| c.plugin_penalties.as_str())
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "plugin_metadata",
            cards
                .iter()
                .map(|c| c.plugin_metadata.as_str())
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "warnings",
            cards
                .iter()
                .map(|c| c.warnings.as_str())
                .collect::<Vec<_>>(),
        ),
    ])
}

fn flush_parquet_buffer(buffer: &mut RenderParquetBuffer) -> Result<(), String> {
    if buffer.pending.is_empty() {
        return Ok(());
    }
    let mut df = cards_to_dataframe(&buffer.pending).map_err(|e| e.to_string())?;
    if buffer.writer.is_none() {
        let file = buffer
            .file
            .take()
            .ok_or_else(|| "parquet file missing".to_string())?;
        let meta = buffer
            .metadata
            .take()
            .ok_or_else(|| "parquet metadata missing".to_string())?;
        let writer = ParquetMetaWriter::new(file, &df, meta)?;
        buffer.writer = Some(writer);
    }
    if let Some(writer) = buffer.writer.as_mut() {
        writer.write_batch(&mut df)?;
    }
    buffer.pending.clear();
    Ok(())
}

fn render_schedule_system(
    config: Res<RenderConfig>,
    mut queue: ResMut<RenderQueue>,
    mut inflight: ResMut<RenderInflight>,
    errors: Res<RenderErrors>,
) {
    if errors.first.is_some() {
        return;
    }
    let pool = AsyncComputeTaskPool::get();
    while inflight.tasks.len() < inflight.max_in_flight {
        let Some((index, metrics)) = queue.0.pop_front() else {
            break;
        };
        let ctx = Arc::clone(&config.ctx);
        let task = pool.spawn(async move { render_gene_record(index, &metrics, &ctx) });
        inflight.tasks.push(task);
    }
}

fn render_collect_system(
    mut inflight: ResMut<RenderInflight>,
    mut pending: ResMut<RenderPending>,
    mut errors: ResMut<RenderErrors>,
) {
    let mut i = 0;
    while i < inflight.tasks.len() {
        if let Some(res) = block_on(poll_once(&mut inflight.tasks[i])) {
            match res {
                Ok(record) => {
                    pending.map.insert(record.index, record);
                }
                Err(err) => {
                    if errors.first.is_none() {
                        errors.first = Some(err);
                    }
                }
            }
            drop(inflight.tasks.swap_remove(i));
        } else {
            i += 1;
        }
    }
}

fn render_flush_system(
    mut pending: ResMut<RenderPending>,
    mut buffers: ResMut<RenderBuffers>,
    mut parquet_buffer: ResMut<RenderParquetBuffer>,
    mut stats: ResMut<NotationStats>,
    mut progress: ResMut<RenderProgress>,
    mut errors: ResMut<RenderErrors>,
) {
    if errors.first.is_some() {
        return;
    }
    while let Some(record) = pending.map.remove(&buffers.next_index) {
        if let Some(writer) = buffers.json.as_mut() {
            if let Err(e) = writeln!(writer, "{}", record.json_line) {
                if errors.first.is_none() {
                    errors.first = Some(format!("write json: {}", e));
                }
                break;
            }
        }
        if let Some(writer) = buffers.csv.as_mut() {
            if let Err(e) = writeln!(writer, "{}", record.csv_line) {
                if errors.first.is_none() {
                    errors.first = Some(format!("write csv: {}", e));
                }
                break;
            }
        }
        if let Some(card) = record.card {
            stats.total += 1;
            if card.taxonomy_status == "disabled" || card.taxonomy_status == "NoResolver" {
                stats.taxonomy_available = false;
            }
            if !card.taxonomy_status.is_empty() {
                *stats
                    .taxonomy_status_counts
                    .entry(card.taxonomy_status.clone())
                    .or_insert(0) += 1;
            }
            let base = classification_base(&card.classification);
            match base {
                "High" => {
                    stats.h_count += 1;
                    let len_ok = card.length_class == "InRange" || card.length_class.is_empty();
                    let orphan_ok = card.orphan_status == "None" || card.orphan_status.is_empty();
                    if len_ok && orphan_ok {
                        stats.h_complete += 1;
                    } else {
                        stats.h_fragmented += 1;
                    }
                }
                "Medium" => {
                    stats.m_count += 1;
                }
                "Low" => {
                    stats.l_count += 1;
                    if card.intrinsic_score >= 0.5 {
                        stats.l_novel += 1;
                    } else {
                        stats.l_artifact += 1;
                    }
                }
                "NoData" => {
                    stats.x_count += 1;
                }
                _ => {
                    stats.l_count += 1;
                    if card.intrinsic_score >= 0.5 {
                        stats.l_novel += 1;
                    } else {
                        stats.l_artifact += 1;
                    }
                }
            }
            // Simple X heuristic for now: Excluded if taxonomy contamination is very high (>0.9)
            if let Some(contam) = card.taxonomy_contamination {
                if contam > 0.9 {
                    stats.x_count += 1;
                }
            }
            if matches!(
                card.taxonomy_status.as_str(),
                "Consensus" | "Borrowed" | "CoarseConsensus"
            ) {
                stats.taxonomy_considered += 1;
                if !card.taxonomy_domain.is_empty() {
                    *stats
                        .taxonomy_domain_counts
                        .entry(card.taxonomy_domain.clone())
                        .or_insert(0) += 1;
                }
                if !card.taxonomy_genus.is_empty() {
                    *stats
                        .taxonomy_genus_counts
                        .entry(card.taxonomy_genus.clone())
                        .or_insert(0) += 1;
                }
            }
            if parquet_buffer.enabled {
                parquet_buffer.pending.push(card);
                if parquet_buffer.pending.len() >= parquet_buffer.row_group_size {
                    if let Err(e) = flush_parquet_buffer(&mut parquet_buffer) {
                        if errors.first.is_none() {
                            errors.first = Some(format!("write parquet: {}", e));
                        }
                        break;
                    }
                }
            }
        }
        buffers.next_index += 1;
        progress.flushed += 1;
    }
}

fn render_finish_system(
    queue: Res<RenderQueue>,
    pending: Res<RenderPending>,
    inflight: Res<RenderInflight>,
    progress: Res<RenderProgress>,
    errors: Res<RenderErrors>,
    mut state: ResMut<CompletionState>,
) {
    let queue_empty = queue.0.is_empty();
    let inflight_empty = inflight.tasks.is_empty();
    let pending_empty = pending.map.is_empty();
    if progress.flushed >= progress.total && queue_empty && inflight_empty && pending_empty {
        state.finished = true;
        return;
    }
    if errors.first.is_some() && inflight_empty {
        state.finished = true;
    }
}
