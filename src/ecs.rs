use std::collections::{BTreeMap, HashMap, VecDeque};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;

use bevy_app::{App, Startup, Update};
use bevy_ecs::prelude::*;
use bevy_tasks::{AsyncComputeTaskPool, Task, TaskPoolBuilder};
use futures_lite::future::{block_on, poll_once};
use needletail::parse_fastx_file;
use std::time::{Duration, Instant};

use crate::hmmer::{run_hmmscan, HmmscanSummary};
use crate::mafft::{run_mafft, AlignmentMetrics};
use crate::render_gene_record;
use crate::{taxonomy, RenderContext, RenderedRecord};

#[derive(Resource, Clone)]
pub struct EcsConfig {
    pub fasta_path: String,
    pub diamond_tsv: String,
    pub threads: usize,
    pub log_json: bool,
}

#[derive(Debug, Clone)]
pub struct GeneMetrics {
    pub gene_id: String,
    #[allow(dead_code)]
    pub length: usize,
    pub hits: usize,
}

#[derive(Resource, Default)]
pub struct CompletionState {
    pub finished: bool,
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
        tasks: Vec::new(),
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
}

#[derive(Resource)]
struct AlignmentConfig {
    mafft_bin: String,
    max_jobs: usize,
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

pub fn run_alignment_pipeline(
    mafft_bin: &str,
    jobs: Vec<(String, Vec<String>)>,
    query_map: HashMap<String, Vec<u8>>,
    ref_map: HashMap<String, Vec<u8>>,
    max_jobs: usize,
) -> HashMap<String, AlignmentMetrics> {
    let total_jobs = jobs.len();
    let job_data = jobs;
    AsyncComputeTaskPool::get_or_init(|| {
        TaskPoolBuilder::new().num_threads(max_jobs.max(1)).build()
    });
    let mut app = App::new();
    app.insert_resource(CompletionState::default())
        .insert_resource(AlignmentConfig {
            mafft_bin: mafft_bin.to_string(),
            max_jobs: max_jobs.max(1),
        })
        .insert_resource(AlignmentStores {
            query_map: Arc::new(query_map),
            ref_map: Arc::new(ref_map),
        })
        .insert_resource(AlignmentResultsRes::default())
        .insert_resource(AlignmentInflight::default())
        .insert_resource(AlignmentProgress {
            total: total_jobs,
            completed: 0,
        })
        .insert_resource(AlignmentDiag {
            total: total_jobs,
            skipped_no_panel: 0,
            skipped_no_query: 0,
            skipped_no_refs: 0,
            scheduled: 0,
            last_log: std::time::Instant::now(),
            next_pct: 5,
        })
        .insert_resource(AlignmentSeedJobs(job_data))
        .add_systems(Startup, alignment_seed_system)
        .add_systems(
            Update,
            (
                alignment_schedule_system,
                alignment_collect_system,
                alignment_finish_system,
            ),
        );

    loop {
        app.update();
        let state = app.world().get_resource::<CompletionState>().unwrap();
        if state.finished {
            break;
        }
    }

    let results = app
        .world_mut()
        .remove_resource::<AlignmentResultsRes>()
        .unwrap_or_default();
    results.map
}

fn alignment_schedule_system(
    mut commands: Commands,
    config: Res<AlignmentConfig>,
    stores: Res<AlignmentStores>,
    mut inflight: ResMut<AlignmentInflight>,
    mut progress: ResMut<AlignmentProgress>,
    mut diag: ResMut<AlignmentDiag>,
    mut query: Query<(Entity, &AlignmentJob), Without<AlignmentHandle>>,
) {
    if inflight.current >= config.max_jobs {
        return;
    }
    let pool = AsyncComputeTaskPool::get();
    for (entity, job) in query.iter_mut() {
        if inflight.current >= config.max_jobs {
            break;
        }
        if job.panel_ids.is_empty() {
            progress.completed += 1;
            if let Ok(mut ecmd) = commands.get_entity(entity) {
                ecmd.insert(AlignmentDone);
            }
            diag.skipped_no_panel += 1;
            continue;
        }
        let Some(qseq) = stores.query_map.get(&job.gene_id).cloned() else {
            log::warn!("mafft skipped {}, query sequence missing", job.gene_id);
            progress.completed += 1;
            if let Ok(mut ecmd) = commands.get_entity(entity) {
                ecmd.insert(AlignmentDone);
            }
            diag.skipped_no_query += 1;
            continue;
        };
        let mut hits: HashMap<String, Vec<u8>> = HashMap::new();
        for id in job.panel_ids.iter() {
            let key = taxonomy::canonical_accession(id);
            if let Some(seq) = stores.ref_map.get(&key) {
                hits.insert(key.clone(), seq.clone());
            }
        }
        if hits.is_empty() {
            progress.completed += 1;
            if let Ok(mut ecmd) = commands.get_entity(entity) {
                ecmd.insert(AlignmentDone);
            }
            diag.skipped_no_refs += 1;
            continue;
        }
        let gene_id = job.gene_id.clone();
        let bin = config.mafft_bin.clone();
        let task = pool.spawn(async move { run_mafft(&bin, &gene_id, &qseq, &hits) });
        if let Ok(mut ecmd) = commands.get_entity(entity) {
            ecmd.insert(AlignmentHandle { task });
        }
        inflight.current += 1;
        diag.scheduled += 1;
    }
}

fn alignment_collect_system(
    mut commands: Commands,
    mut inflight: ResMut<AlignmentInflight>,
    mut results: ResMut<AlignmentResultsRes>,
    mut progress: ResMut<AlignmentProgress>,
    mut diag: ResMut<AlignmentDiag>,
    mut query: Query<(Entity, &AlignmentJob, &mut AlignmentHandle)>,
) {
    for (entity, job, mut handle) in query.iter_mut() {
        if let Some(res) = block_on(poll_once(&mut handle.task)) {
            inflight.current = inflight.current.saturating_sub(1);
            match res {
                Ok(metrics) => {
                    results.map.insert(job.gene_id.clone(), metrics);
                }
                Err(err) => {
                    log::warn!("mafft failed for {}: {}", job.gene_id, err);
                }
            }
            progress.completed += 1;
            if let Ok(mut ecmd) = commands.get_entity(entity) {
                ecmd.remove::<AlignmentHandle>();
                ecmd.insert(AlignmentDone);
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
            "mafft progress: completed={}/{} ({}%), inflight={}, scheduled={}, skipped_no_panel={}, skipped_no_query={}, skipped_no_refs={}",
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

fn alignment_finish_system(mut state: ResMut<CompletionState>, progress: Res<AlignmentProgress>) {
    if progress.completed >= progress.total {
        state.finished = true;
    }
}
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
}

#[derive(Resource, Default)]
struct HmmerResultsRes {
    map: HashMap<String, HmmscanSummary>,
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

pub fn run_hmmer_pipeline(
    hmmscan_bin: &str,
    db_path: &str,
    items: Vec<(String, Vec<u8>)>,
    max_jobs: usize,
    top_n: usize,
    max_ievalue: Option<f64>,
) -> HashMap<String, HmmscanSummary> {
    let total = items.len();
    let seeds = items;
    AsyncComputeTaskPool::get_or_init(|| {
        TaskPoolBuilder::new().num_threads(max_jobs.max(1)).build()
    });
    let mut app = App::new();
    app.insert_resource(CompletionState::default())
        .insert_resource(HmmerConfig {
            bin: hmmscan_bin.to_string(),
            db_path: db_path.to_string(),
            max_jobs: max_jobs.max(1),
            top_n,
            max_ievalue,
        })
        .insert_resource(HmmerResultsRes::default())
        .insert_resource(HmmerInflight::default())
        .insert_resource(HmmerProgress {
            total,
            completed: 0,
        })
        .insert_resource(HmmerSeedJobs(seeds))
        .add_systems(Startup, hmmer_seed_system)
        .add_systems(
            Update,
            (
                hmmer_schedule_system,
                hmmer_collect_system,
                hmmer_finish_system,
            ),
        );

    loop {
        app.update();
        let state = app.world().get_resource::<CompletionState>().unwrap();
        if state.finished {
            break;
        }
    }

    let results = app
        .world_mut()
        .remove_resource::<HmmerResultsRes>()
        .unwrap_or_default();
    results.map
}

fn hmmer_seed_system(mut commands: Commands, mut seeds: ResMut<HmmerSeedJobs>) {
    for (gene_id, seq) in seeds.0.drain(..) {
        commands.spawn(HmmerJob { gene_id, seq });
    }
}

fn hmmer_schedule_system(
    mut commands: Commands,
    config: Res<HmmerConfig>,
    mut inflight: ResMut<HmmerInflight>,
    mut progress: ResMut<HmmerProgress>,
    mut query: Query<(Entity, &mut HmmerJob), Without<HmmerHandle>>,
) {
    if inflight.current >= config.max_jobs {
        return;
    }
    let pool = AsyncComputeTaskPool::get();
    for (entity, mut job) in query.iter_mut() {
        if inflight.current >= config.max_jobs {
            break;
        }
        if job.seq.is_empty() {
            progress.completed += 1;
            if let Ok(mut ecmd) = commands.get_entity(entity) {
                ecmd.insert(HmmerDone);
            }
            continue;
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
            ecmd.insert(HmmerHandle { task });
        }
        inflight.current += 1;
    }
}

fn hmmer_collect_system(
    mut commands: Commands,
    mut inflight: ResMut<HmmerInflight>,
    mut results: ResMut<HmmerResultsRes>,
    mut progress: ResMut<HmmerProgress>,
    mut query: Query<(Entity, &HmmerJob, &mut HmmerHandle)>,
) {
    for (entity, job, mut handle) in query.iter_mut() {
        if let Some(res) = block_on(poll_once(&mut handle.task)) {
            inflight.current = inflight.current.saturating_sub(1);
            match res {
                Ok(summary) => {
                    results.map.insert(job.gene_id.clone(), summary);
                }
                Err(err) => {
                    log::warn!("hmmscan failed for {}: {}", job.gene_id, err);
                }
            }
            progress.completed += 1;
            if let Ok(mut ecmd) = commands.get_entity(entity) {
                ecmd.remove::<HmmerHandle>();
                ecmd.insert(HmmerDone);
            }
        }
    }
}

fn hmmer_finish_system(mut state: ResMut<CompletionState>, progress: Res<HmmerProgress>) {
    if progress.completed >= progress.total {
        state.finished = true;
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
    json: BufWriter<File>,
    csv: BufWriter<File>,
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

const CSV_HEADER_VERBOSE: &str = "gene_id,hits_count,top_hit,top_bitscore,top_evalue,top_qcov,top_scov,bitscore_density,coverage_delta,coverage_ratio,fusion_split,final_score,classification,homology_score,intrinsic_score,taxonomy_score,domains_score,domains_arch_score,orphan_domain_score,length_score,length_ratio,length_class,mafft_enabled,conserved_fraction,pairwise_identity,sequences_aligned,query_gap_fraction,gap_run_count,max_gap_run,missing_exon_run,retained_intron_run,start_concordance,start_class,structvar_class,structvar_gap,structvar_left_len,structvar_right_len,structvar_cov_left,structvar_cov_right,orphan_status,taxonomy_contamination,taxonomy_support,taxonomy_considered,consensus_taxon,taxonomy_status,warnings";

const CSV_HEADER_STANDARD: &str = "gene_id,hits_count,top_hit,top_bitscore,top_evalue,top_qcov,top_scov,bitscore_density,coverage_delta,coverage_ratio,fusion_split,final_score,classification,mafft_enabled,conserved_fraction,pairwise_identity,sequences_aligned,query_gap_fraction,gap_run_count,max_gap_run,domains_score,domains_arch_score,orphan_domain_score,structvar_class,structvar_gap,structvar_left_len,structvar_right_len,structvar_cov_left,structvar_cov_right,orphan_status,taxonomy_score,taxonomy_contamination,taxonomy_support,taxonomy_considered,consensus_taxon,taxonomy_status,warnings";

pub fn run_render_pipeline(
    out_dir: &str,
    metrics: &[GeneMetrics],
    ctx: Arc<RenderContext>,
    max_jobs: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let json_path = Path::new(out_dir).join("qc_report.jsonl");
    let csv_path = Path::new(out_dir).join("qc_summary.csv");
    let json_file = File::create(json_path)?;
    let csv_file = File::create(csv_path)?;
    let mut csv_writer = BufWriter::new(csv_file);
    if ctx.csv_verbose {
        writeln!(csv_writer, "{}", CSV_HEADER_VERBOSE)?;
    } else {
        writeln!(csv_writer, "{}", CSV_HEADER_STANDARD)?;
    }
    let json_writer = BufWriter::new(json_file);
    let seeds: VecDeque<(usize, GeneMetrics)> = metrics.iter().cloned().enumerate().collect();
    AsyncComputeTaskPool::get_or_init(|| {
        TaskPoolBuilder::new().num_threads(max_jobs.max(1)).build()
    });
    let mut app = App::new();
    app.insert_resource(CompletionState::default())
        .insert_resource(RenderConfig { ctx })
        .insert_resource(RenderQueue(seeds))
        .insert_resource(RenderInflight {
            tasks: Vec::new(),
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
        buffers.json.flush()?;
        buffers.csv.flush()?;
    }
    let world = app.world_mut();
    if let Some(errs) = world.remove_resource::<RenderErrors>() {
        if let Some(err) = errs.first {
            return Err(err.into());
        }
    }
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
    mut progress: ResMut<RenderProgress>,
    mut errors: ResMut<RenderErrors>,
) {
    if errors.first.is_some() {
        return;
    }
    while let Some(record) = pending.map.remove(&buffers.next_index) {
        if let Err(e) = writeln!(buffers.json, "{}", record.json_line) {
            if errors.first.is_none() {
                errors.first = Some(format!("write json: {}", e));
            }
            break;
        }
        if let Err(e) = writeln!(buffers.csv, "{}", record.csv_line) {
            if errors.first.is_none() {
                errors.first = Some(format!("write csv: {}", e));
            }
            break;
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
