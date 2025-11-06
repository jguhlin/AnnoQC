use std::collections::{HashMap, VecDeque};

use bevy_app::{App, Startup, Update};
use bevy_ecs::prelude::*;
use bevy_tasks::{AsyncComputeTaskPool, Task, TaskPoolBuilder};
use futures_lite::future::{block_on, poll_once};
use needletail::parse_fastx_file;

#[derive(Resource, Clone)]
pub struct EcsConfig {
    pub fasta_path: String,
    pub diamond_tsv: String,
    pub threads: usize,
}

#[derive(Debug, Clone)]
pub struct GeneMetrics {
    pub gene_id: String,
    pub length: usize,
    pub hits: usize,
}

#[derive(Resource, Default)]
pub struct CompletionState {
    pub finished: bool,
}

#[derive(Resource, Default)]
struct WorkQueue {
    ids: VecDeque<String>,
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

fn intake_system(mut commands: Commands, cfg: Res<EcsConfig>) {
    // Parse FASTA ids and lengths
    let mut ids = VecDeque::new();
    let mut lengths = HashMap::<String, usize>::new();
    let mut reader = parse_fastx_file(&cfg.fasta_path).expect("open fasta");
    while let Some(record) = reader.next() {
        let rec = record.expect("valid record");
        let id = String::from_utf8_lossy(rec.id()).to_string();
        let gene = id.split_whitespace().next().unwrap_or("").to_string();
        if gene.is_empty() {
            continue;
        }
        ids.push_back(gene.clone());
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

    // Convert ids to a queue; lengths needed by tasks
    commands.insert_resource(WorkQueue { ids, hit_counts });
    commands.insert_resource(Results::default());
    commands.insert_resource(InFlight {
        tasks: Vec::new(),
        max_in_flight: cfg.threads.max(1),
    });
    commands.insert_resource(LengthMap(lengths));
}

fn schedule_system(
    mut queue: ResMut<WorkQueue>,
    lengths: Res<LengthMap>,
    mut inflight: ResMut<InFlight>,
) {
    let pool = AsyncComputeTaskPool::get();
    while inflight.tasks.len() < inflight.max_in_flight {
        let Some(gene) = queue.ids.pop_front() else {
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

fn collect_system(mut inflight: ResMut<InFlight>, mut results: ResMut<Results>) {
    let mut i = 0;
    while i < inflight.tasks.len() {
        let done = block_on(poll_once(&mut inflight.tasks[i]));
        if let Some(metrics) = done {
            results.records.push(metrics);
            let _ = inflight.tasks.swap_remove(i);
        } else {
            i += 1;
        }
    }
}

fn finish_system(
    queue: Res<WorkQueue>,
    inflight: Res<InFlight>,
    mut state: ResMut<CompletionState>,
) {
    if queue.ids.is_empty() && inflight.tasks.is_empty() {
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
