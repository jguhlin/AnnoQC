use std::{
    fs::File,
    io::BufWriter,
    io::{BufRead, BufReader, Write},
    path::Path,
    process::Command,
    collections::HashMap,
};

use clap::{Parser, Subcommand};
use curl::easy::Easy;
use libsfasta::prelude::*;
use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Prepare, // { config: Option<String> },
    Analyze { input: String },
}

fn main() {
    env_logger::init();
    let cli = Cli::parse();
    match &cli.command {
        Commands::Prepare => prepare(), // { config } => build(config),
        Commands::Analyze { input } => analyze(input), // { config } => analyze(), // TODO: Add config option
    }
}

fn prepare() {
    let uniprot_swissprot_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz";

    // Check if uniprot_sprot.fasta.gz already exists and is > 0 bytes
    if Path::new("uniprot_sprot.fasta.gz").exists()
        && Path::new("uniprot_sprot.fasta.gz")
            .metadata()
            .unwrap()
            .len()
            > 0
    {
        log::info!("uniprot_sprot.fasta.gz already exists, skipping download");
    } else {
        // Download uniprot_sprot.fasta.gz
        log::info!("Downloading uniprot_sprot.fasta.gz");
        let mut output_buf = BufWriter::new(File::create("uniprot_sprot.fasta.gz").unwrap());

        let mut easy = Easy::new();
        easy.url(uniprot_swissprot_url).unwrap();
        easy.write_function(move |data| {
            output_buf.write_all(data).unwrap();
            Ok(data.len())
        })
        .unwrap();
        easy.perform().unwrap();
        drop(easy);
    }

    // Check if uniprot_sprot.sfasta already exists and is > 0 bytes
    if Path::new("uniprot_sprot.sfasta").exists()
        && Path::new("uniprot_sprot.sfasta").metadata().unwrap().len() > 0
    {
        log::info!("uniprot_sprot.sfasta already exists, skipping conversion");
    } else {
        // Convert uniprot_sprot.fasta.gz to uniprot_sprot.sfasta
        log::info!("Converting uniprot_sprot.fasta.gz to uniprot_sprot.sfasta");
        let converter = Converter::default().with_threads(16);
        let in_buf = BufReader::new(File::open("uniprot_sprot.fasta.gz").unwrap());
        let mut in_buf = Box::new(flate2::read::MultiGzDecoder::new(in_buf));
        let mut out_buf = Box::new(BufWriter::new(
            File::create("uniprot_sprot.sfasta").unwrap(),
        ));
        converter.convert_fasta(&mut in_buf, &mut out_buf);
    }

    log::info!("Clustering uniprot_sprot.fasta.gz");
    log::debug!("Executing diamond makedb --db uniprot_sprot --in uniprot_sprot.fasta.gz");
    let output = Command::new("diamond")
        .arg("makedb")
        .arg("--db")
        .arg("uniprot_sprot.dmnd")
        .arg("--in")
        .arg("uniprot_sprot.fasta.gz")
        .output()
        .expect("Failed to execute diamond makedb");

    log::debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    log::debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    log::debug!("Executing diamond cluster --db uniprot_sprot.dmnd -o clusters");
    let output = Command::new("diamond")
        .arg("cluster")
        .arg("--db")
        .arg("uniprot_sprot.dmnd")
        .arg("--approx-id")
        .arg("50")
        .arg("--threads")
        .arg("16")
        .arg("--member-cover")
        .arg("80")
        .arg("-o")
        .arg("clusters")
        .output()
        .expect("Failed to execute diamond cluster");

    log::debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    log::debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    log::debug!(
        "Executing diamond realign --db uniprot_sprot.dmnd --clusters clusters -o clusters.realign"
    );
    let output = Command::new("diamond")
        .arg("realign")
        .arg("--db")
        .arg("uniprot_sprot.dmnd")
        .arg("--approx-id")
        .arg("50")
        .arg("--clusters")
        .arg("clusters")
        .arg("--threads")
        .arg("16")
        .arg("-o")
        .arg("clusters.realign")
        .output()
        .expect("Failed to execute diamond realign");

    log::debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    log::debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    log::debug!("Executing diamond recluster --db uniprot_sprot.dmnd --clusters clusters.realign --approx-id 30 --member-cover 80 --threads 16 --out clusters.recluster");
    let output = Command::new("diamond")
        .arg("recluster")
        .arg("--db")
        .arg("uniprot_sprot.dmnd")
        .arg("--approx-id")
        .arg("50")
        .arg("--clusters")
        .arg("clusters.realign")
        .arg("--threads")
        .arg("16")
        .arg("--member-cover")
        .arg("80")
        .arg("--out")
        .arg("clusters.recluster")
        .output()
        .expect("Failed to execute diamond recluster");

    log::debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    log::debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));
}

fn analyze(input: &str) {
    // Open clusters.recluster file (2 columns: centroid_id, protein_id)

    let mut clusters = HashMap::new();

    let file = BufReader::new(File::open("clusters.recluster").unwrap());
    for line in file.lines() {
        let line = line.unwrap();
        let mut split = line.split_whitespace();
        let centroid_id = split.next().unwrap();
        let protein_id = split.next().unwrap();
        clusters
            .entry(centroid_id.to_string())
            .or_insert_with(Vec::new)
            .push(protein_id.to_string());
    }

    println!("Found {} clusters", clusters.len());

    // Calculate average cluster size
    let mut total_cluster_size = 0;
    for (_, cluster) in clusters.iter() {
        total_cluster_size += cluster.len();
    }
    let average_cluster_size = total_cluster_size as f64 / clusters.len() as f64;
    println!("Average cluster size: {}", average_cluster_size);

    // Calculate singletons
    let mut singletons = 0;
    for (_, cluster) in clusters.iter() {
        if cluster.len() == 1 {
            singletons += 1;
        }
    }
    println!("Singletons: {}", singletons);

    // Print out number of non-singletons
    let mut non_singletons = 0;
    for (_, cluster) in clusters.iter() {
        if cluster.len() > 1 {
            non_singletons += 1;
        }
    }
    println!("Non-singletons: {}", non_singletons);

    // Get all centroids and export to a file
    let centroids = clusters.keys().cloned().collect::<Vec<String>>();
    
    // Use SFASTA

    // Build diamond db of centroids

    // Blastp ultra sensitive against centroids

    // Parse blastp output

    // Identify those that match a cluster

    







    // Diamond blastp against centroids






}
