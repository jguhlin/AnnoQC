use std::collections::HashMap;

use tempfile::TempDir;

use crate::ecs::GeneMetrics;
use crate::taxonomy::TaxonomyEvidence;
use crate::*;

#[test]
fn export_high_sequences_writes_only_high() {
    let tmp = TempDir::new().unwrap();
    let out_dir = tmp.path().to_str().unwrap();
    let metrics = vec![
        GeneMetrics {
            gene_id: "gene_high".into(),
            length: 10,
            hits: 0,
        },
        GeneMetrics {
            gene_id: "gene_low".into(),
            length: 10,
            hits: 0,
        },
    ];
    let mut seqs = HashMap::new();
    seqs.insert(
        "gene_high".into(),
        (metrics::IntrinsicMetrics::default(), b"MKTAA".to_vec()),
    );
    seqs.insert(
        "gene_low".into(),
        (metrics::IntrinsicMetrics::default(), b"AAAAA".to_vec()),
    );
    let mut scores = HashMap::new();
    scores.insert("gene_high".into(), (0.9, "High".into()));
    scores.insert("gene_low".into(), (0.4, "Low".into()));
    let res = scoring_support::export_high_sequences(out_dir, None, &metrics, &seqs, &scores)
        .expect("export succeeds");
    assert!(res.is_some());
    let (path, count) = res.unwrap();
    assert_eq!(count, 1);
    let fasta = std::fs::read_to_string(path).unwrap();
    assert!(fasta.contains("gene_high"));
    assert!(!fasta.contains("gene_low"));
}

#[test]
fn propagate_transcript_taxonomy_borrows_best() {
    use crate::taxonomy::TaxonomyDetail;
    use crate::taxonomy_support::propagate_transcript_taxonomy;

    let metrics = vec![
        GeneMetrics {
            gene_id: "gene1.t1".into(),
            length: 100,
            hits: 0,
        },
        GeneMetrics {
            gene_id: "gene1.t2".into(),
            length: 100,
            hits: 0,
        },
        GeneMetrics {
            gene_id: "gene2".into(),
            length: 100,
            hits: 0,
        },
    ];
    let mut map: HashMap<String, Option<TaxonomyEvidence>> = HashMap::new();
    let mut ev = TaxonomyEvidence::default();
    ev.detail = TaxonomyDetail::Consensus;
    ev.congruence_score = 0.9;
    ev.contamination_score = 0.1;
    ev.support = 5;
    ev.considered = 5;
    ev.support_fraction = 1.0;
    map.insert("gene1.t1".into(), Some(ev));
    map.insert("gene1.t2".into(), None);
    map.insert("gene2".into(), None);
    propagate_transcript_taxonomy(&mut map, &metrics);
    let borrowed = map.get("gene1.t2").and_then(|v| v.clone()).unwrap();
    assert_eq!(borrowed.detail, TaxonomyDetail::Borrowed);
    assert_eq!(borrowed.support, 0);
    assert!(map.get("gene2").unwrap().is_none());
}

#[test]
fn dbinfo_text_detection_handles_positive_and_negative() {
    let ok = "Database sequences\nTaxon count: 100";
    assert!(prepare::support::dbinfo_text_has_taxonomy(ok));
    let bad = "Database sequences\nNo taxonomy found";
    assert!(!prepare::support::dbinfo_text_has_taxonomy(bad));
}
