use crate::*;

pub(crate) fn render_gene_record(
    index: usize,
    m: &ecs::GeneMetrics,
    ctx: &RenderContext,
) -> Result<RenderedRecord, String> {
    let summary = ctx.stats.get(&m.gene_id);
    let (intrinsic, seq_bytes) = ctx
        .intrinsic_map
        .get(&m.gene_id)
        .map(|t| (t.0.clone(), t.1.clone()))
        .unwrap_or_default();
    let genomic_metrics = ctx.genomic_map.as_ref().and_then(|map| map.get(&m.gene_id));
    let taxonomy_entry = ctx.taxsum_map.get(&m.gene_id).and_then(|x| x.as_ref());
    let panel_prov = ctx
        .panel_prov_map
        .get(&m.gene_id)
        .cloned()
        .unwrap_or_default();

    let mut plugin_results = Vec::new();
    let mut plugin_penalty = 0.0;
    if !ctx.plugins.is_empty() || ctx.rhai_runtime.is_some() {
        let homology = summary.map(|s| plugins::PluginHomology {
            hits_count: s.count,
            top_hit: s.top_sseqid.clone(),
            top_bitscore: s.top_bitscore,
            top_evalue: s.top_evalue.clone(),
            top_qcov: s.top_qcov,
            top_scov: s.top_scov,
            bitscore_density: if s.top_len > 0 {
                s.top_bitscore / s.top_len as f64
            } else {
                0.0
            },
            coverage_delta: s.coverage_delta,
            coverage_ratio: s.coverage_ratio,
        });
        let intrinsic_snapshot = plugins::PluginIntrinsic {
            ambiguous_fraction: intrinsic.ambiguous_fraction,
            max_homopolymer: intrinsic.max_homopolymer,
            low_complexity_fraction: intrinsic.low_complexity_fraction,
            low_complexity_windows: intrinsic.low_complexity_windows,
            orf_start_score: intrinsic.orf_start_score,
        };
        let taxonomy_snapshot = if ctx.taxonomy_enabled {
            taxonomy_entry.map(|ev| plugins::PluginTaxonomy {
                detail: ev.detail.to_string(),
                congruence_score: ev.congruence_score,
                contamination_score: ev.contamination_score,
                support_fraction: ev.support_fraction,
                support: ev.support,
                considered: ev.considered,
                consensus_rank: ev.consensus_rank.clone(),
                consensus_taxid: ev.consensus.as_ref().map(|c| c.taxid),
                consensus_name: ev.consensus.as_ref().and_then(|c| c.name.clone()),
            })
        } else {
            None
        };
        let panel_snapshot = Some(plugins::PluginPanel {
            swissprot: panel_prov.swissprot,
            refprot: panel_prov.refprot,
            cluster: panel_prov.cluster,
        });
        let genomic_snapshot = genomic_metrics.map(|g| plugins::PluginGenomic {
            introns_total: g.introns_total,
            splice_canonical: g.splice_canonical,
            splice_noncanonical: g.splice_major_noncan + g.splice_minor,
            splice_weird: g.splice_weird,
            intron_len_min: g.intron_len_min,
            intron_len_max: g.intron_len_max,
            intron_len_avg: g.intron_len_avg,
        });
        let plugin_input = plugins::PluginInput {
            gene_id: m.gene_id.clone(),
            sequence: String::from_utf8_lossy(&seq_bytes).to_string(),
            homology,
            intrinsic: intrinsic_snapshot,
            taxonomy: taxonomy_snapshot,
            panel: panel_snapshot,
            genomic: genomic_snapshot,
        };
        for p in &ctx.plugins {
            match plugins::run_plugin(p, &plugin_input) {
                Ok(res) => {
                    if let Some(pen) = res.penalty {
                        plugin_penalty += pen;
                    }
                    plugin_results.push(res);
                }
                Err(e) => {
                    log::debug!("plugin {} failed for {}: {}", p.name, m.gene_id, e);
                }
            }
        }
        if let Some(rt) = &ctx.rhai_runtime {
            for res in rt.run(&plugin_input) {
                if let Some(pen) = res.penalty {
                    plugin_penalty += pen;
                }
                plugin_results.push(res);
            }
        }
    }
    let sanitize_csv_field = |value: String| value.replace([',', '\n', '\r'], " ");
    let plugin_names_raw = plugin_results
        .iter()
        .map(|r| r.name.clone())
        .collect::<Vec<_>>()
        .join("|");
    let plugin_scores_raw = plugin_results
        .iter()
        .map(|r| {
            let score = r.score.unwrap_or(0.0);
            format!("{}={:.4}", r.name, score)
        })
        .collect::<Vec<_>>()
        .join("|");
    let plugin_penalties_raw = plugin_results
        .iter()
        .map(|r| {
            let pen = r.penalty.unwrap_or(0.0);
            format!("{}={:.4}", r.name, pen)
        })
        .collect::<Vec<_>>()
        .join("|");
    let plugin_metadata_raw = plugin_results
        .iter()
        .filter_map(|r| r.metadata.as_ref().map(|m| (r.name.as_str(), m)))
        .map(|(name, meta)| format!("{}={}", name, meta))
        .collect::<Vec<_>>()
        .join("|");
    let plugin_names = sanitize_csv_field(plugin_names_raw);
    let plugin_scores = sanitize_csv_field(plugin_scores_raw);
    let plugin_penalties = sanitize_csv_field(plugin_penalties_raw);
    let plugin_metadata = sanitize_csv_field(plugin_metadata_raw);
    let plugin_count = plugin_results.len();

    let aln = ctx.alignment_map.get(&m.gene_id);
    let hmmsum = ctx.hmmsum_map.get(&m.gene_id);
    let comp_entry = ctx.comp_map.get(&m.gene_id);
    let homology_score = comp_entry
        .map(|c| c.homology)
        .unwrap_or_else(|| compute_homology_score(summary));
    let intrinsic_score = comp_entry
        .map(|c| c.intrinsic)
        .unwrap_or_else(|| compute_intrinsic_score(&intrinsic));
    let taxonomy_score = if ctx.taxonomy_enabled {
        if let Some(val) = comp_entry.and_then(|c| c.taxonomy) {
            Some(val)
        } else {
            taxonomy_entry.and_then(|ev| {
                if ev.considered > 0 || ev.top_hit.is_some() {
                    Some(compute_taxonomy_score(Some(ev)))
                } else {
                    None
                }
            })
        }
    } else {
        None
    };
    let domains_arch_value = ctx.arch_map.get(&m.gene_id).copied();
    let domains_arch_score = domains_arch_value.unwrap_or(0.0);
    let (length_score, _len_z, len_ratio, len_class, len_min, len_max, len_in_range, len_panel_n) =
        ctx.len_map.get(&m.gene_id).cloned().unwrap_or((
            0.0,
            0.0,
            0.0,
            String::new(),
            0.0,
            0.0,
            false,
            0,
        ));
    let orphan_score = if ctx.orphan_analysis_enabled {
        comp_entry.map(|c| c.orphan).unwrap_or_else(|| {
            ctx.orphan_map
                .get(&m.gene_id)
                .map(|oa| oa.score)
                .unwrap_or(1.0)
        })
    } else {
        1.0
    };
    let subject_cov_score = comp_entry
        .map(|c| c.subject_cov)
        .unwrap_or_else(|| compute_subject_cov_score(summary));
    let subject_cov_penalty = compute_subject_cov_penalty(summary);
    let rnaseq_score = if ctx.rnaseq_enabled {
        comp_entry.map(|c| c.rnaseq).unwrap_or(0.0)
    } else {
        0.0
    };
    let rnaseq_metrics = ctx.rnaseq_map.get(&m.gene_id);
    let rnaseq_tpm = rnaseq_metrics
        .and_then(|r| r.tpm)
        .map_or(String::new(), |v| v.to_string());
    let rnaseq_num_reads = rnaseq_metrics
        .and_then(|r| r.num_reads)
        .map_or(String::new(), |v| v.to_string());
    let termini_score = comp_entry
        .and_then(|c| c.termini)
        .or_else(|| {
            aln.and_then(|a| {
                if a.mafft_enabled && a.sequences_aligned >= mafft::MIN_PANEL_FOR_TERMINI {
                    Some((a.start_concordance + a.end_concordance) / 2.0)
                } else {
                    None
                }
            })
        })
        .unwrap_or(0.0);
    let conserved_regions_score = comp_entry
        .and_then(|c| c.conserved_regions)
        .or_else(|| {
            aln.and_then(|a| {
                if a.mafft_enabled && a.sequences_aligned >= mafft::MIN_PANEL_FOR_CONSERVED_REGIONS
                {
                    Some(compute_conserved_regions_score(Some(a)))
                } else {
                    None
                }
            })
        })
        .unwrap_or(0.0);
    let block_conservation_score = aln
        .and_then(|a| {
            if a.mafft_enabled && !a.conserved_blocks.is_empty() {
                Some(scoring::compute_block_conservation_score(Some(a)))
            } else {
                None
            }
        })
        .unwrap_or(1.0);
    let genomic_score = comp_entry
        .map(|c| c.genomic)
        .unwrap_or_else(|| compute_genomic_score(genomic_metrics));
    let (base_final_score, classif_base) = ctx
        .scores_map
        .get(&m.gene_id)
        .cloned()
        .unwrap_or((0.0, "Low".to_string()));

    let final_score = (base_final_score - plugin_penalty).clamp(0.0, 1.0);
    let classif_final = if final_score >= ctx.th_high {
        "High"
    } else if final_score >= ctx.th_med {
        "Medium"
    } else {
        "Low"
    }
    .to_string();

    let raw_final_score = ctx
        .raw_scores_map
        .get(&m.gene_id)
        .copied()
        .unwrap_or(base_final_score);
    let fusion_split_flag = summary
        .map(|s| s.coverage_delta > ctx.cov_delta_thresh)
        .unwrap_or(false);
    let mut base_warnings: Vec<String> = Vec::new();
    if m.hits == 0 {
        base_warnings.push("No DIAMOND hits".to_string());
    }
    let mut warn_extra: Vec<String> = Vec::new();
    let mut warning_msgs = base_warnings.clone();
    if let Some(a) = aln {
        if a.missing_exon_run >= ctx.mafft_missing_exon_thresh {
            warn_extra.push("MissingExonPossible".into());
            warning_msgs.push("MissingExonPossible".into());
        }
        if a.retained_intron_run >= ctx.mafft_retained_intron_thresh {
            warn_extra.push("RetainedIntronPossible".into());
            warning_msgs.push("RetainedIntronPossible".into());
        }
        if !a.missing_blocks.is_empty() {
            let missing_count = a.missing_blocks.len();
            warn_extra.push(format!("MissingConservedBlocks(count={})", missing_count));
            warning_msgs.push(format!("MissingConservedBlocks(count={})", missing_count));
        }
        if !a.extra_blocks.is_empty() {
            let extra_count = a.extra_blocks.len();
            warn_extra.push(format!("ExtraConservedBlocks(count={})", extra_count));
            warning_msgs.push(format!("ExtraConservedBlocks(count={})", extra_count));
        }
    }
    let sv_obj = ctx.structvar_map.get(&m.gene_id);
    let structvar_multiplier = if m.hits > 0 {
        compute_structvar_multiplier(sv_obj)
    } else {
        1.0
    };
    if let Some(sv) = sv_obj {
        match sv.classification.as_str() {
            "FusionPossible" => {
                warn_extra.push("FusionPossible".into());
                warning_msgs.push("FusionPossible".into());
            }
            "SplitPossible" => {
                warn_extra.push("SplitPossible".into());
                warning_msgs.push("SplitPossible".into());
            }
            "InternalDuplicationPossible" => {
                warn_extra.push("InternalDuplicationPossible".into());
                warning_msgs.push("InternalDuplicationPossible".into());
            }
            _ => {}
        }
        for w in &sv.warnings {
            warn_extra.push(w.clone());
            warning_msgs.push(w.clone());
        }
    }
    let taxonomy_json = if ctx.taxonomy_enabled {
        if let Some(Some(ev)) = ctx.taxsum_map.get(&m.gene_id) {
            let mut obj = serde_json::Map::new();
            obj.insert("status".into(), serde_json::Value::String("enabled".into()));
            obj.insert(
                "detail".into(),
                serde_json::Value::String(ev.detail.to_string()),
            );
            obj.insert(
                "congruence_score".into(),
                serde_json::json!(ev.congruence_score),
            );
            obj.insert(
                "contamination_score".into(),
                serde_json::json!(ev.contamination_score),
            );
            obj.insert(
                "support_fraction".into(),
                serde_json::json!(ev.support_fraction),
            );
            obj.insert(
                "support_hits".into(),
                serde_json::Value::Number(serde_json::Number::from(ev.support as u64)),
            );
            obj.insert(
                "considered_hits".into(),
                serde_json::Value::Number(serde_json::Number::from(ev.considered as u64)),
            );
            obj.insert(
                "consensus_depth".into(),
                serde_json::Value::Number(serde_json::Number::from(ev.consensus_depth as u64)),
            );
            if let Some(rank) = &ev.consensus_rank {
                obj.insert(
                    "consensus_rank".into(),
                    serde_json::Value::String(rank.clone()),
                );
            }
            if let Some(top) = &ev.top_hit {
                obj.insert(
                    "taxid".into(),
                    serde_json::Value::Number(serde_json::Number::from(top.taxid)),
                );
                if let Some(name) = &top.name {
                    obj.insert("name".into(), serde_json::Value::String(name.clone()));
                }
                obj.insert("lineage".into(), serde_json::json!(top.lineage));
                obj.insert(
                    "top_hit".into(),
                    serde_json::json!({
                        "taxid": top.taxid,
                        "name": top.name,
                        "lineage": top.lineage,
                    }),
                );
            }
            if let Some(consensus) = &ev.consensus {
                obj.insert(
                    "consensus_taxid".into(),
                    serde_json::Value::Number(serde_json::Number::from(consensus.taxid)),
                );
                if let Some(name) = &consensus.name {
                    obj.insert(
                        "consensus_name".into(),
                        serde_json::Value::String(name.clone()),
                    );
                }
                obj.insert(
                    "consensus_lineage".into(),
                    serde_json::json!(consensus.lineage.clone()),
                );
            }
            serde_json::Value::Object(obj)
        } else {
            serde_json::json!({"status":"enabled","detail":"NoResolver"})
        }
    } else {
        serde_json::json!({"status":"disabled"})
    };
    let domains = hmmsum.map(|d| {
        let score = if let Some(ev) = d.top_evalue {
            let le = if ev > 0.0 { -ev.log10() } else { 100.0 };
            (le / 20.0).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let arch_score = ctx.arch_map.get(&m.gene_id).copied().unwrap_or(0.0);
        serde_json::json!({
            "hits_count": d.hits_count,
            "top_accession": d.top_accession,
            "top_evalue": d.top_evalue,
            "domains_score": score,
            "domains_arch_score": arch_score,
            "hits": d.hits.iter().map(|h| serde_json::json!({
                "target_name": h.target_name,
                "accession": h.accession,
                "evalue": h.evalue,
                "score": h.score,
                "bias": h.bias
            })).collect::<Vec<_>>()
        })
    });
    let length_block = ctx
        .len_map
        .get(&m.gene_id)
        .map(|(s, z, r, c, mn, mx, in_rng, n)| {
            serde_json::json!({
                "length_score": s,
                "length_z": z,
                "length_ratio": r,
                "length_class": c,
                "expected_len_min": mn,
                "expected_len_max": mx,
                "in_expected_range": in_rng,
                "panel_n": n,
            })
        });
    let orphan_entry = if ctx.orphan_analysis_enabled {
        ctx.orphan_map.get(&m.gene_id)
    } else {
        None
    };
    let divergence_score = comp_entry
        .and_then(|c| c.divergence)
        .unwrap_or_else(|| scoring::compute_divergence_score(aln));
    let record = serde_json::json!({
        "gene_id": m.gene_id,
        "taxonomy": taxonomy_json,
        "score_components": {
            "taxonomy": taxonomy_score,
            "homology": homology_score,
            "intrinsic": intrinsic_score,
            "domains": domains_arch_score,
            "domains_strength": compute_domains_strength_score(hmmsum),
            "subject_coverage": subject_cov_score,
            "subject_cov_penalty": subject_cov_penalty,
            "length": length_score,
            "orphan": orphan_score,
            "conserved_regions": conserved_regions_score,
            "block_conservation": block_conservation_score,
            "structvar_multiplier": structvar_multiplier,
            "termini": termini_score,
            "genomic": genomic_score,
            "rnaseq": rnaseq_score,
            "divergence": divergence_score
        },
        "final_score_raw": raw_final_score,
        "base_score": base_final_score,
        "final_score": final_score,
        "classification_base": classif_base.clone(),
        "classification_final": classif_final.clone(),
        "homology": {
            "hits_count": m.hits,
            "top_hit": summary.and_then(|s| s.top_sseqid.clone()),
            "top_bitscore": summary.map(|s| s.top_bitscore),
            "top_evalue": summary.as_ref().map(|s| s.top_evalue.clone()),
            "top_qcov": summary.map(|s| s.top_qcov),
            "top_scov": summary.map(|s| s.top_scov),
            "bitscore_density": summary.map(|s| if s.top_len > 0 {
                s.top_bitscore / s.top_len as f64
            } else { 0.0 }),
            "coverage_delta": summary.map(|s| s.coverage_delta),
            "coverage_ratio": summary.map(|s| s.coverage_ratio),
            "fusion_split_flag": fusion_split_flag,
        },
        "panel_provenance": {
            "swissprot": panel_prov.swissprot,
            "refprot": panel_prov.refprot,
            "cluster": panel_prov.cluster,
        },
        "intrinsic": {
            "ambiguous_fraction": intrinsic.ambiguous_fraction,
            "max_homopolymer": intrinsic.max_homopolymer,
            "low_complexity_fraction": intrinsic.low_complexity_fraction,
            "low_complexity_windows": intrinsic.low_complexity_windows,
            "orf_diagnostics": {
                "start_methionine": intrinsic.start_methionine,
                "alt_start_pos": intrinsic.alt_start_pos,
                "internal_stop_count": intrinsic.internal_stop_count,
                "terminal_stop": intrinsic.terminal_stop,
                "orf_start_score": intrinsic.orf_start_score,
                "orf_score": scoring::compute_orf_score(&intrinsic),
            },
        },
        "alignment": aln.as_ref().map(|a| serde_json::json!({
            "mafft_enabled": a.mafft_enabled,
            "strategy_used": a.strategy_used,
            "conserved_fraction": a.conserved_fraction,
            "pairwise_identity": a.pairwise_identity,
            "sequences_aligned": a.sequences_aligned,
            "query_gap_fraction": a.query_gap_fraction,
            "gap_run_count": a.gap_run_count,
            "max_gap_run": a.max_gap_run,
            "motif_mismatch_fraction": a.motif_mismatch_fraction,
            "start_concordance": a.start_concordance,
            "start_class": a.start_class,
            "end_concordance": a.end_concordance,
            "end_class": a.end_class,
            "missing_exon_run": a.missing_exon_run,
            "retained_intron_run": a.retained_intron_run,
        })),
        "domains": domains,
        "length": length_block,
        "orphan_analysis": orphan_entry.map(|oa| serde_json::json!({
            "status": oa.status.as_str(),
            "score": oa.score,
            "details": oa.details.iter().map(|d| serde_json::json!({
                "accession": d.accession,
                "domain_index": d.domain_index,
                "total_domains": d.total_domains,
                "completeness": d.completeness,
                "hmm_from": d.hmm_from,
                "hmm_to": d.hmm_to,
                "hmm_len": d.hmm_len,
            })).collect::<Vec<_>>()
        })),
        "structvar": sv_obj.map(|sv| serde_json::json!({
            "classification": sv.classification,
            "fusion_possible": sv.fusion_possible,
            "split_possible": sv.split_possible,
            "duplication_possible": sv.duplication_possible,
            "spans": sv.spans,
            "fusion_subjects": sv.fusion_subjects,
            "fusion_gap": sv.fusion_gap,
            "fusion_left_len": sv.fusion_left_len,
            "fusion_right_len": sv.fusion_right_len,
            "fusion_cover_fracs": sv.fusion_cover_fracs,
            "subjects": sv.subjects,
            "subject_warnings": sv.warnings,
        })),
        "genomic": genomic_metrics.map(|g| serde_json::json!({
            "introns_total": g.introns_total,
            "splice_canonical": g.splice_canonical,
            "splice_noncanonical": g.splice_major_noncan + g.splice_minor,
            "splice_weird": g.splice_weird,
            "intron_len_min": g.intron_len_min,
            "intron_len_max": g.intron_len_max,
            "intron_len_avg": g.intron_len_avg,
            "genomic_score": genomic_score,
        })),
        "rnaseq": ctx.rnaseq_map.get(&m.gene_id).map(|r| serde_json::json!({
            "tpm": r.tpm,
            "num_reads": r.num_reads,
            "expression_score": r.expression_score,
        })),
        "warnings": if warn_extra.is_empty() { base_warnings.clone() } else { warn_extra.clone() },
        "plugins": serde_json::json!({
            "total_penalty": plugin_penalty,
            "results": plugin_results.clone(),
        }),
    });

    let json_line = serde_json::to_string(&record).map_err(|e| e.to_string())?;

    let (top_hit, top_bitscore, top_evalue, top_qcov, top_scov, bsd, cov_delta, cov_ratio) =
        if let Some(s) = summary {
            (
                s.top_sseqid.clone().unwrap_or_default(),
                s.top_bitscore,
                s.top_evalue.clone(),
                s.top_qcov,
                s.top_scov,
                if s.top_len > 0 {
                    s.top_bitscore / s.top_len as f64
                } else {
                    0.0
                },
                s.coverage_delta,
                s.coverage_ratio,
            )
        } else {
            (String::new(), 0.0, String::new(), 0.0, 0.0, 0.0, 0.0, 0.0)
        };
    let subject_cov_score_str = format!("{:.3}", subject_cov_score);
    let subject_cov_penalty_str = format!("{:.3}", subject_cov_penalty);

    let domains_score_field = if let Some(d) = hmmsum {
        let score = if let Some(ev) = d.top_evalue {
            let le = if ev > 0.0 { -ev.log10() } else { 100.0 };
            (le / 20.0).clamp(0.0, 1.0)
        } else {
            0.0
        };
        format!("{:.4}", score)
    } else {
        String::new()
    };
    let domains_arch_field = domains_arch_value
        .map(|v| format!("{:.4}", v))
        .unwrap_or_default();
    let orphan_status_str = orphan_entry
        .map(|oa| oa.status.as_str().to_string())
        .unwrap_or_default();
    let orphan_score_field = if ctx.orphan_analysis_enabled {
        orphan_entry
            .map(|oa| format!("{:.4}", oa.score))
            .unwrap_or_default()
    } else {
        String::new()
    };
    let warnings_field = if warning_msgs.is_empty() {
        String::new()
    } else {
        warning_msgs.join(";")
    };
    let (
        mafft_enabled,
        conserved,
        pid,
        panel_pid,
        div_ratio,
        seqs_aln,
        qgap,
        gap_runs,
        max_gap,
        missing_run,
        intron_run,
        start_conc,
        start_class,
        end_conc,
        end_class,
    ) = if let Some(a) = aln {
        (
            a.mafft_enabled,
            a.conserved_fraction,
            a.pairwise_identity,
            a.panel_pairwise_identity,
            a.divergence_ratio,
            a.sequences_aligned,
            a.query_gap_fraction,
            a.gap_run_count,
            a.max_gap_run,
            a.missing_exon_run,
            a.retained_intron_run,
            a.start_concordance,
            a.start_class.clone(),
            a.end_concordance,
            a.end_class.clone(),
        )
    } else {
        (
            false,
            0.0,
            0.0,
            0.0,
            0.0,
            0,
            0.0,
            0,
            0,
            0,
            0,
            0.0,
            String::new(),
            0.0,
            String::new(),
        )
    };

    let (
        structvar_class,
        structvar_gap,
        structvar_left_len,
        structvar_right_len,
        structvar_cov_left,
        structvar_cov_right,
    ) = if let Some(sv) = sv_obj {
        (
            sv.classification.clone(),
            sv.fusion_gap.map(|g| g.to_string()).unwrap_or_default(),
            sv.fusion_left_len
                .map(|g| g.to_string())
                .unwrap_or_default(),
            sv.fusion_right_len
                .map(|g| g.to_string())
                .unwrap_or_default(),
            sv.fusion_cover_fracs
                .map(|(a, _)| format!("{:.3}", a))
                .unwrap_or_default(),
            sv.fusion_cover_fracs
                .map(|(_, b)| format!("{:.3}", b))
                .unwrap_or_default(),
        )
    } else {
        (
            String::new(),
            String::new(),
            String::new(),
            String::new(),
            String::new(),
            String::new(),
        )
    };
    let taxonomy_score_value = if ctx.taxonomy_enabled {
        taxonomy_entry
            .and_then(|ev| {
                if ev.considered > 0 || ev.top_hit.is_some() {
                    Some(compute_taxonomy_score(Some(ev)))
                } else {
                    None
                }
            })
            .unwrap_or(0.0)
    } else {
        0.0
    };
    let top_evalue_val = summary
        .and_then(|s| s.top_evalue.parse::<f64>().ok())
        .unwrap_or(0.0);
    let domains_strength_score = compute_domains_strength_score(hmmsum);
    let domains_score = hmmsum.map(|_| domains_strength_score);

    let (
        taxonomy_contamination_field,
        taxonomy_support_field,
        taxonomy_considered_field,
        taxonomy_support_frac_field,
        taxonomy_status_str,
        consensus_label,
        taxonomy_rank_field,
    ) = if ctx.taxonomy_enabled {
        if let Some(ev) = taxonomy_entry {
            let label = ev
                .consensus
                .as_ref()
                .map(|c| {
                    if let Some(name) = &c.name {
                        format!("{} ({})", name, c.taxid)
                    } else {
                        c.taxid.to_string()
                    }
                })
                .unwrap_or_default();
            (
                format!("{:.4}", ev.contamination_score),
                ev.support.to_string(),
                ev.considered.to_string(),
                format!("{:.4}", ev.support_fraction),
                ev.detail.to_string(),
                label,
                ev.consensus_rank.clone().unwrap_or_default(),
            )
        } else {
            (
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                "NoResolver".to_string(),
                String::new(),
                String::new(),
            )
        }
    } else {
        (
            String::new(),
            String::new(),
            String::new(),
            String::new(),
            "disabled".to_string(),
            String::new(),
            String::new(),
        )
    };
    let (taxonomy_domain, taxonomy_genus) = if ctx.taxonomy_enabled {
        if let (Some(ev), Some(resolver)) = (taxonomy_entry, ctx.taxonomy_resolver.as_ref()) {
            if let Some(consensus) = &ev.consensus {
                let domain = resolver
                    .ancestor_at_rank(consensus.taxid, "domain")
                    .and_then(|tid| resolver.name_of(tid).map(|s| s.to_string()))
                    .unwrap_or_else(|| "Unknown".to_string());
                let genus = resolver
                    .ancestor_at_rank(consensus.taxid, "genus")
                    .and_then(|tid| resolver.name_of(tid).map(|s| s.to_string()))
                    .unwrap_or_default();
                (domain, genus)
            } else {
                (String::new(), String::new())
            }
        } else {
            (String::new(), String::new())
        }
    } else {
        (String::new(), String::new())
    };

    let len_ratio_str = if len_ratio > 0.0 {
        format!("{:.3}", len_ratio)
    } else {
        String::new()
    };
    let len_expected_min_str = if len_panel_n > 0 && len_min > 0.0 {
        format!("{:.0}", len_min)
    } else {
        String::new()
    };
    let len_expected_max_str = if len_panel_n > 0 && len_max > 0.0 {
        format!("{:.0}", len_max)
    } else {
        String::new()
    };
    let len_in_range_str = if len_panel_n > 0 {
        len_in_range.to_string()
    } else {
        String::new()
    };
    let len_panel_n_str = if len_panel_n > 0 {
        len_panel_n.to_string()
    } else {
        String::new()
    };
    let panel_swissprot_field = panel_prov.swissprot.to_string();
    let panel_refprot_field = panel_prov.refprot.to_string();
    let panel_cluster_field = panel_prov.cluster.to_string();
    let orphan_component_str = if ctx.orphan_analysis_enabled {
        if orphan_score_field.is_empty() {
            let orphan_component = comp_entry.map(|c| c.orphan).unwrap_or(orphan_score);
            format!("{:.4}", orphan_component)
        } else {
            orphan_score_field.clone()
        }
    } else {
        String::new()
    };
    let genomic_score_str = if ctx.genomic_map.is_some() {
        format!("{:.4}", genomic_score)
    } else {
        String::new()
    };

    let rnaseq_score_str = if ctx.rnaseq_enabled {
        format!("{:.4}", rnaseq_score)
    } else {
        String::new()
    };

    let warnings_field_clone = warnings_field.clone();
    let csv_line = if ctx.csv_verbose {
        let taxonomy_component_str = if ctx.taxonomy_enabled {
            format!("{:.4}", taxonomy_score.unwrap_or(0.0))
        } else {
            String::new()
        };
        vec![
            m.gene_id.clone(),
            m.hits.to_string(),
            panel_swissprot_field.clone(),
            panel_refprot_field.clone(),
            panel_cluster_field.clone(),
            top_hit.clone(),
            format!("{:.3}", top_bitscore),
            top_evalue.clone(),
            format!("{:.3}", top_qcov),
            format!("{:.3}", top_scov),
            format!("{:.3}", bsd),
            format!("{:.3}", cov_delta),
            format!("{:.3}", cov_ratio),
            subject_cov_score_str.clone(),
            subject_cov_penalty_str.clone(),
            fusion_split_flag.to_string(),
            format!("{:.4}", structvar_multiplier),
            format!("{:.3}", final_score),
            classif_base.clone(),
            classif_final.clone(),
            format!("{:.4}", homology_score),
            format!("{:.4}", intrinsic_score),
            genomic_score_str.clone(),
            taxonomy_component_str,
            rnaseq_score_str.clone(),
            rnaseq_tpm.clone(),
            rnaseq_num_reads.clone(),
            domains_score_field.clone(),
            domains_arch_field.clone(),
            orphan_component_str,
            format!("{:.4}", length_score),
            len_ratio_str,
            len_class.clone(),
            len_expected_min_str,
            len_expected_max_str,
            len_in_range_str,
            len_panel_n_str,
            format!("{:.4}", conserved_regions_score),
            format!("{:.4}", termini_score),
            format!("{:.4}", divergence_score),
            intrinsic.start_methionine.to_string(),
            intrinsic.internal_stop_count.to_string(),
            intrinsic.terminal_stop.to_string(),
            format!("{:.4}", scoring::compute_orf_score(&intrinsic)),
            mafft_enabled.to_string(),
            format!("{:.3}", conserved),
            format!("{:.3}", pid),
            format!("{:.3}", panel_pid),
            format!("{:.3}", div_ratio),
            seqs_aln.to_string(),
            format!("{:.3}", qgap),
            gap_runs.to_string(),
            max_gap.to_string(),
            missing_run.to_string(),
            intron_run.to_string(),
            format!("{:.3}", start_conc),
            start_class.clone(),
            format!("{:.3}", end_conc),
            end_class.clone(),
            structvar_class.clone(),
            structvar_gap.clone(),
            structvar_left_len.clone(),
            structvar_right_len.clone(),
            structvar_cov_left.clone(),
            structvar_cov_right.clone(),
            orphan_status_str.clone(),
            taxonomy_contamination_field.clone(),
            taxonomy_support_field.clone(),
            taxonomy_considered_field.clone(),
            taxonomy_support_frac_field.clone(),
            consensus_label.clone(),
            taxonomy_rank_field.clone(),
            taxonomy_status_str.clone(),
            format!("{:.4}", plugin_penalty),
            plugin_names.clone(),
            plugin_scores.clone(),
            plugin_penalties.clone(),
            plugin_metadata.clone(),
            warnings_field,
        ]
        .join(",")
    } else {
        vec![
            m.gene_id.clone(),
            m.hits.to_string(),
            panel_swissprot_field.clone(),
            panel_refprot_field.clone(),
            panel_cluster_field.clone(),
            top_hit.clone(),
            format!("{:.3}", top_bitscore),
            top_evalue.clone(),
            format!("{:.3}", top_qcov),
            format!("{:.3}", top_scov),
            format!("{:.3}", bsd),
            format!("{:.3}", cov_delta),
            format!("{:.3}", cov_ratio),
            subject_cov_score_str.clone(),
            subject_cov_penalty_str.clone(),
            fusion_split_flag.to_string(),
            format!("{:.4}", structvar_multiplier),
            format!("{:.3}", final_score),
            classif_base.clone(),
            classif_final.clone(),
            mafft_enabled.to_string(),
            format!("{:.3}", conserved),
            format!("{:.3}", pid),
            format!("{:.3}", panel_pid),
            format!("{:.3}", div_ratio),
            seqs_aln.to_string(),
            format!("{:.3}", qgap),
            gap_runs.to_string(),
            max_gap.to_string(),
            domains_score_field.clone(),
            domains_arch_field.clone(),
            orphan_score_field.clone(),
            structvar_class.clone(),
            structvar_gap.clone(),
            structvar_left_len.clone(),
            structvar_right_len.clone(),
            structvar_cov_left.clone(),
            structvar_cov_right.clone(),
            orphan_status_str.clone(),
            genomic_score_str.clone(),
            format!("{:.4}", taxonomy_score_value),
            taxonomy_contamination_field.clone(),
            taxonomy_support_field.clone(),
            taxonomy_considered_field.clone(),
            taxonomy_support_frac_field.clone(),
            consensus_label.clone(),
            taxonomy_rank_field.clone(),
            taxonomy_status_str.clone(),
            format!("{:.4}", plugin_penalty),
            plugin_names.clone(),
            plugin_scores.clone(),
            plugin_penalties.clone(),
            plugin_metadata.clone(),
            warnings_field_clone.clone(),
        ]
        .join(",")
    };

    let card = ScoreCard {
        gene_id: m.gene_id.clone(),
        hits_count: m.hits,
        panel_swissprot: panel_prov.swissprot,
        panel_refprot: panel_prov.refprot,
        panel_cluster: panel_prov.cluster,
        top_hit: top_hit.clone(),
        top_bitscore,
        top_evalue: top_evalue_val,
        top_qcov,
        top_scov,
        bitscore_density: bsd,
        coverage_delta: cov_delta,
        coverage_ratio: cov_ratio,
        subject_cov_score,
        subject_cov_penalty,
        fusion_split: fusion_split_flag,
        structvar_multiplier,
        final_score,
        classification_base: classif_base.clone(),
        classification_final: classif_final.clone(),
        homology_score,
        intrinsic_score,
        taxonomy_score,
        domains_score,
        domains_arch_score,
        orphan_domain_score: comp_entry.map(|c| c.orphan).unwrap_or(orphan_score),
        length_score,
        length_ratio: len_ratio,
        length_class: len_class.clone(),
        expected_len_min: if len_panel_n > 0 { Some(len_min) } else { None },
        expected_len_max: if len_panel_n > 0 { Some(len_max) } else { None },
        length_in_expected_range: if len_panel_n > 0 {
            Some(len_in_range)
        } else {
            None
        },
        length_panel_n: if len_panel_n > 0 {
            Some(len_panel_n)
        } else {
            None
        },
        conserved_regions_score,
        termini_score,
        divergence_score,
        mafft_enabled,
        conserved_fraction: conserved,
        pairwise_identity: pid,
        panel_pairwise_identity: panel_pid,
        divergence_ratio: div_ratio,
        sequences_aligned: seqs_aln,
        query_gap_fraction: qgap,
        gap_run_count: gap_runs,
        max_gap_run: max_gap,
        missing_exon_run: missing_run,
        retained_intron_run: intron_run,
        start_concordance: start_conc,
        start_class: start_class.clone(),
        end_concordance: end_conc,
        end_class: end_class.clone(),
        structvar_class: structvar_class.clone(),
        structvar_gap: sv_obj.and_then(|sv| sv.fusion_gap),
        orphan_status: orphan_status_str.clone(),
        taxonomy_contamination: if ctx.taxonomy_enabled {
            taxonomy_entry.map(|ev| ev.contamination_score)
        } else {
            None
        },
        taxonomy_support: if let Some(ev) = taxonomy_entry {
            ev.support
        } else {
            0
        },
        taxonomy_considered: if let Some(ev) = taxonomy_entry {
            ev.considered
        } else {
            0
        },
        taxonomy_support_frac: if let Some(ev) = taxonomy_entry {
            ev.support_fraction
        } else {
            0.0
        },
        consensus_taxon: consensus_label.clone(),
        taxonomy_rank: taxonomy_rank_field.clone(),
        taxonomy_status: taxonomy_status_str.clone(),
        taxonomy_domain,
        taxonomy_genus,
        genomic_introns: genomic_metrics.map(|g| g.introns_total),
        genomic_splice_canonical: genomic_metrics.map(|g| g.splice_canonical),
        genomic_splice_noncanonical: genomic_metrics
            .map(|g| g.splice_major_noncan + g.splice_minor),
        genomic_splice_weird: genomic_metrics.map(|g| g.splice_weird),
        genomic_score,
        rnaseq_score,
        rnaseq_tpm,
        rnaseq_num_reads,
        plugin_penalty,
        plugin_count,
        plugin_names,
        plugin_scores,
        plugin_penalties,
        plugin_metadata,
        warnings: warnings_field_clone,
    };

    Ok(RenderedRecord {
        index,
        json_line,
        csv_line,
        card: Some(card),
    })
}
