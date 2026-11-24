mod algorithms;
mod catalog;
mod data;
mod gui;
mod math;
mod metrics;
mod pdf_report;

use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use algorithms::{run_selected_algorithms, QualityPreset, RunConfig};
use catalog::ALGORITHM_CATALOG;
use data::{prepare_dataset, PreprocessConfig, RawCsvDataset};
use pdf_report::{generate_pdf_report, PdfReportOptions};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();

    // Optional CLI batch mode: cluster_overview --csv input.csv --pdf report.pdf [--k 3] [--max-samples 240] [--compact]
    if args.iter().any(|a| a == "--csv" || a == "--help" || a == "-h") {
        return run_cli(&args);
    }

    // Default: Launch native Rust desktop GUI (eframe/egui — zero HTML/WebView)
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("ClusterOverview — Klyngeanalyse (200 algoritmer)  |  © Robert Larsen")
            .with_inner_size([1240.0, 780.0])
            .with_min_inner_size([840.0, 560.0])
            .with_drag_and_drop(true),
        ..Default::default()
    };

    eframe::run_native(
        "ClusterOverview",
        native_options,
        Box::new(|cc| Ok(Box::new(gui::ClusterOverviewApp::new(cc)))),
    )?;
    Ok(())
}

fn run_cli(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    if args.iter().any(|a| a == "--help" || a == "-h") {
        println!("ClusterOverview v1.0.0 — 200 klyngealgoritmer i ren Rust");
        println!("Copyright (c) 2026 Robert Larsen (MIT License)");
        println!("Bruk:");
        println!("  cluster_overview                                   (Starter native GUI)");
        println!("  cluster_overview --csv data.csv --pdf rapport.pdf [--k 3] [--max-samples 240] [--compact]");
        return Ok(());
    }

    let mut csv_path: Option<PathBuf> = None;
    let mut pdf_path: Option<PathBuf> = None;
    let mut target_k: Option<usize> = None;
    let mut max_samples = 240usize;
    let mut compact = false;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--csv" if i + 1 < args.len() => {
                csv_path = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--pdf" if i + 1 < args.len() => {
                pdf_path = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--k" if i + 1 < args.len() => {
                target_k = args[i + 1].parse().ok();
                i += 2;
            }
            "--max-samples" if i + 1 < args.len() => {
                if let Ok(v) = args[i + 1].parse() {
                    max_samples = v;
                }
                i += 2;
            }
            "--compact" => {
                compact = true;
                i += 1;
            }
            _ => i += 1,
        }
    }

    let csv_p = csv_path.ok_or("Mangler --csv <sti>")?;
    let pdf_p = pdf_path.unwrap_or_else(|| PathBuf::from("clusterOverview_rapport.pdf"));

    println!("Leser CSV: {}", csv_p.display());
    let raw = RawCsvDataset::from_file(&csv_p)?;
    let mut pcfg = PreprocessConfig::default_for(&raw);
    pcfg.max_samples = max_samples;
    let prepared = Arc::new(prepare_dataset(&raw, &pcfg)?);

    let k = target_k.unwrap_or(prepared.suggested_k);
    let rcfg = RunConfig {
        target_k: k,
        seed: 42,
        quality: QualityPreset::Standard,
    };

    let all_ids: Vec<usize> = (1..=ALGORITHM_CATALOG.len()).collect();
    println!(
        "Kjører alle {} klyngealgoritmer (n={}, d={}, k={})...",
        all_ids.len(),
        prepared.n_samples,
        prepared.n_features,
        k
    );
    let t0 = Instant::now();
    let results = run_selected_algorithms(&prepared, &all_ids, &rcfg, None);
    let elapsed = t0.elapsed().as_secs_f64();

    let errors = results.iter().filter(|r| r.error.is_some()).count();
    let valid_multi = results.iter().filter(|r| r.metrics.n_clusters >= 2).count();
    println!(
        "Fullført {} algoritmer på {:.2} sek! (Feil: {}, Flere klynger: {}/{})",
        results.len(),
        elapsed,
        errors,
        valid_multi,
        results.len()
    );

    let mut opts = PdfReportOptions::default();
    opts.subtitle = format!("Datasett: {} ({} rader, {} kolonner)", raw.display_name, raw.rows.len(), raw.columns.len());
    opts.include_scatter_cards = !compact;

    let pages = generate_pdf_report(&pdf_p, &raw, &prepared, &results, &opts)?;
    println!("PDF-rapport skrevet til '{}' ({} sider).", pdf_p.display(), pages);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verify_catalog_has_all_200_algorithms_ordered() {
        assert_eq!(ALGORITHM_CATALOG.len(), 200);
        for (idx, spec) in ALGORITHM_CATALOG.iter().enumerate() {
            assert_eq!(spec.id, idx + 1, "Mismatch at index {}", idx);
            assert!(!spec.name.is_empty());
            assert!(!spec.description_no.is_empty());
        }
    }

    #[test]
    fn verify_all_200_algorithms_run_cleanly_on_iris_and_generate_pdf() {
        let raw = RawCsvDataset::built_in_demo_iris();
        let pcfg = PreprocessConfig::default_for(&raw);
        let prep = Arc::new(prepare_dataset(&raw, &pcfg).expect("prepare iris"));
        let rcfg = RunConfig {
            target_k: 3,
            seed: 42,
            quality: QualityPreset::Fast,
        };
        let ids: Vec<usize> = (1..=200).collect();
        let results = run_selected_algorithms(&prep, &ids, &rcfg, None);
        assert_eq!(results.len(), 200);
        for r in &results {
            assert!(r.error.is_none(), "Algorithm #{} ({}) failed: {:?}", r.id, r.name, r.error);
            assert_eq!(r.labels.len(), prep.n_samples);
        }

        let tmp_pdf = std::env::temp_dir().join("cluster_overview_test_report.pdf");
        let pages = generate_pdf_report(&tmp_pdf, &raw, &prep, &results, &PdfReportOptions::default())
            .expect("generate pdf");
        assert!(pages >= 35, "Expected rich multi-page PDF, got {} pages", pages);
        let meta = std::fs::metadata(&tmp_pdf).expect("pdf exists");
        assert!(meta.len() > 20_000, "PDF too small: {} bytes", meta.len());
        let _ = std::fs::remove_file(&tmp_pdf);
    }
}
