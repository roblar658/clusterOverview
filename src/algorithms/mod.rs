pub mod cat1_centroid;
pub mod cat2_density;
pub mod cat3_hierarchical;
pub mod cat4_distribution;
pub mod cat5_graph;
pub mod cat6_grid_subspace;
pub mod cat7_fuzzy_domain;
pub mod cat8_neural_stream;

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::Instant;

use rayon::prelude::*;

use crate::catalog::{dispatch_algorithm, ALGORITHM_CATALOG};
use crate::data::PreparedData;
use crate::math::relabel_contiguous;
use crate::metrics::{evaluate_clustering, ClusterMetrics};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Category {
    Centroid = 1,
    Density = 2,
    Hierarchical = 3,
    Distribution = 4,
    GraphSpectral = 5,
    GridSubspace = 6,
    FuzzyDomain = 7,
    NeuralStream = 8,
}

impl Category {
    pub const ALL: [Category; 8] = [
        Category::Centroid,
        Category::Density,
        Category::Hierarchical,
        Category::Distribution,
        Category::GraphSpectral,
        Category::GridSubspace,
        Category::FuzzyDomain,
        Category::NeuralStream,
    ];

    pub fn number(&self) -> usize {
        *self as usize
    }

    pub fn title_no(&self) -> &'static str {
        match self {
            Category::Centroid => "1. Sentroide-, medoid- og partisjoneringsbaserte (1–30)",
            Category::Density => "2. Tetthets- og modusbaserte algoritmer (31–55)",
            Category::Hierarchical => "3. Hierarkiske og trebaserte algoritmer (56–80)",
            Category::Distribution => "4. Fordelings-, statistikk- og emnemodeller (81–102)",
            Category::GraphSpectral => "5. Graf-, spektral- og nettverksalgoritmer (103–132)",
            Category::GridSubspace => "6. Rutenett-, underroms- og projeksjonsbaserte (133–150)",
            Category::FuzzyDomain => "7. Fuzzy, biclustering, tidsserie- og domenespesifikke (151–175)",
            Category::NeuralStream => "8. Nevrale, dype, strømmebaserte og evolusjonære (176–200)",
        }
    }

    pub fn short_name_no(&self) -> &'static str {
        match self {
            Category::Centroid => "Sentroide/Partisjon",
            Category::Density => "Tetthet/Modus",
            Category::Hierarchical => "Hierarkisk/Tre",
            Category::Distribution => "Fordeling/Emne",
            Category::GraphSpectral => "Graf/Spektral",
            Category::GridSubspace => "Rutenett/Underrom",
            Category::FuzzyDomain => "Fuzzy/Domene",
            Category::NeuralStream => "Nevral/Strøm/Evo",
        }
    }

    pub fn accent_rgb(&self) -> [u8; 3] {
        match self {
            Category::Centroid => [37, 99, 235],      // Royal Blue
            Category::Density => [13, 148, 136],      // Teal
            Category::Hierarchical => [124, 58, 237], // Violet
            Category::Distribution => [217, 119, 6],  // Amber
            Category::GraphSpectral => [219, 39, 119],// Rose
            Category::GridSubspace => [5, 150, 105],  // Emerald
            Category::FuzzyDomain => [234, 88, 12],   // Orange
            Category::NeuralStream => [79, 70, 229],  // Indigo
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct AlgorithmSpec {
    pub id: usize,
    pub category: Category,
    pub name: &'static str,
    pub description_no: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QualityPreset {
    Fast,
    Standard,
    Thorough,
}

impl QualityPreset {
    pub fn label_no(&self) -> &'static str {
        match self {
            QualityPreset::Fast => "Rask",
            QualityPreset::Standard => "Standard",
            QualityPreset::Thorough => "Grundig",
        }
    }

    pub fn iter_scale(&self) -> f64 {
        match self {
            QualityPreset::Fast => 0.6,
            QualityPreset::Standard => 1.0,
            QualityPreset::Thorough => 1.6,
        }
    }
}

#[derive(Debug, Clone)]
pub struct RunConfig {
    pub target_k: usize,
    pub seed: u64,
    pub quality: QualityPreset,
}

impl RunConfig {
    pub fn max_iters(&self, base: usize) -> usize {
        ((base as f64) * self.quality.iter_scale()).round().max(5.0) as usize
    }

    pub fn effective_k(&self, n: usize) -> usize {
        self.target_k.clamp(2, n.saturating_sub(1).max(2))
    }
}

#[derive(Debug, Clone)]
pub struct AlgorithmOutput {
    pub labels: Vec<i32>,
    pub params_summary: String,
    pub note: String,
}

impl AlgorithmOutput {
    pub fn new(labels: Vec<i32>, params_summary: impl Into<String>, note: impl Into<String>) -> Self {
        Self {
            labels: relabel_contiguous(&labels),
            params_summary: params_summary.into(),
            note: note.into(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct AlgorithmRunResult {
    pub id: usize,
    pub category: Category,
    pub name: String,
    pub description_no: String,
    pub labels: Vec<i32>,
    pub metrics: ClusterMetrics,
    pub duration_ms: f64,
    pub params_summary: String,
    pub note: String,
    pub error: Option<String>,
}

pub struct RunProgress {
    pub total: AtomicUsize,
    pub completed: AtomicUsize,
    pub cancelled: AtomicBool,
}

impl RunProgress {
    pub fn new(total: usize) -> Self {
        Self {
            total: AtomicUsize::new(total),
            completed: AtomicUsize::new(0),
            cancelled: AtomicBool::new(false),
        }
    }
}

pub fn run_selected_algorithms(
    data: &Arc<PreparedData>,
    selected_ids: &[usize],
    cfg: &RunConfig,
    progress: Option<&Arc<RunProgress>>,
) -> Vec<AlgorithmRunResult> {
    let specs: Vec<&AlgorithmSpec> = selected_ids
        .iter()
        .filter_map(|&id| ALGORITHM_CATALOG.get(id.saturating_sub(1)))
        .collect();

    if let Some(p) = progress {
        p.total.store(specs.len(), Ordering::SeqCst);
        p.completed.store(0, Ordering::SeqCst);
    }

    let mut results: Vec<AlgorithmRunResult> = specs
        .par_iter()
        .filter_map(|spec| {
            if let Some(p) = progress {
                if p.cancelled.load(Ordering::Relaxed) {
                    return None;
                }
            }

            let t0 = Instant::now();
            let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                dispatch_algorithm(spec.id, data, cfg)
            }));
            let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;

            let res = match outcome {
                Ok(mut out) => {
                    if out.labels.len() != data.n_samples {
                        out.labels.resize(data.n_samples, -1);
                    }
                    let labels = relabel_contiguous(&out.labels);
                    let metrics = evaluate_clustering(data, &labels);
                    AlgorithmRunResult {
                        id: spec.id,
                        category: spec.category,
                        name: spec.name.to_string(),
                        description_no: spec.description_no.to_string(),
                        labels,
                        metrics,
                        duration_ms: elapsed_ms,
                        params_summary: out.params_summary,
                        note: out.note,
                        error: None,
                    }
                }
                Err(_) => {
                    let labels = vec![-1; data.n_samples];
                    let metrics = evaluate_clustering(data, &labels);
                    AlgorithmRunResult {
                        id: spec.id,
                        category: spec.category,
                        name: spec.name.to_string(),
                        description_no: spec.description_no.to_string(),
                        labels,
                        metrics,
                        duration_ms: elapsed_ms,
                        params_summary: String::new(),
                        note: "Intern feil fanget trygt".to_string(),
                        error: Some("Algoritmen støtte på en numerisk feil.".to_string()),
                    }
                }
            };

            if let Some(p) = progress {
                p.completed.fetch_add(1, Ordering::SeqCst);
            }
            Some(res)
        })
        .collect();

    results.sort_by_key(|r| r.id);
    results
}
