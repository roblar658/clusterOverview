use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

use crate::math::{build_knn_graph, normalize_vec, pairwise_dist_matrix, pca_2d};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColumnKind {
    Numeric,
    Categorical,
    Binary,
}

impl ColumnKind {
    pub fn label_no(&self) -> &'static str {
        match self {
            ColumnKind::Numeric => "Numerisk",
            ColumnKind::Categorical => "Kategorisk",
            ColumnKind::Binary => "Binær",
        }
    }
}

#[derive(Debug, Clone)]
pub struct ColumnProfile {
    pub index: usize,
    pub name: String,
    pub kind: ColumnKind,
    pub missing_count: usize,
    pub unique_count: usize,
    pub min: Option<f64>,
    pub max: Option<f64>,
    pub mean: Option<f64>,
    pub std: Option<f64>,
    pub sample_values: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct RawCsvDataset {
    pub path: Option<PathBuf>,
    pub display_name: String,
    pub encoding: String,
    pub delimiter: u8,
    pub headers: Vec<String>,
    pub rows: Vec<Vec<String>>,
    pub columns: Vec<ColumnProfile>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScalingMethod {
    ZScore,
    MinMax,
    None,
}

impl ScalingMethod {
    pub fn label_no(&self) -> &'static str {
        match self {
            ScalingMethod::ZScore => "Standardisering (Z-score)",
            ScalingMethod::MinMax => "Min-Max [0, 1]",
            ScalingMethod::None => "Ingen skalering",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MissingStrategy {
    MeanMode,
    DropRow,
}

impl MissingStrategy {
    pub fn label_no(&self) -> &'static str {
        match self {
            MissingStrategy::MeanMode => "Imputér (gjennomsnitt / modus)",
            MissingStrategy::DropRow => "Fjern rader med manglende verdier",
        }
    }
}

#[derive(Debug, Clone)]
pub struct PreprocessConfig {
    pub selected_columns: Vec<bool>,
    pub label_column: Option<usize>,
    pub time_column: Option<usize>,
    pub scaling: ScalingMethod,
    pub missing_strategy: MissingStrategy,
    pub max_samples: usize,
    pub seed: u64,
}

impl PreprocessConfig {
    pub fn default_for(raw: &RawCsvDataset) -> Self {
        let mut selected = vec![false; raw.columns.len()];
        let mut label_col = None;
        let mut time_col = None;

        for col in &raw.columns {
            let lower = col.name.to_lowercase();
            let looks_like_id = lower == "id"
                || lower.ends_with("_id")
                || lower == "index"
                || lower.contains("org.nr")
                || lower.contains("orgnr")
                || lower.ends_with(".nr.")
                || lower.ends_with("nr")
                || (col.kind == ColumnKind::Numeric
                    && raw.rows.len() > 20
                    && col.unique_count as f64 / raw.rows.len() as f64 > 0.95
                    && col.min.unwrap_or(0.0) >= 100_000_000.0);
            let is_constant = col.unique_count <= 1;
            let looks_like_label = lower == "label"
                || lower == "class"
                || lower == "species"
                || lower == "segment"
                || lower == "target"
                || lower == "cluster"
                || lower == "category";
            let looks_like_time =
                lower.contains("time") || lower.contains("date") || lower == "timestamp" || lower == "t";

            if looks_like_label && label_col.is_none() && col.unique_count >= 2 && col.unique_count <= 20 {
                label_col = Some(col.index);
            } else if is_constant {
                continue;
            } else if looks_like_time && time_col.is_none() && col.kind == ColumnKind::Numeric {
                time_col = Some(col.index);
                selected[col.index] = true;
            } else if !looks_like_id {
                if col.kind == ColumnKind::Numeric || col.kind == ColumnKind::Binary {
                    selected[col.index] = true;
                } else if col.unique_count <= 15 {
                    selected[col.index] = true;
                }
            }
        }

        // Ensure at least one column is selected
        if !selected.iter().any(|&b| b) {
            for col in &raw.columns {
                if Some(col.index) != label_col {
                    selected[col.index] = true;
                }
            }
        }

        Self {
            selected_columns: selected,
            label_column: label_col,
            time_column: time_col,
            scaling: ScalingMethod::ZScore,
            missing_strategy: MissingStrategy::MeanMode,
            max_samples: 240,
            seed: 42,
        }
    }
}

#[derive(Debug, Clone)]
pub struct PreparedData {
    pub dataset_name: String,
    pub total_raw_rows: usize,
    pub n_samples: usize,
    pub n_features: usize,
    pub feature_names: Vec<String>,
    pub raw_numeric_cols: Vec<Vec<f64>>,
    pub raw_numeric_names: Vec<String>,
    pub x: Vec<Vec<f64>>,
    pub x_unit: Vec<Vec<f64>>,
    pub x_nonneg: Vec<Vec<f64>>,
    pub counts: Vec<Vec<usize>>,
    pub x_binary: Vec<Vec<u8>>,
    pub categorical: Vec<Vec<usize>>,
    pub has_real_categorical: bool,
    pub time_values: Vec<f64>,
    pub pca_2d: Vec<[f64; 2]>,
    pub pca_variance: [f64; 2],
    pub dist_matrix: Vec<Vec<f64>>,
    pub knn_graph: Vec<Vec<(usize, f64)>>,
    pub ground_truth: Option<Vec<i32>>,
    pub ground_truth_name: Option<String>,
    pub source_row_indices: Vec<usize>,
    pub suggested_k: usize,
    pub scaling_label: String,
}

fn detect_delimiter(content: &str) -> u8 {
    let first_lines: Vec<&str> = content.lines().take(5).collect();
    if first_lines.is_empty() {
        return b',';
    }
    let candidates = [b',', b';', b'\t', b'|'];
    let mut best_delim = b',';
    let mut best_score = 0usize;
    for &delim in &candidates {
        let counts: Vec<usize> = first_lines
            .iter()
            .filter(|l| !l.trim().is_empty())
            .map(|l| l.as_bytes().iter().filter(|&&b| b == delim).count())
            .collect();
        if counts.is_empty() || counts[0] == 0 {
            continue;
        }
        let consistent = counts.iter().all(|&c| c == counts[0]);
        let score = counts[0] * if consistent { 10 } else { 1 };
        if score > best_score {
            best_score = score;
            best_delim = delim;
        }
    }
    best_delim
}

/// Decode raw file bytes as UTF-8 (with or without BOM), UTF-16 LE/BE, or Windows-1252 (Norwegian Excel default).
fn decode_csv_bytes(bytes: &[u8]) -> (String, &'static str) {
    // UTF-8 BOM
    if bytes.starts_with(&[0xEF, 0xBB, 0xBF]) {
        if let Ok(s) = std::str::from_utf8(&bytes[3..]) {
            return (s.to_string(), "UTF-8 (BOM)");
        }
    }
    // UTF-16 LE BOM
    if bytes.starts_with(&[0xFF, 0xFE]) && bytes.len() >= 4 {
        let u16s: Vec<u16> = bytes[2..]
            .chunks_exact(2)
            .map(|c| u16::from_le_bytes([c[0], c[1]]))
            .collect();
        return (String::from_utf16_lossy(&u16s), "UTF-16 LE");
    }
    // UTF-16 BE BOM
    if bytes.starts_with(&[0xFE, 0xFF]) && bytes.len() >= 4 {
        let u16s: Vec<u16> = bytes[2..]
            .chunks_exact(2)
            .map(|c| u16::from_be_bytes([c[0], c[1]]))
            .collect();
        return (String::from_utf16_lossy(&u16s), "UTF-16 BE");
    }
    // Standard valid UTF-8
    if let Ok(s) = std::str::from_utf8(bytes) {
        return (s.to_string(), "UTF-8");
    }
    // Fallback: Windows-1252 (CP1252 / ISO-8859-1) — lossless for Norwegian æøåÆØÅ and Excel exports
    const CP1252_80_9F: [char; 32] = [
        '\u{20AC}', '\u{0081}', '\u{201A}', '\u{0192}', '\u{201E}', '\u{2026}', '\u{2020}', '\u{2021}',
        '\u{02C6}', '\u{2030}', '\u{0160}', '\u{2039}', '\u{0152}', '\u{008D}', '\u{017D}', '\u{008F}',
        '\u{0090}', '\u{2018}', '\u{2019}', '\u{201C}', '\u{201D}', '\u{2022}', '\u{2013}', '\u{2014}',
        '\u{02DC}', '\u{2122}', '\u{0161}', '\u{203A}', '\u{0153}', '\u{009D}', '\u{017E}', '\u{0178}',
    ];
    let mut out = String::with_capacity(bytes.len());
    for &b in bytes {
        let ch = match b {
            0x80..=0x9F => CP1252_80_9F[(b - 0x80) as usize],
            _ => b as char,
        };
        out.push(ch);
    }
    (out, "Windows-1252")
}

fn parse_number(s: &str) -> Option<f64> {
    let mut t = s.trim();
    if t.is_empty() {
        return None;
    }
    // Strip common trailing/leading currency or percent symbols
    if let Some(stripped) = t.strip_suffix('%') {
        t = stripped.trim();
    } else if let Some(stripped) = t.strip_suffix("kr").or_else(|| t.strip_suffix("NOK")) {
        t = stripped.trim();
    }
    if t.is_empty() {
        return None;
    }

    // Fast path for clean ASCII float
    if let Ok(v) = t.parse::<f64>() {
        if v.is_finite() {
            return Some(v);
        }
    }

    // Normalize Norwegian/European numbers: remove thousand-separator spaces / NBSP / narrow NBSP
    let compact: String = t
        .chars()
        .filter(|&c| c != ' ' && c != '\u{00A0}' && c != '\u{202F}')
        .collect();
    if compact.is_empty() {
        return None;
    }
    // Reject strings containing letters (except scientific notation e/E)
    if compact.chars().any(|c| c.is_alphabetic() && c != 'e' && c != 'E') {
        return None;
    }

    if let Ok(v) = compact.parse::<f64>() {
        if v.is_finite() {
            return Some(v);
        }
    }
    // European decimal comma (e.g. "12345,67" or "12.345,67")
    if compact.contains(',') {
        let normalized = if compact.contains('.') {
            // Decide which is decimal separator by last occurrence
            let last_comma = compact.rfind(',').unwrap();
            let last_dot = compact.rfind('.').unwrap();
            if last_comma > last_dot {
                compact.replace('.', "").replace(',', ".")
            } else {
                compact.replace(',', "")
            }
        } else {
            compact.replace(',', ".")
        };
        if let Ok(v) = normalized.parse::<f64>() {
            if v.is_finite() {
                return Some(v);
            }
        }
    }
    None
}

impl RawCsvDataset {
    pub fn from_file(path: &Path) -> Result<Self, String> {
        let bytes = fs::read(path)
            .map_err(|e| format!("Kunne ikke lese CSV-fil '{}': {}", path.display(), e))?;
        let (content, enc) = decode_csv_bytes(&bytes);
        let name = path
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "data.csv".to_string());
        Self::from_str_content_with_enc(&content, name, Some(path.to_path_buf()), enc)
    }

    pub fn from_str_content(
        content: &str,
        display_name: String,
        path: Option<PathBuf>,
    ) -> Result<Self, String> {
        Self::from_str_content_with_enc(content, display_name, path, "UTF-8")
    }

    fn from_str_content_with_enc(
        content: &str,
        display_name: String,
        path: Option<PathBuf>,
        encoding: &str,
    ) -> Result<Self, String> {
        let cleaned = content.trim_start_matches('\u{feff}');
        if cleaned.trim().is_empty() {
            return Err("CSV-filen er tom.".to_string());
        }
        let delimiter = detect_delimiter(cleaned);
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(delimiter)
            .flexible(true)
            .has_headers(true)
            .from_reader(cleaned.as_bytes());

        let headers: Vec<String> = rdr
            .headers()
            .map_err(|e| format!("Feil ved lesing av CSV-overskrift: {}", e))?
            .iter()
            .enumerate()
            .map(|(i, h)| {
                let t = h.trim();
                if t.is_empty() {
                    format!("Kolonne_{}", i + 1)
                } else {
                    t.to_string()
                }
            })
            .collect();

        if headers.is_empty() {
            return Err("Ingen kolonner funnet i CSV-filen.".to_string());
        }

        let mut rows = Vec::new();
        for rec in rdr.records() {
            if let Ok(record) = rec {
                let mut row = Vec::with_capacity(headers.len());
                for i in 0..headers.len() {
                    row.push(record.get(i).unwrap_or("").trim().to_string());
                }
                if row.iter().any(|c| !c.is_empty()) {
                    rows.push(row);
                }
            }
        }

        if rows.len() < 4 {
            return Err(format!(
                "For få rader i CSV-filen ({}). Minst 4 datapunkter kreves.",
                rows.len()
            ));
        }

        let columns = Self::profile_columns(&headers, &rows);
        Ok(Self {
            path,
            display_name,
            encoding: encoding.to_string(),
            delimiter,
            headers,
            rows,
            columns,
        })
    }

    fn profile_columns(headers: &[String], rows: &[Vec<String>]) -> Vec<ColumnProfile> {
        let mut profiles = Vec::with_capacity(headers.len());
        for (col_idx, name) in headers.iter().enumerate() {
            let mut missing_count = 0usize;
            let mut nums = Vec::new();
            let mut uniques = HashSet::new();
            let mut samples = Vec::new();
            let mut non_missing = 0usize;

            for row in rows {
                let cell = row.get(col_idx).map(|s| s.as_str()).unwrap_or("").trim();
                let is_na = cell.is_empty()
                    || cell.eq_ignore_ascii_case("na")
                    || cell.eq_ignore_ascii_case("nan")
                    || cell.eq_ignore_ascii_case("null")
                    || cell == "?";
                if is_na {
                    missing_count += 1;
                    continue;
                }
                non_missing += 1;
                if uniques.len() < 5000 {
                    uniques.insert(cell.to_string());
                }
                if samples.len() < 4 && !samples.iter().any(|s: &String| s == cell) {
                    samples.push(cell.to_string());
                }
                if let Some(v) = parse_number(cell) {
                    nums.push(v);
                }
            }

            let numeric_ratio = if non_missing > 0 {
                nums.len() as f64 / non_missing as f64
            } else {
                0.0
            };

            let kind = if numeric_ratio >= 0.85 && !nums.is_empty() {
                let all_01 = nums.iter().all(|&v| (v - 0.0).abs() < 1e-9 || (v - 1.0).abs() < 1e-9);
                if all_01 && uniques.len() <= 2 {
                    ColumnKind::Binary
                } else {
                    ColumnKind::Numeric
                }
            } else if uniques.len() == 2 {
                ColumnKind::Binary
            } else {
                ColumnKind::Categorical
            };

            let (min, max, mean, std) = if !nums.is_empty() && numeric_ratio >= 0.85 {
                let mn = nums.iter().copied().fold(f64::INFINITY, f64::min);
                let mx = nums.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let mu = nums.iter().sum::<f64>() / nums.len() as f64;
                let var = nums.iter().map(|&v| (v - mu).powi(2)).sum::<f64>()
                    / nums.len().max(1) as f64;
                (Some(mn), Some(mx), Some(mu), Some(var.sqrt()))
            } else {
                (None, None, None, None)
            };

            profiles.push(ColumnProfile {
                index: col_idx,
                name: name.clone(),
                kind,
                missing_count,
                unique_count: uniques.len(),
                min,
                max,
                mean,
                std,
                sample_values: samples,
            });
        }
        profiles
    }

    pub fn built_in_demo_iris() -> Self {
        let csv_str = include_str!("../sample_data/iris_demo.csv");
        Self::from_str_content(csv_str, "iris_demo.csv (innebygd)".to_string(), None)
            .expect("Built-in iris demo CSV must parse")
    }

    pub fn built_in_demo_customers() -> Self {
        let csv_str = include_str!("../sample_data/customers_mixed.csv");
        Self::from_str_content(csv_str, "customers_mixed.csv (innebygd)".to_string(), None)
            .expect("Built-in customer demo CSV must parse")
    }

    /// Generates a rich synthetic multi-structure CSV (blobs + crescent + density cluster) for instant benchmarking.
    pub fn built_in_demo_multishape() -> Self {
        let mut rng = ChaCha8Rng::seed_from_u64(2026);
        use rand::Rng;
        let mut s = String::from("x1,x2,x3,x4,timestamp,region,true_group\n");
        // Cluster A: dense Gaussian blob
        for i in 0..35 {
            let x1 = -3.0 + rng.gen_range(-0.6..0.6);
            let x2 = 2.5 + rng.gen_range(-0.6..0.6);
            let x3 = 1.0 + rng.gen_range(-0.4..0.4);
            let x4 = 8.0 + rng.gen_range(-0.8..0.8);
            s.push_str(&format!("{:.3},{:.3},{:.3},{:.3},{},Nord,DenseBlob\n", x1, x2, x3, x4, i));
        }
        // Cluster B: curved crescent arc
        for i in 0..40 {
            let theta = std::f64::consts::PI * (i as f64) / 39.0;
            let x1 = 2.0 + 2.2 * theta.cos() + rng.gen_range(-0.22..0.22);
            let x2 = -1.0 + 2.2 * theta.sin() + rng.gen_range(-0.22..0.22);
            let x3 = 4.5 + 0.5 * theta.cos() + rng.gen_range(-0.3..0.3);
            let x4 = 3.0 + rng.gen_range(-0.5..0.5);
            s.push_str(&format!("{:.3},{:.3},{:.3},{:.3},{},Sor,Crescent\n", x1, x2, x3, x4, 35 + i));
        }
        // Cluster C: elongated diagonal ellipse
        for i in 0..40 {
            let t = rng.gen_range(-1.8..1.8);
            let x1 = -1.0 + t + rng.gen_range(-0.25..0.25);
            let x2 = -3.0 + 0.8 * t + rng.gen_range(-0.25..0.25);
            let x3 = -2.5 + rng.gen_range(-0.4..0.4);
            let x4 = 5.5 + 0.6 * t + rng.gen_range(-0.4..0.4);
            s.push_str(&format!("{:.3},{:.3},{:.3},{:.3},{},Vest,Ellipse\n", x1, x2, x3, x4, 75 + i));
        }
        // 5 mild outliers
        for i in 0..5 {
            let x1 = rng.gen_range(-5.5..5.5);
            let x2 = rng.gen_range(-5.5..5.5);
            let x3 = rng.gen_range(-4.0..6.0);
            let x4 = rng.gen_range(0.0..10.0);
            s.push_str(&format!("{:.3},{:.3},{:.3},{:.3},{},Ost,Outlier\n", x1, x2, x3, x4, 115 + i));
        }
        Self::from_str_content(&s, "multishape_benchmark.csv (syntetisk)".to_string(), None)
            .expect("Built-in synthetic CSV must parse")
    }
}

pub fn prepare_dataset(raw: &RawCsvDataset, cfg: &PreprocessConfig) -> Result<PreparedData, String> {
    let active_cols: Vec<usize> = cfg
        .selected_columns
        .iter()
        .enumerate()
        .filter(|&(idx, &sel)| sel && Some(idx) != cfg.label_column)
        .map(|(idx, _)| idx)
        .collect();

    if active_cols.is_empty() {
        return Err("Velg minst én kolonne for klyngeanalyse.".to_string());
    }

    // Filter rows if DropRow is chosen
    let mut valid_row_indices: Vec<usize> = Vec::with_capacity(raw.rows.len());
    for (r_idx, row) in raw.rows.iter().enumerate() {
        if cfg.missing_strategy == MissingStrategy::DropRow {
            let has_missing = active_cols.iter().any(|&c| {
                let s = row.get(c).map(|v| v.trim()).unwrap_or("");
                s.is_empty()
                    || s.eq_ignore_ascii_case("na")
                    || s.eq_ignore_ascii_case("nan")
                    || s.eq_ignore_ascii_case("null")
            });
            if has_missing {
                continue;
            }
        }
        valid_row_indices.push(r_idx);
    }

    if valid_row_indices.len() < 4 {
        return Err("For få gyldige rader etter filtrering av manglende verdier (minst 4 kreves).".to_string());
    }

    // Deterministic subsampling if dataset exceeds max_samples
    let sample_cap = cfg.max_samples.max(20);
    if valid_row_indices.len() > sample_cap {
        let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed);
        valid_row_indices.shuffle(&mut rng);
        valid_row_indices.truncate(sample_cap);
        valid_row_indices.sort_unstable();
    }

    let n = valid_row_indices.len();
    let mut feature_names = Vec::new();
    let mut raw_numeric_cols: Vec<Vec<f64>> = Vec::new();
    let mut raw_numeric_names: Vec<String> = Vec::new();
    let mut numeric_matrix_cols: Vec<Vec<f64>> = Vec::new();
    let mut real_cat_cols: Vec<Vec<usize>> = Vec::new();

    for &c_idx in &active_cols {
        let profile = &raw.columns[c_idx];
        match profile.kind {
            ColumnKind::Numeric => {
                let fallback_mean = profile.mean.unwrap_or(0.0);
                let vals: Vec<f64> = valid_row_indices
                    .iter()
                    .map(|&r| {
                        raw.rows[r]
                            .get(c_idx)
                            .and_then(|s| parse_number(s))
                            .unwrap_or(fallback_mean)
                    })
                    .collect();
                raw_numeric_names.push(profile.name.clone());
                raw_numeric_cols.push(vals.clone());
                feature_names.push(profile.name.clone());
                numeric_matrix_cols.push(scale_column(&vals, cfg.scaling));
            }
            ColumnKind::Binary => {
                // Map unique non-empty values to 0 / 1
                let mut val_map: HashMap<String, usize> = HashMap::new();
                let mut cat_vals = Vec::with_capacity(n);
                let mut num_vals = Vec::with_capacity(n);
                for &r in &valid_row_indices {
                    let s = raw.rows[r].get(c_idx).map(|v| v.trim()).unwrap_or("");
                    if let Some(v) = parse_number(s) {
                        let b = if v > 0.5 { 1 } else { 0 };
                        cat_vals.push(b);
                        num_vals.push(b as f64);
                    } else {
                        let next_id = val_map.len();
                        let id = *val_map.entry(s.to_string()).or_insert(next_id);
                        let b = if id > 0 { 1 } else { 0 };
                        cat_vals.push(b);
                        num_vals.push(b as f64);
                    }
                }
                real_cat_cols.push(cat_vals);
                raw_numeric_names.push(profile.name.clone());
                raw_numeric_cols.push(num_vals.clone());
                feature_names.push(profile.name.clone());
                numeric_matrix_cols.push(scale_column(&num_vals, cfg.scaling));
            }
            ColumnKind::Categorical => {
                // Compute mode for missing imputation
                let mut freq: HashMap<String, usize> = HashMap::new();
                for &r in &valid_row_indices {
                    let s = raw.rows[r].get(c_idx).map(|v| v.trim()).unwrap_or("");
                    if !s.is_empty() {
                        *freq.entry(s.to_string()).or_insert(0) += 1;
                    }
                }
                let mode = freq
                    .iter()
                    .max_by_key(|(_, &cnt)| cnt)
                    .map(|(k, _)| k.clone())
                    .unwrap_or_else(|| "Ukjent".to_string());

                let mut categories: Vec<String> = freq.keys().cloned().collect();
                categories.sort();
                let cat_map: HashMap<String, usize> = categories
                    .iter()
                    .enumerate()
                    .map(|(i, s)| (s.clone(), i))
                    .collect();

                let cat_vals: Vec<usize> = valid_row_indices
                    .iter()
                    .map(|&r| {
                        let s = raw.rows[r].get(c_idx).map(|v| v.trim()).unwrap_or("");
                        let key = if s.is_empty() { &mode } else { s };
                        *cat_map.get(key).unwrap_or(&0)
                    })
                    .collect();
                real_cat_cols.push(cat_vals.clone());

                // One-hot encode up to 8 categories, otherwise frequency/ordinal encode
                if categories.len() <= 8 && categories.len() >= 2 {
                    for (ci, cname) in categories.iter().enumerate() {
                        let col: Vec<f64> = cat_vals
                            .iter()
                            .map(|&v| if v == ci { 1.0 } else { 0.0 })
                            .collect();
                        feature_names.push(format!("{}={}", profile.name, cname));
                        numeric_matrix_cols.push(scale_column(&col, cfg.scaling));
                    }
                } else {
                    let col: Vec<f64> = cat_vals.iter().map(|&v| v as f64).collect();
                    feature_names.push(profile.name.clone());
                    numeric_matrix_cols.push(scale_column(&col, cfg.scaling));
                }
            }
        }
    }

    if numeric_matrix_cols.is_empty() {
        return Err("Ingen brukbare egenskaper kunne trekkes ut fra valgte kolonner.".to_string());
    }

    let d = numeric_matrix_cols.len();
    let mut x = vec![vec![0.0; d]; n];
    for j in 0..d {
        for i in 0..n {
            x[i][j] = numeric_matrix_cols[j][i];
        }
    }

    let x_unit: Vec<Vec<f64>> = x.iter().map(|r| normalize_vec(r)).collect();

    // Non-negative shifted matrix & integer counts (for NMF, LDA, pLSA, Poisson, DMM, BTM, etc.)
    let mut col_mins = vec![f64::INFINITY; d];
    let mut col_maxs = vec![f64::NEG_INFINITY; d];
    for i in 0..n {
        for j in 0..d {
            col_mins[j] = col_mins[j].min(x[i][j]);
            col_maxs[j] = col_maxs[j].max(x[i][j]);
        }
    }
    let mut x_nonneg = vec![vec![0.0; d]; n];
    let mut counts = vec![vec![0usize; d]; n];
    for i in 0..n {
        for j in 0..d {
            let span = (col_maxs[j] - col_mins[j]).max(1e-9);
            let norm01 = (x[i][j] - col_mins[j]) / span;
            x_nonneg[i][j] = norm01 + 0.05;
            counts[i][j] = ((norm01 * 12.0).round() as usize).max(0);
        }
    }

    // Binary representation (above column median -> 1, else 0)
    let mut medians = vec![0.0; d];
    for j in 0..d {
        let mut col_vals: Vec<f64> = (0..n).map(|i| x[i][j]).collect();
        col_vals.sort_by(|a, b| a.total_cmp(b));
        medians[j] = col_vals[n / 2];
    }
    let mut x_binary = vec![vec![0u8; d]; n];
    for i in 0..n {
        for j in 0..d {
            x_binary[i][j] = if x[i][j] >= medians[j] { 1 } else { 0 };
        }
    }

    // Categorical matrix: use real categorical columns if >= 2, else combine with 4-bin quantile discretized numeric cols
    let has_real_categorical = !real_cat_cols.is_empty();
    let mut cat_cols = real_cat_cols;
    if cat_cols.len() < 2 {
        for j in 0..d.min(8) {
            let mut col_vals: Vec<f64> = (0..n).map(|i| x[i][j]).collect();
            col_vals.sort_by(|a, b| a.total_cmp(b));
            let q1 = col_vals[n / 4];
            let q2 = col_vals[n / 2];
            let q3 = col_vals[(3 * n) / 4];
            let binned: Vec<usize> = (0..n)
                .map(|i| {
                    let v = x[i][j];
                    if v <= q1 {
                        0
                    } else if v <= q2 {
                        1
                    } else if v <= q3 {
                        2
                    } else {
                        3
                    }
                })
                .collect();
            cat_cols.push(binned);
        }
    }
    let c_dim = cat_cols.len();
    let mut categorical = vec![vec![0usize; c_dim]; n];
    for j in 0..c_dim {
        for i in 0..n {
            categorical[i][j] = cat_cols[j][i];
        }
    }

    // Time vector (from time column if selected, else normalized index)
    let time_values: Vec<f64> = if let Some(tc) = cfg.time_column {
        let raw_t: Vec<f64> = valid_row_indices
            .iter()
            .enumerate()
            .map(|(idx, &r)| {
                raw.rows[r]
                    .get(tc)
                    .and_then(|s| parse_number(s))
                    .unwrap_or(idx as f64)
            })
            .collect();
        scale_column(&raw_t, ScalingMethod::MinMax)
    } else {
        (0..n)
            .map(|i| i as f64 / (n.saturating_sub(1)).max(1) as f64)
            .collect()
    };

    // Ground truth labels (optional)
    let (ground_truth, ground_truth_name) = if let Some(lc) = cfg.label_column {
        let mut label_map: HashMap<String, i32> = HashMap::new();
        let mut gt = Vec::with_capacity(n);
        for &r in &valid_row_indices {
            let s = raw.rows[r].get(lc).map(|v| v.trim()).unwrap_or("").to_string();
            let next_id = label_map.len() as i32;
            let id = *label_map.entry(s).or_insert(next_id);
            gt.push(id);
        }
        (Some(gt), Some(raw.columns[lc].name.clone()))
    } else {
        (None, None)
    };

    let (pca_coords, pca_variance) = pca_2d(&x);
    let dist_matrix = pairwise_dist_matrix(&x);
    let knn_k = ((n as f64).sqrt().round() as usize).clamp(5, 14).min(n.saturating_sub(1).max(1));
    let knn_graph = build_knn_graph(&dist_matrix, knn_k);
    let suggested_k = estimate_suggested_k(&x, &dist_matrix, cfg.seed);

    Ok(PreparedData {
        dataset_name: raw.display_name.clone(),
        total_raw_rows: raw.rows.len(),
        n_samples: n,
        n_features: d,
        feature_names,
        raw_numeric_cols,
        raw_numeric_names,
        x,
        x_unit,
        x_nonneg,
        counts,
        x_binary,
        categorical,
        has_real_categorical,
        time_values,
        pca_2d: pca_coords,
        pca_variance,
        dist_matrix,
        knn_graph,
        ground_truth,
        ground_truth_name,
        source_row_indices: valid_row_indices,
        suggested_k,
        scaling_label: cfg.scaling.label_no().to_string(),
    })
}

fn scale_column(vals: &[f64], method: ScalingMethod) -> Vec<f64> {
    if vals.is_empty() {
        return vec![];
    }
    match method {
        ScalingMethod::None => vals.to_vec(),
        ScalingMethod::ZScore => {
            let mean = vals.iter().sum::<f64>() / vals.len() as f64;
            let var =
                vals.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / vals.len().max(1) as f64;
            let std = var.sqrt();
            if std > 1e-12 {
                vals.iter().map(|&v| (v - mean) / std).collect()
            } else {
                vec![0.0; vals.len()]
            }
        }
        ScalingMethod::MinMax => {
            let mn = vals.iter().copied().fold(f64::INFINITY, f64::min);
            let mx = vals.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let span = mx - mn;
            if span > 1e-12 {
                vals.iter().map(|&v| (v - mn) / span).collect()
            } else {
                vec![0.5; vals.len()]
            }
        }
    }
}

fn estimate_suggested_k(x: &[Vec<f64>], dist: &[Vec<f64>], seed: u64) -> usize {
    let n = x.len();
    if n < 6 {
        return 2;
    }
    let max_k = 7.min(n / 3).max(2);
    let mut best_k = 3.min(max_k);
    let mut best_sil = f64::NEG_INFINITY;
    for k in 2..=max_k {
        let labels = quick_kmeans_pp(x, k, seed);
        if let Some(sil) = crate::metrics::silhouette_score(dist, &labels) {
            if sil > best_sil {
                best_sil = sil;
                best_k = k;
            }
        }
    }
    best_k
}

fn quick_kmeans_pp(x: &[Vec<f64>], k: usize, seed: u64) -> Vec<i32> {
    use rand::Rng;
    let n = x.len();
    let d = x[0].len();
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut centers = Vec::with_capacity(k);
    centers.push(x[rng.gen_range(0..n)].clone());
    for _ in 1..k {
        let dists: Vec<f64> = x
            .iter()
            .map(|p| {
                centers
                    .iter()
                    .map(|c| crate::math::euclidean_sq(p, c))
                    .fold(f64::INFINITY, f64::min)
            })
            .collect();
        let total: f64 = dists.iter().sum();
        if total <= 1e-12 {
            centers.push(x[rng.gen_range(0..n)].clone());
        } else {
            let mut r = rng.gen::<f64>() * total;
            let mut chosen = n - 1;
            for (i, &w) in dists.iter().enumerate() {
                r -= w;
                if r <= 0.0 {
                    chosen = i;
                    break;
                }
            }
            centers.push(x[chosen].clone());
        }
    }

    let mut labels = vec![0i32; n];
    for _ in 0..20 {
        let mut changed = false;
        for i in 0..n {
            let mut best_c = 0;
            let mut best_d = f64::INFINITY;
            for (c, center) in centers.iter().enumerate() {
                let d2 = crate::math::euclidean_sq(&x[i], center);
                if d2 < best_d {
                    best_d = d2;
                    best_c = c;
                }
            }
            if labels[i] != best_c as i32 {
                labels[i] = best_c as i32;
                changed = true;
            }
        }
        if !changed {
            break;
        }
        let mut sums = vec![vec![0.0; d]; k];
        let mut cnts = vec![0usize; k];
        for i in 0..n {
            let c = labels[i] as usize;
            cnts[c] += 1;
            for j in 0..d {
                sums[c][j] += x[i][j];
            }
        }
        for c in 0..k {
            if cnts[c] > 0 {
                for j in 0..d {
                    centers[c][j] = sums[c][j] / cnts[c] as f64;
                }
            }
        }
    }
    labels
}
