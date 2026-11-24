use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::Path;

use chrono::Local;
use flate2::write::ZlibEncoder;
use flate2::Compression;
use pdf_writer::{Content, Filter, Finish, Name, Pdf, Rect, Ref, Str, TextStr};

use crate::algorithms::{AlgorithmRunResult, Category};
use crate::data::{PreparedData, RawCsvDataset};

#[derive(Debug, Clone)]
pub struct PdfReportOptions {
    pub title: String,
    pub subtitle: String,
    pub author: String,
    pub notes: String,
    pub include_scatter_cards: bool,
    pub include_appendix_table: bool,
}

impl Default for PdfReportOptions {
    fn default() -> Self {
        Self {
            title: "ClusterOverview — Klyngeanalyserapport".to_string(),
            subtitle: "Komparativ evaluering av klyngealgoritmer".to_string(),
            author: "ClusterOverview (Rust)".to_string(),
            notes: String::new(),
            include_scatter_cards: true,
            include_appendix_table: true,
        }
    }
}

// A4 dimensions in PDF points (1 pt = 1/72 inch)
const PAGE_W: f32 = 595.28;
const PAGE_H: f32 = 841.89;
const MARGIN_X: f32 = 34.0;
const CONTENT_W: f32 = PAGE_W - 2.0 * MARGIN_X;

// Distinct, high-contrast colorblind-friendly cluster palette (RGB 0..1)
pub const CLUSTER_PALETTE_U8: [[u8; 3]; 12] = [
    [37, 99, 235],   // Blue
    [220, 38, 38],   // Red
    [5, 150, 105],   // Emerald
    [217, 119, 6],   // Amber
    [124, 58, 237],  // Violet
    [8, 145, 178],   // Cyan
    [219, 39, 119],  // Pink
    [101, 163, 13],  // Lime
    [234, 88, 12],   // Orange
    [79, 70, 229],   // Indigo
    [13, 148, 136],  // Teal
    [147, 51, 234],  // Purple
];

fn rgb(r: u8, g: u8, b: u8) -> (f32, f32, f32) {
    (r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0)
}

fn cluster_rgb(cid: i32) -> (f32, f32, f32) {
    if cid < 0 {
        return rgb(148, 163, 184); // Noise gray
    }
    let c = CLUSTER_PALETTE_U8[(cid as usize) % CLUSTER_PALETTE_U8.len()];
    rgb(c[0], c[1], c[2])
}

/// Convert UTF-8 string (including Norwegian æøåÆØÅ, en-dash, superscripts) to PDF WinAnsiEncoding bytes.
fn to_winansi(s: &str) -> Vec<u8> {
    let mut out = Vec::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '\u{0020}'..='\u{007E}' => out.push(ch as u8),
            // Latin-1 Supplement (includes æ 0xE6, ø 0xF8, å 0xE5, Æ 0xC6, Ø 0xD8, Å 0xC5, ² 0xB2, ± 0xB1, × 0xD7, µ 0xB5)
            '\u{00A0}'..='\u{00FF}' => out.push(ch as u8),
            '–' => out.push(0x96),
            '—' => out.push(0x97),
            '•' => out.push(0x95),
            '…' => out.push(0x85),
            '‘' => out.push(0x91),
            '’' => out.push(0x92),
            '“' => out.push(0x93),
            '”' => out.push(0x94),
            '€' => out.push(0x80),
            'ε' => out.extend_from_slice(b"eps"),
            'α' => out.extend_from_slice(b"alpha"),
            'β' => out.extend_from_slice(b"beta"),
            'γ' => out.extend_from_slice(b"gamma"),
            'λ' => out.extend_from_slice(b"lambda"),
            'ρ' => out.extend_from_slice(b"rho"),
            'δ' => out.extend_from_slice(b"delta"),
            'σ' => out.extend_from_slice(b"sigma"),
            '≈' => out.push(b'~'),
            '≤' => out.extend_from_slice(b"<="),
            '≥' => out.extend_from_slice(b">="),
            '→' => out.extend_from_slice(b"->"),
            '∩' => out.push(b'&'),
            '∪' => out.push(b'U'),
            _ => out.push(b'?'),
        }
    }
    out
}

/// Approximate width of Helvetica string in points for reliable text wrapping & right alignment.
fn approx_text_width(text: &str, font_size: f32, bold: bool) -> f32 {
    let mut units = 0.0_f32;
    for ch in text.chars() {
        let w = match ch {
            'i' | 'l' | '!' | '.' | ',' | ':' | ';' | '|' | '\'' => 0.28,
            'f' | 'j' | 'r' | 't' | '(' | ')' | '[' | ']' | '-' | '/' => 0.34,
            ' ' => 0.28,
            'm' | 'w' | 'M' | 'W' | 'Æ' | 'æ' => 0.80,
            'A'..='Z' | 'Ø' | 'Å' => 0.68,
            '0'..='9' => 0.56,
            _ => 0.53,
        };
        units += w;
    }
    units * font_size * if bold { 1.05 } else { 1.0 }
}

fn wrap_text(text: &str, max_width: f32, font_size: f32, bold: bool, max_lines: usize) -> Vec<String> {
    let mut lines = Vec::new();
    let mut current = String::new();
    for word in text.split_whitespace() {
        let candidate = if current.is_empty() {
            word.to_string()
        } else {
            format!("{} {}", current, word)
        };
        if approx_text_width(&candidate, font_size, bold) <= max_width || current.is_empty() {
            current = candidate;
        } else {
            lines.push(current);
            if lines.len() >= max_lines {
                // Add ellipsis to last line
                if let Some(last) = lines.last_mut() {
                    if last.len() > 3 {
                        let cut = last.char_indices().rev().nth(2).map(|(i, _)| i).unwrap_or(last.len());
                        last.truncate(cut);
                        last.push_str("...");
                    }
                }
                return lines;
            }
            current = word.to_string();
        }
    }
    if !current.is_empty() && lines.len() < max_lines {
        lines.push(current);
    }
    lines
}

fn truncate_to_width(text: &str, max_width: f32, font_size: f32, bold: bool) -> String {
    if approx_text_width(text, font_size, bold) <= max_width {
        return text.to_string();
    }
    let mut out = String::new();
    for ch in text.chars() {
        let mut cand = out.clone();
        cand.push(ch);
        cand.push_str("...");
        if approx_text_width(&cand, font_size, bold) > max_width {
            break;
        }
        out.push(ch);
    }
    out.push_str("...");
    out
}

#[derive(Clone, Copy)]
enum FontKind {
    Regular,
    Bold,
    Italic,
    Mono,
}

impl FontKind {
    fn pdf_name(&self) -> Name<'static> {
        match self {
            FontKind::Regular => Name(b"F1"),
            FontKind::Bold => Name(b"F2"),
            FontKind::Italic => Name(b"F3"),
            FontKind::Mono => Name(b"F4"),
        }
    }
    fn is_bold(&self) -> bool {
        matches!(self, FontKind::Bold)
    }
}

struct PageBuilder {
    content: Content,
}

impl PageBuilder {
    fn new() -> Self {
        Self { content: Content::new() }
    }

    fn fill_rect(&mut self, x: f32, y: f32, w: f32, h: f32, color: (f32, f32, f32)) {
        self.content.save_state();
        self.content.set_fill_rgb(color.0, color.1, color.2);
        self.content.rect(x, y, w, h);
        self.content.fill_nonzero();
        self.content.restore_state();
    }

    fn stroke_rect(&mut self, x: f32, y: f32, w: f32, h: f32, color: (f32, f32, f32), lw: f32) {
        self.content.save_state();
        self.content.set_stroke_rgb(color.0, color.1, color.2);
        self.content.set_line_width(lw);
        self.content.rect(x, y, w, h);
        self.content.stroke();
        self.content.restore_state();
    }

    fn line(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, color: (f32, f32, f32), lw: f32) {
        self.content.save_state();
        self.content.set_stroke_rgb(color.0, color.1, color.2);
        self.content.set_line_width(lw);
        self.content.move_to(x1, y1);
        self.content.line_to(x2, y2);
        self.content.stroke();
        self.content.restore_state();
    }

    fn fill_circle(&mut self, cx: f32, cy: f32, r: f32, color: (f32, f32, f32)) {
        // 4-segment Bezier circle approximation (k = 0.55228475)
        let k = r * 0.55228475;
        self.content.save_state();
        self.content.set_fill_rgb(color.0, color.1, color.2);
        self.content.move_to(cx + r, cy);
        self.content.cubic_to(cx + r, cy + k, cx + k, cy + r, cx, cy + r);
        self.content.cubic_to(cx - k, cy + r, cx - r, cy + k, cx - r, cy);
        self.content.cubic_to(cx - r, cy - k, cx - k, cy - r, cx, cy - r);
        self.content.cubic_to(cx + k, cy - r, cx + r, cy - k, cx + r, cy);
        self.content.close_path();
        self.content.fill_nonzero();
        self.content.restore_state();
    }

    fn draw_cross(&mut self, cx: f32, cy: f32, r: f32, color: (f32, f32, f32)) {
        self.content.save_state();
        self.content.set_stroke_rgb(color.0, color.1, color.2);
        self.content.set_line_width(0.65);
        self.content.move_to(cx - r, cy - r);
        self.content.line_to(cx + r, cy + r);
        self.content.move_to(cx - r, cy + r);
        self.content.line_to(cx + r, cy - r);
        self.content.stroke();
        self.content.restore_state();
    }

    fn text(&mut self, x: f32, y: f32, size: f32, font: FontKind, color: (f32, f32, f32), s: &str) {
        let bytes = to_winansi(s);
        self.content.save_state();
        self.content.set_fill_rgb(color.0, color.1, color.2);
        self.content.begin_text();
        self.content.set_font(font.pdf_name(), size);
        self.content.next_line(x, y);
        self.content.show(Str(&bytes));
        self.content.end_text();
        self.content.restore_state();
    }

    fn text_right(&mut self, x_right: f32, y: f32, size: f32, font: FontKind, color: (f32, f32, f32), s: &str) {
        let w = approx_text_width(s, size, font.is_bold());
        self.text(x_right - w, y, size, font, color, s);
    }

    fn draw_scatter_plot(
        &mut self,
        x: f32,
        y: f32,
        w: f32,
        h: f32,
        pca_2d: &[[f64; 2]],
        bounds: (f64, f64, f64, f64),
        labels: &[i32],
        point_r: f32,
    ) {
        // Plot frame
        self.fill_rect(x, y, w, h, rgb(252, 253, 255));
        self.stroke_rect(x, y, w, h, rgb(203, 213, 225), 0.6);

        // Subtle gridlines (25%, 50%, 75%)
        for frac in [0.25_f32, 0.5, 0.75] {
            let gx = x + w * frac;
            let gy = y + h * frac;
            self.line(gx, y, gx, y + h, rgb(236, 241, 247), 0.4);
            self.line(x, gy, x + w, gy, rgb(236, 241, 247), 0.4);
        }

        let pad = 6.0_f32;
        let iw = (w - 2.0 * pad).max(4.0);
        let ih = (h - 2.0 * pad).max(4.0);
        let (min_x, max_x, min_y, max_y) = bounds;
        let sx = (max_x - min_x).max(1e-9);
        let sy = (max_y - min_y).max(1e-9);

        // Draw noise first so real clusters render crisply on top
        for (i, pt) in pca_2d.iter().enumerate() {
            let lbl = labels.get(i).copied().unwrap_or(-1);
            if lbl < 0 {
                let px = x + pad + ((pt[0] - min_x) / sx) as f32 * iw;
                let py = y + pad + ((pt[1] - min_y) / sy) as f32 * ih;
                self.draw_cross(px, py, point_r * 0.95, rgb(148, 163, 184));
            }
        }
        for (i, pt) in pca_2d.iter().enumerate() {
            let lbl = labels.get(i).copied().unwrap_or(-1);
            if lbl >= 0 {
                let px = x + pad + ((pt[0] - min_x) / sx) as f32 * iw;
                let py = y + pad + ((pt[1] - min_y) / sy) as f32 * ih;
                self.fill_circle(px, py, point_r, cluster_rgb(lbl));
            }
        }
    }

    fn finish(self) -> Vec<u8> {
        self.content.finish()
    }
}

fn compute_pca_bounds(pca: &[[f64; 2]]) -> (f64, f64, f64, f64) {
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for p in pca {
        min_x = min_x.min(p[0]);
        max_x = max_x.max(p[0]);
        min_y = min_y.min(p[1]);
        max_y = max_y.max(p[1]);
    }
    if !min_x.is_finite() {
        return (-1.0, 1.0, -1.0, 1.0);
    }
    (min_x, max_x, min_y, max_y)
}

fn fmt_opt(v: Option<f64>, decimals: usize) -> String {
    match v {
        Some(x) if x.is_finite() => format!("{:.*}", decimals, x),
        _ => "—".to_string(),
    }
}

fn draw_chrome(pb: &mut PageBuilder, header_right: &str, page_num: usize, total_pages: usize) {
    // Top subtle bar
    pb.line(MARGIN_X, PAGE_H - 34.0, PAGE_W - MARGIN_X, PAGE_H - 34.0, rgb(226, 232, 240), 0.7);
    pb.text(
        MARGIN_X,
        PAGE_H - 27.0,
        7.5,
        FontKind::Bold,
        rgb(51, 65, 85),
        "CLUSTEROVERVIEW  |  KLYNGEANALYSERAPPORT",
    );
    pb.text_right(
        PAGE_W - MARGIN_X,
        PAGE_H - 27.0,
        7.5,
        FontKind::Regular,
        rgb(100, 116, 139),
        header_right,
    );

    // Bottom footer
    pb.line(MARGIN_X, 32.0, PAGE_W - MARGIN_X, 32.0, rgb(226, 232, 240), 0.7);
    pb.text(
        MARGIN_X,
        20.0,
        7.5,
        FontKind::Regular,
        rgb(100, 116, 139),
        "Generert med ClusterOverview (Native Rust) — 200 klyngealgoritmer",
    );
    pb.text_right(
        PAGE_W - MARGIN_X,
        20.0,
        7.5,
        FontKind::Bold,
        rgb(51, 65, 85),
        &format!("Side {} av {}", page_num, total_pages),
    );
}

pub fn estimate_page_count(n_results: usize, opts: &PdfReportOptions) -> usize {
    let base = 3; // Cover + Dataset Profile + Quality Overview
    let cards = if opts.include_scatter_cards {
        // Rough estimate: 6 cards per page across up to 8 categories
        ((n_results as f64) / 5.2).ceil() as usize
    } else {
        0
    };
    let appendix = if opts.include_appendix_table {
        ((n_results as f64) / 34.0).ceil() as usize
    } else {
        0
    };
    base + cards + appendix
}

pub fn generate_pdf_report(
    out_path: &Path,
    raw: &RawCsvDataset,
    data: &PreparedData,
    results: &[AlgorithmRunResult],
    opts: &PdfReportOptions,
) -> Result<usize, String> {
    if results.is_empty() {
        return Err("Ingen algoritmeresultater å eksportere til PDF.".to_string());
    }

    let bounds = compute_pca_bounds(&data.pca_2d);
    let now_str = Local::now().format("%d.%m.%Y kl. %H:%M").to_string();

    // Group results by category (preserving 1..8 order)
    let mut by_cat: HashMap<Category, Vec<&AlgorithmRunResult>> = HashMap::new();
    for r in results {
        by_cat.entry(r.category).or_default().push(r);
    }

    // Build all pages into memory buffers first so we know exact total_pages and bookmark page indices!
    // Each entry: (content_builder_fn_output, optional_header_right, skip_chrome)
    struct BuiltPage {
        content: PageBuilder,
        header_right: String,
        is_cover: bool,
    }

    let mut pages: Vec<BuiltPage> = Vec::new();
    let mut bookmarks: Vec<(String, usize)> = Vec::new(); // (title, 0-based page_index)

    // Precompute aggregate KPIs
    let total_algos = results.len();
    let ok_algos = results.iter().filter(|r| r.error.is_none()).count();
    let total_ms: f64 = results.iter().map(|r| r.duration_ms).sum();
    let avg_ms = total_ms / total_algos.max(1) as f64;

    let best_sil = results
        .iter()
        .filter_map(|r| r.metrics.silhouette.map(|s| (r, s)))
        .max_by(|a, b| a.1.total_cmp(&b.1));
    let best_ch = results
        .iter()
        .filter_map(|r| r.metrics.calinski_harabasz.map(|v| (r, v)))
        .max_by(|a, b| a.1.total_cmp(&b.1));
    let fastest = results
        .iter()
        .min_by(|a, b| a.duration_ms.total_cmp(&b.duration_ms));

    // Reserve Page 0 for Cover (we'll populate TOC page numbers after laying out all pages!)
    bookmarks.push(("Forside & Sammendrag".to_string(), 0));
    pages.push(BuiltPage {
        content: PageBuilder::new(),
        header_right: String::new(),
        is_cover: true,
    });

    // =========================================================================
    // PAGE 2: Dataset Profile & 2D PCA Overview
    // =========================================================================
    bookmarks.push(("Datasettprofil & PCA-oversikt".to_string(), pages.len()));
    {
        let mut pb = PageBuilder::new();
        let mut y = PAGE_H - 64.0;
        pb.text(MARGIN_X, y, 16.0, FontKind::Bold, rgb(15, 23, 42), "Datasettprofil & 2D PCA-projeksjon");
        y -= 16.0;
        pb.text(
            MARGIN_X,
            y,
            9.0,
            FontKind::Regular,
            rgb(71, 85, 105),
            &format!(
                "Kilde: {}  •  {} rader analysert (av {} totalt)  •  {} egenskaper  •  Skalering: {}",
                data.dataset_name, data.n_samples, data.total_raw_rows, data.n_features, data.scaling_label
            ),
        );
        y -= 18.0;

        // Large 2D PCA Scatter Plot + Info Side Panel
        let plot_w = 320.0;
        let plot_h = 220.0;
        let plot_y = y - plot_h;
        let init_labels = data
            .ground_truth
            .clone()
            .unwrap_or_else(|| vec![0i32; data.n_samples]);
        pb.draw_scatter_plot(MARGIN_X, plot_y, plot_w, plot_h, &data.pca_2d, bounds, &init_labels, 2.7);
        pb.text(
            MARGIN_X + 6.0,
            plot_y - 11.0,
            7.5,
            FontKind::Italic,
            rgb(100, 116, 139),
            &format!(
                "PC1 ({:.1}% varians) vs. PC2 ({:.1}% varians) — {}",
                data.pca_variance[0] * 100.0,
                data.pca_variance[1] * 100.0,
                if let Some(gt) = &data.ground_truth_name {
                    format!("fargelagt etter fasitkolonne '{}'", gt)
                } else {
                    "rå datastruktur i redusert 2D-rom".to_string()
                }
            ),
        );

        // Right metadata card
        let rx = MARGIN_X + plot_w + 14.0;
        let rw = CONTENT_W - plot_w - 14.0;
        pb.fill_rect(rx, plot_y, rw, plot_h, rgb(248, 250, 252));
        pb.stroke_rect(rx, plot_y, rw, plot_h, rgb(226, 232, 240), 0.7);
        pb.fill_rect(rx, plot_y + plot_h - 4.0, rw, 4.0, rgb(37, 99, 235));

        let mut iy = plot_y + plot_h - 22.0;
        pb.text(rx + 10.0, iy, 10.0, FontKind::Bold, rgb(15, 23, 42), "Forbehandling & Oppsett");
        iy -= 18.0;
        let info_rows = [
            ("Rå rader i CSV:", format!("{}", data.total_raw_rows)),
            ("Analyserte punkter:", format!("{}", data.n_samples)),
            ("Aktive egenskaper:", format!("{}", data.n_features)),
            ("Skalering:", data.scaling_label.clone()),
            ("Anbefalt k (Silhouette):", format!("k = {}", data.suggested_k)),
            ("PC1 + PC2 forklart:", format!("{:.1}%", (data.pca_variance[0] + data.pca_variance[1]) * 100.0)),
            ("Fasitkolonne (ARI):", data.ground_truth_name.clone().unwrap_or_else(|| "Ingen".to_string())),
            ("Kategoriske kolonner:", if data.has_real_categorical { "Ja (ekte)".into() } else { "Kvantil-diskretisert".into() }),
        ];
        for (k_lbl, v_lbl) in info_rows {
            pb.text(rx + 10.0, iy, 8.0, FontKind::Regular, rgb(100, 116, 139), k_lbl);
            let val_short = truncate_to_width(&v_lbl, rw - 20.0, 8.2, true);
            pb.text(rx + 10.0, iy - 10.0, 8.2, FontKind::Bold, rgb(30, 41, 59), &val_short);
            iy -= 23.0;
        }

        // Column Profiling Table
        y = plot_y - 34.0;
        pb.text(MARGIN_X, y, 12.0, FontKind::Bold, rgb(15, 23, 42), "Kolonneprofil (CSV-struktur)");
        y -= 16.0;

        // Table Header
        pb.fill_rect(MARGIN_X, y - 4.0, CONTENT_W, 18.0, rgb(30, 41, 59));
        let cols_x = [MARGIN_X + 6.0, MARGIN_X + 135.0, MARGIN_X + 205.0, MARGIN_X + 260.0, MARGIN_X + 325.0, MARGIN_X + 390.0, MARGIN_X + 455.0];
        let hdrs = ["Kolonne", "Type", "Unike", "Mangl.", "Min", "Maks", "Snitt (Std)"];
        for (cx, h) in cols_x.iter().zip(hdrs.iter()) {
            pb.text(*cx, y + 1.0, 8.0, FontKind::Bold, rgb(255, 255, 255), h);
        }
        y -= 20.0;

        for (idx, col) in raw.columns.iter().take(18).enumerate() {
            if y < 60.0 {
                break;
            }
            if idx % 2 == 1 {
                pb.fill_rect(MARGIN_X, y - 4.0, CONTENT_W, 16.0, rgb(248, 250, 252));
            }
            let cname = truncate_to_width(&col.name, 122.0, 8.0, true);
            pb.text(cols_x[0], y, 8.0, FontKind::Bold, rgb(30, 41, 59), &cname);
            pb.text(cols_x[1], y, 8.0, FontKind::Regular, rgb(71, 85, 105), col.kind.label_no());
            pb.text(cols_x[2], y, 8.0, FontKind::Regular, rgb(51, 65, 85), &col.unique_count.to_string());
            pb.text(cols_x[3], y, 8.0, FontKind::Regular, rgb(51, 65, 85), &col.missing_count.to_string());
            pb.text(cols_x[4], y, 8.0, FontKind::Regular, rgb(51, 65, 85), &fmt_opt(col.min, 2));
            pb.text(cols_x[5], y, 8.0, FontKind::Regular, rgb(51, 65, 85), &fmt_opt(col.max, 2));
            let ms = match (col.mean, col.std) {
                (Some(m), Some(s)) => format!("{:.2} ({:.2})", m, s),
                _ => truncate_to_width(&col.sample_values.join(", "), 68.0, 7.5, false),
            };
            pb.text(cols_x[6], y, 7.8, FontKind::Regular, rgb(51, 65, 85), &ms);
            y -= 16.0;
        }
        if raw.columns.len() > 18 {
            pb.text(
                MARGIN_X,
                y - 4.0,
                7.8,
                FontKind::Italic,
                rgb(100, 116, 139),
                &format!("... og {} flere kolonner i CSV-filen.", raw.columns.len() - 18),
            );
        }

        pages.push(BuiltPage {
            content: pb,
            header_right: "Datasettprofil & PCA".to_string(),
            is_cover: false,
        });
    }

    // =========================================================================
    // PAGE 3: Category Comparison Bar Chart + Top 15 Leaderboard
    // =========================================================================
    bookmarks.push(("Kvalitetsoversikt & Toppliste".to_string(), pages.len()));
    {
        let mut pb = PageBuilder::new();
        let mut y = PAGE_H - 64.0;
        pb.text(MARGIN_X, y, 16.0, FontKind::Bold, rgb(15, 23, 42), "Kvalitetsoversikt & Topprangerte Algoritmer");
        y -= 16.0;
        pb.text(
            MARGIN_X,
            y,
            8.5,
            FontKind::Regular,
            rgb(71, 85, 105),
            "Sammenligning av Silhouette-kvalitet per algoritmekategori (snitt vs. beste i hver kategori) og topp 20 totalt.",
        );
        y -= 18.0;

        // Category Performance Horizontal Bar Chart
        let chart_h = 205.0;
        let chart_y = y - chart_h;
        pb.fill_rect(MARGIN_X, chart_y, CONTENT_W, chart_h, rgb(248, 250, 252));
        pb.stroke_rect(MARGIN_X, chart_y, CONTENT_W, chart_h, rgb(226, 232, 240), 0.7);

        pb.text(MARGIN_X + 10.0, y - 16.0, 9.5, FontKind::Bold, rgb(15, 23, 42), "Silhouette-score per kategori (Mørk = Snitt, Lys = Beste)");
        let label_w = 138.0;
        let bar_x0 = MARGIN_X + label_w;
        let bar_max_w = CONTENT_W - label_w - 86.0;
        let mut row_y = y - 38.0;

        for cat in Category::ALL {
            let list = by_cat.get(&cat).map(|v| v.as_slice()).unwrap_or(&[]);
            let sils: Vec<f64> = list.iter().filter_map(|r| r.metrics.silhouette).collect();
            let avg_s = if sils.is_empty() { 0.0 } else { sils.iter().sum::<f64>() / sils.len() as f64 };
            let max_s = sils.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let max_s_clean = if max_s.is_finite() { max_s } else { 0.0 };

            let acc = cat.accent_rgb();
            let acc_c = rgb(acc[0], acc[1], acc[2]);
            let acc_light = (
                (acc[0] as f32 * 0.45 + 255.0 * 0.55) / 255.0,
                (acc[1] as f32 * 0.45 + 255.0 * 0.55) / 255.0,
                (acc[2] as f32 * 0.45 + 255.0 * 0.55) / 255.0,
            );

            pb.text(
                MARGIN_X + 10.0,
                row_y + 2.0,
                8.0,
                FontKind::Bold,
                rgb(30, 41, 59),
                &format!("{}. {}", cat.number(), cat.short_name_no()),
            );

            // Track background
            pb.fill_rect(bar_x0, row_y, bar_max_w, 12.0, rgb(226, 232, 240));
            // Max bar
            let w_max = (max_s_clean.clamp(0.0, 1.0) as f32) * bar_max_w;
            if w_max > 0.5 {
                pb.fill_rect(bar_x0, row_y, w_max, 12.0, acc_light);
            }
            // Avg bar
            let w_avg = (avg_s.clamp(0.0, 1.0) as f32) * bar_max_w;
            if w_avg > 0.5 {
                pb.fill_rect(bar_x0, row_y + 2.5, w_avg, 7.0, acc_c);
            }

            pb.text(
                bar_x0 + bar_max_w + 8.0,
                row_y + 2.5,
                7.6,
                FontKind::Mono,
                rgb(51, 65, 85),
                &format!("{:.2} / {:.2}", avg_s, max_s_clean),
            );
            row_y -= 20.5;
        }

        // Top 20 Leaderboard Table
        y = chart_y - 24.0;
        pb.text(MARGIN_X, y, 12.0, FontKind::Bold, rgb(15, 23, 42), "Topp 20 Algoritmer (rangert etter Silhouette-score)");
        y -= 18.0;

        let mut ranked: Vec<&AlgorithmRunResult> = results.iter().collect();
        ranked.sort_by(|a, b| {
            let sa = a.metrics.silhouette.unwrap_or(-2.0);
            let sb = b.metrics.silhouette.unwrap_or(-2.0);
            sb.total_cmp(&sa)
        });

        pb.fill_rect(MARGIN_X, y - 4.0, CONTENT_W, 18.0, rgb(15, 23, 42));
        let tx = [
            MARGIN_X + 6.0,   // Rank
            MARGIN_X + 28.0,  // Name
            MARGIN_X + 215.0, // Category
            MARGIN_X + 305.0, // k
            MARGIN_X + 333.0, // Noise
            MARGIN_X + 372.0, // Sil
            MARGIN_X + 415.0, // DBI
            MARGIN_X + 455.0, // CHI
            MARGIN_X + 498.0, // Time
        ];
        let th = ["#", "Algoritme", "Kategori", "k", "Støy", "Sil.", "DBI", "CHI", "ms"];
        for (cx, h) in tx.iter().zip(th.iter()) {
            pb.text(*cx, y + 1.0, 7.8, FontKind::Bold, rgb(255, 255, 255), h);
        }
        y -= 19.0;

        for (rank, r) in ranked.iter().take(20).enumerate() {
            if y < 52.0 {
                break;
            }
            if rank % 2 == 1 {
                pb.fill_rect(MARGIN_X, y - 4.0, CONTENT_W, 17.5, rgb(248, 250, 252));
            }
            let acc = r.category.accent_rgb();
            pb.fill_rect(MARGIN_X, y - 4.0, 3.0, 17.5, rgb(acc[0], acc[1], acc[2]));

            pb.text(tx[0], y, 8.0, FontKind::Bold, rgb(37, 99, 235), &format!("{}", rank + 1));
            let nm = truncate_to_width(&format!("#{:03} {}", r.id, r.name), 182.0, 7.8, true);
            pb.text(tx[1], y, 7.8, FontKind::Bold, rgb(15, 23, 42), &nm);
            pb.text(tx[2], y, 7.4, FontKind::Regular, rgb(71, 85, 105), r.category.short_name_no());
            pb.text(tx[3], y, 7.8, FontKind::Regular, rgb(30, 41, 59), &r.metrics.n_clusters.to_string());
            pb.text(tx[4], y, 7.6, FontKind::Regular, rgb(100, 116, 139), &format!("{:.0}%", r.metrics.noise_ratio * 100.0));
            pb.text(tx[5], y, 7.8, FontKind::Bold, rgb(5, 150, 105), &fmt_opt(r.metrics.silhouette, 3));
            pb.text(tx[6], y, 7.6, FontKind::Regular, rgb(51, 65, 85), &fmt_opt(r.metrics.davies_bouldin, 2));
            pb.text(tx[7], y, 7.6, FontKind::Regular, rgb(51, 65, 85), &fmt_opt(r.metrics.calinski_harabasz, 1));
            pb.text(tx[8], y, 7.6, FontKind::Mono, rgb(100, 116, 139), &format!("{:.1}", r.duration_ms));
            y -= 17.5;
        }

        pages.push(BuiltPage {
            content: pb,
            header_right: "Kvalitetsoversikt & Toppliste".to_string(),
            is_cover: false,
        });
    }

    // =========================================================================
    // CATEGORY PAGES WITH 2D VECTOR SCATTER CARDS (6 cards / A4 page)
    // =========================================================================
    let mut category_start_pages: Vec<(Category, usize, usize)> = Vec::new();
    if opts.include_scatter_cards {
        for cat in Category::ALL {
            let Some(cat_list) = by_cat.get(&cat) else { continue };
            if cat_list.is_empty() {
                continue;
            }
            let start_pg = pages.len();
            bookmarks.push((cat.title_no().to_string(), start_pg));
            category_start_pages.push((cat, start_pg + 1, cat_list.len()));

            let acc = cat.accent_rgb();
            let acc_c = rgb(acc[0], acc[1], acc[2]);

            // Category summary stats
            let sils: Vec<f64> = cat_list.iter().filter_map(|r| r.metrics.silhouette).collect();
            let avg_sil = if sils.is_empty() { None } else { Some(sils.iter().sum::<f64>() / sils.len() as f64) };
            let best_in_cat = cat_list
                .iter()
                .filter_map(|r| r.metrics.silhouette.map(|s| (r, s)))
                .max_by(|a, b| a.1.total_cmp(&b.1));

            // Layout: 6 cards per page (2 cols x 3 rows). First page of category includes a compact 46pt category banner!
            let mut idx = 0usize;
            let mut first_in_cat = true;
            while idx < cat_list.len() {
                let mut pb = PageBuilder::new();
                let mut top_y = PAGE_H - 46.0;

                if first_in_cat {
                    // Category Hero Strip
                    let ban_h = 46.0;
                    let ban_y = top_y - ban_h;
                    pb.fill_rect(MARGIN_X, ban_y, CONTENT_W, ban_h, rgb(248, 250, 252));
                    pb.stroke_rect(MARGIN_X, ban_y, CONTENT_W, ban_h, rgb(203, 213, 225), 0.7);
                    pb.fill_rect(MARGIN_X, ban_y, 6.0, ban_h, acc_c);

                    pb.text(MARGIN_X + 14.0, ban_y + 28.0, 12.0, FontKind::Bold, rgb(15, 23, 42), cat.title_no());
                    let best_str = match best_in_cat {
                        Some((br, bs)) => format!("Beste: {} (Sil. {:.3})", br.name, bs),
                        None => "Ingen fler-klyngeløsning".to_string(),
                    };
                    let sub = format!(
                        "{} algoritmer  •  Snitt Silhouette: {}  •  {}",
                        cat_list.len(),
                        fmt_opt(avg_sil, 3),
                        best_str
                    );
                    pb.text(
                        MARGIN_X + 14.0,
                        ban_y + 11.0,
                        8.0,
                        FontKind::Regular,
                        rgb(71, 85, 105),
                        &truncate_to_width(&sub, CONTENT_W - 24.0, 8.0, false),
                    );
                    top_y = ban_y - 10.0;
                    first_in_cat = false;
                }

                // Compute 2x3 grid geometry fitting cleanly above footer (y = 42.0)
                let avail_h = top_y - 42.0;
                let gap_x = 11.0_f32;
                let gap_y = 9.0_f32;
                let card_w = (CONTENT_W - gap_x) * 0.5;
                let card_h = (avail_h - 2.0 * gap_y) / 3.0;

                for slot in 0..6 {
                    if idx >= cat_list.len() {
                        break;
                    }
                    let r = cat_list[idx];
                    idx += 1;
                    let col = slot % 2;
                    let row = slot / 2;
                    let cx = MARGIN_X + (col as f32) * (card_w + gap_x);
                    let cy = top_y - ((row + 1) as f32) * card_h - (row as f32) * gap_y;

                    draw_algorithm_card(&mut pb, cx, cy, card_w, card_h, r, data, bounds, acc_c);
                }

                pages.push(BuiltPage {
                    content: pb,
                    header_right: cat.title_no().to_string(),
                    is_cover: false,
                });
            }
        }
    }

    // =========================================================================
    // APPENDIX: COMPLETE COMPARISON TABLE (ALL EXECUTED ALGORITHMS)
    // =========================================================================
    let appendix_start_page = pages.len() + 1;
    if opts.include_appendix_table {
        bookmarks.push(("Komplett sammenligningstabell".to_string(), pages.len()));
        let rows_per_page = 36usize;
        for (chunk_idx, chunk) in results.chunks(rows_per_page).enumerate() {
            let mut pb = PageBuilder::new();
            let mut y = PAGE_H - 60.0;
            if chunk_idx == 0 {
                pb.text(
                    MARGIN_X,
                    y,
                    14.0,
                    FontKind::Bold,
                    rgb(15, 23, 42),
                    &format!("Vedlegg: Komplett sammenligningstabell ({} algoritmer)", results.len()),
                );
                y -= 18.0;
            }

            // Header
            pb.fill_rect(MARGIN_X, y - 4.0, CONTENT_W, 17.0, rgb(15, 23, 42));
            let tx = [
                MARGIN_X + 5.0,
                MARGIN_X + 29.0,
                MARGIN_X + 206.0,
                MARGIN_X + 292.0,
                MARGIN_X + 318.0,
                MARGIN_X + 356.0,
                MARGIN_X + 398.0,
                MARGIN_X + 438.0,
                MARGIN_X + 485.0,
            ];
            let hdrs = ["ID", "Algoritme", "Kategori", "k", "Støy", "Sil.", "DBI", "CHI", "Tid ms"];
            for (cx, h) in tx.iter().zip(hdrs.iter()) {
                pb.text(*cx, y + 1.0, 7.5, FontKind::Bold, rgb(255, 255, 255), h);
            }
            y -= 18.0;

            for (ri, r) in chunk.iter().enumerate() {
                if ri % 2 == 1 {
                    pb.fill_rect(MARGIN_X, y - 4.0, CONTENT_W, 17.5, rgb(248, 250, 252));
                }
                let acc = r.category.accent_rgb();
                pb.fill_rect(MARGIN_X, y - 4.0, 2.5, 17.5, rgb(acc[0], acc[1], acc[2]));

                pb.text(tx[0], y, 7.4, FontKind::Mono, rgb(71, 85, 105), &format!("{:03}", r.id));
                let nm = truncate_to_width(&r.name, 172.0, 7.4, true);
                pb.text(tx[1], y, 7.4, FontKind::Bold, rgb(15, 23, 42), &nm);
                pb.text(tx[2], y, 7.1, FontKind::Regular, rgb(100, 116, 139), r.category.short_name_no());
                pb.text(tx[3], y, 7.4, FontKind::Regular, rgb(30, 41, 59), &r.metrics.n_clusters.to_string());
                pb.text(tx[4], y, 7.2, FontKind::Regular, rgb(100, 116, 139), &format!("{:.0}%", r.metrics.noise_ratio * 100.0));

                let sil_col = match r.metrics.silhouette {
                    Some(s) if s >= 0.45 => rgb(5, 150, 105),
                    Some(s) if s >= 0.20 => rgb(217, 119, 6),
                    _ => rgb(100, 116, 139),
                };
                pb.text(tx[5], y, 7.4, FontKind::Bold, sil_col, &fmt_opt(r.metrics.silhouette, 3));
                pb.text(tx[6], y, 7.2, FontKind::Regular, rgb(51, 65, 85), &fmt_opt(r.metrics.davies_bouldin, 2));
                pb.text(tx[7], y, 7.2, FontKind::Regular, rgb(51, 65, 85), &fmt_opt(r.metrics.calinski_harabasz, 1));
                pb.text(tx[8], y, 7.2, FontKind::Mono, rgb(100, 116, 139), &format!("{:.1}", r.duration_ms));
                y -= 17.5;
            }

            pages.push(BuiltPage {
                content: pb,
                header_right: "Vedlegg — Komplett tabell".to_string(),
                is_cover: false,
            });
        }
    }

    // =========================================================================
    // NOW POPULATE PAGE 0 (COVER PAGE WITH EXACT TOC PAGE NUMBERS!)
    // =========================================================================
    let total_pages = pages.len();
    {
        let mut pb = PageBuilder::new();
        // Full-width dark executive hero banner
        let hero_h = 250.0;
        pb.fill_rect(0.0, PAGE_H - hero_h, PAGE_W, hero_h, rgb(15, 23, 42));
        // Accent diagonal/stripe bar
        pb.fill_rect(0.0, PAGE_H - hero_h, PAGE_W, 6.0, rgb(37, 99, 235));

        // Decorative mini cluster constellation on top right of cover banner
        let deco_pts = [
            (465.0, 770.0, 5.0, CLUSTER_PALETTE_U8[0]),
            (482.0, 782.0, 4.0, CLUSTER_PALETTE_U8[0]),
            (474.0, 756.0, 3.5, CLUSTER_PALETTE_U8[0]),
            (496.0, 768.0, 4.5, CLUSTER_PALETTE_U8[0]),
            (520.0, 710.0, 5.0, CLUSTER_PALETTE_U8[2]),
            (536.0, 722.0, 4.0, CLUSTER_PALETTE_U8[2]),
            (528.0, 696.0, 3.5, CLUSTER_PALETTE_U8[2]),
            (548.0, 706.0, 4.2, CLUSTER_PALETTE_U8[2]),
            (448.0, 684.0, 4.8, CLUSTER_PALETTE_U8[3]),
            (464.0, 672.0, 3.8, CLUSTER_PALETTE_U8[3]),
            (436.0, 668.0, 3.6, CLUSTER_PALETTE_U8[3]),
            (510.0, 790.0, 3.2, CLUSTER_PALETTE_U8[4]),
            (532.0, 778.0, 4.5, CLUSTER_PALETTE_U8[4]),
            (544.0, 792.0, 3.6, CLUSTER_PALETTE_U8[4]),
        ];
        for &(px, py, r, c) in &deco_pts {
            pb.fill_circle(px, py, r, rgb(c[0], c[1], c[2]));
        }

        pb.text(MARGIN_X, PAGE_H - 58.0, 9.5, FontKind::Bold, rgb(96, 165, 250), "CLUSTEROVERVIEW  •  ANALYSERAPPORT");
        let title_lines = wrap_text(&opts.title, 395.0, 24.0, true, 2);
        let mut ty = PAGE_H - 96.0;
        for line in &title_lines {
            pb.text(MARGIN_X, ty, 24.0, FontKind::Bold, rgb(255, 255, 255), line);
            ty -= 28.0;
        }
        pb.text(
            MARGIN_X,
            ty - 6.0,
            11.5,
            FontKind::Regular,
            rgb(203, 213, 225),
            &truncate_to_width(&opts.subtitle, 400.0, 11.5, false),
        );

        // Cover metadata chips inside dark hero
        let meta_y = PAGE_H - hero_h + 26.0;
        pb.text(
            MARGIN_X,
            meta_y + 14.0,
            8.5,
            FontKind::Regular,
            rgb(148, 163, 184),
            &format!("Datasett: {}   |   Forfatter: {}   |   Dato: {}", data.dataset_name, opts.author, now_str),
        );
        pb.text(
            MARGIN_X,
            meta_y,
            8.5,
            FontKind::Regular,
            rgb(148, 163, 184),
            &format!(
                "Omfang: {} algoritmer over 8 hovedkategorier   |   {} sider",
                total_algos, total_pages
            ),
        );

        // 4 Executive KPI Cards right below hero
        let kpi_top = PAGE_H - hero_h - 24.0;
        let kpi_h = 74.0;
        let kpi_gap = 10.0;
        let kpi_w = (CONTENT_W - 3.0 * kpi_gap) / 4.0;

        let best_sil_val = best_sil.map(|(_, s)| format!("{:.3}", s)).unwrap_or_else(|| "—".into());
        let best_sil_sub = best_sil.map(|(r, _)| truncate_to_width(&r.name, kpi_w - 16.0, 7.4, false)).unwrap_or_default();
        let best_ch_val = best_ch.map(|(_, v)| format!("{:.1}", v)).unwrap_or_else(|| "—".into());
        let best_ch_sub = best_ch.map(|(r, _)| truncate_to_width(&r.name, kpi_w - 16.0, 7.4, false)).unwrap_or_default();
        let fast_val = fastest.map(|r| format!("{:.2} ms", r.duration_ms)).unwrap_or_else(|| "—".into());
        let fast_sub = fastest.map(|r| truncate_to_width(&r.name, kpi_w - 16.0, 7.4, false)).unwrap_or_default();

        let kpis = [
            ("ALGORITMER KJØRT", format!("{} / {}", ok_algos, total_algos), format!("Snitt {:.1} ms / alg.", avg_ms), rgb(37, 99, 235)),
            ("BESTE SILHOUETTE", best_sil_val, best_sil_sub, rgb(5, 150, 105)),
            ("BESTE CALINSKI-H.", best_ch_val, best_ch_sub, rgb(124, 58, 237)),
            ("RASKESTE ALGORITME", fast_val, fast_sub, rgb(217, 119, 6)),
        ];

        for (i, (label, big, sub, accent)) in kpis.iter().enumerate() {
            let kx = MARGIN_X + i as f32 * (kpi_w + kpi_gap);
            let ky = kpi_top - kpi_h;
            pb.fill_rect(kx, ky, kpi_w, kpi_h, rgb(248, 250, 252));
            pb.stroke_rect(kx, ky, kpi_w, kpi_h, rgb(226, 232, 240), 0.7);
            pb.fill_rect(kx, ky + kpi_h - 3.5, kpi_w, 3.5, *accent);
            pb.text(kx + 8.0, ky + kpi_h - 17.0, 7.0, FontKind::Bold, rgb(100, 116, 139), label);
            pb.text(kx + 8.0, ky + kpi_h - 42.0, 16.5, FontKind::Bold, rgb(15, 23, 42), big);
            pb.text(kx + 8.0, ky + 10.0, 7.4, FontKind::Regular, rgb(71, 85, 105), sub);
        }

        // Executive Summary + Optional User Notes
        let mut y = kpi_top - kpi_h - 28.0;
        pb.text(MARGIN_X, y, 12.5, FontKind::Bold, rgb(15, 23, 42), "Innholdsfortegnelse & Rapportstruktur");
        y -= 16.0;

        pb.fill_rect(MARGIN_X, y - 270.0, CONTENT_W, 276.0, rgb(252, 253, 255));
        pb.stroke_rect(MARGIN_X, y - 270.0, CONTENT_W, 276.0, rgb(226, 232, 240), 0.7);

        let mut ty = y - 16.0;
        let mut toc_entries: Vec<(String, String, (f32, f32, f32))> = vec![
            ("Datasettprofil & 2D PCA-projeksjon".to_string(), "Side 2".to_string(), rgb(37, 99, 235)),
            ("Kvalitetsoversikt, Kategoridiagram & Topp-20".to_string(), "Side 3".to_string(), rgb(37, 99, 235)),
        ];
        for (cat, pg, cnt) in &category_start_pages {
            let acc = cat.accent_rgb();
            toc_entries.push((
                format!("{}  [{} alg.]", cat.title_no(), cnt),
                format!("Side {}", pg),
                rgb(acc[0], acc[1], acc[2]),
            ));
        }
        if opts.include_appendix_table {
            toc_entries.push((
                format!("Vedlegg: Komplett sammenligningstabell ({} algoritmer)", total_algos),
                format!("Side {}", appendix_start_page),
                rgb(51, 65, 85),
            ));
        }

        for (title, pg_str, dot_c) in &toc_entries {
            pb.fill_circle(MARGIN_X + 14.0, ty + 3.0, 3.0, *dot_c);
            pb.text(MARGIN_X + 24.0, ty, 9.0, FontKind::Bold, rgb(30, 41, 59), title);
            pb.text_right(MARGIN_X + CONTENT_W - 14.0, ty, 9.0, FontKind::Mono, rgb(71, 85, 105), pg_str);
            pb.line(MARGIN_X + 24.0, ty - 5.0, MARGIN_X + CONTENT_W - 14.0, ty - 5.0, rgb(241, 245, 249), 0.5);
            ty -= 23.0;
        }

        // User notes box at bottom of cover if provided, else methodology box
        let box_y = 54.0;
        let box_h = 84.0;
        pb.fill_rect(MARGIN_X, box_y, CONTENT_W, box_h, rgb(248, 250, 252));
        pb.stroke_rect(MARGIN_X, box_y, CONTENT_W, box_h, rgb(226, 232, 240), 0.7);
        pb.fill_rect(MARGIN_X, box_y, 4.0, box_h, rgb(15, 23, 42));

        let (box_title, box_body) = if !opts.notes.trim().is_empty() {
            ("Merknader fra analytiker", opts.notes.as_str())
        } else {
            (
                "Metodikk & Evalueringsmetrikker",
                "Alle algoritmer er kjørt på identisk forbehandlet datamatrise for direkte sammenlignbarhet. Silhouette [-1, 1] (høyere er bedre) måler klyngeseparasjon og kohesjon; Davies-Bouldin (lavere er bedre) måler overlapp mellom klynger; Calinski-Harabasz (høyere er bedre) måler variansforholdet mellom og innad i klynger. Støypunkter (-1) visualiseres som grå kryss i 2D PCA-projeksjonen.",
            )
        };
        pb.text(MARGIN_X + 12.0, box_y + box_h - 17.0, 9.0, FontKind::Bold, rgb(15, 23, 42), box_title);
        let note_lines = wrap_text(box_body, CONTENT_W - 24.0, 8.2, false, 4);
        let mut ny = box_y + box_h - 32.0;
        for ln in note_lines {
            pb.text(MARGIN_X + 12.0, ny, 8.2, FontKind::Regular, rgb(71, 85, 105), &ln);
            ny -= 12.5;
        }

        pages[0].content = pb;
    }

    // Add consistent headers and footers to all non-cover pages
    for (idx, pg) in pages.iter_mut().enumerate() {
        if !pg.is_cover {
            draw_chrome(&mut pg.content, &pg.header_right, idx + 1, total_pages);
        }
    }

    // =========================================================================
    // WRITE COMPRESSED PDF 1.7 WITH BOOKMARKS (OUTLINE) VIA pdf-writer
    // =========================================================================
    let mut pdf = Pdf::new();
    let mut next_id = 1i32;
    let mut alloc = || {
        let r = Ref::new(next_id);
        next_id += 1;
        r
    };

    let catalog_id = alloc();
    let pages_id = alloc();
    let outline_id = alloc();
    let font_reg_id = alloc();
    let font_bold_id = alloc();
    let font_ital_id = alloc();
    let font_mono_id = alloc();

    let page_ids: Vec<Ref> = (0..total_pages).map(|_| alloc()).collect();
    let content_ids: Vec<Ref> = (0..total_pages).map(|_| alloc()).collect();
    let bm_ids: Vec<Ref> = (0..bookmarks.len()).map(|_| alloc()).collect();

    // Catalog
    {
        let mut cat = pdf.catalog(catalog_id);
        cat.pages(pages_id);
        if !bm_ids.is_empty() {
            cat.outlines(outline_id);
        }
    }

    // Fonts (Standard PDF Type1 with WinAnsiEncoding -> full Norwegian æøåÆØÅ!)
    for (fid, base) in [
        (font_reg_id, Name(b"Helvetica")),
        (font_bold_id, Name(b"Helvetica-Bold")),
        (font_ital_id, Name(b"Helvetica-Oblique")),
        (font_mono_id, Name(b"Courier")),
    ] {
        let mut f = pdf.type1_font(fid);
        f.base_font(base);
        f.encoding_predefined(Name(b"WinAnsiEncoding"));
    }

    // Pages root
    pdf.pages(pages_id)
        .kids(page_ids.iter().copied())
        .count(total_pages as i32);

    // Individual pages + FlateDecode-compressed content streams
    for (i, built) in pages.into_iter().enumerate() {
        let raw_stream = built.content.finish();
        let mut encoder = ZlibEncoder::new(Vec::new(), Compression::fast());
        encoder
            .write_all(&raw_stream)
            .map_err(|e| format!("PDF komprimeringsfeil: {}", e))?;
        let compressed = encoder
            .finish()
            .map_err(|e| format!("PDF komprimeringsfeil: {}", e))?;

        let mut pg = pdf.page(page_ids[i]);
        pg.parent(pages_id);
        pg.media_box(Rect::new(0.0, 0.0, PAGE_W, PAGE_H));
        pg.contents(content_ids[i]);
        let mut res = pg.resources();
        let mut fonts = res.fonts();
        fonts.pair(FontKind::Regular.pdf_name(), font_reg_id);
        fonts.pair(FontKind::Bold.pdf_name(), font_bold_id);
        fonts.pair(FontKind::Italic.pdf_name(), font_ital_id);
        fonts.pair(FontKind::Mono.pdf_name(), font_mono_id);
        fonts.finish();
        res.finish();
        pg.finish();

        pdf.stream(content_ids[i], &compressed)
            .filter(Filter::FlateDecode);
    }

    // PDF Bookmarks (Outline)
    if !bm_ids.is_empty() {
        pdf.outline(outline_id)
            .first(bm_ids[0])
            .last(*bm_ids.last().unwrap())
            .count(bm_ids.len() as i32);
        for (i, (title, pg_idx)) in bookmarks.iter().enumerate() {
            let mut item = pdf.outline_item(bm_ids[i]);
            item.parent(outline_id);
            item.title(TextStr(title));
            if i > 0 {
                item.prev(bm_ids[i - 1]);
            }
            if i + 1 < bm_ids.len() {
                item.next(bm_ids[i + 1]);
            }
            item.dest().page(page_ids[*pg_idx]).xyz(0.0, PAGE_H, None);
        }
    }

    let bytes = pdf.finish();
    fs::write(out_path, bytes)
        .map_err(|e| format!("Kunne ikke lagre PDF til '{}': {}", out_path.display(), e))?;

    Ok(total_pages)
}

fn draw_algorithm_card(
    pb: &mut PageBuilder,
    x: f32,
    y: f32,
    w: f32,
    h: f32,
    r: &AlgorithmRunResult,
    data: &PreparedData,
    bounds: (f64, f64, f64, f64),
    accent: (f32, f32, f32),
) {
    // Card container
    pb.fill_rect(x, y, w, h, rgb(255, 255, 255));
    pb.stroke_rect(x, y, w, h, rgb(203, 213, 225), 0.65);
    pb.fill_rect(x, y + h - 3.5, w, 3.5, accent);

    // Header: #ID badge + Name + k badge
    let hdr_y = y + h - 18.0;
    pb.fill_rect(x + 7.0, hdr_y - 3.0, 26.0, 12.5, accent);
    pb.text(x + 9.5, hdr_y, 7.2, FontKind::Bold, rgb(255, 255, 255), &format!("#{:03}", r.id));

    let k_pill = format!("k={}", r.metrics.n_clusters);
    let k_pill_w = 24.0;
    pb.fill_rect(x + w - k_pill_w - 7.0, hdr_y - 3.0, k_pill_w, 12.5, rgb(241, 245, 249));
    pb.stroke_rect(x + w - k_pill_w - 7.0, hdr_y - 3.0, k_pill_w, 12.5, rgb(203, 213, 225), 0.45);
    pb.text(x + w - k_pill_w - 3.5, hdr_y, 7.2, FontKind::Bold, rgb(30, 41, 59), &k_pill);

    let max_name_w = w - 74.0;
    let short_name = truncate_to_width(&r.name, max_name_w, 8.6, true);
    pb.text(x + 37.0, hdr_y, 8.6, FontKind::Bold, rgb(15, 23, 42), &short_name);

    // 2-line Norwegian description directly under header
    let desc_lines = wrap_text(&r.description_no, w - 14.0, 7.0, false, 2);
    let mut dy = hdr_y - 12.5;
    for ln in &desc_lines {
        pb.text(x + 7.0, dy, 7.0, FontKind::Regular, rgb(71, 85, 105), ln);
        dy -= 8.8;
    }

    // Middle body: Left = 2D PCA Scatter Plot, Right = Metrics table
    let body_top = hdr_y - 31.0;
    let plot_w = (w * 0.49).min(126.0);
    let plot_h = (body_top - (y + 34.0)).clamp(68.0, 118.0);
    let plot_y = body_top - plot_h;
    pb.draw_scatter_plot(x + 7.0, plot_y, plot_w, plot_h, &data.pca_2d, bounds, &r.labels, 1.85);

    // Right metrics box
    let mx = x + 7.0 + plot_w + 7.0;
    let mw = (x + w - 7.0) - mx;
    pb.fill_rect(mx, plot_y, mw, plot_h, rgb(248, 250, 252));
    pb.stroke_rect(mx, plot_y, mw, plot_h, rgb(226, 232, 240), 0.45);

    let mut my = plot_y + plot_h - 11.5;
    let metric_rows = [
        ("Silhouette:", fmt_opt(r.metrics.silhouette, 3)),
        ("Davies-B.:", fmt_opt(r.metrics.davies_bouldin, 2)),
        ("Calinski-H.:", fmt_opt(r.metrics.calinski_harabasz, 1)),
        ("Dunn-ind.:", fmt_opt(r.metrics.dunn_index, 3)),
        ("Støy:", format!("{} ({:.0}%)", r.metrics.n_noise, r.metrics.noise_ratio * 100.0)),
        ("Tid:", format!("{:.1} ms", r.duration_ms)),
    ];
    for (idx_m, (lbl, val)) in metric_rows.iter().enumerate() {
        if my < plot_y + 12.0 {
            break;
        }
        pb.text(mx + 4.5, my, 6.9, FontKind::Regular, rgb(100, 116, 139), lbl);
        let val_c = if idx_m == 0 {
            match r.metrics.silhouette {
                Some(s) if s >= 0.45 => rgb(5, 150, 105),
                Some(s) if s >= 0.20 => rgb(217, 119, 6),
                _ => rgb(71, 85, 105),
            }
        } else {
            rgb(15, 23, 42)
        };
        pb.text_right(mx + mw - 4.5, my, 7.0, FontKind::Bold, val_c, val);
        my -= 12.2;
    }

    // ARI line if ground truth exists and space permits
    if let Some(ari) = r.metrics.adjusted_rand_index {
        if my >= plot_y + 4.0 {
            pb.text(mx + 4.5, my, 6.9, FontKind::Regular, rgb(100, 116, 139), "ARI (fasit):");
            pb.text_right(mx + mw - 4.5, my, 7.0, FontKind::Bold, rgb(37, 99, 235), &format!("{:.3}", ari));
        }
    } else if let Some(sil) = r.metrics.silhouette {
        // Silhouette visual quality bar at bottom of metrics box
        let bar_y = plot_y + 5.0;
        let bar_w = mw - 9.0;
        pb.fill_rect(mx + 4.5, bar_y, bar_w, 4.5, rgb(226, 232, 240));
        let fill_frac = ((sil + 1.0) * 0.5).clamp(0.0, 1.0) as f32;
        let bar_c = if sil >= 0.45 { rgb(5, 150, 105) } else if sil >= 0.2 { rgb(217, 119, 6) } else { rgb(220, 38, 38) };
        pb.fill_rect(mx + 4.5, bar_y, bar_w * fill_frac, 4.5, bar_c);
    }

    // Cluster size proportion stacked bar directly below scatter plot
    let prop_y = y + 21.0;
    let prop_w = w - 14.0;
    let total_pts = data.n_samples.max(1) as f32;
    let mut cur_x = x + 7.0;
    for (cid, &sz) in r.metrics.cluster_sizes.iter().enumerate() {
        let seg_w = (sz as f32 / total_pts) * prop_w;
        if seg_w > 0.4 {
            pb.fill_rect(cur_x, prop_y, seg_w, 5.0, cluster_rgb(cid as i32));
            cur_x += seg_w;
        }
    }
    if r.metrics.n_noise > 0 {
        let seg_w = (r.metrics.n_noise as f32 / total_pts) * prop_w;
        pb.fill_rect(cur_x, prop_y, seg_w, 5.0, rgb(148, 163, 184));
    }

    // Bottom diagnostic note & parameter line
    let foot = if let Some(err) = &r.error {
        format!("Feil: {}", err)
    } else if !r.note.is_empty() {
        format!("{}  |  {}", r.params_summary, r.note)
    } else {
        r.params_summary.clone()
    };
    let foot_s = truncate_to_width(&foot, w - 14.0, 6.7, false);
    pb.text(x + 7.0, y + 8.5, 6.7, FontKind::Italic, rgb(100, 116, 139), &foot_s);
}
