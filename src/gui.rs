use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::atomic::Ordering;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use eframe::egui;
use egui::{Color32, CornerRadius, Margin, Pos2, Rect, RichText, Stroke, StrokeKind, Vec2};
use egui_extras::{Column, TableBuilder};

use crate::algorithms::{
    run_selected_algorithms, AlgorithmRunResult, Category, QualityPreset, RunConfig, RunProgress,
};
use crate::catalog::ALGORITHM_CATALOG;
use crate::data::{
    prepare_dataset, MissingStrategy, PreprocessConfig, PreparedData, RawCsvDataset, ScalingMethod,
};
use crate::metrics::adjusted_rand_index;
use crate::pdf_report::{
    estimate_page_count, generate_pdf_report, PdfReportOptions, CLUSTER_PALETTE_U8,
};

#[derive(Clone, Copy, PartialEq, Eq)]
enum ActiveTab {
    Overview,
    Algorithms,
    Compare,
    PdfExport,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SortColumn {
    Id,
    Name,
    Clusters,
    Noise,
    Silhouette,
    DaviesBouldin,
    CalinskiHarabasz,
    Duration,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum PlotAxisMode {
    Pca2D,
    RawColumns(usize, usize),
}

struct BackgroundRun {
    progress: Arc<RunProgress>,
    started: Instant,
    output: Arc<Mutex<Option<Vec<AlgorithmRunResult>>>>,
    prepared: Arc<PreparedData>,
}

pub struct ClusterOverviewApp {
    active_tab: ActiveTab,

    raw_dataset: Option<RawCsvDataset>,
    preprocess_cfg: PreprocessConfig,
    prepared_preview: Option<Arc<PreparedData>>,

    // Run parameters
    run_cfg: RunConfig,
    selected_algos: Vec<bool>, // length 200
    algo_sidebar_search: String,

    // Running state & results
    bg_run: Option<BackgroundRun>,
    results: Vec<AlgorithmRunResult>,
    last_prepared: Option<Arc<PreparedData>>,
    last_run_wall_secs: f64,

    // Results table UI state
    result_search: String,
    result_category_filter: Option<Category>,
    sort_col: SortColumn,
    sort_asc: bool,
    selected_result_id: Option<usize>,
    plot_axis_mode: PlotAxisMode,

    // Compare tab state
    compare_a_id: usize,
    compare_b_id: usize,

    // PDF Export state
    pdf_opts: PdfReportOptions,
    last_pdf_path: Option<PathBuf>,
    last_pdf_pages: Option<usize>,

    // Status banner / feedback
    status_msg: Option<(String, bool)>, // (message, is_error)
}

impl ClusterOverviewApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let mut app = Self {
            active_tab: ActiveTab::Overview,
            raw_dataset: None,
            preprocess_cfg: PreprocessConfig {
                selected_columns: vec![],
                label_column: None,
                time_column: None,
                scaling: ScalingMethod::ZScore,
                missing_strategy: MissingStrategy::MeanMode,
                max_samples: 240,
                seed: 42,
            },
            prepared_preview: None,
            run_cfg: RunConfig {
                target_k: 3,
                seed: 42,
                quality: QualityPreset::Standard,
            },
            selected_algos: vec![true; ALGORITHM_CATALOG.len()],
            algo_sidebar_search: String::new(),
            bg_run: None,
            results: Vec::new(),
            last_prepared: None,
            last_run_wall_secs: 0.0,
            result_search: String::new(),
            result_category_filter: None,
            sort_col: SortColumn::Id,
            sort_asc: true,
            selected_result_id: None,
            plot_axis_mode: PlotAxisMode::Pca2D,
            compare_a_id: 2,  // K-Means++
            compare_b_id: 31, // DBSCAN
            pdf_opts: PdfReportOptions::default(),
            last_pdf_path: None,
            last_pdf_pages: None,
            status_msg: Some((
                "Velkommen! Åpne en CSV-fil (eller dra den inn i vinduet), eller last et innebygd eksempeldatasett.".to_string(),
                false,
            )),
        };
        app.apply_theme(&cc.egui_ctx);
        // Load built-in Iris demo right away so the user can explore immediately or replace with their CSV!
        app.load_raw_dataset(RawCsvDataset::built_in_demo_iris());
        app
    }

    fn apply_theme(&self, ctx: &egui::Context) {
        let mut visuals = egui::Visuals::light();
        visuals.panel_fill = Color32::from_rgb(248, 250, 252);
        visuals.window_fill = Color32::from_rgb(255, 255, 255);
        visuals.extreme_bg_color = Color32::from_rgb(255, 255, 255);
        visuals.faint_bg_color = Color32::from_rgb(241, 245, 249);
        visuals.selection.bg_fill = Color32::from_rgb(37, 99, 235);
        visuals.widgets.noninteractive.bg_stroke = Stroke::new(1.0_f32, Color32::from_rgb(226, 232, 240));
        visuals.widgets.inactive.bg_fill = Color32::from_rgb(241, 245, 249);
        visuals.widgets.hovered.bg_fill = Color32::from_rgb(226, 232, 240);
        visuals.widgets.noninteractive.corner_radius = CornerRadius::same(6);
        visuals.widgets.inactive.corner_radius = CornerRadius::same(6);
        visuals.widgets.hovered.corner_radius = CornerRadius::same(6);
        visuals.widgets.active.corner_radius = CornerRadius::same(6);
        ctx.set_visuals(visuals);
        // Ensure all labels and text wrap cleanly inside available bounds instead of overflowing off-screen
        ctx.style_mut(|s| {
            s.wrap_mode = Some(egui::TextWrapMode::Wrap);
        });
    }

    fn load_raw_dataset(&mut self, raw: RawCsvDataset) {
        self.preprocess_cfg = PreprocessConfig::default_for(&raw);
        let name = raw.display_name.clone();
        let n_rows = raw.rows.len();
        let n_cols = raw.columns.len();
        self.raw_dataset = Some(raw);
        self.results.clear();
        self.selected_result_id = None;
        self.last_pdf_path = None;
        self.plot_axis_mode = PlotAxisMode::Pca2D;
        self.rebuild_preview();
        if let Some(prep) = &self.prepared_preview {
            self.run_cfg.target_k = prep.suggested_k;
        }
        self.pdf_opts.subtitle = format!("Datasett: {} ({} rader, {} kolonner)", name, n_rows, n_cols);
        self.status_msg = Some((
            format!("Lastet '{}' ({} rader, {} kolonner). Klar for klyngeanalyse!", name, n_rows, n_cols),
            false,
        ));
    }

    fn rebuild_preview(&mut self) {
        if let Some(raw) = &self.raw_dataset {
            match prepare_dataset(raw, &self.preprocess_cfg) {
                Ok(prep) => {
                    self.prepared_preview = Some(Arc::new(prep));
                }
                Err(e) => {
                    self.prepared_preview = None;
                    self.status_msg = Some((e, true));
                }
            }
        }
    }

    fn open_csv_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .set_title("Velg CSV-fil for klyngeanalyse")
            .add_filter("CSV-filer", &["csv", "tsv", "txt"])
            .pick_file()
        {
            match RawCsvDataset::from_file(&path) {
                Ok(ds) => self.load_raw_dataset(ds),
                Err(e) => self.status_msg = Some((e, true)),
            }
        }
    }

    fn selected_algo_ids(&self) -> Vec<usize> {
        self.selected_algos
            .iter()
            .enumerate()
            .filter(|&(_, &on)| on)
            .map(|(i, _)| i + 1)
            .collect()
    }

    fn start_analysis(&mut self) {
        if self.bg_run.is_some() {
            return;
        }
        let ids = self.selected_algo_ids();
        if ids.is_empty() {
            self.status_msg = Some(("Velg minst én algoritme i venstre panel.".to_string(), true));
            return;
        }
        let Some(raw) = &self.raw_dataset else {
            self.status_msg = Some(("Last inn en CSV-fil først.".to_string(), true));
            return;
        };
        let prepared = match prepare_dataset(raw, &self.preprocess_cfg) {
            Ok(p) => Arc::new(p),
            Err(e) => {
                self.status_msg = Some((e, true));
                return;
            }
        };
        self.prepared_preview = Some(prepared.clone());
        let n_ids = ids.len();
        let progress = Arc::new(RunProgress::new(n_ids));
        let output = Arc::new(Mutex::new(None));
        let cfg = self.run_cfg.clone();

        let p_clone = progress.clone();
        let o_clone = output.clone();
        let d_clone = prepared.clone();
        thread::spawn(move || {
            let res = run_selected_algorithms(&d_clone, &ids, &cfg, Some(&p_clone));
            if let Ok(mut guard) = o_clone.lock() {
                *guard = Some(res);
            }
        });

        self.bg_run = Some(BackgroundRun {
            progress,
            started: Instant::now(),
            output,
            prepared,
        });
        self.status_msg = Some((
            format!("Kjører {} klyngealgoritmer parallelt på alle CPU-kjerner...", n_ids),
            false,
        ));
    }

    fn export_pdf_dialog(&mut self) {
        let (Some(raw), Some(prep)) = (&self.raw_dataset, &self.last_prepared) else {
            self.status_msg = Some(("Kjør analysen før du eksporterer PDF-rapport.".to_string(), true));
            return;
        };
        if self.results.is_empty() {
            self.status_msg = Some(("Ingen resultater å eksportere ennå.".to_string(), true));
            return;
        }
        let default_name = format!(
            "clusterOverview_rapport_{}.pdf",
            chrono::Local::now().format("%Y%m%d_%H%M")
        );
        if let Some(path) = rfd::FileDialog::new()
            .set_title("Lagre klyngeanalyserapport som PDF")
            .set_file_name(&default_name)
            .add_filter("PDF-dokument", &["pdf"])
            .save_file()
        {
            match generate_pdf_report(&path, raw, prep, &self.results, &self.pdf_opts) {
                Ok(n_pages) => {
                    self.last_pdf_path = Some(path.clone());
                    self.last_pdf_pages = Some(n_pages);
                    self.status_msg = Some((
                        format!(
                            "PDF-rapport lagret til '{}' ({} sider, {} algoritmer).",
                            path.display(),
                            n_pages,
                            self.results.len()
                        ),
                        false,
                    ));
                }
                Err(e) => {
                    self.status_msg = Some((e, true));
                }
            }
        }
    }

    fn poll_background_run(&mut self, ctx: &egui::Context) {
        if let Some(bg) = &self.bg_run {
            ctx.request_repaint_after(std::time::Duration::from_millis(60));
            let done_res = bg.output.lock().ok().and_then(|mut g| g.take());
            if let Some(res) = done_res {
                let elapsed = bg.started.elapsed().as_secs_f64();
                let prep = bg.prepared.clone();
                let was_cancelled = bg.progress.cancelled.load(Ordering::Relaxed);
                self.bg_run = None;
                self.last_run_wall_secs = elapsed;
                self.last_prepared = Some(prep);
                self.results = res;
                if let Some(first) = self.results.first() {
                    self.selected_result_id = Some(first.id);
                    self.compare_a_id = first.id;
                }
                if let Some(second) = self.results.get(1) {
                    self.compare_b_id = second.id;
                }
                self.status_msg = Some((
                    if was_cancelled {
                        format!(
                            "Analyse avbrutt etter {:.2} sek ({} algoritmer fullført).",
                            elapsed,
                            self.results.len()
                        )
                    } else {
                        format!(
                            "Ferdig! {} klyngealgoritmer fullført på {:.2} sek. Klar for inspeksjon og PDF-eksport.",
                            self.results.len(),
                            elapsed
                        )
                    },
                    false,
                ));
            }
        }
    }

    fn handle_file_drops(&mut self, ctx: &egui::Context) {
        let dropped: Vec<PathBuf> = ctx.input(|i| {
            i.raw
                .dropped_files
                .iter()
                .filter_map(|f| f.path.clone())
                .collect()
        });
        if let Some(path) = dropped.first() {
            match RawCsvDataset::from_file(path) {
                Ok(ds) => self.load_raw_dataset(ds),
                Err(e) => self.status_msg = Some((e, true)),
            }
        }
    }

    fn card_frame(&self) -> egui::Frame {
        egui::Frame::new()
            .fill(Color32::WHITE)
            .stroke(Stroke::new(1.0_f32, Color32::from_rgb(226, 232, 240)))
            .corner_radius(CornerRadius::same(8))
            .inner_margin(Margin::same(12))
    }

    fn muted_color(&self) -> Color32 {
        Color32::from_rgb(100, 116, 139)
    }

    fn heading_color(&self) -> Color32 {
        Color32::from_rgb(15, 23, 42)
    }
}

impl eframe::App for ClusterOverviewApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.handle_file_drops(ctx);
        self.poll_background_run(ctx);

        // =====================================================================
        // TOP TOOLBAR
        // =====================================================================
        egui::TopBottomPanel::top("top_bar")
            .frame(
                egui::Frame::new()
                    .fill(Color32::WHITE)
                    .stroke(Stroke::new(1.0_f32, Color32::from_rgb(226, 232, 240)))
                    .inner_margin(Margin::symmetric(12, 8)),
            )
            .show(ctx, |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.spacing_mut().item_spacing = Vec2::new(8.0, 6.0);
                    // Brand badge
                    let badge_bg = Color32::from_rgb(37, 99, 235);
                    egui::Frame::new()
                        .fill(badge_bg)
                        .corner_radius(CornerRadius::same(6))
                        .inner_margin(Margin::symmetric(8, 4))
                        .show(ui, |ui| {
                            ui.label(RichText::new("CO").strong().size(13.0).color(Color32::WHITE));
                        });
                    ui.vertical(|ui| {
                        ui.spacing_mut().item_spacing.y = 0.0;
                        ui.label(RichText::new("ClusterOverview").strong().size(15.0).color(self.heading_color()));
                        ui.label(
                            RichText::new("200 Klyngealgoritmer  •  Native Rust  •  PDF-Rapport")
                                .size(10.5)
                                .color(self.muted_color()),
                        );
                    });

                    ui.add_space(6.0);

                    // Primary toolbar buttons
                    if ui.button(RichText::new("Åpne CSV...").strong()).clicked() {
                        self.open_csv_dialog();
                    }

                    ui.menu_button("Eksempeldata", |ui| {
                        if ui.button("Iris blomsterdata (120 rader, 4 num + fasit)").clicked() {
                            self.load_raw_dataset(RawCsvDataset::built_in_demo_iris());
                            ui.close_menu();
                        }
                        if ui.button("Kundesegmentering (blandet num + kategorisk)").clicked() {
                            self.load_raw_dataset(RawCsvDataset::built_in_demo_customers());
                            ui.close_menu();
                        }
                        if ui.button("Syntetisk flergeometri (blobs + bue + ellipse)").clicked() {
                            self.load_raw_dataset(RawCsvDataset::built_in_demo_multishape());
                            ui.close_menu();
                        }
                    });

                    let sel_cnt = self.selected_algos.iter().filter(|&&b| b).count();
                    let is_running = self.bg_run.is_some();
                    let run_btn = egui::Button::new(
                        RichText::new(if is_running {
                            "Kjører...".to_string()
                        } else {
                            format!("Kjør analyse ({})", sel_cnt)
                        })
                        .strong()
                        .color(Color32::WHITE),
                    )
                    .fill(if is_running {
                        Color32::from_rgb(100, 116, 139)
                    } else {
                        Color32::from_rgb(37, 99, 235)
                    });
                    if ui.add_enabled(!is_running && self.raw_dataset.is_some(), run_btn).clicked() {
                        self.start_analysis();
                    }

                    let pdf_enabled = !self.results.is_empty() && !is_running;
                    let pdf_btn = egui::Button::new(
                        RichText::new("Eksporter PDF")
                            .strong()
                            .color(if pdf_enabled { Color32::WHITE } else { Color32::from_rgb(148, 163, 184) }),
                    )
                    .fill(if pdf_enabled {
                        Color32::from_rgb(5, 150, 105)
                    } else {
                        Color32::from_rgb(226, 232, 240)
                    });
                    if ui.add_enabled(pdf_enabled, pdf_btn).clicked() {
                        self.export_pdf_dialog();
                    }
                });
            });

        // =====================================================================
        // BOTTOM STATUS BAR
        // =====================================================================
        egui::TopBottomPanel::bottom("status_bar")
            .frame(
                egui::Frame::new()
                    .fill(Color32::from_rgb(241, 245, 249))
                    .inner_margin(Margin::symmetric(12, 5)),
            )
            .show(ctx, |ui| {
                // Allocate right-aligned status items & copyright first so they never clip,
                // then render the status message truncated cleanly inside the remaining width.
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(
                        RichText::new("© 2026 Robert Larsen")
                            .size(11.0)
                            .color(self.muted_color()),
                    );
                    if let Some(path) = &self.last_pdf_path {
                        ui.separator();
                        if ui.small_button("Åpne siste PDF").clicked() {
                            open_file_in_os(path);
                        }
                    }
                    if !self.results.is_empty() {
                        ui.separator();
                        ui.label(
                            RichText::new(format!(
                                "{} klar • {:.2}s",
                                self.results.len(),
                                self.last_run_wall_secs
                            ))
                            .size(11.0)
                            .color(Color32::from_rgb(5, 150, 105)),
                        );
                    }
                    ui.separator();
                    ui.with_layout(egui::Layout::left_to_right(egui::Align::Center), |ui| {
                        if let Some((msg, is_err)) = &self.status_msg {
                            let c = if *is_err {
                                Color32::from_rgb(220, 38, 38)
                            } else {
                                self.muted_color()
                            };
                            ui.add(egui::Label::new(RichText::new(msg).size(11.5).color(c)).truncate())
                                .on_hover_text(msg);
                        }
                    });
                });
            });

        // =====================================================================
        // LEFT CONFIGURATION SIDEBAR
        // =====================================================================
        egui::SidePanel::left("left_sidebar")
            .resizable(true)
            .default_width(295.0)
            .min_width(240.0)
            .max_width(400.0)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().auto_shrink([false; 2]).show(ui, |ui| {
                    ui.add_space(6.0);
                    self.render_sidebar_dataset_config(ui);
                    ui.add_space(8.0);
                    self.render_sidebar_run_params(ui);
                    ui.add_space(8.0);
                    self.render_sidebar_algo_picker(ui);
                    ui.add_space(10.0);
                });
            });

        // =====================================================================
        // CENTRAL WORKSPACE WITH 4 TABS
        // =====================================================================
        egui::CentralPanel::default().show(ctx, |ui| {
            // Tab Bar (wraps cleanly on narrower screens)
            ui.horizontal_wrapped(|ui| {
                let tabs = [
                    (ActiveTab::Overview, "Oversikt & Data"),
                    (ActiveTab::Algorithms, "Algoritmeresultater"),
                    (ActiveTab::Compare, "Sammenligning (A vs B)"),
                    (ActiveTab::PdfExport, "PDF-Rapport"),
                ];
                for (tab, label) in tabs {
                    let active = self.active_tab == tab;
                    let text = if tab == ActiveTab::Algorithms && !self.results.is_empty() {
                        format!("{} ({})", label, self.results.len())
                    } else {
                        label.to_string()
                    };
                    let btn = egui::Button::new(
                        RichText::new(text).strong().size(13.0).color(if active {
                            Color32::WHITE
                        } else {
                            self.heading_color()
                        }),
                    )
                    .fill(if active {
                        Color32::from_rgb(37, 99, 235)
                    } else {
                        Color32::TRANSPARENT
                    });
                    if ui.add(btn).clicked() {
                        self.active_tab = tab;
                    }
                }
            });
            ui.separator();

            // Live progress bar banner while algorithms are executing in parallel
            if let Some(bg) = &self.bg_run {
                let done = bg.progress.completed.load(Ordering::Relaxed);
                let total = bg.progress.total.load(Ordering::Relaxed).max(1);
                let frac = done as f32 / total as f32;
                let elapsed = bg.started.elapsed().as_secs_f64();

                self.card_frame().show(ui, |ui| {
                    ui.horizontal_wrapped(|ui| {
                        ui.label(
                            RichText::new(format!("Kjører klyngeanalyse: {} / {} algoritmer", done, total))
                                .strong()
                                .color(Color32::from_rgb(37, 99, 235)),
                        );
                        ui.label(RichText::new(format!("({:.1} sek)", elapsed)).color(self.muted_color()));
                        if ui.button("Avbryt").clicked() {
                            bg.progress.cancelled.store(true, Ordering::SeqCst);
                        }
                    });
                    ui.add(
                        egui::ProgressBar::new(frac)
                            .show_percentage()
                            .animate(true),
                    );
                });
                ui.add_space(6.0);
            }

            match self.active_tab {
                ActiveTab::Overview => self.render_overview_tab(ui),
                ActiveTab::Algorithms => self.render_algorithms_tab(ui),
                ActiveTab::Compare => self.render_compare_tab(ui),
                ActiveTab::PdfExport => self.render_pdf_tab(ui),
            }
        });
    }
}

impl ClusterOverviewApp {
    fn render_sidebar_dataset_config(&mut self, ui: &mut egui::Ui) {
        let mut config_changed = false;
        let heading = self.heading_color();
        let muted = self.muted_color();
        self.card_frame().show(ui, |ui| {
            ui.label(RichText::new("1. Datakilde & Forbehandling").strong().size(13.0).color(heading));
            ui.add_space(4.0);

            if let Some(raw) = &self.raw_dataset {
                let name_resp = ui.label(RichText::new(&raw.display_name).strong().size(12.0));
                if let Some(p) = &raw.path {
                    name_resp.on_hover_text(p.display().to_string());
                }
                ui.label(
                    RichText::new(format!(
                        "{} rader  •  {} kolonner  •  {}  •  '{}'",
                        raw.rows.len(),
                        raw.columns.len(),
                        raw.encoding,
                        if raw.delimiter == b'\t' { "\\t".to_string() } else { (raw.delimiter as char).to_string() }
                    ))
                    .size(11.0)
                    .color(muted),
                );
                ui.add_space(6.0);

                let combo_w = (ui.available_width() - 88.0).clamp(95.0, 210.0);

                // Scaling combo
                ui.horizontal(|ui| {
                    ui.add_sized([76.0, 18.0], egui::Label::new("Skalering:"));
                    egui::ComboBox::from_id_salt("scaling_combo")
                        .width(combo_w)
                        .truncate()
                        .selected_text(self.preprocess_cfg.scaling.label_no())
                        .show_ui(ui, |ui| {
                            for m in [ScalingMethod::ZScore, ScalingMethod::MinMax, ScalingMethod::None] {
                                if ui.selectable_value(&mut self.preprocess_cfg.scaling, m, m.label_no()).changed() {
                                    config_changed = true;
                                }
                            }
                        });
                });

                // Missing value strategy
                ui.horizontal(|ui| {
                    ui.add_sized([76.0, 18.0], egui::Label::new("Manglende:"));
                    egui::ComboBox::from_id_salt("missing_combo")
                        .width(combo_w)
                        .truncate()
                        .selected_text(self.preprocess_cfg.missing_strategy.label_no())
                        .show_ui(ui, |ui| {
                            for m in [MissingStrategy::MeanMode, MissingStrategy::DropRow] {
                                if ui.selectable_value(&mut self.preprocess_cfg.missing_strategy, m, m.label_no()).changed() {
                                    config_changed = true;
                                }
                            }
                        });
                });

                // Ground-truth label column (for ARI)
                ui.horizontal(|ui| {
                    ui.add_sized([76.0, 18.0], egui::Label::new("Fasit (ARI):"));
                    let cur_lbl = self
                        .preprocess_cfg
                        .label_column
                        .and_then(|i| raw.columns.get(i))
                        .map(|c| c.name.as_str())
                        .unwrap_or("Ingen");
                    egui::ComboBox::from_id_salt("label_col_combo")
                        .width(combo_w)
                        .truncate()
                        .selected_text(cur_lbl)
                        .show_ui(ui, |ui| {
                            if ui.selectable_value(&mut self.preprocess_cfg.label_column, None, "Ingen").changed() {
                                config_changed = true;
                            }
                            for c in &raw.columns {
                                if ui
                                    .selectable_value(
                                        &mut self.preprocess_cfg.label_column,
                                        Some(c.index),
                                        format!("{} ({} unike)", c.name, c.unique_count),
                                    )
                                    .changed()
                                {
                                    config_changed = true;
                                }
                            }
                        });
                });

                // Max sample cap slider (stacked cleanly so slider + value never clip)
                ui.add_space(2.0);
                ui.label("Maks utvalg:");
                ui.scope(|ui| {
                    ui.spacing_mut().slider_width = (ui.available_width() - 84.0).max(70.0);
                    if ui
                        .add(egui::Slider::new(&mut self.preprocess_cfg.max_samples, 50..=800).suffix(" rader"))
                        .changed()
                    {
                        config_changed = true;
                    }
                });

                // Collapsible column selector
                ui.add_space(4.0);
                egui::CollapsingHeader::new(format!(
                    "Aktive kolonner ({}/{})",
                    self.preprocess_cfg.selected_columns.iter().filter(|&&b| b).count(),
                    raw.columns.len()
                ))
                .default_open(false)
                .show(ui, |ui| {
                    ui.horizontal_wrapped(|ui| {
                        if ui.small_button("Alle").clicked() {
                            self.preprocess_cfg.selected_columns.fill(true);
                            config_changed = true;
                        }
                        if ui.small_button("Kun numeriske").clicked() {
                            for (i, c) in raw.columns.iter().enumerate() {
                                self.preprocess_cfg.selected_columns[i] =
                                    c.kind == crate::data::ColumnKind::Numeric;
                            }
                            config_changed = true;
                        }
                    });
                    for col in &raw.columns {
                        if Some(col.index) == self.preprocess_cfg.label_column {
                            ui.label(
                                RichText::new(format!("• {} [Fasitkolonne]", col.name))
                                    .size(11.0)
                                    .color(muted),
                            );
                            continue;
                        }
                        if let Some(sel) = self.preprocess_cfg.selected_columns.get_mut(col.index) {
                            if ui
                                .checkbox(sel, format!("{} ({})", col.name, col.kind.label_no()))
                                .changed()
                            {
                                config_changed = true;
                            }
                        }
                    }
                });
            } else {
                ui.label(RichText::new("Ingen CSV-fil lastet ennå.").color(muted));
            }
        });
        if config_changed {
            self.rebuild_preview();
        }
    }

    fn render_sidebar_run_params(&mut self, ui: &mut egui::Ui) {
        let heading = self.heading_color();
        let muted = self.muted_color();
        let suggested_k = self.prepared_preview.as_ref().map(|p| p.suggested_k);
        self.card_frame().show(ui, |ui| {
            ui.label(RichText::new("2. Analyseparametere").strong().size(13.0).color(heading));
            ui.add_space(4.0);

            ui.label("Antall klynger (k):");
            ui.scope(|ui| {
                ui.spacing_mut().slider_width = (ui.available_width() - 56.0).max(70.0);
                ui.add(egui::Slider::new(&mut self.run_cfg.target_k, 2..=12));
            });
            if let Some(sk) = suggested_k {
                ui.horizontal_wrapped(|ui| {
                    ui.label(
                        RichText::new(format!("Anbefalt k: {}", sk))
                            .size(11.0)
                            .color(muted),
                    );
                    if ui.small_button("Bruk anbefalt").clicked() {
                        self.run_cfg.target_k = sk;
                    }
                });
            }

            let combo_w = (ui.available_width() - 88.0).clamp(95.0, 210.0);
            ui.horizontal(|ui| {
                ui.add_sized([76.0, 18.0], egui::Label::new("Kvalitet:"));
                egui::ComboBox::from_id_salt("quality_combo")
                    .width(combo_w)
                    .truncate()
                    .selected_text(self.run_cfg.quality.label_no())
                    .show_ui(ui, |ui| {
                        for q in [QualityPreset::Fast, QualityPreset::Standard, QualityPreset::Thorough] {
                            ui.selectable_value(&mut self.run_cfg.quality, q, q.label_no());
                        }
                    });
            });

            ui.horizontal_wrapped(|ui| {
                ui.add_sized([76.0, 18.0], egui::Label::new("Frø (seed):"));
                let mut seed_i = self.run_cfg.seed as i64;
                if ui.add(egui::DragValue::new(&mut seed_i).range(0..=999999)).changed() {
                    self.run_cfg.seed = seed_i.max(0) as u64;
                    self.preprocess_cfg.seed = self.run_cfg.seed;
                }
            });
        });
    }

    fn render_sidebar_algo_picker(&mut self, ui: &mut egui::Ui) {
        self.card_frame().show(ui, |ui| {
            let total_sel = self.selected_algos.iter().filter(|&&b| b).count();
            ui.label(
                RichText::new(format!("3. Algoritmeutvalg ({} / 200)", total_sel))
                    .strong()
                    .size(13.0)
                    .color(self.heading_color()),
            );
            ui.add_space(4.0);

            ui.horizontal_wrapped(|ui| {
                if ui.small_button("Velg alle (200)").clicked() {
                    self.selected_algos.fill(true);
                }
                if ui.small_button("Fjern alle").clicked() {
                    self.selected_algos.fill(false);
                }
                if ui.small_button("Rask kjerne (18)").clicked() {
                    self.selected_algos.fill(false);
                    for id in [2, 4, 10, 12, 31, 32, 34, 51, 59, 64, 80, 81, 93, 103, 107, 151, 176, 181] {
                        self.selected_algos[id - 1] = true;
                    }
                }
            });
            ui.add_space(4.0);
            ui.add(
                egui::TextEdit::singleline(&mut self.algo_sidebar_search)
                    .hint_text("Søk blant 200 algoritmer...")
                    .desired_width(f32::INFINITY),
            );
            ui.add_space(4.0);

            let query = self.algo_sidebar_search.trim().to_lowercase();

            for cat in Category::ALL {
                let specs: Vec<_> = ALGORITHM_CATALOG
                    .iter()
                    .filter(|s| s.category == cat)
                    .collect();
                let matching: Vec<_> = if query.is_empty() {
                    specs.clone()
                } else {
                    specs
                        .iter()
                        .copied()
                        .filter(|s| {
                            s.name.to_lowercase().contains(&query)
                                || s.description_no.to_lowercase().contains(&query)
                                || s.id.to_string() == query
                        })
                        .collect()
                };
                if matching.is_empty() {
                    continue;
                }

                let cat_sel = specs.iter().filter(|s| self.selected_algos[s.id - 1]).count();
                let mut all_on = cat_sel == specs.len();
                let acc = cat.accent_rgb();

                let id = ui.make_persistent_id(format!("cat_pick_{}", cat.number()));
                let mut state = egui::collapsing_header::CollapsingState::load_with_default_open(
                    ui.ctx(),
                    id,
                    !query.is_empty(),
                );
                if !query.is_empty() {
                    state.set_open(true);
                }

                let header_txt = RichText::new(format!(
                    "{}. {} ({}/{})",
                    cat.number(),
                    cat.short_name_no(),
                    cat_sel,
                    specs.len()
                ))
                .size(11.6)
                .strong()
                .color(Color32::from_rgb(acc[0], acc[1], acc[2]));

                let header_res = ui.horizontal(|ui| {
                    state.show_toggle_button(ui, egui::collapsing_header::paint_default_icon);
                    if ui.checkbox(&mut all_on, "").changed() {
                        for s in &specs {
                            self.selected_algos[s.id - 1] = all_on;
                        }
                    }
                    let lbl = ui.add(
                        egui::Label::new(header_txt)
                            .truncate()
                            .sense(egui::Sense::click()),
                    );
                    if lbl.clicked() {
                        state.toggle(ui);
                    }
                    lbl.on_hover_text(cat.title_no());
                });

                state.show_body_indented(&header_res.response, ui, |ui| {
                    for s in matching {
                        let on = &mut self.selected_algos[s.id - 1];
                        ui.checkbox(on, RichText::new(format!("#{:03} {}", s.id, s.name)).size(11.0))
                            .on_hover_text(s.description_no);
                    }
                });
            }
        });
    }

    // =========================================================================
    // TAB 1: OVERVIEW & DATA PREVIEW
    // =========================================================================
    fn render_overview_tab(&mut self, ui: &mut egui::Ui) {
        egui::ScrollArea::vertical().auto_shrink([false; 2]).show(ui, |ui| {
            // If results exist, display Executive KPI Dashboard & Category Chart first!
            if !self.results.is_empty() {
                self.render_executive_kpis(ui);
                ui.add_space(10.0);
                self.render_category_bar_chart(ui);
                ui.add_space(10.0);
            }

            // Data preview + PCA 2D preview (side-by-side on wide screens, stacked on narrow screens)
            if let Some(raw) = &self.raw_dataset {
                let render_pca_card = |ui: &mut egui::Ui| {
                    self.card_frame().show(ui, |ui| {
                        ui.label(RichText::new("2D PCA-forhåndsvisning av datasett").strong().size(13.5).color(self.heading_color()));
                        if let Some(prep) = &self.prepared_preview {
                            ui.label(
                                RichText::new(format!(
                                    "PC1 ({:.1}%) vs PC2 ({:.1}%)  •  {} punkter  •  {} egenskaper",
                                    prep.pca_variance[0] * 100.0,
                                    prep.pca_variance[1] * 100.0,
                                    prep.n_samples,
                                    prep.n_features
                                ))
                                .size(11.0)
                                .color(self.muted_color()),
                            );
                            let lbls = prep.ground_truth.clone().unwrap_or_else(|| vec![0i32; prep.n_samples]);
                            draw_egui_scatter(ui, prep, &lbls, PlotAxisMode::Pca2D, 250.0);
                        }
                    });
                };

                let render_cols_card = |ui: &mut egui::Ui| {
                    self.card_frame().show(ui, |ui| {
                        ui.label(RichText::new("Kolonneoversikt (CSV)").strong().size(13.5).color(self.heading_color()));
                        ui.add_space(4.0);
                        egui::ScrollArea::both().max_height(265.0).id_salt("col_stats_scroll").show(ui, |ui| {
                            egui::Grid::new("col_stats_grid").striped(true).spacing([10.0, 4.0]).show(ui, |ui| {
                                ui.label(RichText::new("Kolonne").strong().size(11.5));
                                ui.label(RichText::new("Type").strong().size(11.5));
                                ui.label(RichText::new("Unike").strong().size(11.5));
                                ui.label(RichText::new("Min / Maks").strong().size(11.5));
                                ui.label(RichText::new("Snitt").strong().size(11.5));
                                ui.end_row();
                                for c in &raw.columns {
                                    ui.label(RichText::new(&c.name).size(11.0));
                                    ui.label(RichText::new(c.kind.label_no()).size(11.0).color(self.muted_color()));
                                    ui.label(RichText::new(c.unique_count.to_string()).size(11.0));
                                    let mm = match (c.min, c.max) {
                                        (Some(a), Some(b)) => format!("{:.2} .. {:.2}", a, b),
                                        _ => c.sample_values.iter().take(2).cloned().collect::<Vec<_>>().join(", "),
                                    };
                                    ui.label(RichText::new(mm).size(11.0));
                                    let mu = c.mean.map(|v| format!("{:.2}", v)).unwrap_or_else(|| "—".into());
                                    ui.label(RichText::new(mu).size(11.0));
                                    ui.end_row();
                                }
                            });
                        });
                    });
                };

                if ui.available_width() >= 720.0 {
                    ui.columns(2, |cols| {
                        render_pca_card(&mut cols[0]);
                        render_cols_card(&mut cols[1]);
                    });
                } else {
                    render_pca_card(ui);
                    ui.add_space(10.0);
                    render_cols_card(ui);
                }

                ui.add_space(10.0);

                // CSV Data Sample Table (first 30 rows)
                self.card_frame().show(ui, |ui| {
                    ui.label(
                        RichText::new(format!("Forhåndsvisning av CSV-data (viser første 30 av {} rader)", raw.rows.len()))
                            .strong()
                            .size(13.0)
                            .color(self.heading_color()),
                    );
                    ui.add_space(4.0);
                    egui::ScrollArea::both().max_height(320.0).id_salt("csv_preview_h").show(ui, |ui| {
                        egui::Grid::new("csv_preview_grid").striped(true).spacing([12.0, 3.0]).show(ui, |ui| {
                            ui.label(RichText::new("#").strong().size(11.0).color(self.muted_color()));
                            for h in &raw.headers {
                                ui.label(RichText::new(h).strong().size(11.0));
                            }
                            ui.end_row();
                            for (ri, row) in raw.rows.iter().take(30).enumerate() {
                                ui.label(RichText::new((ri + 1).to_string()).size(10.5).color(self.muted_color()));
                                for cell in row {
                                    ui.label(RichText::new(cell).size(10.5));
                                }
                                ui.end_row();
                            }
                        });
                    });
                });
            }
        });
    }

    fn render_executive_kpis(&self, ui: &mut egui::Ui) {
        let best_sil = self
            .results
            .iter()
            .filter_map(|r| r.metrics.silhouette.map(|s| (r, s)))
            .max_by(|a, b| a.1.total_cmp(&b.1));
        let best_ch = self
            .results
            .iter()
            .filter_map(|r| r.metrics.calinski_harabasz.map(|v| (r, v)))
            .max_by(|a, b| a.1.total_cmp(&b.1));
        let fastest = self
            .results
            .iter()
            .min_by(|a, b| a.duration_ms.total_cmp(&b.duration_ms));
        let avg_ms = self.results.iter().map(|r| r.duration_ms).sum::<f64>()
            / self.results.len().max(1) as f64;

        let cards = [
            (
                "ALGORITMER FULLFØRT",
                format!("{}", self.results.len()),
                format!("Total veggklokke: {:.2} sek (snitt {:.1} ms)", self.last_run_wall_secs, avg_ms),
                Color32::from_rgb(37, 99, 235),
            ),
            (
                "BESTE SILHOUETTE",
                best_sil.map(|(_, s)| format!("{:.3}", s)).unwrap_or_else(|| "—".into()),
                best_sil.map(|(r, _)| format!("#{:03} {}", r.id, r.name)).unwrap_or_default(),
                Color32::from_rgb(5, 150, 105),
            ),
            (
                "BESTE CALINSKI-HARABASZ",
                best_ch.map(|(_, v)| format!("{:.1}", v)).unwrap_or_else(|| "—".into()),
                best_ch.map(|(r, _)| format!("#{:03} {}", r.id, r.name)).unwrap_or_default(),
                Color32::from_rgb(124, 58, 237),
            ),
            (
                "RASKESTE ALGORITME",
                fastest.map(|r| format!("{:.2} ms", r.duration_ms)).unwrap_or_else(|| "—".into()),
                fastest.map(|r| format!("#{:03} {}", r.id, r.name)).unwrap_or_default(),
                Color32::from_rgb(217, 119, 6),
            ),
        ];

        let draw_kpi = |ui: &mut egui::Ui, title: &str, val: &str, sub: &str, accent: Color32| {
            self.card_frame().show(ui, |ui| {
                ui.label(RichText::new(title).strong().size(10.0).color(accent));
                ui.label(RichText::new(val).strong().size(21.0).color(self.heading_color()));
                ui.label(RichText::new(sub).size(11.0).color(self.muted_color()));
            });
        };

        if ui.available_width() >= 760.0 {
            ui.columns(4, |cols| {
                for (i, (title, val, sub, accent)) in cards.iter().enumerate() {
                    draw_kpi(&mut cols[i], title, val, sub, *accent);
                }
            });
        } else {
            ui.columns(2, |cols| {
                draw_kpi(&mut cols[0], cards[0].0, &cards[0].1, &cards[0].2, cards[0].3);
                draw_kpi(&mut cols[1], cards[1].0, &cards[1].1, &cards[1].2, cards[1].3);
            });
            ui.add_space(6.0);
            ui.columns(2, |cols| {
                draw_kpi(&mut cols[0], cards[2].0, &cards[2].1, &cards[2].2, cards[2].3);
                draw_kpi(&mut cols[1], cards[3].0, &cards[3].1, &cards[3].2, cards[3].3);
            });
        }
    }

    fn render_category_bar_chart(&self, ui: &mut egui::Ui) {
        self.card_frame().show(ui, |ui| {
            ui.label(
                RichText::new("Kategorisammenligning (Snitt og Beste Silhouette-score)")
                    .strong()
                    .size(13.5)
                    .color(self.heading_color()),
            );
            ui.add_space(6.0);
            let row_w = ui.available_width();
            for cat in Category::ALL {
                let list: Vec<_> = self.results.iter().filter(|r| r.category == cat).collect();
                if list.is_empty() {
                    continue;
                }
                let sils: Vec<f64> = list.iter().filter_map(|r| r.metrics.silhouette).collect();
                let avg_s = if sils.is_empty() { 0.0 } else { sils.iter().sum::<f64>() / sils.len() as f64 };
                let max_s = sils.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let max_clean = if max_s.is_finite() { max_s } else { 0.0 };
                let acc = cat.accent_rgb();
                let c_main = Color32::from_rgb(acc[0], acc[1], acc[2]);
                let cat_lbl = format!("{}. {} ({})", cat.number(), cat.short_name_no(), list.len());
                let score_lbl = format!("Snitt {:.3} | Maks {:.3}", avg_s, max_clean);

                let paint_bar = |ui: &mut egui::Ui, width: f32, height: f32| {
                    let (rect, _) = ui.allocate_exact_size(Vec2::new(width.max(40.0), height), egui::Sense::hover());
                    let painter = ui.painter();
                    let track = Color32::from_rgb(226, 232, 240);
                    painter.rect_filled(rect, 4.0, track);
                    let max_w = (max_clean.clamp(0.0, 1.0) as f32) * rect.width();
                    painter.rect_filled(
                        Rect::from_min_size(rect.min, Vec2::new(max_w, rect.height())),
                        4.0,
                        c_main.gamma_multiply(0.35),
                    );
                    let avg_w = (avg_s.clamp(0.0, 1.0) as f32) * rect.width();
                    painter.rect_filled(
                        Rect::from_min_size(Pos2::new(rect.min.x, rect.min.y + 3.0), Vec2::new(avg_w, rect.height() - 6.0)),
                        3.0,
                        c_main,
                    );
                };

                if row_w < 460.0 {
                    // Compact stacked row for narrow windows so right-hand scores never clip
                    ui.horizontal_wrapped(|ui| {
                        ui.label(RichText::new(&cat_lbl).size(11.5).strong());
                        ui.label(RichText::new(&score_lbl).monospace().size(10.5).color(self.muted_color()));
                    });
                    paint_bar(ui, ui.available_width(), 14.0);
                    ui.add_space(3.0);
                } else {
                    let label_w = (row_w * 0.27).clamp(120.0, 165.0);
                    let score_w = 156.0_f32;
                    let bar_w = (row_w - label_w - score_w - 22.0).max(50.0);
                    ui.horizontal(|ui| {
                        ui.add_sized(
                            [label_w, 18.0],
                            egui::Label::new(RichText::new(&cat_lbl).size(11.5).strong()).truncate(),
                        )
                        .on_hover_text(cat.title_no());
                        paint_bar(ui, bar_w, 16.0);
                        ui.add_sized(
                            [score_w, 18.0],
                            egui::Label::new(
                                RichText::new(&score_lbl)
                                    .monospace()
                                    .size(10.5)
                                    .color(self.muted_color()),
                            )
                            .truncate(),
                        );
                    });
                }
            }
        });
    }

    // =========================================================================
    // TAB 2: ALGORITHMS RESULTS TABLE + INTERACTIVE INSPECTOR
    // =========================================================================
    fn render_algorithms_tab(&mut self, ui: &mut egui::Ui) {
        if self.results.is_empty() {
            self.card_frame().show(ui, |ui| {
                ui.vertical_centered(|ui| {
                    ui.add_space(40.0);
                    ui.label(RichText::new("Ingen resultater ennå").strong().size(18.0));
                    ui.label(RichText::new("Klikk på 'Kjør analyse' i topplinjen for å kjøre valgte klyngealgoritmer.").color(self.muted_color()));
                    ui.add_space(40.0);
                });
            });
            return;
        }

        // Filter bar (wraps cleanly so buttons never overflow off-screen)
        ui.horizontal_wrapped(|ui| {
            ui.label("Søk:");
            let search_w = (ui.available_width() * 0.28).clamp(130.0, 210.0);
            ui.add(
                egui::TextEdit::singleline(&mut self.result_search)
                    .hint_text("Navn, nummer eller beskrivelse...")
                    .desired_width(search_w),
            );
            ui.label("Kategori:");
            let cat_txt = self
                .result_category_filter
                .map(|c| format!("{}. {}", c.number(), c.short_name_no()))
                .unwrap_or_else(|| "Alle 8 kategorier".to_string());
            egui::ComboBox::from_id_salt("res_cat_filter")
                .width(155.0)
                .truncate()
                .selected_text(cat_txt)
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.result_category_filter, None, "Alle 8 kategorier");
                    for c in Category::ALL {
                        ui.selectable_value(
                            &mut self.result_category_filter,
                            Some(c),
                            c.title_no(),
                        );
                    }
                });
            if ui.button("Sorter etter Silhouette").clicked() {
                self.sort_col = SortColumn::Silhouette;
                self.sort_asc = false;
            }
            if ui.button("Sorter etter ID").clicked() {
                self.sort_col = SortColumn::Id;
                self.sort_asc = true;
            }
        });
        ui.add_space(6.0);

        // Filtered & sorted indices
        let q = self.result_search.trim().to_lowercase();
        let mut indices: Vec<usize> = self
            .results
            .iter()
            .enumerate()
            .filter(|(_, r)| {
                if let Some(cf) = self.result_category_filter {
                    if r.category != cf {
                        return false;
                    }
                }
                if !q.is_empty() {
                    return r.name.to_lowercase().contains(&q)
                        || r.description_no.to_lowercase().contains(&q)
                        || r.id.to_string().contains(&q);
                }
                true
            })
            .map(|(i, _)| i)
            .collect();

        indices.sort_by(|&a, &b| {
            let ra = &self.results[a];
            let rb = &self.results[b];
            let ord = match self.sort_col {
                SortColumn::Id => ra.id.cmp(&rb.id),
                SortColumn::Name => ra.name.cmp(&rb.name),
                SortColumn::Clusters => ra.metrics.n_clusters.cmp(&rb.metrics.n_clusters),
                SortColumn::Noise => ra.metrics.noise_ratio.total_cmp(&rb.metrics.noise_ratio),
                SortColumn::Silhouette => ra
                    .metrics
                    .silhouette
                    .unwrap_or(-2.0)
                    .total_cmp(&rb.metrics.silhouette.unwrap_or(-2.0)),
                SortColumn::DaviesBouldin => ra
                    .metrics
                    .davies_bouldin
                    .unwrap_or(9999.0)
                    .total_cmp(&rb.metrics.davies_bouldin.unwrap_or(9999.0)),
                SortColumn::CalinskiHarabasz => ra
                    .metrics
                    .calinski_harabasz
                    .unwrap_or(-1.0)
                    .total_cmp(&rb.metrics.calinski_harabasz.unwrap_or(-1.0)),
                SortColumn::Duration => ra.duration_ms.total_cmp(&rb.duration_ms),
            };
            if self.sort_asc { ord } else { ord.reverse() }
        });

        // Responsive two-pane layout: side-by-side when wide enough, stacked vertically on narrow screens
        let total_w = ui.available_width();
        let avail_h = ui.available_height().max(260.0);

        if total_w >= 720.0 {
            let gap = 8.0_f32;
            let left_w = ((total_w - gap) * 0.55).floor().max(320.0);
            let right_w = (total_w - gap - left_w).max(240.0);

            ui.horizontal_top(|ui| {
                ui.allocate_ui_with_layout(
                    Vec2::new(left_w, avail_h),
                    egui::Layout::top_down(egui::Align::Min),
                    |ui| self.render_results_table_card(ui, &indices),
                );
                ui.add_space(gap);
                ui.allocate_ui_with_layout(
                    Vec2::new(right_w, avail_h),
                    egui::Layout::top_down(egui::Align::Min),
                    |ui| {
                        egui::ScrollArea::vertical().id_salt("inspector_scroll").show(ui, |ui| {
                            self.render_algorithm_inspector_card(ui);
                        });
                    },
                );
            });
        } else {
            // Narrow screen fallback: stack table above inspector so nothing ever goes off-screen
            let table_h = (avail_h * 0.45).clamp(180.0, 320.0);
            ui.allocate_ui_with_layout(
                Vec2::new(total_w, table_h),
                egui::Layout::top_down(egui::Align::Min),
                |ui| self.render_results_table_card(ui, &indices),
            );
            ui.add_space(8.0);
            egui::ScrollArea::vertical().id_salt("inspector_scroll_stacked").show(ui, |ui| {
                self.render_algorithm_inspector_card(ui);
            });
        }
    }

    fn render_results_table_card(&mut self, ui: &mut egui::Ui, indices: &[usize]) {
        self.card_frame().show(ui, |ui| {
            let mut clicked_id = None;
            let mut sort_toggle = None;

            TableBuilder::new(ui)
                .striped(true)
                .resizable(true)
                .cell_layout(egui::Layout::left_to_right(egui::Align::Center))
                .column(Column::initial(34.0).at_least(28.0).clip(true))
                .column(Column::remainder().at_least(95.0).clip(true))
                .column(Column::initial(28.0).at_least(24.0).clip(true))
                .column(Column::initial(40.0).at_least(34.0).clip(true))
                .column(Column::initial(50.0).at_least(42.0).clip(true))
                .column(Column::initial(46.0).at_least(38.0).clip(true))
                .column(Column::initial(50.0).at_least(42.0).clip(true))
                .column(Column::initial(46.0).at_least(38.0).clip(true))
                .sense(egui::Sense::click())
                .header(22.0, |mut header| {
                    let hdr_defs = [
                        ("#", SortColumn::Id),
                        ("Algoritme", SortColumn::Name),
                        ("k", SortColumn::Clusters),
                        ("Støy", SortColumn::Noise),
                        ("Sil.", SortColumn::Silhouette),
                        ("DBI", SortColumn::DaviesBouldin),
                        ("CHI", SortColumn::CalinskiHarabasz),
                        ("ms", SortColumn::Duration),
                    ];
                    for (txt, sc) in hdr_defs {
                        header.col(|ui| {
                            let arrow = if self.sort_col == sc {
                                if self.sort_asc { " ▲" } else { " ▼" }
                            } else {
                                ""
                            };
                            if ui
                                .selectable_label(false, RichText::new(format!("{}{}", txt, arrow)).strong().size(11.0))
                                .clicked()
                            {
                                sort_toggle = Some(sc);
                            }
                        });
                    }
                })
                .body(|body| {
                    body.rows(22.0, indices.len(), |mut row| {
                        let r = &self.results[indices[row.index()]];
                        let is_sel = self.selected_result_id == Some(r.id);
                        row.set_selected(is_sel);

                        row.col(|ui| {
                            let acc = r.category.accent_rgb();
                            ui.label(
                                RichText::new(format!("{:03}", r.id))
                                    .monospace()
                                    .size(11.0)
                                    .color(Color32::from_rgb(acc[0], acc[1], acc[2])),
                            );
                        });
                        row.col(|ui| {
                            ui.add(egui::Label::new(RichText::new(&r.name).size(11.5).strong()).truncate())
                                .on_hover_text(&r.name);
                        });
                        row.col(|ui| {
                            ui.label(RichText::new(r.metrics.n_clusters.to_string()).size(11.0));
                        });
                        row.col(|ui| {
                            ui.label(RichText::new(format!("{:.0}%", r.metrics.noise_ratio * 100.0)).size(11.0));
                        });
                        row.col(|ui| {
                            let txt = r.metrics.silhouette.map(|v| format!("{:.3}", v)).unwrap_or_else(|| "—".into());
                            let col = match r.metrics.silhouette {
                                Some(s) if s >= 0.45 => Color32::from_rgb(5, 150, 105),
                                Some(s) if s >= 0.20 => Color32::from_rgb(217, 119, 6),
                                _ => self.muted_color(),
                            };
                            ui.label(RichText::new(txt).strong().size(11.0).color(col));
                        });
                        row.col(|ui| {
                            let txt = r.metrics.davies_bouldin.map(|v| format!("{:.2}", v)).unwrap_or_else(|| "—".into());
                            ui.label(RichText::new(txt).size(11.0));
                        });
                        row.col(|ui| {
                            let txt = r.metrics.calinski_harabasz.map(|v| format!("{:.1}", v)).unwrap_or_else(|| "—".into());
                            ui.label(RichText::new(txt).size(11.0));
                        });
                        row.col(|ui| {
                            ui.label(RichText::new(format!("{:.1}", r.duration_ms)).monospace().size(10.5));
                        });

                        if row.response().clicked() {
                            clicked_id = Some(r.id);
                        }
                    });
                });

            if let Some(id) = clicked_id {
                self.selected_result_id = Some(id);
            }
            if let Some(sc) = sort_toggle {
                if self.sort_col == sc {
                    self.sort_asc = !self.sort_asc;
                } else {
                    self.sort_col = sc;
                    self.sort_asc = matches!(sc, SortColumn::Id | SortColumn::Name | SortColumn::DaviesBouldin | SortColumn::Duration);
                }
            }
        });
    }

    fn render_algorithm_inspector_card(&mut self, ui: &mut egui::Ui) {
        let sel_res = self
            .selected_result_id
            .and_then(|id| self.results.iter().find(|r| r.id == id))
            .or_else(|| self.results.first())
            .cloned();

        if let (Some(r), Some(prep)) = (sel_res, self.last_prepared.clone()) {
            self.card_frame().show(ui, |ui| {
                let acc = r.category.accent_rgb();
                ui.horizontal_wrapped(|ui| {
                    egui::Frame::new()
                        .fill(Color32::from_rgb(acc[0], acc[1], acc[2]))
                        .corner_radius(CornerRadius::same(4))
                        .inner_margin(Margin::symmetric(6, 2))
                        .show(ui, |ui| {
                            ui.label(RichText::new(format!("#{:03}", r.id)).strong().size(11.0).color(Color32::WHITE));
                        });
                    ui.label(RichText::new(&r.name).strong().size(15.0).color(self.heading_color()));
                });
                ui.label(RichText::new(r.category.title_no()).size(11.0).color(Color32::from_rgb(acc[0], acc[1], acc[2])));
                ui.add_space(4.0);
                ui.label(RichText::new(&r.description_no).size(12.0));
                ui.add_space(6.0);

                // Axis selector for interactive scatter plot
                ui.horizontal_wrapped(|ui| {
                    ui.label(RichText::new("Projeksjon:").size(11.5));
                    let axis_lbl = match self.plot_axis_mode {
                        PlotAxisMode::Pca2D => "PCA 2D (PC1 vs PC2)".to_string(),
                        PlotAxisMode::RawColumns(x, y) => format!(
                            "{} vs {}",
                            prep.raw_numeric_names.get(x).map(|s| s.as_str()).unwrap_or("X"),
                            prep.raw_numeric_names.get(y).map(|s| s.as_str()).unwrap_or("Y"),
                        ),
                    };
                    let axis_w = (ui.available_width() - 80.0).clamp(130.0, 240.0);
                    egui::ComboBox::from_id_salt("plot_axis_combo")
                        .width(axis_w)
                        .truncate()
                        .selected_text(axis_lbl)
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut self.plot_axis_mode, PlotAxisMode::Pca2D, "PCA 2D (PC1 vs PC2)");
                            let m = prep.raw_numeric_names.len().min(6);
                            for a in 0..m {
                                for b in (a + 1)..m {
                                    ui.selectable_value(
                                        &mut self.plot_axis_mode,
                                        PlotAxisMode::RawColumns(a, b),
                                        format!("{} vs {}", prep.raw_numeric_names[a], prep.raw_numeric_names[b]),
                                    );
                                }
                            }
                        });
                });

                ui.add_space(4.0);
                draw_egui_scatter(ui, &prep, &r.labels, self.plot_axis_mode, 250.0);
                ui.add_space(8.0);

                // Metrics grid (constrained column width so values never overflow right edge)
                let col_max = (ui.available_width() * 0.62).max(110.0);
                egui::Grid::new("inspector_metrics")
                    .striped(true)
                    .num_columns(2)
                    .max_col_width(col_max)
                    .spacing([12.0, 4.0])
                    .show(ui, |ui| {
                        let rows = [
                            ("Antall klynger (k)", format!("{}", r.metrics.n_clusters)),
                            ("Støypunkter (-1)", format!("{} ({:.1}%)", r.metrics.n_noise, r.metrics.noise_ratio * 100.0)),
                            ("Silhouette-score [-1..1]", r.metrics.silhouette.map(|v| format!("{:.4}", v)).unwrap_or_else(|| "—".into())),
                            ("Davies-Bouldin (lavere=bedre)", r.metrics.davies_bouldin.map(|v| format!("{:.4}", v)).unwrap_or_else(|| "—".into())),
                            ("Calinski-Harabasz (høyere=bedre)", r.metrics.calinski_harabasz.map(|v| format!("{:.2}", v)).unwrap_or_else(|| "—".into())),
                            ("Dunn-indeks", r.metrics.dunn_index.map(|v| format!("{:.4}", v)).unwrap_or_else(|| "—".into())),
                            ("Inertia (SSE)", format!("{:.2}", r.metrics.inertia_sse)),
                            ("Balanse-entropi [0..1]", format!("{:.3}", r.metrics.balance_entropy)),
                            ("Kjøretid", format!("{:.2} ms", r.duration_ms)),
                        ];
                        for (k, v) in rows {
                            ui.label(RichText::new(k).size(11.5).color(self.muted_color()));
                            ui.label(RichText::new(v).strong().size(11.5));
                            ui.end_row();
                        }
                        if let Some(ari) = r.metrics.adjusted_rand_index {
                            ui.label(RichText::new("Adjusted Rand Index (fasit)").size(11.5).color(self.muted_color()));
                            ui.label(RichText::new(format!("{:.4}", ari)).strong().size(11.5).color(Color32::from_rgb(37, 99, 235)));
                            ui.end_row();
                        }
                    });

                ui.add_space(6.0);
                ui.label(RichText::new("Klyngestørrelser:").strong().size(11.5));
                ui.horizontal_wrapped(|ui| {
                    for (cid, &sz) in r.metrics.cluster_sizes.iter().enumerate() {
                        let c = CLUSTER_PALETTE_U8[cid % CLUSTER_PALETTE_U8.len()];
                        egui::Frame::new()
                            .fill(Color32::from_rgb(c[0], c[1], c[2]).gamma_multiply(0.18))
                            .stroke(Stroke::new(1.0_f32, Color32::from_rgb(c[0], c[1], c[2])))
                            .corner_radius(CornerRadius::same(4))
                            .inner_margin(Margin::symmetric(6, 2))
                            .show(ui, |ui| {
                                ui.label(RichText::new(format!("C{}: {}", cid, sz)).size(10.5).strong());
                            });
                    }
                });

                ui.add_space(6.0);
                ui.label(RichText::new(format!("Parametere: {}", r.params_summary)).size(11.0).color(self.muted_color()));
                if !r.note.is_empty() {
                    ui.label(RichText::new(format!("Merknad: {}", r.note)).italics().size(11.0).color(self.muted_color()));
                }
            });
        }
    }

    // =========================================================================
    // TAB 3: SIDE-BY-SIDE COMPARISON (A vs B)
    // =========================================================================
    fn render_compare_tab(&mut self, ui: &mut egui::Ui) {
        if self.results.len() < 2 {
            self.card_frame().show(ui, |ui| {
                ui.vertical_centered(|ui| {
                    ui.add_space(36.0);
                    ui.label(RichText::new("Kjør minst 2 algoritmer for å sammenligne side-om-side.").color(self.muted_color()));
                    ui.add_space(36.0);
                });
            });
            return;
        }
        let Some(prep) = self.last_prepared.clone() else { return };

        egui::ScrollArea::vertical().show(ui, |ui| {
            // Algorithm A & B pickers (wraps cleanly on narrow screens)
            self.card_frame().show(ui, |ui| {
                let cb_w = (ui.available_width() * 0.34).clamp(150.0, 260.0);
                ui.horizontal_wrapped(|ui| {
                    ui.label(RichText::new("Algoritme A:").strong());
                    let cur_a = self.results.iter().find(|r| r.id == self.compare_a_id).unwrap_or(&self.results[0]);
                    egui::ComboBox::from_id_salt("cmp_a")
                        .width(cb_w)
                        .truncate()
                        .selected_text(format!("#{:03} {}", cur_a.id, cur_a.name))
                        .show_ui(ui, |ui| {
                            for r in &self.results {
                                ui.selectable_value(&mut self.compare_a_id, r.id, format!("#{:03} {}", r.id, r.name));
                            }
                        });

                    ui.add_space(10.0);
                    ui.label(RichText::new("Algoritme B:").strong());
                    let cur_b = self.results.iter().find(|r| r.id == self.compare_b_id).unwrap_or(&self.results[1]);
                    egui::ComboBox::from_id_salt("cmp_b")
                        .width(cb_w)
                        .truncate()
                        .selected_text(format!("#{:03} {}", cur_b.id, cur_b.name))
                        .show_ui(ui, |ui| {
                            for r in &self.results {
                                ui.selectable_value(&mut self.compare_b_id, r.id, format!("#{:03} {}", r.id, r.name));
                            }
                        });
                });
            });
            ui.add_space(8.0);

            let ra = self.results.iter().find(|r| r.id == self.compare_a_id).unwrap_or(&self.results[0]).clone();
            let rb = self.results.iter().find(|r| r.id == self.compare_b_id).unwrap_or(&self.results[1]).clone();
            let ari_ab = adjusted_rand_index(&ra.labels, &rb.labels);

            let draw_cmp_card = |ui: &mut egui::Ui, r: &AlgorithmRunResult| {
                self.card_frame().show(ui, |ui| {
                    ui.label(RichText::new(format!("#{:03} {}", r.id, r.name)).strong().size(14.0).color(self.heading_color()));
                    ui.label(RichText::new(&r.description_no).size(11.0).color(self.muted_color()));
                    ui.add_space(4.0);
                    draw_egui_scatter(ui, &prep, &r.labels, self.plot_axis_mode, 240.0);
                    ui.add_space(4.0);
                    ui.label(
                        RichText::new(format!(
                            "k={}  •  Silhouette: {}  •  DBI: {}  •  Tid: {:.1} ms",
                            r.metrics.n_clusters,
                            r.metrics.silhouette.map(|v| format!("{:.3}", v)).unwrap_or_else(|| "—".into()),
                            r.metrics.davies_bouldin.map(|v| format!("{:.2}", v)).unwrap_or_else(|| "—".into()),
                            r.duration_ms
                        ))
                        .strong()
                        .size(11.5),
                    );
                });
            };

            if ui.available_width() >= 640.0 {
                ui.columns(2, |cols| {
                    draw_cmp_card(&mut cols[0], &ra);
                    draw_cmp_card(&mut cols[1], &rb);
                });
            } else {
                draw_cmp_card(ui, &ra);
                ui.add_space(8.0);
                draw_cmp_card(ui, &rb);
            }

            ui.add_space(8.0);
            self.card_frame().show(ui, |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.label(RichText::new("Enighet mellom A og B (Adjusted Rand Index):").strong().size(13.0));
                    ui.label(
                        RichText::new(format!("{:.4}", ari_ab))
                            .strong()
                            .size(15.0)
                            .color(Color32::from_rgb(37, 99, 235)),
                    );
                    ui.label(
                        RichText::new(if ari_ab > 0.75 {
                            "(Svært høy enighet)"
                        } else if ari_ab > 0.40 {
                            "(Moderat enighet)"
                        } else {
                            "(Ulik klyngestruktur)"
                        })
                        .color(self.muted_color()),
                    );
                });
            });
        });
    }

    // =========================================================================
    // TAB 4: PDF REPORT EXPORT STUDIO
    // =========================================================================
    fn render_pdf_tab(&mut self, ui: &mut egui::Ui) {
        egui::ScrollArea::vertical().show(ui, |ui| {
            self.card_frame().show(ui, |ui| {
                ui.label(RichText::new("PDF-Rapportkonfigurasjon").strong().size(15.0).color(self.heading_color()));
                ui.label(
                    RichText::new("Eksporter en profesjonell A4 vektor-PDF med forside, datasettprofil, kategoridiagrammer, 2D PCA-spredningsplott per algoritme og klikkbare bokmerker.")
                        .size(11.5)
                        .color(self.muted_color()),
                );
                ui.add_space(10.0);

                let field_w = (ui.available_width() - 155.0).clamp(160.0, 460.0);
                egui::Grid::new("pdf_cfg_grid").num_columns(2).spacing([12.0, 8.0]).show(ui, |ui| {
                    ui.label(RichText::new("Rapporttittel:").strong());
                    ui.add(egui::TextEdit::singleline(&mut self.pdf_opts.title).desired_width(field_w));
                    ui.end_row();

                    ui.label(RichText::new("Undertittel:").strong());
                    ui.add(egui::TextEdit::singleline(&mut self.pdf_opts.subtitle).desired_width(field_w));
                    ui.end_row();

                    ui.label(RichText::new("Forfatter / Org.:").strong());
                    ui.add(egui::TextEdit::singleline(&mut self.pdf_opts.author).desired_width(field_w));
                    ui.end_row();

                    ui.label(RichText::new("Merknader (forside):").strong());
                    ui.add(
                        egui::TextEdit::multiline(&mut self.pdf_opts.notes)
                            .desired_width(field_w)
                            .desired_rows(3)
                            .hint_text("Valgfrie kommentarer eller metodebeskrivelse som vises nederst på forsiden..."),
                    );
                    ui.end_row();
                });

                ui.add_space(8.0);
                ui.checkbox(
                    &mut self.pdf_opts.include_scatter_cards,
                    "Inkluder individuelle algoritmekort med 2D PCA vektor-spredningsplott (6 kort per side — Anbefalt)",
                );
                ui.checkbox(
                    &mut self.pdf_opts.include_appendix_table,
                    "Inkluder komplett sammenligningstabell som vedlegg bakerst i rapporten",
                );

                ui.add_space(12.0);
                let est_pages = estimate_page_count(self.results.len(), &self.pdf_opts);
                ui.label(
                    RichText::new(format!(
                        "Estimert omfang: ~{} A4-sider  •  {} analyserte algoritmer  •  Vektorgrafikk + PDF-bokmerker",
                        est_pages,
                        self.results.len()
                    ))
                    .strong()
                    .color(self.muted_color()),
                );
                ui.add_space(10.0);

                ui.horizontal_wrapped(|ui| {
                    let enabled = !self.results.is_empty() && self.bg_run.is_none();
                    let btn = egui::Button::new(
                        RichText::new("Generer og lagre PDF-rapport...")
                            .strong()
                            .size(14.0)
                            .color(Color32::WHITE),
                    )
                    .fill(if enabled {
                        Color32::from_rgb(5, 150, 105)
                    } else {
                        Color32::from_rgb(148, 163, 184)
                    })
                    .min_size(Vec2::new(210.0, 36.0));
                    if ui.add_enabled(enabled, btn).clicked() {
                        self.export_pdf_dialog();
                    }

                    if let Some(path) = self.last_pdf_path.clone() {
                        if ui.button(RichText::new("Åpne generert PDF").strong()).clicked() {
                            open_file_in_os(&path);
                        }
                        if let Some(parent) = path.parent() {
                            if ui.button("Åpne mappe").clicked() {
                                open_file_in_os(parent);
                            }
                        }
                    }
                });

                if let (Some(path), Some(pgs)) = (&self.last_pdf_path, self.last_pdf_pages) {
                    ui.add_space(10.0);
                    ui.label(
                        RichText::new(format!("Siste rapport lagret: {} ({} sider)", path.display(), pgs))
                            .color(Color32::from_rgb(5, 150, 105))
                            .strong(),
                    );
                }
            });
        });
    }
}

/// Interactive native egui 2D scatter plot with gridlines, noise markers, and hover tooltips.
fn draw_egui_scatter(
    ui: &mut egui::Ui,
    prep: &PreparedData,
    labels: &[i32],
    mode: PlotAxisMode,
    height: f32,
) {
    let coords: Vec<[f64; 2]> = match mode {
        PlotAxisMode::Pca2D => prep.pca_2d.clone(),
        PlotAxisMode::RawColumns(cx, cy) => (0..prep.n_samples)
            .map(|i| {
                [
                    prep.raw_numeric_cols.get(cx).and_then(|c| c.get(i)).copied().unwrap_or(0.0),
                    prep.raw_numeric_cols.get(cy).and_then(|c| c.get(i)).copied().unwrap_or(0.0),
                ]
            })
            .collect(),
    };

    let (rect, response) = ui.allocate_exact_size(
        Vec2::new(ui.available_width().max(120.0), height),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);

    let bg = Color32::from_rgb(252, 253, 255);
    let border = Color32::from_rgb(203, 213, 225);
    let grid_c = Color32::from_rgb(236, 241, 247);

    painter.rect_filled(rect, 6.0, bg);
    painter.rect_stroke(rect, 6.0, Stroke::new(1.0_f32, border), StrokeKind::Inside);

    for f in [0.25_f32, 0.5, 0.75] {
        let gx = rect.left() + rect.width() * f;
        let gy = rect.top() + rect.height() * f;
        painter.line_segment([Pos2::new(gx, rect.top()), Pos2::new(gx, rect.bottom())], Stroke::new(1.0_f32, grid_c));
        painter.line_segment([Pos2::new(rect.left(), gy), Pos2::new(rect.right(), gy)], Stroke::new(1.0_f32, grid_c));
    }

    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for p in &coords {
        min_x = min_x.min(p[0]);
        max_x = max_x.max(p[0]);
        min_y = min_y.min(p[1]);
        max_y = max_y.max(p[1]);
    }
    let sx = (max_x - min_x).max(1e-9);
    let sy = (max_y - min_y).max(1e-9);
    let pad = 16.0_f32;
    let iw = (rect.width() - 2.0 * pad).max(10.0);
    let ih = (rect.height() - 2.0 * pad).max(10.0);

    let mut screen_pts = Vec::with_capacity(coords.len());
    for (i, p) in coords.iter().enumerate() {
        let px = rect.left() + pad + ((p[0] - min_x) / sx) as f32 * iw;
        let py = rect.bottom() - pad - ((p[1] - min_y) / sy) as f32 * ih;
        let pos = Pos2::new(px, py);
        screen_pts.push(pos);
        let lbl = labels.get(i).copied().unwrap_or(-1);
        if lbl < 0 {
            let r = 3.2;
            let st = Stroke::new(1.2_f32, Color32::from_rgb(148, 163, 184));
            painter.line_segment([Pos2::new(px - r, py - r), Pos2::new(px + r, py + r)], st);
            painter.line_segment([Pos2::new(px - r, py + r), Pos2::new(px + r, py - r)], st);
        } else {
            let c = CLUSTER_PALETTE_U8[(lbl as usize) % CLUSTER_PALETTE_U8.len()];
            painter.circle_filled(pos, 3.8, Color32::from_rgb(c[0], c[1], c[2]));
            painter.circle_stroke(pos, 3.8, Stroke::new(0.6_f32, Color32::from_black_alpha(45)));
        }
    }

    // Hover point tooltip
    if let Some(pointer) = response.hover_pos() {
        let nearest = screen_pts
            .iter()
            .enumerate()
            .map(|(i, &pt)| (i, pt.distance(pointer)))
            .min_by(|a, b| a.1.total_cmp(&b.1));
        if let Some((idx, dist)) = nearest {
            if dist <= 10.0 {
                let lbl = labels.get(idx).copied().unwrap_or(-1);
                painter.circle_stroke(screen_pts[idx], 6.5, Stroke::new(2.0_f32, Color32::from_rgb(245, 158, 11)));
                let cl_str = if lbl < 0 { "Støy (-1)".to_string() } else { format!("Klynge {}", lbl) };
                egui::show_tooltip_at_pointer(ui.ctx(), ui.layer_id(), egui::Id::new("scatter_tip"), |ui| {
                    ui.label(RichText::new(format!("Rad #{}  •  {}", prep.source_row_indices[idx] + 1, cl_str)).strong());
                    ui.label(format!("X={:.3}, Y={:.3}", coords[idx][0], coords[idx][1]));
                });
            }
        }
    }

    // Distinct cluster count watermark
    let n_cl = labels.iter().filter(|&&l| l >= 0).collect::<HashSet<_>>().len();
    let n_noise = labels.iter().filter(|&&l| l < 0).count();
    painter.text(
        Pos2::new(rect.right() - 8.0, rect.top() + 8.0),
        egui::Align2::RIGHT_TOP,
        format!("k={} | støy={}", n_cl, n_noise),
        egui::FontId::proportional(10.5),
        Color32::from_rgb(100, 116, 139),
    );
}

fn open_file_in_os(path: &std::path::Path) {
    #[cfg(target_os = "windows")]
    {
        let _ = std::process::Command::new("cmd")
            .args(["/C", "start", "", &path.to_string_lossy()])
            .spawn();
    }
    #[cfg(target_os = "macos")]
    {
        let _ = std::process::Command::new("open").arg(path).spawn();
    }
    #[cfg(target_os = "linux")]
    {
        let _ = std::process::Command::new("xdg-open").arg(path).spawn();
    }
}
