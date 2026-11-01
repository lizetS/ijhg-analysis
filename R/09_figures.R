# ===========================================================================
# 09_figures.R -- manuscript Figures 3 and 4, plus the MCMC diagnostic
# figures requested by Reviewer 2 (point 6).
#
# Figure 1 (epidemic curve) and Figure 2 (descriptive maps) are produced from
# the same objects but are not recreated here; see docs/MANUSCRIPT_MAP.md.
#
# Outputs: output/figures/fig3_RR_<outcome>.pdf
#          output/figures/fig4_hotspots.pdf
#          output/figures/figS2_trace_acf.pdf
#          output/figures/figS3_diagnostics.pdf
# ===========================================================================


source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

log_msg("--- 09_figures ---")
dat <- read_rds("analysis_dataset.rds")
ar  <- utils::read.csv(file.path(DIR$tables, "area_rr_summary.csv"),
                       stringsAsFactors = FALSE)
dg  <- utils::read.csv(file.path(DIR$tables, "diagnostics.csv"),
                       stringsAsFactors = FALSE)

geom_only <- dat %>%
  filter(covid_variant == WAVES[1]) %>%
  select(area_code = Etiqueta, ASalud, geometry)

# ---- Figure 3: adjusted RR maps -------------------------------------------
for (oc in unique(ar$outcome)) {
  d <- geom_only %>% left_join(ar %>% filter(outcome == oc), by = "area_code")
  p <- ggplot(d) +
    geom_sf(aes(fill = rr_mean), colour = "white", linewidth = 0.15) +
    facet_wrap(~ factor(wave, levels = WAVE_LABELS[WAVES]), nrow = 2) +
    scale_fill_viridis_c(option = "magma", direction = -1,
                         name = "Adjusted\nRR") +
    labs(title = sprintf("Adjusted relative risk — %s", oc)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold"))
  path <- file.path(DIR$figures, sprintf("fig3_RR_%s.pdf", oc))
  ggsave(path, p, width = 10, height = 6.5)
  log_msg("wrote ", path)
}

# ---- Figure 4: hotspot persistence ----------------------------------------
hs <- geom_only %>%
  left_join(utils::read.csv(file.path(DIR$tables, "hotspot_classification.csv"),
                            stringsAsFactors = FALSE) %>%
              filter(rule == "posterior mean RR > 1.5"),
            by = "area_code")

if (nrow(hs)) {
  p <- ggplot(hs) +
    geom_sf(aes(fill = factor(category,
                              levels = c("Never", "Occasional", "Frequent", "Persistent"))),
            colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = c(Never = "#2166ac", Occasional = "#92c5de",
                                 Frequent = "#f4a582", Persistent = "#b2182b"),
                      name = "Hotspot type") +
    labs(title = "Hotspot persistence across waves") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank())
  path <- file.path(DIR$figures, "fig4_hotspots.pdf")
  ggsave(path, p, width = 7, height = 5.5)
  log_msg("wrote ", path)
}

# ---- Supplement: trace and ACF --------------------------------------------
models <- list.files(DIR$models, pattern = "\\.rds$", full.names = TRUE)
models <- models[grepl("^incidence_", basename(models))]

if (length(models)) {
  traces <- bind_rows(lapply(models, function(f) {
    o <- readRDS(f)
    tr <- as.data.frame(o$trace)
    tr$iteration <- seq_len(nrow(tr))
    tr$model <- paste(o$outcome, WAVE_LABELS[[o$wave]])
    tr %>% tidyr::pivot_longer(-c(iteration, model, chain),
                               names_to = "parameter", values_to = "value")
  }))

  sel <- traces %>% filter(parameter %in% c("density", "hacin"))

  p <- ggplot(sel, aes(iteration, value, colour = factor(chain))) +
    geom_line(linewidth = 0.25, alpha = 0.7) +
    facet_grid(parameter ~ model, scales = "free_y") +
    labs(title = "Trace plots (thinned)", x = "iteration", y = "value") +
    theme_minimal(base_size = 7) +
    theme(legend.position = "none")
  path <- file.path(DIR$figures, "figS2_trace_acf.pdf")
  ggsave(path, p, width = 12, height = 7)
  log_msg("wrote ", path)
}

# ---- Supplement: Rhat and ESS ---------------------------------------------
dg_ok <- dg %>% filter(status == "ok")
if (nrow(dg_ok)) {
  p1 <- ggplot(dg_ok, aes(x = rhat_max)) +
    geom_histogram(bins = 15, fill = "#2166ac") +
    geom_vline(xintercept = 1.05, linetype = 2, colour = "red") +
    labs(title = "Maximum R-hat per model", x = "max R-hat", y = "models") +
    theme_minimal(base_size = 9)
  p2 <- ggplot(dg_ok, aes(x = ess_min)) +
    geom_histogram(bins = 15, fill = "#1b7837") +
    geom_vline(xintercept = 1000, linetype = 2, colour = "red") +
    labs(title = "Minimum effective sample size per model",
         x = "min ESS", y = "models") +
    theme_minimal(base_size = 9)
  path <- file.path(DIR$figures, "figS3_diagnostics.pdf")
  ggsave(path, gridExtra::arrangeGrob(p1, p2, ncol = 2), width = 9, height = 3.8)
  log_msg("wrote ", path)
}

log_msg("09_figures complete")
