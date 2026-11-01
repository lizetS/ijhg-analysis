# ===========================================================================
# 12_figures_descriptive.R -- manuscript Figures 1 and 2.
#
#   Figure 1  daily incidence and 14-day moving averages, with 14-day linear
#             trend segments and the epidemic-wave boundaries. Requires the
#             daily case series, which is part of the restricted data; the
#             figure is skipped if the file is absent.
#
#   Figure 2  descriptive maps of incidence rate, mortality rate and case
#             fatality by health area and wave, with quintile intervals
#             standardised across waves for each indicator.
#
# Outputs: output/figures/fig1_epidemic_curve.pdf
#          output/figures/fig2_descriptive_maps.pdf
# ===========================================================================

source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

log_msg("--- figures 2 and S1 ---")
dat <- read_rds("analysis_dataset.rds")
pop_total <- sum(dat$poblacion[dat$covid_variant == WAVES[1]])

# ---------------------------------------------------------------------------
# Figure 2 -- descriptive maps
# ---------------------------------------------------------------------------
map_block <- function(outcome, var, title, palette) {
  d <- dat %>%
    mutate(value = .data[[var]]) %>%
    mutate(value = ifelse(is.infinite(value), NA, value))
  br <- quantile(d$value, probs = seq(0, 1, 0.2), na.rm = TRUE)
  br <- unique(br)
  d <- d %>% mutate(cat = cut(value, breaks = br, include.lowest = TRUE,
                              labels = FALSE))
  ggplot(d) +
    geom_sf(aes(fill = factor(cat)), colour = "white", linewidth = 0.1) +
    facet_wrap(~ factor(covid_variant, levels = WAVES, labels = WAVE_LABELS[WAVES]),
               nrow = 1) +
    scale_fill_manual(values = palette, na.value = "grey90", name = "Quintile",
                      labels = if (length(palette) == length(br) - 1) 1:(length(br) - 1) else NULL) +
    labs(title = title) +
    theme_void(base_size = 7) +
    theme(strip.text = element_text(size = 6.5, face = "bold"),
          legend.position = "right",
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6.5),
          plot.title = element_text(size = 8, face = "bold", margin = margin(b = 10)))
}

pal <- function(cols, n) colorRampPalette(cols)(n)
nq <- 5
blocks <- list(
  map_block("incidence", "incidence_rate", "Cumulative incidence rate (per 100,000)",
            pal(c("#deebf7", "#08519c"), nq)),
  map_block("mortality", "mortality_rate", "Cumulative mortality rate (per 100,000)",
            pal(c("#fee6ce", "#a63603"), nq)),
  map_block("cfr", "lethality", "Case fatality rate (%)",
            pal(c("#fcbba1", "#a50f15"), nq))
)
path <- file.path(DIR$figures, "fig2_descriptive_maps.pdf")
pdf(path, width = 14, height = 8.5)
gridExtra::grid.arrange(grobs = blocks, ncol = 1)
dev.off()
log_msg("wrote ", path)

# ---------------------------------------------------------------------------
# Figure S1 -- distribution of the key area-level risk factors
# ---------------------------------------------------------------------------
map_single <- function(var, title, cols) {
  d <- dat %>% filter(covid_variant == WAVES[1])
  br <- unique(quantile(d[[var]], probs = seq(0, 1, 0.2), na.rm = TRUE))
  d$cat <- cut(d[[var]], breaks = br, include.lowest = TRUE, labels = FALSE)
  ggplot(d) +
    geom_sf(aes(fill = factor(cat)), colour = "white", linewidth = 0.15) +
    scale_fill_manual(values = pal(cols, length(br) - 1), na.value = "grey90",
                      name = "Quintile") +
    labs(title = title) +
    theme_void(base_size = 8) +
    theme(legend.position = "right",
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 6.5),
          legend.title = element_text(size = 7),
          plot.title = element_text(size = 8.5, face = "bold"))
}

s1 <- list(
  map_single("casasAlq", "A  Rental housing density (establishments per km2)",
             c("#e5f5e0", "#00441b")),
  map_single("ageing_index_cub", "B  Ageing index (% aged 60 or over)",
             c("#deebf7", "#08306b")),
  map_single("hacin", "C  Overcrowding (residents per bedroom)",
             c("#fee6ce", "#7f2704")),
  map_single("Dis_x_CMF_A", "D  Distance to the nearest CMF (metres)",
             c("#f2f0f7", "#3f007d"))
)
path <- file.path(DIR$figures, "figS1b_continuous_maps.pdf")
pdf(path, width = 11, height = 9)
gridExtra::grid.arrange(grobs = s1, ncol = 2)
dev.off()
log_msg("wrote ", path)

log_msg("12_figures_descriptive complete")
