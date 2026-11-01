# ===========================================================================
# 19_figure_S1.R -- Figure S1: maps of the area-level covariates.
#
# Nine categorical maps in the layout of the submitted figure. Classes are
# quantile classes (quintiles for five-class maps, terciles for three-class
# maps), computed from the 82 area values with a rank-based split so that ties
# cannot produce empty classes. Two maps keep fixed cut-points, as in the
# submitted legends, because quantiles are meaningless there:
#   vaccination (Delta wave): <30 %, 30-49 %, 50-70 %  (almost binary: it
#                             follows the vaccination blocks)
#   population served by CMF: <1250, 1250-2500, >2500
# Note: the submitted classes were not exactly equal-sized for education,
# ageing, crowding and shanty towns, so quantile classes will not reproduce
# the class membership of the submitted figure for those maps.
# Shanty towns are mapped as a count per area (CantDBprec_N), as in the
# submitted figure.
#
# Outputs: output/figures/figS1_covariate_maps.pdf / .png
#          output/tables/figS1_class_cutpoints.csv
# ===========================================================================

source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))
log_msg("--- 19_figure_S1 ---")

dat <- read_rds("analysis_dataset.rds")
vacc <- dat %>% st_drop_geometry() %>%
  filter(covid_variant == "Delta") %>%
  select(Etiqueta, vacc_delta = cobertura_vac)
base <- dat %>% filter(covid_variant == WAVES[1]) %>%
  left_join(vacc, by = "Etiqueta")
assert(nrow(base) == 82, "expected 82 areas")
assert(!anyNA(base$vacc_delta), "vaccination coverage missing for some areas")

lab5 <- c("Very Low", "Low", "Medium", "High", "Very High")
lab3 <- c("Low", "Medium", "High")
pal  <- function(name, n, rev = FALSE) {
  p <- scales::brewer_pal(palette = name)(max(n, 3))[seq_len(n)]
  if (rev) rev(p) else p
}

specs <- list(
  list(var = "density", title = "Population Density", labels = lab5, colours = pal("RdPu", 5)),
  list(var = "GradE", title = "Grade of Education", labels = lab3, colours = pal("Reds", 3)),
  list(var = "ageing_index_cub", title = "Aging Index", labels = lab5, colours = pal("YlOrRd", 5)),
  list(var = "Dis_x_CMF_A", title = "Distance to CMF", labels = lab5, colours = pal("BuGn", 5, rev = TRUE)),
  list(var = "PobAt_CMF_A", title = "Population served by CMF",
       labels = c("<1250", "1250-2500", ">2500"), colours = pal("Greens", 3),
       breaks = c(-Inf, 1250, 2500, Inf)),
  list(var = "vacc_delta", title = "Vaccination: Delta Wave",
       labels = c("<30%", "30%-49%", "50%-70%"), colours = pal("YlGn", 3),
       breaks = c(-Inf, 30, 50, Inf)),
  list(var = "hacin", title = "Overcrowding Index", labels = lab3, colours = pal("YlGnBu", 3)),
  list(var = "CantDBprec_N", title = "Shanty Towns", labels = lab3,
       colours = c("#CCEBC5", "#B3CDE3", "#FBB4AE")),
  list(var = "casasAlq", title = "Rental Houses", labels = lab5, colours = pal("Blues", 5))
)
for (s in specs) assert(s$var %in% names(base), paste("missing variable:", s$var))

classify <- function(x, s) {
  if (!is.null(s$breaks)) {
    cut(x, breaks = s$breaks, labels = s$labels, include.lowest = TRUE, right = FALSE)
  } else {
    cut(rank(x, ties.method = "first"), breaks = length(s$labels), labels = s$labels)
  }
}

cuts <- list(); maps <- list()
for (s in specs) {
  d <- base
  d$cls <- classify(d[[s$var]], s)
  cuts[[s$var]] <- d %>% st_drop_geometry() %>% group_by(class = cls) %>%
    summarise(n_areas = n(), min = min(.data[[s$var]]), max = max(.data[[s$var]]),
              .groups = "drop") %>% mutate(variable = s$var, .before = 1)
  maps[[s$var]] <- ggplot(d) +
    geom_sf(aes(fill = cls), colour = "white", linewidth = 0.15) +
    scale_fill_manual(values = s$colours, name = "Category", drop = FALSE) +
    labs(title = s$title) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(hjust = 0, size = 10),
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 7), legend.title = element_text(size = 7))
}

write_csv_audit(bind_rows(cuts), "figS1_class_cutpoints.csv")

for (ext in c("pdf", "png")) {
  path <- file.path(DIR$figures, paste0("figS1_covariate_maps.", ext))
  if (ext == "pdf") pdf(path, width = 11, height = 8.5) else png(path, width = 11, height = 8.5, units = "in", res = 300)
  gridExtra::grid.arrange(grobs = unname(maps), ncol = 3)
  dev.off()
  log_msg("wrote ", path)
}
log_msg("19_figure_S1 complete")
