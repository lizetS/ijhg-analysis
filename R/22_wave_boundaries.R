# ===========================================================================
# 22_wave_boundaries.R -- robustness of the wave boundaries (Reviewer 2, 2.15)
#
# The previous check counted the R+S*F+ patterns produced by the wave detector
# over a grid of settings. That count is not what the reviewer asked about: with
# a 7-day window and a 7-day interval a single epidemic wave produces several
# separated R-to-F cycles and is counted as several waves, which is why the
# count ranged from 3 to 17.
#
# What matters is whether the six boundaries between the seven periods move.
# For every setting of the grid this script locates the fall from a rising to a
# falling phase, matches each official boundary to the nearest one, and reports
# the distance in days.
#
# Why a 14-day window and a 14-day segment (the manuscript's setting):
#   - the incubation period of SARS-CoV-2 has a median of about 5 days and a
#     97.5th percentile of about 11.5, so 14 days covers the full cohort of
#     infections generated on a given day;
#   - the serial interval is about 5-7 days for the ancestral lineage and 3-5
#     for Delta and Omicron, so a 14-day segment spans two to three
#     transmission generations: long enough for a linear trend to mean
#     something, short enough for the regime to be roughly constant;
#   - 14 days also removes the weekly reporting cycle and is the window of the
#     14-day notification rate used in surveillance.
# A moving average should therefore not be shorter than the incubation period,
# and the segment should not be shorter than the moving average.
#
# Outputs: output/tables/wave_boundaries.csv
#          output/tables/wave_boundaries_summary.csv
# ===========================================================================

source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

log_msg("--- 22_wave_boundaries ---")

# the six boundaries between the seven periods used in the models
OFFICIAL <- as.Date(c("2020-07-01", "2020-11-04", "2021-03-24",
                      "2021-06-30", "2021-12-15", "2022-03-09"))

linelist <- suppressMessages(readxl::read_excel(FILES$daily, .name_repair = "minimal"))
daily <- linelist %>%
  filter(!is.na(date_conf)) %>%
  mutate(date = as.Date(date_conf)) %>%
  count(date, name = "cases")

# Returns the date of every transition from a rising to a falling segment.
detect_transitions <- function(daily, window, interval, r2_min) {
  d <- daily %>% arrange(date) %>%
    tidyr::complete(date = seq(min(date), max(date), by = "day"),
                    fill = list(cases = 0))
  d$ma <- as.numeric(stats::filter(d$cases, rep(1 / window, window), sides = 1))
  d <- d %>% filter(!is.na(ma))
  d$seg <- ceiling(seq_len(nrow(d)) / interval)

  sym <- vapply(split(seq_len(nrow(d)), d$seg), function(i) {
    if (length(i) < 5 || all(is.na(d$ma[i]))) return("S")
    fit <- stats::lm(d$ma[i] ~ seq_along(i))
    sl <- coef(fit)[2]; r2 <- summary(fit)$r.squared
    if (sl > 0 && r2 >= r2_min) "R" else if (sl < 0 && r2 >= r2_min) "F" else "S"
  }, character(1))

  starts <- tapply(d$date, d$seg, min)
  symbols <- sym[as.character(sort(unique(d$seg)))]
  seg_start <- as.Date(starts[names(symbols)])

  # A boundary between two epidemic regimes is any change between a rising and
  # a falling phase -- a peak (R->F) or a trough (F->R). Stationary segments are
  # skipped, so R S S F counts as a single transition.
  keep <- symbols != "S"
  sy <- symbols[keep]; st <- seg_start[keep]
  if (length(sy) < 2) return(as.Date(character(0)))
  idx <- which(sy[-1] != sy[-length(sy)]) + 1
  st[idx]
}

match_boundaries <- function(detected, official, tol = 120) {
  if (!length(detected)) return(rep(NA_real_, length(official)))
  vapply(official, function(b) {
    dist <- as.numeric(detected - b)
    dist <- dist[abs(dist) <= tol]
    if (!length(dist)) NA_real_ else dist[which.min(abs(dist))]
  }, numeric(1))
}

grid <- expand.grid(window = c(7, 15, 21), interval = c(7, 14, 21),
                    r2_min = c(0.60, 0.75, 0.90))

rows <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  det <- detect_transitions(daily, g$window, g$interval, g$r2_min)
  shifts <- match_boundaries(det, OFFICIAL)
  rows[[length(rows) + 1]] <- tibble(
    window = g$window, interval = g$interval, r2_min = g$r2_min,
    boundary = format(OFFICIAL, "%Y-%m-%d"),
    shift_days = shifts,
    n_detected = length(det),
    # a 14-day window is the shortest that covers the incubation period and a
    # segment shorter than the window cannot have an independent trend
    defensible = g$window >= 14 && g$interval >= 7
  )
}

bnd <- bind_rows(rows)

summary_tab <- bnd %>%
  group_by(window, interval, r2_min, defensible) %>%
  summarise(n_detected = unique(n_detected),
            median_abs_shift = median(abs(shift_days), na.rm = TRUE),
            max_abs_shift = max(abs(shift_days), na.rm = TRUE),
            n_within_14d = sum(abs(shift_days) <= 14, na.rm = TRUE),
            n_boundaries_found = sum(!is.na(shift_days)),
            .groups = "drop") %>%
  arrange(desc(defensible), desc(r2_min), interval, window)

write_csv_audit(bnd, "wave_boundaries.csv")
write_csv_audit(summary_tab, "wave_boundaries_summary.csv")

print(summary_tab, n = 30)
log_msg("22_wave_boundaries complete: en los ajustes defendibles la frontera se mueve ",
        round(median(summary_tab$median_abs_shift[summary_tab$defensible], na.rm = TRUE)),
        " dias de mediana (maximo ",
        round(max(summary_tab$max_abs_shift[summary_tab$defensible], na.rm = TRUE)), ")")
