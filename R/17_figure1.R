# Sourced by run_all.R after config.R; when run on its own, load config so that
# FILES and DIR resolve, including any RUN_TAG sub-folder.
if (!exists("DIR")) source("config.R")

# Instalar/cargar paquetes
# The individual packages are loaded instead of the tidyverse meta-package: the
# meta-package pins a dplyr version that conflicts with the one shipped by some
# conda builds, and this script only needs dplyr and tibble from it.
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(lubridate)
library(zoo)
library(readxl)
library(patchwork)
library(ggrepel)

# ============================================================
# Leer y procesar datos
# ============================================================
datos <- read_excel(if (exists("FILES")) FILES$serie else file.path("data", "raw", "serie.xlsx"))

datos <- datos %>%
  rename(Fecha = 1, Incidencia = Incidence, Mortalidad = Mortality) %>%
  mutate(Fecha = as.Date(Fecha, origin = "1899-12-30")) %>%
  arrange(Fecha) %>%
  mutate(
    MA14_Incidencia = rollmean(Incidencia, k = 14, fill = NA, align = "right"),
    MA14_Mortalidad = rollmean(Mortalidad, k = 14, fill = NA, align = "right"),
    Periodo = floor(as.numeric(difftime(Fecha, min(Fecha), units = "days")) / 14)
  )

resumen_incidencia <- datos %>%
  group_by(Periodo) %>%
  summarise(
    FechaInicio = min(Fecha), FechaFin = max(Fecha),
    modelo = list(lm(Incidencia ~ as.numeric(Fecha))),
    .groups = "drop"
  ) %>%
  mutate(
    slope = map_dbl(modelo, ~ coef(.x)[2]),
    r2 = map_dbl(modelo, ~ summary(.x)$r.squared),
    etiqueta = case_when(slope > 0 & r2 > 0.75 ~ "R",
                         slope < 0 & r2 > 0.75 ~ "F",
                         TRUE ~ "S"),
    FechaEtiqueta = FechaInicio + 15
  )

fechas_verticales <- as.Date(c("2020-07-01","2020-11-04","2021-03-24","2021-06-30","2021-12-15","2022-03-09"))
factor_escala <- max(datos$Incidencia, na.rm = TRUE) / max(datos$Mortalidad, na.rm = TRUE)
xlim_compartido <- range(datos$Fecha)
max_inc <- max(datos$Incidencia, na.rm = TRUE)

y_variant  <- 1.28
y_letters  <- 1.13
y_dates    <- 0.99
y_top_lim  <- 1.34

variant_periods <- tibble(
  variant = c("Early wave", "Wuhan", "D614G", "Beta", "Delta", "Omicron BA.1", "Omicron BA.2"),
  inicio  = c(min(datos$Fecha), fechas_verticales),
  fin     = c(fechas_verticales, max(datos$Fecha))
) %>%
  mutate(mid = inicio + floor(as.numeric(fin - inicio) / 2))

# ============================================================
# BASE DE DATOS DE MEDIDAS
# ============================================================
medidas <- tribble(
  ~categoria,           ~subfase,                   ~inicio,               ~fin,                  ~nivel, ~end_type,
  
  "Border closure",     "Initial closure",           as.Date("2020-03-20"), as.Date("2020-04-01"), 1,      "closed",
  "Border closure",     "Full shutdown",              as.Date("2020-04-01"), as.Date("2020-10-31"), 1,      "closed",
  "Border closure",     "Restricted reopening",       as.Date("2021-01-01"), as.Date("2021-11-14"), 2,      "closed",
  "Border closure",     "Expanded reopening",         as.Date("2021-11-15"), as.Date("2022-05-31"), 3,      "closed",
  
  "Localised lockdown", "First lockdown",             as.Date("2020-03-24"), as.Date("2020-04-30"), 1,      "closed",
  "Localised lockdown", "curfew",            as.Date("2020-09-01"), as.Date("2020-09-15"), 1,      "closed",
  "Localised lockdown", "prolonged restrictions",             as.Date("2021-01-01"), as.Date("2021-07-31"), 1,      "closed",
  
  "Isolation",          "Institutional",              as.Date("2020-03-11"), as.Date("2020-10-13"), 1,      "closed",
  "Isolation",          "Home + PHC follow-up",       as.Date("2020-10-13"), as.Date("2021-03-02"), 2,      "closed",
  "Isolation",          "Home, reinforced",           as.Date("2021-03-02"), as.Date("2021-07-17"), 1.5,    "closed",
  "Isolation",          "Home as default",            as.Date("2021-07-17"), as.Date("2022-01-13"), 2,      "closed",
  "Isolation",          "relaxed",                    as.Date("2022-01-13"), as.Date("2022-05-31"), 3,      "closed",
  
  "Vaccination",        "Health workers",             as.Date("2021-03-01"), as.Date("2021-04-30"), 1,      "closed",
  "Vaccination",        "Vulnerable groups",          as.Date("2021-04-01"), as.Date("2021-05-31"), 1,      "closed",
  "Vaccination",        "General adults",             as.Date("2021-06-01"), as.Date("2021-08-31"), 1,      "closed",
  "Vaccination",        "Children",                   as.Date("2021-09-01"), as.Date("2021-12-31"), 1,      "closed",
  "Vaccination",        "Boosters",                   as.Date("2021-11-01"), as.Date("2022-05-31"), 1,      "closed",
  
  "Mask mandate",       "Mandatory mask use",         as.Date("2020-03-24"), as.Date("2022-05-30"), 1,      "closed"
) %>%
  mutate(
    categoria = factor(categoria, levels = c("Mask mandate","Vaccination","Isolation",
                                             "Localised lockdown","Border closure")),
    ypos = as.numeric(categoria),
    punto_medio = inicio + floor(as.numeric(fin - inicio) / 2)
  )

colores_medidas <- c(
  "Border closure"     = "#B71C1C",
  "Localised lockdown" = "#EF6C00",
  "Isolation"          = "#1A239E",
  "Vaccination"        = "#1B5E20",
  "Mask mandate"       = "#616161"
)

medidas <- medidas %>%
  mutate(
    nivel_txt = case_when(
      nivel <= 1   ~ "Strict",
      nivel <= 2   ~ "Intermediate",
      TRUE         ~ "Relaxed"
    ),
    nivel_txt = factor(nivel_txt, levels = c("Strict", "Intermediate", "Relaxed")),
    alpha_val = case_when(
      nivel_txt == "Strict"       ~ 1.00,
      nivel_txt == "Intermediate" ~ 0.62,
      TRUE                        ~ 0.30
    ),
    linetype_borde = if_else(end_type == "closed", "solid", "22")
  )

# ============================================================
# PANEL PRINCIPAL
# ============================================================

# --- FIX: cuadro de leyenda de curvas más ancho, para que el texto no se corte ---
legend_xmin <- as.Date("2021-10-25")   # antes: 2022-01-01 (se recorre a la izquierda para dar más ancho)
legend_xmax <- as.Date("2022-05-31")
legend_text_x <- legend_xmin + 45      # inicio del texto, con margen desde el borde del cuadro
legend_seg_x0 <- legend_xmin + 12
legend_seg_x1 <- legend_xmin + 40

p_principal <- ggplot(datos, aes(x = Fecha)) +
  geom_col(aes(y = Incidencia), fill = "#D1D3D4", width = 1, alpha = 0.8) +
  geom_line(aes(y = MA14_Incidencia), color = "#1A239E", linewidth = 1.2, na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, aes(y = Incidencia, group = Periodo),
              color = "#B71C1C", linewidth = 0.8) +
  geom_line(aes(y = MA14_Mortalidad * factor_escala * 0.3),
            color = "#1B5E20", linewidth = 1.2, na.rm = TRUE, linetype = "solid") +
  geom_vline(xintercept = resumen_incidencia$FechaInicio,
             linetype = "dashed", color = "#616161", alpha = 0.5, linewidth = 0.4) +
  geom_vline(xintercept = fechas_verticales, linetype = "solid", color = "#212121") +
  
  geom_text(data = variant_periods,
            aes(x = mid, y = max_inc * y_variant, label = variant),
            color = "#212121", size = 3.1, fontface = "bold", vjust = 0) +
  
  geom_label(data = resumen_incidencia,
             aes(x = FechaEtiqueta, y = max_inc * y_letters, label = etiqueta),
             color = "#B71C1C", size = 2.6, fontface = "bold",
             label.size = 0, label.padding = unit(0.08, "lines"),
             fill = alpha("white", 0.75)) +
  
  geom_text(data = tibble(fecha = fechas_verticales),
            aes(x = fecha + 4, y = max_inc * y_dates, label = format(fecha, "%Y-%m-%d")),
            color = "#212121", size = 2.2, angle = 90, hjust = 0, vjust = 0.5,
            fontface = "italic", inherit.aes = FALSE) +
  
  scale_y_continuous(
    name = "Incidence rate (cases per 100,000 population)",
    sec.axis = sec_axis(trans = ~./(factor_escala * 0.3),
                        name = "Mortality rate (deaths per 100,000 population)",
                        breaks = scales::pretty_breaks(n = 6)),
    limits = c(0, max_inc * y_top_lim),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_date(limits = xlim_compartido, date_breaks = "3 months",
               date_labels = "%b %Y", expand = expansion(mult = c(0.02, 0.02))) +
  labs(title = "COVID-19 incidence and mortality rates in Havana, Cuba, 2020-2022",
       subtitle = "14-day moving averages (MA) with 14-day linear trend segments; labels indicate predominant SARS-CoV-2 variant/lineage per period",
       x = NULL) +
  theme_minimal(base_size = 13, base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 15, margin = margin(b = 6)),
    plot.subtitle = element_text(hjust = 0, color = "#424242", size = 10.5, margin = margin(b = 14)),
    # --- FIX: título del eje Y izquierdo pegado a los valores ---
    axis.title.y.left = element_text(color = "#1A239E", face = "bold", size = 10.5, margin = margin(r = 1)),
    axis.title.y.right = element_text(color = "#1B5E20", face = "bold", size = 10.5, margin = margin(l = 1)),
    axis.text.y.left = element_text(color = "#1A239E", size = 10),
    axis.text.y.right = element_text(color = "#1B5E20", size = 10),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    plot.margin = margin(25, 25, 5, 8),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  # --- Cuadro de leyenda de curvas, ahora más ancho ---
  annotate("rect",
           xmin = legend_xmin, xmax = legend_xmax,
           ymin = max_inc * 0.60, ymax = max_inc * 0.82,
           fill = "white", alpha = 0.85, color = "#BDBDBD", linewidth = 0.3) +
  annotate("segment", x = legend_seg_x0, xend = legend_seg_x1,
           y = max_inc * 0.79, yend = max_inc * 0.79,
           color = "#1A239E", linewidth = 1.2) +
  annotate("text", x = legend_text_x, y = max_inc * 0.79,
           label = "14-day MA incidence", hjust = 0, size = 3.4, color = "#212121") +
  annotate("segment", x = legend_seg_x0, xend = legend_seg_x1,
           y = max_inc * 0.72, yend = max_inc * 0.72,
           color = "#1B5E20", linewidth = 1.2) +
  annotate("text", x = legend_text_x, y = max_inc * 0.72,
           label = "14-day MA mortality", hjust = 0, size = 3.4, color = "#212121") +
  annotate("segment", x = legend_seg_x0, xend = legend_seg_x1,
           y = max_inc * 0.65, yend = max_inc * 0.65,
           color = "#B71C1C", linewidth = 0.8) +
  annotate("text", x = legend_text_x, y = max_inc * 0.65,
           label = "linear trends", hjust = 0, size = 3.4, color = "#212121")

# ============================================================
# PANEL DE MEDIDAS
# ============================================================
p_medidas <- ggplot(medidas) +
  geom_vline(xintercept = fechas_verticales, linetype = "solid", color = "#212121", alpha = 0.3) +
  geom_rect(
    aes(xmin = inicio, xmax = fin, ymin = ypos - 0.30, ymax = ypos + 0.30,
        fill = categoria, alpha = alpha_val, linetype = linetype_borde),
    color = "black", linewidth = 0.3
  ) +
  geom_text_repel(
    aes(x = punto_medio, y = ypos, label = subfase),
    size = 2.4, color = "black", fontface = "bold", lineheight = 0.85,
    nudge_y = 0.50,
    direction = "x",
    segment.size = 0.3, segment.color = "grey40",
    min.segment.length = 0,
    box.padding = 0.1,
    point.padding = 0.05,
    force = 2,
    max.overlaps = Inf,
    seed = 42
  ) +
  scale_fill_manual(values = colores_medidas, guide = "none") +
  scale_alpha_continuous(range = c(0.30, 1.00), guide = "none") +
  scale_linetype_identity(guide = "none") +
  scale_x_date(limits = xlim_compartido, date_breaks = "3 months",
               date_labels = "%b %Y", expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(
    breaks = 1:5, labels = levels(medidas$categoria),
    limits = c(0.4, 5.9),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Date", y = NULL) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(size = 10.5, color = "#212121", face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 6)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#F0F0F0", linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(5, 25, 10, 8),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# ============================================================
# COMBINAR PANELES
# ============================================================
figura_final <- p_principal / p_medidas +
  plot_layout(heights = c(2.7, 1.7)) +
  plot_annotation(
    caption = paste(
      "Data source: Ministry of Public Health (MINSAP) COVID-19 registries. ",
      "R: Rising trend (slope > 0, R2 > 0.75) | F: Falling | S: Stationary/unclear.\n",
      "Vertical lines mark transitions between periods of predominant circulating SARS-CoV-2 variant/lineage (dates shown above the curve)."
    )
  ) &
  theme(plot.caption = element_text(hjust = 0, color = "#757575", size = 8, lineheight = 1.15))

figura_final

dir.create(fig_dir <- if (exists("DIR")) DIR$figures else file.path("output", "figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(fig_dir, "fig1_epidemic_curve_measures.pdf"), figura_final, width = 12.5, height = 9.8, units = "in")
ggsave(file.path(fig_dir, "fig1_epidemic_curve_measures.png"), figura_final, width = 12.5, height = 9.8, units = "in", dpi = 600)
