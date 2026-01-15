# ============================================================
# Exploratory plots for mineral-extracted FTICR-MS (MAOM) data
# Focus: mineralogy (Kaolinite vs Smectite), depth (TOP/BTM),
# and hypotheses about NOSC / compound-class differences.
# ============================================================

# ---- Libraries ----
library(tidyverse)
library(vegan)
library(ggrepel)
library(patchwork)
library(scales)
library(glue)
library(dplyr)

# ---- File paths (edit if needed) ----
fticr_path <- "processed_data/Processed_MAOM_60880_Data.csv"
mol_path   <- "processed_data/Processed_MAOM_60880_Mol.csv"
meta_path  <- "processed_data/Metadata.csv"

# ---- Output folder ----
out_dir <- "plots_exploratory"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Helper: safe %||% ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Read data ----
ft_raw  <- readr::read_csv(fticr_path, show_col_types = FALSE)
mol_raw <- readr::read_csv(mol_path, show_col_types = FALSE)
meta    <- readr::read_csv(meta_path, show_col_types = FALSE)

# ---- Standardize key column names ----
# FTICR data: first column is molecular formula (often "Unnamed: 0")
formula_col <- names(ft_raw)[1]
ft <- ft_raw %>%
  rename(molecular_formula = all_of(formula_col)) %>%
  mutate(molecular_formula = str_squish(molecular_formula))

mol <- mol_raw %>%
  # sometimes has an extra unnamed first column; keep molecular_formula + properties
  mutate(molecular_formula = str_squish(molecular_formula)) %>%
  select(-any_of("Unnamed: 0"))

# Metadata cleanup
meta <- meta %>%
  mutate(
    Proposal_ID   = as.character(Proposal_ID),
    Sampling_Set  = as.character(Sampling_Set),
    Core_ID       = as.character(Core_ID),
    Clay_Mineralogy = as.character(Clay_Mineralogy),
    Site_Name     = as.character(Site_Name)
  )

# ---- Parse sample names from columns like 60880_1_TOP ----
sample_cols <- setdiff(names(ft), "molecular_formula")

sample_key <- tibble(sample_name = sample_cols) %>%
  tidyr::extract(
    sample_name,
    into = c("Proposal_ID", "Sampling_Set", "Depth"),
    regex = "^([^_]+)_([^_]+)_([^_]+)$",
    remove = FALSE
  ) %>%
  mutate(
    Proposal_ID  = as.character(Proposal_ID),
    Sampling_Set = as.character(Sampling_Set),
    Depth        = as.character(Depth)
  ) %>%
  left_join(meta, by = c("Proposal_ID", "Sampling_Set"))

# Flag missing metadata (e.g., Sampling_Set 6)
if (any(is.na(sample_key$Clay_Mineralogy))) {
  missing_sets <- sample_key %>%
    filter(is.na(Clay_Mineralogy)) %>%
    distinct(Proposal_ID, Sampling_Set) %>%
    mutate(id = paste0(Proposal_ID, "_", Sampling_Set)) %>%
    pull(id)
  
  message("WARNING: Missing metadata for: ", paste(missing_sets, collapse = ", "),
          ". These will be labeled as Unknown in plots.")
}

sample_key <- sample_key %>%
  mutate(
    Clay_Mineralogy = replace_na(Clay_Mineralogy, "Unknown"),
    Site_Name       = replace_na(Site_Name, "Unknown"),
    # FIX: vectorized NA replacement
    Core_ID         = if_else(is.na(Core_ID),
                              paste0(Proposal_ID, "_", Sampling_Set),
                              Core_ID),
    Depth           = factor(Depth, levels = c("TOP", "BTM"), ordered = TRUE)
  )
# ---- Long format: formula presence per sample ----
# Assumes 0/1 (presence/absence). If your matrix is intensities, this still works,
# but you may want to treat >0 as present and/or weight summaries by intensity.
ft_long <- ft %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to = "sample_name",
               values_to = "value") %>%
  mutate(
    present = as.integer(value > 0)
  ) %>%
  filter(present == 1) %>%
  select(molecular_formula, sample_name, value, present) %>%
  left_join(sample_key, by = "sample_name") %>%
  left_join(mol, by = "molecular_formula")

# ---- Sample x formula matrix for ordination ----
ft_mat <- ft %>%
  column_to_rownames("molecular_formula") %>%
  t()
ft_mat <- as.matrix(ft_mat)
# convert to presence/absence for Jaccard ordination if not already 0/1
ft_mat_pa <- (ft_mat > 0) * 1

# ============================================================
# 1) Composition differences: NMDS + PERMANOVA
# ============================================================

set.seed(123)
nmds <- metaMDS(ft_mat_pa, distance = "jaccard", trymax = 200, autotransform = FALSE)

nmds_sites <- scores(nmds, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  left_join(sample_key, by = "sample_name")

p_nmds <- ggplot(nmds_sites, aes(NMDS1, NMDS2, color = Clay_Mineralogy, shape = Depth)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = Clay_Mineralogy), level = 0.68, linewidth = 0.6, show.legend = FALSE) +
  ggrepel::geom_text_repel(aes(label = Site_Name), size = 3, max.overlaps = 30, show.legend = FALSE) +
  coord_fixed() +
  theme_bw() +
  labs(
    title = "NMDS (Jaccard) of MAOM FTICR-MS composition",
    subtitle = glue("Stress = {round(nmds$stress, 3)}"),
    color = "Clay mineralogy",
    shape = "Depth"
  )

ggsave(file.path(out_dir, "01_NMDS_Jaccard_mineralogy_depth.png"),
       p_nmds, width = 9, height = 6, dpi = 300)

# PERMANOVA: does mineralogy explain composition? (and depth?)
# Use Core_ID as a strata for paired TOP/BTM permutations when possible
perm <- adonis2(ft_mat_pa ~ Clay_Mineralogy * Depth,
                data = sample_key %>% mutate(Depth = as.character(Depth)),
                method = "jaccard",
                strata = sample_key$Core_ID,
                permutations = 999)

writeLines(capture.output(perm), con = file.path(out_dir, "01_PERMANOVA_results.txt"))

# ============================================================
# 2) Richness (formula count) and key molecular indices by group
#    (Targets: kaolinite higher NOSC; smectite more reduced / aliphatic)
# ============================================================

sample_metrics <- ft_long %>%
  group_by(sample_name, Clay_Mineralogy, Depth, Site_Name, Core_ID) %>%
  summarise(
    richness = n_distinct(molecular_formula),
    frac_N   = mean(N > 0, na.rm = TRUE),
    frac_S   = mean(S > 0, na.rm = TRUE),
    frac_P   = mean(P, na.rm = TRUE),
    mean_NOSC = mean(NOSC, na.rm = TRUE),
    mean_AImod = mean(AImod, na.rm = TRUE),
    mean_GFE = mean(GFE, na.rm = TRUE),
    mean_OC = mean(OC, na.rm = TRUE),
    mean_HC = mean(HC, na.rm = TRUE),
    .groups = "drop"
  )

p_rich <- ggplot(sample_metrics, aes(Clay_Mineralogy, richness, fill = Depth)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12), size = 2, alpha = 0.9) +
  theme_bw() +
  labs(title = "Formula richness by mineralogy and depth", x = "", y = "# formulas (present)")

p_nosc <- ggplot(sample_metrics, aes(Clay_Mineralogy, mean_NOSC, fill = Depth)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12), size = 2, alpha = 0.9) +
  theme_bw() +
  labs(title = "Mean NOSC by mineralogy and depth", x = "", y = "Mean NOSC")

p_fracN <- ggplot(sample_metrics, aes(Clay_Mineralogy, frac_N, fill = Depth)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12), size = 2, alpha = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_bw() +
  labs(title = "Fraction N-containing formulas by mineralogy and depth", x = "", y = "Fraction (N > 0)")

(p_rich / p_nosc / p_fracN) +
  plot_annotation(
    caption = "Exploratory summaries from presence/absence formulas; aligns to hypotheses about oxidation state & N-rich compounds."
  )

ggsave(file.path(out_dir, "02_Richness_NOSC_FracN_boxplots.png"),
       (p_rich / p_nosc / p_fracN), width = 9, height = 10, dpi = 300)

# ============================================================
# 3) Van Krevelen space as density (less overplotting than points)
#    (Targets: smectite higher aliphatic / carb-protein-like; kaolinite more oxidized)
# ============================================================

# Ensure Class1 exists; if not, fall back to any available class column
class_col <- c("Class1", "Class2", "Class3") %>% keep(~ .x %in% names(ft_long)) %>% .[1] %||% "Class1"

# VK density faceted by mineralogy x depth
vk_base <- ft_long %>%
  filter(!is.na(OC), !is.na(HC)) %>%
  mutate(Depth = factor(Depth, levels = c("TOP", "BTM"), ordered = TRUE))

p_vk <- ggplot(vk_base, aes(x = OC, y = HC)) +
  stat_bin2d(bins = 45) +
  facet_grid(Depth ~ Clay_Mineralogy) +
  theme_bw() +
  labs(
    title = "Van Krevelen density (O:C vs H:C) by mineralogy and depth",
    x = "O:C", y = "H:C", fill = "Count"
  ) +
  coord_cartesian(xlim = c(0, 1.2), ylim = c(0, 2.5))

ggsave(file.path(out_dir, "03_VK_density_mineralogy_depth.png"),
       p_vk, width = 10, height = 6, dpi = 300)

# Optional: VK colored by class, but sampled down for readability
set.seed(1)


vk_sample <- vk_base %>%
  group_by(sample_name) %>%
  mutate(.rn = row_number()) %>%
  # shuffle within group
  slice_sample(prop = 1) %>%
  # keep first 800 (or fewer if group smaller)
  slice_head(n = 800) %>%
  ungroup() %>%
  select(-.rn)

p_vk_class <- ggplot(vk_sample, aes(x = OC, y = HC, color = .data[[class_col]])) +
  geom_point(alpha = 0.35, size = 0.8) +
  facet_grid(Depth ~ Clay_Mineralogy) +
  theme_bw() +
  labs(
    title = glue("Van Krevelen scatter (sampled) colored by {class_col}"),
    x = "O:C", y = "H:C", color = class_col
  ) +
  coord_cartesian(xlim = c(0, 1.2), ylim = c(0, 2.5))

ggsave(file.path(out_dir, "03_VK_scatter_class_sampled.png"),
       p_vk_class, width = 11, height = 6, dpi = 300)

# ============================================================
# 4) Compound-class composition by mineralogy/depth (stacked bars)
#    (Targets: more carbohydrate/protein-like and N-containing in smectite)
# ============================================================

class_comp <- ft_long %>%
  filter(!is.na(.data[[class_col]])) %>%
  count(sample_name, Clay_Mineralogy, Depth, .data[[class_col]], name = "n") %>%
  group_by(sample_name) %>%
  mutate(rel = n / sum(n)) %>%
  ungroup()

p_class <- ggplot(class_comp,
                  aes(x = sample_name, y = rel, fill = .data[[class_col]])) +
  geom_col(width = 0.9) +
  facet_grid(Depth ~ Clay_Mineralogy, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = glue("Relative compound-class composition ({class_col})"),
    x = "", y = "Relative abundance (by formula count)", fill = class_col
  )

ggsave(file.path(out_dir, "04_Class_composition_stackedbars.png"),
       p_class, width = 12, height = 7, dpi = 300)

# ============================================================
# 5) Formula-level property distributions (NOSC, AImod) by group
#    (Directly targets the kaolinite-vs-smectite oxidation hypothesis)
# ============================================================

props_long <- ft_long %>%
  select(sample_name, Clay_Mineralogy, Depth, Core_ID, Site_Name,
         molecular_formula, NOSC, AImod, GFE, OC, HC, N, S, P) %>%
  pivot_longer(cols = c(NOSC, AImod, GFE, OC, HC),
               names_to = "property",
               values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(
    property = recode(property,
                      NOSC = "NOSC",
                      AImod = "AImod",
                      GFE = "Gibbs Free Energy",
                      OC = "O:C",
                      HC = "H:C")
  )

p_props <- ggplot(props_long, aes(x = Clay_Mineralogy, y = value, fill = Depth)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~property, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(title = "Formula-level molecular property distributions by mineralogy and depth",
       x = "", y = "Value") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(out_dir, "05_Property_distributions_boxplots.png"),
       p_props, width = 12, height = 7, dpi = 300)

# ============================================================
# 6) “Shared vs unique” formulas: mineralogy contrast
#    (What formulas appear only in kaolinite vs only in smectite vs both?)
# ============================================================

# Determine presence by mineralogy (collapsing across samples within mineralogy)
min_presence <- ft_long %>%
  filter(Clay_Mineralogy %in% c("Kaolinite", "Smectite")) %>%
  group_by(molecular_formula) %>%
  summarise(
    in_kaolinite = any(Clay_Mineralogy == "Kaolinite"),
    in_smectite  = any(Clay_Mineralogy == "Smectite"),
    category = case_when(
      in_kaolinite & in_smectite ~ "Shared",
      in_kaolinite ~ "Kaolinite only",
      in_smectite  ~ "Smectite only",
      TRUE ~ "Other"
    ),
    .groups = "drop"
  ) %>%
  left_join(mol, by = "molecular_formula")

cat_summary <- min_presence %>%
  count(category) %>%
  mutate(pct = n / sum(n))

p_shared <- ggplot(cat_summary, aes(x = category, y = n, fill = category)) +
  geom_col() +
  geom_text(aes(label = paste0(n, "\n(", percent(pct, accuracy = 0.1), ")")),
            vjust = -0.2, size = 3.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none") +
  labs(title = "Shared vs mineralogy-unique formulas", x = "", y = "# formulas")

ggsave(file.path(out_dir, "06_Shared_vs_unique_formulas_mineralogy.png"),
       p_shared, width = 8, height = 5, dpi = 300)

# Compare NOSC (and other props) for those categories
p_cat_props <- min_presence %>%
  filter(category %in% c("Shared", "Kaolinite only", "Smectite only")) %>%
  pivot_longer(cols = c(NOSC, AImod, GFE, OC, HC),
               names_to = "property",
               values_to = "value") %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = category, y = value, fill = category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~property, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none") +
  labs(title = "Molecular properties of shared vs mineralogy-unique formulas",
       x = "", y = "Value")

ggsave(file.path(out_dir, "06_Category_property_boxplots.png"),
       p_cat_props, width = 12, height = 7, dpi = 300)

# ============================================================
# 7) Simple “hypothesis check” plots (targeted summaries)
#    - Kaolinite: higher NOSC + higher O:C (more oxidized)
#    - Smectite: higher fraction of N-containing + more aliphatic (higher H:C) etc.
# ============================================================

targeted <- sample_metrics %>%
  mutate(
    # crude "aliphatic tendency" proxy: mean H:C and lower AImod
    aliphatic_proxy = mean_HC - mean_AImod
  ) %>%
  pivot_longer(cols = c(mean_NOSC, mean_OC, mean_HC, mean_AImod, frac_N, aliphatic_proxy),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = recode(metric,
                         mean_NOSC = "Mean NOSC",
                         mean_OC = "Mean O:C",
                         mean_HC = "Mean H:C",
                         mean_AImod = "Mean AImod",
                         frac_N = "Fraction N-containing",
                         aliphatic_proxy = "Aliphatic proxy (H:C - AImod)"))

p_target <- ggplot(targeted, aes(x = Clay_Mineralogy, y = value, color = Depth)) +
  geom_point(size = 2.5, position = position_jitter(width = 0.12)) +
  geom_boxplot(aes(fill = Depth), alpha = 0.25, outlier.shape = NA) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(title = "Targeted metrics aligned to mineralogy hypotheses", x = "", y = "Value")

ggsave(file.path(out_dir, "07_Targeted_metrics_facets.png"),
       p_target, width = 12, height = 7, dpi = 300)

# ---- Done ----
message("All plots saved to: ", out_dir)
message("PERMANOVA results saved to: ", file.path(out_dir, "01_PERMANOVA_results.txt"))
