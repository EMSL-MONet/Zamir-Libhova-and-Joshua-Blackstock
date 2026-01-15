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

# ---- Read data ----
ft_raw  <- readr::read_csv(fticr_path, show_col_types = FALSE)
mol_raw <- readr::read_csv(mol_path, show_col_types = FALSE)
meta    <- readr::read_csv(meta_path, show_col_types = FALSE)

# ---- Standardize key column names ----
# FTICR data: first column is molecular formula (sometimes "Unnamed: 0")
formula_col <- names(ft_raw)[1]
ft <- ft_raw %>%
  rename(molecular_formula = all_of(formula_col)) %>%
  mutate(molecular_formula = str_squish(molecular_formula))

mol <- mol_raw %>%
  mutate(molecular_formula = str_squish(molecular_formula)) %>%
  select(-any_of("Unnamed: 0"))

# ---- Metadata cleanup ----
meta <- meta %>%
  mutate(
    Proposal_ID     = as.character(Proposal_ID),
    Sampling_Set    = as.character(Sampling_Set),
    Core_ID         = as.character(Core_ID),
    Clay_Mineralogy = as.character(Clay_Mineralogy),
    Site_Name       = as.character(Site_Name)
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

# ---- Robust NA/blank handling + mineralogy normalization (incl Carbonatic) ----
sample_key <- sample_key %>%
  mutate(
    Clay_Mineralogy = replace_na(Clay_Mineralogy, "Unknown"),
    Clay_Mineralogy = str_squish(Clay_Mineralogy),
    Site_Name       = replace_na(Site_Name, "Unknown"),
    
    # treat blanks as NA then fill Core_ID per row
    Core_ID         = na_if(Core_ID, ""),
    Core_ID         = coalesce(Core_ID, paste0(Proposal_ID, "_", Sampling_Set)),
    
    Depth           = factor(Depth, levels = c("TOP", "BTM"), ordered = TRUE),
    
    # normalized mineralogy class
    mineralogy_class = case_when(
      str_detect(str_to_lower(Clay_Mineralogy), "kaolin") ~ "Kaolinite",
      str_detect(str_to_lower(Clay_Mineralogy), "smect")  ~ "Smectite",
      str_detect(str_to_lower(Clay_Mineralogy), "carbonat|calcit|dolomit|carbonate") ~ "Carbonatic",
      str_detect(str_to_lower(Clay_Mineralogy), "mixed")  ~ "Mixed",
      TRUE ~ "Other/Unknown"
    ),
    
    mineralogy_class = factor(
      mineralogy_class,
      levels = c("Kaolinite", "Smectite", "Carbonatic", "Mixed", "Other/Unknown")
    )
  )

# Optional sanity check
print(sample_key %>% count(Clay_Mineralogy, mineralogy_class, sort = TRUE))

# ---- Long format: formula presence per sample ----
# Assumes 0/1 presence/absence; if intensities, present := value > 0
ft_long <- ft %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to = "sample_name",
               values_to = "value") %>%
  mutate(present = as.integer(value > 0)) %>%
  filter(present == 1) %>%
  select(molecular_formula, sample_name, value, present) %>%
  left_join(sample_key, by = "sample_name") %>%
  left_join(mol, by = "molecular_formula")

# ---- Sample x formula matrix for ordination ----
ft_mat <- ft %>%
  column_to_rownames("molecular_formula") %>%
  t() %>%
  as.matrix()

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

p_nmds <- ggplot(nmds_sites, aes(NMDS1, NMDS2, color = mineralogy_class, shape = Depth)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = mineralogy_class), level = 0.68, linewidth = 0.6, show.legend = FALSE) +
  ggrepel::geom_text_repel(aes(label = Site_Name), size = 3, max.overlaps = 30, show.legend = FALSE) +
  coord_fixed() +
  theme_bw() +
  labs(
    title = "NMDS (Jaccard) of MAOM FTICR-MS composition",
    subtitle = glue("Stress = {round(nmds$stress, 3)}"),
    color = "Mineralogy (normalized)",
    shape = "Depth"
  )

ggsave(file.path(out_dir, "01_NMDS_Jaccard_mineralogy_depth.png"),
       p_nmds, width = 9, height = 6, dpi = 300)

# PERMANOVA: mineralogy + depth (paired within Core_ID where possible)
perm <- adonis2(
  ft_mat_pa ~ mineralogy_class * Depth,
  data = sample_key %>% mutate(Depth = as.character(Depth)),
  method = "jaccard",
  strata = sample_key$Core_ID,
  permutations = 999
)

writeLines(capture.output(perm), con = file.path(out_dir, "01_PERMANOVA_results.txt"))

# ============================================================
# 2) Richness and key molecular indices by group
# ============================================================

sample_metrics <- ft_long %>%
  group_by(sample_name, mineralogy_class, Depth, Site_Name, Core_ID) %>%
  summarise(
    richness   = n_distinct(molecular_formula),
    frac_N     = mean(N > 0, na.rm = TRUE),
    frac_S     = mean(S > 0, na.rm = TRUE),
    frac_P     = mean(P > 0, na.rm = TRUE),
    mean_NOSC  = mean(NOSC, na.rm = TRUE),
    mean_AImod = mean(AImod, na.rm = TRUE),
    mean_GFE   = mean(GFE, na.rm = TRUE),
    mean_OC    = mean(OC, na.rm = TRUE),
    mean_HC    = mean(HC, na.rm = TRUE),
    .groups = "drop"
  )

p_rich <- ggplot(sample_metrics, aes(mineralogy_class, richness, fill = Depth)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12), size = 2, alpha = 0.9) +
  theme_bw() +
  labs(title = "Formula richness by mineralogy and depth", x = "", y = "# formulas (present)")

p_nosc <- ggplot(sample_metrics, aes(mineralogy_class, mean_NOSC, fill = Depth)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12), size = 2, alpha = 0.9) +
  theme_bw() +
  labs(title = "Mean NOSC by mineralogy and depth", x = "", y = "Mean NOSC")

p_fracN <- ggplot(sample_metrics, aes(mineralogy_class, frac_N, fill = Depth)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12), size = 2, alpha = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_bw() +
  labs(title = "Fraction N-containing formulas by mineralogy and depth", x = "", y = "Fraction (N > 0)")

ggsave(file.path(out_dir, "02_Richness_NOSC_FracN_boxplots.png"),
       (p_rich / p_nosc / p_fracN), width = 9, height = 10, dpi = 300)

# ============================================================
# 3) Van Krevelen density + sampled scatter (by class)
# ============================================================

# Pick a class column if present
class_col <- c("Class1", "Class2", "Class3") %>% keep(~ .x %in% names(ft_long)) %>% .[1]
if (is.na(class_col)) class_col <- "Class1"  # fallback

vk_base <- ft_long %>%
  filter(!is.na(OC), !is.na(HC)) %>%
  mutate(Depth = factor(Depth, levels = c("TOP", "BTM"), ordered = TRUE))

p_vk <- ggplot(vk_base, aes(x = OC, y = HC)) +
  stat_bin2d(bins = 45) +
  facet_grid(Depth ~ mineralogy_class) +
  theme_bw() +
  labs(
    title = "Van Krevelen density (O:C vs H:C) by mineralogy and depth",
    x = "O:C", y = "H:C", fill = "Count"
  ) +
  coord_cartesian(xlim = c(0, 1.2), ylim = c(0, 2.5))

ggsave(file.path(out_dir, "03_VK_density_mineralogy_depth.png"),
       p_vk, width = 10, height = 6, dpi = 300)

# Sample for readability (safe slice_sample)
set.seed(1)
set.seed(1)

vk_sample <- vk_base %>%
  group_by(sample_name) %>%
  mutate(.rand = runif(n())) %>%
  arrange(.rand, .by_group = TRUE) %>%
  slice_head(n = 800) %>%
  ungroup() %>%
  select(-.rand)


p_vk_class <- ggplot(vk_sample, aes(x = OC, y = HC, color = .data[[class_col]])) +
  geom_point(alpha = 0.35, size = 0.8) +
  facet_grid(Depth ~ mineralogy_class) +
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
# ============================================================

class_comp <- ft_long %>%
  filter(!is.na(.data[[class_col]])) %>%
  count(sample_name, mineralogy_class, Depth, .data[[class_col]], name = "n") %>%
  group_by(sample_name) %>%
  mutate(rel = n / sum(n)) %>%
  ungroup()

p_class <- ggplot(class_comp,
                  aes(x = sample_name, y = rel, fill = .data[[class_col]])) +
  geom_col(width = 0.9) +
  facet_grid(Depth ~ mineralogy_class, scales = "free_x", space = "free_x") +
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
# 5) Formula-level property distributions by mineralogy/depth
# ============================================================

props_long <- ft_long %>%
  select(sample_name, mineralogy_class, Depth, Core_ID, Site_Name,
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

p_props <- ggplot(props_long, aes(x = mineralogy_class, y = value, fill = Depth)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~property, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(title = "Formula-level molecular property distributions by mineralogy and depth",
       x = "", y = "Value") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(out_dir, "05_Property_distributions_boxplots.png"),
       p_props, width = 12, height = 7, dpi = 300)

# ============================================================
# 6) Shared/unique formulas: 2-class AND 3-class options
# ============================================================

# ---- 6A) Two-class shared/unique (Kaolinite vs Smectite) ----
min_presence_ks <- ft_long %>%
  filter(mineralogy_class %in% c("Kaolinite", "Smectite")) %>%
  group_by(molecular_formula) %>%
  summarise(
    in_kaolinite = any(mineralogy_class == "Kaolinite"),
    in_smectite  = any(mineralogy_class == "Smectite"),
    category = case_when(
      in_kaolinite & in_smectite ~ "Shared",
      in_kaolinite ~ "Kaolinite only",
      in_smectite  ~ "Smectite only",
      TRUE ~ "Other"
    ),
    .groups = "drop"
  ) %>%
  left_join(mol, by = "molecular_formula")

cat_summary <- min_presence_ks %>%
  count(category) %>%
  mutate(pct = n / sum(n))

p_shared <- ggplot(cat_summary, aes(x = category, y = n, fill = category)) +
  geom_col() +
  geom_text(aes(label = paste0(n, "\n(", percent(pct, accuracy = 0.1), ")")),
            vjust = -0.2, size = 3.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none") +
  labs(title = "Shared vs mineralogy-unique formulas (Kaolinite vs Smectite)", x = "", y = "# formulas")

ggsave(file.path(out_dir, "06A_Shared_unique_Kaolinite_vs_Smectite.png"),
       p_shared, width = 9, height = 5, dpi = 300)

cat_data <- min_presence_ks %>%
  filter(category %in% c("Shared", "Kaolinite only", "Smectite only")) %>%
  pivot_longer(cols = c(NOSC, AImod, GFE, OC, HC),
               names_to = "property",
               values_to = "value") %>%
  filter(!is.na(value))

if (nrow(cat_data) == 0) {
  message("Skipping 06A_Category_property_boxplots: no rows after filtering.\n",
          "Check metadata mineralogy labels and mol property columns.")
} else {
  p_cat_props <- ggplot(cat_data, aes(x = category, y = value, fill = category)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    facet_wrap(~property, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          legend.position = "none") +
    labs(title = "Molecular properties of shared vs mineralogy-unique formulas (K vs S)",
         x = "", y = "Value")
  
  ggsave(file.path(out_dir, "06A_Category_property_boxplots_K_vs_S.png"),
         p_cat_props, width = 12, height = 7, dpi = 300)
}

# ---- 6B) Three-class overlap (Kaolinite / Smectite / Carbonatic) ----
# This is better shown with an UpSet plot if you enable ComplexUpset.
# If you don't want extra packages, you can skip this section.

# Uncomment if you installed ComplexUpset:
# if (requireNamespace("ComplexUpset", quietly = TRUE)) {
#   library(ComplexUpset)
#
#   min_presence_3 <- ft_long %>%
#     filter(mineralogy_class %in% c("Kaolinite", "Smectite", "Carbonatic")) %>%
#     group_by(molecular_formula) %>%
#     summarise(
#       Kaolinite  = any(mineralogy_class == "Kaolinite"),
#       Smectite   = any(mineralogy_class == "Smectite"),
#       Carbonatic = any(mineralogy_class == "Carbonatic"),
#       .groups = "drop"
#     )
#
#   p_upset <- ComplexUpset::upset(
#     min_presence_3,
#     intersect = c("Kaolinite", "Smectite", "Carbonatic"),
#     name = "Formula overlap",
#     width_ratio = 0.2
#   )
#
#   ggsave(file.path(out_dir, "06B_UpSet_overlap_K_S_C.png"),
#          p_upset, width = 10, height = 6, dpi = 300)
# } else {
#   message("ComplexUpset not installed; skipping 06B UpSet overlap plot.")
# }

# ============================================================
# 7) Targeted “hypothesis check” metrics
# ============================================================

targeted <- sample_metrics %>%
  mutate(
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

p_target <- ggplot(targeted, aes(x = mineralogy_class, y = value, color = Depth)) +
  geom_point(size = 2.5, position = position_jitter(width = 0.12)) +
  geom_boxplot(aes(fill = Depth), alpha = 0.25, outlier.shape = NA) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(title = "Targeted metrics aligned to mineralogy hypotheses", x = "", y = "Value") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(out_dir, "07_Targeted_metrics_facets.png"),
       p_target, width = 12, height = 7, dpi = 300)

message("All plots saved to: ", out_dir)
message("PERMANOVA results saved to: ", file.path(out_dir, "01_PERMANOVA_results.txt"))
