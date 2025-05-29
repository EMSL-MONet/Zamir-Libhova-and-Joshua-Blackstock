# Load required libraries
library(tidyverse)
library(vegan)
library(ggrepel)
library(patchwork)
library(readr)

# ==== Read data files =====
# From processed_data folder
fticr_data <- read_csv("processed_data/Processed_FTICR_60880_Data.csv")
names(fticr_data)[1] = 'molecular_formula'
mol_properties <- read_csv("processed_data/Processed_FTICR_60880_Mol.csv")
mol_properties = mol_properties[,-1]
vk_abundance <- read_csv("processed_data/Relative_Abundance_FTICR_60880_VK_Class1.csv")
summary_properties <- read_csv("processed_data/Summary_FTICR_60880_Properties.csv")

# From data_from_EMSL_data_Central folder
coordinates <- read_csv("data_from_EMSL_data_Central/Coordinates.csv")%>%
  mutate(
    proposal_id = as.character(proposal_id),
    sampling_set = as.character(sampling_set)
  )

# ===== Data processing for visualization =====
sample_info <- tibble(
  sample_name = colnames(fticr_data)[colnames(fticr_data) != "molecular_formula"]
) %>%
  separate(sample_name, into = c("proposal_id", "sample_set", "depth"), 
           sep = "_", remove = FALSE) %>%
  mutate(
    proposal_id = as.character(proposal_id),
    sample_set = as.character(sample_set)
  ) %>%
  left_join(coordinates, by = c("proposal_id", "sample_set" = "sampling_set"))


# Prepare data for NMDS analysis
fticr_matrix <- fticr_data %>%
  column_to_rownames("molecular_formula") %>%
  t()

# Create metadata for plots - directly from column names
metadata <- tibble(
  sample_name = colnames(fticr_data)[colnames(fticr_data) != "molecular_formula"]
) %>%
  separate(sample_name, into = c("proposal_id", "sample_set", "depth"), 
           sep = "_", remove = FALSE) %>%
  mutate(
    proposal_id = as.character(proposal_id),
    sample_set = as.character(sample_set)
  ) %>%
  left_join(coordinates, by = c("proposal_id", "sample_set" = "sampling_set"))


# Join molecular formula data with properties for VK plots
formula_presence <- fticr_data %>%
  pivot_longer(
    cols = -molecular_formula,
    names_to = "sample_name",
    values_to = "present"
  ) %>%
  filter(present == 1) %>%
  left_join(mol_properties, by = "molecular_formula")

# Color palette for van Krevelen classes
vk_colors <- c(
  "Amino Sugar" = "#FF9AA2",
  "Carbohydrate" = "#FFDAC1",
  "Cond Hydrocarbon" = "#E2F0CB",
  "Lignin" = "#B5EAD7",
  "Lipid" = "#C7CEEA",
  "Other" = "#9AA6C4",
  "Protein" = "#F5B5FC",
  "Tannin" = "#BDB2FF",
  "Unsat Hydrocarbon" = "#FFCAB0"
)

# ==== 1. NMDS Plot ====

set.seed(123) # For reproducibility
nmds_result <- metaMDS(fticr_matrix, distance = "jaccard", trymax = 100)

# Extract NMDS coordinates and join with metadata
nmds_coords <- as.data.frame(scores(nmds_result, display = "sites")) %>%
  rownames_to_column("sample_name") %>%
  left_join(metadata, by = "sample_name")

# Create NMDS plot
nmds_plot <- ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = geo_loc_name, shape = depth)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  labs(title = "NMDS Analysis of FTICR Data",
       subtitle = paste("Stress =", round(nmds_result$stress, 3)),
       color = "Location",
       shape = "Depth") +
  theme(legend.position = "right") +
  coord_fixed()

print(nmds_plot)
ggsave("plots/nmds_plot.png", nmds_plot, width = 8, height = 6, dpi = 300)

# ==== 2. Van Krevelen Plots ====

# Function to draw VK class regions
add_vk_regions <- function() {
  # Define class boundaries based on H:C and O:C ratios
  regions <- list(
    "Lipid" = data.frame(x = c(0, 0.3, 0.3, 0), y = c(1.5, 1.5, 2.5, 2.5)),
    "Protein" = data.frame(x = c(0.3, 0.6, 0.6, 0.3), y = c(1.5, 1.5, 2.2, 2.2)),
    "Amino Sugar" = data.frame(x = c(0.55, 0.7, 0.7, 0.55), y = c(1.5, 1.5, 2.2, 2.2)),
    "Carbohydrate" = data.frame(x = c(0.6, 1.2, 1.2, 0.6), y = c(1.5, 1.5, 2.2, 2.2)),
    "Lignin" = data.frame(x = c(0.1, 0.6, 0.6, 0.1), y = c(0.7, 0.7, 1.5, 1.5)),
    "Tannin" = data.frame(x = c(0.6, 1.2, 1.2, 0.6), y = c(0.5, 0.5, 1.5, 1.5)),
    "Unsat Hydrocarbon" = data.frame(x = c(0, 0.1, 0.1, 0), y = c(0.7, 0.7, 2.0, 2.0)),
    "Cond Hydrocarbon" = data.frame(x = c(0, 0.7, 0.7, 0), y = c(0.2, 0.2, 0.7, 0.7))
  )
  
  region_polygons <- lapply(names(regions), function(name) {
    geom_polygon(data = regions[[name]], aes(x = x, y = y), 
                 fill = NA, color = "darkgrey", linetype = "dashed", inherit.aes = FALSE)
  })
  
  labels <- lapply(names(regions), function(name) {
    mid_x <- mean(regions[[name]]$x)
    mid_y <- mean(regions[[name]]$y)
    geom_text(aes(x = mid_x, y = mid_y, label = name), 
              inherit.aes = FALSE, size = 3, fontface = "italic", color = "grey30")
  })
  
  c(region_polygons, labels)
}

# Create Van Krevelen plots for each sample
vk_plots <- formula_presence %>%
  group_by(sample_name) %>%
  group_map(~ {
    sample <- .y$sample_name
    ggplot(.x, aes(x = OC, y = HC, color = Class1)) +
      geom_point(alpha = 0.7, size = 1) +
      scale_color_manual(values = vk_colors) +
      add_vk_regions() +
      xlim(0, 1.2) +
      ylim(0, 2.5) +
      labs(title = paste("Van Krevelen Plot -", sample),
           x = "O:C Ratio",
           y = "H:C Ratio",
           color = "Compound Class") +
      theme_minimal() +
      theme(legend.position = "bottom",
            legend.box = "horizontal",
            legend.key.size = unit(0.5, "cm"))
  })

# Save all VK plots
for(i in seq_along(vk_plots)) {
  sample_name <- formula_presence %>% 
    distinct(sample_name) %>% 
    pull(sample_name) %>% 
    .[i]
  
  ggsave(
    filename = paste0("plots/vk_plot_", sample_name, ".png"),
    plot = vk_plots[[i]],
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ==== 3. Bar plot of relative abundance for TOP samples ====
# Filter for TOP samples and prepare data
top_vk_abundance <- vk_abundance %>%
  filter(str_detect(sample_name, "_TOP$")) %>%
  pivot_longer(
    cols = -sample_name,
    names_to = "class",
    values_to = "abundance"
  ) %>%
  mutate(
    class = factor(class),
    sample_short = str_extract(sample_name, "\\d+_\\d+") # Extract sample identifiers
  )

# Create bar plot
bar_plot <- ggplot(top_vk_abundance, 
                   aes(x = sample_short, y = abundance, fill = class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = vk_colors) +
  labs(title = "Relative Abundance of Van Krevelen Classes in TOP Samples",
       x = " ",
       y = "Relative Abundance (%)",
       fill = "Compound Class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

print(bar_plot)
ggsave("plots/vk_abundance_top_bar_plot.png", bar_plot, width = 10, height = 6, dpi = 300)

# Filter for TOP samples and prepare data
btm_vk_abundance <- vk_abundance %>%
  filter(str_detect(sample_name, "_BTM$")) %>%
  pivot_longer(
    cols = -sample_name,
    names_to = "class",
    values_to = "abundance"
  ) %>%
  mutate(
    class = factor(class),
    sample_short = str_extract(sample_name, "\\d+_\\d+") # Extract sample identifiers
  )

# Create bar plot
bar_plot <- ggplot(btm_vk_abundance, 
                   aes(x = sample_short, y = abundance, fill = class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = vk_colors) +
  labs(title = "Relative Abundance of Van Krevelen Classes in BTM Samples",
       x = " ",
       y = "Relative Abundance (%)",
       fill = "Compound Class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

print(bar_plot)
ggsave("plots/vk_abundance_btm_bar_plot.png", bar_plot, width = 10, height = 6, dpi = 300)

# ===== 4. Visualize AImod, NOSC, and GFE across samples =====
# Join sample info with summary properties
sample_properties <- summary_properties %>%
  mutate(
    sample_name = paste(Proposal_ID, Sample_number, Depth, sep = "_"),
    Proposal_ID = as.character(Proposal_ID),
    Sample_number = as.character(Sample_number)
  ) %>%
  left_join(coordinates, by = c("Proposal_ID" = "proposal_id", 
                                "Sample_number" = "sampling_set"))

# Prepare data for visualization
properties_long <- sample_properties %>%
  select(sample_name, Mean_AImod, Mean_NOSC, Mean_GFE, geo_loc_name, Depth) %>%
  pivot_longer(
    cols = c(Mean_AImod, Mean_NOSC, Mean_GFE),
    names_to = "property",
    values_to = "value"
  ) %>%
  mutate(
    property = factor(property, 
                      levels = c("Mean_AImod", "Mean_NOSC", "Mean_GFE"),
                      labels = c("Mean Modified Aromaticity Index", "Mean NOSC", "Mean Gibbs Free Energy"))
  )

# Create a faceted plot for molecular properties
properties_plot <- ggplot(properties_long, 
                          aes(x = sample_name, y = value, fill = geo_loc_name, shape = Depth)) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~property, scales = "free_y", ncol = 1) +
  theme_bw() +
  labs( x = "",
       y = "Value",
       fill = "Location",
       shape = "Depth") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

print(properties_plot)
ggsave("plots/Mean_molecular_properties_plot.png", properties_plot, width = 10, height = 8, dpi = 300)

#Boxplot visualization for molecular properties by location and depth
properties_boxplot <- ggplot(properties_long, 
                             aes(x = geo_loc_name, y = value, fill = Depth)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~property, scales = "free_y") +
  theme_bw() +
  labs( x = "Location",
       y = "Value",
       fill = "Depth") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(properties_boxplot)
ggsave("plots/Mean_molecular_properties_boxplot.png", properties_boxplot, width = 10, height = 8, dpi = 300)

# === Unique and shared molecular formulas ===
# Analyze unique and shared formulas between TOP and BTM samples
formula_depth_analysis <- fticr_data %>%
  pivot_longer(
    cols = -molecular_formula,
    names_to = "sample_name",
    values_to = "present"
  ) %>%
  filter(present == 1) %>%
  mutate(depth = str_extract(sample_name, "(TOP|BTM)$")) %>%
  group_by(molecular_formula) %>%
  summarize(
    n_samples = n(),
    in_top = any(depth == "TOP"),
    in_btm = any(depth == "BTM"),
    category = case_when(
      in_top & in_btm ~ "Shared",
      in_top ~ "TOP only",
      in_btm ~ "BTM only"
    )
  )

# Join with molecular properties
formula_depth_properties <- formula_depth_analysis %>%
  left_join(mol_properties, by = "molecular_formula")

# Create summary of unique/shared formulas
formula_summary <- formula_depth_analysis %>%
  count(category) %>%
  mutate(percentage = n / sum(n) * 100)

# Plot Venn diagram-like visualization
venn_data <- formula_summary %>%
  mutate(label = paste0(category, "\n", n, " (", round(percentage, 1), "%)"))

ggplot(venn_data, aes(x = "", y = n, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Set2") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  theme_void() +
  labs(title = "Distribution of Molecular Formulas by Depth",
       fill = "Category") +
  theme(legend.position = "bottom")

ggsave("plots/formula_distribution_by_depth.png", width = 8, height = 8, dpi = 300)

# Analyze molecular properties by category
category_properties <- formula_depth_properties %>%
  group_by(category) %>%
  summarize(
    mean_AImod = mean(AImod, na.rm = TRUE),
    mean_NOSC = mean(NOSC, na.rm = TRUE),
    mean_GFE = mean(GFE, na.rm = TRUE),
    mean_OC = mean(OC, na.rm = TRUE),
    mean_HC = mean(HC, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Visualize properties by category
category_plot <- formula_depth_properties %>%
  pivot_longer(cols = c(AImod, NOSC, GFE, OC, HC),
               names_to = "property",
               values_to = "value") %>%
  ggplot(aes(x = category, y = value, fill = category)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~property, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(title = "Molecular Properties by Depth Category",
       x = "",
       y = "Value",
       fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/molecular_properties_by_category.png", category_plot, 
       width = 10, height = 8, dpi = 300)

# === VK ====
# Compare van Krevelen classes between TOP and BTM
vk_class_depth <- vk_abundance %>%
  mutate(depth = str_extract(sample_name, "(TOP|BTM)$")) %>%
  pivot_longer(
    cols = -c(sample_name, depth),
    names_to = "class",
    values_to = "abundance"
  )

# Statistical test for each class
vk_stats <- vk_class_depth %>%
  group_by(class) %>%
  do(test_result = wilcox.test(abundance ~ depth, data = ., paired = TRUE)) %>%
  mutate(p_value = test_result$p.value,
         significance = if_else(p_value < 0.05, "*", ""))

# Calculate mean differences
vk_diff <- vk_class_depth %>%
  group_by(class, depth) %>%
  summarize(mean_abundance = mean(abundance, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = depth, values_from = mean_abundance) %>%
  mutate(difference = TOP - BTM) %>%
  left_join(select(vk_stats, class, p_value, significance), by = "class")

# Plot differences with significance
diff_plot <- ggplot(vk_diff, aes(x = reorder(class, difference), y = difference, fill = difference > 0)) +
  geom_col() +
  geom_text(aes(label = significance, y = ifelse(difference > 0, difference + 1, difference - 1)), 
            size = 5) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#66c2a5", "FALSE" = "#fc8d62"),
                    labels = c("TRUE" = "Higher in TOP", "FALSE" = "Higher in BTM"),
                    name = "") +
  labs(title = "Difference in Van Krevelen Class Abundance (TOP - BTM)",
       subtitle = "* indicates p < 0.05",
       x = "Compound Class",
       y = "Difference in Relative Abundance (%)") +
  theme_bw()

ggsave("plots/vk_class_difference.png", diff_plot, width = 9, height = 7, dpi = 300)
