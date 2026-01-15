rm(list=ls(all=T))
# ===== Load Libraries =====
# Load required libraries
library(readr)
library(here)
library(tidyverse)
library(readxl)
library(tidyverse)
library(vegan)
library(ggrepel)
library(stringr)
# ===== Load MAOM data from Google Drive =====

import_files = function(FILEPATH){
 
  # Filter for only samples that belong to Joshua's proposal or blanks
  
  zip_filePaths <- list.files(path = FILEPATH, pattern = "60880|Blank", full.names = TRUE, recursive = TRUE)

  
  # now, identify all the MAOM csv files that we just extracted 
  csv_filePaths <- zip_filePaths[grep(pattern = '02Jan25|27Dec24|26Feb25', zip_filePaths)]
  csv_filePaths <- csv_filePaths[!grepl(pattern = 'MeOH', csv_filePaths)]
  csv_filePaths <- csv_filePaths[grep(pattern = '.csv', csv_filePaths)]
  
  # read and combine the MAOM csv files
  icr_dat <- do.call(bind_rows, lapply(csv_filePaths, function(path) {
    
    data = read.csv(path) %>% 
      mutate(source = basename(path)) # add file name
    
  }))
  
  # finally, delete the temporary files from the parent directory
  file.remove(csv_filePaths)
  
  # this is our final output
  icr_dat
}


icr_report = import_files("C:/Users/gara009/Downloads/FY23_MAOM_ProcessedData_extracted/FY23_MAOM_ProcessedData/")



# ====== Initial Processing  ============
# ===== Create the Mol File ======
# Step 1 is to split this file into two: (a) `mol` file, which has info for each molecule/peak; and (b) `dat` file, which has data about the samples


### Create molecular metadata file

# This dataframe contains information pertaining to each peak, including the molecular formula and other molecular indices.  
# We first select only the columns with atomic composition.

mol = 
  icr_report %>% 
  dplyr::select(`Molecular.Formula`, C,H,O,N,S) %>% # NOTE: NO P in the formula assignment on 09/29/25
  rename(molecular_formula = `Molecular.Formula`) %>% 
  distinct()

# Now, we want to process the `mol` file and calculate various indices
# 
# (a) indices
# 
# - AImod (aromatic index), calculated based on [Koch and Dittmar (2016)](https://doi.org/10.1002/rcm.7433)
# - NOSC (nominal oxidation state of carbon), calculated based on [Riedel et al. 2012](https://doi.org/10.1021/es203901u)
# - GFE (Gibbs Free Energy of carbon oxidation), calculated from NOSC, as per [LaRowe & Van Cappellen 2011](https://doi.org/10.1016/j.gca.2011.01.020)
# - H/C, or the ratio of hydrogen to carbon in the molecule
# - O/C, or the ratio of oxygen to carbon in the molecule

mol$P = NA

mol = 
  mol %>% 
  mutate(across(c("N","S","P"), ~replace_na(.,0)),
         AImod = round((1 + C - (0.5*O) - S - (0.5 * (N+P+H)))/(C - (0.5*O) - S - N - P), 4),
         AImod = ifelse(is.na(AImod), 0, AImod),
         AImod = ifelse(AImod == "Inf", 0, AImod),
         AImod = ifelse(AImod == "-Inf", 0, AImod),
         NOSC = round(4 - ((4*C + H - 3*N - 2*O + 5*P - 2*S)/C),4),
         GFE = 60.3-(28.5 * NOSC),
         HC = round(H/C, 2),
         OC = round(O/C, 2)
  )

# (b) Elemental class
# 
# We can also group the molecules based on the elemental composition, i.e. "CHO", "CHONS", "CHONP", "CHONPS" classes


mol = 
  mol %>% 
  mutate(
    El = str_remove_all(molecular_formula, "13C"),
    El = str_remove_all(El, "34S"),
    El = str_remove_all(El, "17O"),
    El = str_remove_all(El, "18O"),
    El = str_remove_all(El, "15N"),
    El = str_remove_all(El, "[0-9]"),
    El = str_remove_all(El, " "))


#(c) molecular class

# Next, we assign classes (aromatic, aliphatic, etc.). These are Van Krevelen classes, typically assigned based on the H/C, O/C, and AImod indices. 
# We calculate three sets of classes; users can use any of these, or assign their own classes as appropriate.
# 
# 1. `Class1` from [Kim et al. 2003](https://doi.org/10.1021/ac034415p) uses H/C and O/C to classify molecules into "lipid", "unsaturated hydrocarbon", "protein", "lignin", "carbohydrate", amino sugar", "tannin", and "condensed hydrocarbon"
# 2. `Class2` from [Seidel et al. 2014](https://doi.org/10.1016/j.gca.2014.05.038) uses H/C, O/C, and AImod to classify molecules into "aromatic", "condensed aromatic", "highly unsaturated compounds including polyphenols/lignins", and "aliphatic"
# 3. `Class3` from [Seidel et al. 2017](https://doi.org/10.3389/feart.2017.00031) includes classes "aromatic", "condensed aromatic", "highly unsaturated compounds including polyphenols/lignins", "carbohydrate", "lipid", "aliphatic", and "aliphatic containing N".

mol = 
  mol %>% 
  mutate(
    Class1 = case_when(HC >= 1.55 & HC <= 2.25 & OC >= 0 & OC <= 0.3 ~ "Lipid",
                       HC >= 0.7 & HC <= 1.5 & OC >= 0.05 & OC <= 0.15 ~ "Unsat Hydrocarbon",
                       HC >= 1.45 & HC <= 2 & OC >= 0.3 & OC <= 0.55 ~ "Protein",
                       HC >= 0.81 & HC <= 1.45 & OC >= 0.28 & OC <= 0.65 ~ "Lignin",
                       HC >= 1.48 & HC <= 2.15 & OC >= 0.68 & OC <= 1 ~ "Carbohydrate",
                       HC >= 1.34 & HC <= 1.8 & OC >= 0.54 & OC <= 0.71 ~ "Amino Sugar",
                       HC >= 0.7 & HC <= 1.3 & OC >= 0.65 & OC <= 1.05 ~ "Tannin",
                       HC >= 0.3 & HC <= 0.81 & OC >= 0.12 & OC <= 0.7 ~ "Cond Hydrocarbon",
                       TRUE ~ "Other"),
    Class2 = case_when(AImod > 0.66 ~ "condensed aromatic",
                       AImod <= 0.66 & AImod > 0.50 ~ "aromatic",
                       AImod <= 0.50 & HC < 1.5 ~ "unsaturated/lignin",
                       HC >= 1.5 ~ "aliphatic"),
    Class2 = replace_na(Class2, "other"),
    Class3 = case_when(AImod > 0.66 ~ "condensed aromatic",
                       AImod <= 0.66 & AImod > 0.50 ~ "aromatic",
                       AImod <= 0.50 & HC < 1.5 ~ "unsaturated/lignin",
                       HC >= 2.0 & OC >= 0.9 ~ "carbohydrate",
                       HC >= 2.0 & OC < 0.9 ~ "lipid",
                       HC < 2.0 & HC >= 1.5 & N == 0 ~ "aliphatic",
                       HC < 2.0 & HC >= 1.5 & N > 0 ~ "aliphatic+N")
  )
         

# ==== Create `dat` file ======
icr_report$source = gsub('Blank','Blank_',icr_report$source)
icr_report$source = gsub('Blank__','Blank_',icr_report$source)
icr_report$source = gsub('QC_','',icr_report$source)

dat <-
  icr_report %>%
  dplyr::select(source, Molecular.Formula, Calculated.m.z, contains("Peak.Area")) %>%
  janitor::clean_names() %>%
  separate(source, sep = "_", into = c("Proposal_ID", "Sampling_Set", "Core_Section", "Rep", "extra1", "acquisition_raw"), extra = "merge") %>%
  mutate(
    Rep = parse_number(Rep),
    acquisition = parse_number(acquisition_raw),  # Extract number from "r1" -> 1
    sample_name = paste0(Proposal_ID, "_", Sampling_Set, "_", Core_Section)
  ) %>%
  select(-extra1, -acquisition_raw)  # Remove temporary columns

# ====== Clean samples ======
# ==== Remove peaks present in 50% of blanks =====
number_of_blanks <- length(unique(icr_report$source[grep('Blank', icr_report$source)]))

peaks_to_remove <-
  dat %>%
  filter(Proposal_ID == 'Blank') %>%  # Only blank samples
  mutate(peak_area = case_when(peak_area > 0 ~ 1)) %>%
  group_by(molecular_formula, calculated_m_z) %>%  # Group by unique peak identifiers
  summarise(
    blank_count = sum(peak_area > 0, na.rm = TRUE),  # Count blanks where peak is present
    .groups = 'drop'
  ) %>%
  filter(blank_count >= number_of_blanks/2) %>%  # Keep peaks in ≥50% of blanks
  select(molecular_formula, calculated_m_z)

# Remove these peaks from main dataset
dat_cleaned <- dat %>%
  anti_join(peaks_to_remove, by = c("molecular_formula", "calculated_m_z")) %>%
  filter(Proposal_ID != 'Blank')

# ===== Merge aquisition ======
#### 3 acquisitions

#Each sample contains a signal from 3 different instrument acquisitions, these need to be merged before we address extraction replicates. 

##  Merge 3 acquisitions, keeping peaks present in ≥2 out of 3 acquisitions

dat_acq = 
  dat_cleaned %>%
  group_by(Proposal_ID, Sampling_Set, Core_Section, Rep, sample_name, 
           molecular_formula, calculated_m_z) %>%
  summarise(
    peak_area_total = sum(peak_area, na.rm = TRUE),
    acquisitions_with_peak = sum(peak_area > 0, na.rm = TRUE),  # Count acquisitions with peak
    n_acquisitions = n(),  # Total acquisitions (should be 3)
    .groups = 'drop'
  ) %>%
  # Keep only peaks present in ≥2 out of 3 acquisitions
  filter(acquisitions_with_peak >= 2) %>%
  # Clean up columns
  select(-acquisitions_with_peak) %>%
  rename(peak_area = peak_area_total)


#### 3 replicates
#We have 3 extractions from each core section and site, the data from these extractions needs to be merged such that you can have one sample per site. Merge replicate extractions where you keep a molecular formula if it was present in at least 2 out of the 3 extraction replicates for that site and core location. 

## Now, we only select molecules that were identified in 2/3 of the total reps

dat_reps_keep = dat_acq %>%
  group_by(Proposal_ID, Sampling_Set, Core_Section, sample_name, 
           molecular_formula, calculated_m_z) %>%
  summarise(
    peak_area_total = sum(peak_area, na.rm = TRUE),
    reps_with_peak = sum(peak_area > 0, na.rm = TRUE),  # Count reps with peak
    n_reps = n(),  # Total reps (should be 3)
    .groups = 'drop'
  ) %>%
  # Keep only peaks present in ≥2 out of 3 reps
  filter(reps_with_peak >= 2) %>%
  # Clean up columns
  select(-reps_with_peak) %>%
  rename(peak_area = peak_area_total)



# This is the list of all the peaks "present" (identified) in our samples.

# ====== Combine this file with the `mol` file =====

icr_processed = 
  dat_reps_keep %>% 
  left_join(mol)

names(icr_processed)

wide_data = icr_processed %>%
  dplyr::select(-Proposal_ID,-Sampling_Set,-Core_Section,-calculated_m_z,
                -C,-H,-O,-N,-S,-AImod,-NOSC,-GFE,
                -El,-Class1,-Class2,-Class3,-HC,-OC)

wide_data <- icr_processed %>%
  distinct(molecular_formula, sample_name) %>%  # Ensure unique combinations
   mutate(value = 1) %>%  # Create an auxiliary column for 1s
  pivot_wider(
    names_from = sample_name,
    values_from = value,
    values_fill = list(value = 0)  # Fill with 0s where there's no association
  )

merged_fticr_data = as.data.frame(wide_data)
row.names(merged_fticr_data) = merged_fticr_data$molecular_formula
row.names(mol) = mol$molecular_formula

merged_fticr_data = merged_fticr_data %>%
  dplyr::select(-molecular_formula)


# =====  Additional Clean up for Data analysis ====

# Extra cleaning if needed 
# clean up missing peaks across all samples
#merged_fticr_data = merged_fticr_data[-which(rowSums(merged_fticr_data) == 0),]

# removing singletons (formulas found only in one site)
singletons = apply(merged_fticr_data, 1, function(x) length(which(x > 0))) # identify
merged_fticr_data = merged_fticr_data[-which(singletons == 1),]

# making sure the columns match in the mol file
mol <- mol[rownames(mol) %in% rownames(merged_fticr_data), ]


# ====== Calculate Average Molecular Indices in the merged file ====
# Calculate averages and total formulas
result <- lapply(1:ncol(merged_fticr_data), function(i) {
  # Create a subset of mol where the corresponding merged_fticr_data column is 1
  subset_mol <- mol[merged_fticr_data[,i] == 1, ]
  
  list(
    Mean_C = mean(subset_mol$C, na.rm = T),
    Mean_H = mean(subset_mol$H, na.rm = T),
    Mean_O = mean(subset_mol$O, na.rm = T),
    Mean_N = mean(subset_mol$N, na.rm = T),
    Mean_S = mean(subset_mol$S, na.rm = T),
    Mean_P = mean(subset_mol$P, na.rm = T),
    Mean_AImod = mean(subset_mol$AImod, na.rm = T),
    Mean_NOSC = mean(subset_mol$NOSC, na.rm = T),
    Mean_GFE = mean(subset_mol$GFE, na.rm = T),
    total_formulas = sum(merged_fticr_data[,i] != 0)
  )
})

# Convert the list to a dataframe
mol_properties_average <- do.call(rbind, lapply(result, as.data.frame))
rownames(mol_properties_average) <- colnames(merged_fticr_data)

# Add Sample_ID and separate into Site_Code and Layer
mol_properties_average <- mol_properties_average %>%
  tibble::rownames_to_column("Sample_ID") %>%
  dplyr::mutate(Sample_ID = as.character(Sample_ID)) %>%
  tidyr::separate(Sample_ID, into = c("Proposal_ID","Sample_number", "Depth"), sep = "_", remove = FALSE) 


# ===== Calculate Percent Relative Abundance of VK classes using class 1 ====
# Calculate the percent relative abundance
vk_classes = icr_processed %>%
  # Group by sample and Class1 (van Krevelen class)
  group_by(sample_name, Class1) %>%
  # Count the number of peaks in each class
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  # Calculate the total number of peaks per sample
  group_by(sample_name) %>%
  mutate(total_peaks = sum(count),
         # Calculate the percent relative abundance
         percent_abundance = (count / total_peaks) * 100) %>%
  # Create a wide format table of the percent abundances
  select(sample_name, Class1, percent_abundance) %>%
  pivot_wider(names_from = "Class1", values_from = "percent_abundance", values_fill = 0)


# ==== Export data =====

write.csv(merged_fticr_data,"processed_data/Processed_MAOM_60880_Data.csv")
write.csv(mol,"processed_data/Processed_MAOM_60880_Mol.csv")

write.csv(mol_properties_average,"processed_data/Summary_MAOM_60880_Properties.csv", row.names = F)
write.csv(vk_classes,"processed_data/Relative_Abundance_FTICR_MAOM_VK_Class1.csv", row.names = F)
