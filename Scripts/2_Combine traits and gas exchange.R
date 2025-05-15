# Getting a handle on the datafile

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

#data <- read_excel("Dataset/French_Guiana/20241014_Dataset_Phenoflux.xlsx", sheet = "for_r") # old - this had incorrections in unique_code
data <- read_excel("Dataset/French_Guiana/20250402_Dataset_Phenoflux.xlsx", sheet = "for_r") # this has now been corrected

#ACi curve parameters from Julien's github.
load("Dataset/French_Guiana/Outputs from Julien Lamour/2_Fitted_ACi_data.rdata")

#fix the class of the columns
data <- data %>%
  mutate(across(c(Age_leaf_max, Age_leaf_min,leaf_thickness:LMA_leaf, nb_tablets_LMA:LMA_lamina, `fv_/_fm`:PTLP, `Vcmax(25)`, `Jmax(25)`), 
                ~ as.numeric(as.character(.))))

#inspect columns
unique(data$Phenoflux_ID) #no longer exists
unique(data$sp) #should only be 8 species so fix the 2 typos

data$sp[data$sp == "Lecythis_poiteauii"] <- "Lecythis_poiteaui"
data$sp[data$sp == "Pradosia_cochleria"] <- "Pradosia_cochlearia"
# now it is all good here.

#now make sure all code_unique... columns say the exact same thing and then clean up by having just one of those columns.
code_unique_cols <- grep("^code_unique", colnames(data), value = TRUE) # Select columns starting with "code_unique"

data$code_unique_match <- apply(data[, code_unique_cols], 1, function(row) {  # Check if all values match for each row
  all(row == row[1], na.rm = TRUE)
})

data$code_unique...2[data$code_unique_match=="FALSE"] #these ones may not be the same
data$code_unique...10[data$code_unique_match=="FALSE"] 
data$code_unique...22[data$code_unique_match=="FALSE"] 
#now correct but note comments in the excel file

# extract date from unique code so that I have it in correct format since other column is being annoying
data <- data %>%
  mutate(date = sub(".*_(.*)", "\\1", code_unique...2))

data <- data %>%
  mutate(date = as.POSIXct(date, format = "%Y%m%d"))

#now extract info leaf_code so I have unique trees and cohorts and leaves
data <- data %>%
  separate(leaf_code, into = c("tree", "cohort", "leaf"), sep = "\\.", remove = FALSE)

#rename the unique code to get rid of the number
data <- rename(data,  "code_unique" = "code_unique...2")

#now get the midpoint between the min and max leaf age
data$Age_leaf_mid <- data$Age_leaf_min + ((data$Age_leaf_max - data$Age_leaf_min) / 2)
data$Age_leaf_range <- data$Age_leaf_max - data$Age_leaf_min

hist(data$Age_leaf_range) #note that range should be 14 days (2 weeks) for all right?

# now add the Bilan to the dataset
#first rename the sample id column in Bilan to match our dataset
Bilan <- rename(Bilan, "code_unique" = "SampleID")

#next check what values are unique to each dataframe
setdiff(unique(Bilan$code_unique), unique(data$code_unique)) #unique to Bilan
setdiff(unique(data$code_unique), unique(Bilan$code_unique)) #none unique to data

#fix the unique code attention_surface part of Bilan
Bilan$code_unique[Bilan$code_unique == "10.1.1_20221011_attention_surface"] <- "10.1.1_20221011"

data <- merge(data, Bilan, all = TRUE, by = "code_unique")

#### NOW ADD THE LEAF ELEMENTAL DATA HERE
nutrients <- read_delim("Dataset/French_Guiana/2303_IC_IM_Resultats_corrige.csv", 
                                            delim = ";", escape_double = FALSE, trim_ws = TRUE)

nutrients_clean <- nutrients %>%
  # Replace ',' with '.' and convert to numeric
  mutate(`%N` = as.numeric(str_replace(`%N`, ",", ".")),
         `%S` = as.numeric(str_replace(`%S`, ",", "."))) %>%
  # Fix missing decimal in %C
  mutate(`%C` = as.numeric(if_else(
    `%C` > 100,  # Threshold to catch unscaled values
    str_replace(as.character(`%C`), "^([0-9]{2})([0-9]+)$", "\\1.\\2"),
    as.character(`%C`))))

nutrients_subset <- nutrients_clean %>%
  select(`ID Demandeur corrige`, `%N`, `%C`, `%S`)

#check that the names match in hte datasets before merging
setdiff(unique(nutrients_subset$`ID Demandeur corrige`), unique(data$code_unique)) #unique to nutrients
setdiff(unique(data$code_unique), unique(nutrients_subset$`ID Demandeur corrige`)) #unique to data

nutrients_subset <- rename(nutrients_subset, "code_unique" = "ID Demandeur corrige")

# Then, perform a left join (or inner_join, if you only want matching rows)
data <- merge(data, nutrients_subset, by = "code_unique") #removes the 1 that was unique to data that was anatomical

data <- rename(data, "Nmass" = "%N", "Cmass" = "%C", "Smass" = "%S")

##### REPLICATION CHECK
unique(data$sp)
unique(data$tree)
unique(data$cohort)
unique(data$leaf)

#how many replicate trees per species?
data %>%
  subset(!is.na(`Vcmax(25)`)) %>% # two missing
  #subset(!is.na(Vcmax25)) %>% # two missing
  #subset(!is.na(LMA_leaf)) %>% # one missing
  #subset(!is.na(leaf_thickness)) %>% # none missing
  #subset(!is.na(PTLP)) %>% # none missing
  group_by(sp) %>%
  summarise(unique_trees = n_distinct(tree)) #3 trees x 8 species = 24 trees total


# how many visits per tree?
#how many cohorts per tree or species?
data %>%
  group_by(sp, tree) %>%
  summarise(unique_visits = n_distinct(date)) %>% print(n = 24) #between 1-6 visits per tree

#how many cohorts per tree or species? Trying to see if I can guess the phenological strategy from this data.
data %>%
  group_by(sp, tree, date) %>%
  summarise(unique_cohorts = n_distinct(cohort)) %>% print(n = 104)
# Caryocar glabrum: 1 cohort at a time except on 2022-11-28
# Eschweilera coriaceae: 1-3 cohorts at a time
# Lecythis poiteaui: only ever 1 cohort at a time.
# Pouteria eurenifolia: 1-3 cohorts at a time
# Pradosia cohclearia: Usually only 1 cohort at a time except 2022-02-01
# Recordoxylon speciosum: Usually only 1 cohort at a time except 2022-10-14
# Sextonia rubra: often 2 cohorts at a time.
# Sterculia pruriens: 1-2 Cohorts at a time.

### NOW ADD ACI CURVE PARAMATER FITS FROM PLANTECOPHYS WITHOUT TPU LIMITATION TO COMPARE
#Aci <- read_csv("Outputs/Aci_NOTPU_plantecophys.csv")
##### OUR FINAL ESTIMATES SELECT WHETHER TPU IMPROVES THE FIT OR NOT AND USES THE BEST MODEL
Aci <- read_csv("Outputs/Processed licor data/final_Aci_data.csv")

setdiff(unique(Aci$SampleID), unique(data$code_unique)) #unique to Aci 
setdiff(unique(data$code_unique), unique(Aci$SampleID)) #none unique to data

#fix the unique code attention_surface part of Bilan
Aci$SampleID[Aci$SampleID == "10.1.1_20221011_attention_surface"] <- "10.1.1_20221011"

Aci$code_unique <- Aci$SampleID
Aci$date <- NULL
Aci$obs_date <- NULL

## NOTE SOME I STILL NEED TO WORK OUT WHAT THE ISSUE IS.
data <- merge(data, Aci, by = c("code_unique", "leaf_code", "tree", "cohort", "leaf"), all = TRUE)

#now get the gasx spot measurements extract artound 415ppm
spots <- read_csv("Outputs/Processed licor data/Gasx_spot.csv")

setdiff(unique(spots$SampleID), unique(data$code_unique)) #unique to spots 
setdiff(unique(data$code_unique), unique(spots$SampleID)) #none unique to data

#fix the unique code attention_surface
spots$SampleID[spots$SampleID == "10.1.1_20221011_attention_surface"] <- "10.1.1_20221011"
spots$code_unique <- spots$SampleID
spots$date <- NULL
spots$obs_date <- NULL

## NOTE SOME I STILL NEED TO WORK OUT WHAT THE ISSUE IS.
data <- merge(data, spots, by = c("code_unique", "leaf_code", "tree", "cohort", "leaf"), all = TRUE)


# CALCULATE N ON AN AREA BASIS
data$Narea <- data$Nmass * data$LMA_leaf
data$Carea <- data$Cmass * data$LMA_leaf
data$CNratio <- data$Cmass / data$Nmass
data$Vcmax25mass <- data$Vcmax_Estimate_final / data$LMA_leaf #using MY data

#get month column
data$month <- lubridate::month(data$date, label = TRUE, abbr = TRUE)

#calculate g1 from the spot measurements.
data$cica <- data$CO2_s/data$CO2_r
data$vpd <- plantecophys::RHtoVPD(data$RHs, data$Tleaf.y)
data$g1.spots <- (data$cica * sqrt(1.6))/(1-data$cica) #ASSUMING vpd OF 1.6 WHICH SHOULD BE THE AVERAGE DURING MEASUREMENTS
data$g1.spots <- (data$cica * sqrt(data$vpd))/(1-data$cica) 


################## 
# ARE LEAF AGE AND SAMPLING DATE CORRELATED IN EACH SPECIES?
#################

# Convert sampling_date to numeric
library(lubridate)

data <- data %>%
  mutate(julian_date = yday(date))  # Or as.numeric(sampling_date) for days since 1970-01-01

# Calculate Pearson correlation within each species
cor_by_species <- data %>%
  group_by(sp) %>%
  summarise(correlation = cor(Age_leaf_max, julian_date, use = "complete.obs"))

cor_by_species

ggplot(data, aes(x = julian_date, y =Age_leaf_max)) +
  geom_point(alpha = 0.4, size = 3, aes(col = cohort)) +
  geom_smooth(method = "lm", se = FALSE, col = "blue") +
  facet_wrap(~sp, scales = "free") +
  theme_bw() +
  labs(x = "Julian Sampling Date", y = "Leaf Age")

# Basic linear model with interaction
lm_model <- lm(Age_leaf_max ~ julian_date * sp, data = data)
summary(lm_model)

library(lme4)
mixed_model <- lmer(Age_leaf_max ~ julian_date * sp + (1|tree), data = data)
summary(mixed_model)
confint(mixed_model)

library(lmerTest)

# Refit model
mod <- lmer(Age_leaf_max ~ julian_date * sp + (1 | tree), data = data)

# View summary with p-values
summary(mod)



######

library(dplyr)
library(purrr)
library(tidyr)
library(broom)

# Run cor.test per species
cor_results <- data %>%
  filter(!is.na(Age_leaf_mid), !is.na(julian_date), !is.na(sp)) %>%
  group_by(sp) %>%
  nest() %>%
  mutate(
    test = map(data, ~ cor.test(.x$julian_date, .x$Age_leaf_mid, method = "pearson")),
    tidied = map(test, broom::tidy)
  ) %>%
  unnest(tidied) %>%
  select(sp, estimate, p.value, conf.low, conf.high)

# View significant results
cor_results %>%
  #ilter(p.value < 0.05) %>%
  arrange(p.value)

# Join correlation results back to the main data
plot_data <- data %>%
  filter(!is.na(Age_leaf_mid), !is.na(julian_date), !is.na(sp)) %>%
  left_join(cor_results, by = "sp") %>%
  mutate(significant = ifelse(p.value < 0.05, "p < 0.05", "ns"))

# Plot
ggplot(plot_data, aes(x = julian_date, y = Age_leaf_mid)) +
  geom_point(alpha = 0.6, aes(fill = cohort), shape = 21) +
  geom_smooth(method = "lm", se = TRUE, col = "black", 
              aes(linetype = significant)) +
  facet_wrap(~ sp, ncol = 4, nrow = 2) +
  geom_text(
    aes(x = -Inf, y = Inf,
        label = paste0("r = ", round(estimate, 2), "\np = ", signif(p.value, 2))),
    hjust = -0.5, vjust = 1.1, size = 3) +
  theme_bw(base_size = 15) +
  labs(x = "Julian sampling date", y = "Leaf age mid (days)") +
  theme(strip.text = element_text(size = 9), legend.position = "top")+
  scale_linetype_manual(values = c("p < 0.05" = "solid", "ns" = "dashed"))

ggsave("./Outputs/Plots/age_vs_season_correlation.png",
       type = "Cairo", bg = "white", dpi = 300, width = 10, height = 6)

###### NOW SHOW HOW LEAF AGE RANGE VARIES
ggplot(data, aes(Age_leaf_range))+
  geom_histogram()+
  theme_bw(base_size = 15)+
  labs(x = "Leaf Age Uncertainty (days)")

ggsave("./Outputs/Plots/Leafage_range.png",
       type = "Cairo", bg = "white", dpi = 300, width = 4, height = 4)

#View(data[data$Age_leaf_range > 23,])

##### FLAG SUSPICIOUS DATA POINTS
# Note this is based off our expected trends with leaf age.
# So don't remove outliers without checking raw data to see if we can explain the weirdness.
# Also - in some cases the leaf age estimate may be wrong. So should double check that as well.
# I am immediately suspicious of Age_leaf_range that are too small or large
# If it was a 3 week not 2 week window then should be around 21 days.
#data$thickness_QC <- 0
#data$LMA_leaf_QC <- 0
#data$Nmass_QC <- 0
#data$PTLP_QC <- 0
#data$QC_totals <- 0

#Here if something looks weird I will give it a 1. 
#If a leaf has multiple weird looking data maybe 
#it is a leaf age calc error thing or a real weird leaf?
# 
# ggplot(data = subset(data, !is.na(Age_leaf_mid)), aes(Age_leaf_mid, leaf_thickness))+
#   geom_point(aes(col = Age_leaf_range))+
#   geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
#   scale_color_viridis_c()+
#   facet_wrap(~sp, ncol = 4, nrow = 2)+
#   theme_bw()
# 
# #thickness suspicions are: 
# #highest for caryocar, 
# #lowest for pradosia
# # and 2 highest for sterculia
# 
# ggplot(data, aes(leaf_thickness, LMA_leaf, col=sp))+
#   geom_point()+
#   facet_wrap(~sp, scales="free", ncol = 4)+
#   theme(legend.position = "top")
# 
# # could identify outliers based on trait coordination?
# 
# # Filter for non-missing values
# data_filtered <- data %>%
#   filter(!is.na(leaf_thickness), !is.na(LMA_leaf), !is.na(sp))
# 
# # Fit models per species and identify outliers
# outlier_data <- data_filtered %>%
#   group_by(sp) %>%
#   nest() %>%
#   mutate(
#     model = map(data, ~ lm(LMA_leaf ~ leaf_thickness, data = .x)),
#     augmented = map2(model, data, ~ augment(.x, data = .y))
#   ) %>%
#   unnest(augmented) %>%
#   group_by(sp) %>%
#   mutate(
#     residual_sd = sd(.resid, na.rm = TRUE),
#     is_outlier = abs(.resid) > 3 * residual_sd
#   )
# 
# 
# ggplot(outlier_data, aes(x = leaf_thickness, y = LMA_leaf)) +
#   geom_point(aes(col = is_outlier), alpha = 0.8) +
#   facet_wrap(~ sp, scales = "free", ncol = 4) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   theme_bw() +
#   theme(legend.position = "top") +
#   labs(colour = "Outlier", title = "Outlier detection in LMA vs leaf thickness")
# 
# ggsave("./Outputs/Plots/Outliers_LMA.png",
#        type = "Cairo", bg = "white", dpi = 300, width = 8, height = 6)

library(purrr)
library(broom)
# 
# flag_outliers <- function(data, xvar, yvar, qc_col, threshold = 2.5) {
#   data %>%
#     filter(!is.na(.data[[xvar]]), !is.na(.data[[yvar]]), !is.na(sp)) %>%
#     group_by(sp) %>%
#     nest() %>%
#     mutate(
#       model = purrr::map(data, ~ lm(as.formula(paste(yvar, "~", xvar)), data = .x)),
#       augmented = purrr::map2(model, data, ~ augment(.x, data = .y))
#     ) %>%
#     unnest(augmented) %>%
#     group_by(sp) %>%
#     mutate(
#       residual_sd = sd(.resid, na.rm = TRUE),
#       is_outlier = abs(.resid) > threshold * residual_sd,
#       !!qc_col := ifelse(is_outlier, 1, 0)
#     ) %>%
#     ungroup() %>%
#     select(-any_of(c(".fitted", ".resid", ".hat", ".sigma", ".cooksd", ".std.resid", "residual_sd", "is_outlier", "data", "model", "augmented")))
# }
# 
# # LMA outliers from leaf_thickness
# data_qc1 <- flag_outliers(data, xvar = "leaf_thickness", yvar = "LMA_leaf", qc_col = "LMA_leaf_QC")
# 
# # thickness outliers from LMA
# data_qc2 <- flag_outliers(data_qc1, xvar = "LMA_leaf", yvar = "leaf_thickness", qc_col = "thickness_QC")
# 
# # TLP outliers from %N
# data_qc3 <- flag_outliers(data_qc2, xvar = "Nmass", yvar = "PTLP", qc_col = "PTLP_QC")
# 
# # N% outliers from TLP
# data_qc4 <- flag_outliers(data_qc3, xvar ="PTLP", yvar = "Nmass", qc_col = "Nmass_QC")
# 
# 
# data_qc <- data_qc4 %>%
#   mutate(QC_totals = rowSums(across(c(LMA_leaf_QC, thickness_QC, PTLP_QC, Nmass_QC)), na.rm = TRUE))
# 
# 
# ggplot(data_qc, aes(x = leaf_thickness, y = LMA_leaf)) +
#   geom_point(aes(col = factor(QC_totals)), alpha = 0.8) +
#   facet_wrap(~ sp, ncol = 4) +
#   theme_bw() +
#   theme(legend.position = "top") +
#   scale_color_brewer(palette = "Dark2")+
#   labs(colour = "Outlier", title = "Outlier detection in LMA vs leaf thickness")
# 
# ggsave("./Outputs/Plots/Outliers_LMA_andthickness.png",
#        type = "Cairo", bg = "white", dpi = 300, width = 8, height = 6)
# 
# ggplot(data_qc, aes(x = Nmass, y = PTLP)) +
#   geom_point(aes(col = factor(QC_totals)), alpha = 0.8) +
#   facet_wrap(~ sp, ncol = 4) +
#   theme_bw() +
#   theme(legend.position = "top") +
#   scale_color_brewer(palette = "Dark2")

#View(data_qc[data_qc$QC_totals > 0,])

### NOW JUST FLAG BASED ON DISTRIBUTION NOT ON TRAIT CORRELATIONS SINCE THERE ARE ASSUMPTIONS INVOLVED THERE
flag_trait_outliers <- function(data, trait_col, qc_col, threshold = 3) {
  data %>%
    group_by(sp) %>%
    mutate(
      trait_mean = mean(.data[[trait_col]], na.rm = TRUE),
      trait_sd = sd(.data[[trait_col]], na.rm = TRUE),
      !!qc_col := case_when(
        is.na(.data[[trait_col]]) ~ NA_real_,  # If the trait is NA, set QC to NA
        abs(.data[[trait_col]] - trait_mean) > threshold * trait_sd ~ 1,  # Flag as outlier if beyond threshold
        TRUE ~ 0  # Otherwise, set QC to 0
      )
    ) %>%
    ungroup() %>%
    select(-trait_mean, -trait_sd)
}

# flagging outliers for each trait
data_qc <- flag_trait_outliers(data, trait_col = "LMA_leaf", qc_col = "LMA_leaf_QC_dist")
data_qc <- flag_trait_outliers(data_qc, trait_col = "leaf_thickness", qc_col = "thickness_QC_dist")
data_qc <- flag_trait_outliers(data_qc, trait_col = "leaf_area", qc_col = "leaf_area_QC_dist")
data_qc <- flag_trait_outliers(data_qc, trait_col = "leaf_drymass", qc_col = "leaf_drymass_QC_dist")
data_qc <- flag_trait_outliers(data_qc, trait_col = "Narea", qc_col = "Narea_QC_dist")
data_qc <- flag_trait_outliers(data_qc, trait_col = "Nmass", qc_col = "Nmass_QC_dist")
data_qc <- flag_trait_outliers(data_qc, trait_col = "Cmass", qc_col = "Cmass_QC_dist")
data_qc <- flag_trait_outliers(data_qc, trait_col = "PTLP", qc_col = "PTLP_QC_dist")
data_qc <- flag_trait_outliers(data_qc, trait_col = "Vcmax_Estimate_final", qc_col = "Vcmax_Estimate_final_QC_dist")
data_qc <- flag_trait_outliers(data_qc, trait_col = "Jmax_Estimate_final", qc_col = "Jmax_Estimate_final_QC_dist")
data_qc <- flag_trait_outliers(data_qc, trait_col = "TPU_Estimate_final", qc_col = "TPU_Estimate_final_QC_dist")

#ADD UP ALL OF THE QC FLAGS
data_qc <- data_qc %>%
  mutate(QC_totals = rowSums(across(c(LMA_leaf_QC_dist:PTLP_QC_dist)), na.rm = TRUE))

####### save this combined data for other stuff
write.csv(data_qc, "./Outputs/Combined_data.csv", row.names = FALSE)




### PRELIMINARY PLOTS
ggplot(data = subset(data, !is.na(Age_leaf_mid)), aes(Age_leaf_mid, leaf_thickness))+
  geom_point(aes(col = Age_leaf_range))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_viridis_c()+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()

ggsave("./Outputs/Plots/frenchguiana_thickness_leafagerange.png", type = "Cairo", bg = "white", dpi = 300, width = 12, height = 5)

ggplot(data, aes(julian_date, leaf_thickness))+
  #geom_point(aes(col = Age_leaf_range))+
  geom_point(aes(col = thickness_QC))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cc"))+
  #scale_color_viridis_c()+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()


#normal functional traits
ggplot(data = subset(data, !is.na(Age_leaf_mid)), aes(Age_leaf_mid, leaf_thickness, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")
#Pradosia cochlearia gets thicker, others don't change

ggsave("./Outputs/Plots/frenchguiana_leafthickness_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)

ggplot(data = subset(data, !is.na(Age_leaf_mid)), aes(Age_leaf_mid, leaf_drymass, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_leafdrymass_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)

ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, LMA_leaf, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")
#LMA somewhat increases for all because stabilising but particularl for Lecythis poiteaui and Pradosia cochlearia (two highest LMA species)

ggsave("./Outputs/Plots/frenchguiana_lma_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)


ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, leaf_area, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_leafarea_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)

#for Pouteria_eugenifolia and Recordoxylon speciosum leaf area declines with age - but that shouldn't be possible.
#Therefore should we account for leaf area when analysing data? see Nara's paper again.

ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, Narea, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_Narea_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)

ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, Nmass, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_Nmass_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)

ggplot(data = subset(data, !is.na(Age_leaf_mid)), aes(Age_leaf_mid, Nmass))+
  geom_point(aes(col=Age_leaf_range))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_viridis_c()


#spot gasx data - not not under steady state
ggplot(data = subset(data, !is.na(Age_leaf_mid)), aes(Age_leaf_mid, A))+
  geom_point(aes(col=Age_leaf_range))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_viridis_c()

ggplot(data = subset(data, !is.na(Age_leaf_mid)), aes(Age_leaf_mid, gsw))+
  geom_point(aes(col=Age_leaf_range))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_viridis_c()

ggplot(data = subset(data, !is.na(Age_leaf_mid)), aes(Age_leaf_mid, A/gsw))+
  geom_point(aes(col=Age_leaf_range))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_viridis_c()

ggplot(data = subset(data, !is.na(Age_leaf_mid)), aes(Age_leaf_mid, g1.spots))+
  geom_point(aes(col=Age_leaf_range))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_viridis_c()


#photosynthetic traits

#compare Vcmax from Julien vs the other data: note "Vcmax25" is Julien Lamour and Vcmax(25) is the original from shiny app
A <- ggplot(data = subset(data, Vcmax_QC != "bad" & Jmax_QC != "bad"), aes(Vcmax25, `Vcmax(25)`, fill = month))+
  geom_point(shape = 21)+
  geom_abline()+
  theme_bw()+
  scale_x_continuous(name = "Vcmax25 (Lamour)")+
  scale_y_continuous(name = "Vcmax25 (Shiny)")
#The two approaches are definitely different between October - December and the rest of the year.
A
B <- ggplot(data = subset(data, Vcmax_QC != "bad" & Jmax_QC != "bad"), aes(Vcmax25, Vcmax_Estimate_final, fill = month))+
  geom_point(shape = 21)+
  geom_abline()+
  theme_bw()+
  scale_x_continuous(name = "Vcmax25 (Lamour)")+
  scale_y_continuous(name = "Vcmax25 (Plantecophys - varying TPU)")

C <- ggplot(data = subset(data, Vcmax_QC != "bad" & Jmax_QC != "bad"), aes(Vcmax_Estimate_final, `Vcmax(25)`, fill = month))+
  geom_point(shape = 21)+
  geom_abline()+
  theme_bw()+
  scale_x_continuous(name = "Vcmax25 (Plantecophys - varying TPU)")+
  scale_y_continuous(name = "Vcmax25 (Shiny)")

ggpubr::ggarrange(A, B, C, common.legend = TRUE, labels = "AUTO")

ggsave("./Outputs/Plots/frenchguiana_Vcmax_fitcomparison.png", type = "Cairo", bg = "white", dpi = 300, width = 7, height = 7)

#now plot these with leaf age
ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, TPU25, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_TPUlamour_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)

ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, `Vcmax(25)`, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_Vcmaxshiny_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)

ggplot(data = subset(data, !is.na(Age_leaf_mid)), aes(Age_leaf_mid, Vcmax25, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_Vcmaxlamour_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)

ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, Vcmax_Estimate_final, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_Vcmaxkalinotpu_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)

ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, Jmax_Estimate_final/Vcmax_Estimate_final, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_JVratiokalinotpu_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)


ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, `Jmax(25)`/`Vcmax(25)`, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_JVratioshiny_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)

ggplot(data = subset(data, !is.na(Age_leaf_max) & Vcmax_QC != "bad" & Jmax_QC != "bad"), aes(Age_leaf_max, Jmax25/Vcmax25, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_JVratiolamour_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)

#note JVratio is relatively consistent across the different fits.

#now plot TPU to Vcmax ratio
ggplot(data = subset(data, !is.na(Age_leaf_max) & Vcmax_QC != "bad" & Jmax_QC != "bad"), aes(Age_leaf_max, TPU25/Vcmax25, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_TPUVratiolamour_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)


ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, `Jmax(25)`, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_Jmaxshiny_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)

ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, Jmax25, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_Jmaxslamour_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)

ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, Jmax_Estimate, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_Jmaxkalinotpu_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)


# OTHER TRAITS

ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, `fv_/_fm`, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_fvfm_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)


#TLP traits
ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_mid, PTLP))+
  geom_point(aes(fill = Age_leaf_range), shape = 21, size = 2)+
  geom_smooth(col = "black", method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  facet_wrap(~sp, ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "top")+
  scale_fill_viridis_c(option = "C")+
  labs(fill = "Leaf age uncertainty (days)", y = "TLP", x = "Estimated Leaf age (days)")

#examploe of effect of leaf range
ggsave("./Outputs/Plots/frenchguiana_TLP_leafagerange.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)


ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_mid, PTLP, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_TLP_leafage.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)
#for many it seems PTLP declines and then either stabilises or increases - except for Pouteria_eugenifolia and Sextonia rubra

ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(julian_date, Age_leaf_range, col = cohort))+
         geom_point()+
  facet_wrap(sp~tree)


##### OUT OF CURIOUSITY LET US PLOT THE OBSERVATION DATE INSTEAD OF LEAF AGE
ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(date, PTLP, col = sp, fill = sp, group = sp))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme_bw()+
  theme(legend.position = "none")

ggsave("./Outputs/Plots/frenchguiana_TLP_date.png", type = "Cairo", bg = "white", dpi = 300, width = 10, height = 5)


ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(date, PTLP, group = sp))+
  geom_point(aes(col = Age_leaf_max))+
  geom_smooth(col = "black", method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  facet_wrap(~sp, scales = "free_y")+
  theme(legend.position = "top")+
  theme_bw()+
  scale_color_viridis_c()

ggplot(data = subset(data, !is.na(Age_leaf_max)), aes(Age_leaf_max, PTLP, group = sp))+
  geom_point(aes(col = date), size = 3)+
  geom_smooth(col = "black", method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  facet_wrap(~sp, scales = "free_y")+
  theme(legend.position = "top")+
  theme_bw()+
  scale_color_viridis_c()







##### Calculate mean and se for each collection (different leaves for each group)

library(purrr)

# Identify numeric variables
numeric_vars <- data %>%
  select(where(is.numeric)) %>%
  names()

# Group and summarise
summary_stats <- data %>%
  group_by(sp, date, tree, cohort) %>%
  summarise(across(all_of(numeric_vars),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
                        n = ~sum(!is.na(.x))),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")

ggplot(data = subset(summary_stats, !is.na(Age_leaf_max_mean)), aes(Age_leaf_max_mean, PTLP_mean, group = sp))+
  geom_point(size = 1.5, shape = 21)+
  geom_errorbar(aes(ymin = PTLP_mean - PTLP_se, ymax = PTLP_mean + PTLP_se))+
  geom_smooth(col = "black", method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  facet_wrap(~sp, scales = "free_y")+
  theme(legend.position = "top")+
  theme_bw()

ggplot(data = subset(summary_stats, !is.na(Age_leaf_max_mean)), aes(Age_leaf_max_mean, `Vcmax(25)_mean`, group = sp))+
  geom_point(size = 1.5, shape = 21)+
  geom_errorbar(aes(ymin = `Vcmax(25)_mean` - `Vcmax(25)_se`, ymax = `Vcmax(25)_mean` + `Vcmax(25)_se`))+
  geom_smooth(col = "black", method = "gam", formula = y ~ s(x, k = 4, bs = "cs"))+
  facet_wrap(~sp, scales = "free_y")+
  theme(legend.position = "top")+
  theme_bw()


library(ggplot2)

##### LET US FIT THE GAM WITH RANDOM EFFECTS I.E. LEAF NESTED WITHIN TREE

library(mgcv)

#data <- rename(data, "Vcmax25" = `Vcmax(25)`, "Jmax25" = `Jmax(25)`)

# Make sure 'tree' and 'leaf' are factors
data$tree <- as.factor(data$tree)
data$leaf <- as.factor(data$leaf)
data$sp <- as.factor(data$sp)

# Fit the GAM model
gam_vcmax <- gam(Vcmax25 ~ s(Age_leaf_max, by = sp, bs = "cs", k = 6) +
                   sp +
                   s(tree, bs = "re"),
                 data = data,
                 method = "REML")


# Create a prediction dataset
newdata <- data %>%
  filter(!is.na(Age_leaf_max)) %>%
  group_by(sp) %>%
  summarise(Age_leaf_max = seq(min(Age_leaf_max, na.rm = TRUE),
                               max(Age_leaf_max, na.rm = TRUE),
                               length.out = 100),
            .groups = "drop") %>%
  # Add dummy tree and leaf levels (needed for random effect structure)
  mutate(tree = factor("dummy_tree"),
         leaf = factor("dummy_leaf"))

# Predict from the fitted model
newdata$pred <- predict(gam_vcmax, newdata = newdata, exclude = c("s(tree)", "s(leaf,tree)"))

library(ggplot2)

ggplot(data = subset(data, !is.na(Age_leaf_max)), 
       aes(x = Age_leaf_max, y = Vcmax25, col = sp, fill = sp)) +
  geom_point(alpha = 0.5) +
  geom_line(data = newdata, aes(x = Age_leaf_max, y = pred), inherit.aes = FALSE, colour = "black") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~sp, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "top")


library(mgcv)
library(ggplot2)
library(dplyr)

# Rename columns if needed
#data <- rename(data, "Vcmax25" = `Vcmax(25)`, "Jmax25" = `Jmax(25)`)

# Ensure grouping variables are factors
data$tree <- as.factor(data$tree)
data$leaf <- as.factor(data$leaf)
data$sp <- as.factor(data$sp)

# Fit the GAM model
gam_vcmax <- gam(Vcmax25 ~ s(Age_leaf_max, by = sp, bs = "cs", k = 6) +
                   sp +
                   s(tree, bs = "re"),
                 data = data,
                 method = "REML")
summary(gam_vcmax)

# Create a new prediction dataset
newdata <- data %>%
  filter(!is.na(Age_leaf_max)) %>%
  group_by(sp) %>%
  summarise(Age_leaf_max = seq(min(Age_leaf_max, na.rm = TRUE),
                               max(Age_leaf_max, na.rm = TRUE),
                               length.out = 100),
            .groups = "drop") %>%
  mutate(tree = factor("dummy_tree"),
         leaf = factor("dummy_leaf"))

# Predict with standard errors
predictions <- predict(gam_vcmax, newdata = newdata,
                       exclude = c("s(tree)", "s(leaf,tree)"),
                       se.fit = TRUE)

# Add predictions and confidence intervals to newdata
newdata$pred <- predictions$fit
newdata$se <- predictions$se.fit
newdata$lower <- newdata$pred - 1.96 * newdata$se
newdata$upper <- newdata$pred + 1.96 * newdata$se

# Plot with confidence ribbons
ggplot(data = subset(data, !is.na(Age_leaf_max)), 
       aes(x = Age_leaf_max, y = Vcmax25, col = sp, fill = sp)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(data = newdata, aes(x = Age_leaf_max, ymin = lower, ymax = upper, group = sp),
              inherit.aes = FALSE, alpha = 0.2, fill = "grey7") +
  geom_line(data = newdata, aes(x = Age_leaf_max, y = pred), inherit.aes = FALSE, colour = "black") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~sp, ncol = 4, nrow = 2) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")


###

# Fit the GAM model
gam_LMA_leaf <- gam(LMA_leaf ~ s(Age_leaf_max, by = sp, bs = "cs", k = 8) +
                   sp +
                   s(tree, bs = "re"),
                 data = data,
                 method = "REML")

summary(gam_LMA_leaf)

# Create a new prediction dataset
newdata <- data %>%
  filter(!is.na(Age_leaf_max)) %>%
  group_by(sp) %>%
  summarise(Age_leaf_max = seq(min(Age_leaf_max, na.rm = TRUE),
                               max(Age_leaf_max, na.rm = TRUE),
                               length.out = 100),
            .groups = "drop") %>%
  mutate(tree = factor("dummy_tree"),
         leaf = factor("dummy_leaf"))

# Predict with standard errors
predictions <- predict(gam_LMA_leaf, newdata = newdata,
                       exclude = c("s(tree)", "s(leaf,tree)"),
                       se.fit = TRUE)

# Add predictions and confidence intervals to newdata
newdata$pred <- predictions$fit
newdata$se <- predictions$se.fit
newdata$lower <- newdata$pred - 1.96 * newdata$se
newdata$upper <- newdata$pred + 1.96 * newdata$se

# Plot with confidence ribbons
ggplot(data = subset(data, !is.na(Age_leaf_max)), 
       aes(x = Age_leaf_max, y = LMA_leaf, col = sp, fill = sp)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(data = newdata, aes(x = Age_leaf_max, ymin = lower, ymax = upper, group = sp),
              inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_line(data = newdata, aes(x = Age_leaf_max, y = pred), inherit.aes = FALSE, colour = "black") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~sp, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "top")

### INCLUDE DAY OF YEAR AS A COVARIATE IN THE GAM

# Fit the GAM model including seasonal effect using Julian date
gam_vcmax <- gam(Vcmax25 ~ s(Age_leaf_max, by = sp, bs = "cs", k = 5) +
                   s(julian_date, bs = "cc", k = 5) +  # cyclic spline for seasonality
                   sp +
                   s(tree, bs = "re"),
                 data = data,
                 method = "REML")

summary(gam_vcmax)

# Create new prediction data
newdata <- data %>%
  filter(!is.na(Age_leaf_max)) %>%
  group_by(sp) %>%
  summarise(Age_leaf_max = seq(min(Age_leaf_max, na.rm = TRUE),
                               max(Age_leaf_max, na.rm = TRUE),
                               length.out = 100),
            .groups = "drop") %>%
  mutate(tree = factor("dummy_tree"),
         julian_date = median(data$julian_date, na.rm = TRUE))  # fixed seasonality for prediction

# Predict with standard errors
predictions <- predict(gam_vcmax, newdata = newdata,
                       exclude = c("s(tree)", "s(leaf,tree)"),
                       se.fit = TRUE)

# Add predictions and confidence intervals
newdata$pred <- predictions$fit
newdata$se <- predictions$se.fit
newdata$lower <- newdata$pred - 1.96 * newdata$se
newdata$upper <- newdata$pred + 1.96 * newdata$se

# Plot
ggplot(data = subset(data, !is.na(Age_leaf_max)), 
       aes(x = Age_leaf_max, y = Vcmax25, col = sp, fill = sp)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(data = newdata, aes(x = Age_leaf_max, ymin = lower, ymax = upper, group = sp),
              inherit.aes = FALSE, alpha = 0.2, fill = "grey7") +
  geom_line(data = newdata, aes(x = Age_leaf_max, y = pred), inherit.aes = FALSE, colour = "black") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~sp, ncol = 4, nrow = 2) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")

