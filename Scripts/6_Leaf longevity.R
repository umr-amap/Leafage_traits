# Some Code to Try to Estimate leaf longevity
# Uses the optimality-based approach to predict leaf longevity as the result of a cost-benefit trade-off between
# Construction and maintenance costs of the leaf (e.g. building and supporting tissues)
# Carbon gain through photosynthesis over the leaf's lifespan.
#The central assumption is that a leaf should be retained only as long as it yields net carbon gain


# For each tree we want the peak for the traits LMA, Vcmax

library(readr)
library(ggplot2)
library(dplyr)
library(mgcv)
library(purrr)
library(tibble)

data <- read_csv("Outputs/Combined_data.csv")

#quick plot to look at peak Vcmax and LMA

ggplot(data, aes(Age_leaf_mid, Vcmax25, col = sp))+
  geom_point()+
  facet_wrap(~sp)

ggplot(data, aes(Age_leaf_mid, LMA_leaf, col = sp))+
  geom_point()+
  facet_wrap(~sp)

################################################################################
# IDENTIFY THE TRAIT VALUES AT MATURATION AND THE LEAF AGE AT WHICH THEY PEAK

#for now let us just use the 90th % for the mature vals.
data_species <- data %>%
  group_by(sp) %>%
  summarise(
    LMA_leaf_90 = quantile(LMA_leaf, probs = 0.9, na.rm = TRUE),
    Vcmax25_90 = quantile(Vcmax25, probs = 0.9, na.rm = TRUE),
    Vcmax25mass_90 = Vcmax25_90 / LMA_leaf_90,
    max_age_leaf_mid = max(Age_leaf_mid, na.rm = TRUE),
    .groups = "drop"
  )

#Don't get the age using this method

################################################################################

# Test doing this using GAMS

# Function to fit GAM and extract peak age and trait value
find_trait_peak_gam <- function(df, trait_col, age_col = "Age_leaf_mid") {
  df <- df %>% filter(!is.na(.data[[trait_col]]), !is.na(.data[[age_col]]))
  
  if (nrow(df) < 5 || length(unique(df[[age_col]])) < 4) {
    return(data.frame(peak_age = NA, peak_value = NA))
  }
  
  fit <- gam(as.formula(paste(trait_col, "~ s(", age_col, ", k = 4, bs = 'cs')")), data = df)
  
  age_seq <- seq(min(df[[age_col]], na.rm = TRUE),
                 max(df[[age_col]], na.rm = TRUE),
                 length.out = 200)
  preds <- predict(fit, newdata = data.frame(Age_leaf_mid = age_seq))
  
  peak_index <- which.max(preds)
  peak_age <- age_seq[peak_index]
  peak_value <- preds[peak_index]
  
  return(data.frame(peak_age = peak_age, peak_value = peak_value))
}


# Apply to both traits for each species
trait_peaks <- data %>%
  filter(!is.na(sp)) %>%
  group_by(sp) %>%
  group_split() %>%
  map_df(function(df) {
    lma_peak <- find_trait_peak_gam(df, trait_col = "LMA_leaf")
    vcmax_peak <- find_trait_peak_gam(df, trait_col = "Vcmax25")
    
    tibble(
      sp = unique(df$sp),
      peak_LMA_age = lma_peak$peak_age,
      peak_LMA_value = lma_peak$peak_value,
      peak_Vcmax25_age = vcmax_peak$peak_age,
      peak_Vcmax25_value = vcmax_peak$peak_value
    )
  })


# View the summary
trait_peaks

hist(trait_peaks$peak_LMA_age)
hist(trait_peaks$peak_LMA_value)
hist(trait_peaks$peak_Vcmax25_age)
hist(trait_peaks$peak_Vcmax25_value)

#add these values to main data for plotting
data <- merge(data, trait_peaks, by = "sp")

# add these values to the plots
#PLOTTING LMA
ggplot(data, aes(x = Age_leaf_mid, y = LMA_leaf)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~sp)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"), se = FALSE, colour = "darkgreen") +
  geom_vline(aes(xintercept = peak_LMA_age), linetype = "dashed", colour = "blue") +
  geom_hline(aes(yintercept = peak_LMA_value), linetype = "dashed", colour = "blue") +
  geom_point(aes(x = peak_LMA_age, y = peak_LMA_value),
             colour = "red", size = 3) +
  labs(title = "LMA vs Age",
       x = "Leaf Age (mid, days)",
       y = "LMA (g m⁻²)")


#PLOTTING VCMAX25
ggplot(data, aes(x = Age_leaf_mid, y = Vcmax25)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~sp)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"), se = FALSE, colour = "darkgreen") +
  geom_vline(aes(xintercept = peak_Vcmax25_age), linetype = "dashed", colour = "blue") +
  geom_hline(aes(yintercept = peak_Vcmax25_value), linetype = "dashed", colour = "blue") +
  geom_point(aes(x = peak_Vcmax25_age, y = peak_Vcmax25_value),
             colour = "red", size = 3) +
  labs(title = "Vcmax25 vs Age",
       x = "Leaf Age (mid, days)",
       y = "Vcmax25")


#### Ok updating method to find the first maximum using derivative instead of the absoluate maximum. 
# Approximates the slope and look for where the slope shifts from positive to negative.
# Finds first stabilisation point (plateau onset).
# Falls back to which.max if no local max if found.

# Logging skipped species
skipped_species <- c()

# Function to find first local maximum in smoothed GAM curve
# Note this finds the first minimum in the curve for specified traitsd i.e. TLP

find_first_local_peak <- function(df, trait_col, age_col = "Age_leaf_mid", sp_name = NULL, find_min = FALSE) {
  df <- df %>% filter(!is.na(.data[[trait_col]]), !is.na(.data[[age_col]]))
  
  if (nrow(df) < 5 || length(unique(df[[age_col]])) < 4) {
    if (!is.null(sp_name)) skipped_species <<- c(skipped_species, sp_name)
    return(data.frame(peak_age = NA, peak_value = NA))
  }
  
  fit <- gam(as.formula(paste(trait_col, "~ s(", age_col, ", k = 4, bs = 'cs')")), data = df)
  age_seq <- seq(min(df[[age_col]], na.rm = TRUE),
                 max(df[[age_col]], na.rm = TRUE),
                 length.out = 200)
  preds <- predict(fit, newdata = data.frame(Age_leaf_mid = age_seq))
  slope <- diff(preds) / diff(age_seq)
  
  if (find_min) {
    peak_index <- which(diff(sign(slope)) == 2)[1]  # slope goes from negative to positive
    if (is.na(peak_index)) peak_index <- which.min(preds)
  } else {
    peak_index <- which(diff(sign(slope)) == -2)[1]  # slope goes from positive to negative
    if (is.na(peak_index)) peak_index <- which.max(preds)
  }
  
  peak_age <- age_seq[peak_index]
  peak_value <- preds[peak_index]
  return(data.frame(peak_age = peak_age, peak_value = peak_value))
}


# Apply to all species and traits

selected_traits <- c("LMA_leaf", "Vcmax25", "Nmass", "leaf_thickness", "leaf_drymass", "leaf_area","PTLP")

invert_peak_traits <- c("PTLP")  # you can add more if needed

trait_peaks <- data %>%
  filter(!is.na(sp)) %>%
  group_by(sp) %>%
  group_split() %>%
  map_df(function(df) {
    sp_name <- unique(df$sp)
    out <- tibble(sp = sp_name)
    
    for (trait in selected_traits) {
      find_min <- trait %in% invert_peak_traits
      peak <- find_first_local_peak(df, trait_col = trait, sp_name = sp_name, find_min = find_min)
      
      out[[paste0("peak_", trait, "_age")]] <- peak$peak_age
      out[[paste0("peak_", trait, "_value")]] <- peak$peak_value
    }
    
    if ("LMA_leaf" %in% selected_traits && "Vcmax25" %in% selected_traits) {
      df_lma <- df %>% filter(!is.na(LMA_leaf), !is.na(Age_leaf_mid))
      vcmax_peak_age <- out$`peak_Vcmax25_age`
      
      lma_at_vcmax_age <- NA
      if (nrow(df_lma) >= 5 && length(unique(df_lma$Age_leaf_mid)) >= 4 && !is.na(vcmax_peak_age)) {
        lma_gam <- gam(LMA_leaf ~ s(Age_leaf_mid, k = 4, bs = "cs"), data = df_lma)
        lma_at_vcmax_age <- predict(lma_gam, newdata = data.frame(Age_leaf_mid = vcmax_peak_age))
      }
      
      out$LMA_at_Vcmax25_peak <- lma_at_vcmax_age
      out$Vcmax25mass_peak <- out$`peak_Vcmax25_value` / lma_at_vcmax_age
    }
    
    out$max_Age_leaf_mid <- max(df$Age_leaf_mid, na.rm = TRUE)
    out$max_Age_leaf_max <- max(df$Age_leaf_max, na.rm = TRUE)
    
    return(out)
  })


# Output summary
trait_peaks
skipped_species

# Estimate b from mass-based Vcmax25 (Vcmax25m in µmol g⁻¹ s⁻¹)
trait_peaks$estimate_b <- 10^(2.35 - 1.36 * log10(trait_peaks$Vcmax25mass_peak)) 
#based on tropical data from metaanalysis in Xu et al 2017

summary(trait_peaks$estimate_b) #ranges 410 - 1475

ggplot(trait_peaks, aes(Vcmax25mass_peak, estimate_b, fill = sp))+
  geom_point(size = 5, shape = 21)

ggplot(trait_peaks, aes(max_Age_leaf_max, estimate_b, fill = sp))+
  geom_point(size = 5, shape = 21)


#add these values to main data for plotting
data <- merge(data, trait_peaks, by = "sp")

# add these values to the plots
#PLOTTING LMA
ggplot(data, aes(x = Age_leaf_mid, y = LMA_leaf)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~sp)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"), se = FALSE, colour = "darkgreen") +
  geom_vline(aes(xintercept = peak_LMA_leaf_age), linetype = "dashed", colour = "blue") +
  geom_hline(aes(yintercept = peak_LMA_leaf_value), linetype = "dashed", colour = "blue") +
  geom_point(aes(x = peak_LMA_leaf_age, y = peak_LMA_leaf_value),
             colour = "red", size = 3) +
  geom_point(aes(x = peak_Vcmax25_age, y = LMA_at_Vcmax25_peak),
             colour = "hotpink", size = 3) + #LMA when Vcmax peaked.
  labs(title = "LMA vs Age",
       x = "Leaf Age (mid, days)",
       y = "LMA (g m⁻²)")


#PLOTTING VCMAX25
ggplot(data, aes(x = Age_leaf_mid, y = Vcmax25)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~sp)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"), se = FALSE, colour = "darkgreen") +
  geom_vline(aes(xintercept = peak_Vcmax25_age), linetype = "dashed", colour = "blue") +
  geom_hline(aes(yintercept = peak_Vcmax25_value), linetype = "dashed", colour = "blue") +
  geom_point(aes(x = peak_Vcmax25_age, y = peak_Vcmax25_value),
             colour = "red", size = 3) +
  labs(title = "Vcmax25 vs Age",
       x = "Leaf Age (mid, days)",
       y = "Vcmax25")

# PLOT DRY MASS
ggplot(data, aes(x = Age_leaf_mid, y = leaf_drymass)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~sp, scales = "free_y")+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"), se = FALSE, colour = "darkgreen") +
  geom_vline(aes(xintercept = peak_leaf_drymass_age), linetype = "dashed", colour = "blue") +
  geom_hline(aes(yintercept = peak_leaf_drymass_value), linetype = "dashed", colour = "blue") +
  geom_point(aes(x = peak_leaf_drymass_age, y = peak_leaf_drymass_value),
             colour = "red", size = 3) +
  labs(title = "leaf_drymass vs Age",
       x = "Leaf Age (mid, days)",
       y = "leaf_drymass")

# PLOT TLP
ggplot(data, aes(x = Age_leaf_mid, y = PTLP)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~sp, scales = "free_y")+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"), se = FALSE, colour = "darkgreen") +
  geom_vline(aes(xintercept = peak_PTLP_age), linetype = "dashed", colour = "blue") +
  geom_hline(aes(yintercept = peak_PTLP_value), linetype = "dashed", colour = "blue") +
  geom_point(aes(x = peak_PTLP_age, y = peak_PTLP_value),
             colour = "red", size = 3) +
  labs(title = "TLP vs Age",
       x = "Leaf Age (mid, days)",
       y = "TLP")

#lets compare the leaf maturation ages for these traits
ggplot(trait_peaks, aes(peak_LMA_leaf_age, peak_Vcmax25_age, col = sp))+
  geom_point(size = 3)+
  geom_abline()

ggplot(trait_peaks, aes(peak_LMA_leaf_age, peak_PTLP_age, col = sp))+
  geom_point(size = 3)+
  geom_abline()

ggplot(trait_peaks, aes(peak_leaf_drymass_age, peak_PTLP_age, col = sp))+
  geom_point(size = 3)+
  geom_abline()

ggplot(trait_peaks, aes(peak_leaf_thickness_age, peak_PTLP_age, col = sp))+
  geom_point(size = 3)+
  geom_abline()

ggplot(trait_peaks, aes(peak_Vcmax25_age, peak_PTLP_age, col = sp))+
  geom_point(size = 3)+
  geom_abline()

#note that this predicts that Vcmax peaks earlier than LMA. is this what we expect? yes
#also, when calculating Vcmaxmass basis for each species, these will be based off different ages, or should we standardise it
#and use the LMA from when Vcmax peaked?

ggplot(trait_peaks, aes(peak_LMA_leaf_value, peak_Vcmax25_value, col = sp))+
  geom_point(size = 3)
#note no correlation between the peak LMA and peak Vcmax value for these species, similar in Xu paper.


#now that we have an estimate for b (which could roughly be equivalent to leaf longevity)
#lets use this to calcualte relative leaf age - though note this is a simplistic way of doing it.
data$Age_leaf_mid_rel <- data$Age_leaf_mid / data$estimate_b

#as a quick comparison - lets do this also for trait peaks to see if the range of peaks for Vcmax or LMA is reduced when using relative age
trait_peaks$peak_LMA_leaf_age_rel <- trait_peaks$peak_LMA_leaf_age / trait_peaks$estimate_b
trait_peaks$peak_Vcmax25_age_rel <- trait_peaks$peak_Vcmax25_age / trait_peaks$estimate_b


#note that 0 = newly mature leaf and 1 = end of functional lifespan.
#note that our age to reach maturation is based on our data not assuming 2 months like they have

summary(trait_peaks$peak_LMA_age)
summary(trait_peaks$peak_LMA_age_rel)

summary(trait_peaks$peak_Vcmax25_age)
summary(trait_peaks$peak_Vcmax25_age_rel)

#PLOTTING LMA
ggplot(data, aes(x = Age_leaf_mid_rel, y = LMA_leaf)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~sp)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"), se = FALSE, colour = "darkgreen") +
  #geom_vline(aes(xintercept = peak_LMA_age), linetype = "dashed", colour = "blue") +
  #geom_hline(aes(yintercept = peak_LMA_value), linetype = "dashed", colour = "blue") +
  #geom_point(aes(x = peak_LMA_age, y = peak_LMA_value),
  #           colour = "red", size = 3) +
  #geom_point(aes(x = peak_Vcmax25_age, y = LMA_at_Vcmax25_peak),
  #           colour = "hotpink", size = 3) + #LMA when Vcmax peaked.
  labs(title = "LMA vs Age",
       x = "Relative Leaf Age (mid, days)",
       y = "LMA (g m⁻²)")


#PLOTTING VCMAX25
ggplot(data, aes(x = Age_leaf_mid_rel, y = Vcmax25)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~sp)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"), se = FALSE, colour = "darkgreen") +
  #geom_vline(aes(xintercept = peak_Vcmax25_age), linetype = "dashed", colour = "blue") +
  #geom_hline(aes(yintercept = peak_Vcmax25_value), linetype = "dashed", colour = "blue") +
  #geom_point(aes(x = peak_Vcmax25_age, y = peak_Vcmax25_value),
  #           colour = "red", size = 3) +
  labs(title = "Vcmax25 vs Age",
       x = "Relative Leaf Age (mid, days)",
       y = "Vcmax25")

#plot another trait
#PLOTTING VCMAX25
ggplot(data, aes(x = Age_leaf_mid_rel, y = PTLP)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~sp, scales = "free_y")+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"), se = FALSE, colour = "darkgreen") +
  geom_vline(aes(xintercept = peak_PTLP_age/estimate_b), linetype = "dashed", colour = "blue") +
  geom_hline(aes(yintercept = peak_PTLP_value/estimate_b), linetype = "dashed", colour = "blue") +
  geom_point(aes(x = peak_PTLP_age/estimate_b, y = peak_PTLP_value),
             colour = "red", size = 3) +
  labs(title = "TLP vs Age",
       x = "Relative Leaf Age (mid, days)",
       y = "TLP")

ggplot(data, aes(x = Age_leaf_mid, y = PTLP)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~sp, scales = "free_y")+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"), se = FALSE, colour = "darkgreen") +
  geom_vline(aes(xintercept = peak_PTLP_age), linetype = "dashed", colour = "blue") +
  geom_hline(aes(yintercept = peak_PTLP_value), linetype = "dashed", colour = "blue") +
  geom_point(aes(x = peak_PTLP_age, y = peak_PTLP_value),
             colour = "red", size = 3) +
  labs(title = "TLP vs Age",
       x = "Absolute Leaf Age (mid, days)",
       y = "TLP")

ggplot(data, aes(x = Age_leaf_max, y = leaf_thickness)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~sp, scales = "free_y")+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"), se = FALSE, colour = "darkgreen") +
  #geom_vline(aes(xintercept = peak_Vcmax25_age), linetype = "dashed", colour = "blue") +
  #geom_hline(aes(yintercept = peak_Vcmax25_value), linetype = "dashed", colour = "blue") +
  #geom_point(aes(x = peak_Vcmax25_age, y = peak_Vcmax25_value),
  #           colour = "red", size = 3) +
  labs(title = "area vs Age",
       x = "Absolute Leaf Age (mid, days)",
       y = "area")


