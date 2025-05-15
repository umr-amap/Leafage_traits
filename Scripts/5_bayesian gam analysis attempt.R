# PACKAGES
library(dplyr)
library(tidyr)
library(readr)
library(brms)
library(tidyverse)
library(ggplot2)
library(mgcv)

#DATA
data <- read_csv("Outputs/Combined_data.csv")

# OPTIONAL STEP TO EXCLUDE OUTLIERS HERE
data <- data %>%
  mutate(
    LMA_leaf = ifelse(LMA_leaf_QC_dist > 0, NA_real_, LMA_leaf),
    leaf_thickness = ifelse(thickness_QC_dist > 0, NA_real_, leaf_thickness),
    leaf_area = ifelse(leaf_area_QC_dist > 0, NA_real_, leaf_area),
    leaf_drymass = ifelse(leaf_drymass_QC_dist > 0, NA_real_, leaf_drymass),
    Narea = ifelse(Narea_QC_dist > 0, NA_real_, Narea),
    Nmass = ifelse(Nmass_QC_dist > 0, NA_real_, Nmass),
    Cmass = ifelse(Cmass_QC_dist > 0, NA_real_, Cmass),
    PTLP = ifelse(PTLP_QC_dist > 0, NA_real_, PTLP),
    Vcmax_Estimate_final = ifelse(Vcmax_Estimate_final_QC_dist > 0, NA_real_, Vcmax_Estimate_final),
    Jmax_Estimate_final = ifelse(Jmax_Estimate_final_QC_dist > 0, NA_real_, Jmax_Estimate_final),
    TPU_Estimate_final = ifelse(TPU_Estimate_final_QC_dist > 0, NA_real_, TPU_Estimate_final)
  )


# -----------------------------------------------
# First try using GAMs - not bayesian, to disentangle season vs leaf age and identify species effects.
# Use leaf_age_mid for now
# -----------------------------------------------

# ------------------------------------------
# Step 1: Load and Prepare the Data
# ------------------------------------------

# Assume data is in a dataframe called 'data'
# You should load your data here (use read.csv() or similar for your dataset)
# For example:
# data <- read.csv("your_data_file.csv")

# Simplified model: Use leaf_age_mid for leaf age and julian_date for seasonality
# Ensure that leaf_age_mid and julian_date are numeric and available in the data.

data$sp <- as.factor(data$sp)
data$tree <- as.factor(data$tree)

data <- subset(data, !is.na(Age_leaf_mid))

# ------------------------------------------
# Step 2: Define the Model Formula
# ------------------------------------------

# First check for correlation

fit0 <- mgcv::gam(Age_leaf_mid ~ s(julian_date, by = sp, bs = "cc", k = 4) + 
                    s(tree, bs = "re"), 
                  data = data, family = gaussian())

summary(fit0)
plot(fit0, pages = 1)  # Plot smooth terms for both predictors (leaf age and julian date)

# We include species-specific smooth terms for leaf age (Age_leaf_mid) and seasonality (julian_date).
# We also include a random intercept for tree to account for within-tree variability.

# Full Model: Species-specific smooth terms for leaf age and julian_date
# The random intercept for tree captures tree-level variability.

# Fit the model using mgcv::gam
fit1 <- mgcv::gam(LMA_leaf ~ s(Age_leaf_mid, by = sp, k = 4) + 
                   s(julian_date, by = sp, bs = "cc", k = 4) + 
                   s(tree, bs = "re"), 
                  data = data, family = gaussian()) 
#species specific smooths for both leaf age and seasonality and tree included as random effect

fit2 <- mgcv::gam(LMA_leaf ~ s(Age_leaf_mid, k = 4) + 
                    s(julian_date, bs = "cc", k = 4) + 
                    s(tree, bs = "re") + 
                    s(sp, bs = "re"), 
                  data = data, family = gaussian())
#species included as a fixed effect but doesnt have species specific smooths

fit3 <- mgcv::gam(LMA_leaf ~ s(Age_leaf_mid, k = 4) + 
                    s(julian_date, bs = "cc", k = 4) + 
                    s(tree, bs = "re"), 
                  data = data, family = gaussian())
#species not included at all

fit4 <- mgcv::gam(LMA_leaf ~ s(Age_leaf_mid, by = sp, k = 4) + 
                    s(julian_date, bs = "cc", k = 4) + 
                    s(tree, bs = "re"), 
                  data = data, family = gaussian()) 
#species-specific smooths only for leaf age

fit5 <- mgcv::gam(LMA_leaf ~ s(Age_leaf_mid, k = 4) + 
                    s(julian_date, by = sp, bs = "cc", k = 4) + 
                    s(tree, bs = "re"), 
                  data = data, family = gaussian()) 
#species-specific smooths only for seasonality

fit6 <- mgcv::gam(
  LMA_leaf ~ s(Age_leaf_mid, by = sp, k = 3) + 
    s(julian_date, bs = "cc", k = 3) + 
    s(Age_leaf_mid, julian_date, by = sp, k = 3), 
  data = data, 
  family = gaussian()
)


fit7 <-   mgcv::gam(
  LMA_leaf ~ s(Age_leaf_mid, by = sp, k = 3) + 
    s(julian_date, by = sp, bs = "cc", k = 3) + 
    s(Age_leaf_mid, julian_date, by = sp, k = 3), 
  data = data, 
  family = gaussian()
)

#includes the interaction between leaf age and seasonality since they can be correlated due to the timing of leaf expansion

#compare models
anova(fit1, fit4, fit5, fit3, fit6, fit7, test = "Chisq") #model indicates species specific smooths are important for seasonality but not for leaf age.
#strong evidence for a species-specific seasonality effect.

#indicates species should be included as a fixed effect. let us check if it improves it to have by= sp just for lewaf age or just for seasonality



# ------------------------------------------
# Step 4: Model Diagnostics
# ------------------------------------------

# Check the model summary
summary(fit1)
summary(fit2)
summary(fit3)
summary(fit6)
summary(fit7)

fit <- fit7

# Visualize the smooth terms
plot(fit, pages = 1)  # Plot smooth terms for both predictors (leaf age and julian date)

# Check model residuals
par(mfrow = c(1, 2))
gam.check(fit)  # Check diagnostics: residuals, fitted vs residuals plot

# ------------------------------------------
# Step 5: Visualize the Effect of Leaf Age and Seasonality on Vcmax
# ------------------------------------------

# Create a data frame for predictions based on the fitted model
# Generate a sequence of leaf age and julian_date for plotting
age_seq <- seq(min(data$Age_leaf_mid), max(data$Age_leaf_mid), length.out = 100)
julian_seq <- seq(1, 365, length.out = 100)

# Create prediction grid
newdata <- expand.grid(Age_leaf_mid = age_seq, julian_date = julian_seq, sp = unique(data$sp))

# Predict values from the fitted GAM
preds <- predict(fit, newdata = newdata, type = "response", se.fit = TRUE)

# Add the predictions to the new data frame
newdata$pred <- preds$fit
newdata$se <- preds$se.fit

# Create a plot of the effect of leaf age on Vcmax
ggplot(newdata, aes(x = Age_leaf_mid, y = pred)) +
  geom_line(aes(color = sp), size = 1) +
  geom_ribbon(aes(ymin = pred - 1.96 * se, ymax = pred + 1.96 * se, fill = sp), alpha = 0.2) +
  labs(x = "Leaf Age (days)", y = "Vcmax Estimate", title = "Effect of Leaf Age on Vcmax") +
  theme_minimal() +
  theme(legend.position = "none")

# Create a plot of the effect of seasonality (julian_date) on Vcmax
ggplot(newdata, aes(x = julian_date, y = pred)) +
  geom_line(aes(color = sp), size = 1) +
  geom_ribbon(aes(ymin = pred - 1.96 * se, ymax = pred + 1.96 * se, fill = sp), alpha = 0.2) +
  labs(x = "Julian Date", y = "Vcmax Estimate", title = "Effect of Seasonality on Vcmax") +
  theme_minimal() +
  theme(legend.position = "none")

# ------------------------------------------
# Step 6: Extract Species-Specific Smooth Effects and Check for Significance
# ------------------------------------------

# Extract the smooth terms for species
species_smooths <- fit$smooth[[1]]  # Smooth term for Age_leaf_mid by species
print(species_smooths)

# Check if the species-specific effects are significant
# Look at the 95% CI (if the CI does not include zero, the effect is likely significant)
# You can examine these using summary(fit) or by plotting the smooth terms directly.

# ------------------------------------------
# Step 7: Save the Model Results
# ------------------------------------------

# Save the model output
saveRDS(fit, "species_vcmax_model.rds")

# Save the predictions and summary data
write.csv(newdata, "predictions_with_vcmax.csv")



#################################################################################
#################################################################################

# -----------------------------------------------
# Bayesian GAMM for trait-age response (Vcmax25)
# -----------------------------------------------

# STEP 1: Simulate multiple plausible leaf ages
# For each observation, we sample 10 possible leaf ages uniformly from [min, max]
set.seed(123)  # for reproducibility

data_long <- data %>%
  filter(!is.na(Vcmax_Estimate_final), !is.na(Age_leaf_min), !is.na(Age_leaf_max)) %>%
  rowwise() %>%
  mutate(age_draws = list(runif(10, Age_leaf_min, Age_leaf_max))) %>%
  unnest(cols = c(age_draws))

# -----------------------------------------------
# STEP 2: Define the model formula
# We now include species-specific smooth terms for both leaf age and seasonality (julian_date).
# We use species-specific smooth terms via s(age_draws, by = sp) and s(julian_date, by = sp) for seasonality
# and a random intercept for tree
# -----------------------------------------------

formula <- bf(
  Vcmax_Estimate_final ~ s(age_draws, by = sp, k = 4) + s(julian_date, by = sp, bs = "cc", k = 4) + (1 | tree)
)

# -----------------------------------------------
# STEP 3: Set weakly informative priors
# These allow reasonable flexibility while avoiding overfitting
# -----------------------------------------------

prior <- c(
  prior(normal(50, 30), class = "Intercept"),
  prior(student_t(3, 0, 20), class = "sigma"),
  prior(student_t(3, 0, 20), class = "sd", group = "tree"),
  prior(student_t(3, 0, 20), class = "sds")  # for spline smoothness
)

# -----------------------------------------------
# STEP 4: Fit the model
# This may take some time depending on the size of your data
# -----------------------------------------------

#initial settings just to test the code works
fit <- brm(
  formula = formula,
  data = data_long,
  prior = prior,
  chains = 2,         # Use 2 chains for a quick test (instead of 4)
  iter = 1000,        # Use fewer iterations for a quick test (instead of 3000)
  warmup = 500,       # Use fewer warmup iterations for a quick test (instead of 1000)
  cores = 4,          # Use 4 cores for parallel processing (adjust based on your system)
  control = list(adapt_delta = 0.9),  # Use a slightly lower adapt_delta for testing
  seed = 123
)

#took 741 seconds

 fit <- brm(
   formula = formula,
   data = data_long,
   prior = prior,
   chains = 4,
   iter = 3000,
   warmup = 1000,
   cores = 4,
   control = list(adapt_delta = 0.95),  # stabilise for splines
   seed = 123
 )

# -----------------------------------------------
# STEP 5: Check diagnostics
# -----------------------------------------------

#what do i need to do and show here that the model does or does not work?
summary(fit)
library(bayesplot)
mcmc_trace(fit)
pairs(fit)

# -----------------------------------------------
# STEP 6: Visualise the independent effects of leaf age and seasonality on the trait
# -----------------------------------------------

# Extract the conditional effects of the model (this gives you smooth terms)
cond_effects <- conditional_effects(fit, robust = TRUE)

# Plot the effect of leaf age
ggplot(cond_effects$`age_draws:sp`, aes(x = age_draws, y = estimate__, ymin = lower__, ymax = upper__)) +
  geom_line(color = "blue") +
  geom_ribbon(alpha = 0.2) +
  facet_wrap(~sp, scales = "free", ncol = 4, nrow = 2) + 
  labs(x = "Leaf Age (days)", y = "Effect on Vcmax", title = "Effect of Leaf Age on Vcmax") +
  theme_minimal()

# Plot the effect of seasonality (julian_date)
ggplot(cond_effects$`julian_date:sp`, aes(x = julian_date, y = estimate__, ymin = lower__, ymax = upper__)) +
  geom_line(color = "green") +
  geom_ribbon(alpha = 0.2) +
  facet_wrap(~sp, ncol = 4, nrow = 2) + 
  labs(x = "Julian Date", y = "Effect on Vcmax", title = "Effect of Seasonality (Julian Date) on Vcmax") +
  theme_minimal()

#compare with the raw data
ggplot(data, aes(Age_leaf_mid, Vcmax_Estimate_final))+
  geom_point(aes(col = julian_date))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"))+
  facet_wrap(~sp, scales = "free", ncol = 4, nrow = 2)+
  scale_color_viridis_c()

#compare with the raw data
ggplot(data, aes(julian_date, Vcmax_Estimate_final))+
  geom_point(aes(col = Age_leaf_mid))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cc"))+
  #facet_wrap(~sp, scales = "free", ncol = 4, nrow = 2)+
  scale_color_viridis_c()

ggplot(data, aes(Age_leaf_mid, PTLP))+
  geom_point(aes(col = julian_date))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"))+
  facet_wrap(~sp, scales = "free", ncol = 4, nrow = 2)+
  scale_color_viridis_c()

#compare with the raw data
ggplot(data, aes(julian_date, PTLP))+
  geom_point(aes(col = Age_leaf_mid))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cc"))+
 # facet_wrap(~sp, scales = "free", ncol = 4, nrow = 2)+
  scale_color_viridis_c()

ggplot(data, aes(Age_leaf_mid, gsw))+
  geom_point(aes(col = julian_date))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"))+
  facet_wrap(~sp, scales = "free", ncol = 4, nrow = 2)+
  scale_color_viridis_c()

#compare with the raw data
ggplot(data, aes(julian_date, gsw))+
  geom_point(aes(col = Age_leaf_mid))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cc"))+
  # facet_wrap(~sp, scales = "free", ncol = 4, nrow = 2)+
  scale_color_viridis_c()

 # -----------------------------------------------
# STEP ?: FIT MODEL WITH NO SPECIES-SPECIFIC SEASONALITY TERM
# -----------------------------------------------
# Reduced Model (species as a random effect)
reduced_formula <- bf(
  Vcmax_Estimate_final ~ s(age_draws, k = 4) + s(julian_date, bs = "cc", k = 4) + (1 | sp/tree)
)

# Fit the reduced model
reduced_fit <- brm(
  formula = reduced_formula,
  data = data_long,
  prior = prior,
  chains = 2,
  iter = 1000,
  warmup = 500,
  cores = 4,
  control = list(adapt_delta = 0.9),
  seed = 123
)

# Compare models using Bayes Factor
bayes_factor <- bayes_factor(fit, reduced_fit)
print(bayes_factor)

#Interpretation: If the Bayes Factor is significantly greater than 1, 
#it suggests that the full model (with species-specific smooth terms) 
#fits the data significantly better than the reduced model (with species as a random effect).
#If the Bayes Factor is close to 1, it suggests that there is little difference 
#between the two models, and you might prefer the simpler reduced model.

# Compute LOO-CV for both models
full_loo <- loo(fit)
reduced_loo <- loo(reduced_fit)

# Compare LOO-CV results
loo_compare(full_loo, reduced_loo)
# Interpretation: The model with the lower LOO (lower expected log predictive density) 
# is generally preferred. If the full model has a significantly better LOO score, 
# it suggests that the additional complexity of modeling species-specific smooths 
# improves the predictive power. If the difference is small or in favor of the reduced model,
# it suggests that species-specific smooths are not necessary, and a random effect for species is sufficient.


# -----------------------------------------------
# STEP 6: Extract species-specific smooths and peaks
# We'll use conditional_effects() or posterior_epred() to do this
# That part comes after fitting completes
# -----------------------------------------------

# STEP 1: Define a sequence of leaf ages and a sequence of julian dates to predict over
age_seq <- seq(0, 300, length.out = 100)
julian_seq <- seq(1, 365, length.out = 300)  # assuming 365 days in a year

# STEP 2: Create a new data frame for predictions for each species
newdata <- expand_grid(
  age_draws = age_seq,
  julian_date = julian_seq,
  sp = unique(data_long$sp),
  tree = NA  # set to NA so predictions use population-level effect
)

# STEP 3: Get posterior predicted values
# posterior_epred returns the posterior mean prediction for each draw
# THIS STEP TAKES AGES

library(future)

# Set up parallel processing
plan(multisession, workers = 4)  # Adjust workers to the number of cores you have

# Run posterior_epred in parallel
preds <- posterior_epred(fit, newdata = newdata, allow_new_levels = TRUE)

# Reset to default sequential plan after use
plan(sequential)

# STEP 4: Convert to tidy data frame and compute mean & CI
pred_summary <- preds %>%
  as_tibble() %>%
  bind_cols(newdata) %>%
  group_by(sp, age_draws, julian_date) %>%
  summarise(
    mean = mean(c_across(starts_with("V"))),
    lower = quantile(c_across(starts_with("V")), 0.025),
    upper = quantile(c_across(starts_with("V")), 0.975),
    .groups = "drop"
  )

# STEP 5: Find the peak (maximum) Vcmax_Estimate_final for each species
peaks <- pred_summary %>%
  group_by(sp) %>%
  filter(mean == max(mean)) %>%
  slice(1) %>%  # in case of ties
  ungroup() %>%
  select(sp, peak_age = age_draws, peak_value = mean)

# STEP 6: Plot the results
ggplot(pred_summary, aes(x = age_draws, y = mean)) +
  geom_line(color = "darkgreen") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_point(data = peaks, aes(x = peak_age, y = peak_value), colour = "red", size = 3) +
  facet_wrap(~sp) +
  labs(x = "Leaf Age (days)", y = "Vcmax_Estimate_final (µmol m⁻² s⁻¹)",
       title = "Species-specific Vcmax_Estimate_final trajectories with peak points") +
  theme_minimal()

# Define the age sequence and prediction grid
age_seq <- seq(0, 400, length.out = 200)
newdata <- expand_grid(
  age_draws = age_seq,
  julian_date = julian_seq,  # also include seasonal variation
  sp = unique(data_long$sp),
  tree = NA  # population-level predictions
)

# Get posterior predictions: rows = posterior draws, columns = newdata rows
posterior_samples <- posterior_epred(fit, newdata = newdata)

# Combine posterior draws and metadata
nd_combined <- bind_cols(as_tibble(posterior_samples), newdata)

# For each species and draw, find the age where Vcmax_Estimate_final is maximal
peak_draws <- nd_combined %>%
  pivot_longer(cols = starts_with("V"), names_to = "draw", values_to = "pred") %>%
  group_by(draw, sp) %>%
  slice_max(order_by = pred, n = 1, with_ties = FALSE) %>%
  ungroup()

# Summarise to get posterior mean and credible intervals of peak age and value
peaks_summary <- peak_draws %>%
  group_by(sp) %>%
  summarise(
    peak_age_mean = mean(age_draws),
    peak_age_lower = quantile(age_draws, 0.025),
    peak_age_upper = quantile(age_draws, 0.975),
    peak_value_mean = mean(pred),
    peak_value_lower = quantile(pred, 0.025),
    peak_value_upper = quantile(pred, 0.975),
    .groups = "drop"
  )

# Save to CSV
write_csv(peaks_summary, "species_Vcmax_Estimate_final_peaks_with_CI.csv")
