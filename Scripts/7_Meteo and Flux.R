# Meteo and Flux data

library(readxl)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

data <- read_excel("Dataset/French_Guiana/Meteo and Eddy data/GX-METEO+EDDY-2004-2024E GF Propre - IM.xlsx")
data <- subset(data, Year >= 2020) #subset for time period we have phenology/trait data for.

#data cleaning for vpd
data <- subset(data, VPD55 > 0 )

#clean up column names
data <- data %>%
  mutate(julian_date = as.numeric(`Julian Day`),
         hour = as.numeric(`Hour/min`))  # Assuming this column contains hour and minute as decimal

#calculate daytime daily summaries
data.daily <- data %>%
  subset(Rg > 25) %>%
  group_by(Year, Month, julian_date) %>%
  summarise(GPP.max = max(`GPP-L`, na.rm = TRUE),
            GPP.sum = sum(`GPP-L`, na.rm = TRUE),
            GPP.mean = mean(`GPP-L`, na.rm = TRUE),
            #LE.max = max(LE_f, na.rm = TRUE),
            #LE.mean = mean(LE_f, na.rm = TRUE),
            Reco.max = max(`Reco-L`, na.rm = TRUE),
            Reco.mean = mean(`Reco-L`, na.rm = TRUE),
            Tair.max = max(`Temp(55)`, na.rm = TRUE),
            Tair.mean = mean(`Temp(55)`, na.rm = TRUE),
            VPD.max = max(VPD55, na.rm = TRUE),
            VPD.mean = mean(VPD55, na.rm = TRUE),
            SWC.max = max(SWC615, na.rm = TRUE),
            SWC.mean = mean(SWC615, na.rm = TRUE),
            RG.max = max(Rg, na.rm = TRUE),
            RG.mean = mean(Rg, na.rm = TRUE))

ggplot(data.daily, aes(julian_date, GPP.mean))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cc"), col = "black")+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank())

data.daily$Year <- as.factor(data.daily$Year)

a <- ggplot(data.daily, aes(julian_date, GPP.mean, col = Year, fill = Year))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cc"),
              aes(group = Year))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank())+
  scale_color_viridis_d(option = "A")+
  scale_fill_viridis_d(option = "A")
#same seasonal pattern, but timing of the peak GPP different across years. 

b <- ggplot(data.daily, aes(julian_date, Reco.mean, col = Year, fill = Year))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cc"),
              aes(group = Year))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank())+
  scale_color_viridis_d(option = "A")+
  scale_fill_viridis_d(option = "A")

c <- ggplot(data.daily, aes(julian_date, Tair.mean, col = Year, fill = Year))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cc"),
              aes(group = Year))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank())+
  scale_color_viridis_d(option = "A")+
  scale_fill_viridis_d(option = "A")

d <- ggplot(data.daily, aes(julian_date, VPD.mean, col = Year, fill = Year))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cc"),
              aes(group = Year))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank())+
  scale_color_viridis_d(option = "A")+
  scale_fill_viridis_d(option = "A")

e <- ggplot(data.daily, aes(julian_date, RG.mean, col = Year, fill = Year))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cc"),
              aes(group = Year))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank())+
  scale_color_viridis_d(option = "A")+
  scale_fill_viridis_d(option = "A")

f <- ggplot(data.daily, aes(julian_date, SWC.mean, col = Year, fill = Year))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cc"),
              aes(group = Year))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank())+
  scale_color_viridis_d(option = "A")+
  scale_fill_viridis_d(option = "A")

ggpubr::ggarrange(a, b, c, d, e, f, common.legend = TRUE, align = "hv",
                  ncol = 2, nrow = 3)

ggsave("./Outputs/Plots/Meteo_Fluxes/Meteoflux_doy.png", type = "Cairo",
       bg = "white", dpi = 300, width = 8, height = 10)

# check how correlated the meteo data is
ggplot(data.daily, aes(Tair.mean, VPD.mean))+
  geom_point()

ggplot(data.daily, aes(Tair.mean, SWC.mean))+
  geom_point()

ggplot(data.daily, aes(VPD.mean, SWC.mean))+
  geom_point()

ggplot(data.daily, aes(Tair.mean, RG.mean))+
  geom_point()

#could use Tair and SWC - since other traits not only Vcmax may be affected e.g. TLP

#### Now let us read in the trait measurement dates, so we can calculate what
# The mean conditions were over the 30 days prior to measurements.
traits <- read_csv("Outputs/Combined_data.csv")

# Function to calculate rolling means over the 2 weeks before
calc_rolling_means <- function(traits, data, var) {
  sapply(traits$julian_date, function(jd) {
    mean(data %>% filter(julian_date < jd, julian_date >= jd - 14) %>% pull({{ var }}), na.rm = TRUE)
  })
}

# Function to calculate mid-day mean on sampling day (9:30am–2pm)
calc_midday_means <- function(traits, data, var) {
  sapply(traits$julian_date, function(jd) {
    mean(data %>%
           filter(julian_date == jd, hour >= 9.5, hour <= 14) %>%
           pull({{ var }}), na.rm = TRUE)
  })
}

# Apply to both Tair and SWC615
traits$Tair_14day_mean <- calc_rolling_means(traits, data, `Temp(55)`)
traits$SWC615_14day_mean <- calc_rolling_means(traits, data, SWC615)

traits$Tair_midday_mean <- calc_midday_means(traits, data, `Temp(55)`)
traits$SWC615_midday_mean <- calc_midday_means(traits, data, SWC615)

#Vcmax25
ggplot(traits, aes(Tair_14day_mean, Vcmax_Estimate_final, fill = Age_leaf_mid))+
  geom_point(shape = 21, size = 3)+
  scale_fill_viridis_c()

ggplot(traits, aes(Tair_midday_mean, Vcmax_Estimate_final, fill = Age_leaf_mid))+
  geom_point(shape = 21, size = 3)+
  scale_fill_viridis_c()

ggplot(traits, aes(SWC615_14day_mean, Vcmax_Estimate_final, fill = Age_leaf_mid))+
  geom_point(shape = 21, size = 3)+
  scale_fill_viridis_c()

ggplot(traits, aes(SWC615_midday_mean, Vcmax_Estimate_final, fill = Age_leaf_mid))+
  geom_point(shape = 21, size = 3)+
  scale_fill_viridis_c()

#TLP
ggplot(traits, aes(Tair_14day_mean, PTLP, fill = Age_leaf_mid))+
  geom_point(shape = 21, size = 3)+
  scale_fill_viridis_c()

ggplot(traits, aes(Tair_midday_mean, PTLP, fill = Age_leaf_mid))+
  geom_point(shape = 21, size = 3)+
  scale_fill_viridis_c()

ggplot(traits, aes(SWC615_14day_mean, PTLP, fill = Age_leaf_mid))+
  geom_point(shape = 21, size = 3)+
  scale_fill_viridis_c()+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme(legend.position = "top")

ggplot(traits, aes(SWC615_midday_mean, PTLP, fill = Age_leaf_mid))+
  geom_point(shape = 21, size = 3)+
  scale_fill_viridis_c()+
  facet_wrap(~sp, scales = "free_y", ncol = 4, nrow = 2)+
  theme(legend.position = "top")

ggplot(data = subset(traits, !is.na(Age_leaf_mid) & !is.na(Vcmax_Estimate_final)), aes(Tair_14day_mean, Age_leaf_mid, fill = Vcmax_Estimate_final))+
  geom_point(shape = 21, size = 3, alpha = 0.6)+
  scale_fill_viridis_c()+
  theme(legend.position = "top")


#### curious about the timing of leaf flushing and senescence in the dataset.
pheno <- read_excel("Dataset/French_Guiana/20250402_Dataset_Phenoflux.xlsx", 
                    sheet = "Fluch Phénoflux (2)")

pheno <- pheno[1:24,]

# Convert to long format
pheno_long <- pheno %>%
  pivot_longer(
    cols = `23.10.2020`:`20.04.2023`,  # Replace with the actual range of your date columns
    names_to = "Date",
    values_to = "Phenology"
  ) %>%
  mutate(Date = as.Date(Date, format = "%d.%m.%Y"))  # Convert date column to proper Date format

pheno_long <- pheno_long %>%
  mutate(
    Year     = year(Date),
    Month    = month(Date, label = TRUE, abbr = TRUE),  # Use `label = FALSE` if you want numeric
    Julian   = yday(Date),
    pheno.L  = if_else(str_detect(Phenology, "\\bL;"), 1, 0),
    pheno.D  = if_else(str_detect(Phenology, "\\bD;"), 1, 0),
    pheno.F  = if_else(str_detect(Phenology, "\\bF;"), 1, 0),    # Flushing
    pheno.Fl = if_else(str_detect(Phenology, "\\bFl;"), 1, 0)    # Flowering
  )

# show proportion of trees flushing across all trees
# 1. Aggregate to get daily flushing probability
pheno_summary <- pheno_long %>%
  group_by(Year, Month, Julian) %>%
  summarise(
    n_obs     = sum(!is.na(pheno.F) | !is.na(pheno.D) | !is.na(pheno.Fl)),
    p_flush   = sum(pheno.F, na.rm = TRUE) / sum(!is.na(pheno.F)),
    p_senesce = sum(pheno.D, na.rm = TRUE) / sum(!is.na(pheno.D))
    #,p_flower  = sum(pheno.Fl, na.rm = TRUE) / sum(!is.na(pheno.Fl))
  )

pheno_long_summary <- pheno_summary %>%
  pivot_longer(
    cols = c(p_flush, p_senesce),
    names_to = "Phenophase",
    values_to = "Proportion"
  )

ggplot(pheno_long_summary, aes(x = Julian, y = Proportion, colour = Phenophase)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10), se = TRUE) +
  labs(
    x = "Julian Day",
    y = "Proportion of Trees",
    colour = "Phenophase"
  )+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank())

ggsave("./Outputs/Plots/Meteo_Fluxes/Phenology_doy.png", type = "Cairo",
       bg = "white", dpi = 300, width = 7, height = 7)


pheno_long_summary$Year <- as.factor(pheno_long_summary$Year)

ggplot(data = subset(pheno_long_summary, Phenophase == "p_flush"), aes(x = Julian, y = Proportion, colour = Year)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cc"),
              aes(group = Year))+
  labs(
    x = "Julian Day",
    y = "Proportion of Trees",
    colour = "Phenophase",
    title = "Seasonal Phenology Patterns"
  )+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank())+
  scale_color_viridis_d(option = "A")

# NOW FOR EACH SPECIES
pheno_specsummary <- pheno_long %>%
  group_by(Year, Month, Julian, Genus, Species) %>%
  summarise(
    n_obs     = sum(!is.na(pheno.F) | !is.na(pheno.D) | !is.na(pheno.Fl)),
    p_flush   = sum(pheno.F, na.rm = TRUE) / sum(!is.na(pheno.F)),
    p_senesce = sum(pheno.D, na.rm = TRUE) / sum(!is.na(pheno.D))
    #,p_flower  = sum(pheno.Fl, na.rm = TRUE) / sum(!is.na(pheno.Fl))
  )

pheno_long_specsummary <- pheno_specsummary %>%
  pivot_longer(
    cols = c(p_flush, p_senesce),
    names_to = "Phenophase",
    values_to = "Proportion"
  )

ggplot(data = subset(pheno_long_specsummary, Phenophase == "p_flush"), aes(x = Julian, y = Proportion, colour = Genus)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cc"),
              aes(group = Genus))+
  labs(
    x = "Julian Day",
    y = "Proportion of Trees",
    colour = "Phenophase",
    title = "Seasonal Phenology Patterns"
  )+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank())+
  scale_color_viridis_d(option = "A")

ggsave("./Outputs/Plots/Meteo_Fluxes/Phenology_species_doy.png", type = "Cairo",
       bg = "white", dpi = 300, width = 7, height = 7)

