# Reprocess A-Ci curves my own way to see how it compares
# I noticed that Vcmax and Jmax determined from Julien Lamour's processing vs. the original shiny app gives very different results.
# I am reprocessing them how I would to see which it closely resembles.

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# Desired column names
desired_cols <- c("obs", "A", "Ci", "CO2_s", "CO2_r", "gsw", "Pa", "Qin", "RHcham", "Tleaf", "Tleaf2", "file")

# List all Excel files
file_list <- list.files(path = "./Dataset/French_Guiana/LI6800_A_Ci/", 
                        pattern = "\\.xlsx$", full.names = TRUE)

# Process each file
combined_df <- file_list %>%
  set_names(nm = basename(.)) %>%
  map_df(~ {
    df <- read_excel(.x, skip = 13) %>%
      select(any_of(desired_cols)) %>%
      slice(-1) %>%  # remove first row
      mutate(across(everything(), as.numeric)) %>%
      mutate(filename = basename(.x))
    return(df)
  })

#tleaf not always plugged into the same port
combined_df[combined_df$Tleaf > 100, "Tleaf"] = combined_df[combined_df$Tleaf > 100, "Tleaf2"]
combined_df$Tleaf2 <- NULL

#now just checking range of data generally
hist(combined_df$Tleaf)
hist(combined_df$Ci)
hist(combined_df$A)
hist(combined_df$gsw)

#get SampleID
combined_df$SampleID=sub("\\.xlsx$", "", basename(combined_df$filename))
combined_df$SampleID=sub("ACi_", "", combined_df$SampleID)

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
combined_df$SampleID_num=as.numeric(as.factor(combined_df$SampleID))

#now extract info leaf_code so I have unique trees and cohorts and leaves
combined_df <- combined_df %>%
  separate(SampleID, into = c("leaf_code","date"), sep = "\\_", remove = FALSE)

#remove the xlsx

#now extract info leaf_code so I have unique trees and cohorts and leaves
combined_df <- combined_df %>%
  separate(leaf_code, into = c("tree", "cohort", "leaf"), sep = "\\.", remove = FALSE)

# quick check to see what gs is doing during these curves.
ggplot(combined_df, aes(CO2_s, gsw))+
  geom_point(alpha=0.4)+
  facet_wrap(~tree, scales = "free")

#before curating for aci curves - let us get the values at 415 ppm hoping they are similar to that measured (but not logged) under steady state
unique(combined_df$SampleID)

spot_df <- combined_df %>%
  subset(CO2_r >= 410 & CO2_r <= 420) %>%
  subset(Ci >= 150) %>%
  group_by(across(where(is.character))) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

hist(spot_df$A)
hist(spot_df$gsw)
hist(spot_df$A/spot_df$gsw)
hist(spot_df$Ci)

write.csv(spot_df, "./Outputs/Gasx_spot.csv", row.names = FALSE)

#a little piece from juliens processing
curated_data <- combined_df

# Adding a Quality Check column with two values: "ok" or "bad", where "bad" correspond to bad data 
# that is not used to fit the ACi curves in the next step.
curated_data$QC='ok'

# Displaying Ci values to spot bad data
hist(curated_data$Ci)
hist(curated_data$gsw)
hist(curated_data$Tleaf)

# Flagging the below 0 Ci as bad data
curated_data[curated_data$Ci<0,'QC']='bad'
curated_data[curated_data$Tleaf>100,'QC']='bad'
# I create a "Quality check" table where I manually inform bad points. 
# When I start the quality check the QC table is empty: QC_table=cbind.data.frame(SampleID_num=c(),Record=c()) 
# I plot the data in the next step (line starting with pdf below). 
# I then manually check the plots and write down the SampleID_num and the Record corresponding to bad points
# This works well if there are few points to remove.
# You can use another method if you need.

QC_table=cbind.data.frame(filename=c(55),
                          Record=c(219)) 
# Here I flag the bad curves by writing down the filename of the bad curves.
ls_bad_curve=c(20,53,68,142)

# Here I flag all the bad curves and bad points
curated_data[paste(curated_data$filename,curated_data$obs)%in%paste(QC_table$filename,QC_table$obs),'QC']='bad'
curated_data[curated_data$filename%in%ls_bad_curve,'QC']='bad'

# For this dataset, I also remove the Ci part above 500 ppm of many curves for which there is an oscillation at high CO2.
# The reason for the oscillation is probably a too fast ramp of CO2
ls_oscillate=c(15,17,18,21,27,29,31,37,59,61,63,67,83,87,89,91,99,110,117,122,142,143)

curated_data[curated_data$filename%in%ls_oscillate&curated_data$Ci>500,'QC']='bad'

# I also remove data with Ci > 1500 which often have weird behaviour
#curated_data[curated_data$Ci>1500,'QC']='bad' #note this is what Lamour does - but I am taking it down to 1000 for the weird behaviour
curated_data[curated_data$Ci>1000,'QC']='bad' #note this is what Lamour does - but I am taking it down to 1000 for the weird behaviour


# PLOT THIS
QC_plot <- ggplot(curated_data, aes(Ci, A, col = QC))+
  geom_point()+
  facet_wrap(~SampleID_num)+
  theme_bw()

# Keeping only the good points and curves
curated_data=curated_data[curated_data$QC=='ok',]

ggplot(curated_data, aes(Ci, A, col = QC))+
  geom_point()+
  facet_wrap(~SampleID_num)+
  theme_bw()

ggplot(curated_data, aes(Ci, A, col = factor(date)))+
  geom_point()+
  facet_wrap(~tree)+
  theme_bw()

################################################################################

# Now fit A-Ci curves using plantecophys package

library(plantecophys)

curated_data <- data.frame(curated_data)

#first fit without TPU much faster - fit with bilinear method so it can be comparable to the TPU fits
aci <- fitacis(data = curated_data,  group = ("SampleID"), id = c("leaf_code", "tree", "cohort", "leaf", "date", "SampleID_num"),
        varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
        Tcorrect = TRUE,
        Patm = 101,
        fitTPU = FALSE,
        fitmethod = "bilinear"
       )
#aci.res <- coef(aci) %>% data.frame() #usually i just do this but doesn't work if some of the curves are fit with different method

plot(aci[[1]])

for (i in seq_along(aci)) {
  sample_id <- names(aci)[i]                      # Name in list (e.g., "1.1.1_20220712")
  sample_num <- aci[[i]]$df$SampleID_num[1]       # Numeric ID from the data frame
  
  # Filename includes both identifiers
  filename <- paste0("./Outputs/Plots/Inspect ACi curves/Plantecophys-notpu/aci_plot_", sample_num, "_", sample_id, ".pdf")
  
  pdf(filename, width = 6, height = 5)
  plot(aci[[i]],
       main = paste("SampleID:", sample_id),
       sub = paste("SampleID_num:", sample_num, "   RMSE:", round(aci[[i]]$RMSE, 2)))
  dev.off()
}

#now look at saved curves and identify things that don't have the best fit.
#these notes are when i only excluded vals >1500 not >1000

# no ? - 6.1.2_20230202 - curve looks bad, maybe need to subset for Ci up to 800?
# no 2 - 1.1.1_20221008 - bad fit, but why? doesnt properly saturated maybe that is why
# no 3 - 1.1.1_20221128 - bad fit, TPU affecting from 1000 on?
# no 6 - 1.1.2_20221008 - high RMSE due to noisy data at low ci but fit itself looks fine
# no 9 - 10.1.1_20221011 - highish RMSE, fit not the worst but looks a litte off to me.
# no 11 - 10.2.1_20221011 - noisy but fit looks ok
# no 12 - 12.1.1_20220630 - fit looks a bit wrong. maybe need to remove points after 1100
# no 14 - 12.1.1_20220725 - noisy but fit looks good.
# no 15 - 12.1.1_20230306 - bad fit due to TPU limitation - remove points above 1000
# no 17 - 12.1.2_20221205 - bad fit due to TPU limitation - try removing points above 1000?
# no 18 - 12.1.2_20230306 - bad fit due to TPU limitation - try removing points above 1000?
# no 19 - 12.1.3_20220630 - bad fit, random increase from 800-900. Remove poitns above 800?? maybe it wont affect it
# no 20 - 12.1.3_20220725 - terrible curve. remove?? fit might be fine but its honestly not great
# no 21 - 21.1.3_202221205 - bad fit due to TPU limitation or perhaps an overshoot? maybe remove points above 600?
# no 22 - 12.1.3_20230306 - bad fit due to TPU limitation - try removing points above 1000?
# no 23 - 12.2.1_20221205 - bad fit due to TPU limitation - but dont think removing points will help
# no 24 - 12.2.2_20221205 - bad fit due to TPU limitation - try removing vals above 1000?
# no 27 - 13.1.3_20221018 - bad fit due to TPU limitation - try removing vals above 1000? - perhaps the point isnt actually terrible
# no 29 - 14.1.1_20220720 - bad fit due to TPU limitation or overshoot - try removing vals above 1000
# no 31 - 14.1.1_20230209 - bad fit due to TPU limitation - try removing vals above 1000
# no 32 - 14.1.2_20220720 - bad fit due to TPU limitation - try removing vals above 1000
# no 34 - 14.1.2_20230209 - bad fit due to TPU limitation - try removing vals above 1000
# no 35 - 14.1.3_20220720 - TPU limitation above 800 but fit not terrible
# no 37 - 14.1.3_20230209 - bad fit due to TPU limitation or an overshoot - try removing vals above 1000 but this one probably wont help
# no 42 - 15.1.1_20230209 - bad fit due to TPU limitation - try removing vals above 1000
# no 43 - 15.1.2_20220720 - bad fit due to TPU limitation - try removing vals above 1000
# no 53 - 16.1.2_20230124 - bad fit due to noisy curve, this should be removed super underestimated.
# no 55 - 16.1.3_20230124 - weird curve, fit is okish - does the overshoot hing and doesnt go up to a high Ci. 
# no 56 - 16.2.1_20221014 - probably will improve after removing >1000 vals
# no 58 - 17.1.1_20230124 - bad fit due to TPU limitation and overshoot - probably wont improve
# no 59 - 17.1.2_20220712 - bad fit due to TPU and maybe overshoot - removing als >900 might help?
# no 61 - 17.1.2_20230124 - bad fit due to TPU and overshoot but actually the estimate may be correct
# no 63 - 17.1.3_20230124 - bad fit due to TPU and overshoot - only goes up to 800 so wont improve with subsetting
# no 65 - 18.1.2_20220707 - not the best fit due to TPU - removing about 1000 may help but actually estimate may be ok
# no 67 - 18.1.3_20230124 - terrible fit - overshoot and/or TPU? An is also generally low. Removing vals above 1000 might help
# no 68 - 19.1.1_20221006 - bad fit - need to remove vals above 1000 due to weird increase - wont be fixed by fitting TPU
# no 69 - 19.1.1_20221130 - bad fit but unclear why. maybe TPU in this case would fix the estimate. weird cyclic thing at bottom
# no 72 - 19.1.2_20230216 - bad fit due to TPU - dont think subsetting will help
# no 74 - 19.2.2_20221130 - bad fit due to TPU- dont think subsetting will help but can try.
# no 78 - 2.1.2_20220725 - bad fit - can try subsetting but might be due to the bottom part of the curve. 
# no 83 - 20.1.3_20221010 - bad fit due to TPU and overshoot. TPU would help subsetting probably wont
# no 86 - 21.1.1_20221010 - bad fit, but not sure what would help
# no 87 - 21.1.3_20230117 - very bad fit, TPU and overshoot, dont think subsetting will help here.
# no 88 - 21.2.1_20221010 - bad fit, overshoot and tpu. probably needs TPU to be better.
# no 89 - 21.2.1_20230117 - bad fit but estimate may be ok
# no 97 - 22.1.3_20221012 - bad fit TPU and overshoot, subsetting > 1000 may help
# no 99 - 22.2.1_20221205 - bad fit, TPU and overshoot, subsetting probably wont help here
# no 106 - 23.1.1_20221011 - bad fit, TPU and overshoot and outliers, needs TPU fitting here to be correct
# no 107 - 23.1.1_20221124 - bad fit, needs TPU
# no 110 - 23.1.2_20221124 - bad fit, ovrshoot and TPU, remove vals after 800 may help
# no 111 - 23.1.3_20221124 - bad fit, needs TPU or removal of vals > 800
# no 117 - 3.1.1_20221012 - terrible fit, TPU and overshoot. Subsetting wont help
# no 118 - 3.1.1_20221130 - bad fit, needs TPU
# no 119 - 3.1.1_20230308 - bad fit, TPU afffected but subsetting > 800 may help
# no 122 - 3.1.2_20221130 - bad fit - overshoot and TPU - subsetting wont help.
# no 125 - 4.1.2_20220707 - bad fit - TPU, subsetting > 1000 may help
#no 126 - same as above
#etc etc

# basically - we need to both subset and fit TPU. 

####################
library(tibble)

aci.res <- map_dfr(aci, function(fit) {
  rmse <- if (!is.null(fit$RMSE)) fit$RMSE else NA_real_
  
  if (!is.null(fit$pars) && nrow(fit$pars) == 3) {
    pars <- as.data.frame(fit$pars) %>%
      rownames_to_column("parameter") %>%
      pivot_wider(names_from = parameter,
                  values_from = c(Estimate, `Std. Error`),
                  names_glue = "{parameter}_{.value}") %>%
      mutate(RMSE = rmse)
  } else {
    pars <- tibble(
      Vcmax_Estimate = NA, Vcmax_Std.Error = NA,
      Jmax_Estimate = NA, Jmax_Std.Error = NA,
      Rd_Estimate = NA, Rd_Std.Error = NA,
      RMSE = rmse
    )
  }
})

# add id columns back in.
aci.res$SampleID <- unique(curated_data$SampleID)

hist(aci.res$Vcmax_Estimate) #seems to be at least one high Vcmax outlier
hist(aci.res$Jmax_Estimate)

#checking plot of Vcmax vs. Jmax to see if it looks appropriate
ggplot(aci.res, aes(Vcmax_Estimate, Jmax_Estimate))+
  geom_point(size = 3, aes(col = `Vcmax_Std. Error`))+
  scale_color_viridis_c(limits = c(0,6))+
  theme_minimal(base_size = 15)+
  theme(legend.position = "top")
#can exclude outliers wit SE > 6

ggplot(aci.res, aes(Vcmax_Estimate, Jmax_Estimate))+
  geom_point(size = 3, aes(col = `Jmax_Std. Error`))+
  scale_color_viridis_c(limits = c(0,6))+
  theme_minimal(base_size = 15)

#get info from sample ID
aci.res <- aci.res %>%
  separate(SampleID, into = c("leaf_code", "obs_date"), sep = "_", remove = FALSE)

# extract date from unique code so that I have it in correct format since other column is being annoying
aci.res <- aci.res %>%
  mutate(date = as.POSIXct(obs_date, format = "%Y%m%d"))

#now extract info leaf_code so I have unique trees and cohorts and leaves
aci.res <- aci.res %>%
  separate(leaf_code, into = c("tree", "cohort", "leaf"), sep = "\\.", remove = FALSE)

#now make a quality check based on the SE or RMSE
aci.res$Vcmax_QC <- "good"
aci.res$Vcmax_QC[aci.res$`Vcmax_Std. Error` >= 6] <- "bad"

aci.res$Jmax_QC <- "good"
aci.res$Jmax_QC[aci.res$`Jmax_Std. Error` >= 6] <- "bad"

#reorder columns
aci.res <- aci.res %>% 
  select(SampleID, leaf_code, tree, cohort, leaf, obs_date, date, everything())

#save this data
write.csv(aci.res, "./Outputs/Aci_NOTPU_plantecophys.csv", row.names = FALSE)

#plot
ggplot(aci.res, aes(obs_date, Vcmax_Estimate))+
  geom_point(aes(col = Vcmax_QC))+
  facet_wrap(~tree)

ggplot(aci.res, aes(obs_date, Jmax_Estimate))+
  geom_point(aes(col = Jmax_QC))+
  facet_wrap(~tree)

ggplot(data = subset(aci.res, Vcmax_QC != "bad" & Jmax_QC != "bad"), aes(obs_date, Vcmax_Estimate/Jmax_Estimate))+
  geom_point()+
  facet_wrap(~tree)

ggplot(data = subset(aci.res, Vcmax_QC != "bad" & Jmax_QC != "bad"), aes(obs_date, Vcmax_Estimate))+
  geom_point()+
  facet_wrap(~tree)

################################################################################

################################################################################
#now fit with TPU to compare - it takes AGES. Licor manual suggests fitting to a smaller dataset e.g. up to 1000 Ci
#however the TPU limitation sometimes occurs at much lower Ci values I notice.

# fit a subset for now to see if the values look closer to julien lamours

curated_data <- curated_data %>%
  mutate(date = as.POSIXct(date, format = "%Y%m%d"))

curated_data$month <- lubridate::month(curated_data$date, label = TRUE, abbr = TRUE)

test <- subset(curated_data, month == "Oct")

test.aci.tpu <- fitacis(data = test,  group = ("SampleID"), id = c("leaf_code", "tree", "cohort", "leaf", "date", "month", "SampleID_num"),
                   varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
                   Tcorrect = TRUE,
                   Patm = 101,
                   fitTPU = TRUE,
                   quiet = FALSE
)


test.aci.tpu.res <- map_dfr(test.aci.tpu, function(fit) {
  rmse <- if (!is.null(fit$RMSE)) fit$RMSE else NA_real_
  
  if (!is.null(fit$pars) && nrow(fit$pars) == 4) {
    pars <- as.data.frame(fit$pars) %>%
      rownames_to_column("parameter") %>%
      pivot_wider(names_from = parameter,
                  values_from = c(Estimate, `Std. Error`),
                  names_glue = "{parameter}_{.value}") %>%
      mutate(RMSE = rmse)
  } else {
    pars <- tibble(
      Vcmax_Estimate = NA, Vcmax_Std.Error = NA,
      Jmax_Estimate = NA, Jmax_Std.Error = NA,
      Rd_Estimate = NA, Rd_Std.Error = NA,
      TPU_Estimate = NA, TPU_Std.Error = NA,
      RMSE = rmse
    )
  }
})

# add id columns back in.
test.aci.tpu.res$SampleID <- unique(test$SampleID)

hist(test.aci.tpu.res$Vcmax_Estimate)
hist(test.aci.tpu.res$Jmax_Estimate)

summary(test.aci.tpu.res$Vcmax_Estimate)
summary(test.aci.tpu.res$Jmax_Estimate) #one jmax bad

#checking plot of Vcmax vs. Jmax to see if it looks appropriate
ggplot(test.aci.tpu.res, aes(Vcmax_Estimate, Jmax_Estimate))+
  geom_point(size = 3, aes(col = `Vcmax_Std. Error`))+
  scale_color_viridis_c(limits = c(0,6))+
  theme_minimal(base_size = 15)

#get info from sample ID
test.aci.tpu.res <- test.aci.tpu.res %>%
  separate(SampleID, into = c("leaf_code", "obs_date"), sep = "_", remove = FALSE)

# extract date from unique code so that I have it in correct format since other column is being annoying
test.aci.tpu.res <- test.aci.tpu.res %>%
  mutate(date = as.POSIXct(obs_date, format = "%Y%m%d"))

#now extract info leaf_code so I have unique trees and cohorts and leaves
test.aci.tpu.res <- test.aci.tpu.res %>%
  separate(leaf_code, into = c("tree", "cohort", "leaf"), sep = "\\.", remove = FALSE)

#now make a quality check based on the SE or RMSE
test.aci.tpu.res$Vcmax_QC <- "good"
#test.aci.tpu.res$Vcmax_QC[test.aci.tpu.res$`Vcmax_Std. Error` >= 6] <- "bad"

test.aci.tpu.res$Jmax_QC <- "good"
test.aci.tpu.res$Jmax_QC[test.aci.tpu.res$Jmax_Estimate >= 500] <- "bad"
#test.aci.tpu.res$Jmax_QC[test.aci.tpu.res$`Jmax_Std. Error` >= 6] <- "bad"

ggplot(data = subset(test.aci.tpu.res, Jmax_QC == "good"), aes(Vcmax_Estimate, Jmax_Estimate))+
  geom_point(size = 3, aes(col = `Vcmax_Std. Error`))+
  scale_color_viridis_c(limits = c(0,6))+
  theme_minimal(base_size = 15)

#reorder columns
test.aci.tpu.res <- test.aci.tpu.res %>% 
  select(SampleID, leaf_code, tree, cohort, leaf, obs_date, date, everything())

#now combine fits with and without tpu
# Define the key columns for merging
key_cols <- c("leaf_code", "tree", "cohort", "leaf", "date")

# For the first dataframe, add suffix "_tpufit" to all non-key columns
tpufit_renamed <- test.aci.tpu.res %>%
  rename_with(~ paste0(., "_tpufit"), .cols = setdiff(names(.), key_cols))

# For the second dataframe, add suffix "_notpufit" to all non-key columns
notpufit_renamed <- aci.res %>%
  rename_with(~ paste0(., "_notpufit"), .cols = setdiff(names(.), key_cols))

# Merge the two dataframes by the key columns
merged_df <- full_join(tpufit_renamed, notpufit_renamed, by = key_cols)


ggplot(data = subset(merged_df, Jmax_QC_tpufit == "good"), 
       aes(Vcmax_Estimate_tpufit, Vcmax_Estimate_notpufit, fill = `Vcmax_Std. Error_notpufit`))+
  geom_point(size = 3, shape = 21)+
  geom_abline()+
  scale_fill_viridis_c()

ggplot(data = subset(merged_df, Jmax_QC_tpufit == "good"), 
       aes(Jmax_Estimate_tpufit, Jmax_Estimate_notpufit, fill = `Jmax_Std. Error_tpufit`))+
  geom_point(size = 3, shape = 21)+
  geom_abline()+
  scale_fill_viridis_c()

#ok data become basically the same as lamours when fitting TPU.

#################################################################################
#create 5 subsets of the data since this takes so long.

# Unique SampleIDs
sample_ids <- unique(curated_data$SampleID)

# Split into 5 roughly equal groups
sample_chunks <- split(sample_ids, cut(seq_along(sample_ids), 5, labels = FALSE))

# Create 5 data frame subsets
curated_data_sub1 <- curated_data %>% filter(SampleID %in% sample_chunks[[1]])
curated_data_sub2 <- curated_data %>% filter(SampleID %in% sample_chunks[[2]])
curated_data_sub3 <- curated_data %>% filter(SampleID %in% sample_chunks[[3]])
curated_data_sub4 <- curated_data %>% filter(SampleID %in% sample_chunks[[4]])
curated_data_sub5 <- curated_data %>% filter(SampleID %in% sample_chunks[[5]])

#now run on each chunk one by one
aci_sub1 <- fitacis(data = curated_data_sub1, group = "SampleID",
                    id = c("leaf_code", "tree", "cohort", "leaf", "date", "SampleID_num"),
                    varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
                    Tcorrect = TRUE,
                    Patm = 101,
                    fitTPU = TRUE,
                    quiet = FALSE)
beepr::beep()

aci_sub2 <- fitacis(data = curated_data_sub2, group = "SampleID",
                    id = c("leaf_code", "tree", "cohort", "leaf", "date", "SampleID_num"),
                    varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
                    Tcorrect = TRUE,
                    Patm = 101,
                    fitTPU = TRUE,
                    quiet = FALSE)
beepr::beep()

aci_sub3 <- fitacis(data = curated_data_sub3, group = "SampleID",
                    id = c("leaf_code", "tree", "cohort", "leaf", "date", "SampleID_num"),
                    varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
                    Tcorrect = TRUE,
                    Patm = 101,
                    fitTPU = TRUE,
                    quiet = FALSE)
beepr::beep()

aci_sub4 <- fitacis(data = curated_data_sub4, group = "SampleID",
                    id = c("leaf_code", "tree", "cohort", "leaf", "date", "SampleID_num"),
                    varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
                    Tcorrect = TRUE,
                    Patm = 101,
                    fitTPU = TRUE,
                    quiet = FALSE)
beepr::beep()

aci_sub5 <- fitacis(data = curated_data_sub5, group = "SampleID",
                    id = c("leaf_code", "tree", "cohort", "leaf", "date", "SampleID_num"),
                    varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
                    Tcorrect = TRUE,
                    Patm = 101,
                    fitTPU = TRUE,
                    quiet = FALSE)
beepr::beep()

#### NOW COMBINE THE DATA AND SAVE SO I NEVER HAVE TO DO THIS AGAIN
#aci.tpu.res <- coef(aci.tpu) %>% data.frame()


##### FUNCTION TO SAVE ALL THE PLOTS FROM THE TPU FITTING

# Function to plot and save the ACi curves for a given dataset
plot_save_aci <- function(aci_data, output_dir = "./Outputs/Plots/Inspect ACi curves/Plantecophys-tpufit/") {
  for (i in seq_along(aci_data)) {
    sample_id <- names(aci_data)[i]  # Name in list (e.g., "1.1.1_20220712")
    sample_num <- aci_data[[i]]$df$SampleID_num[1]  # Numeric ID from the data frame
    
    # Filename includes both identifiers
    filename <- paste0(output_dir, "aci_plot_", sample_num, "_", sample_id, ".pdf")
    
    # Create the output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Plot and save the figure
    pdf(filename, width = 6, height = 5)
    plot(aci_data[[i]],
         main = paste("SampleID:", sample_id),
         sub = paste("SampleID_num:", sample_num, "   RMSE:", round(aci_data[[i]]$RMSE, 2)))
    dev.off()
  }
}

# Plot and save results for subsets
plot_save_aci(aci_sub1)
plot_save_aci(aci_sub2)
plot_save_aci(aci_sub3)
plot_save_aci(aci_sub4)
plot_save_aci(aci_sub5)

##### NOW SAVE THE FIT RESULTS IN A DATAFRAME AND SAVE THE FILE
library(tibble)

# Function to extract parameters and RMSE from a single fit
extract_aci_results <- function(fit) {
  rmse <- if (!is.null(fit$RMSE)) fit$RMSE else NA_real_
  
  # Extract parameters and standard errors
  pars <- if (!is.null(fit$pars)) {
    as.data.frame(fit$pars) %>%
      rownames_to_column("parameter") %>%
      pivot_wider(names_from = parameter,
                  values_from = c(Estimate, `Std. Error`),
                  names_glue = "{parameter}_{.value}")
  } else {
    tibble(
      Vcmax_Estimate = NA, Vcmax_Std.Error = NA,
      Jmax_Estimate = NA, Jmax_Std.Error = NA,
      Rd_Estimate = NA, Rd_Std.Error = NA,
      TPU_Estimate = NA, TPU_Std.Error = NA
    )
  }
  
  # Add RMSE to the results
  pars <- pars %>% mutate(RMSE = rmse)
  
  return(pars)
}

# Combine results for each subset (aci_sub1 to aci_sub5)
aci_tpu_res <- map_dfr(list(aci_sub1, aci_sub2, aci_sub3, aci_sub4, aci_sub5), 
                       function(subset_aci) {
                         # Extract results from each curve in the subset
                         map_dfr(subset_aci, extract_aci_results)
                       })

# Add SampleID to the results (ensure each row corresponds to a unique SampleID)
aci_tpu_res$SampleID <- rep(unique(curated_data$SampleID), each = nrow(aci_tpu_res) / length(unique(curated_data$SampleID)))

# View the final result
# You can write the results to a CSV file or examine it directly
write.csv(aci_tpu_res, "Outputs/Aci_fits_TPU_plantecophys.csv", row.names = FALSE)

ggplot(aci_tpu_res, aes(Vcmax_Estimate, Jmax_Estimate))+
  geom_point()

hist(aci_tpu_res$Vcmax_Estimate)
hist(aci_tpu_res$Jmax_Estimate)
summary(aci_tpu_res$Jmax_Estimate)

hist(aci_tpu_res$TPU_Estimate)

ggplot(data = subset(aci_tpu_res, Jmax_Estimate < 400), aes(Vcmax_Estimate, Jmax_Estimate))+
  geom_point()

ggplot(aci_tpu_res, aes(Vcmax_Estimate, TPU_Estimate))+
  geom_point()

#### NOW COMPARE THE MODELS WITH WITH VS WITHOUT TPU LIMITATION
names(aci.res)
names(aci_tpu_res)

library(dplyr)

# Rename columns in aci_tpu_res to include _tpu suffix
aci_tpu_res_renamed <- aci_tpu_res %>%
  rename(
    Vcmax_Estimate_tpu = Vcmax_Estimate,
    Jmax_Estimate_tpu = Jmax_Estimate,
    Rd_Estimate_tpu = Rd_Estimate,
    TPU_Estimate = TPU_Estimate,
    Vcmax_Std.Error_tpu = `Vcmax_Std. Error`,
    Jmax_Std.Error_tpu = `Jmax_Std. Error`,
    Rd_Std.Error_tpu = `Rd_Std. Error`,
    TPU_Std.Error = `TPU_Std. Error`,
    RMSE_tpu = RMSE
  )

# Merge aci.res and aci_tpu_res_renamed by SampleID
merged_data <- left_join(aci.res, aci_tpu_res_renamed, by = "SampleID")

# View the merged dataframe
head(merged_data)

##### NOW PLOT TO SEE WHICH FITS ARE BETTER WITH VS WITHOUT TPU
ggplot(merged_data, aes(RMSE, RMSE_tpu))+
  geom_point()+
  geom_abline()

ggplot(merged_data, aes(`Vcmax_Std. Error`, Vcmax_Std.Error_tpu))+
  geom_point()+
  geom_abline()

#now let us make a column that indicates which model is supposedly better
merged_data$TPU_improved_vcmax <- 0
merged_data$TPU_improved_rmse <- 0
merged_data$TPU_improved_vcmax[merged_data$Vcmax_Std.Error_tpu < merged_data$`Vcmax_Std. Error`] <- 1
merged_data$TPU_improved_rmse[merged_data$RMSE_tpu < merged_data$RMSE] <- 1
merged_data$TPU_improved_both <- merged_data$TPU_improved_rmse + merged_data$TPU_improved_vcmax

#now when did it make weird jmax estimates???
merged_data$TPU_created_spuriousjmax <- 0
merged_data$TPU_created_spuriousjmax[merged_data$Jmax_Estimate_tpu > 500] <- 1

merged_data$useTPU <- merged_data$TPU_improved_both
merged_data$useTPU[merged_data$TPU_created_spuriousjmax > 0] <- 0

length(which(merged_data$TPU_improved_both > 0))
which(merged_data$TPU_improved_both > 0)

length(which(merged_data$useTPU > 0))
which(merged_data$useTPU > 0)

# Optionally, write to CSV
write.csv(merged_data, "Outputs/merged_aci_data.csv", row.names = FALSE)

#### NOW MAKE A DATAFRAME THAT GIVES THE SELECTED PARAMETERS AND INDICATES WHICH FIT THEY CAME FROM

# Create a new dataframe where parameters default to the non-TPU values unless useTPU > 0
merged_with_tpu <- merged_data %>%
  mutate(
    Vcmax_Estimate_final = ifelse(useTPU > 0, Vcmax_Estimate_tpu, Vcmax_Estimate),
    Jmax_Estimate_final = ifelse(useTPU > 0, Jmax_Estimate_tpu, Jmax_Estimate),
    Rd_Estimate_final = ifelse(useTPU > 0, Rd_Estimate_tpu, Rd_Estimate),
    TPU_Estimate_final = ifelse(useTPU > 0, TPU_Estimate, NA_real_),
    Vcmax_Std.Error_final = ifelse(useTPU > 0, Vcmax_Std.Error_tpu, `Vcmax_Std. Error`),
    Jmax_Std.Error_final = ifelse(useTPU > 0, Jmax_Std.Error_tpu, `Jmax_Std. Error`),
    Rd_Std.Error_final = ifelse(useTPU > 0, Rd_Std.Error_tpu, `Rd_Std. Error`),
    TPU_Std.Error_final = ifelse(useTPU > 0, TPU_Std.Error, NA_real_),
    RMSE_final = ifelse(useTPU > 0, RMSE_tpu, RMSE)
  )

# View the resulting dataframe
head(merged_with_tpu)

ggplot(merged_with_tpu, aes(Vcmax_Estimate_final, Jmax_Estimate_final))+
  geom_point()

ggplot(merged_with_tpu, aes(obs_date, Vcmax_Estimate_final))+
  geom_point()

ggplot(merged_with_tpu, aes(obs_date, Jmax_Estimate_final/Vcmax_Estimate_final))+
  geom_point()

ggplot(merged_with_tpu, aes(obs_date, TPU_Estimate_final/Vcmax_Estimate_final))+
  geom_point()

# so much better!

# Optionally, write the resulting dataframe to a CSV
write.csv(merged_with_tpu, "Outputs/final_Aci_data.csv", row.names = FALSE)


