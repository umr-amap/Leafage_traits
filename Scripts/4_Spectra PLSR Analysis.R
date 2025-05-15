# PLSR Analysis to predict Vcmax from Spectra.

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)

##### READ IN REFLECTANCE
data <- read_table("Outputs/Processed reflectance data/Reflectance_Paracou_Processed.csv")
data <- data.frame("Wavelength" = data$`"wavelength"`,"Ref_smooth" =data$`"smooth_value_sg"`,
                   "file_name" = data$`"file_name"`,
                   "leaf_code" = data$`"leaf_code"`, "obs_date" = data$`"obs_date"`)

data <- data %>% separate(leaf_code, into = c("tree", "cohort", "leaf"), sep = "\\.", remove = FALSE)

##### READ IN TRAIT DATA
traits<- read_csv("Outputs/Combined_data.csv")

#### PREP DATASETS FOR MERGING
names(data)
head(data)

head(traits$leaf_code) 
head(traits$tree)
head(traits$cohort)
head(traits$leaf)

#need to fix the spectra data id columns
data$obs_date <- gsub("[^0-9]", "", data$obs_date)
data$tree <- as.numeric(gsub('\"', '', data$tree))
data$leaf <- as.numeric(gsub('\"', '', data$leaf))
data$leaf_code <- gsub('\"', '', data$leaf_code)

data <- data %>%
  mutate(date = ymd(obs_date))

data$code_unique <- paste(data$leaf_code, data$obs_date, sep = "_")

##### CHECK AND CLEAN REFLECTANCE DATA
p <- ggplot(data = subset(data, tree == 16), aes(Wavelength, Ref_smooth, col = file_name))+
  geom_point()+
  facet_wrap(~tree)
p
plotly::ggplotly(p)

#for now remove the start and end but need to work out why that is happening and fix in original script. also remove that one weir dindividual 
data <- subset(data, Wavelength >= 500 & Wavelength <= 940)
data <- subset(data, file_name != "\"reflection2_16.1.2_20230124.txt\"")

# get average for each wavelength and code_unique

mean_ref_smooth <- data %>%
  group_by(code_unique, leaf_code, tree, cohort, leaf, obs_date, date, Wavelength) %>%
  summarise(Ref_mean = mean(Ref_smooth, na.rm = TRUE),
            Ref_sd = sd(Ref_smooth, na.rm = TRUE))

#also get average reflectance across all wavelengths
avg_ref <- data %>%
  group_by(code_unique, leaf_code, tree, cohort, leaf, obs_date, date) %>%
  summarise(Avg_Ref_tot = mean(Ref_smooth, na.rm = TRUE))

# merge
traits <- rename(traits, "code_unique" = "code_unique...1")
traits$obs_date <- NULL

#merge
merged_data <- merge(mean_ref_smooth, traits, by = c("code_unique", "leaf_code", "tree", "leaf", "date"))
merged_data <- merge(merged_data, avg_ref, by = c("code_unique", "leaf_code", "tree", "leaf", "date"))


#quick plot
ggplot(merged_data, aes(Wavelength, Ref_mean, col = Age_leaf_mid))+
  geom_point()+
  facet_wrap(~sp)+
  scale_color_viridis_c()

ggplot(merged_data, aes(Age_leaf_mid, Avg_Ref_tot))+
  geom_point()+
  facet_wrap(~sp)
#generally younger leaves have higher reflectance, than matureleaves, perhaps goes up again for old leaves??

ggplot(merged_data, aes(julian_date, Avg_Ref_tot))+
  geom_point()+
  facet_wrap(~sp)



#################################################################################
# NOW WE WANT TO SEE IF WE CAN PREDICT VCMAX USING THE SPECTRA

#follow this guide: https://academic.oup.com/jxb/article/72/18/6175/6299948#304472874

###############################################################################

outdir <- "./Outputs/PLSR_out"

### Load libraries
list.of.packages <- c("pls","dplyr","here","plotrix","ggplot2","gridExtra", "tidyverse", "spectratrait")
invisible(lapply(list.of.packages, library, character.only = TRUE))

### Setup options

# Script options
pls::pls.options(plsralg = "oscorespls")
pls::pls.options("plsralg")

# Default par options
opar <- par(no.readonly = T)

# What is the target variable?
inVar <- "Vcmax_Estimate_final"




# #just keep the unique id and the trait
# trait <- trait[,c(38, 36, 11)] #Vcmax_Estimate_final as an example
# unique(duplicated(trait$unique_id))
# 
# ### Create plsr dataset
# Start.wave <- 500
# End.wave <- 940
# wv <- seq(Start.wave,End.wave,1)
# Spectra <- as.matrix(dat_raw[,names(dat_raw) %in% wv])
# colnames(Spectra) <- c(paste0("Wave_",wv))
# head(Spectra)[1:6,1:10]


# Step 1: Select the relevant columns
plsr_data <- merged_data %>%
  dplyr::select(code_unique, Vcmax_Estimate_final, Wavelength, Ref_mean)

# Step 2: Reshape the data to ensure Wavelength is in columns
# Assuming each observation is associated with multiple wavelengths
plsr_data <- plsr_data %>%
  gather(key = "Wavelength", value = "Ref_mean", starts_with("Wave_"))





plsr_data <- plsr_data[complete.cases(plsr_data[,names(plsr_data) %in% 
                                                  c(inVar,paste0("Wave_",wv))]),]

#check if trait data is normally distributed
ggplot(plsr_data, aes(Vcmax_Estimate_final))+
  geom_histogram() #data is skewed so should be transformed

#transform data if needed
#plsr_data$cci_log <- log10(plsr_data$cci_mean)

#or exclude outlier
plsr_data <- subset(plsr_data, Vcmax_Estimate_final <= 80)

###############################################################################
### Create cal/val datasets
## Make a stratified random sampling in the strata USDA_Species_Code and Domain

method <- "dplyr" #base/dplyr
# base R - a bit slow
# dplyr - much faster
split_data <- spectratrait::create_data_split(dataset=plsr_data, approach=method, split_seed=2356812, 
                                              prop=0.8, group_variables="code_unique")
names(split_data)
cal.plsr.data <- split_data$cal_data
head(cal.plsr.data)[1:8]
val.plsr.data <- split_data$val_data
head(val.plsr.data)[1:8]
rm(split_data)

# Datasets:
print(paste("Cal observations: ",dim(cal.plsr.data)[1],sep=""))
print(paste("Val observations: ",dim(val.plsr.data)[1],sep=""))

cal_hist_plot <- qplot(cal.plsr.data[,paste0(inVar)],geom="histogram",
                       main = paste0("Cal. Histogram for ",inVar),
                       xlab = paste0(inVar),ylab = "Count",fill=I("grey50"),col=I("black"),alpha=I(.7))
val_hist_plot <- qplot(val.plsr.data[,paste0(inVar)],geom="histogram",
                       main = paste0("Val. Histogram for ",inVar),
                       xlab = paste0(inVar),ylab = "Count",fill=I("grey50"),col=I("black"),alpha=I(.7))
histograms <- grid.arrange(cal_hist_plot, val_hist_plot, ncol=2)

################################################################################
### Format PLSR data for model fitting 
cal_spec <- as.matrix(cal.plsr.data[, which(names(cal.plsr.data) %in% paste0("Wave_",wv))])
cal.plsr.data <- data.frame(cal.plsr.data[, which(names(cal.plsr.data) %notin% paste0("Wave_",wv))],
                            Spectra=I(cal_spec))
head(cal.plsr.data)

val_spec <- as.matrix(val.plsr.data[, which(names(val.plsr.data) %in% paste0("Wave_",wv))])
val.plsr.data <- data.frame(val.plsr.data[, which(names(val.plsr.data) %notin% paste0("Wave_",wv))],
                            Spectra=I(val_spec))

################################################################################
### plot cal and val spectra
par(mfrow=c(1,2)) # B, L, T, R
spectratrait::f.plot.spec(Z=cal.plsr.data$Spectra,wv=wv,plot_label="Calibration")
spectratrait::f.plot.spec(Z=val.plsr.data$Spectra,wv=wv,plot_label="Validation")


dev.copy(png,file.path(outdir,paste0(inVar,'_Cal_Val_Spectra.png')), 
         height=2500,width=4900, res=340)
dev.off();
par(mfrow=c(1,1))


################################################################################
### Use permutation to determine the optimal number of components
if(grepl("Windows", sessionInfo()$running)){
  pls.options(parallel = NULL)
} else {
  pls.options(parallel = parallel::detectCores()-1)
}

method <- "firstPlateau" #pls, firstPlateau, firstMin
random_seed <- 2356812
seg <- 40 #originally had this at 100 but need to have less segments than observations this only matters if method = "pls"
maxComps <- 18
iterations <- 1000
prop <- 0.60
if (method=="pls") {
  nComps <- spectratrait::find_optimal_components(dataset=cal.plsr.data, targetVariable=inVar,
                                                  method=method, 
                                                  maxComps=maxComps, seg=seg, 
                                                  random_seed=random_seed)
  print(paste0("*** Optimal number of components: ", nComps))
} else {
  nComps <- spectratrait::find_optimal_components(dataset=cal.plsr.data, targetVariable=inVar,
                                                  method=method, 
                                                  maxComps=maxComps, iterations=iterations, 
                                                  seg=seg, prop=prop, 
                                                  random_seed=random_seed)
}
dev.copy(png,file.path(outdir,paste0(paste0(inVar,"_PLSR_Component_Selection.png"))), 
         height=2800, width=3400,  res=340)
dev.off();


################################################################################
### Fit final model
segs <-40
plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")),scale=FALSE,ncomp=nComps,validation="CV",
                 segments=segs, segment.type="interleaved",trace=FALSE,data=cal.plsr.data)
fit <- plsr.out$fitted.values[,1,nComps]
pls.options(parallel = NULL)

# External validation fit stats
par(mfrow=c(1,2)) # B, L, T, R
pls::RMSEP(plsr.out, newdata = val.plsr.data)
plot(pls::RMSEP(plsr.out,estimate=c("test"),newdata = val.plsr.data), main="MODEL RMSEP",
     xlab="Number of Components",ylab="Model Validation RMSEP",lty=1,col="black",cex=1.5,lwd=2)
box(lwd=2.2)

pls::R2(plsr.out, newdata = val.plsr.data)
plot(pls::R2(plsr.out,estimate=c("test"),newdata = val.plsr.data), main="MODEL R2",
     xlab="Number of Components",ylab="Model Validation R2",lty=1,col="black",cex=1.5,lwd=2)
box(lwd=2.2)
dev.copy(png,file.path(outdir,paste0(paste0(inVar,"_Validation_RMSEP_R2_by_Component.png"))), 
         height=2800, width=4800,  res=340)
dev.off();
par(opar)

#################################################################################
### PLSR fit observed vs. predicted plot data
#calibration
cal.plsr.output <- data.frame(cal.plsr.data[, which(names(cal.plsr.data) %notin% "Spectra")], 
                              PLSR_Predicted=fit,
                              PLSR_CV_Predicted=as.vector(plsr.out$validation$pred[,,nComps]))
cal.plsr.output <- cal.plsr.output %>%
  mutate(PLSR_CV_Residuals = PLSR_CV_Predicted-get(inVar))
head(cal.plsr.output)
cal.R2 <- round(pls::R2(plsr.out,intercept=F)[[1]][nComps],2)
cal.RMSEP <- round(sqrt(mean(cal.plsr.output$PLSR_CV_Residuals^2)),2)

val.plsr.output <- data.frame(val.plsr.data[, which(names(val.plsr.data) %notin% "Spectra")],
                              PLSR_Predicted=as.vector(predict(plsr.out, 
                                                               newdata = val.plsr.data, 
                                                               ncomp=nComps, type="response")[,,1]))
val.plsr.output <- val.plsr.output %>%
  mutate(PLSR_Residuals = PLSR_Predicted-get(inVar))
head(val.plsr.output)
val.R2 <- round(pls::R2(plsr.out,newdata=val.plsr.data,intercept=F)[[1]][nComps],2)
val.RMSEP <- round(sqrt(mean(val.plsr.output$PLSR_Residuals^2)),2)

rng_quant <- quantile(cal.plsr.output[,inVar], probs = c(0.001, 0.999))
cal_scatter_plot <- ggplot(cal.plsr.output, aes(x=PLSR_CV_Predicted, y=get(inVar))) + 
  theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1, color="dark grey", 
                                          linetype="dashed", size=1.5) + xlim(rng_quant[1], rng_quant[2]) + 
  ylim(rng_quant[1], rng_quant[2]) +
  labs(x=paste0("Predicted ", paste(inVar), " (units)"),
       y=paste0("Observed ", paste(inVar), " (units)"),
       title=paste0("Calibration: ", paste0("Rsq = ", cal.R2), "; ", paste0("RMSEP = ", cal.RMSEP))) +
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

cal_resid_histogram <- ggplot(cal.plsr.output, aes(x=PLSR_CV_Residuals)) +
  geom_histogram(alpha=.5, position="identity") + 
  geom_vline(xintercept = 0, color="black", 
             linetype="dashed", size=1) + theme_bw() + 
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

rng_quant <- quantile(val.plsr.output[,inVar], probs = c(0.001, 0.999))
val_scatter_plot <- ggplot(val.plsr.output, aes(x=PLSR_Predicted, y=get(inVar))) + 
  theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1, color="dark grey", 
                                          linetype="dashed", size=1.5) + xlim(rng_quant[1], rng_quant[2]) + 
  ylim(rng_quant[1], rng_quant[2]) +
  labs(x=paste0("Predicted ", paste(inVar), " (units)"),
       y=paste0("Observed ", paste(inVar), " (units)"),
       title=paste0("Validation: ", paste0("Rsq = ", val.R2), "; ", paste0("RMSEP = ", val.RMSEP))) +
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

val_resid_histogram <- ggplot(val.plsr.output, aes(x=PLSR_Residuals)) +
  geom_histogram(alpha=.5, position="identity") + 
  geom_vline(xintercept = 0, color="black", 
             linetype="dashed", size=1) + theme_bw() + 
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

# plot cal/val side-by-side
scatterplots <- grid.arrange(cal_scatter_plot, val_scatter_plot, cal_resid_histogram, 
                             val_resid_histogram, nrow=2, ncol=2)
ggsave(filename = file.path(outdir,paste0(inVar,"_Cal_Val_Scatterplots.png")), 
       plot = scatterplots, device="png", width = 32, height = 30, units = "cm",
       dpi = 300)

################################################################################
### Generate Coefficient and VIP plots
vips <- spectratrait::VIP(plsr.out)[nComps,]

par(mfrow=c(2,1))
plot(plsr.out$coefficients[,,nComps], x=wv,xlab="Wavelength (nm)",
     ylab="Regression coefficients",lwd=2,type='l')
box(lwd=2.2)
plot(seq(Start.wave,End.wave,1),vips,xlab="Wavelength (nm)",ylab="VIP",cex=0.01)
lines(seq(Start.wave,End.wave,1),vips,lwd=3)
abline(h=0.8,lty=2,col="dark grey")
box(lwd=2.2)
dev.copy(png,file.path(outdir,paste0(inVar,'_Coefficient_VIP_plot.png')), 
         height=3100, width=4100, res=340)
dev.off();
par(opar)

################################################################################
if(grepl("Windows", sessionInfo()$running)){
  pls.options(parallel =NULL)
} else {
  pls.options(parallel = parallel::detectCores()-1)
}

seg <- 40
jk.plsr.out <- pls::plsr(as.formula(paste(inVar,"~","Spectra")), scale=FALSE, 
                         center=TRUE, ncomp=nComps, validation="CV", 
                         segments = seg, segment.type="interleaved", trace=FALSE, 
                         jackknife=TRUE, data=cal.plsr.data)
pls.options(parallel = NULL)

Jackknife_coef <- spectratrait::f.coef.valid(plsr.out = jk.plsr.out, 
                                             data_plsr = cal.plsr.data, 
                                             ncomp = nComps, inVar=inVar)
Jackknife_intercept <- Jackknife_coef[1,,,]
Jackknife_coef <- Jackknife_coef[2:dim(Jackknife_coef)[1],,,]

interval <- c(0.025,0.975)
Jackknife_Pred <- val.plsr.data$Spectra %*% Jackknife_coef + 
  matrix(rep(Jackknife_intercept, length(val.plsr.data[,inVar])), byrow=TRUE, 
         ncol=length(Jackknife_intercept))
Interval_Conf <- apply(X = Jackknife_Pred, MARGIN = 1, FUN = quantile, 
                       probs=c(interval[1], interval[2]))
sd_mean <- apply(X = Jackknife_Pred, MARGIN=1, FUN=sd)
sd_res <- sd(val.plsr.output$PLSR_Residuals)
sd_tot <- sqrt(sd_mean^2+sd_res^2)
val.plsr.output$LCI <- Interval_Conf[1,]
val.plsr.output$UCI <- Interval_Conf[2,]
val.plsr.output$LPI <- val.plsr.output$PLSR_Predicted-1.96*sd_tot
val.plsr.output$UPI <- val.plsr.output$PLSR_Predicted+1.96*sd_tot
head(val.plsr.output)

# JK regression coefficient plot
spectratrait::f.plot.coef(Z = t(Jackknife_coef), wv = wv, 
                          plot_label="Jackknife regression coefficients",position = 'bottomleft')
abline(h=0,lty=2,col="grey50")
box(lwd=2.2)
dev.copy(png,file.path(outdir,paste0(inVar,'_Jackknife_Regression_Coefficients.png')), 
         height=2100, width=3800, res=340)
dev.off();

# JK validation plot
rmsep_percrmsep <- spectratrait::percent_rmse(plsr_dataset = val.plsr.output, 
                                              inVar = inVar, 
                                              residuals = val.plsr.output$PLSR_Residuals, 
                                              range="full")
RMSEP <- rmsep_percrmsep$rmse
perc_RMSEP <- rmsep_percrmsep$perc_rmse
r2 <- round(pls::R2(plsr.out, newdata = val.plsr.data, intercept=F)$val[nComps],2)
expr <- vector("expression", 3)
expr[[1]] <- bquote(R^2==.(r2))
expr[[2]] <- bquote(RMSEP==.(round(RMSEP,2)))
expr[[3]] <- bquote("%RMSEP"==.(round(perc_RMSEP,2)))
rng_vals <- c(min(val.plsr.output$LPI), max(val.plsr.output$UPI))
par(mfrow=c(1,1), mar=c(4.2,5.3,1,0.4), oma=c(0, 0.1, 0, 0.2))
plotrix::plotCI(val.plsr.output$PLSR_Predicted,val.plsr.output[,inVar], 
                li=val.plsr.output$LPI, ui=val.plsr.output$UPI, gap=0.009,sfrac=0.004, 
                lwd=1.6, xlim=c(rng_vals[1], rng_vals[2]), ylim=c(rng_vals[1], rng_vals[2]), 
                err="x", pch=21, col="black", pt.bg=scales::alpha("grey70",0.7), scol="grey50",
                cex=2, xlab=paste0("Predicted ", paste(inVar), " (units)"),
                ylab=paste0("Observed ", paste(inVar), " (units)"),
                cex.axis=1.5,cex.lab=1.8)
abline(0,1,lty=2,lw=2)
legend("topleft", legend=expr, bty="n", cex=1.5)
box(lwd=2.2)
dev.copy(png,file.path(outdir,paste0(inVar,"_PLSR_Validation_Scatterplot.png")), 
         height=2800, width=3200,  res=340)
dev.off();

################################################################################
#output jackknife results

# JK Coefficents
out.jk.coefs <- data.frame(Iteration=seq(1,seg,1),Intercept=Jackknife_intercept,t(Jackknife_coef))
head(out.jk.coefs)[1:6]
write.csv(out.jk.coefs,file=file.path(outdir,paste0(inVar,'_Jackkife_PLSR_Coefficients.csv')),
          row.names=FALSE)

################################################################################
#export model output

print(paste("Output directory: ", getwd()))

# Observed versus predicted
write.csv(cal.plsr.output,file=file.path(outdir,paste0(inVar,'_Observed_PLSR_CV_Pred_',nComps,
                                                       'comp.csv')),row.names=FALSE)

# Validation data
write.csv(val.plsr.output,file=file.path(outdir,paste0(inVar,'_Validation_PLSR_Pred_',nComps,
                                                       'comp.csv')),row.names=FALSE)

# Model coefficients
coefs <- coef(plsr.out,ncomp=nComps,intercept=TRUE)
write.csv(coefs,file=file.path(outdir,paste0(inVar,'_PLSR_Coefficients_',nComps,'comp.csv')),
          row.names=TRUE)

# PLSR VIP
write.csv(vips,file=file.path(outdir,paste0(inVar,'_PLSR_VIPs_',nComps,'comp.csv')))

# confirm files were written to temp space
print("**** PLSR output files: ")
print(list.files(getwd())[grep(pattern = inVar, list.files(getwd()))])
################################################################################



