# Absorptance data processing

# Load the tidyr package
library(tidyr)
library(ggplot2)
library(dplyr)

# List all .txt files in the "Spectral_data" folder
file_list <- list.files(path = "./Dataset/French_Guiana/Spectral_data/", pattern = "\\.txt$", full.names = TRUE)

# Read each file and add the file name as a column
data_list <- lapply(file_list, function(file) {
  # Read the comma-delimited file
  df <- read.csv(file, header = FALSE) # Assuming no header in files
  # Rename the columns to "wavelength" and "value"
  colnames(df) <- c("wavelength", "value")
  # Add the file name as a new column
  df$file_name <- basename(file)
  return(df)
})

# Combine all the data frames into one
combined_data <- do.call(rbind, data_list)

#do some cleaning of the dataframe

# Separate the 'file_name' column into 'measurement', 'leaf_code', and 'obs_date'
combined_data <- combined_data %>%
  separate(file_name, into = c("measurement", "leaf_code", "obs_date"), sep = "_", remove = FALSE)

##### Julien's spectra processing steps include:

# "Removing observations in the range 891.4 : 891.7 which exhibit a clear bias"
combined_data = combined_data[combined_data$wavelength < 891.4 | combined_data$wavelength > 897.7, ]
#combined_data$value[combined_data$wavelength >= 891.4 & combined_data$wavelength <= 897.7] <- NA

#subset to only include wavelengths 490 - 950
#combined_data$wavelength <- round(combined_data$wavelength, digits = 0)
combined_data <- combined_data[combined_data$wavelength >= 490 & combined_data$wavelength <= 950,]

# since we only need reflectance data for now let us subset for that too
ref_data <- combined_data %>%
  filter(grepl("reflection", measurement))

# try their version but my way also
# Apply Savitzky-Golay filtering with a window size of 11 and polynomial degree 2
ref_data$smooth_value_sg <- signal::sgolayfilt(ref_data$value, p = 1, n = 51)

# Generate a continuous sequence of wavelengths from the minimum to the maximum wavelength in the dataset
wavelength_range <- seq(min(ref_data$wavelength), max(ref_data$wavelength), by = 1)

# Now perform interpolation for the missing wavelengths
ref_data_interpolated <- ref_data %>%
  group_by(file_name, measurement, leaf_code, obs_date) %>%
  do({
    # Interpolate over the continuous range of wavelengths
    interp <- approx(x = .$wavelength, y = .$smooth_value_sg, xout = wavelength_range, method = "linear", rule = 2)
    
    # Return interpolated data
    data.frame(wavelength = interp$x, smooth_value_sg = interp$y, file_name = unique(.$file_name), measurement = unique(.$measurement), leaf_code = unique(.$leaf_code), obs_date = unique(.$obs_date))
  }) %>%
  ungroup()

#now subset for one individual just to test out what the different files are since they are actually not what I expected.
sub_data <- subset(ref_data_interpolated, leaf_code == "1.1.1")

ggplot(sub_data, aes(wavelength, smooth_value_sg, col = leaf_code))+
  geom_point(col = "black")

#now plot all
data <- data.frame(ref_data_interpolated)

#save this dataset
write.table(data, "./Outputs/Reflectance_Paracou_Processed.csv", row.names = FALSE)


