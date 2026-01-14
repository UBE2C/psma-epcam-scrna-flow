#!user/bin/Rscript

#=======================================#
#|                                     |#
#|    Smart-seq2 - Quality control II  |#
#|                                     |#
#=======================================#


#####################################################################################################################################################
#################################################################   Start Script  ###################################################################
#####################################################################################################################################################




#################################################################  Package management  ##################################################################
                                                                #                      #


# Define the necessary packages
cran_packages <- c("tidyverse", "rebus", "this.path")


# Install and load necessary packages
install_required_packages = function(packages) {
  for(i in seq_along(packages)) {
    cran_item <-  packages[i]

    if (!is.element(cran_item, installed.packages())) {
      message("The package: ", cran_item, " is not installed. Would you like to install these? \n1: YES\n2:NO")
      answer <- scan(file = "", what = "numeric", n = 1)

      if (answer == 1) {
        install.packages(cran_item)

      } else if (answer == 2) {
        message("The required packages will not be automatically installed. If you want to run this script install them manually.")

      } else {
        simpleError("Only the numbers 1 and 2 are accepted, please re-run the function.")

      }
      
    } else {
      message("The required package ", cran_item, " is already installed, you can continue with the script.")

    }

  }

}

install_required_packages(cran_packages)
lapply(cran_packages, library, character.only = TRUE)

#################################################################   End Section  ##################################################################






#################################################################         Set the working directory       ##################################################################
                                                                # then load and process the analysis files  #


# Set the working directory to the one the script is in
local_dir <- this.path::this.dir()
setwd(local_dir)


# List the smear analysis files
smears <- list.files(path = local_dir, pattern = "Smear Analysis Result.csv", )


# Read in the files 
Smear_lst <- list()
for (e in seq_along(smears)) {
  Smear_lst[[e]] <- read.csv(file = paste0(path, "/", smears[e]), sep = ";")
}
names(Smear_lst) <- smears
print(Smear_lst)


# Build a data frame and replace the comma separators with dots to properly represent floats
smear_df <- data.frame(do.call(rbind, smear_matrix_lst))
colnames(smear_df) <- str_split_fixed(colnames(Smear_lst[[1]]), ";", n = 10)
for (e in 1:nrow(smear_df)){
  smear_df[e, ] <- str_replace_all(smear_df[e, ], ",", ".")
}


# Convert the numeric columns into the proper type
for(e in 4:ncol(smear_df)){
  smear_df[, e] <- as.numeric(smear_df[, e])
}
glimpse(smear_df)


# Remove the ladder and blank rows
smear_df <- subset(smear_df, smear_df$`Sample ID` != "ladder" & smear_df$`Sample ID` != "Blank")
smear_df <- subset(smear_df, !is.nan(smear_df$`ng/uL`))

smear_df_clean <- mutate(smear_df, `fmole/ul` = `nmole/L`, `total fmole in 5 ul` = 5 * `fmole/ul`)
print(smear_df_clean)

# Save the resulting clean.csv
write_csv(smear_df_clean, file = paste0(path, "/", "Smear_analysis_clean_merged.csv"))

#################################################################   End Section  ##################################################################






#################################################################         Calculate the required dilutions       ##################################################################
                                                                #        and the required pipetting volumes      #


#NOTE: The strategy is to split the samples into two -> samples above the median concentration ([C]), will all be diluted to the median [C] and ones below the median [C]
#will all be concentrated on beads and if possible the multiplex will be diluted to the median [C]

# Split the data frame based on sample [C]
print(median(smear_df_clean$`nmole/L`))

below_med <- smear_df_clean[smear_df_clean$`nmole/L` < median(smear_df_clean$`nmole/L`), ]
print(below_med)

above_med <- smear_df_clean[smear_df_clean$`nmole/L` > median(smear_df_clean$`nmole/L`), ]
print(above_med)


# Further fraction low [C] samples so samples below 1nM can be excluded from the library
reseq_low_C <- below_med[below_med$`nmole/L` > 1, ]
print(reseq_low_C)

cut_samples <- below_med[below_med$`nmole/L` < 1, ]
print(cut_samples)


# Calculate the dilution factor for the high [C] samples (dilution factor = current[C]/ desired[C]) using the median [C] as desired[C]
reseq_high_C <- mutate(above_med, 
                       `dilution factor` = above_med$`nmole/L`/ median(smear_df_clean$`nmole/L`))
print(reseq_high_C)


# Save these dfs as .csv files
write_csv(reseq_low_C, file = paste0(path, "/", "Smear_analysis_below_median_nM.csv"))
write_csv(cut_samples, file = paste0(path, "/", "Smear_analysis_cut_samples.csv"))
write_csv(reseq_high_C, file = paste0(path, "/", "Smear_analysis_above_median_nM.csv"))

#################################################################   End Section  ##################################################################




###################################################################################################################################################
#################################################################   End Script  ###################################################################
###################################################################################################################################################