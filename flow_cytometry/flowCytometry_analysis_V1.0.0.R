#!user/bin/Rscript

#==================================#
#|                                |#
#|    FACS results - analysis     |#
#|                                |#
#==================================#






#################################################################  Setting up, formatting and  ##################################################################
                                                                #  cleaning the relevant data  #




## Prepare for the analysis by loading the necessary libraries and reading in the .fcs files
cran_packages <- c("tidyverse", "stringi", "BiocManager",
              "scales", "RCurl", "cowplot", "rebus", "ggsci", "ggsignif",
              "progress", "metap", "doSNOW", "foreach", "scCustomize",
              "ggpubr", "R.utils", "devtools", "remotes", "rstudioapi", "scales")

bioconductor_packages <- c("flowCore", "flowStats", "flowViz", "ggcyto")


# This function will check if the required packages are installed and if not, it installs them
# NOTE: by default the function will only check for CRAN packages, so if you want to install bioconductor
# or github packages you will have to feed those in separately
install_required_packages = function (CRAN_package_lst = NULL, BioC_package_lst = NULL, github_package_lst = NULL, github_install_param = NULL) {
    if (is.null(CRAN_package_lst) == FALSE) {
        for (cran_item in CRAN_package_lst) {
            if (is.element(cran_item, installed.packages()) == FALSE) {
                message("The package: ", cran_item, " is not installed. Installing package...")
                install.packages(cran_item)
            } else {
                message("The package: ", cran_item, " is installed. Load the package using the library function.")
            }
        }

    } else {
        message("No CRAN packages were requested.")
    }
      
    if (is.null(BioC_package_lst) == FALSE) {
        for (bioc_item in BioC_package_lst) {
            if (is.element(bioc_item, installed.packages()) == FALSE) {
                message("The package: ", bioc_item, " is not installed. Installing package...")
                BiocManager::install(bioc_item)
            } else {
                message("The package: ", bioc_item, " is installed. Load the package using the library function.")
            }
        }

    } else {
        message("\n", "No Bioconductor packages were requested.")

    }

    if (is.null(github_package_lst) == FALSE) {
        for (github_item in github_package_lst) {
            if (is.element(github_item, installed.packages()) == FALSE) {
                message("The package: ", github_item, " is not installed. Installing package...")
                for (install_param in github_install_param) {
                    devtools::install_github(install_param)
                }
            } else {
                message("The package: ", github_item, " is installed. Load the package using the library function.")
            }
        }

    } else {
        message("\n", "No github packages were requested.")
    }  
}
install_required_packages(CRAN_package_lst = cran_packages, BioC_package_lst = bioconductor_packages)


# Used library packages
packages <- c("tidyverse", "rebus", "progress", "ggpubr", "ggsci", "rstudioapi", "doSNOW", "foreach", "flowCore", "flowStats", "flowViz", "ggcyto")

lapply(packages, library, character.only = TRUE)

# Define the working directory and the directory containing the .fcs and other necessary files
path <- this.path::here()
setwd(path)
FACS_files <- c("/FACS_fcs_csv_files")


# Define the script version for output files, for easier version control, and make a dedicated folder in the WD
define_version_create_folder = function(directory_name = "FACS") {
    #define the current script version
    #NOTE: the version numbers must match the defined pattern. I  updated the function to accept sub-versions below two digits.
    if (nchar(stringr::str_extract(rstudioapi::getActiveDocumentContext()$path, "_V.*")) == 7) {
        script_version <- paste0(directory_name, stringr::str_extract(rstudioapi::getActiveDocumentContext()$path, "_V" %R% rebus::one_or_more(ASCII_ALNUM) %R% one_or_more(PUNCT) %R%
                                      one_or_more(ASCII_ALNUM)))
    } else {
        script_version <- paste0(directory_name, stringr::str_extract(rstudioapi::getActiveDocumentContext()$path, "_V" %R% rebus::one_or_more(ASCII_ALNUM) %R% one_or_more(PUNCT) %R%
                                      one_or_more(ASCII_ALNUM) %R% one_or_more(PUNCT) %R% one_or_more(ASCII_ALNUM)))
    }
    
    #check if the current version already has a folder or not and if not it creates one
    if (dir.exists(paths = paste0(script_version, "/", "Plots")) == TRUE) {
        message("The current script directory is already present. No new directory will be created.")
    } else {
        message("The current script version has no directory yet.", "\n", "Creating one now using the path:", paste0(getwd(), script_version, "/", "Plots"))
        dir.create(path = paste0(script_version, "/", "Plots"), recursive = TRUE)
        
    }
    assign("script_version", script_version, envir = .GlobalEnv)
}
define_version_create_folder()


# Load fcs files from the wd in a recursive manner
load_fcs = function(path = getwd(), file_pattern = ".fcs$", test_set = FALSE, full_path = TRUE, recursive = TRUE) {
  files <- list.files(path = path, pattern = file_pattern, full.names = full_path, recursive = recursive)
  sNames <- paste0(stringr::str_extract(files, "[0-9]+[-]+[0-9]+[-]+[0-9]+"), " ",
                   stringr::str_extract(files, "Specimen.*."))
  
  if (test_set == TRUE) {
      test_lst <- list()
      for (e in 1:5) {
          test_lst[[e]] <- flowCore::read.FCS(files[e])
      }
      names(test_lst) <- sNames[1:5]
      message("A test set was requested, compiling a 5 element set into the object test_lst and the SortList")
      
  } else {
      message("No test_set was requested, Compiling the FCS files into the SortList")
  }
  
  
  progressBar <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                  total = length(files),
                                  complete = "=",
                                  incomplete = "-",
                                  current = ">",
                                  clear = FALSE,
                                  width = 100)
  
  sList <- list()
  for (e in seq_along(files)) {
    progressBar$tick()
      
    sList[[e]] <- flowCore::read.FCS(files[e])
  }
  names(sList) <- sNames
  
  if (test_set == TRUE)  {
      assign("test_lst", test_lst, envir = .GlobalEnv)
      assign("SortList", sList, envir = .GlobalEnv)
      message("Both the test_set and the SortList have been compiled into the following objects: test_list and SortList")
  } else {
      assign("SortList", sList, envir = .GlobalEnv)
      message("Task complete, the FCS files were compiled into the following object: SortList")
  }
  
}
load_fcs(path = paste0(path, FACS_files))


# Clean up unnecessary variables from this section
rm(bioconductor_packages, cran_packages, packages, FACS_files)

#################################################################   End Section  ##################################################################






#################################################################  Prepare the control FlowFrames for  ##################################################################
                                                                #         data normalization           #




## Load the saved unified dataframe if not present in the .GlobalEnv
if (!exists("ctrlMedian_df", envir = .GlobalEnv)) {
  
  if (file.exists(paste0(script_version, "/", "C4-2_EpCAM-PSMA-median_MFI_for_normalization", "_", script_version, ".csv"))) {
    
    ctrlMedian_df <- read_csv(paste0(script_version, "/", "C4-2_EpCAM-PSMA-median_MFI_for_normalization", "_", script_version, ".csv"))
    ctrlMedian_df <- as.data.frame(ctrlMedian_df)
    message("The C4-2 control median dataframe was successfully loaded into the .GlobalEnv.")
    
  } else {
    
    message("The C4-2 control median dataframe is not loaded in the .GlobalEnvironment and the save file: \n",
            paste0("C4-2_EpCAM-PSMA-median_MFI_for_normalization", "_", script_version, ".csv"), "\n",
            "was not found in the current script version directory. Please run the code block below to generate it.")
    
  } 
  
} else {
  
  message("ctrlMedian_df is already in the .Global.Env, no further actions are necessary.")
  
}




## Extracting the control experiments done with C4-2 cells at the beginning every sort session


# Subset the control samples for normalization (FF stands for FlowFrame)
sCtrl <- names(SortList)[grep(pattern = "C4-2", x = names(SortList))]
sCtrl <- sCtrl[!grepl(pattern = "INX", x = sCtrl)]

ctrlFFs <- SortList[c(sCtrl)]
ctrlSumm <- lapply(ctrlFFs, summary)
channel_names <- colnames(ctrlFFs[[1]]@exprs)
names(channel_names) <- NULL
names(ctrlFFs)[30] <- "2022-05-23 Specimen_001_C4-2 pos ctr_003.fcs"




## Extracting the dataframes containing the MFI values from the ctrlFlowFrames (ctrlFFs) to check the distribution


# This loop will re-assign the MFI values from the ctrlFFs to a new list, and makes a unified df
ctrlFF_dfs <- vector(mode = "list", length = length(ctrlFFs))
for(e in seq_along(ctrlFFs)){
    ctrlFF_dfs[[e]] <- ctrlFFs[[e]]@exprs
}
ctrlFF_dfs <- do.call(rbind, ctrlFF_dfs)
ctrlFF_dfs <- as.data.frame(ctrlFF_dfs)
head(ctrlFF_dfs)


# I will clean the colnames and subset the ctrlFF_dfs to get only the fluorophore channels
ctrlFF_dfs <- mutate(ctrlFF_dfs, "Time" = NULL)
colnames(ctrlFF_dfs) <- c("FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H",
                           "SSC.W", "EpCAM", "DAPI_CD45", "Caspase_3.7", "PSMA")

ctrlFF_dfs_Fluo <- ctrlFF_dfs[, c(7, 10)]
head(ctrlFF_dfs_Fluo)


# To see the data distribution and spot if there are some values which are strong outliers I transform the dataframe to a long version
# and I make a density plot.

# Note: it is not recommended to use the mean values for MFI data analysis since it normally does not follow normal distribution
# instead one should use the median values, so all the calculations and normalization will be done by the median values
long_ctrlFF_dfs_Fluo <- pivot_longer(ctrlFF_dfs_Fluo, cols = c(EpCAM, PSMA),
                             names_to = "Channels", values_to = "MFI_values")

# Density plot to see if the data is normally distributed or not (probably not)
ggplot(data = long_ctrlFF_dfs_Fluo,
       aes(x = log10(MFI_values), fill = Channels)) +
    #scale_x_continuous(limits = c(0, 6000)) +
    geom_density(alpha = 0.5) +
    theme_classic()
rm(long_ctrlFF_dfs_Fluo)




## Calculating the medians of the control samples for normalization and QC


# Calculate the median MFI values for the control samples which will be used for normalization and QC
medianFlowFrame = function(FlowFrameList) {
    tmp_lst <- vector(mode = "list", length = length(FlowFrameList))
    for (i in seq_along(FlowFrameList)) {
        
        tmp_lst[[i]] <- FlowFrameList[[i]]@exprs
    }
    
    channel_names <- colnames(FlowFrameList[[1]]@exprs)
    names(channel_names) <- NULL
    df <- data.frame(matrix(nrow = 1, ncol = length(channel_names)))
    colnames(df) <- channel_names
    
    median_lst <- vector(mode = "list", length = length(tmp_lst))
    for (e in seq_along(tmp_lst)){
        median_lst[[e]] <- df
        for (i in seq_len(ncol(tmp_lst[[1]]))) {
            m <- median(tmp_lst[[e]][, i])
            median_lst[[e]][i] <- m
        }
    }
    names(median_lst) <- names(FlowFrameList)
    
    return(median_lst)  
}
ctrlMedian <- medianFlowFrame(ctrlFFs)
head(ctrlMedian)


# Due to the fact that Both Martina and Me ran/used the same control set on the same day, there might be duplicates in the dataframe
# I write a function which check if elements of the list are redundant/equal. If yes, remove them.
equal_elements = function(data, remove = FALSE) {
  tmp <- vector(mode = "character", length = length(data))
  for (i in seq_along(data)) {
    if (i < length(data)) {
      if (all(data[[i]] == data[[i + 1]])) {
        tmp[i] <- paste0("Elements ", i, " and ", i + 1, " are equal")
      } else {
        tmp[i] <- paste0("Elements ", i, " and ", i + 1, " are not equal")
      }
    } else {
      break
    }  
  }
  
  if(remove == TRUE){
    to_remove <- which(grepl("are equal", tmp)) # which(...): This function then finds the indices of TRUE values in the logical vector returned by grepl(...). These indices correspond to the elements that match the pattern.
  
    if (length(to_remove) > 0) {
        message(paste("Elements", paste(to_remove, collapse = ", "), "are duplicated so one of the elements will be removed"))
        data_2 <- data[-to_remove] # this -to_remove will remove the selected values
        } else {
        message("no duplicated elements found, returning the unmodified object")
        }
    old_object_name <- deparse(substitute(data))
    new_object_name <- paste0(old_object_name, "_clean")
    assign(new_object_name, data_2, envir = .GlobalEnv)
    message(paste("Returning a new, clean object with the name:", paste0(old_object_name, "_clean")))

  } else {
    return(tmp)

  }
  
}
equal_elements(ctrlMedian, remove = TRUE)  
head(ctrlMedian_clean)


# This function merges the ctrl dataframes (MFI values) for easier visualization
ctrlMedian_df <- do.call(rbind, ctrlMedian_clean)


# This snippet will remove the Time col. and changes the channel names to the easier to understand ones
# NOTE: here are the annotated channel names:
# 640-670/30 - Incucyte Caspase 3/7
# 405-450/50 - DAPI + CD45
# 488-530/30 - EpCAM
# 561-780/60 - PSMA
ctrlMedian_df <- mutate(ctrlMedian_df, "Time" = NULL)
df_colnames <- c("FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H",
                 "SSC.W", "EpCAM", "DAPI_CD45", "Caspase_3.7", "PSMA")
colnames(ctrlMedian_df) <- df_colnames
head(ctrlMedian_df)

ctrlMedian_Fluo <- ctrlMedian_df[,c(7, 10)]
head(ctrlMedian_Fluo)

rm(ctrlSumm, ctrlMedian, ctrlMedian_clean, ctrlFF_dfs)




## Calculating the distribution values like MAD (Median Absolute Deviation) and MADM (Median Absolute Deviation of the medians) for QC

# 1st calculating QC values
calculating_grand_median = function(data) {
    gMedian <- data.frame(matrix(nrow = 1, ncol = ncol(data)))
    for (e in seq_len(ncol(data))) {
        gMedian[, e] <- median(data[, e])
    }
    colnames(gMedian) <- colnames(data)
    message("The grand median has been calculated from the median values of the supplied dataframe, and a new object has been created: ctrlFluo_GrMedian")
    assign("ctrlFluo_GrMedian", gMedian, envir = .GlobalEnv)
}
calculating_grand_median(ctrlMedian_Fluo)
head(ctrlFluo_GrMedian)


calculating_individual_MAD = function(data, col_names = NULL) {
    #Extract the MFI dfs from the FFs into a temp list
    FF_df_lst <- vector(mode = "list", length = length(data))
    for (e in seq_along(data)) {
        FF_df_lst[[e]] <- data[[e]]@exprs
    }
    
    #Calculate the MAD values of each column of each df in the list
    df <- data.frame(matrix(nrow = 1, ncol = ncol(FF_df_lst[[1]])))
    MAD_lst <- vector(mode = "list", length = length(FF_df_lst))
    for (e in seq_along(FF_df_lst)) {
        MAD_lst[[e]] <- df
        for (i in seq_len(ncol(FF_df_lst[[1]]))) {
            mad <- mad(FF_df_lst[[e]][, i])
            MAD_lst[[e]][[i]] <- mad
        
        }
         
    }
    names(MAD_lst) <- names(data)
    ctrlFluo_ind_MAD <- do.call(rbind, MAD_lst)
    
    if (is.null(col_names)) {
        message("No colname was assigned.")
    } else {
        if (length(col_names) == 10) {
            message("The intended colname vector is shorter than the number of columns in the new object.")
            message("I assume you intended to remove the last column: Time, and add the desired channel names.")
            message("Therefore, removing last column and adding proper colnames.")
            ctrlFluo_ind_MAD <- ctrlFluo_ind_MAD[, -ncol(ctrlFluo_ind_MAD)]
            colnames(ctrlFluo_ind_MAD) <- col_names
        } else if (length(col_names) == 11) {
            message("The intended colname vector does match the column number in the new df so I assume you want to keep the Time column")
            message("Therefore I will just assign the new colnames")
        } else {
            message("The intended colname does not match the number of columns (11 with Time 10 without Time) so please check the intended colnames.")
            message("Nothing will be assigned.")
        }
               
                
    }
    #return(ctrlFluo_ind_MAD)
    assign("ctrlFF_ind_MAD", ctrlFluo_ind_MAD, envir = .GlobalEnv)
}
calculating_individual_MAD(ctrlFFs, col_names = df_colnames)
print(ctrlFF_ind_MAD)

calculating_total_MAD = function(data) {
    MAD <- numeric(length = ncol(data))    
    ctrlFluo_MAD <- data.frame(matrix(ncol = ncol(data), nrow = 1))
    
    for (i in seq_len(ncol(data))) {
        MAD[i] <- mad(data[, i])
        ctrlFluo_MAD[, i] <- MAD[i]
        
    }  
    
    colnames(ctrlFluo_MAD) <- colnames(data)
    assign("ctrlFluo_total_MAD", ctrlFluo_MAD, envir = .GlobalEnv)
} 
calculating_total_MAD(ctrlFF_dfs_Fluo)
print(ctrlFluo_total_MAD)

calculating_MADM = function(data) {
    
    MADM <- data.frame(matrix(nrow = 1, ncol = ncol(data)))
    for (i in seq_len(ncol(data))) {
        MADM[, i] <- mad(data[, i])
        
    }  
    colnames(MADM) <- colnames(data)
    message("A MADM calculation complete, a new object is created in your workspace: ctrlFluo_MADM")
    assign("ctrlFluo_MADM", MADM, envir = .GlobalEnv)
} 
calculating_MADM(ctrlFF_ind_MAD[, c(7, 10)])
print(ctrlFluo_MADM)



# Visualizing the Median distribution and the MADM values
# I reformat the means/medians df to better visualize fluorophore data
## NOTE: for final figure printing disable the title
long_ctrlMed_F <- tidyr::pivot_longer(ctrlMedian_Fluo, cols = c(EpCAM, PSMA),
                               names_to = "Channels", values_to = "MFI_medians")
head(long_ctrlMed_F)

ctrlMedian_dist_box_p <- ggplot2::ggplot(data = long_ctrlMed_F,
                                    aes(fill = Channels, x = Channels, 
                                        y = log10(MFI_medians))) +
                                    stat_boxplot() +
                                    geom_dotplot(show.legend = FALSE, stackdir = "center", binaxis = "y",
                                                 binwidth = 0.05, alpha = 0.5) +
                                    geom_pwc(method = "wilcox_test", label = "p = {p}", label.size = 6) +
                                    scale_fill_tron() +
                                    #ggtitle("control EpCAM and PSMA median MFI values") +
                                    theme_classic() +
                                    theme(plot.title = element_text(hjust = 0.5),
                                          text = element_text(size = 24))
print(ctrlMedian_dist_box_p)

ggsave(filename = paste0("Ctrl_samples_medians", "_", script_version, ".png"), ctrlMedian_dist_box_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)

rm(ctrlMedian_dist_box_p)


ctrlMedian_dist_dens_p <- ggplot2::ggplot(data = long_ctrlMed_F,
                                aes(x = log10(MFI_medians), fill = Channels)) +
                                #scale_x_continuous(limits = c(0, 6000)) +
                                geom_density(alpha = 0.5) +
                                scale_fill_tron() +
                                #ggtitle("control EpCAM and PSMA median MFI values - density plot") +
                                labs(x = expression("log" [10] *" median MFI values")) +
                                theme_classic() +
                                theme(plot.title = element_text(hjust = 0.5),
                                      text = element_text(size = 24))
print(ctrlMedian_dist_dens_p)

ggsave(filename = paste0("Ctrl_samples_medians_density", "_", script_version, ".png"), ctrlMedian_dist_dens_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(ctrlMedian_dist_dens_p)


# Creating a summary stat df for error bar visualization
long_ctrlFluo_GrMedian <- tidyr::pivot_longer(ctrlFluo_GrMedian, cols = c(EpCAM, PSMA),
                               names_to = "Channels", values_to = "MFI_grand_medians")
long_ctrlFluo_MADM <- tidyr::pivot_longer(ctrlFluo_MADM, cols = c(EpCAM, PSMA),
                                   names_to = "Channels", values_to = "MADM")
ctrlFluo_sum_stat <- cbind(long_ctrlFluo_GrMedian, long_ctrlFluo_MADM[, 2])
print(ctrlFluo_sum_stat)
rm(long_ctrlFluo_GrMedian, long_ctrlFluo_MADM)

# Plotting the Grand Median and MADM values on a bar chart with the error bars
MADM_ctrlFluo_p <- ggplot2::ggplot(data = ctrlFluo_sum_stat,
                         aes(x = Channels, y = log10(MFI_grand_medians), fill = Channels)) +
    geom_col(alpha = 0.5, color = "black") +
    geom_errorbar(aes(x = Channels, ymin = log10(MFI_grand_medians - MADM), ymax = log10(MFI_grand_medians + MADM),
                      width = 0.1, color = "red")) +
    scale_fill_tron() +
    labs(x = expression("Channels"), y = expression("log"[10]*" Grand median MFI values")) +
    #ggtitle("MADM plot of multiple control runs") +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(MADM_ctrlFluo_p)

ggsave(filename = paste0("Ctrl_samples_medians_MADM", "_", script_version, ".png"), MADM_ctrlFluo_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(MADM_ctrlFluo_p)


# Add the sample names as a column and save the created control-medians df
ctrlMedian_df <- mutate(ctrlMedian_df, Names = rownames(ctrlMedian_df), .before = FSC.A)
head(ctrlMedian_df)


# Save the created control-medians df
write_csv(x = ctrlMedian_df, file = paste0(script_version, "/", "C4-2_EpCAM-PSMA-median_MFI_for_normalization", "_", script_version, ".csv"))


# Removing the unnecessary QC objects to remove clutter and free memory
rm(ctrlFF_dfs_Fluo, ctrlFF_ind_MAD, ctrlFluo_total_MAD, ctrlFluo_GrMedian, ctrlFluo_sum_stat, ctrlMedian_Fluo, long_ctrlMed_F, ctrlFluo_MADM, e,
df_colnames, channel_names)

#################################################################   End Section  ##################################################################






#################################################################  Prepare the sample FlowFrames for  ##################################################################
                                                                #          data normalization         #




## Load the saved unified dataframe if not present in the .GlobalEnv
if (!exists("iSort_unif", envir = .GlobalEnv)) {
  
  if (file.exists(paste0(script_version, "/", "Unified_index_sorts_with_metadata", "_", script_version, ".csv"))) {
    
    iSort_unif <- read_csv(paste0(script_version, "/", "Unified_index_sorts_with_metadata", "_", script_version, ".csv"))
    iSort_unif <- as.data.frame(iSort_unif)
    message("The iSort_unif dataframe was successfully loaded into the .GlobalEnv.")
  } else {
    
    message("The iSort_unif dataframe is not loaded in the .GlobalEnvironment and the save file: \n",
            paste0("Unified_index_sorts_with_metadata", "_", script_version, ".csv"), "\n",
            "was not found in the current script version directory. Please run the code block below to generate it.")
    
  } 
  
} else {
  
  message("iSort_unif is already in the .Global.Env, no further actions are necessary.")
  
}




## Preparing the real samples for normalization


# I will subset the sort list to an index sort list iSorts
iSorts <- SortList[grepl("INX", names(SortList))]

#I will correct some indexing mistakes in iSorts
grep("INX_ALM80C8-004-1", names(iSorts))
names(iSorts)[[18]]
grep("INX_ALM76C6-003-1", names(iSorts))
names(iSorts)[[16]]
names(iSorts)[[17]]

keyword(iSorts[[16]])$"$FIL" #Correct
keyword(iSorts[[17]])$"$FIL"
keyword(iSorts[[18]])$"$FIL"

keyword(iSorts[[16]])$"$FIL" <- "Specimen_001_INX_ALM76C6-003-1_strip_1_001_006.fcs"

names(iSorts)[[18]] <- "2022-04-13 Specimen_001_INX_ALM80C8-004-2 strip 1_001_005.fcs"
keyword(iSorts[[18]])$"$FIL" <- "Specimen_001_INX_ALM80C8-004-2 strip 1_001_005.fcs"

grep("ALM83C8", names(iSorts))
names(iSorts)[[48]]
keyword(iSorts[[48]])$"$FIL" <- "Specimen_001_INX_ALM83C8-001-3_str_1_001_005.fcs"

grep("LM78C7-005-1", names(iSorts))
names(iSorts)[[153]]
names(iSorts)[[154]]
keyword(iSorts[[153]])$"$FIL" <- "Specimen_001_INX_ALM78C7-005-1_001_005.fcs"
keyword(iSorts[[154]])$"$FIL" <- "Specimen_001_INX_ALM78C7-005-1_002_006.fcs"


# I will clean the list by getting rid of the C4-2 and leukocyte control sort samples as they are unnecessary
iSorts <- iSorts[!grepl("C4-2", names(iSorts))]
iSorts <- iSorts[!grepl("leukocyte", names(iSorts))]


# This function will correct the incorrect sample names in the main FlowFrames (.fcs)
fix_sample_name_in_FFs = function(data, incorrect_name_pattern) {
    # first we list out the previously corrected sample names
    specimen_names <- list()
    for (e in seq_along(data)) {
        specimen_names[[e]] <- stringr::str_extract(flowCore::keyword(data[[e]])$FIL, "Specimen*.*") 
    }
    
    # then this loop and if statement will check if the sample name is incorrect or not
    # if it is, it will be fixed based on the correct sample name list
    incorrect_lst_positions <- c()
    modif_data <- data
    for (e in seq_along(data)) {
        if(stringr::str_detect(flowCore::keyword(data[[e]])$"$FIL", incorrect_name_pattern) == TRUE) {
            incorrect_lst_positions[e] <- e
            flowCore::keyword(modif_data[[e]])$"$FIL" <- specimen_names[[e]]
        } else {
            #do nothing
        }
    }
    incorrect_lst_positions <- incorrect_lst_positions[!is.na(incorrect_lst_positions)]
    message(paste0("The following list elements were corrected:", deparse(substitute(incorrect_lst_positions))))
    return(modif_data)
}
iSorts <- fix_sample_name_in_FFs(iSorts, "leukocyte")
iSorts <- fix_sample_name_in_FFs(iSorts, "LM66C3-003")
iSorts <- fix_sample_name_in_FFs(iSorts, "LM68C6-005")
iSorts <- fix_sample_name_in_FFs(iSorts, "LM67C6-011")
iSorts <- fix_sample_name_in_FFs(iSorts, "LM81C6-017")
iSorts <- fix_sample_name_in_FFs(iSorts, "LM75C7-006")

# I noticed that the colname order for the ctrl flow frames are different
# than the colname order of the index flow frames so I made this function
# to check if all the index flowframes have the same colname order or not
check_colnames = function(FFs, colname_order) {
    tmp <- c()
    for (e in seq_along(FFs)) {
        if (length(summary(unique(colnames(FFs[[e]]@exprs) == colname_order))) == 2) {
            tmp[e] <- c("it is fine")
        } else {
            tmp[e] <- c("this one is off")
        }
    }
    return(tmp)
}

index_colnames <- colnames(iSorts[[1]]@exprs)
ctrl_colnames <- colnames(ctrlFFs[[1]]@exprs)
names(index_colnames) <- NULL
names(ctrl_colnames) <- NULL
check_colnames(iSorts, ctrl_colnames)
rm(index_colnames, ctrl_colnames)


# It seems that for each dataframe the colnames are off between the ctrl
# and the sample flowframes, therefore I will check the order manually
colnames(ctrlFFs[[1]]) == colnames(iSorts[[1]])

# NOTE: after checking it seems that it is only the scatter channels which are in different order, and the fluorophore channels are fine
# therefore no re-ordering is necessary as I will only look at the fluorophore channels


# This snippet allows the extraction of index sort tables from the index FFs and to 
# unify them in one main df
iSort_dfs <- list()
for (e in seq_along(iSorts)) {
    iSort_dfs[[e]] <- getIndexSort(iSorts[[e]])
}
names(iSort_dfs) <- names(iSorts)

iSort_df <- do.call(rbind, iSort_dfs)


# I will reformat the iSort_df and remove columns which are not important for the normalization (will also remove scatter columns)
colnames(iSort_df)
iSort_df <- iSort_df[, -c(2, 3, 4,5, 6, 7, 12, 13, 14)]
iSort_colnames <- c("Name", "EpCAM", "DAPI_CD45", "Caspase_3.7", "PSMA")
colnames(iSort_df) <- iSort_colnames
rm(iSort_colnames)




## Here I'm gathering additional data for my iSort_df to make it more comprehensive

# This gets the sort dates
iDates <- vector(mode = "character", length = length(iSorts))
for (i in seq_along(iSorts)) {
  FF_dates <- flowCore::keyword(iSorts[[i]])$"$DATE"
  conv_date <- as.Date(FF_dates, format = "%d-%b-%Y")
  iDates[i] <- format(conv_date, format = "%Y-%m-%d") 
}
rm(FF_dates, conv_date)


# This snippet subtracts the sample IDs 
# important note: there a too many naming mistakes in here, so I had to write a
# super complicated string matching pattern. Leave the pattern order, otherwise it won't work
sampleID <- str_extract(iSort_df$Name,
                        "INX_" %R% one_or_more(PUNCT) %R% one_or_more(ALNUM) %R%
                            one_or_more(PUNCT) %R% one_or_more(ALNUM) %R%
                            one_or_more(PUNCT) %R% one_or_more(ALNUM) %R%
                            one_or_more(PUNCT) %R% one_or_more(ALNUM) %|%
                            "INX_" %R% one_or_more(ALNUM) %R%
                            SPACE %R% one_or_more(ALNUM) %R%
                            one_or_more(PUNCT) %R% one_or_more(ALNUM) %R%
                            one_or_more(PUNCT) %R% one_or_more(ALNUM) %|%
                            "INX_" %R% one_or_more(ALNUM) %R%
                            "-" %R% one_or_more(ALNUM) %R%
                            "-" %R% one_or_more(ALNUM) %R%
                            "-" %R% one_or_more(ALNUM) %|%
                            "INX_" %R% one_or_more(ALNUM) %R%
                            one_or_more(PUNCT) %R% one_or_more(ALNUM) %R%
                            one_or_more(PUNCT) %R% one_or_more(ALNUM))
sampleID <- str_remove(sampleID, "INX_")


# This part loads the file names separated by response group so the response can be subtracted
# Non-responders
nonResp <- list.files(path = paste0(path, "/", "FACS_fcs_csv_files", "/", "Non-responder indexes"))
nonResp <- nonResp[str_detect(nonResp, "INX")]
nonResp <- str_extract(nonResp,
                       "INX_" %R% one_or_more(PUNCT) %R% one_or_more(ALNUM) %R%
                           one_or_more(PUNCT) %R% one_or_more(ALNUM) %R%
                           "-" %R% one_or_more(ALNUM) %R%
                           "-" %R% one_or_more(ALNUM) %|%
                           "INX_" %R% one_or_more(ALNUM) %R%
                           "-" %R% one_or_more(ALNUM) %R%
                           "-" %R% one_or_more(ALNUM) %|%
                           "INX_" %R% one_or_more(ALNUM) %R%
                           SPACE %R% one_or_more(ALNUM) %R%
                           "-" %R% one_or_more(ALNUM) %R%
                           "-" %R% one_or_more(ALNUM))
nonResp <- str_remove(nonResp, "INX_")
nonResp <- unique(nonResp)


# Responders
Resp <- list.files(path = paste0(path, "/", "FACS_fcs_csv_files", "/", "Responder indexes"))
Resp <- Resp[str_detect(Resp, "INX")]
Resp <- str_extract(Resp,
                    "INX_" %R% one_or_more(PUNCT) %R% one_or_more(ALNUM) %R%
                        one_or_more(PUNCT) %R% one_or_more(ALNUM) %R%
                        "-" %R% one_or_more(ALNUM) %R%
                        "-" %R% one_or_more(ALNUM) %|%
                        "INX_" %R% one_or_more(ALNUM) %R%
                        "-" %R% one_or_more(ALNUM) %R%
                        "-" %R% one_or_more(ALNUM) %|%
                        "INX_" %R% one_or_more(ALNUM) %R%
                        SPACE %R% one_or_more(ALNUM) %R%
                        "-" %R% one_or_more(ALNUM) %R%
                        "-" %R% one_or_more(ALNUM))
Resp <- str_remove(Resp, "INX_")
Resp <- unique(Resp)


# These samples are not in the final seq data, therefore I will not use them for the FACS analysis
Unkn <- list.files(path = paste0(path, "/", "FACS_fcs_csv_files", "/", "Unclassified indexes"))
Unkn <- Unkn[str_detect(Unkn, "INX")]
Unkn <- str_extract(Unkn,
                    "INX_" %R% one_or_more(PUNCT) %R% one_or_more(ALNUM) %R%
                        one_or_more(PUNCT) %R% one_or_more(ALNUM) %R%
                        "-" %R% one_or_more(ALNUM) %R%
                        "-" %R% one_or_more(ALNUM) %|%
                        "INX_" %R% one_or_more(ALNUM) %R%
                        "-" %R% one_or_more(ALNUM) %R%
                        "-" %R% one_or_more(ALNUM) %|%
                        "INX_" %R% one_or_more(ALNUM) %R%
                        SPACE %R% one_or_more(ALNUM) %R%
                        "-" %R% one_or_more(ALNUM) %R%
                        "-" %R% one_or_more(ALNUM))
Unkn <- str_remove(Unkn, "INX_")
Unkn <- unique(Unkn)

# This loop basically does the heavy lifting, by checking if the sampleID is
# part of a response group or not. If yes it assigns the proper group into a string

# NOTE: additionally, I realized that the sample naming is fixed for the sorted index files
# from which I borrowed the names, so I fixed the naming for the sample IDs
# however, this fixed sample ID should not be used outside of this snippet....
fixed_ID <- str_replace_all(sampleID, pattern = "ALM-70C8-006-2",
                            replacement = "ALM70C8-006-2")
fixed_ID <- str_replace_all(fixed_ID, pattern = "ALM69-C5-003-3",
                            replacement = "ALM69C5-003-3")
fixed_ID <- str_replace_all(fixed_ID, pattern = "ALM71C9-04-1",
                            replacement = "ALM71C9-004-1")
fixed_ID <- str_replace_all(fixed_ID, pattern = SPACE,
                            replacement = "")
fixed_ID <- str_replace_all(fixed_ID, pattern = "\\(A\\)",
                            replacement = "A")
# IMPORTANT NOTE: this typo happened during sorting so it is screwed up!
fixed_ID <- str_replace_all(fixed_ID, pattern = "LM86C6-005-2",
                            replacement = "LM68C6-005-2")

resp.Group <- c()
for (e in seq_along(fixed_ID)){
    if (fixed_ID[e] %in% Resp == TRUE) {
        resp.Group[e] <- c("Responder")
    } else if (fixed_ID[e] %in% nonResp == TRUE) {
        resp.Group[e] <- c("Nonresponder")
    } else if (fixed_ID[e] %in% Unkn == TRUE) {
        resp.Group[e] <- c("Unknown")
    } else {
        resp.Group[e] <- NA
    }
}
summary(is.na(resp.Group))


# Get the matching dates
# NOTE: some samples have (A) in their name which causes issues with the pattern matching,
# as brackets are special characters which has to be escaped. Therefore the following loop
# will produce NAs where the bracketed sample names are located.
# This split loop is much more efficient at doing the same thing as the main loop
dates_df <- data.frame(matrix(ncol = 2, nrow = length(iSorts)))
colnames(dates_df) <- c("Sample", "Date")
for (e in seq_along(iSorts)) {
    dates_df[e, 1] <- keyword(iSorts[[e]])$"$FIL" # this is a special inherent function to FF objects
    dates_df[e, 2] <- str_extract(names(iSorts)[e], START %R% one_or_more(ALNUM) %R%
                                      "-" %R% one_or_more(ALNUM) %R%
                                      "-" %R% one_or_more(ALNUM))
}


# This loop will match the dates of the samples
Dates <- vector(mode = "character", length = nrow(iSort_df))
for (e in seq_along(iSort_df$Name)) {
    Dates[e] <- dates_df$Date[iSort_df$Name[e] == dates_df$Sample]
}
# NOTE: don't use grep, grepl or str_detect, because for some reason it does not work well here


# I believe I have the necessary metadata to unify everything in the cleaned up
# database
metadata <- list(sampleID, fixed_ID, resp.Group, Dates)

iSort_unif <- mutate(iSort_df, dates = metadata[[4]], resp.Group = metadata[[3]], oSampleID = metadata[[1]],
                     fSampleID = metadata[[2]], .before = "EpCAM") #oSampleID = original sample ID, fSampleID = fixed sample ID


# I will save the generated dataframe
write_csv(x = iSort_unif, file = paste0(script_version, "/", "Unified_index_sorts_with_metadata", "_", script_version, ".csv"))


# Create a dataframe containing the sample specific CTC counts, IDs and respGroups
count_CTCs = function(dataframe) {
    # I will bind the script version to a local variable to stop VS Codium from complaining :)
    script_version <- script_version

    # This line will subset the unique sampleIDs
    uniqueIDs <- unique(dataframe$fSampleID)
    
    # This line will create a new df containing the unique sampleIDs, pure patient IDs the corresponding CTC numbers and the response group
    cellNumbers <- data.frame(fSampleID = uniqueIDs, patientID = NA, cellNumbers = NA, responseGroup = NA)
    
    # I will create a temp df to hold the unique sample info for the new df, and run a loop to fill up the cellNumbers df
    filteredSample <- c()
    for (i in seq_along(uniqueIDs)) {
        
        filteredSample <- dplyr::filter(dataframe, fSampleID == uniqueIDs[i])
        
        cellNumbers$cellNumbers[i] <- nrow(filteredSample)
        cellNumbers$responseGroup [i] <- unique(filteredSample$resp.Group)
    }

    # This line will convert the sampleIDs to pure patient IDs (without the treatment cycle part)
    # I will import the %R% operator from the rebus package to this function so VSCodium will recognize it
    #' @importFrom rebus %R%
    cellNumbers$patientID <- stringr::str_remove_all(uniqueIDs, rebus::one_or_more("-") %R% rebus::one_or_more(ALNUM) %R% END)
    
    readr::write_csv(x = cellNumbers, file = paste0(script_version, "/", "CTC_numbers_per_patient", "_", script_version, ".csv"))
    message("The CTC numbers were saved as a .csv under the following name:", "\n", "CTC_numbers_per_patient", "_", script_version, ".csv")
    
    message("\n", "Here are the CTC numbers and additional information:")
    print(cellNumbers)
    return(cellNumbers)
    
}
CTC_counts <- count_CTCs(iSort_unif)


# Create a dataframe containing the sum of CTCs/patient, IDs and respGroups
summ_counts = function(dataframe) {
  ## Define variables
  
  unique_IDs <- vector(mode = "character")
  output_df <- data.frame(matrix(nrow = length(unique_IDs), ncol = 3))
  
  ## Calculate the sum of CTCs per patient
  unique_IDs <- unique(dataframe$patientID)
  for (i in seq_along(unique_IDs)) {
    output_df[i, 1] <- unique_IDs[i]
    output_df[i, 2] <- sum(dataframe[grep(x = dataframe$patientID, pattern = unique_IDs[i]), ][, 3])
    output_df[i, 3] <- unique(dataframe[grep(x = dataframe$patientID, pattern = unique_IDs[i]), ][, 4])
  }
  
  ## Manage the output dataframe
  colnames(output_df) <- c("patient_ID", "cell_number", "response_group")
  
  ## Save and print the output
  readr::write_csv(x = output_df, file = paste0(script_version, "/", "aggregated_CTC_numbers_per_patient", "_", script_version, ".csv"))
  message("The CTC numbers were saved as a .csv under the following name:", "\n", "aggregated_CTC_numbers_per_patient", "_", script_version, ".csv")
  message("\n", "Here are the CTC numbers and additional information:")
  print(output_df)
  
  ## Return the output
  return(output_df)
  
}
sum_CTC_counts <- summ_counts(CTC_counts)


# Remove unnecessary variables from this section
rm(dates_df, CTC_counts, sum_CTC_counts, fixed_ID, metadata, sampleID, resp.Group, Dates, Resp, nonResp, Unkn)

#################################################################   End Section  ##################################################################






#################################################################  Finalize the control and sample FlowFrames  ##################################################################
                                                                #          and do data normalization           #




## Load the saved unified dataframe if not present in the .GlobalEnv
if (!exists("niSortU_RvNR", envir = .GlobalEnv)) {

    if (file.exists(paste0(script_version, "/", "RespvNonresp_NUnif_index_sorts_with_metadata", "_", script_version, ".csv"))) {

      niSortU_RvNR <- read_csv(paste0(script_version, "/", "RespvNonresp_NUnif_index_sorts_with_metadata", "_", script_version, ".csv"))
      niSortU_RvNR <- as.data.frame(niSortU_RvNR)
        message("The niSortU_RvNR dataframe was successfully loaded into the .GlobalEnv.")
    } else {

        message("The niSortU_RvNR dataframe is not loaded in the .GlobalEnvironment and the save file: \n",
        paste0("RespvNonresp_NUnif_index_sorts_with_metadata", "_", script_version, ".csv"), "\n",
        "was not found in the current script version directory. Please run the code block below to generate it.")

    } 

} else {
  
  message("niSortU_RvNR is already in the .Global.Env, no further actions are necessary.")
  
}




## Prep the ctrlMean_df for the normalization by adding the sort dates to the
## df in a separate column
ctrlDates <- str_extract(ctrlMedian_df$Names, START %R% one_or_more(ASCII_ALNUM) %R% "-" 
                        %R% one_or_more(ASCII_ALNUM) %R% "-" %R% one_or_more(ASCII_ALNUM))
ctrlMedian_df <- mutate(ctrlMedian_df, Dates = ctrlDates)




## Run the normalization based on date matching 
norm_iSortu <- iSort_unif
for (i in seq_along(ctrlMedian_df$Dates)) {
  tmp <- grep(ctrlMedian_df$Dates[i], iSort_unif$dates)
  for (e in seq_along(tmp)) {
    norm_iSortu[tmp[e], 6:9] <- iSort_unif[tmp[e], 6:9]/ctrlMedian_df[i, 8:11]
  }
  
}
head(norm_iSortu)
norm_iSortu <- norm_iSortu[, -c(7, 8)]
head(norm_iSortu)


# I will save the generated  normalized dataframe
write_csv(x = norm_iSortu, file = paste0(script_version, "/", "Normalized_unif_index_sorts_with_metadata", "_", script_version, ".csv"))




## As I'm mostly interested in the resp vs nonresp, I will create a df where the unknown
## resp.Group is removed
niSortU_RvNR <- norm_iSortu[!grepl("Unknown", norm_iSortu$resp.Group), ]
head(niSortU_RvNR)


# I will save the generated clean normalized dataframe
write_csv(x = niSortU_RvNR, file = paste0(script_version, "/", "RespvNonresp_NUnif_index_sorts_with_metadata", "_", script_version, ".csv"))


# Removing unnecessary variables from this section
rm(ctrlMedian_df, ctrlDates, ctrlFFs, e, i, iDates, sCtrl, tmp)

#################################################################   End Section  ##################################################################






#################################################################  Normalized data visualization  ##################################################################
                                                                #          and filtering          #



## Visualize some parameters with the normalized samples
niSortU_RvNR <- dplyr::select(niSortU_RvNR,
                               fSampleID, resp.Group, EpCAM, PSMA)
long_niSortU_RvNR <- pivot_longer(niSortU_RvNR, cols = c(EpCAM, PSMA),
                               names_to = "Channels", values_to = "MFI_values")


# Scatter plot of the normalized EpCAM and PSMA channel values
QC_scatter_p <- ggplot(niSortU_RvNR, 
                       aes(x = log10(EpCAM), y = log10(PSMA),
                           fill = resp.Group, color = resp.Group)) +
    geom_point(size = 3, alpha = 0.5) +
    scale_fill_tron() +
    scale_color_tron() +
    geom_hline(yintercept = 2, color = "red", linetype = "dashed") +
    geom_vline(xintercept = 2, color = "red", linetype = "dashed") +
    geom_hline(yintercept = -1.5, color = "orange", linetype = "dashed") +
    geom_vline(xintercept = -1, color = "orange", linetype = "dashed") +
    labs(x = expression("log"[10]*" EpCAM MFI"), y = expression("log"[10]*" PSMA MFI"),
         fill = expression("Response groups"), color = expression("Response groups")) +
    #ggtitle("EpCAM vs PSMA MFI values - gating strategy for analysis") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(QC_scatter_p)

ggsave(filename = paste0("QC_gating_scatterplot", "_", script_version, ".png"), QC_scatter_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(QC_scatter_p)




## Data QC and filtering 


# If I want to use the log transformed values for any analysis
# I will have to filter out the PSMA and EpCAM values below 0
# or I have to input low MFI values instead
# This function will impute the mean of 2 randomly selected MFI values from the top10 lowest MFI values
imput_low_MFI_value = function(data, channel1, channel2, replacement_for_random_sampling = FALSE) {
    # selecting the above 0 MFI values from the dataset based on the two channels
    channel1_values_above_zero <- subset(data[, channel1], data[, channel1] > 0)
    channel2_values_above_zero <- subset(data[, channel2], data[, channel2] > 0)
    
    # taking the 10 lowest above 0 MFI values from each channels
    top10_low_channel1_values <- head(channel1_values_above_zero[order(channel1_values_above_zero)], 10)
    top10_low_channel2_values <- head(channel2_values_above_zero[order(channel2_values_above_zero)], 10)
    
    #imputing values into channels 1 and 2
    imputed_df <- data
    for (e in seq_len(nrow(imputed_df))) {
        #the if statement for channel 1
        if (imputed_df[e, channel1] <= 0) {
            imputed_df[e, channel1] <- mean(sample(top10_low_channel1_values,
                                                   size = 2,
                                                   replace = replacement_for_random_sampling))
            
        } else {
            #do nothing
        }
        
        # The if statement for channel 2
        if (imputed_df[e, channel2] <= 0) {
            imputed_df[e, channel2] <- mean(sample(top10_low_channel2_values,
                                                     size = 2,
                                                     replace = replacement_for_random_sampling))
        } else {
            #do nothing
        }
        
    }
    
    
    # creating a name for the new object based on the name of the original dataframe
    old_object_name <- deparse(substitute(data))
    new_object_name <- paste0(old_object_name, "_imputed")
    
    # creating the new object
    assign(new_object_name, imputed_df, envir = .GlobalEnv)
    message(paste("New values were imputed into:", old_object_name))
    message(paste("A new imputed dataframe was created withe the name:", new_object_name))

}

imput_low_MFI_value(niSortU_RvNR, "EpCAM", "PSMA")
head(niSortU_RvNR_imputed)
min(niSortU_RvNR_imputed$EpCAM)    
min(niSortU_RvNR_imputed$PSMA)


# Filter the data to remove extremely low and high values
niSortFilt <- niSortU_RvNR_imputed[log10(niSortU_RvNR_imputed$EpCAM) >= -1 | log10(niSortU_RvNR_imputed$PSMA) >= -1.5, ]
niSortFilt <- niSortFilt[log10(niSortFilt$EpCAM) < 2 & log10(niSortFilt$PSMA) < 2, ]
niSortFilt <- niSortFilt[!log10(niSortFilt$EpCAM) < -5, ]
niSortFilt <- niSortFilt[!log10(niSortFilt$PSMA) < -5, ]


# Add patient IDs to the normalized, imputed, filtered dataset
niSortFilt <- mutate(niSortFilt,
    patientID = stringr::str_extract(niSortFilt$fSampleID, pattern = "[A-Z0-9]+[-]+[0-9]+"), .after = fSampleID)
head(niSortFilt)





## Create new dataframes for some summary stats after filtering

# This function creates and saves a new cell number dataframe
count_filtered_cell_numbers = function(filtered_data_frame) {
    # define the global script version as a local variable
    script_version <- script_version
    
    # create a sample ID string
    sampleIDs <- filtered_data_frame$fSampleID
    uniqueSampleIDs <- unique(sampleIDs)

    # create a new cell numbers dataframe
    filtered_cell_numbers <- data.frame(matrix(ncol = 3, nrow = length(uniqueSampleIDs)))
    colnames(filtered_cell_numbers) <- c("sampleID", "responseGroup", "totalCellNumber")

    # fill the dataframe using the for loop
    for (i in seq_along(uniqueSampleIDs)) {
        
        filtered_cell_numbers[i, 1] <- uniqueSampleIDs[i]
        filtered_cell_numbers[i, 2] <- unique(dplyr::filter(filtered_data_frame, fSampleID == uniqueSampleIDs[i])$resp.Group)
        filtered_cell_numbers[i, 3] <- nrow(dplyr::filter(filtered_data_frame, fSampleID == uniqueSampleIDs[i]))
        
    }

    # save the resulting dataframe
    readr::write_csv(x = filtered_cell_numbers,
                    file = paste0(script_version, "/", "Filtered_data_cell_numbers", "_", script_version, ".csv"))
    message("The new cell numbers dataframe was saved to the current script directory: \n",
    paste0(getwd(), "/", script_version), "as: ", "\n",
    paste0("Filtered_data_cell_numbers", "_", script_version, ".csv"))

    # return the resulting dataframe for assignment
    return(filtered_cell_numbers)

}

filtered_cell_numbers_df <- count_filtered_cell_numbers(niSortFilt)
head(filtered_cell_numbers_df)

# This function creates and saves a new general stats dataframe
create_filtered_data_stats = function(cell_numbers_df) {
    # define the script version as a local variable
    script_version <- script_version
    
    # create a new stat dataframe
    filtered_data_stats <- data.frame(matrix(nrow = 3, ncol = 4))
    rownames(filtered_data_stats) <- c("Responder", "Nonresponder", "All_samples")
    colnames(filtered_data_stats) <- c("totalSampleNumber", "totalCellNumber", "mean", "median")

    # fill the stat dataframe with the relevant data
    filtered_data_stats[1, 1] <- nrow(dplyr::filter(cell_numbers_df, responseGroup == "Responder"))
    filtered_data_stats[1, 2] <- sum(dplyr::filter(cell_numbers_df, responseGroup == "Responder")$totalCellNumber)
    filtered_data_stats[1, 3] <- mean(dplyr::filter(cell_numbers_df, responseGroup == "Responder")$totalCellNumber)
    filtered_data_stats[1, 4] <- median(dplyr::filter(cell_numbers_df, responseGroup == "Responder")$totalCellNumber)
    filtered_data_stats[2, 1] <- nrow(dplyr::filter(cell_numbers_df, responseGroup == "Nonresponder"))
    filtered_data_stats[2, 2] <- sum(dplyr::filter(cell_numbers_df, responseGroup == "Nonresponder")$totalCellNumber)
    filtered_data_stats[2, 3] <- mean(dplyr::filter(cell_numbers_df, responseGroup == "Nonresponder")$totalCellNumber)
    filtered_data_stats[2, 4] <- median(dplyr::filter(cell_numbers_df, responseGroup == "Nonresponder")$totalCellNumber)
    filtered_data_stats[3, 1] <- nrow(cell_numbers_df)
    filtered_data_stats[3, 2] <- sum(cell_numbers_df$totalCellNumber)
    filtered_data_stats[3, 3] <- mean(cell_numbers_df$totalCellNumber)
    filtered_data_stats[3, 4] <- median(cell_numbers_df$totalCellNumber)
    
    # save the resulting stat dataframe 
    readr::write_csv(x = filtered_data_stats,
        file = paste0(script_version, "/", "General_cell_number_stats_after_data_filtering", "_", script_version, ".csv"))
    message("The created new stat dataframe was saved at the following location: \n",
    paste0(getwd(), "/", script_version),"as: \n",
    paste0("General_cell_number_stats_after_data_filtering", "_", script_version, ".csv"))

    # return the created stat dataframe
    return(filtered_data_stats)

}

filtered_cell_number_stats_df <- create_filtered_data_stats(filtered_cell_numbers_df)
head(filtered_cell_number_stats_df)

# Calculate and save the p-values for the cell number means per resp. group
total_sorted_RvNR_means_p_val <- t.test(dplyr::filter(filtered_cell_numbers_df, responseGroup == "Responder")$totalCellNumber,
    dplyr::filter(filtered_cell_numbers_df, responseGroup == "Nonresponder")$totalCellNumber)

p_val_to_save <- capture.output(print(total_sorted_RvNR_means_p_val))
readr::write_lines(x = p_val_to_save,
    file = paste0(script_version, "/", "total_sorted_RvNR_means_p-value", "_", script_version, ".txt"))

# Calculate and save the p-values for the total cell numbers per resp. group
total_sorted_RvNR_p_val <- wilcox.test(dplyr::filter(filtered_cell_numbers_df, responseGroup == "Responder")$totalCellNumber,
    dplyr::filter(filtered_cell_numbers_df, responseGroup == "Nonresponder")$totalCellNumber)

p_val_to_save <- capture.output(print(total_sorted_RvNR_p_val))
readr::write_lines(x = p_val_to_save,
    file = paste0(script_version, "/", "total_sorted_RvNR_p-value", "_", script_version, ".txt"))




## Visualize various parameters of the filtered data


# EpCAM vs PSMA scatter plot following filtering
after_QC_scatter_p <- ggplot(niSortFilt, 
                       aes(x = log10(EpCAM), y = log10(PSMA),
                           fill = resp.Group, color = resp.Group)) +
    geom_point(size = 3, alpha = 0.5) +
    scale_fill_tron(labels = c("Nonresponder", "Responder")) +
    scale_color_tron(labels = c("Nonresponder", "Responder")) +
    geom_hline(yintercept = 2, color = "red", linetype = "dashed") +
    geom_vline(xintercept = 2, color = "red", linetype = "dashed") +
    geom_hline(yintercept = -1.5, color = "orange", linetype = "dashed") +
    geom_vline(xintercept = -1, color = "orange", linetype = "dashed") +
    labs(x = expression("log"[10]*" EpCAM MFI"), y = expression("log"[10]*" PSMA MFI"),
         fill = expression("Response groups"), color = expression("Response groups")) +
    #ggtitle("EpCAM vs PSMA gating strategy") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))   
print(after_QC_scatter_p)

ggsave(filename = paste0("After_QC_gating_scatter_plot", "_", script_version, ".png"), after_QC_scatter_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(after_QC_scatter_p)


# Bar chart showing the total cell numbers sorted / response group
nResp <- nrow(dplyr::filter(niSortFilt, resp.Group == "Responder"))
nNonresp <- nrow(dplyr::filter(niSortFilt, resp.Group == "Nonresponder"))

total_cell_number_df <- data.frame("nResp" = nResp, "nNonresp" = nNonresp)
print(total_cell_number_df)

long_total_cell_number_df <-  pivot_longer(total_cell_number_df, cols = c(nResp, nNonresp),
                                           names_to = "resp.Group", values_to = "Cell_numbers")


raw_cell_numbers_per_resp_group_p <- ggplot(data = long_total_cell_number_df,
                                        aes(x = resp.Group, y = Cell_numbers, fill = resp.Group)) +
    geom_col(show.legend = FALSE, position = "dodge", color = "black") +
    geom_text(aes(label = Cell_numbers), vjust = -0.5, position = position_dodge(width = 0.9)) +
    #geom_pwc(method = "wilcox_test", label = "p = {p}") + #this geom is for pairwise comparison and is amazing!
    #scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, 10)) +
    scale_x_discrete(labels = c("Nonresponder", "Responder")) +
    scale_fill_tron(labels = c("Nonresponder", "Responder")) +
    #scale_fill_manual(values = c("#00BFC4", "#F8766D"), labels = c("Responder", "Nonresponder")) +
    labs(fill = "Response groups", x = expression("Response groups"), y = expression("Cell numbers")) +
    #ggtitle("Total cell numbers/response groups") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
          strip.text = element_blank(), text = element_text(size = 24))
print(raw_cell_numbers_per_resp_group_p)

ggsave(filename = paste0("Total_raw_cell_numbers_per_response_groups", "_", script_version, ".png"),
       raw_cell_numbers_per_resp_group_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(raw_cell_numbers_per_resp_group_p)
rm(nResp, nNonresp ,total_cell_number_df, long_total_cell_number_df)


# Bar charts showing the mean cell numbers per response groups with error bars
summary_cell_number_stats_p <- ggplot(data = filtered_cell_numbers_df, aes(x = responseGroup, y = totalCellNumber, fill = responseGroup)) +
                                stat_summary(fun = "mean", geom = "col", position = "dodge", color = "black") +
                                stat_summary(fun.data = "mean_se", geom = "errorbar", position = "dodge", width = 0.25, color = "black") +
                                geom_pwc(method = "t_test", label = "p = {p}", y.position = 14, label.size = 6) + #this geom is for pairwise comparison and is amazing!
                                #scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 10)) +
                                scale_x_discrete(labels = c("Nonresponder", "Responder")) +
                                scale_fill_tron(labels = c("Nonresponder", "Responder")) +
                                #geom_text(aes(label = round(after_stat(y), 2)), stat = "summary", fun = "mean", nudge_y = 3) +
                                labs(x = "Response groups", y = "Mean Cell numbers", fill = "Response groups") +
                                #ggtitle("Mean cell numbers per response group") +
                                theme_classic() +
                                theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                                                                    strip.text = element_blank(), text = element_text(size = 24))
print(summary_cell_number_stats_p)

ggsave(filename = paste0("Total_raw_cell_numbers_per_response_groups_summary_stat", "_", script_version, ".png"),
       summary_cell_number_stats_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(summary_cell_number_stats_p)


# Box plot showing the median cell numbers per response groups with ranges and a wilcox test p value
# NOTE: this is more appropriate to use in my opinion as the sample distribution seems very much skewed and not normal
summary_cell_number_stat_median_p <- ggplot(data = filtered_cell_numbers_df, aes(x = responseGroup, y = totalCellNumber, fill = responseGroup)) +
                                    geom_boxplot() +
                                    geom_dotplot(stackdir = "center", binaxis = "y", binwidth = 1, alpha = 0.5) +
                                    geom_pwc(method = "wilcox_test", label = "p = {p}", y.position = 45, label.size = 6) + #this geom is for pairwise comparison and is amazing!
                                    scale_x_discrete(labels = c("Nonresponder", "Responder")) +
                                    scale_fill_tron(labels = c("Nonresponder", "Responder")) +
                                    labs(x = "Response groups", y = "Median Cell numbers", fill = "Response groups") +
                                    #ggtitle("Median cell numbers per response group") +
                                    theme_classic() +
                                    theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                                                                        strip.text = element_blank(), text = element_text(size = 24))
print(summary_cell_number_stat_median_p)

ggsave(filename = paste0("Total_raw_cell_numbers_per_response_groups_summary_stat_MEDIANS", "_", script_version, ".png"),
       summary_cell_number_stat_median_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(summary_cell_number_stat_median_p)





# Based on the minimum thresholds set here I will additionally label each sample to be single EpCAM positive (SE)
# single PSMA positive (SP) or double positive (DP)
niSortFilt$EvP_status <- NA
for (e in seq_len(nrow(niSortFilt))) {
    if (log10(niSortFilt$EpCAM[e]) < -1) {
        niSortFilt$EvP_status[e] <- "SP"
    } else if (log10(niSortFilt$PSMA[e]) < -1.5) {
        niSortFilt$EvP_status[e] <- "SE"
    } else {
        niSortFilt$EvP_status[e] <- "DP"
    }
}


# Creating additional palettes for plotting (as the reponder, non-responder palette should not be used)
plotting_colors_st <- pal_startrek(palette = "uniform")(7)
plotting_colors_tron <- pal_tron(palette = "legacy")(7)


# Creating a plot to check if the EpCAM vs PSMA status assignment was successful
QC_EvP_status_p <- ggplot(niSortFilt,
                          aes(x = log10(EpCAM), y = log10(PSMA),
                              fill = EvP_status, color = EvP_status)) +
    geom_point(size = 3, alpha = 0.5) +
    #geom_hline(yintercept = 2, color = "red", linetype = "dashed") +
    #geom_vline(xintercept = 2, color = "red", linetype = "dashed") +
    geom_hline(yintercept = -1.5, color = "orange", linetype = "dashed") +
    geom_vline(xintercept = -1, color = "orange", linetype = "dashed") +
    scale_fill_manual(labels = c("EpCAM+PSMA+", "EpCAM+PSMA-", "EpCAM-PSMA+"), values = plotting_colors_tron[c(7, 4, 6)]) +
    scale_color_manual(labels = c("EpCAM+PSMA+", "EpCAM+PSMA-", "EpCAM-PSMA+"), values = plotting_colors_tron[c(7, 4, 6)]) +
    labs(x = expression("log"[10] *" EpCAM MFI"), y = expression("log"[10]*" PSMA MFI"),
         fill = expression("Marker status"), color = expression("Marker status")) +
    #ggtitle("EpCAM vs PSMA gating strategy - marker distribution") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(QC_EvP_status_p)

ggsave(filename = paste0("Marker_status_scatter_plot", "_", script_version, ".png"), QC_EvP_status_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(QC_EvP_status_p)


# Save the normalized, imputed and filtered dataframe
write_csv(x = niSortFilt, file = paste0(script_version, "/", "Total_norm_input_filt_indexDataframe", "_", script_version, ".csv"))


# Remove unnecessary variables from this section
rm(p_val_to_save, total_sorted_RvNR_p_val, total_sorted_RvNR_means_p_val)

################################################################   End Section  ##################################################################






#################################################################  Marker distribution analysis  ##################################################################
                                                                #        on both channels        #




## A quick exploration of the surface marker distribution


# I will transform the filtered dataframe to a long format
long_niSortFilt <- pivot_longer(niSortFilt, cols = c(EpCAM, PSMA),
                                names_to = "Channels", values_to = "MFI_values")
head(long_niSortFilt)


# First, visualize the filtered dataframe in terms of total EpCAM and PSMA signal per response group
EpCAM_and_PSMA_RvNR_p <- ggplot(data = long_niSortFilt,
                                aes(fill = resp.Group, x = Channels, 
                                    y = log10(MFI_values))) +
    #geom_boxplot() +
    geom_dotplot(stackdir = "center", binaxis = "y",
                 binwidth = 0.05, alpha = 0.5) +
    scale_fill_tron() +
    facet_wrap(vars(resp.Group), strip.position = "bottom") +
    labs(x = expression(Channels), y = expression("log"[10]*" MFI values"),
         fill = expression("Response groups")) +
    ggtitle("Total MFI values per response groups") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
          strip.text = element_blank(), text = element_text(size = 24))
print(EpCAM_and_PSMA_RvNR_p)

ggsave(filename = paste0("Total_MFI_values_per_response_group_facet_wrap", "_", script_version, ".png"), EpCAM_and_PSMA_RvNR_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(EpCAM_and_PSMA_RvNR_p)




## Creating summary stats for for the marker distribution data


# I will look at some basic stats for the response groups regarding the total cell numbers
surface_marker_stats = function(filtered_dataframe) {
    total_cell_number <- nrow(filtered_dataframe)
    print(total_cell_number)


    # Working on the responders
    niSortFilt_resp <- dplyr::filter(filtered_dataframe, resp.Group == "Responder")

    resp_total_cell_number <- nrow(niSortFilt_resp)
    print(resp_total_cell_number)

    resp_single_positive_EpCAM <- length(grep("SE", niSortFilt_resp$EvP_status))
    print(resp_single_positive_EpCAM)

    resp_single_positive_PSMA <- length(grep("SP", niSortFilt_resp$EvP_status))
    print(resp_single_positive_PSMA)

    resp_double_positive <- length(grep("DP", niSortFilt_resp$EvP_status))
    print(resp_double_positive)

    resp_cell_number_df <- data.frame("Total_cell_number" = resp_total_cell_number, "Double_positives" = resp_double_positive, "EpCAM_positives" = resp_single_positive_EpCAM,
                                    "PSMA_positives" = resp_single_positive_PSMA)
    rownames(resp_cell_number_df) <- c("Responder")


    # Working on the Nonresponders
    niSortFilt_nonresp <- dplyr::filter(filtered_dataframe, resp.Group == "Nonresponder")

    nonresp_total_cell_number <- nrow(niSortFilt_nonresp)
    print(nonresp_total_cell_number)

    nonresp_single_positive_EpCAM <- length(grep("SE", niSortFilt_nonresp$EvP_status))
    print(nonresp_single_positive_EpCAM)

    nonresp_single_positive_PSMA <- length(grep("SP", niSortFilt_nonresp$EvP_status))
    print(nonresp_single_positive_PSMA)

    nonresp_double_positive <- length(grep("DP", niSortFilt_nonresp$EvP_status))
    print(nonresp_double_positive)

    nonresp_cell_number_df <- data.frame("Total_cell_number" = nonresp_total_cell_number, "Double_positives" = nonresp_double_positive, "EpCAM_positives" = nonresp_single_positive_EpCAM,
                                        "PSMA_positives" = nonresp_single_positive_PSMA)
    rownames(nonresp_cell_number_df) <- c("Nonresponder")


    # creating the unified summary stat df
    unif_cell_numbers_per_marker <- rbind(resp_cell_number_df, nonresp_cell_number_df)
    print(unif_cell_numbers_per_marker)


    # Return the summary stat df
    return(unif_cell_numbers_per_marker)
}
unif_cell_numbers_per_marker <- surface_marker_stats(niSortFilt)


# Looking at the total distribution of the response groups among the surface marker populations
unif_cell_numbers_per_marker_percent <- unif_cell_numbers_per_marker
unif_cell_numbers_per_marker_percent <- unif_cell_numbers_per_marker_percent / unif_cell_numbers_per_marker_percent$Total_cell_number * 100
unif_cell_numbers_per_marker_percent <- mutate(.data = unif_cell_numbers_per_marker_percent, resp.Group = c("Responder", "Nonresponder"))
print(unif_cell_numbers_per_marker_percent)


# Looking at the distribution of the response groups inside the surface marker populations
unif_cell_numbers_per_marker_inpop <- unif_cell_numbers_per_marker
for (ci in 2:4) {
    for (ri in 1:2) {
        unif_cell_numbers_per_marker_inpop[ri, ci] <- unif_cell_numbers_per_marker[ri, ci] / sum(unif_cell_numbers_per_marker[, ci]) * 100
    }
}
print(unif_cell_numbers_per_marker_inpop)
rm(ci, ri)

unif_cell_numbers_per_marker_inpop <- mutate(.data = unif_cell_numbers_per_marker_inpop, resp.Group = c("Responder", "Nonresponder"))
print(unif_cell_numbers_per_marker_inpop)



# Saving the summary data
write_csv(unif_cell_numbers_per_marker, file = paste0(script_version, "/", "Total_raw_cell_numbers_with_marker_dist", "_", script_version, ".csv"))
write_csv(unif_cell_numbers_per_marker_percent, file = paste0(script_version, "/", "Total_raw_cell_number_percent_with_marker_dist.csv", "_", script_version, ".csv"))
write_csv(unif_cell_numbers_per_marker_inpop, file = paste0(script_version, "/", "Cell_number_percent_dist_within_marker_groups.csv", "_", script_version, ".csv"))




## Visualizing the basic marker distribution data


# Visualizing the raw cell numbers based on marker distribution (thx bing!)
unif_cell_numbers_per_marker <- mutate(.data = unif_cell_numbers_per_marker, resp.Group = c("Responder", "Nonresponder"))
long_unif_cell_numbers_per_marker <- pivot_longer(unif_cell_numbers_per_marker, cols = c(Double_positives, EpCAM_positives, PSMA_positives),
                                                  names_to = "Marker_status", values_to = "Cell_numbers")


raw_cell_numbers_per_marker_p <- ggplot(data = long_unif_cell_numbers_per_marker,
                                        aes(x = Marker_status, y = Cell_numbers, fill = resp.Group)) +
    geom_col(position = "dodge", color = "black") +
    geom_text(aes(label = Cell_numbers), vjust = -0.5, position = position_dodge(width = 0.9)) +
    #geom_pwc(method = "wilcox_test", label = "p = {p}") + #this geom is for pairwise comparison and is amazing!
    #scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, 10)) +
    scale_x_discrete(labels = c("EpCAM+\nPSMA+", "EpCAM+\nPSMA-", "EpCAM-\nPSMA+")) +
    scale_fill_tron(labels = c("Nonresponder", "Responder")) +
    #scale_fill_manual(values = c("#00BFC4", "#F8766D"), labels = c("Responder", "Nonresponder")) +
    labs(fill = "Response groups", x = expression("Marker status"), y = expression("Cell numbers")) +
    #ggtitle("Total cell numbers/marker groups") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
          strip.text = element_blank(), text = element_text(size = 24))
print(raw_cell_numbers_per_marker_p)

ggsave(filename = paste0("Total_raw_cell_numbers_per_staining", "_", script_version, ".png"), raw_cell_numbers_per_marker_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(raw_cell_numbers_per_marker_p)
rm(long_unif_cell_numbers_per_marker)


# Visualizing the total distribution of the response groups among the surface marker populations (thx bing!)
long_unif_cell_numbers_per_marker_percent <- pivot_longer(unif_cell_numbers_per_marker_percent, cols = c(Double_positives, EpCAM_positives, PSMA_positives),
                                                          names_to = "Marker_status", values_to = "Percentages")

marker_percentages_p <- ggplot(data = long_unif_cell_numbers_per_marker_percent,
                               aes(x = Marker_status, y = Percentages, fill = resp.Group)) +
    geom_col(position = "dodge", color = "black") +
    #geom_pwc(method = "wilcox_test", label = "p = {p}") + #this geom is for pairwise comparison and is amazing!
    scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, 10)) +
    scale_x_discrete(labels = c("EpCAM+\nPSMA+", "EpCAM+\nPSMA-", "EpCAM-\nPSMA+")) +
    scale_fill_tron(labels = c("Nonresponder", "Responder")) +
    labs(fill = "Response groups", x = "Marker status", y = "Cell number percentages") +
    #scale_fill_manual(values = c("#00BFC4", "#F8766D"), labels = c("Responder", "Nonresponder")) +
    geom_text(aes(label = round(Percentages, 2)), vjust = -0.5, position = position_dodge2(0.9)) +
    #ggtitle("CTC distribution among marker groups") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))

print(marker_percentages_p)

ggsave(filename = paste0("Total_cell_numbers_percenteges_per_staining", "_", script_version, ".png"), marker_percentages_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(marker_percentages_p)
rm(long_unif_cell_numbers_per_marker_percent)


# Visualizing the distribution of the response groups inside the surface marker populations (thx bing!)
long_unif_cell_numbers_per_marker_inpop <- pivot_longer(unif_cell_numbers_per_marker_inpop, cols = c(Double_positives, EpCAM_positives, PSMA_positives),
                                                          names_to = "Marker_status", values_to = "Percentages")
print(long_unif_cell_numbers_per_marker_inpop)

marker_percentages_p <- ggplot(data = long_unif_cell_numbers_per_marker_inpop,
                               aes(x = Marker_status, y = Percentages, fill = resp.Group)) +
    geom_col(position = "dodge", color = "black") +
    #geom_pwc(method = "wilcox_test", label = "p = {p}") + #this geom is for pairwise comparison and is amazing!
    scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, 10)) +
    scale_x_discrete(labels = c("EpCAM+\nPSMA+", "EpCAM+\nPSMA-", "EpCAM-\nPSMA+")) +
    scale_fill_tron(labels = c("Nonresponder", "Responder")) +
    labs(fill = "Response groups", x = "Marker status", y = "Cell number percentages") +
    #scale_fill_manual(values = c("#00BFC4", "#F8766D"), labels = c("Responder", "Nonresponder")) +
    geom_text(aes(label = round(Percentages, 2)), vjust = -0.5, position = position_dodge2(0.9)) +
    #ggtitle("CTC distribution inside marker groups") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(marker_percentages_p)

ggsave(filename = paste0("Response_group_percentage_distrib_per_surface_marker_pop", "_", script_version, ".png"), marker_percentages_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(marker_percentages_p)
rm(long_unif_cell_numbers_per_marker_inpop)



## Preparing the stats for a Chi-square goodness of fit test for the marker distribution data
## NOTE: the idea here is to see if the surface marker distribution differs between response groups
## in each surface marker population


# Contingency table for the EpCAM+PSMA+ population
# NOTE: the table will contain the count data and the expected count data
# if the distribution is 50-50, so there is no difference between the response groups
marker_cont_table_DP <- data.frame(observed_count = unif_cell_numbers_per_marker[, 2],
                                expected_count = rbind(c(sum(unif_cell_numbers_per_marker[, 2]) * 0.5), c(sum(unif_cell_numbers_per_marker[, 2]) * 0.5)))


# Contingency table for the EpCAM+PSMA- population
# NOTE: the table will contain the count data and the expected count data
# if the distribution is 50-50, so there is no difference between the response groups
marker_cont_table_E <- data.frame(observed_count = unif_cell_numbers_per_marker[, 3],
                                expected_count = rbind(c(sum(unif_cell_numbers_per_marker[, 3]) * 0.5), c(sum(unif_cell_numbers_per_marker[, 3]) * 0.5)))


# Contingency table for the EpCAM-PSMA+ population
# NOTE: the table will contain the count data and the expected count data
# if the distribution is 50-50, so there is no difference between the response groups
marker_cont_table_P <- data.frame(observed_count = unif_cell_numbers_per_marker[, 4],
                                expected_count = rbind(c(sum(unif_cell_numbers_per_marker[, 4]) * 0.5), c(sum(unif_cell_numbers_per_marker[, 4]) * 0.5)))




## Conducting a Chi-square goodness of fit test 


# Chi-square for the EpCAM+PSMA+ population and saving the result
chi_sqr_test_DP <- chisq.test(marker_cont_table_DP)
chi_sqr_test_DP_p_val <- capture.output(print(chi_sqr_test_DP))
write_lines(x = chi_sqr_test_DP_p_val,
            file = paste0(script_version, "/", "Chi_sqr_goodness_of_fit_EpCAM+PSMA+_marker_distribution_p-value", "_", script_version, ".txt"))

# Chi-square for the EpCAM+PSMA- population
chi_sqr_test_E <- chisq.test(marker_cont_table_E)
chi_sqr_test_E_p_val <- capture.output(print(chi_sqr_test_E))
write_lines(x = chi_sqr_test_E_p_val,
            file = paste0(script_version, "/", "Chi_sqr_goodness_of_fit_EpCAM+PSMA-_marker_distribution_p-value", "_", script_version, ".txt"))


# Chi-square for the EpCAM-PSMA+ population
chi_sqr_test_P <- chisq.test(marker_cont_table_P)
chi_sqr_test_P_p_val <- capture.output(print(chi_sqr_test_P))
write_lines(x = chi_sqr_test_P_p_val,
            file = paste0(script_version, "/", "Chi_sqr_goodness_of_fit_EpCAM-PSMA+_marker_distribution_p-value", "_", script_version, ".txt"))


# Visualizing the raw cell numbers based on marker distribution (thx bing!) with the corresponding p-values based on the chi-sqr test
unif_cell_numbers_per_marker <- mutate(.data = unif_cell_numbers_per_marker, resp.Group = c("Responder", "Nonresponder"))
long_unif_cell_numbers_per_marker <- pivot_longer(unif_cell_numbers_per_marker, cols = c(Double_positives, EpCAM_positives, PSMA_positives),
                                                  names_to = "Marker_status", values_to = "Cell_numbers")

raw_cell_numbers_per_marker_p <- ggplot(data = long_unif_cell_numbers_per_marker,
                                        aes(x = Marker_status, y = Cell_numbers, fill = resp.Group)) +
    geom_col(position = "dodge", color = "black") +
    geom_text(aes(label = Cell_numbers), vjust = -0.5, position = position_dodge(width = 0.9)) +
    geom_signif(y_position = rep(120, 3), xmin = 1:3-0.2, xmax = 1:3+0.2, 
    annotation = c(paste("p =", round(chi_sqr_test_DP$p.value, 4)), paste("p =", round(chi_sqr_test_E$p.value, 4)), paste("p =", round(chi_sqr_test_P$p.value, 4))),
    tip_length = 0.03) +
    #geom_pwc(method = "wilcox_test", label = "p = {p}") + #this geom is for pairwise comparison and is amazing!
    #scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, 10)) +
    scale_x_discrete(labels = c("EpCAM+\nPSMA+", "EpCAM+\nPSMA-", "EpCAM-\nPSMA+")) +
    scale_fill_tron(labels = c("Nonresponder", "Responder")) +
    #scale_fill_manual(values = c("#00BFC4", "#F8766D"), labels = c("Responder", "Nonresponder")) +
    labs(fill = "Response groups", x = expression("Marker status"), y = expression("Cell numbers")) +
    #ggtitle("Total cell numbers/marker groups") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
          strip.text = element_blank(), text = element_text(size = 24))
print(raw_cell_numbers_per_marker_p)

ggsave(filename = paste0("Total_raw_cell_numbers_per_staining_with_p-vals", "_", script_version, ".png"), raw_cell_numbers_per_marker_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(raw_cell_numbers_per_marker_p)
rm(long_unif_cell_numbers_per_marker)



# Remove unnecessary variables from this section
rm(chi_sqr_test_DP, chi_sqr_test_DP_p_val, chi_sqr_test_E, chi_sqr_test_E_p_val, chi_sqr_test_P, chi_sqr_test_P_p_val,
long_niSortFilt_EpCAM, long_niSortFilt_nonresp, long_niSortFilt_resp, long_niSortU_RvNR, marker_cont_table_DP,
marker_cont_table_E, marker_cont_table_P, unif_cell_numbers_per_marker, unif_cell_numbers_per_marker_percent, unif_cell_numbers_per_marker_inpop)

################################################################   End Section  ##################################################################






#################################################################  A more in-depth channel  ##################################################################
                                                                # and statistical analysis  #




## Analyzing the EpCAM channel


# I create an EpCAM only df
long_niSortFilt_EpCAM <- dplyr::filter(long_niSortFilt, Channels == "EpCAM")


#testing for normality and p values for each response groups to see if they are
#significantly different or not (thx GPT)
shapiro.test(dplyr::filter(long_niSortFilt_EpCAM,
                           long_niSortFilt_EpCAM$resp.Group == "Responder")$MFI_values) 
qqnorm(dplyr::filter(long_niSortFilt_EpCAM,
                     long_niSortFilt_EpCAM$resp.Group == "Responder")$MFI_values)
qqline(dplyr::filter(long_niSortFilt_EpCAM,
                     long_niSortFilt_EpCAM$resp.Group == "Responder")$MFI_values)

shapiro.test(dplyr::filter(long_niSortFilt_EpCAM,
                           long_niSortFilt_EpCAM$resp.Group == "Nonresponder")$MFI_values) 
qqnorm(dplyr::filter(long_niSortFilt_EpCAM,
                     long_niSortFilt_EpCAM$resp.Group == "Nonresponder")$MFI_values)
qqline(dplyr::filter(long_niSortFilt_EpCAM,
                     long_niSortFilt_EpCAM$resp.Group == "Nonresponder")$MFI_values)

# As the data is clearly not normally distributed (in theory) I can't use t-test,
# therefore I will do a log10 transformation (just like for the figures)
# which will give the data a more normal distribution)
# alternatively I can use a wilcox.test (Mann-Whitney U test)
ttest_NRvR_EpCAM <- t.test(log10(dplyr::filter(long_niSortFilt_EpCAM,
                                               long_niSortFilt_EpCAM$resp.Group == "Responder")$MFI_values),
                           log10(dplyr::filter(long_niSortFilt_EpCAM,
                                           long_niSortFilt_EpCAM$resp.Group == "Nonresponder")$MFI_values))
print(ttest_NRvR_EpCAM)
rm(ttest_NRvR_EpCAM)

MWUt_NRvR_EpCAM <- wilcox.test(dplyr::filter(long_niSortFilt_EpCAM,
                                            long_niSortFilt_EpCAM$resp.Group == "Responder")$MFI_values,
                          dplyr::filter(long_niSortFilt_EpCAM,
                                        long_niSortFilt_EpCAM$resp.Group == "Nonresponder")$MFI_values)
print(MWUt_NRvR_EpCAM)
rm(MWUt_NRvR_EpCAM)

# Visualize the EpCAM channel between response groups
long_niSortFilt_EpCAM <- dplyr::filter(long_niSortFilt, Channels == "EpCAM")

EpCAM_RvNR_p <- ggplot(long_niSortFilt_EpCAM,
                       aes(x = resp.Group, y = log10(MFI_values), fill = resp.Group)) +
    geom_boxplot(show.legend = FALSE) +
    geom_pwc(method = "wilcox_test", label = "p = {p}", label.size = 6) +
    scale_fill_tron() +
    #geom_segment(aes(x = "Nonresponder", xend = "Responder", y = 1.5, yend = 1.5),
    #             col = "black", linewidth = 1) +
    #stat_compare_means(label.y = 1.6, label.x = 1.4) +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5) +
    labs(x = "Response Groups", y = expression("EpCAM log"[10]*" MFI values"), fill = "Response Group") +
    #ggtitle("EpCAM MFI comparison between response groups") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(EpCAM_RvNR_p)

ggsave(filename = paste0("Total_EpCAM_MFI_Resp_v_Nonresp", "_", script_version, ".png"), EpCAM_RvNR_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(EpCAM_RvNR_p)
rm(long_niSortFilt_EpCAM)


# I create a PSMA only df
long_niSortFilt_PSMA <- dplyr::filter(long_niSortFilt, Channels == "PSMA")


#testing for normality and p values for each response groups to see if they are
#significantly different or not (thx GPT)
shapiro.test(dplyr::filter(long_niSortFilt_PSMA,
                           long_niSortFilt_PSMA$resp.Group == "Responder")$MFI_values) 
qqnorm(dplyr::filter(long_niSortFilt_PSMA,
                     long_niSortFilt_PSMA$resp.Group == "Responder")$MFI_values)
qqline(dplyr::filter(long_niSortFilt_PSMA,
                     long_niSortFilt_PSMA$resp.Group == "Responder")$MFI_values)

shapiro.test(dplyr::filter(long_niSortFilt_PSMA,
                           long_niSortFilt_PSMA$resp.Group == "Nonresponder")$MFI_values) 
qqnorm(dplyr::filter(long_niSortFilt_PSMA,
                     long_niSortFilt_PSMA$resp.Group == "Nonresponder")$MFI_values)
qqline(dplyr::filter(long_niSortFilt_PSMA,
                     long_niSortFilt_PSMA$resp.Group == "Nonresponder")$MFI_values)

# As the data is clearly not normally distributed (in theory) I can't use t-test,
# therefore I will do a log10 transformation (just like for the figures)
# which will give the data a more normal distribution)
# alternatively I can use a wilcox.test (Mann-Whitney U test)
ttest_NRvR_PSMA <- t.test(log10(dplyr::filter(long_niSortFilt_PSMA,
                                               long_niSortFilt_PSMA$resp.Group == "Responder")$MFI_values),
                           log10(dplyr::filter(long_niSortFilt_PSMA,
                                               long_niSortFilt_PSMA$resp.Group == "Nonresponder")$MFI_values))
print(ttest_NRvR_PSMA)

MWUt_NRvR_PSMA <- wilcox.test(dplyr::filter(long_niSortFilt_PSMA,
                                            long_niSortFilt_PSMA$resp.Group == "Responder")$MFI_values,
                              dplyr::filter(long_niSortFilt_PSMA,
                                            long_niSortFilt_PSMA$resp.Group == "Nonresponder")$MFI_values)
print(MWUt_NRvR_PSMA)


# Visualize the PSMA channel between response groups
long_niSortFilt_PSMA <- dplyr::filter(long_niSortFilt, Channels == "PSMA")

PSMA_RvNR_p <- ggplot(long_niSortFilt_PSMA,
                       aes(x = resp.Group, y = log10(MFI_values), fill = resp.Group)) +
    geom_boxplot(show.legend = FALSE) +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5) +
    geom_pwc(method = "wilcox_test", label = "p = {p}", label.size = 6) +
    scale_fill_tron() +
    #geom_segment(aes(x = "Nonresponder", xend = "Responder", y = 1.5, yend = 1.5),
    #             col = "black", linewidth = 1) +
    #stat_compare_means(label.y = 1.6, label.x = 1.4) +
    labs(x = "Response Groups", y = expression("PSMA log"[10]*" MFI values"), fill = "Response Group") +
    #ggtitle("PSMA MFI comparison between response groups") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(PSMA_RvNR_p)

ggsave(filename = paste0("Total_PSMA_MFI_Resp_v_Nonresp", "_", script_version, ".png"), PSMA_RvNR_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(PSMA_RvNR_p)

################################################################   End Section  ##################################################################






#################################################################  Surface marker analysis  ##################################################################
                                                                #    from cycle to cycle    #




## Exploring and visualizing the EpCAM and PSMA MFI changes over the cycles


# Preparing the full dataset for the cycle wise analysis (cCycle = collection cycle, tCycle = treatment cycle)
cycles <- mutate(.data = niSortFilt,
                      cCycle = str_extract(niSortFilt$fSampleID, one_or_more(ALNUM) %R% END),
                      tCycle = as.numeric((str_extract(niSortFilt$fSampleID, one_or_more(ALNUM) %R% END))) - 1)
long_cycles <- pivot_longer(cycles, cols = c(EpCAM, PSMA),
                                 names_to = "Channels", values_to = "MFI_values")

# Save the df with the added cycle data
readr::write_csv(x = cycles, file = paste0(path, "/", script_version, "/", "NormImpFilt_data_with_cycles_", script_version,  ".csv"))


# Looking at the EpCAM data first in separated response groups
# I will prepare the EpCAM data first, by filtering the long_cycles df
long_CEpCAM <- dplyr::filter(long_cycles, Channels == "EpCAM")
head(long_CEpCAM)


# I will prep the data and run this loop, which will calculate the median EpCAM values 
# for responders
long_CEpCAM_res <- dplyr::filter(.data = long_CEpCAM, resp.Group == "Responder")
head(long_CEpCAM_res)

tCycle_REpCAM_medians <- c()
for (i in min(long_CEpCAM_res$tCycle):max(long_CEpCAM_res$tCycle)) {
    tCycle_REpCAM_medians[i + 1] <- median(long_CEpCAM_res$MFI_values[long_CEpCAM_res$tCycle == i])
}
# Note: in the i + 1, the +1 is crucial, as the cycle numbers are starting with 0 so the vector indexing would also start with a 0, and that
# is not possible, so I loose an element if I don't ad the +1
print(tCycle_REpCAM_medians)


# Repeating the calculations for nonresponders
long_CEpCAM_nonres <- dplyr::filter(.data = long_CEpCAM, resp.Group == "Nonresponder")
head(long_CEpCAM_nonres)

tCycle_NREpCAM_medians <- c()
for (i in min(long_CEpCAM_nonres$tCycle):max(long_CEpCAM_nonres$tCycle)) {
    tCycle_NREpCAM_medians[i + 1] <- median(long_CEpCAM_nonres$MFI_values[long_CEpCAM_nonres$tCycle == i])
}
print(tCycle_NREpCAM_medians)


# I will also calculate the MAD values for each cycle
# for responders
tCycle_REpCAM_mad <- c()
for (i in min(long_CEpCAM_res$tCycle):max(long_CEpCAM_res$tCycle)) {
    tCycle_REpCAM_mad[i + 1] <- mad(long_CEpCAM_res$MFI_values[long_CEpCAM_res$tCycle == i])
}
tCycle_REpCAM_mad


# Repeat the MAD calculation for nonresponders
tCycle_NREpCAM_mad <- c()
for (i in min(long_CEpCAM_nonres$tCycle):max(long_CEpCAM_nonres$tCycle)) {
    tCycle_NREpCAM_mad[i + 1] <- mad(long_CEpCAM_nonres$MFI_values[long_CEpCAM_nonres$tCycle == i])
}
print(tCycle_NREpCAM_mad)


# This loop will calculate the number of cells for each treatment cycle
# I decided, not to look at treatment cycles which did not have at least 10 cells
# as mean or median for these are unreliable
# I will start with the responders
respCellperCyc <- c()
for (i in min(long_CEpCAM_res$tCycle):max(long_CEpCAM_res$tCycle)) {
  respCellperCyc[i + 1] <- length(long_CEpCAM_res$tCycle[long_CEpCAM_res$tCycle == i])
}
print(respCellperCyc)

# I will repeat the calculations for nonresponders
nonrespCellperCyc <- c()
for (i in min(long_CEpCAM_nonres$tCycle):max(long_CEpCAM_nonres$tCycle)) {
    nonrespCellperCyc[i + 1] <- length(long_CEpCAM_nonres$tCycle[long_CEpCAM_nonres$tCycle == i])
}
print(nonrespCellperCyc)


# I will create a new df with the EpCAM medians and cycle numbers
# for responders
tCycle_REpCAM_medians_df <- data.frame(matrix(nrow = length(tCycle_REpCAM_medians), ncol = 5))
colnames(tCycle_REpCAM_medians_df) <- c("tCycle", "cell_number", "EpCAM_medians", "EpCAM_MAD", "resp.Group")
tCycle_REpCAM_medians_df [, 1] <- c(min(long_CEpCAM_res$tCycle):max(long_CEpCAM_res$tCycle))
tCycle_REpCAM_medians_df [, 2] <- respCellperCyc
tCycle_REpCAM_medians_df [, 3] <- tCycle_REpCAM_medians
tCycle_REpCAM_medians_df [, 4] <- tCycle_REpCAM_mad
tCycle_REpCAM_medians_df [, 5] <- "Responder"
print(tCycle_REpCAM_medians_df)

# I will also create the new df for nonresponders
tCycle_NREpCAM_medians_df <- data.frame(matrix(nrow = length(tCycle_NREpCAM_medians), ncol = 5))
colnames(tCycle_NREpCAM_medians_df) <- c("tCycle", "cell_number", "EpCAM_medians", "EpCAM_MAD", "resp.Group")
tCycle_NREpCAM_medians_df [, 1] <- c(min(long_CEpCAM_nonres$tCycle):max(long_CEpCAM_nonres$tCycle))
tCycle_NREpCAM_medians_df [, 2] <- nonrespCellperCyc
tCycle_NREpCAM_medians_df [, 3] <- tCycle_NREpCAM_medians
tCycle_NREpCAM_medians_df [, 4] <- tCycle_NREpCAM_mad
tCycle_NREpCAM_medians_df [, 5] <- "Nonresponder"
print(tCycle_NREpCAM_medians_df)


# Unifying the EpCAM means and MADs into one df and saving it
tCycle_EpCAM_medians_unif <- rbind(tCycle_REpCAM_medians_df, tCycle_NREpCAM_medians_df)
write_csv(x = tCycle_EpCAM_medians_unif, file = paste0(script_version, "/", "Treatment_cycles_EpCAM_MFI_medians_MADs", "_", script_version, ".csv"))


# Plot the cycles data  for EpCAM responders (the responder df is filtered as the last cycle had less than 10 cells)
long_CEpCAM_res_filt <- dplyr::filter(.data = long_CEpCAM_res, tCycle != 4)

# Plot with pairwise comparison
tCycle_REpCAM_wilcox_p <- ggplot(data = long_CEpCAM_res_filt,
                                 aes(x = factor(tCycle), y = log10(MFI_values), fill = factor(tCycle), group = factor(tCycle))) +
    geom_pwc(method = "wilcox_test", label = "p = {p}", ref.group = "0", label.size = 6) + #this geom is for pairwise comparison and is amazing!
    geom_boxplot(show.legend = FALSE, position = "dodge") +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5,
                 position = "dodge") +
    scale_fill_manual(labels = c("EpCAM+PSMA+", "EpCAM+PSMA-", "EpCAM-PSMA+"), values = plotting_colors_tron[c(7, 3, 4, 6)]) +
    scale_color_manual(labels = c("EpCAM+PSMA+", "EpCAM+PSMA-", "EpCAM-PSMA+"), values = plotting_colors_tron[c(7, 3, 4, 6)]) +
    labs(x = "Treatment cycles", y = expression("EpCAM log"[10] * " MFI values"), fill = "Treatment cycles") +
    #ggtitle("EpCAM MFI comparison between treatment cycles - Responders") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(tCycle_REpCAM_wilcox_p)

ggsave(filename = paste0("EpCAM_MFI_over_treatment_cycles_Responders_wilcox", "_", script_version, ".png"), tCycle_REpCAM_wilcox_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(tCycle_REpCAM_wilcox_p)


# Plot the cycles data for EpCAM nonresponders

# Pairwise comparison
tCycle_NREpCAM_wilcox_p <- ggplot(data = long_CEpCAM_nonres,
                                 aes(x = factor(tCycle), y = log10(MFI_values), fill = factor(tCycle), group = factor(tCycle))) +
    geom_pwc(method = "wilcox_test", label = "p = {p}", ref.group = "0", label.size = 6) + #this geom is for pairwise comparison and is amazing!
    geom_boxplot(show.legend = FALSE, position = "dodge") +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5,
                 position = "dodge") +
    scale_fill_manual(labels = c("EpCAM+PSMA+", "EpCAM+PSMA-", "EpCAM-PSMA+"), values = plotting_colors_tron[c(7, 3, 4, 6)]) +
    scale_color_manual(labels = c("EpCAM+PSMA+", "EpCAM+PSMA-", "EpCAM-PSMA+"), values = plotting_colors_tron[c(7, 3, 4, 6)]) +
    labs(x = "Treatment cycles", y = expression("EpCAM log"[10] * " MFI values"), fill = "Treatment cycles") +
    #ggtitle("EpCAM MFI comparison between treatment cycles - Nonresponders") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(tCycle_NREpCAM_wilcox_p)

ggsave(filename = paste0("EpCAM_MFI_over_treatment_cycles_Nonresponders_wilcox", "_", script_version, ".png"), tCycle_NREpCAM_wilcox_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(tCycle_NREpCAM_wilcox_p)



# Looking at the PSMA data first in both response groups


# I will prepare the PSMA data first, by filtering the long_cycles df
long_CPSMA <- dplyr::filter(long_cycles, Channels == "PSMA")
head(long_CPSMA)

# I will prep the data and run this loop, which will calculate the median EpCAM values 
# for responders
long_CPSMA_res <- dplyr::filter(.data = long_CPSMA, resp.Group == "Responder")
head(long_CPSMA_res)

tCycle_RPSMA_medians <- c()
for(i in min(long_CPSMA_res$tCycle):max(long_CPSMA_res$tCycle)) {
    tCycle_RPSMA_medians[i + 1] <- median(long_CPSMA_res$MFI_values[long_CPSMA_res$tCycle == i])
}
# Note: in the i + 1, the +1 is crucial, as the cycle numbers are starting with 0 so the vector indexing would also start with a 0, and that
# is not possible, so I loose an element if I don't ad the +1
print(tCycle_RPSMA_medians)

# For nonresponders
long_CPSMA_nonres <- dplyr::filter(.data = long_CPSMA, resp.Group == "Nonresponder")
head(long_CPSMA_nonres)

tCycle_NRPSMA_medians <- c()
for (i in min(long_CPSMA_nonres$tCycle):max(long_CPSMA_nonres$tCycle)) {
    tCycle_NRPSMA_medians[i + 1] <- median(long_CPSMA_nonres$MFI_values[long_CPSMA_nonres$tCycle == i])
}
print(tCycle_NRPSMA_medians)


# I will also calculate the MAD values for each cycle
# for responders
tCycle_RPSMA_mad <- c()
for(i in min(long_CPSMA_res$tCycle):max(long_CPSMA_res$tCycle)) {
    tCycle_RPSMA_mad[i + 1] <- mad(long_CPSMA_res$MFI_values[long_CPSMA_res$tCycle == i])
}
tCycle_RPSMA_mad

# For nonresponders
tCycle_NRPSMA_mad <- c()
for (i in min(long_CPSMA_nonres$tCycle):max(long_CPSMA_nonres$tCycle)) {
    tCycle_NRPSMA_mad[i + 1] <- mad(long_CPSMA_nonres$MFI_values[long_CPSMA_nonres$tCycle == i])
}
print(tCycle_NRPSMA_mad)


# This loop will calculate the number of cells for each treatment cycle
# I decided not look at treatment cycles which did not have at least 10 cells
# as mean or median for these are unreliable
# for responders
respCellperCyc_PSMA <- c()
for (i in min(long_CPSMA_res$tCycle):max(long_CPSMA_res$tCycle)) {
    respCellperCyc_PSMA[i + 1] <- length(long_CPSMA_res$tCycle[long_CPSMA_res$tCycle == i])
}
print(respCellperCyc_PSMA)

# For nonresponders
nonrespCellperCyc_PSMA <- c()
for (i in min(long_CPSMA_nonres$tCycle):max(long_CPSMA_nonres$tCycle)) {
    nonrespCellperCyc_PSMA[i + 1] <- length(long_CPSMA_nonres$tCycle[long_CPSMA_nonres$tCycle == i])
}
print(nonrespCellperCyc_PSMA)


# I will create a new df with the PSMA medians and cycle numbers
# for responders
tCycle_RPSMA_medians_df <- data.frame(matrix(nrow = length(tCycle_RPSMA_medians), ncol = 5))
colnames(tCycle_RPSMA_medians_df) <- c("tCycle", "cell_number", "PSMA_medians", "PSMA_MAD", "resp.Group")
tCycle_RPSMA_medians_df [, 1] <- c(min(long_CPSMA_res$tCycle):max(long_CPSMA_res$tCycle))
tCycle_RPSMA_medians_df [, 2] <- respCellperCyc_PSMA
tCycle_RPSMA_medians_df [, 3] <- tCycle_RPSMA_medians
tCycle_RPSMA_medians_df [, 4] <- tCycle_RPSMA_mad
tCycle_RPSMA_medians_df [, 5] <- "Responder"
print(tCycle_RPSMA_medians_df)

# For nonresponders
tCycle_NRPSMA_medians_df <- data.frame(matrix(nrow = length(tCycle_NRPSMA_medians), ncol = 5))
colnames(tCycle_NRPSMA_medians_df) <- c("tCycle", "cell_number", "PSMA_medians", "PSMA_MAD", "resp.Group")
tCycle_NRPSMA_medians_df [, 1] <- c(min(long_CPSMA_nonres$tCycle):max(long_CPSMA_nonres$tCycle))
tCycle_NRPSMA_medians_df [, 2] <- nonrespCellperCyc_PSMA
tCycle_NRPSMA_medians_df [, 3] <- tCycle_NRPSMA_medians
tCycle_NRPSMA_medians_df [, 4] <- tCycle_NRPSMA_mad
tCycle_NRPSMA_medians_df [, 5] <- "Nonresponder"
print(tCycle_NRPSMA_medians_df)


# Unifying the PSMA means and MADs into one df and saving it
tCycle_PSMA_medians_unif <- rbind(tCycle_RPSMA_medians_df, tCycle_NRPSMA_medians_df)
write_csv(x = tCycle_PSMA_medians_unif, file = paste0(script_version, "/", "Treatment_cycles_PSMA_MFI_medians_MADs", "_", script_version, ".csv"))


# Plot the cycles data  for PSMA responders (In the responder group the 4th cycle only have 3 cells so it will be removed)
long_CPSMA_res_filt <- dplyr::filter(.data = long_CPSMA_res, tCycle != 4)

# Pairwise comparison
tCycle_RPSMA_wilcox_p <- ggplot(data = long_CPSMA_res_filt,
                                 aes(x = factor(tCycle), y = log10(MFI_values), fill = factor(tCycle), group = factor(tCycle))) +
    geom_pwc(method = "wilcox_test", label = "p = {p}", ref.group = "0", label.size = 6) + #this geom is for pairwise comparison and is amazing!
    geom_boxplot(show.legend = FALSE, position = "dodge") +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5,
                 position = "dodge") +
    scale_fill_manual(labels = c("EpCAM+PSMA+", "EpCAM+PSMA-", "EpCAM-PSMA+"), values = plotting_colors_tron[c(7, 3, 4, 6)]) +
    scale_color_manual(labels = c("EpCAM+PSMA+", "EpCAM+PSMA-", "EpCAM-PSMA+"), values = plotting_colors_tron[c(7, 3, 4, 6)]) +
    labs(x = "Treatment cycles", y = expression("PSMA log"[10] * " MFI values"), fill = "Treatment cycles") +
    #ggtitle("PSMA MFI comparison between treatment cycles - Responders") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(tCycle_RPSMA_wilcox_p)

ggsave(filename = paste0("PSMA_MFI_over_treatment_cycles_Responders_wilcox", "_", script_version, ".png"), tCycle_RPSMA_wilcox_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(tCycle_RPSMA_wilcox_p)


# Plot the cycles data for PSMA nonresponders

# Pairwise comparison
tCycle_NRPSMA_wilcox_p <- ggplot(data = long_CPSMA_nonres,
                                  aes(x = factor(tCycle), y = log10(MFI_values), fill = factor(tCycle), group = factor(tCycle))) +
    geom_pwc(method = "wilcox_test", label = "p = {p}", ref.group = "0", label.size = 6) + #this geom is for pairwise comparison and is amazing!
    geom_boxplot(show.legend = FALSE, position = "dodge") +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5,
                 position = "dodge") +
    scale_fill_manual(labels = c("EpCAM+PSMA+", "EpCAM+PSMA-", "EpCAM-PSMA+"), values = plotting_colors_tron[c(7, 3, 4, 6)]) +
    scale_color_manual(labels = c("EpCAM+PSMA+", "EpCAM+PSMA-", "EpCAM-PSMA+"), values = plotting_colors_tron[c(7, 3, 4, 6)]) +
    labs(x = "Treatment cycles", y = expression("PSMA log"[10] * " MFI values"), fill = "Treatment cycles") +
    #ggtitle("PSMA MFI comparison between treatment cycles - Nonresponders") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(tCycle_NRPSMA_wilcox_p)

ggsave(filename = paste0("PSMA_MFI_over_treatment_cycles_Nonresponders_wilcox", "_", script_version, ".png"), tCycle_NRPSMA_wilcox_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(tCycle_NRPSMA_wilcox_p)




## Response group based comparisons inside a single channel


# For EpCAM
long_CEpCAM_filt <- dplyr::filter(.data = long_CEpCAM, tCycle != 4)

tCycle_RvNREpCAM_wilcox_p <- ggplot(data = long_CEpCAM_filt,
                                    aes(x = factor(tCycle), y = log10(MFI_values), fill = resp.Group)) +
    geom_boxplot(show.legend = TRUE, position = position_dodge(width = 1, preserve = "single")) +
    geom_dotplot(show.legend = TRUE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5,
                 position = position_dodge(width = 1, preserve = "total")) +
    geom_pwc(method = "wilcox_test", label = "p = {p}", label.size = 6) + #this geom is for pairwise comparison and is amazing!
    scale_fill_tron() +
    #scale_fill_manual(values = brewer.pal(5, "RdYlBu")) +
    #stat_compare_means(label.y = 1.2, label.x = 1.2, method = "wilcox.test", paired = TRUE) +
    labs(x = "Treatment cycles", y = expression("EpCAM log"[10] * " MFI values"), fill = "Response groups") +
    #ggtitle("EpCAM MFI comparison between treatment cycles and response groups") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(tCycle_RvNREpCAM_wilcox_p)

ggsave(filename = paste0("EpCAM_MFI_over_treatment_cycles_Resp_vs_Nonresp_wilcox", "_", script_version, ".png"), tCycle_RvNREpCAM_wilcox_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(tCycle_RvNREpCAM_wilcox_p)

# For PSMA
long_CPSMA_filt <- dplyr::filter(.data = long_CPSMA, tCycle != 4)

tCycle_RvNRPSMA_wilcox_p <- ggplot(data = long_CPSMA_filt,
                                   aes(x = factor(tCycle), y = log10(MFI_values), fill = resp.Group)) +
    geom_boxplot(show.legend = TRUE, position = position_dodge(width = 1, preserve = "single")) +
    geom_dotplot(show.legend = TRUE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5,
                 position = position_dodge(width = 1, preserve = "total")) +
    geom_pwc(method = "wilcox_test", label = "p = {p}", label.size = 6) + #this geom is for pairwise comparison and is amazing!
    scale_fill_tron() +
    #scale_fill_manual(values = brewer.pal(5, "RdYlBu")) +
    #stat_compare_means(label.y = 1.2, label.x = 1.2, method = "wilcox.test", paired = TRUE) +
    labs(x = "Treatment cycles", y = expression("PSMA log"[10] * "MFI values"), fill = "Response groups") +
    #ggtitle("PSMA MFI comparison between treatment cycles and response groups") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(tCycle_RvNRPSMA_wilcox_p)

ggsave(filename = paste0("PSMA_MFI_over_treatment_cycles_Resp_vs_Nonresp_wilcox", "_", script_version, ".png"), tCycle_RvNRPSMA_wilcox_p,
       device = "png", path = paste0(script_version, "/", "Plots"),
       width = 3000, height = 1800, units = "px", dpi = 320)
rm(tCycle_RvNRPSMA_wilcox_p)

################################################################   End Section  ##################################################################


##############################################################################################################################################################################
################################################################   This is the end of the analysis for now  ##################################################################
##############################################################################################################################################################################













































