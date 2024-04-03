options(timeout = max(600, getOption("timeout")))
options(download.file.method.GEOquery = "curl")
options(scipen = 9)

library(GEOquery)
library(data.table)

#####################################################################################################
#                                                                                                   #
# Global variables
#                                                                                                   #
#####################################################################################################

#====================================================================================================
# Define allowable platforms
#====================================================================================================
platforms = c("illumina", "affymetrix", "agilent")

#====================================================================================================
# Define working directory
#====================================================================================================
setwd("C:/Users/weraz/Pictures/Sem6/R-bioconductor/tests")

#====================================================================================================
# Define datasets to be downloaded
#====================================================================================================
datasets = fread("./utils/Datasets.csv")
accession_numbers <- datasets[grep("^\\+", datasets$Info), `Accession number`]

#####################################################################################################
#                                                                                                   #
# Global functions
#                                                                                                   #
#####################################################################################################

getPlatform <- function(descr, accesion) {
  #============================================================================================================================
  # Search for platform name
  #     @descr: string
  #============================================================================================================================
  descr <- tolower(descr)
  result <- sapply(platforms, function(platform) grepl(platform, descr))
  
  # Check if any platform was found
  if (any(colSums(result) > 0)) {
    return(platforms[colSums(result) > 0])
  } else {
    stop(paste0("[", accession, "] Error: Platform not supported."))
  }
}

prepareRawFolder <- function(accession, baseDir) {
  #============================================================================================================================
  # Prepare /raw folder
  #     
  #============================================================================================================================
  # Create /raw
  if (!file.exists(paste0(baseDir, accession))) {
    suppressWarnings(dir.create(file.path(baseDir, accession, "raw"), recursive = TRUE))
  }
  # Check if TAR file exists
  tar_file <- paste0(file.path(baseDir, accession), "/", accession, "_RAW.tar")
  if (file.exists(tar_file)) {
    # Extract TAR file
    untar(tar_file, exdir = paste0(file.path(baseDir, accession), "/raw/"))
  } else {
    # Probably Illumina
    non_normalized_files <- c(list.files(path = file.path(baseDir, accession), pattern = "*raw*", full.names = TRUE), list.files(path = file.path(baseDir, accession), pattern = "*non-normalized*", full.names = TRUE))
    if (length(non_normalized_files) == 0){
      unlink(file.path(baseDir, accession), recursive = TRUE)
      stop(paste0("[", accession, "] Error: Supplementary raw data files not provided."))
    }
    
    file.rename(non_normalized_files[1], file.path(baseDir, accession, "raw", basename(non_normalized_files[1])))
  }
}

moveFolders <- function(originalDir, newDir) {
  #============================================================================================================================
  # Clean folders from data/accesion to data/platform/accesion
  #     
  #============================================================================================================================
  if (!file.exists(newDir)) {
    suppressWarnings(dir.create(newDir, recursive = TRUE))
  }
  
  files <- list.files(originalDir, full.names = TRUE)
  for (file in files) {
    file.copy(file, newDir)
  }
  unlink(originalDir, recursive = TRUE)
}

downloadGEO <- function(accession, baseDir = "./data"){
  #============================================================================================================================
  # Download both raw files and preprocessed data from GEO
  #
  #============================================================================================================================
  message(paste0("[", accession, "] Starting installation..."))
  
  # If path does not exist, create it
  if (!file.exists(paste0(baseDir, accession))) {
    suppressWarnings(dir.create(file.path(baseDir, accession), recursive = TRUE))
  }
  
  # Processed data
  geo <- getGEO(accession, destdir=file.path(baseDir, accession), GSEMatrix = TRUE)
  geo <- geo[[1]]
  platform <- getPlatform(paste0(geo$data_processing, geo$scan_protocol, geo$hyb_protocol, geo$description), accession)

  # Move to directory platform
  moveFolders(file.path(baseDir, accession), file.path(baseDir, platform, accession))


  # Raw files
  getGEOSuppFiles(accession, baseDir=file.path(baseDir, platform))
  prepareRawFolder(accession, file.path(baseDir, platform))

  message(paste0("[", accession, "] Dataset installed successfully"))
  return(TRUE)
}

#####################################################################################################
#                                                                                                   #
# Main
#                                                                                                   #
#####################################################################################################

downloaded <- 0
for (accession_number in accession_numbers){
  skip_to_next <- FALSE
  tryCatch(downloadGEO(accession_number), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { 
    next 
  }
  downloaded <- downloaded + 1
}
message(paste0(downloaded, "/", length(accession_numbers), " installed successfully"))