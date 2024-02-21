## ---------------------------
##
## Copyright (c) Rhodes Ford, 2024
## Author: Dr. B. Rhodes Ford
## Email: rhodes.ford@ucsf.edu
## Date Created: 2024-02-13
## 
## 
## 
## Script name: Seurat scRNAseq processing
##
## 
##
## ---------------------------
## Notes
## 
## Purpose of script: To take publicly available or internal scRNAseq data and
## write a script
## 
## README detailing any pre-R processing steps can be found within parent directory
## Inputs can either be downloaded directly from GEO or processed w/ cellranger
## Inputs should be placed in a directory called "rawdata" within project dir
##
## Input (option 1): Directory with multiple directories, 
##   each with *matrix.mtx.gz, *barcodes.tsv.gz, and *genes.tsv.gz files 
## Input (option 2): A single directory with multiple samples each having 
##   a corresponding *matrix.mtx.gz, *barcodes.tsv.gz, and *genes.tsv.gz file 
##
## Outputs: Found in the following directories
##   
##    
##
## ---------------------------



## load in the packages we will need

library(dplyr)
library(Seurat)
library(patchwork)

# load in any necessary functions I've written functions
#source("~/general_scripts/")       # loads up all the packages we need

## ---------------------------

## Define Functions

FileParsing <- function(directory) {
  # directory can include either other directories, with each of the three necessary 10x files or it can have multiple named versions of these files
  
  dir.contents <- list.files(directory)
  
  dir.contents <- lapply(dir.contents, 
                         function(x) paste(directory, x, sep = "/"))
  
  dir.bool <- unlist(lapply(dir.contents, 
                            function(x) file_test("-d", x)))
  
  file.bool <- unlist(lapply(dir.contents, 
                             function(x) file_test("-f", x)))
  
  print(dir.contents)
  
  if (length(dir.contents[dir.bool]) == length(dir.contents)) {

    print("There are only directories in this parent directory. Analysis begun.")

    seuratObjList <- list()

    for (i in 1:length(dir.contents)) {
      # Load the PBMC dataset
      obj.10x <- Read10X(data.dir = paste(dir.contents[i], "/", sep = ""))
      # Initialize the Seurat object with the raw (non-normalized data).
      obj.seurat <- CreateSeuratObject(counts = obj.10x, project = filelist[i], min.cells = 3, min.features = 200)
      seurat.obj.list[[i]] <- obj.seurat
      }
  } else if (length(dir.contents[file.bool]) == length(dir.contents)) {
      print("There are only files in this parent directory. Analysis begun.")

      samplelist <- dir.contents[grep(pattern = "matrix.mtx.gz", x = dir.contents)]
      filepattern <- gsub("_matrix.mtx.gz", "", samplelist)

      print("The following samples are present in the provided directory:")
      print(filepattern)

      for (i in 1:length(filepattern)) {

      }

    seuratObjList <- list()
    } else {
      stop("There are both files and directories in this parent. Please make sure file organization is consistent with either mtx, barcodes, and genes files in distinct directories or named and found in parent directory with no other directories present")
    }
    return(seuratObjList)
}

## ---------------------------

## Main (Executed statements)

# Set working directory
setwd("/c4/home/brford/breastCancerProject/")


# Run functions
gse206638 <- seuratProcessing("GSE206638_RAW/")
gse114724 <- seuratProcessing("GSE114724_RAW")

