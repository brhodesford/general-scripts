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
## Outputs: Seurat objects for each sample included in 10x dataset
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

FileParsing10xDirs <- function(directory){
  
  file.list <- list.files(directory)
  
  dir.contents <- lapply(file.list, 
                         function(x) paste(directory, x, sep = "/"))
  
  # Make list to hold Seurat objects
  seurat.obj.list <- list()
  
  for (i in 1:length(dir.contents)) {
    # Make the 10x object
    obj.10x <- Read10X(data.dir = paste(dir.contents[i], "/", sep = ""))
    
    # Initialize the Seurat object with the raw (non-normalized data).
    obj.seurat <- CreateSeuratObject(counts = obj.10x, project = file.list[i], min.cells = 3, min.features = 200)
    
    # put seurat object into a list
    seurat.obj.list[[i]] <- obj.seurat
  }
  return(seurat.obj.list)
}


FileParsing <- function(directory) {
  # directory can include either other directories, with each of the three necessary 10x files or it can have multiple named versions of these files
  
  file.list <- list.files(directory)
  
  dir.contents <- lapply(file.list, 
                         function(x) paste(directory, x, sep = "/"))
  
  dir.bool <- unlist(lapply(dir.contents, 
                            function(x) file_test("-d", x)))
  
  file.bool <- unlist(lapply(dir.contents, 
                             function(x) file_test("-f", x)))
  
  print(dir.contents)
  
  if (length(dir.contents[dir.bool]) == length(dir.contents)) {
    
    print("There are only directories in this parent directory. Directory analysis begun.")
    seurat.obj.list <- FileParsing10xDirs(directory)
    
  } else if (length(dir.contents[file.bool]) == length(dir.contents)) {
    print("There are only files in this parent directory. Analysis begun.")
    
    # Get list of matrix, barcode, and genes files
    sample.list.mat <- file.list[grep(pattern = "matrix.mtx.gz", x = file.list)]
    sample.list.mat <- gsub("_matrix.mtx.gz", "", sample.list.mat)
    
    sample.list.bar <- file.list[grep(pattern = "barcodes.tsv.gz", x = file.list)]
    sample.list.bar <- gsub("barcodes.tsv.gz", "", sample.list.bar)
    
    if(length(file.list[grep(pattern = "genes.tsv.gz", x = file.list)]) != 0){
      sample.list.genes <- file.list[grep(pattern = "genes.tsv.gz", x = file.list)]
      sample.list.genes <- gsub("genes.tsv.gz", "", sample.list.genes)
    }else if(length(file.list[grep(pattern = "features.tsv.gz", x = file.list)]) != 0){
      sample.list.genes <- file.list[grep(pattern = "features.tsv.gz", x = file.list)]
      sample.list.genes <- gsub("features.tsv.gz", "", sample.list.genes)
    }
    
    # Make sure lengths of the lists of each file type are the same length (all files provided)
    if(length(sample.list.mat) == length(sample.list.bar) & length(sample.list.mat) == length(sample.list.genes)){
      sample.list <- sample.list.mat
    }else{
      stop("Not all samples included in this directory have all three necessary files. Please make sure to include all three necessary files for each sample")
    }
    
    print("The following samples are present in the provided directory with all three 10x files:")
    print(sample.list)
    
    # Make new parent directory
    new.parent.dir <- paste(directory, '_matched_files', sep ="")
    dir.create(new.parent.dir)
    
    for (i in 1:length(sample.list)) {
      # Get the three 10x files for sample i
      sample.files <- file.list[grepl(pattern = sample.list[i],x = file.list)]
      
      # Make new directory containing the three 10x files for that sample
      new.dir <- paste(new.parent.dir,sample.list[i],sep = "/")
      dir.create(new.dir)
      
      # Move three 10x files into the new directory 
      file.rename(paste(directory,sample.files, sep = "/"), paste(new.dir, sample.files, sep='/'))
      
      newdir.file.list <- list.files(new.dir)
      if(length(newdir.file.list) != 3){
        stop("Something has gone wrong. The three necessary 10x files were not transferred properly")
      }else{
        newdir.file.path <- paste(new.dir,newdir.file.list,sep = '/')
        for(j in 1:length(newdir.file.path)){
          if(grepl('barcodes.tsv.gz',newdir.file.path[j])){
            file.rename(newdir.file.path[j],paste(new.dir,'barcodes.tsv.gz', sep = "/"))
          }else if(grepl('matrix.mtx.gz',newdir.file.path[j])){
            file.rename(newdir.file.path[j],paste(new.dir,'matrix.mtx.gz', sep = "/"))
          }else if(grepl('genes.tsv.gz',newdir.file.path[j])){
            file.rename(newdir.file.path[j],paste(new.dir,'features.tsv.gz', sep = "/"))
          }else if(grepl('features.tsv.gz',newdir.file.path[j])){
            file.rename(newdir.file.path[j],paste(new.dir,'features.tsv.gz', sep = "/"))
          }
        }
      }
    }
    # Run directory analysis with new directory-separated data
    seurat.obj.list <- FileParsing10xDirs(new.parent.dir)
    
  } else {
    stop("There are both files and directories in this parent. Please make sure file organization is consistent with either mtx, barcodes, and genes files in distinct directories or named and found in parent directory with no other directories present")
  }
  return(seurat.obj.list)
}

## ---------------------------