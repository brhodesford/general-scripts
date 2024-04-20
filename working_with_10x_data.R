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
<<<<<<< HEAD
=======
library(celldex)
library(SingleR)
library(SingleCellExperiment)
library(pheatmap)
library(ggplot2)
library(harmony)
library(CellChat)
library(dittoSeq)
library(clusterProfiler)

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
# This function is called by FileParsing and makes Seurat Objects from matrix,
# features, and barcode files 
#
# Args : 
#   directory: Parent directory containing directories containing each of the 
#   files mentioned above
# 
# Returns: List of seurat objects that can be used in RunAndVizQC


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
# This function goes through given directory and 1) sorts the files into 
# subdirectories that contain one of each of the necessary files to be used 
# to make 10x object and Seurat object, 2) calls FileParsing10xDirs to make 
# a Seurat objext 
#
# Args : 
#   directory: Parent directory containing either 1) many labeled matrix, 
#   barcode, and feature files that can be moved into distinct subdirectories 
#   or 2) subdirectories each containing each of the files mentioned above
# 
# Returns: List of seurat objects that can be used in RunAndVizQC


RunAndVizQC <- function(list.of.seurat.objs, outputdir, info.file.path = NULL){
  
  # Check if seurat objects in list have percent.mt and add if not
  for (i in 1:length(list.of.seurat.objs)){
    orig.ident.name <- as.character(unique(list.of.seurat.objs[[i]]$orig.ident))
    if (!("percent.mt" %in% colnames(list.of.seurat.objs[[i]][[]]))){
      cat("percent.mt is being added to", orig.ident.name, "\n")
      list.of.seurat.objs[[i]][["percent.mt"]] <- PercentageFeatureSet(list.of.seurat.objs[[i]],
                                                                       pattern = "^MT-")
    }
  }
  # Check if info.file.path provided and run preQC version if not
  if (is.null(info.file.path)){
    print('Running script with no info sheet. Outputs include violin plots and table giving sample names and preQC cell counts, which can be edited to include nFeatures, nCounts, and percent.mt filter values and provided to RunAndVizQC as info.file')

    # Write output path for the table with QC information
    table.path.pre <- paste(outputdir,"/FiltersAndFeatures_preQC.csv", sep = "")
    # Make empty table for important information
    table.var <- list()

    # for loop to make violin plots and make first column of table for export

    for (i in 1:length(list.of.seurat.objs)){
      table.var[[i]] <- unique(Idents(list.of.seurat.objs[[i]]))
      if (file.exists(outputdir)) {

        filename <- paste(outputdir, "/", list.of.seurat.objs[[i]]$orig.ident,
                          "_qc_violinplot.pdf", sep = "")

        # Visualize with Violin Plots
        pdf(filename)
        print(VlnPlot(list.of.seurat.objs[[i]],
                      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                      ncol = 3))
        dev.off()
        print(paste("PDF can be found at ", unique(filename), sep = ""))

      }else{
        stop("That output directory does not exist")
      }
    }

    # Convert list from above to data.frame and name columns appropriately
    table.var <- as.data.frame(unlist(table.var))
    colnames(table.var)[1] <- "Sample"

    # Add Pre QC number of Cells Value to table
    table.var$"Pre_QC_Cells" <- unlist(lapply(list.of.seurat.objs, function(seurat.obj) {
      length(colnames(seurat.obj))
    }))

    # Add Pre QC average of nFeature_RNA to table
    table.var$"PreQC_Average_nFeature_RNA" <- unlist(lapply(list.of.seurat.objs, function(seurat.obj) {
      mean(seurat.obj$nFeature_RNA)
    }))

    # Add Pre QC average of nCount_RNA to table
    table.var$"PreQC_Average_nCount_RNA" <- unlist(lapply(list.of.seurat.objs, function(seurat.obj) {
      mean(seurat.obj$nCount_RNA)
    }))

    # Add Pre QC average of percent.mt to table
    table.var$"PreQC_Average_percent.mt" <- unlist(lapply(list.of.seurat.objs, function(seurat.obj) {
      mean(seurat.obj$percent.mt)
    }))

    # Add Pre QC number of Cells Value to table
    table.var$"QC_Filter_nFeature_RNA" <- NA
    table.var$"QC_Filter_nCount_RNA" <- NA
    table.var$"QC_Filter_percent.mt" <- NA

    # Print information about csv and write out file
    print("Structure of CSV can be found below")
    print(str(table.var))
    write.csv(table.var,table.path.pre,row.names = F)
    print(paste("CSV can be found at ", table.path.pre, sep = ""))
<<<<<<< HEAD
=======
    return(list.of.seurat.objs)
>>>>>>> 2650a4b (added functions, including wrappers to chunk Seurat workflow)

  }else if(file.exists(info.file.path)){
    # ^ Check if info.file.path provided and exists and run postQC version if so
    
    print("Running script with info sheet. Outputs include table adding postQC cell count and cells lost columns. Returns list of filtered Seurat objects")

    table.path.post <- paste(outputdir,"/FiltersAndFeatures_postQC.csv", sep = "")

    # Read in table found at <info.file.path>
    table.var <- read.csv(info.file.path)

    if(nrow(table.var) == length(list.of.seurat.objs)){
      # Get and print out pertinent filtering information
      nfeature <- table.var$"QC_Filter_nFeature_RNA"
      ncount <- table.var$"QC_Filter_nCount_RNA"
      percent.mt <- table.var$"QC_Filter_percent.mt"
      cat("nFeatures filter is: ", nfeature, "\n")
      cat("nCounts filter is: ", ncount, "\n")
      cat("percent.mt filter is: ", percent.mt, "\n")

      # Apply each filter to the appro
      list.of.seurat.objs.filtered <- list()
      for (i in 1:nrow(table.var)){
        temp.seurat.obj <- list.of.seurat.objs[[i]]
        
        list.of.seurat.objs.filtered[[i]] <- subset(x = temp.seurat.obj,(nFeature_RNA > table.var[i, "QC_Filter_nFeature_RNA"]) & (nCount_RNA < table.var[i, "QC_Filter_nCount_RNA"]) & (percent.mt < table.var[i, "QC_Filter_percent.mt"]))
      }
      
      # Add Post_QC_Cells Column
      table.var["Post_QC_Cells"] <- unlist(lapply(list.of.seurat.objs.filtered,
                                                    function(seurat.obj) {
                                                      length(colnames(seurat.obj))
                                                      }))
      
      # Add Cells_Lost Column
      table.var["Cells_Lost"] <- table.var["Pre_QC_Cells"] - table.var["Post_QC_Cells"]
      
      # Output Path for CSV
      table.path.post <- paste(outputdir,"/FiltersAndFeatures_postQC.csv", sep = "")
      
      # Print information about csv and write out file
      print("Structure of CSV can be found below")
      print(str(table.var))
      write.csv(table.var,table.path.post,row.names = F)
      print(paste("CSV can be found at ", table.path.post, sep = ""))
      
<<<<<<< HEAD
      
=======
      return(list.of.seurat.objs.filtered)
>>>>>>> 2650a4b (added functions, including wrappers to chunk Seurat workflow)
    }else{
      stop("The number of rows in the info file provided does not match the number of samples provided in the list of seurat objects. Please check these objects and try again. If you continute to have trouble, a starting template can be provided if this script is run without the info.file.path argument.")
    }
  }else{
    stop("The info file provided is not valid. Check file path and try again.")
  }
<<<<<<< HEAD
  return(list.of.seurat.objs.filtered)
=======
>>>>>>> 2650a4b (added functions, including wrappers to chunk Seurat workflow)
}
# This function has two main functions 1) perform QC filtering as indicated in 
# the info file and 2) make a table with QC information pre and post QC 
# filtering. It can be run sequentially without and then with an info file. The 
# first run will make PDFs showing violin plots for the three filtering 
# parameters and create a template file that can be used to fill in and then 
# provided as the info file with filtering information. 
#
# Args : 
#   list.of.seurat.objs: Can provide list made in a few different ways: 
#     1) manually create a list of already made Seurat objects, 2) provide 
#     returned value of FileParsing function above, 3) provide returned value of 
#     first of two sequential runs of this function (will have percent.mt added to 
#     meta.data)
#   outputdir: Specify output directory for PDFs of Violin Plots and CSV with QC 
#     filtering information
#   info.file.path: Default is NULL (used to generate template, which then can 
#     be filled out for use in second sequential run). Otherwise, provide path 
#     to file with QC filtering information. The column titles for 
#     filtering information must be as follows:
#       QC_Filter_nFeature_RNA
#       QC_Filter_nCount_RNA
#       QC_Filter_percent.mt
#    For nFeature_RNA values > (greater than provided) will be filtered
#    For nCount_RNA values < (less than provided) will be filtered
#    For precent.mt values < (less than provided) will be filtered
# 
# Files Created: PDFs for Violin Plots for each Seurat Object, FiltersAndFeatures 
# CSV Template (pre_QC) or FiltersAndFeatures postQC CSV
# Returns: List of seurat objects that can be used for downstream Seurat analysis

## ---------------------------


ListMerge <- function(list.of.seurat.objs, project.id){
  
  list.cell.ids <- as.character(unlist(lapply(list.of.seurat.objs,
                                                   function(x){
                                                     unique(x$orig.ident)
                                                   })))
  
  # merge all seurat objects in list into one large seurat object
  large.seurat <- Merge_Seurat_List(
    list.of.seurat.objs,
    add.cell.ids = list.cell.ids,
    merge.data = TRUE,
    project = project.id
  )
  
  saveRDS(large.seurat, paste(project.id, "_mergedSeurat.rds", sep = ""))
  
  return(large.seurat)
}
# This function takes a list of Seurat objects, adds a new parameter called 
# orig.dataset.name, and merges each object into one large Seurat object that 
# can be used for downstream analysis. 
#
# Args: 
#   list.of.seurat.objs: Self-explanatory list of Seurat objects to be merged
#
#   project.id: name for the project to name files
#
# Files Created: RDS for large.seurat object
# Returns: merged large.seurat object


NormAndPCA <- function(seurat.obj, obj.name, n.var.features = 2000){
  # Normalize using LogNormalize and scale data by 10000 (Seurat defaults)
  seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find Variable Features
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = n.var.features)
  
  # Scale the data
  all.genes <- rownames(seurat.obj)
  seurat.obj <- ScaleData(seurat.obj, features = all.genes)
  
  # Run principle component analysis
  seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
  print("Saving Pre-Clustering RDS")
  saveRDS(seurat.obj,paste(obj.name, "_preClustering.rds", sep = ""))
  
  print("Saving Dimmension Heatmap")
  # Visualize those components
  pdf(paste(obj.name, '_dimHeatmap.pdf', sep = ""))
  print(DimHeatmap(seurat.obj, dims = 1:15, cells = 500, balanced = TRUE))
  dev.off()
  
  print("Saving Elbow Plot")
  pdf(paste(obj.name, "_ElbowPlot.pdf", sep = ""))
  print(ElbowPlot(seurat.obj, ndims = 50))
  dev.off()
  
  return(seurat.obj)
}
# This function takes a seurat object and performs the initial steps of Seurat 
# analysis, including Normalization, Variable Feature identification, Scaling,
# and Principle Component Analysis. Dimension heatmaps and an elbow plot are 
# outputs in order to identify dimensionality for clustering.
#
# Args:
#   seurat.obj: Self-explanatory Seurat object on which to perform the normalization, 
#   scaling, and principle component analyses
#   
#   obj.name: name of the Seurat object to be used to name output files.
#   
#   n.var.features (optional): number of features to include for the find 
#   variable features. Default is 2000.
#
# Files Created: *_dimHeatmap.pdf is heatmaps for 15 first dimensions of variable
#                 features from PCA
#                 
#                *_ElbowPlot is an elbow plot to identify the inflection point 
#                 value, which can be used for clustering
#
# Returns: seurat.obj with these new values included

PerformSeuratClustering <- function(seurat.obj, obj.name, n.dims, resolution.vector, gene.lis){
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:n.dims)
  
  if (length(resolution.vector == 1)){
    cluster.res <- resolution.vector[1]
    seurat.obj <- FindClusters(seurat.obj, resolution = cluster.res)
    # Look at cluster IDs of the first 5 cells
    print(paste("Final resolution provided as ", cluster.res))
    print(head(Idents(seurat.obj), 5))
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:n.dims)
    
    # Plot UMAP in PDF
    print(paste("Saving UMAP with ", n.dims, " dimensions and ", cluster.res, " resolution"))
    pdf(paste(obj.name, "_UMAP_", n.dims, "dimensions", cluster.res, "resolution.pdf", sep = ""))
    print(DimPlot(seurat.obj, reduction = "umap"))
    dev.off()
  }else{
    for (i in 1:length(resolution.vector)){
      cluster.res <- resolution.vector[i]
      print(cluster.res)
      # make temporary seurat object because no values will be saved, just PDFs of each UMAP made
      temp.seurat.obj <- FindClusters(seurat.obj, resolution = cluster.res)
      # Look at cluster IDs of the first 5 cells
      print(head(Idents(temp.seurat.obj), 5))
      temp.seurat.obj <- RunUMAP(temp.seurat.obj, dims = 1:n.dims)
      
      # Plot UMAP in PDF
      print(paste("Saving UMAP with ", n.dims, " dimensions and ", cluster.res, " resolution", sep = ""))
      pdf(paste(obj.name, "_UMAP_", n.dims, "dimensions_", cluster.res, "resolution.pdf", sep = ""))
      print(DimPlot(seurat.obj, reduction = "umap"))
      dev.off()
    }
    # if more than one resolution is provided, then no changes are made to the
    # seurat.obj itself. Will need to be run with chosen final resolution 
    # (length(resolution.vector) == 1)
    print("No final resolution was provided, so no changes to seurat.obj were 
        saved. Rerun with final resolution provided (length == 1) as
        resolution.vector")
  }
  print("Saving Post Clustering RDS")
  saveRDS(seurat.obj, paste(obj.name, "_postClustering.rds", sep = ""))
  return(seurat.obj)
}

# This function takes a seurat object and performs the FindNeighbors, Clustering, 
# and UMAP steps of Seurat pipeline. It takes a vector to try a variety of 
# resolution values for clustering. Once the resolution value is decided, it can 
# be provided, and the seurat.obj with these values added will be returned.
#
# Args:
#   seurat.obj: Self-explanatory Seurat object on which to perform the 
#   clustering and UMAP
#   
#   obj.name: name of the Seurat object to be used to name output files.
#   
#   n.dims: number of dimensions to use for FindNeighbors and UMAP functions
#   
#   resolution.vector: vector of values to be tried or single value if final 
#   resolution decided
#
# Files Created: 
#                
#                 
#
# Returns: seurat.obj with these new values included

GetSeuratTop10 <- function(seurat.obj, obj.name){
  # find markers for every cluster compared to all remaining cells, report only the positive
  # ones
  seurat.markers <- FindAllMarkers(seurat.obj, only.pos = TRUE)
  seurat.markers.filtered <- seurat.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)
  write.csv(seurat.markers.filtered, paste(obj.name, "_markersDGE_filtered.csv", sep = ""))
  
  seurat.markers.filtered %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  
  pdf(paste(obj.name, "_top10Markers.pdf", sep = ""), width = 20, height = 10)
  print(DoHeatmap(seurat.obj, features = top10$gene) + theme(axis.text=element_text(size=3)))
  dev.off()
  
  return(seurat.markers.filtered)
}

# This function takes a Seurat object and gets the top 10 markers of each 
# cluster in order to identify which cell types these clusters represent. It is
# used in the Annotation SingleR function but can also be used independently
#
#
#
# Args:
#   seurat.obj: Self-explanatory Seurat object on which to perform 
#   
#   
#   obj.name: name of the Seurat object to be used to name output files.
#   
# Files Created: *_markersDGE_filtered.csv: differential gene expression table
#                
#               *_top10Markers.pdf: Heatmap showing top 10 features of each 
#                                   cluster
#
# Returns: seurat.obj with these new values included


DifferentialAndAnnotation <- function(seurat.obj, obj.name, gene.list, anno.type = "manual"){
  if(anno.type == "manual"){
    print("Performing Differential Gene Expression Analysis and Visualization for Manual Annotation")
  }else if (anno.type == "automated"){
    print("Performing Differential Gene Expression Analysis, Visualization, and Automated Annotation with SingleR")
  }else
    stop("Only manual or automated are acceptable inputs for anno.type. Please try again")
  # Hard coded gene lists for cell types in human tumors
  
  # Get differentially expressed genes
  seurat.markers.filtered <- GetSeuratTop10(seurat.obj, obj.name)
  
  # Make heatmap of hard coded cell type markers above
  pdf(paste(obj.name, "_cannonicalCellMarkers_Heatmap.pdf", sep = ""), width = 20, height = 10)
  print(DoHeatmap(seurat.obj, features = gene.list) + NoLegend())
  dev.off()
  
  # Only perform SingleR if 
  if(anno.type == "automated"){
    ref <- celldex::HumanPrimaryCellAtlasData()
    # convert Seurat object to SingleCellExperiment object
    sce.obj <- as.SingleCellExperiment(seurat.obj)
    
    # Save SCE object as RDS
    saveRDS(sce.obj,paste(obj.name,"_postclustering_sce.rds", sep=""))
    
    # Run SingleR Annotation on SCE object using previously described Seurat clusters
    obj.anno <- SingleR(sce.obj,ref = ref,  labels =ref$label.main,clusters = sce.obj$seurat_clusters)
    
    table(obj.anno$labels)
    
    tab <- table(Assigned=obj.anno$pruned.labels, Cluster=levels(sce.obj$seurat_clusters))
    
    pdf(paste(obj.name,"_SingleR_Anno.pdf",sep = ""))
    print(plotScoreHeatmap(obj.anno))
    # Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
    print(pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101)))
    dev.off()
    dev
    
    tab.df <- as.data.frame(tab)
    tab.df <- tab.df[tab.df$Freq==1,]
    
    ident.labels <- tab.df$Assigned[match(seurat.obj$seurat_clusters, tab.df$Cluster)]
    
    # Add the new metadata column 'ident.labels' to the Seurat object
    seurat.obj <- AddMetaData(seurat.obj, metadata = ident.labels, col.name = "ident.labels")
  }
  return(seurat.obj)
} 
# This function takes a Seurat object, performs an automated annotation and
# provides outputs that may be helpful in performing a manual annotation, such
# as heatmaps and feature plots with canonical markers, as well as a heatmap 
# showing top 10 markers for each cluster.
#
#
# Args:
#   seurat.obj: Self-explanatory Seurat object on which to perform 
#   
#   
#   obj.name: name of the Seurat object to be used to name output files.
#   
# Files Created: 
#   *_cannonicalCellMarkers_Heatmap.pdf: Heatmap of hardcoded canonical cell 
#                                        markers, to inform manual clustering or
#                                        confirm automated clustering
#
#   *_cannonicalCellMarkers_featurePlot.pdf: Feature Plot of hardcoded canonical  
#                                            cellmarkers, to inform manual 
#                                            clustering or confirm automated 
#                                            clustering
#
#   *_postclustering_sce.rds: RDS of the object converted ot the 
#                             SingleCellExperiment object type
#
#   *_SingleR_Anno.pdf: PDF of two outputs of annotation program, including 
#                 
#
# Returns: seurat.obj with these new values (ident.labels metadata object from 
#          annotation) included

RunCellChat <- function(seurat.obj, obj.name){
  ##### Cell Chat ##### 
  ## Create Cell Chat Object
  # transform seurat obj to cellchat obj
  cellchat <- createCellChat(object = seurat.obj, # this is my seurat obj
                             group.by = "ident.labels")
  # set the database#
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)
  # use a subset of CellChatDB for cell-cell communication analysis
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
  
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  
  ## Preprocessing the expression data for cell-cell communication analysis
  # subset the expression data of signaling genes for saving computation cost
  print(cellchat@idents)
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multisession", workers = 4) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  cellchat <- projectData(cellchat, PPI.human)
  print("83")
  
  # Compute the communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat, type = "triMean")

  cellchat <- filterCommunication(cellchat, min.cells = 10)

  # Extract the inferred cellular communication network as a data frame
  df.net <- subsetCommunication(cellchat)

  write.csv(df.net, paste(obj.name,"_cellChat.csv"))
  
  return(df.net)
}
# This function takes a seurat object and 
#
#
#
#
#
# Args:
#   seurat.obj: Self-explanatory Seurat object on which to perform 
#   
#   
#   obj.name: name of the Seurat object to be used to name output files.
#   
#
# Files Created: 
#                
#                 
#
# Returns: seurat.obj with these new values included



CellTypeSpecificProccessing <- function(seurat.obj, obj.name, cluster.numbers, cluster.name){ 
  # multiple cluster numbers that refer to a single cell type can be provided, 
  # only one cluster.name can be provided
  
  file.output.name <- paste(obj.name, cluster.name, sep = "_")
  
  celltype.obj <- subset(seurat.obj, subset = seurat_clusters %in% cluster.numbers)
  
  # Find Variable Features
  celltype.obj <- FindVariableFeatures(celltype.obj, selection.method = "vst", nfeatures = 2000)
  
  # Scale the data
  all.genes <- rownames(celltype.obj)
  celltype.obj <- ScaleData(celltype.obj, features = all.genes)
  
  # Run principal component analysis
  celltype.obj <- RunPCA(celltype.obj, features = VariableFeatures(object = seurat.obj))
  saveRDS(celltype.obj, paste(file.output.name, "_preClustering.rds", sep = ""))
  
  # Visualize principal components
  pdf(paste(file.output.name, '_dimHeatmap.pdf', sep = ""))
  print(DimHeatmap(celltype.obj, dims = 1:15, cells = 500, balanced = TRUE))
  dev.off()
  
  pdf(paste(file.output.name, "_ElbowPlot.pdf", sep = ""))
  print(ElbowPlot(celltype.obj))
  dev.off()
  
  return(celltype.obj)
}


GetCellTypeGeneList <- function(celltype.obj, obj.name){
  # Get differentially expressed genes and visualize top 10
  cell.markers <- FindAllMarkers(celltype.obj, only.pos = TRUE)
  cell.markers.filtered <- cell.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)
  write.csv(cell.markers.filtered, paste(obj.name, "_markersDGE_filtered.csv", sep = ""))
  
  cell.markers.filtered %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  write.csv(top10, paste(obj.name, "_markersDGE_filtered_top10.csv", sep = ""))
  
  pdf(paste(obj.name, "_top10Markers.pdf", sep = ""), width = 20, height = 10)
  DoHeatmap(GSE181919.Tcells, features = top10$gene) + theme(axis.text=element_text(size=3))
  dev.off()
  
}


## ---------------------------

