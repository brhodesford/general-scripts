## ---------------------------
##
## Copyright (c) Rhodes Ford, 2024
## Author: Dr. B. Rhodes Ford
## Email: rhodes.ford@ucsf.edu
## Date Created: 2024-02-20
## Version: v1.0
## 
## 
## 
## Script name: Populate Directory
##
## 
##
## ---------------------------
## Notes
## 
## Purpose of script: To add my prefered subdirectories to provided directory 
## 
## Input: parent directory
##
## Outputs: parent directory populated with parent subdirectories
##   
##    
## ---------------------------

subdir.list = c('data', 'scripts', 'figs', 'tables',  'prose')

## Define Functions
PopulateReadMe <- function(directory){
  
  # Specify the README file path (directory + README.md)
  readme.file <- paste(directory, 'README.md')
  
  current.year <- format(Sys.Date(), "%Y")
  current.date <- format(Sys.Date())
  
  
  # Open the README file for writing or append if it already exists
  file.conn <- file(readme.file, 
                    open = ifelse(file.exists(readme.file), "a", "w"))
  
  # Write content to the README file
  cat("## ---------------------------\n\n", 
      file = file.conn)
  
  cat(paste("## Copyright (c) Rhodes Ford, ", 
            current.year, "\n", 
            sep = ""), 
      file = file.conn)
  
  cat("## Author: Dr. B. Rhodes Ford\n",file = file.conn)
  
  cat("## Email: rhodes.ford@ucsf.edu\n", file = file.conn)
  
  cat(paste("## Date Created: ", 
            current.date, 
            "\n\n", 
            sep = ""), 
      file = file.conn)
  
  cat("## ---------------------------\n\n", file = file.conn)
  
  cat("## Project Title: n\n", file = file.conn)
  
  cat("##Project Description: \n\n", file = file.conn)
  
  cat("## ---------------------------\n\n", file = file.conn)
  
  cat(paste("## Parent Directory: ", 
            directory, 
            "\n\n", 
            sep = ""), 
      file = file.conn)
  
  cat("## Subdirectories: ", file = file.conn)
  
  cat(paste(subdir.list, 
            collapse = ", "), 
      file = file.conn)
  
  cat("\n\n", file = file.conn)
  
  cat("## ---------------------------\n\n", file = file.conn)

}
# This function is run in the overall RunPopulate command to make and populate 
# a README file to be filled out for a specific project.
#
# Args : 
#   directory: Parent directory for the project
# 
# Returns: No return value. The files are written into the directory.


PopulateSubDirs <- function(directory){
  
  # Make subdir.list with parent directory
  subdir.list <- paste(directory, subdir.list, sep = "/")
  
  # Check if each file in subdir.list exists
  sapply(subdir.list, function(dir) {
    if (!file.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  })
}
# This function is run in the overall RunPopulate command to make and populate 
# the necessary subdirectories (listed above in line 26)
#
# Args : 
#   directory: Parent directory for the project
# 
# Returns: No return value. The files are written into the directory.

RunPopulate <- function(directory){
  
  # Run PopulateReadMe function to create and populate README file within 
  # parent directory
  PopulateReadMe(directory)
  
  PopulateSubDirs(directory)

}
# This function executes PopulateReadMe and PopulateSubDirs functions
#
# Args : 
#   directory: Parent directory for the project
# 
# Returns: No return value. The files are written into the directory.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  
  stop("No trailing arguments provided.\n")

} else if (length(args) == 1){
  
  cat("One trailing argument provided. ")
  cat("Proceeding with analysis with default list of subdirectories\n")
  
} else {
  
  stop("Two trailing arguments provided. This script cannot use provided directory list.")
  
}

RunPopulate(args[1])

## ---------------------------
## Notes:
## 
## Improvements to be made:
##   Take second argument of provided subdirectory list.
##   Possibility of automated reorganization?