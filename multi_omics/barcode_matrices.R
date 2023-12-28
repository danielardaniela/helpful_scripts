# Script to Analyze Intersection of Barcodes Across Arrays, HTO Libraries, and ADT Libraries
# Author: Daniela D. Russo
# Date: 12/27/2023
# Usage: Rscript barcode_matrices.R "array1,array2,..." "HTO1,HTO2,..." "ADT1,ADT2,..."

# Load necessary libraries
library(Seurat)

#Purpose: if you would like to ensure you have good sensitivity for your HTO/ADT libraries or if have mixed up/forgotten which of your HTO/ADT libraries belong to a given gene expression library, this will help.

# this assumes that you have you gene expression libraries in a folder named "starsolo_out", and your HTO and ADT libraries in a folder named citeseq_counts
# the script will give you a matrix for each of your gene expression libraries. The matrix will show you the number of intersecting barcodes for each combination of HTO and ADT libraries 


# Function to check if files exist
check_files_exist <- function(paths) {
  for (path in paths) {
    if (!file.exists(path)) {
      return(FALSE)
    }
  }
  return(TRUE)
}

# Function to process an array
process_array <- function(array, hto_libs, adt_libs) {
  # Initialize a matrix to store the number of common barcodes
  common_barcode_matrix <- matrix(nrow = length(hto_libs), ncol = length(adt_libs), dimnames = list(hto_libs, adt_libs))

  # Load UMI matrix for the array
  mtx_path <- paste0("starsolo_out/", array, "/Solo.out/Gene/raw/matrix.mtx")
  features_path <- paste0("starsolo_out/", array, "/Solo.out/Gene/raw/features.tsv")
  cells_path <- paste0("starsolo_out/", array, "/Solo.out/Gene/raw/barcodes.tsv")
  umis_filtered <- ReadMtx(mtx = mtx_path, features = features_path, cells = cells_path)

  for (hto_lib in hto_libs) {
    for (adt_lib in adt_libs) {
      # Construct file paths for HTO and ADT
      hto_path <- paste0("citeseq_counts/", hto_lib, "/umi_count/")
      adt_path <- paste0("citeseq_counts/", adt_lib, "/umi_count/")

      # Check if HTO and ADT files exist
      if (!check_files_exist(c(hto_path, adt_path))) {
        next
      }

      # Load HTO and ADT reads matrix
      umis_hto <- Read10X(data.dir = hto_path, gene.column = 1)
      umis_adt <- Read10X(data.dir = adt_path, gene.column = 1)

      # Intersect barcodes
      joint_bcs_hto <- intersect(colnames(umis_filtered), colnames(umis_hto))
      joint_bcs_adt <- intersect(colnames(umis_filtered), colnames(umis_adt))
      common_barcodes <- intersect(joint_bcs_hto, joint_bcs_adt)

      # Store the number of common barcodes in the matrix
      common_barcode_matrix[hto_lib, adt_lib] <- length(common_barcodes)
    }
  }

  # Print or write the matrix to a file
  print(paste("Common Barcodes Matrix for Array:", array))
  print(common_barcode_matrix)
  
  # Optional: write to file
  # write.csv(common_barcode_matrix, paste0("common_barcodes_matrix_", array, ".csv"))
}

# Main Script

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript script_name.R \"array1,array2,...\" \"HTO1,HTO2,...\" \"ADT1,ADT2,...\"", call. = FALSE)
}

# Split the arguments into arrays, HTO libraries, and ADT libraries
arrays <- strsplit(args[1], ",")[[1]]
hto_libs <- strsplit(args[2], ",")[[1]]
adt_libs <- strsplit(args[3], ",")[[1]]

# Process each array
for (array in arrays) {
  process_array(array, hto_libs, adt_libs)
}
