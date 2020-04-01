read_kb_counts <- function(dir, name) {
  # Loading scRNA-seq count matrix
  #
  # Generates a combined scRNA-seq matrix from the ouput of 
  # the Kallisto | Bustools wrapper. Needs the location of
  # the count tables to combined. Returned matrix consists of 
  # genes (with gene names) in the rows, and cells (with 
  # barcode identifier) in the columns. 
  #
  # dir: location of the "_output" folders generated 
  #   by the kb-wrapper
  # name: the name of count assay you want to load (e.g. 
  #   spliced, unspliced or cells_x_genes, the latter when 
  #   running without velocity)
  
  ## Loading packages & files ##
  library(Matrix)
  library(tidyr)
  library(dplyr)
  # Location of barcode file
  barcode_file <- "/home/snabel/scrna/barcode_384.tab"
  
  ## Generate matrix ##
  
  # Iterate over all _output folders, generating matrix per plate
  # combining matrices by gene matches
  dir <- normalizePath(dir, mustWork = TRUE)
  output_folders <- list.files(dir, pattern = "_output", 
                               recursive = TRUE, include.dirs = TRUE)
  i <- 1
  for (folder in output_folders){
    print(paste("Reading:",folder))
    plate <- paste0(dir, "/", folder, "/counts_unfiltered/", name)
    m <- readMM(paste0(plate, ".mtx"))
    m <- Matrix::t(m)
    m <- as(m, "dgCMatrix")
    # the matrix has genes in columns, cells in rows, 
    # stored in .genes.txt and .barcodes.txt
    genes <- as.vector(read.table(file(paste0(plate, ".genes.txt")))[,1])
    barcodes <- as.vector(read.table(file(paste0(plate, ".barcodes.txt")))[,1])
    # retrieve unique plate-id from folder
    platename <- gsub("_output", "", folder)
    colnames(m) <- paste(barcodes, platename, sep = "_")
    rownames(m) <- genes
    # create a combined matrix for all plates in the folder
    if (i == 1) {
      combined <- m
    } else if (identical(rownames(combined),rownames(m)) == TRUE){
      # Only binds the matrices if genes are identical and in the same order
      combined <- cbind(combined, m)
    }
    i <- i + 1
  }
  
  ## Replace gene identifiers with gene names ##
  
  t2g <- read.table("/home/snabel/refs/genomes/GRCh38/kb_genome/erccrep_index/GRCh38_transcripts_to_genes.txt", col.names = c("transcript", "gene", "gene_symbol")) %>%
    select(-transcript) %>%
    distinct()
  rownames(combined) <- t2g$gene_symbol[match(rownames(combined), t2g$gene)]
  
  ## Replace cell barcodes for well identifier ##
  
  # barcode file contains the well identifier and corresponding DNA barcode
  plate_order <- read.table(barcode_file, sep = "\t", col.names = c("well","barcode"))
  # generate a data.frame to match barcode and wellid 
  cells <- data.frame("cell" = colnames(combined))
  cells$barcode <- gsub("_.*", "", cells$cell)
  cells$well <- plate_order$well[match(cells$barcode, plate_order$barcode)]
  # Remove DNA barcode and add wellid 
  cells$cell_id <- paste(gsub("^.*?_", "", cells$cell), cells$well, sep = "_")
  cells$cell_id <- gsub("-", "_", cells$cell_id)
  # replace cell names of the count matrix
  colnames(combined) <- cells$cell_id
  
  return(combined)
}