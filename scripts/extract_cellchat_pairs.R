#!/usr/bin/env Rscript

# Script to extract all ligand-receptor pairs from CellChatDB.human.rda

# Set working directory to the project root
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  setwd(args[1])
} else {
  # Default to the current directory
  setwd(getwd())
}

# Load the CellChat database
cat("Loading CellChatDB.human.rda...\n")
load("cell-chat-data/CellChatDB.human.rda")

# Check if data was loaded correctly
if (!exists("CellChatDB.human")) {
  stop("Failed to load CellChatDB.human variable")
}

# Create a directory for extracted data if it doesn't exist
output_dir <- "output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Extract the interactions data frame
cat("Extracting interaction data...\n")
interactions <- CellChatDB.human$interaction

# Check the structure and content
cat("Database contains", nrow(interactions), "interaction records\n")
cat("Columns in the interaction data:", paste(names(interactions), collapse=", "), "\n")

# Create an expanded list of ligand-receptor pairs
# Some entries may have multiple ligands or receptors separated by "_"
lr_pairs <- list()
pair_count <- 0

cat("Extracting individual ligand-receptor pairs...\n")
for (i in 1:nrow(interactions)) {
  row <- interactions[i,]
  
  # Get the ligand and receptor entries
  ligand_entry <- row$ligand
  receptor_entry <- row$receptor
  
  # Split into individual ligands and receptors if they contain "_"
  ligands <- unlist(strsplit(ligand_entry, "_"))
  receptors <- unlist(strsplit(receptor_entry, "_"))
  
  # Create all possible combinations
  for (lig in ligands) {
    for (rec in receptors) {
      pair_count <- pair_count + 1
      lr_pairs[[pair_count]] <- list(
        ligand = lig,
        receptor = rec,
        interaction_name = row$interaction_name,
        pathway = row$pathway_name,
        annotation = as.character(row$annotation)
      )
    }
  }
}

# Convert the list to a data frame
cat("Creating data frame from", length(lr_pairs), "ligand-receptor pairs...\n")
lr_df <- do.call(rbind, lapply(lr_pairs, as.data.frame))

# Write the data to a CSV file
output_file <- file.path(output_dir, "cellchat_lr_pairs.csv")
cat("Writing ligand-receptor pairs to", output_file, "...\n")
write.csv(lr_df, output_file, row.names = FALSE)

# Output some statistics
cat("\nSummary Statistics:\n")
cat("Total ligand-receptor pairs extracted:", nrow(lr_df), "\n")
cat("Unique ligands:", length(unique(lr_df$ligand)), "\n")
cat("Unique receptors:", length(unique(lr_df$receptor)), "\n")

# Display a few example pairs
cat("\nExample ligand-receptor pairs:\n")
print(head(lr_df, 10))

cat("\nExtraction complete. File saved to", output_file, "\n") 