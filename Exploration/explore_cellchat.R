#!/usr/bin/env Rscript

# Script to explore CellChat data

# Set working directory
setwd("/Users/jacob/Desktop/astrocytes-ontology-research/cell-chat-data")

# Load the data
load("CellChatDB.human.rda")
load("PPI.human.rda")

# Print summary information about the CellChatDB
cat("\n===== CellChatDB.human Summary =====\n")
print(names(CellChatDB.human))

# Explore the structure of CellChatDB
cat("\n===== Structure of CellChatDB.human =====\n")
for (name in names(CellChatDB.human)) {
  cat("\nComponent:", name, "\n")
  if (is.data.frame(CellChatDB.human[[name]])) {
    cat("  Dimensions:", dim(CellChatDB.human[[name]]), "\n")
    cat("  Columns:", names(CellChatDB.human[[name]]), "\n")
    cat("  Sample rows:\n")
    print(head(CellChatDB.human[[name]], 3))
  } else {
    cat("  Type:", class(CellChatDB.human[[name]]), "\n")
  }
}

# Print summary information about the PPI data
cat("\n===== PPI.human Summary =====\n")
print(names(PPI.human))

# Explore the structure of PPI.human
cat("\n===== Structure of PPI.human =====\n")
for (name in names(PPI.human)) {
  cat("\nComponent:", name, "\n")
  if (is.data.frame(PPI.human[[name]])) {
    cat("  Dimensions:", dim(PPI.human[[name]]), "\n")
    cat("  Columns:", names(PPI.human[[name]]), "\n")
    cat("  Sample rows:\n")
    print(head(PPI.human[[name]], 3))
  } else {
    cat("  Type:", class(PPI.human[[name]]), "\n")
  }
}

# Extract ligand-receptor pairs
cat("\n===== Ligand-Receptor Pairs =====\n")
if ("interaction" %in% names(CellChatDB.human)) {
  interactions <- CellChatDB.human$interaction
  cat("Total number of interactions:", nrow(interactions), "\n")
  
  # Count unique ligands and receptors
  if ("ligand" %in% names(interactions) && "receptor" %in% names(interactions)) {
    unique_ligands <- unique(interactions$ligand)
    unique_receptors <- unique(interactions$receptor)
    cat("Number of unique ligands:", length(unique_ligands), "\n")
    cat("Number of unique receptors:", length(unique_receptors), "\n")
    
    # Show a few examples
    cat("\nSample ligands:", head(unique_ligands, 10), "\n")
    cat("\nSample receptors:", head(unique_receptors, 10), "\n")
    
    # Show a few example interactions
    cat("\nSample interactions:\n")
    sample_interactions <- head(interactions[, c("ligand", "receptor", "annotation")], 10)
    print(sample_interactions)
  }
}

# Save a summary to a text file
sink("cellchat_summary.txt")
cat("===== CellChatDB.human Summary =====\n")
print(names(CellChatDB.human))
if ("interaction" %in% names(CellChatDB.human)) {
  interactions <- CellChatDB.human$interaction
  cat("\nTotal interactions:", nrow(interactions), "\n")
  cat("Interaction columns:", names(interactions), "\n")
}
cat("\n===== PPI.human Summary =====\n")
print(names(PPI.human))
if (length(names(PPI.human)) > 0 && is.data.frame(PPI.human[[1]])) {
  cat("\nPPI dimensions:", dim(PPI.human[[1]]), "\n")
  cat("PPI columns:", names(PPI.human[[1]]), "\n")
}
sink()

cat("\nSummary saved to cellchat_summary.txt\n") 