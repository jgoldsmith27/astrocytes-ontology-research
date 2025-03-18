#!/usr/bin/env Rscript

# Script to explore PPI.human.rda data

# Set working directory
setwd("/Users/jacob/Desktop/astrocytes-ontology-research/cell-chat-data")

# Load the data
load("PPI.human.rda")

# Basic examination
cat("\n===== PPI.human Basic Structure =====\n")
cat("Class of PPI.human:", class(PPI.human), "\n")
cat("Object type:", typeof(PPI.human), "\n")
cat("Object mode:", mode(PPI.human), "\n")
cat("Object length:", length(PPI.human), "\n")

# Try to determine the structure of the object
cat("\n===== PPI.human Content Structure =====\n")

if (is.list(PPI.human)) {
  cat("PPI.human is a list with", length(PPI.human), "elements\n")
  cat("List names:", names(PPI.human), "\n")
} else if (is.data.frame(PPI.human)) {
  cat("PPI.human is a data frame with dimensions:", dim(PPI.human), "\n")
  cat("Column names:", names(PPI.human), "\n")
  cat("Sample rows:\n")
  print(head(PPI.human, 3))
} else if (is.matrix(PPI.human)) {
  cat("PPI.human is a matrix with dimensions:", dim(PPI.human), "\n")
  cat("Sample rows/cols:\n")
  print(PPI.human[1:3, 1:3])
} else {
  cat("PPI.human is a different type of object\n")
  # Try to dump the structure
  cat("Structure summary:\n")
  str(PPI.human)
}

# If it's a matrix or data frame, check dimensions
if (is.matrix(PPI.human) || is.data.frame(PPI.human)) {
  cat("\nDimensions:", dim(PPI.human), "\n")
}

# Try to examine the first few rows or elements
cat("\n===== PPI.human Content Preview =====\n")
tryCatch({
  if (is.list(PPI.human) && length(PPI.human) > 0) {
    cat("First element preview:\n")
    print(head(PPI.human[[1]], 3))
  } else {
    cat("First rows preview:\n")
    print(head(PPI.human, 3))
  }
}, error = function(e) {
  cat("Error previewing data:", e$message, "\n")
})

# Try to list all values if it's small enough
cat("\n===== Full PPI.human Content (if small enough) =====\n")
tryCatch({
  if (length(PPI.human) < 100) {
    print(PPI.human)
  } else {
    cat("Object too large to display completely\n")
  }
}, error = function(e) {
  cat("Error displaying full data:", e$message, "\n")
})

# Create network statistics if it represents a network
cat("\n===== Network Analysis (if applicable) =====\n")
tryCatch({
  if (is.matrix(PPI.human) && nrow(PPI.human) == ncol(PPI.human)) {
    cat("Matrix appears to be a square adjacency matrix with", nrow(PPI.human), "nodes\n")
    # Count non-zero entries which represent edges
    edge_count <- sum(PPI.human != 0)
    cat("Number of edges (non-zero entries):", edge_count, "\n")
  } else if (is.data.frame(PPI.human) && "from" %in% names(PPI.human) && "to" %in% names(PPI.human)) {
    cat("Data frame appears to be an edge list with", nrow(PPI.human), "edges\n")
    # Count unique nodes
    nodes <- unique(c(PPI.human$from, PPI.human$to))
    cat("Number of unique nodes:", length(nodes), "\n")
  }
}, error = function(e) {
  cat("Error analyzing network structure:", e$message, "\n")
})

# Save a summary to a text file
sink("ppi_summary.txt")
cat("===== PPI.human Summary =====\n")
cat("Class:", class(PPI.human), "\n")
cat("Type:", typeof(PPI.human), "\n")
if (is.list(PPI.human)) {
  cat("List length:", length(PPI.human), "\n")
} else if (is.matrix(PPI.human) || is.data.frame(PPI.human)) {
  cat("Dimensions:", dim(PPI.human), "\n")
}
sink()

cat("\nSummary saved to ppi_summary.txt\n") 