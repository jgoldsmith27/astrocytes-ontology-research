# Load required libraries
library(Matrix)

# Load the data files
load("data/cellchat/raw/CellChatDB.human.rda")
print("CellChatDB.human contents:")
print(ls())
print("\nStructure of CellChatDB.human:")
str(CellChatDB.human)

# Export each component of CellChatDB.human
if (exists("CellChatDB.human")) {
  if (!dir.exists("data/cellchat/processed")) {
    dir.create("data/cellchat/processed", recursive = TRUE)
  }
  
  # Export interaction data
  write.csv(CellChatDB.human$interaction, "data/cellchat/processed/interaction.csv", row.names = FALSE)
  
  # Export complex data
  write.csv(CellChatDB.human$complex, "data/cellchat/processed/complex.csv", row.names = FALSE)
  
  # Export cofactor data
  write.csv(CellChatDB.human$cofactor, "data/cellchat/processed/cofactor.csv", row.names = FALSE)
  
  # Export gene info data
  write.csv(CellChatDB.human$geneInfo, "data/cellchat/processed/gene_info.csv", row.names = FALSE)
}

# Try to load PPI data
tryCatch({
  load("data/cellchat/raw/PPI.human.rda")
  print("\nPPI.human contents:")
  print(ls())
  print("\nStructure of PPI.human:")
  str(PPI.human)
  
  # Export PPI data as edge list
  if (exists("PPI.human")) {
    # Convert sparse matrix to edge list
    edges <- summary(PPI.human)
    edge_df <- data.frame(
      protein1 = rownames(PPI.human)[edges$i],
      protein2 = colnames(PPI.human)[edges$j],
      weight = edges$x
    )
    
    # Export edge list
    write.csv(edge_df, "data/cellchat/processed/ppi.csv", row.names = FALSE)
    
    # Print summary
    print("\nPPI network summary:")
    print(paste("Number of proteins:", nrow(PPI.human)))
    print(paste("Number of interactions:", nrow(edge_df)))
  }
}, error = function(e) {
  print(paste("Error loading PPI data:", e$message))
}) 