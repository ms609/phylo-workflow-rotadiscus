# Print current working directory, which should contain the scripts,
# matrix and trees.
getwd()

# If getwd() does not contain the relevant files, set wd to working directory
wd <- "./"
outgroup <- c("") # Specify taxa on which to root tree

source(paste0(wd, "../R/common.R"))
source(paste0(wd, "../R/plot.R"))

treeFiles <- list.files(wd, "*.tre")
trees <- lapply(treeFiles, ReadTntTree)
if (outgroup == "") {
  outgroup <- TipLabels(trees[[1]])[1]
}

for (i in seq_along(treeFiles)) {
  tr <- trees[[i]]
  
  # Ignore outgroup taxa that aren't in tree
  og <- intersect(outgroup, TipLabels(tr)[[1]])
  if (length(og)) {
    # Root trees on outgroup
    tr <- RootTree(tr, og)
  }
  rogues <- Rogue::QuickRogue(tr, p = 1)
  cons <- SortTree(ConsensusWithout(tr, rogues[-1, "taxon"]))
  
  pdf(gsub(".tre", ".pdf", treeFiles[i], fixed = TRUE), 
      width = 8, height = 10)
  
  ColPlot(cons, ec = "black")
  if (nrow(rogues) > 1) {
    legend("topleft", rogues[-1, "taxon"], bty = "n", lty = 2)
  }
  k <- KValue(treeFiles[i])
  legend("topright", treeFiles[i], bty = "n")
  
  if (length(tr) > 5) {
    distances <- TreeDist::ClusteringInfoDistance(tr)
    map <- cmdscale(distances, k = 3)
    
    # Prepare plotting area
    par(mar = rep(0, 4))
    plot(map, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)
    
    # Add minimum spanning tree
    TreeTools::MSTEdges(distances, plot = TRUE, map[, 1], map[, 2],
                        col = '#00000030', lty = 2)
    
    # Connect trees by order found
    lines(map[, 1], map[, 2], col = "#ffccaa", lty = 1)
    
    # Add points
    TreeDist::Plot3(map,
    #                col = treeCols,
                    pch = 16, cex = 2,
                    add = TRUE)
    
    # Add legends
    legend("topleft", 
           c("Minimum spanning tree (mapping distortion)",
             "Order in which trees found"),
           lty = c(2, 1),
           col = c("#00000030", "#ffccaaaa"),
           bty = "n")
  }
  dev.off()
}

message(" # # # Complete # # # ")
