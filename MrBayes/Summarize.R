library("TreeTools")
dataset <- "Rotadiscus_Data_S2.nex"
burninF <- 0.2
outgroup <- "Cnidaria"

escapeDataset <- gsub("([\\-\\.])", "\\\\\\1", dataset)
mbFiles <- list.files(pattern = paste0(escapeDataset, "\\.run\\d+\\.t$"),
                      full.names = TRUE)

mbTrees <- setNames(
  lapply(mbFiles,
         function (mbFile) tryCatch(
           read.nexus(mbFile),
           error = function (err) if (err$message == "NA/NaN argument") {
             # Unterminated tree block, perhaps because a search is ongoing
             withEnd <- tempfile()
             on.exit(unlink(withEnd))
             writeLines(c(readLines(mbFile), "\nEND;"), withEnd)
             read.nexus(withEnd)
           })),
  paste0("mb", seq_along(mbFiles))
)

# Remove burnin
burnOff <- lapply(mbTrees, function (trees) {
  nTree <- length(trees)
  trees[(nTree * burninF):nTree]
})
burntOff <- RootTree(do.call(c, burnOff), outgroup)
cons <- Consensus(burntOff, 0.5)

oPar <- par(mar = rep(0, 4), cex = 0.8)
plot(cons)
par(oPar)

Summary <- function(trees) {
  # Reduce to key relationships
  myFive <- c("Cnidaria", "Vertebrata", "Rotadiscus",
              "Echinoidea", "Saccoglossus")
  shrunk <- KeepTip(trees, myFive) |>
    RenumberTips(myFive)
  topols <- vapply(as.TreeNumber(shrunk), as.integer, integer(1))
  probs <- 100 * table(topols, dnn = NULL) / length(topols)
  cat(paste(vapply(
    names(probs),
    function(x) write.tree(as.phylo(as.integer(x), tipLabels = myFive)),
    ""), signif(probs), "%\n"
  ))
  
  text(barplot(probs), probs, labels = signif(probs), xpd = NA, pos = 3)
}

Summary(burntOff)
