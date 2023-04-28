library("TreeTools")
dataset <- "Rotadiscus_Data_S2.nex"
burninF <- 0.25
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

# Reduce to key relationships
myFive <- c("Cnidaria", "Vertebrata", "Rotadiscus", "Echinoidea", "Saccoglossus")
shrunk <- KeepTip(burntOff, myFive) |>
  RenumberTips(myFive)
topols <- vapply(as.TreeNumber(shrunk), as.integer, integer(1))
probs <- 100 * table(topols, dnn = NULL) / length(topols)
paste(vapply(
  names(probs),
  function(x) write.tree(as.phylo(as.integer(x), tipLabels = myFive)),
  ""), signif(probs), "%"
)

text(barplot(probs), probs, labels = signif(probs), xpd = NA, pos = 3)

