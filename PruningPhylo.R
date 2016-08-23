#####################################################
########## PRUNING ALGORITHM ########################
#####################################################

# Load appropriate R libraries
library(ape)
library(picante)

# Upload megaphylogeny to global environment
tree <- read.tree(file.choose())
spp <- tree$tip.label

# NOTE: Hedtke phylo confirms each bee genera is a monophyletic clade 

# Use the new vector to prune the existing mega-phylogeny using prune.sample 
# From the species of interest, find the appropriate species. If not found, use the closest related species from existing phylogenies

# Genera of interest from our species pool
# Megaphylogeny contains all of our genera of interest 
gen <- c("Anthidium", "Anthophora", "Chelostoma", "Coelioxys", "Heriades", "Hoplitis", "Hylaeus", "Megachile", "Osmia", "Stelis")

# custom fucntion: prune.tree("D99_phylo.txt", gen)
prune.tree <- function (newick, genera) {
  # Read newick file into a phylo object 
  tree <- read.tree(newick)
  # Extract each genera (and its corresponding species) using REGEX into a vector 
  t <- unlist(lapply(gen, FUN = function(x) grep(paste0(x, "+"), tree$tip.label, perl = TRUE, value = TRUE)))
  # Create a pseudo community matrix (species names as columns)
  n <- rbind(t)
  colnames(n) <- n
  # Prune the megaphylogeny
  hi <- prune.sample(n, tree)
  # Export pruned tree into a newick file 
  write.tree(hi, "pruned_tree.new")
  print("created newick file! check working directory")
  # return an image of the pruned phylogeny, with bootstrap support
  return(plot.phylo(hi, font = 1, cex = 0.5, show.node.label = TRUE))
}
