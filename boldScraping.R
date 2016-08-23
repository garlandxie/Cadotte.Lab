#####################################################
############## Install Packages #####################
#####################################################

# BOLD Package
install.packages("bold")

# MSA Package
source("http://www.bioconductor.org/biocLite.R")
biocLite("msa")

# Load appropriate libraries
library(bold)
library(msa)

######################################################
############### BOLD SCRAPING ########################
######################################################

# If you have your own list of taxon ID's from BOLD
sample <- c("JSHYN895-11","BCHYM1500-13", "HDBNS068-03","HCBNS158-03")

# Gather sequence data from BOLD for the corresponding taxon ID's 
# Concatenate taxon name and sequence data into a single fasta file
bold.fasta <- function (sample, file.name) {
  seqs <- bold_seq(ids = sample)
  tt <- unlist(lapply(seqs, function(y) paste(">", gsub(" ", "_", y$name), "\n", y$sequence, "\n")))
  write(tt, file = file.name)
  print("done!, check the fasta file in your working directory")
} 

#####################################################
########## MULTIPLE SEQUENCE ALIGNMENT ##############
#####################################################

# From an existing fasta file, create a multiple sequence alignment.
# Avoids using external software and streamlines bioinformatic pipeline within R.

align.Fasta <- function (file, algo) {
mySeqfile <- readAAStringSet(file, format="fasta")
myFirstAlign <- msa(mySeqfile, method = algo)
return (print(myFirstAlign, show="complete"))
}
