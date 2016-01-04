library(XML)
library(rentrez)
library(data.table)
# querying the database
# getting the id's for the accessions


xml_parse <- function(query, n_terms) {
  q <- toString(query)
  Tt <- entrez_search(db="nuccore", term=q, retmax = n_terms)
  if (length(Tt$ids) == 0) {
    stop("No items were found, Please try again")
  }
  Tt
} 

read.seq <- function (idx) {
  tr <- entrez_fetch(db="nuccore", id=Tt$ids, rettype="xml", parsed=T)
  tr1 <- getNodeSet(tr, paste0("//Seq-entry[Seq-entry_seq/Bioseq/Bioseq_id/Seq-id=",idx,"]"),fun=xmlClone)
  ls_acc <- xpathApply(tr1[[1]], "//Textseq-id_accession", xmlValue)[[1]]
  ls_acc_ver <- xpathApply(tr1[[1]], "//Textseq-id_version", xmlValue)[[1]]
  ls_seqs <- xpathApply(tr1[[1]], "//IUPACna", xmlValue)[[1]]
  ls_titles <- xpathApply(tr1[[1]], "//Seqdesc_title", xmlValue)[[1]]
  ls_names <- xpathApply(tr1[[1]], "//descendant::SeqFeatData", xmlValue)
  ls_from <- as.integer(xpathApply(tr1[[1]],"//Seq-interval_from", xmlValue))+1
  ls_to <- as.integer(xpathSApply(tr1[[1]], "//Seq-interval_to", xmlValue))+1
  if (any(c(length(ls_seqs)  != 1, 
            length(ls_titles)!= 1, 
            length(ls_names) != length(ls_from), 
            length(ls_names) != length(ls_to), 
            length(ls_names) < 1, 
            length(ls_acc)   != 1))) {
    stop(paste0("Error in seq ", idx, " (", length(ls_seqs), length(ls_titles), 
                length(ls_names),length(ls_from), length(ls_to), ")"))
  }
  return(data.frame(whole.sequence  = toString(ls_seqs),
                    length    = nchar(ls_seqs),
                    title     = toString(ls_titles),
                    accession = toString(paste0(ls_acc,".",ls_acc_ver)),
                    name      = toString(ls_names), 
                    from      = as.integer(ls_from), 
                    to        = as.integer(ls_to)))
}

lres <- lapply(Tt$ids, read.seq)
df.res <- rbindlist(lres, fill = T)
df.res$sub.sequence <- apply(df.res, 1, function(x) substr(df.res$whole.sequence, df.res$from, df.res$to))
df.res$sub.seq.length <- nchar(df.res$sub.sequence)
