library(XML)
library(rentrez)
library(data.table)

# calling functions
xml_parse <- function(query, silent = F, retmax = 10) {
  # querying the database to retrieve a list of unique ID's for each hit
  # n_terms indicates maximum number of hits returned by the search
  Tt <- entrez_search(db="nuccore", term=query, retmax = retmax)
  if (length(Tt$ids) == 0) {
    if (!silent) cat("No results for", query, "\n")
    return(paste0("[",query, "]"))
  }
  return(Tt$ids)
} 

read.seq <- function (idx) {
  
  tr <- entrez_fetch(db="nuccore", id=idx, rettype="xml", parsed=T)
  tr1 <- getNodeSet(tr, paste0("//Seq-entry[Seq-entry_seq/Bioseq/Bioseq_id/Seq-id=",idx,"]"),fun=xmlClone)
  ls_acc <- xpathApply(tr1[[1]], "//Textseq-id_accession", xmlValue)[[1]]
  ls_acc_ver <- xpathApply(tr1[[1]], "//Textseq-id_version", xmlValue)[[1]]
  ls_seqs <- xpathApply(tr1[[1]], "//IUPACna", xmlValue)[[1]]
  ls_titles <- xpathApply(tr1[[1]], "//Seqdesc_title", xmlValue)[[1]]
  ls_org <- try(xpathApply(tr1[[1]], "//Org-ref_taxname", xmlValue)[[1]], T)
  if (class(ls_org) == "try-error") ls_org <- NA_character_
  ls_lineage <- try(xpathApply(tr1[[1]], "//OrgName_lineage", xmlValue)[[1]], T)
  if (class(ls_lineage) == "try-error") {
    ls_lineage <- NA_character_
    ls_kingdom <- NA_character_
  } else {
    ls_kingdom <- strsplit(ls_lineage, "; ")[[1]][2]
  }
  ls_names <- c(xpathApply(tr1[[1]], "//descendant::SeqFeatData", xmlValue), recursive = T)
  ls_from <- as.integer(xpathApply(tr1[[1]],"//Seq-interval_from", xmlValue)) + 1
  ls_to <- as.integer(xpathSApply(tr1[[1]], "//Seq-interval_to", xmlValue)) + 1

  if (any(c(length(ls_seqs)    != 1, 
            length(ls_titles)  != 1, 
            length(ls_names)   != length(ls_from), 
            length(ls_names)   != length(ls_to), 
            length(ls_names)   < 1, 
            length(ls_acc)     != 1,
            length(ls_acc_ver) != 1, 
            length(ls_org)     != 1,
            length(ls_kingdom) != 1))) {
    
    stop(paste0("Error in seq ", idx, 
                "(", 
                length(ls_seqs), 
                length(ls_titles), 
                length(ls_names),
                length(ls_from), 
                length(ls_to), 
                length(ls_org),
                length(ls_kingdom),")"
                )
    )
  }
  
  return(data.table(kingdom         = ls_kingdom,
                    spp             = ls_org,
                    length          = nchar(ls_seqs),
                    title           = ls_titles,
                    accession       = paste0(ls_acc, ".", ls_acc_ver),
                    name            = ls_names, 
                    from            = as.integer(ls_from), 
                    to              = as.integer(ls_to),
                    whole.sequence  = ls_seqs))
  
}

scrape.gen <- function(spp, sequences, retmax = 10) {
  
  queryList <- as.character(outer(paste(spp, "[ORGN]"), sequences, paste))
  Tt <- c(lapply(queryList, xml_parse, retmax = retmax), recursive = T)
  lerrors <- grepl("^\\[", Tt)
  errors <- Tt[lerrors]
  Tt <- Tt[!lerrors]
  if (length(Tt) > 0) {
    lres <- lapply(Tt, read.seq)
    df.res <- rbindlist(lres, fill = T)
    df.res$sub.sequence <- substr(df.res$whole.sequence, df.res$from, df.res$to)
    df.res$sub.seq.length <- (df.res$to - df.res$from) + 1 
    return(list(df = df.res, errors = errors))
  }
  return(NULL)
}

