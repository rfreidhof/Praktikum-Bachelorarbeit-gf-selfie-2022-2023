# converts specifically a featureCounts result to TPM. The fiven file is the featureCounts table.

library(readr)

#get commandline arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "./out.tsv"
}

featureCounts = read_tsv(args[1], col_names = TRUE, comment = '#')

counts = as.data.frame(featureCounts[,c(7)])
featureLength = featureCounts[,6]
geneids = as.data.frame(featureCounts[,c(1)])
meanFragmentLength = mean(as.matrix(featureLength))
featureLength = t(featureLength)


#' Source: https://gist.github.com/slowkow/c6ab0348747f86e2748b
#' Slightly changed by me.
#' Convert counts to transcripts per million (TPM).
#' 
#' Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to transcripts per million.
#' 
#'    Lior Pachter. Models for transcript quantification from RNA-Seq.
#'    arXiv:1104.3889v2 
#'    
#'    Wagner, et al. Measurement of mRNA abundance using RNA-seq data:
#'    RPKM measure is inconsistent among samples. Theory Biosci. 24 July 2012.
#'    doi:10.1007/s12064-012-0162-3
#'    
#' @param counts A numeric matrix of raw feature counts i.e.
#'  fragments assigned to each gene.
#' @param featureLength A numeric vector with feature lengths.
#' @param meanFragmentLength A numeric vector with mean fragment lengths.
#' @return tpm A numeric matrix normalized by library size and feature length.
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  featureLength = t(featureLength)
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- as.data.frame(counts[idx,])
  effLen <- as.data.frame(effLen[idx,])
  featureLength <- featureLength[idx]
  geneids <- as.data.frame(geneids[idx,])
  

  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix. Deprecated. Add the gene ids and restore names.
  tpm <- cbind(geneids,tpm)
  colnames(tpm) <- colnames(featureCounts[,c(1,7)])
  return(tpm)
}

tpm=counts_to_tpm(counts, featureLength, meanFragmentLength)

write.table(tpm, file=args[2], quote=FALSE, sep='\t', col.names = NA)
