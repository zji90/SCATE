#' Preprocess scATAC-seq
#'
#' Preprocessing scATAC-seq samples
#'
#' This function filters out scATAC-seq with low library size and transforms the reads into middle points of the reads.
#' @param input Either a character vector of locations to bam files (when type is bam) or list of GRanges object of scATAC-seq reads (when type is gr).
#' @param type Either 'bam' or 'gr'. 'bam' if the input is locations to bam files. 'gr' if the input is a list of GRanges object.
#' @param pairsingle Either 'paired' or 'single'. Specify whether the bam file is single- or paired-end. Only works when type='bam'.
#' @param libsizefilter Numeric variable giving the minimum library size. scATAC-seq samples with library size smaller than this cutoff will be discarded.
#' @return GRanges object of list of GRanges object after preprocessing.
#' @export
#' @import GenomicRanges
#' @author Zhicheng Ji, Weiqiang Zhou, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' satacprocess(list(cell1=GRanges(seqnames="chr1",IRanges(start=1:100+1e6,end=1:100+1e6)),cell2=GRanges(seqnames="chr2",IRanges(start=1:100+1e6,end=1:100+1e6))),type='gr',libsizefilter=10)

satacprocess <- function(input,type='bam',pairsingle='paired',libsizefilter=1000) {
      if (type=='bam') {
            if (pairsingle=='paired') {
                  satac <- sapply(sapply(input,readGAlignmentPairs),GRanges,on.discordant.seqnames="drop")
            } else {
                  satac <- sapply(sapply(input,readGAlignments),GRanges)      
            }
      } else {
            satac <- input
      }
      satac <- satac[sapply(satac,length) >= libsizefilter]
      n <- names(satac)
      satac <- lapply(satac,function(i) {
            start(i) <- end(i) <- round((start(i) + end(i))/2)
            i
      })
      names(satac) <- n
      satac
}

