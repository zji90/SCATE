#' Cell clustering
#'
#' Perform Cell Clustering
#'
#' This function generates averaged signals for CRE clusters and cluster cells.
#' @param satac If type='reads', satac should be a list of GRanges object of scATAC-seq reads. Each element corresponds to one single cell. The GRanges should be the middle point of the reads with length of 1 base pair. Use 'satacprocess' to preprocess raw reads. If type='peaks', satac should be a list of data frames of scATAC-seq peaks. For each data frame, first column is chromsome name, second column is start site, third column is end site, and fourth column is the number of reads of the peak.
#' @param type Character variable of either 'reads' or 'peaks'.
#' @param peakOverlapMethod Character variable of either 'full' or 'middle'. Only effective when type = 'peaks'. If peakOverlapMethod='full', then the full range of the peak will be used to find overlap with bins, and all bins overlapping with this peak will be assigned the read counts of this peak. If peakOverlapMethod='middle', only the middle base pair of the peak will be used to find overlap with bins.
#' @param genome Character variable of either "hg19" or "mm10".
#' @param clunum Numeric variable giving the number of clusters. If NULL the cluster number will be determined automatically
#' @param perplexity Numeric variable specifying perplexity for tSNE. Reduce perplexity when sample size is small.
#' @param filtervar If TRUE, filter out features with low variability.
#' @param datapath Character variable of the path to the customized database (eg myfolder/database.rds). The database can be made using 'makedatabase' function. If not null, 'genome' is ignored.
#' @return A list of three components: tsne results, clustering results and aggregated signal for CRE cluster.
#' @export
#' @import GenomicRanges mclust splines Rtsne preprocessCore
#' @author Zhicheng Ji, Weiqiang Zhou, Wenpin Hou, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' set.seed(12345)
#' celldata <- sapply(1:50,function(i) {pos <- sample(1:1e9,50000); GRanges(seqnames=sample(paste0("chr",1:20),50000,replace=T),IRanges(start=pos,end=pos))})
#' names(celldata) <- paste0('cell',1:50)
#' cellcluster(celldata,type='reads',genome="hg19",filtervar=FALSE,perplexity=1,clunum=3) # reads as input

cellcluster <- function(satac,type='reads',peakOverlapMethod = 'full',genome='hg19',clunum=NULL,perplexity=30,filtervar=TRUE,datapath=NULL) {
      if (!is.null(datapath)) {
            loaddata <- readRDS(datapath)
      } else {
            loaddata <- readRDS(paste0(system.file(package="SCATE"),"/extdata/",genome,".rds"))
      }
      allclunum <- loaddata$allclunum
      clu <- loaddata$cluster[,allclunum==5000]
      tabclu <- as.vector(table(clu))
      gr <- loaddata$gr
      id <- loaddata$id
      gr <- gr[id]
      clusteraggregate <- sapply(names(satac),function(sid) {
            if (type=='reads') {
                  len <- length(satac[[sid]])/1e8
                  rowsum(log2(countOverlaps(gr,satac[[sid]],ignore.strand=TRUE)/len+1),clu)/tabclu
            } else if (type=='peaks') {
                  len <- sum(satac[[sid]][,4])/1e8
                  if (peakOverlapMethod == 'full'){
                        peak <- GRanges(seqnames=satac[[sid]][,1],IRanges(start=satac[[sid]][,2],end=satac[[sid]][,3]))
                  } else {
                        tmp <- floor((satac[[sid]][,2]+satac[[sid]][,3])/2)
                        peak <- GRanges(seqnames=satac[[sid]][,1],IRanges(start=tmp,end=tmp))
                  }
                  o <- as.matrix(findOverlaps(gr,peak))
                  count <- rep(0,length(gr))
                  count[o[,1]] <- satac[[sid]][o[,2],4]
                  rowsum(log2(count/len+1),clu)/tabclu
            } else {
                  stop('Wrong type')
            }
      })
      
      row.names(clusteraggregate) <- 1:nrow(clusteraggregate)
      clusteraggregate <- clusteraggregate[rowMeans(clusteraggregate > 0) > 0.1,]
      dn <- dimnames(clusteraggregate)
      clusteraggregate <- normalize.quantiles(clusteraggregate)
      dimnames(clusteraggregate) <- dn
      if (filtervar) {
            rsd <- apply(clusteraggregate,1,sd)
            rm <- rowMeans(clusteraggregate)
            clusteraggregate <- clusteraggregate[rsd > fitted(lm(rsd~bs(rm))),,drop=F]
      }
      prres <- prcomp(t(clusteraggregate),scale. = T)$x
      
      tsne <- Rtsne(prres[,1:min(50,ncol(prres))],pca=F,perplexity=perplexity)$Y
      row.names(tsne) <- row.names(prres)
      
      if (is.null(clunum)) {
            cluster <- Mclust(tsne,G=1:20,prior = priorControl(),verbose=F)
      } else {
            cluster <- Mclust(tsne,G=clunum,prior = priorControl(),verbose=F)
      }
      cluster <- apply(cluster$z,1,which.max)
      list(tsne=tsne,cluster=cluster,clusteraggregate=clusteraggregate)
}

