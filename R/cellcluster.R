#' Cell clustering
#'
#' Perform Cell Clustering
#'
#' This function generates averaged signals for CRE clusters and cluster cells.
#' @param satac GRanges object or list of GRanges object of scATAC-seq reads. The GRanges should be the middle point of the reads with length of 1 base pair. Use 'satacprocess' to preprocess raw reads.
#' @param genome Character variable of either "hg19" or "mm10".
#' @param clunum Numeric variable giving the number of clusters. If NULL the cluster number will be determined automatically
#' @param perplexity Numeric variable specifying perplexity for tSNE. Reduce perplexity when sample size is small.
#' @param datapath Character variable of the path to the customized database (eg myfolder/database.rds). The database can be made using 'makedatabase' function. If not null, 'genome' is ignored.
#' @return A list of three components: tsne results, clustering results and aggregated signal for CRE cluster.
#' @export
#' @import GenomicRanges mclust splines Rtsne preprocessCore
#' @author Zhicheng Ji, Weiqiang Zhou, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' cellcluster(list(cell1=GRanges(seqnames="chr1",IRanges(start=1:100+1e6,end=1:100+1e6)),cell2=GRanges(seqnames="chr2",IRanges(start=1:100+1e6,end=1:100+1e6))),genome="hg19")

cellcluster <- function(satac,genome,clunum=NULL,perplexity=30,datapath=NULL) {
      set.seed(12345)
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
            len <- length(satac[[sid]])/1e8
            rowsum(log2(countOverlaps(gr,satac[[sid]],ignore.strand=T)/len+1),clu)/tabclu
      })
      
      row.names(clusteraggregate) <- 1:nrow(clusteraggregate)
      clusteraggregate <- clusteraggregate[rowMeans(clusteraggregate > 0) > 0.1,]
      dn <- dimnames(clusteraggregate)
      clusteraggregate <- normalize.quantiles(clusteraggregate)
      dimnames(clusteraggregate) <- dn
      rsd <- apply(clusteraggregate,1,sd)
      rm <- rowMeans(clusteraggregate)
      clusteraggregate <- clusteraggregate[rsd > fitted(lm(rsd~bs(rm))),]
      
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

