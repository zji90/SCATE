#' Make customized database
#'
#' Make customized database with new bulk DNase-seq data
#'
#' This function makes a new customized database if users have new bulk DNase-seq data and such information can be contribued to the model building of SCATE.
#' @param datapath path to the data package folder (e.g. myfolder/hg19/). User must first download the data package to use this function. The data package for hg19 and mm10 can be downloaded from http://jilab.biostat.jhsph.edu/projects/scate/hg19.zip or http://jilab.biostat.jhsph.edu/projects/scate/mm10.zip. The compressed file should be unzipped. If users do not want to use existing data compendium (e.g. to build a database in a new species), datapath should be set NULL, and 'genome' will be ignored.
#' @param savepath path to save the generated database. e.g. myfolder/database.rds.
#' @param bamfile location of bulk DNase-seq bamfiles.
#' @param cre dataframe of new CRE sites to be added to the database. First column: chromosome name. Second column: start position. Third column: end position.
#' @param genome Character variable of either "hg19" or "mm10". Default is 'hg19'. Ignored when datapath is NULL.
#' @param genomerange Data frame with two columns. First column is the chromosome and second column is the length of the genome. Only useful when datapath is NULL. Example is https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
#' @export
#' @import GenomicAlignments
#' @author Zhicheng Ji, Weiqiang Zhou, Hongkai Ji <zji4@@zji4.edu>

makedatabase <- function(datapath,savepath,bamfile=NULL,cre=NULL,genome='hg19',genomerange=NULL) {
      
      suppressMessages(library(GenomicAlignments))
      
      if (is.null(datapath)) {
            
            chrlen <- genomerange
            allres <- NULL
            #sort chromosome
            for (i in 1:nrow(chrlen)) {
                  #keep all chromosome
                  start <- 200 * (0:floor((chrlen[i,2])/200))
                  end <- start + 199
                  allres <- rbind(allres,data.frame(chr=chrlen[i,1],start=start,end=end))
            }
            gr <- GRanges(seqnames=allres$chr,IRanges(start=allres$start,end=allres$end))
            over <- as.matrix(findOverlaps(gr,blacklist))[,1]
            gr <- gr[setdiff(1:length(gr),over),]

            if (!is.null(cre)) {
                  cre <- GRanges(seqnames=cre[,1],IRanges(start=cre[,2],end=cre[,3]))
                  creid <- as.matrix(findOverlaps(gr,cre))[,1]	
            } else {
                  creid <- NULL
            }
            
            len <- rep(0,length(bamfile))
            names(len) <- bamfile
            count <- matrix(0,nrow=length(gr),ncol=length(bamfile),dimnames = list(NULL,bamfile))
            for (f in bamfile) {
                  tmp <- readGAlignmentPairs(f)
                  if (length(tmp) == 0) tmp <- readGAlignments(f)
                  tmp <- GRanges(tmp)
                  mid <- round((start(tmp) + end(tmp))/2)
                  tmp <- GRanges(seqnames=as.character(seqnames(tmp)),IRanges(start=mid,end=mid))
                  count[,f] <- countOverlaps(gr,tmp)
                  len[f] <- length(tmp)
            }
            rm('tmp')
            len <- len/1e8
            norm <- log2(t(t(count)/len)+1)
            
            
            start <- end <- rep(0,length(gr))
            sn <- as.character(seqnames(gr))
            for (s in seqlevels(gr)) {
                  id <- which(sn==s)
                  start[id] <- pmax(min(id),id - 250)
                  end[id] <- pmin(max(id),id + 250)
            }
            
            pos <- cbind(start,end)
            
            newid <- lapply(1:ncol(count),function(i) {  
                  id <- which(count[,i] >= 10 & norm[,i] >= 5)
                  backnorm <- sapply(id,function(j) mean(norm[pos[j,1]:pos[j,2],i]))
                  id[which(norm[id,i]/backnorm >= 5)]	
            })
            
            id <- unique(c(unlist(newid),creid))
      } else {
            packdata <- readRDS(paste0(system.file(package="SCATE"),"/extdata/",genome,".rds"))
            gr <- packdata$gr
            
            if (!is.null(cre)) {
                  cre <- GRanges(seqnames=cre[,1],IRanges(start=cre[,2],end=cre[,3]))
                  creid <- as.matrix(findOverlaps(gr,cre))[,1]	
            } else {
                  creid <- NULL
            }
            
            len <- rep(0,length(bamfile))
            names(len) <- bamfile
            count <- matrix(0,nrow=length(gr),ncol=length(bamfile),dimnames = list(NULL,bamfile))
            for (f in bamfile) {
                  tmp <- readGAlignmentPairs(f)
                  if (length(tmp) == 0) tmp <- readGAlignments(f)
                  tmp <- GRanges(tmp)
                  mid <- round((start(tmp) + end(tmp))/2)
                  tmp <- GRanges(seqnames=as.character(seqnames(tmp)),IRanges(start=mid,end=mid))
                  count[,f] <- countOverlaps(gr,tmp)
                  len[f] <- length(tmp)
            }
            rm('tmp')
            len <- len/1e8
            norm <- log2(t(t(count)/len)+1)
            
            pos <- readRDS(paste0(datapath,"/backgrpos.rds"))
            newid <- lapply(1:ncol(count),function(i) {  
                  id <- which(count[,i] >= 10 & norm[,i] >= 5)
                  backnorm <- sapply(id,function(j) mean(norm[pos[j,1]:pos[j,2],i]))
                  id[which(norm[id,i]/backnorm >= 5)]	
            })
            
            id <- unique(c(packdata$id,unlist(newid),creid))
            
            dcount <- readRDS(paste0(datapath,'/countmat.rds'))
            dlen <- readRDS(paste0(datapath,'/len.rds'))/1e8
            dnorm <- log2(t(t(dcount)/dlen)+1)
            count <- cbind(count,dcount)
            len <- c(len,dlen)
            norm <- cbind(norm,dnorm)
            rm('dcount','dlen','dnorm')      
      }
      
      zid <- which(rowSums(norm) == 0)
      id <- setdiff(id,zid)
      excid <- setdiff(1:length(gr),c(id,zid))
      
      orisum <- sum1 <- sum2 <- sum3 <- rep(0,length(gr))
      for (f in colnames(count)) {
            orisum <- orisum + count[,f]
            sum1 <- sum1 + count[,f]/len[f]
            sum2 <- sum2 + count[,f]^2/len[f]^2
            sum3 <- sum3 + count[,f]/len[f]^2
      }
      mean1 <- sum1/ncol(count)
      mean2 <- sum2/ncol(count)
      mean3 <- sum3/ncol(count)
      alls <- sqrt(log(pmax(1,(mean2-mean3)/mean1^2)))
      allm <- log(mean1)-alls^2/2
      allm <- round(allm,digits=4)
      alls <- round(alls,digits=4)
      s <- alls[id]
      m <- allm[id]
      excs <- alls[excid]
      excm <- allm[excid]
      
      norm <- norm[id,]
      rm('alls','allm','count','sum1','sum2','sum3','orisum','mean1','mean2','mean3')
      scalematrix <- function(data) {
            cm <- rowMeans(data)
            csd <- apply(data, 1, sd)
            (data - cm) / csd
      }
      
      data <- scalematrix(norm)
      data <- round(data,digits=4)
      
      library(ClusterR)
      oriclu <- kmeans(data,5000)$cluster
      
      splitclu <- lapply(1:max(oriclu),function(cid) {
            distmat <- dist(norm[oriclu==cid,])
            hclust(distmat)
      })
      
      clulist <- list()
      clulist[['5000']] <- oriclu
      
      options(scipen=999)
      spclu <- split(1:length(oriclu),oriclu)
      names(spclu) <- NULL
      
      for (split in c(2, 4, 8, 16, 32, 64)) {
            cluster <- oriclu
            curclu <- 1
            for (i in 1:max(oriclu)) {
                  if (sum(oriclu==i) == 1) {
                        cluster[oriclu==i] <- curclu
                        curclu <- curclu + 1
                  } else {
                        tmpclunum <- min(sum(oriclu==i),split)
                        tmpcluster <- cutree(splitclu[[i]],tmpclunum)
                        for (j in 1:tmpclunum) {
                              cluster[spclu[[i]][tmpcluster==j]] <- curclu
                              curclu <- curclu + 1
                        }
                  }
            }
            clulist[[as.character(max(cluster))]] <- cluster
      }
      
      options(scipen=999)
      clucenter <- matrix(0,nrow=max(oriclu),ncol=ncol(norm))
      for (i in 1:max(oriclu)) {
            clucenter[i,] <- colMeans(norm[which(oriclu==i),,drop=F])
      }
      
      hclu <- hclust(dist(scalematrix(clucenter)))
      
      for (clufac in c(2, 4, 8, 16, 32)) {
            clunum <- round(max(oriclu) / clufac)
            clu <- cutree(hclu,k = clunum)
            cluster <- rep(0,length(oriclu))
            for (i in 1:max(oriclu)) {
                  cluster[which(oriclu==i)] <- clu[i]
            }
            clulist[[as.character(clunum)]] <- cluster
      }
      
      cluster <- do.call(cbind,clulist)
      allclunum <- as.numeric(colnames(cluster))
      data <- list(ms=list(m=m,s=s),id=id,gr=gr,cluster=cluster,excid=excid,excms=list(m=excm,s=excs),allclunum=allclunum)
      saveRDS(data,file=savepath)
}

