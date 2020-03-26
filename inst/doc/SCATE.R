## ------------------------------------------------------------------------
# load in SCATE
options(warn=-1)
suppressMessages(library(SCATE))
set.seed(12345)
# set up locations to bam files
bamlist <- list.files(paste0(system.file(package="SCATE")
,"/extdata/example"),full.names = T)
head(bamlist)

## ------------------------------------------------------------------------
satac <- satacprocess(input=bamlist,type='bam',libsizefilter=1000)
# Number of elements in satac
length(satac)
# Content in the first element as an example
satac[[1]]

## ----eval=FALSE----------------------------------------------------------
#  satac <- satacprocess(input=grlist,type='gr',libsizefilter=1000)

## ------------------------------------------------------------------------
names(satac) <- sub('.*/','',names(satac))

## ------------------------------------------------------------------------
clusterres <- cellcluster(satac,genome="hg19",clunum=2,perplexity=5)

## ------------------------------------------------------------------------
# tSNE results
tsne <- clusterres[[1]]
head(tsne)
# clustering results
cluster <- clusterres[[2]]
cluster
# cell 'GSM1596840.bam' belongs to cluster 1, and cell 'SRR1779746.bam' belongs to cluster 2.
# aggregated signal for CRE cluster
aggsig <- clusterres[[3]]
aggsig[1:3,1:3]

## ------------------------------------------------------------------------
library(ggplot2)
plotdata <- data.frame(tSNE1=tsne[,1],tSNE2=tsne[,2],Cluster=as.factor(cluster))
ggplot(plotdata,aes(x=tSNE1,y=tSNE2,col=Cluster)) + geom_point()

## ------------------------------------------------------------------------
res <- SCATE(satac,genome="hg19",cluster=cluster,clusterid=NULL,ncores=10,verbose=TRUE)
# check the 10000-10005th row of the matrix
res[10000:10005,]

## ------------------------------------------------------------------------
# use similar ways to construct the cluster
usercellcluster <- rep(1:2,each=9)
names(usercellcluster) <- names(satac)
# check the contents of the cluster
usercellcluster

## ----eval=FALSE----------------------------------------------------------
#  userclusterres <- SCATE(satac,genome="hg19",cluster=usercellcluster)

## ------------------------------------------------------------------------
region <- data.frame(chr=c('chr5','chr5'),start=c(50000,50700),end=c(50300,51000))
region

## ------------------------------------------------------------------------
extractres <- extractfeature(res,region,mode='overlap')
extractres

## ----eval=FALSE----------------------------------------------------------
#  extractres <- extractfeature(res,region,mode='overlap',folder='destination folder')

## ------------------------------------------------------------------------
peakres <- peakcall(res)
# check the result for the first cluster
head(peakres[[1]])

## ----eval=FALSE----------------------------------------------------------
#  write.table(peakres[[1]],file='your file.bed',sep='\t',quote=F,col.names = F,row.names = F)

## ----eval=FALSE----------------------------------------------------------
#  piperes <- SCATEpipeline(bamlist,genome="hg19",clunum=2,perplexity=5,ncores=10)
#  # get the cell cluster results, same as calling 'cellcluster' function.
#  cluster <- piperes[['cellcluster']]
#  # get the SCATE outputs, same as calling 'SCATE' function.
#  SCATEres <- piperes[['SCATE']]
#  # get the peak calling results, same as calling 'peakcall' function.
#  peakres <- piperes[['peak']]

## ----eval=FALSE----------------------------------------------------------
#  extractres <- extractfeature(piperes[['SCATE']],region,mode='overlap',folder='destination folder')

## ----eval=FALSE----------------------------------------------------------
#  makedatabase(datapath,savepath,bamfile=bamfile,cre=cre,genome='hg19')

## ----eval=FALSE----------------------------------------------------------
#  makedatabase(datapath=NULL,savepath,bamfile=bamfile,cre=cre,genomerange=genomerange)

## ----eval=FALSE----------------------------------------------------------
#  piperes <- SCATEpipeline(bamlist,datapath='path to new database')

## ----eval=FALSE----------------------------------------------------------
#  clusterres <- cellcluster(satac,datapath='path to new database',clunum=2,perplexity=5)
#  cluster <- clusterres[[2]]
#  res <- SCATE(satac,datapath='path to new database',cluster=cluster,clusterid=NULL)

## ----eval=FALSE----------------------------------------------------------
#  sessionInfo()

