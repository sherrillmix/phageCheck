library(dnar)
library(parallel)
library(IRanges)

condenseBlast<-function(xx){
  out<-by(xx[,c('qStart','qEnd')],xx$qName,function(yy){
    width(reduce(IRanges(yy$qStart,yy$qEnd)))
  })
  return(unlist(out))
}


dataDir<-"/media/THING1/alexandra/4Gut_Phage_and_Bacteria_DNA_Modifications/014Blast16s/"
faFiles<-list.files(dataDir,'fasta$')
faFiles<-faFiles[!grepl('^gg',faFiles)]
blastFiles<-file.path('work',basename(sub('\\.fasta','.blast.gz',faFiles)))
blastAllFiles<-file.path('work',basename(sub('\\.fasta','.virus.blast.gz',faFiles)))
names(faFiles)<-names(blastFiles)<-sub('.fasta','',sub('_100bp_.*','',faFiles))
seqs<-mclapply(file.path(dataDir,faFiles),read.fa,mc.cores=10)
blasts<-mclapply(blastFiles,read.table,header=FALSE,mc.cores=10,col.names=c('qName','tName','percIdent','alignLength','mismatch','gapOpens','qStart','qEnd','tStart','tEnd','eVal','bitScore'),stringsAsFactors=FALSE)
blastAlls<-mclapply(blastAllFiles,read.table,header=FALSE,mc.cores=10,col.names=c('qName','tName','percIdent','alignLength','mismatch','gapOpens','qStart','qEnd','tStart','tEnd','eVal','bitScore'),stringsAsFactors=FALSE)
names(seqs)<-names(blasts)<-names(blastAlls)<-names(faFiles)
results<-data.frame('name'=unlist(lapply(seqs,function(x)x$name)),'nchar'=unlist(lapply(seqs,function(x)nchar(x$seq))),'set'=rep(names(seqs),sapply(seqs,nrow)),stringsAsFactors=FALSE)
blastCover<-lapply(blasts,condenseBlast)
blastAllCover<-lapply(blastAlls,condenseBlast)
results$blastAllCover<-NA
results$blastCover<-NA
for(ii in names(blastCover)){
  message(ii)
  results[results$set==ii,'blastCover']<-blastCover[[ii]][results[results$set==ii,'name']]
  results[results$set==ii,'blastAllCover']<-blastAllCover[[ii]][results[results$set==ii,'name']]
}
results[is.na(results$blastCover),'blastCover']<-0
results[is.na(results$blastAllCover),'blastAllCover']<-0
results$blastProp<-results$blastCover/results$nchar
results$blastAllProp<-results$blastAllCover/results$nchar
blastSums<-tapply(results$blastCover,results$set,sum)
blastAllSums<-tapply(results$blastAllCover,results$set,sum)
nBases<-tapply(results$nchar,results$set,sum)
notPhage<-nBases-blastSums
notVirus<-nBases-blastAllSums
print(select<-c('bacteria_fresh','phage25'))
print(fisher.test(matrix(c(notPhage[select],blastSums[select]),nrow=2)))
print(fisher.test(matrix(c(notVirus[select],blastAllSums[select]),nrow=2)))
print(select<-c('bacteria_total','phage25'))
print(fisher.test(matrix(c(notPhage[select],blastSums[select]),nrow=2)))
print(select<-c('bacteria_total','phage44'))
print(fisher.test(matrix(c(notPhage[select],blastSums[select]),nrow=2)))
data.frame('phage'=blastSums,'notPhage'=notPhage)
print(select<-c('bacteria_fresh','phage44'))
print(fisher.test(matrix(c(notPhage[select],blastSums[select]),nrow=2)))
print(data.frame('phage'=blastSums,'notPhage'=notPhage))

