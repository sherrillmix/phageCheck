library(dnar)
library(parallel)
library(IRanges)

dataDir<-"/media/THING1/alexandra/4Gut_Phage_and_Bacteria_DNA_Modifications/014Blast16s/"
faFiles<-list.files(dataDir,'fasta$')
blastFiles<-file.path('work',basename(sub('\\.fasta','\\.blast.gz',faFiles)))
names(faFiles)<-names(blastFiles)<-sub('.fasta','',sub('_100bp_.*','',faFiles))
seqs<-mclapply(file.path(dataDir,faFiles),read.fa,mc.cores=10)
blasts<-mclapply(blastFiles,read.table,header=FALSE,mc.cores=10,col.names=c('qName','tName','percIdent','alignLength','mismatch','gapOpens','qStart','qEnd','tStart','tEnd','eVal','bitScore'),stringsAsFactors=FALSE)
names(seqs)<-names(blasts)<-names(faFiles)
results<-data.frame('name'=unlist(lapply(seqs,function(x)x$name)),'nchar'=unlist(lapply(seqs,function(x)nchar(x$seq))),'set'=rep(names(seqs),sapply(seqs,nrow)),stringsAsFactors=FALSE)
blastCover<-lapply(blasts,function(xx){
  out<-by(xx[,c('qStart','qEnd')],xx$qName,function(yy){
    width(reduce(IRanges(yy$qStart,yy$qEnd)))
  })
  return(unlist(out))
})

results$blastCover<-NA
for(ii in names(blastCover)){
  message(ii)
  results[results$set==ii,'blastCover']<-blastCover[[ii]][results[results$set==ii,'name']]
}
results[is.na(results$blastCover),'blastCover']<-0
