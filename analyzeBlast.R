library(dnar)
library(parallel)
library(IRanges)

condenseBlast<-function(xx,subtract=NULL){
  flip<-xx$qStart>xx$qEnd
  tmp<-xx$qEnd
  xx$qEnd[flip]<-xx$qStart[flip]
  xx$qStart[flip]<-tmp[flip]
  if(!is.null(subtract)){
    subRanges<-out<-by(subtract[,c('qStart','qEnd')],subtract$qName,function(yy){
      return(reduce(IRanges(yy$qStart,yy$qEnd)))
    })
  }
  out<-by(xx[,c('qStart','qEnd','qName')],xx$qName,function(yy){
    reduced<-reduce(IRanges(yy$qStart,yy$qEnd))
    if(!is.null(subtract)&&!is.null(subRanges[[yy$qName[1]]])){
      reduced<-setdiff(reduced,subRanges[[yy$qName[1]]])
    }
    return(sum(width(reduced)))
  })
  return(unlist(out))
}

getBestBlastHit<-function(xx){
  best<-ave(xx$bitScore,xx$qName,FUN=max)
  hits<-tapply(xx[xx$bitScore==best,'tName'],xx[xx$bitScore==best,'qName'],FUN=function(x)paste(x,collapse='|'))
  return(hits)
}
getBestBlastRow<-function(xx){
  best<-ave(xx$bitScore,xx$qName,FUN=max)
  hits<-split(xx[xx$bitScore==best,],xx[xx$bitScore==best,'qName'])
  return(hits)
}


dataDir<-"/media/THING1/alexandra/4Gut_Phage_and_Bacteria_DNA_Modifications/014Blast16s/"
faFiles<-list.files(dataDir,'fasta$')
faFiles<-faFiles[!grepl('^gg',faFiles)]
blastFiles<-file.path('work',basename(sub('\\.fasta','.blast.gz',faFiles)))
blastAllFiles<-file.path('work',basename(sub('\\.fasta','.virus.blast.gz',faFiles)))
blastBactFiles<-file.path('work',basename(sub('\\.fasta','.bacteria.blast.gz',faFiles)))
blastVirginFiles<-file.path('work',basename(sub('\\.fasta','.virgin.blast.gz',faFiles)))
names(faFiles)<-names(blastFiles)<-names(blastAllFiles)<-names(blastBactFiles)<-names(blastVirginFiles)<-sub('.fasta','',sub('_100bp_.*','',faFiles))

seqs<-mclapply(file.path(dataDir,faFiles),read.fa,mc.cores=10)
blasts<-mclapply(blastFiles,read.table,header=FALSE,mc.cores=10,col.names=c('qName','tName','percIdent','alignLength','mismatch','gapOpens','qStart','qEnd','tStart','tEnd','eVal','bitScore'),stringsAsFactors=FALSE)
blastAlls<-mclapply(blastAllFiles,read.table,header=FALSE,mc.cores=10,col.names=c('qName','tName','percIdent','alignLength','mismatch','gapOpens','qStart','qEnd','tStart','tEnd','eVal','bitScore'),stringsAsFactors=FALSE)
blastBacts<-mclapply(blastBactFiles,read.table,header=FALSE,mc.cores=10,col.names=c('qName','tName','percIdent','alignLength','mismatch','gapOpens','qStart','qEnd','tStart','tEnd','eVal','bitScore'),stringsAsFactors=FALSE)
blastVirgin<-mclapply(blastVirginFiles,read.table,header=FALSE,mc.cores=10,col.names=c('qName','tName','percIdent','alignLength','mismatch','gapOpens','qStart','qEnd','tStart','tEnd','eVal','bitScore'),stringsAsFactors=FALSE)
names(seqs)<-names(blasts)<-names(blastAlls)<-names(blastBacts)<-names(blastVirgin)<-names(faFiles)
results<-data.frame('name'=unlist(lapply(seqs,function(x)x$name)),'nchar'=unlist(lapply(seqs,function(x)nchar(x$seq))),'set'=rep(names(seqs),sapply(seqs,nrow)),stringsAsFactors=FALSE)
blastCover<-lapply(blasts,condenseBlast)
blastAllCover<-lapply(blastAlls,condenseBlast)
blastBactCover<-lapply(blastBacts,condenseBlast)
blastVirginCover<-lapply(blastVirgin,condenseBlast)
blastBactCoverMinusPhage<-mapply(function(bac,phage1,phage2)condenseBlast(bac,rbind(phage1,phage2)),blastBacts,blasts,blastAlls)
bactHits<-lapply(blastBacts,getBestBlastHit)
bactHitRows<-lapply(blastBacts,getBestBlastRow)
results$blastAllCover<-NA
results$blastCover<-NA
results$bactCover<-NA
results$virginCover<-NA
results$bactMinusPhage<-NA
results$bactHit<-NA
for(ii in names(blastCover)){
  message(ii)
  results[results$set==ii,'blastCover']<-blastCover[[ii]][results[results$set==ii,'name']]
  results[results$set==ii,'blastAllCover']<-blastAllCover[[ii]][results[results$set==ii,'name']]
  results[results$set==ii,'bactCover']<-blastBactCover[[ii]][results[results$set==ii,'name']]
  results[results$set==ii,'virginCover']<-blastVirginCover[[ii]][results[results$set==ii,'name']]
  results[results$set==ii,'bactMinusPhage']<-blastBactCoverMinusPhage[[ii]][results[results$set==ii,'name']]
  results[results$set==ii,'bactHit']<-bactHits[[ii]][results[results$set==ii,'name']]
}
results[is.na(results$blastCover),'blastCover']<-0
results[is.na(results$blastAllCover),'blastAllCover']<-0
results[is.na(results$bactCover),'bactCover']<-0
results[is.na(results$virginCover),'virginCover']<-0
results[is.na(results$bactMinusPhage),'bactMinusPhage']<-0
results$blastProp<-results$blastCover/results$nchar
results$blastAllProp<-results$blastAllCover/results$nchar
results$bactProp<-results$bactCover/results$nchar
results$virginProp<-results$virginCover/results$nchar
results$bactMinusPhageProp<-results$bactMinusPhage/results$nchar
if(any(results$bactCover-results$bactMinusPhage>results$blastCover+results$blastAllCover))stop(simpleError('Something wrong with phage subtraction'))
results$simpleHit<-sapply(lapply(strsplit(sub('^_','',results[,'bactHit']),'\\|'),function(x)sapply(unique(lapply(strsplit(x,'_'),'[',1:2)),paste,collapse=' ')),paste,collapse='|')

blastSums<-tapply(results$blastCover,results$set,sum)
blastAllSums<-tapply(results$blastAllCover,results$set,sum)
blastBactSums<-tapply(results$bactMinusPhage,results$set,sum)
nBases<-tapply(results$nchar,results$set,sum)
notPhage<-nBases-blastSums
notVirus<-nBases-blastAllSums
notBact<-nBases-blastBactSums
print(select<-c('bacteria_fresh','phage25'))
print(fisher.test(matrix(c(notPhage[select],blastSums[select]),nrow=2)))
print(fisher.test(matrix(c(notVirus[select],blastAllSums[select]),nrow=2)))
print(fisher.test(matrix(c(notBact[select],blastBactSums[select]),nrow=2)))
print(select<-c('bacteria_total','phage25'))
print(fisher.test(matrix(c(notPhage[select],blastSums[select]),nrow=2)))
print(fisher.test(matrix(c(notVirus[select],blastAllSums[select]),nrow=2)))
print(select<-c('bacteria_total','phage44'))
print(fisher.test(matrix(c(notPhage[select],blastSums[select]),nrow=2)))
print(fisher.test(matrix(c(notVirus[select],blastAllSums[select]),nrow=2)))
print(fisher.test(matrix(c(notBact[select],blastBactSums[select]),nrow=2)))
data.frame('phage'=blastSums,'notPhage'=notPhage)
print(select<-c('bacteria_fresh','phage44'))
print(fisher.test(matrix(c(notPhage[select],blastSums[select]),nrow=2)))
print(data.frame('phage'=blastSums,'notPhage'=notPhage))
print(fisher.test(matrix(c(notBact[select],blastBactSums[select]),nrow=2)))

cbind('bact'=blastBactSums,'phage'=blastSums,nBases,'propPhage'=round(blastSums/nBases,3),'propBact'=round(blastBactSums/nBases,3))

phage25BactHits<-sapply(lapply(strsplit(sub('^_','',results[results$set=='phage25'&results$bactMinusPhageProp==1,][,'bactHit']),'\\|'),function(x)sapply(unique(lapply(strsplit(x,'_'),'[',1:2)),paste,collapse=' ')),paste,collapse='|')
sort(table(phage25BactHits))
phage25GenusHits<-sapply(lapply(strsplit(sub('^_','',results[results$set=='phage25'&results$bactMinusPhageProp==1,][,'bactHit']),'\\|'),function(x)sapply(unique(lapply(strsplit(x,'_'),'[',1)),paste,collapse=' ')),paste,collapse='|')
sort(table(phage25GenusHits[phage25BactHits!='Faecalibacterium prausnitzii']))

bacteriaBactHits<-sapply(lapply(strsplit(sub('^_','',results[grepl('bacteria',results$set)&results$bactMinusPhageProp==1,][,'bactHit']),'\\|'),function(x)sapply(unique(lapply(strsplit(x,'_'),'[',1:2)),paste,collapse=' ')),paste,collapse='|')
head(t(t(sort(table(bacteriaBactHits),decreasing=TRUE))),10)


#check location of Faecalibacterium prausnitzii hits
fpraus<-do.call(rbind,bactHitRows[['phage25']][results[results$bactMinusPhageProp==1&results$simpleHit=='Faecalibacterium prausnitzii'&results$set=='phage25','name']])
table(fpraus$tName)
pdf('fpraus.pdf')
for(ii in sort(unique(fpraus$tName))){
  thisData<-fpraus[fpraus$tName==ii,]
  thisData<-thisData[order(thisData$tStart),]
  plot(1,1,type='n',xlim=range(c(thisData$tStart,thisData$tEnd)),main=ii,ylim=c(.5,nrow(thisData)+.5),yaxt='n',ylab='Contigs',xlab='Genomic position')
  segments(thisData$tStart,1:nrow(thisData),thisData$tEnd,1:nrow(thisData),lwd=3)
}
dev.off()


means<-do.call(rbind,lapply(
  split(results[,c('blastProp','blastAllProp','bactProp','virginProp')],results$set),
  function(xx)apply(xx,2,mean)
))
baseProps<-do.call(rbind,mapply(
  function(xx,yy)apply(xx,2,sum)/apply(yy,2,sum),
  split(results[,c('blastCover','blastAllCover','bactCover','virginCover')],results$set),
  split(results[,c('nchar','nchar','nchar','nchar')],results$set),
  SIMPLIFY=FALSE
))
colnames(means)<-colnames(baseProps)<-c('Phage','Virus','Bacteria','Virgin')

cols<-rainbow.lab(nrow(means))

basedBar<-function(datMatrix,base=1,col='white',...){
  ylim<-range(datMatrix)
  pos<-matrix(1:((nrow(datMatrix)+1)*ncol(datMatrix)),ncol=ncol(datMatrix))[1:nrow(datMatrix),]
  xlim<-range(pos)+c(-.5,.5)
  plot(1,1,type='n',...,xlim=xlim,ylim=ylim,bty='n')
  rect(as.vector(pos)-.5,base,as.vector(pos)+.5,as.vector(datMatrix),col=col)
  return(pos)
}

pdf('blastSummary.pdf',height=4,width=4)
  par(mar=c(1.2,4,.3,.1))
  #contig mean
  #log relative to Fresh
  info<-basedBar(apply(means,2,function(x)x[-1]/x[1]),xaxt='n',log='y',yaxt='n',ylab='Mean enrichment in contig coverage\n(relative to bacteria_fresh)',col=cols[-1],mgp=c(2.1,1,0))
  logAxis(2,las=1,mgp=c(1,.7,0))
  axis(1,apply(info,2,mean),colnames(means),padj=1,mgp=c(0,-.5,0),lwd=NA)
  legend('topright',fill=cols[-1],rownames(means)[-1],bty='n')
  abline(h=1,lty=2)
  #nonlog
  info<-barplot(means,beside=TRUE,xaxt='n',ylab='Mean proportion of contigs matching',col=cols,las=1)
  axis(1,apply(info,2,mean),colnames(means),padj=1,mgp=c(0,-.5,0),lwd=NA)
  legend('topleft',fill=cols,rownames(means),bty='n')
  #base props
  #log
  info<-basedBar(apply(baseProps,2,function(x)x[-1]/x[1]),xaxt='n',log='y',yaxt='n',ylab='Mean enrichment in bases matching\n(relative to bacteria_fresh)',col=cols[-1],mgp=c(2.1,1,0))
  logAxis(2,las=1,mgp=c(1,.7,0))
  axis(1,apply(info,2,mean),colnames(means),padj=1,mgp=c(0,-.5,0),lwd=NA)
  legend('topright',fill=cols[-1],rownames(means)[-1],bty='n')
  abline(h=1,lty=2)
  #nonlog
  info<-barplot(baseProps,beside=TRUE,xaxt='n',ylab='Proportion of bases matching',col=cols,las=1)
  axis(1,apply(info,2,mean),colnames(means),padj=1,mgp=c(0,-.5,0),lwd=NA)
  legend('topleft',fill=cols,rownames(means),bty='n')
dev.off()
