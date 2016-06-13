
load('../../DATA/CTTV015_Landers_Sabatini_Library_analysis//R//LS_raw_guide_counts.rdata')
load('../../DATA/CTTV015_Landers_Sabatini_Library_analysis//R//LS_geneWise_guide.rdata')
load('../../DATA/CTTV015_Landers_Sabatini_Library_analysis//R//LS_guide_annotations.rdata')
load('../../DATA/CTTV015_Landers_Sabatini_Library_analysis//R//LS_crispr_scores.rdata')

load('../../DATA/CTTV015_Landers_Sabatini_Library_analysis//R//geneToAdd.rdata')

### function_library.R

## Data Manipulation
averageReplicates_and_finalize<-function(countsObj){
    
    cn<-str_split(colnames(countsObj),'..Rep.')
    duplicated_idx<-which(unlist(lapply(cn,FUN = 'length'))>1)
    
    cn<-unlist(cn)
    
    namesFirstPart<-cn[duplicated_idx+0:(length(duplicated_idx)-1)]
    namesSecondPart<-cn[setdiff(1:length(cn),1:(length(duplicated_idx)*2))]
    
    sampleNames<-c(namesFirstPart,namesSecondPart)
    usampleNames<-unique(sampleNames)
    
    ncellLines<-length(usampleNames)
    
    finalSet<-matrix(NA,nrow(countsObj),ncellLines,dimnames = list(rownames(countsObj),usampleNames))
    for (i in 1:ncellLines){
        idx<-which(sampleNames==usampleNames[i])
        if(length(idx)>1){
            finalSet[,i]<-rowMeans(countsObj[,idx])    
        }else{
            finalSet[,i]<-countsObj[,i]
        }
    }
    
    return(finalSet)
}
buildGeneWise_guideSets<-function(guideAnnotObj){
    ugenes<-unique(guideAnnotObj$Symbol)
    ugenes<-ugenes[!is.na(ugenes)]
    ngenes<-length(ugenes)
    guideIndexes<-matrix(NA,length(ugenes),10,dimnames = list(ugenes,paste('g',1:10,'_idx',sep=''))) 
    
    pb<-txtProgressBar(min = 1,max = ngenes,style = 3)
    
    for (i in 1:ngenes){
        setTxtProgressBar(pb = pb,value = i)
        currentGids<-which(!is.na(match(guideAnnotObj$Symbol,ugenes[i])))
        currentGids<-c(currentGids,rep(NA,10-length(currentGids)))
        guideIndexes[i,]<-currentGids
    }
    
    close(pb)
    
    return(guideIndexes)
}
getGuides<-function(gsymbol,what='names'){
    idxSet<-LS_geneWise_guides[gsymbol,]
    idxSet<-idxSet[!is.na(idxSet)]
    
    
    if (what=='names'){
        GS<-rownames(LS_guide_annotations[idxSet,]) 
    }else{
        if(what!='all'){
            GS<-LS_guide_annotations[idxSet,what]    
        }
    }
    
    GS<-GS[which(rowSums(LS_raw_guide_counts[GS,seq(1,10,2)])>=400)]
    
    return(GS)
    
}

## Statistics
RNAcountsFCs<-function(countsObj){
    
    counts<-countsObj
         counts<-countsObj[which(rowSums(countsObj)>=400),]
    TOTALreads_per_sample<-t(matrix(rep(colSums(counts),nrow(counts)),ncol(counts),nrow(counts)))
         
    counts<-(counts+1)
         
    BpM_FC<-log2((counts[,seq(2,ncol(counts),2)]/TOTALreads_per_sample[,seq(2,ncol(counts),2)])/
                     rowMeans(counts[,seq(1,ncol(counts),2)]/rowMeans(TOTALreads_per_sample[,seq(1,ncol(counts),2)])))
    
    return(BpM_FC)
}

CRISPR_Score<-function(GUIDE_FCs,gene,sampleName,guides=NULL){
    sampleName<-paste(sampleName,'.final',sep='')
    
    GUIDE_FCs<-cbind(rowMeans(GUIDE_FCs[,1:2]),GUIDE_FCs[,3:5])
    colnames(GUIDE_FCs)[1]<-'KBM7.final'
    
    if(length(guides)==0){
        targetingGuides<-getGuides(gene)    
    }else{
        targetingGuides<-guides
    }
    
    
    AA<-colMeans(GUIDE_FCs[targetingGuides,])
    
    nsamples<-length(sampleName)
     PVALS<-rep(NA,nsamples)
     for (i in 1:nsamples){
        KST<-ks.test(GUIDE_FCs[targetingGuides,sampleName[i]],GUIDE_FCs[,sampleName[i]])
        PVALS[i]<-KST$p.value
     }
     names(PVALS)<-sampleName
    
    RES<-cbind(AA[sampleName],PVALS)
    colnames(RES)<-c('CS','p')
    return(RES)
}

GREEDY_optimal_guide_search<-function(GUIDE_FCs,gene){
    
    print(gene)
    
    options(warn=-1)
    sampleName=c('KBM7','K562','Jiyoye','Raji')
    
    originalGuides<-getGuides(gene)
    
    noriginalGuides<-length(originalGuides)
    
    if(noriginalGuides>=5){
        guideCombos<-combn(x = noriginalGuides,5)    
    }else{
        guideCombos<-matrix(1:length(noriginalGuides),length(originalGuides),1)
    }
    
    ncomb<-ncol(guideCombos)

    
    avgDist<-rep(NA,ncomb+1)
    cCS<-matrix(NA,ncomb+1,4)
    cP<-matrix(NA,ncomb+1,4)
    
    
    
    original_CS_pattern<-CRISPR_Score(GUIDE_FCs,gene,sampleName)
    
    avgDist[1]<-0
    cCS[1,]<-t(original_CS_pattern[,1])
    cP[1,]<-t(original_CS_pattern[,2])
    
    original_CS_pattern[,1]<-abs(original_CS_pattern[,1])
    original_CS_pattern[,2]<- -log10(original_CS_pattern[,2])
    
    names(avgDist)[1]<-paste(paste('all',noriginalGuides,'guides'))
    
    
    
    pb<-txtProgressBar(min = 0,max = ncomb,style = 3)  
    
    for (i in 1:ncomb){
  
        setTxtProgressBar(pb,i)
        currentGuideSet<-originalGuides[guideCombos[,i]]
        current_CS_pattern<-CRISPR_Score(GUIDE_FCs,gene,sampleName,guides = currentGuideSet)
        
        cCS[i+1,]<-current_CS_pattern[,1]
        cP[i+1,]<-current_CS_pattern[,2]
        current_CS_pattern[,1]<-abs(current_CS_pattern[,1])
        current_CS_pattern[,2]<- -log10(current_CS_pattern[,2])
        
        distMatrix<-as.matrix(dist(rbind(original_CS_pattern,current_CS_pattern)))
        avgDist[i+1]<-mean(c(distMatrix[5,1],distMatrix[6,2],distMatrix[7,3],distMatrix[8,4]))
        names(avgDist)[i+1]<-paste(currentGuideSet,collapse=',')
    }
    
    close(pb)
    colnames(cCS)<-paste(sampleName,'CS')
    colnames(cP)<-paste(sampleName,'p')
    
    toSave<-cbind(avgDist,cCS,cP)
    toSave<-toSave[order(toSave[,1]),]
    
    toSave<-rbind(colnames(toSave),toSave)
    
    bestSubset<-originalGuides[guideCombos[,order(avgDist)[1]]]
    
    best_CS_pattern<-CRISPR_Score(GUIDE_FCs,gene,sampleName,guides = bestSubset)
    best_CS_pattern[,1]<-abs(best_CS_pattern[,1])
    best_CS_pattern[,2]<- -log10(best_CS_pattern[,2])
    
#     
#     par(mfrow=c(1,2))
    
    
    DIRECTORY<-paste('../../RESULTS/CTTV015_Landers_Sabatini_Library_analysis/outComes/',str_sub(gene,start = 1,end = 1),'/',sep = '')
    
    if(!file.exists(DIRECTORY)){
        dir.create(DIRECTORY)
    }
    
    
    
    write.table(toSave,quote=FALSE,sep='\t',row.names = TRUE,col.names = FALSE,
                file=paste(DIRECTORY,
                           gene,'.txt',sep=''))
#     pdf(paste('../../RESULTS/CTTV015_Landers_Sabatini_Library_analysis/pdfs/',gene,'.pdf',sep=''),
#         width = 15,height = 7)
#     layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths = c(5,2))
#     plot(c(original_CS_pattern[,1],best_CS_pattern[,1]),
#          c(original_CS_pattern[,2],best_CS_pattern[,2]),
#          pch=c(rep(16,4),rep(17,4)),bg='white',col=c('red','blue','green','purple','red','blue','green','purple'),
#          cex=2,xlab='absolute CRISPR score',ylab='-log10 p-value',main=gene)
# 
#     COLS<-c('red','blue','green','purple')
#     for (i in 1:4){
#         lines(c(original_CS_pattern[i,1],best_CS_pattern[i,1]),
#               c(original_CS_pattern[i,2],best_CS_pattern[i,2]),
#               col=COLS[i],lwd=2)            
#     }
#          
#     abline(h=-log10(0.05),lty=2)
#     abline(v=1,lty=2)
#     par(mar=c(0,0,4,0))
#     plot(0,0,xaxt='n',yaxt='n',frame.plot = FALSE,col=NA,xlab='',ylab='')
#     legend('top',pch=c(16,17,15,15,15,15),legend=c(paste('all',noriginalGuides,'guides'),
#                                              'best 5 guide subset',
#                                              'KBM7','K562','Jiyoye','Raji'),
#            col=c('black','black',COLS))
#     dev.off()
    options(warn=0)
}

GENES<-rownames(LS_crispr_scores)

GUIDE_FCs<-RNAcountsFCs(LS_raw_guide_counts)

nGENES<-length(GENES)
for (i in 1:nGENES){
    GREEDY_optimal_guide_search(GUIDE_FCs,GENES[i])
}