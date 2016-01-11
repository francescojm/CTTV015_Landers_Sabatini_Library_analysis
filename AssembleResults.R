
library(stringr)

Letter_dirs<-dir('../../RESULTS/CTTV015_Landers_Sabatini_Library_analysis/outComes/')
nLetters<-length(Letter_dirs)

# includedGuides<-read.table('../../DATA/CTTV015_Landers_Sabatini_Library_analysis/internal/Copy of aac7041_SM_Table_S1-toFrancesco2.txt',sep='\t',header = TRUE,stringsAsFactors = FALSE)
# includedGuides<-includedGuides[which(includedGuides$Present.in.the.current.library.=='Yes' | includedGuides$BbsI.sites.=='Yes'),1]
# save(includedGuides,file='../../DATA/CTTV015_Landers_Sabatini_Library_analysis/alreadyincludedGuides.rdata')

for (i in 1:nLetters){
    currentFC<-dir(paste('../../RESULTS/CTTV015_Landers_Sabatini_Library_analysis/outComes/',Letter_dirs[i],sep=''))

    nfiles<-length(currentFC)
    
    for (j in 1:nfiles){
            print(c(Letter_dirs[i],currentFC[j]))
            currentFN<-paste('../../RESULTS/CTTV015_Landers_Sabatini_Library_analysis/outComes/',Letter_dirs[i],'/',currentFC[j],sep='')
            fileContent<-read.table(currentFN,sep='\t',header = TRUE,stringsAsFactors = FALSE)
            
            flag<-2
            while(flag>0){
             currentList<-fileContent[flag,1]   
             currentList<-unlist(str_split(currentList,","))
             if(sum(is.element(currentList,includedGuides))==0){flag<-0}
             else{flag<-flag+1}
            }
            
            
            if (length(currentList)==5){
                currentGene<-unlist(str_split(currentFC[j],'.txt'))[1]
                currentBunch<-cbind(rep(currentGene,5),currentList)
            }else{
                    currentBunch<-NULL
                }
    
    }
    
    if (i == 1){
        TOTRES<-currentBunch
    }else{
        TOTRES<-rbind(TOTRES,currentBunch)
    }
    }