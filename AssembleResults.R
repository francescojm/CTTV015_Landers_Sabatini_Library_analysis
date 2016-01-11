Letter_dirs<-dir('../../RESULTS/CTTV015_Landers_Sabatini_Library_analysis/outComes/')

nLetters<-length(Letter_dirs)

for (i in 1:nLetters){
    
    
    currentFC<-dir(paste('../../RESULTS/CTTV015_Landers_Sabatini_Library_analysis/outComes/',Letter_dirs[i],sep=''))

    nfiles<-length(currentFC)
    
    for (j in 1:nfiles){
            print(c(Letter_dirs[i],currentFC[j]))
            currentFN<-paste('../../RESULTS/CTTV015_Landers_Sabatini_Library_analysis/outComes/',Letter_dirs[i],'/',currentFC[j],sep='')
            
            
            fileContent<-read.table(currentFN,sep='\t',header = TRUE)
        }
    }