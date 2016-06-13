library(stringr)
fc<-dir('../../DATA/CTTV015_Landers_Sabatini_Library_analysis//raw//txts')

nf<-length(fc)


membmat<-matrix(0,nf,10)

for (i in 1:nf){
    print(i)
    TMP<-read.table(paste('../../DATA/CTTV015_Landers_Sabatini_Library_analysis//raw//txts/',fc[i],sep=''),sep='\t',stringsAsFactors = FALSE)
    guides<-TMP[3,1]
    
    guides<-unlist(str_split(guides,','))
    guides<-unlist(str_split(guides,'_'))[seq(2,10,2)]

    membmat[i,as.numeric(guides)]<-1
    
}