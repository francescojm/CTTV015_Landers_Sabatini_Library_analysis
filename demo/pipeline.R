plasmid<-load.file('data/HCT116_d14_C1.txt')
library1<-load.file('data/HCT116_d14_T1.txt')
library2<-load.file('data/HCT116_d14_T2.txt')

knitr::kable(stats.data(dataset=plasmid, namecolumn = 1, fullmatchcolumn = 3,
                        extractpattern=expression("^(.+?)_.+"), readcount.unmapped.total = 1786217,
                        type="dataset")[1:10,1:5])


carpools.read.distribution(plasmid, fullmatchcolumn=3, breaks=200,
                           title='plasmid', xlab="log2 Readcount", ylab="# sgRNAs",statistics=TRUE) 
carpools.read.distribution(library1, fullmatchcolumn=3, breaks=200,
                           title='library1', xlab="log2 Readcount", ylab="# sgRNAs",statistics=TRUE) 
carpools.read.distribution(library2, fullmatchcolumn=3, breaks=200,
                           title='library2', xlab="log2 Readcount", ylab="# sgRNAs",statistics=TRUE) 


carpools.read.distribution(library2, fullmatchcolumn=3, breaks=200,
                           title='library2', xlab="log2 Readcount", ylab="# sgRNAs",statistics=TRUE,plotgene = 'CASP8') 

carpools.read.count.vs(dataset=list(library1,plasmid),
                       pairs=FALSE, namecolumn=1, fullmatchcolumn=3, title="", pch=16,
                       normalize=FALSE, norm.function="median",  labelcolor="blue",
                       center=TRUE, aggregated=TRUE)

data.wilcox = stat.DESeq(untreated.list = list(plasmid),agg.function = 'mean',
                          treated.list = list(library2), namecolumn=1, fullmatchcolumn=3,
                         extractpattern=expression("^(.+?)(_.+)"),
                         sorting=FALSE, filename.deseq = "ANALYSIS-DESeq2-sgRNA.tab",
                         fitType="parametric")