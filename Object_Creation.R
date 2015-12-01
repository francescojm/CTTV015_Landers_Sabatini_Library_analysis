## Object_Creation.R

source('Preamble.R')



cat('\n- Creating guide annotation object...')
LS_guide_annotations<-read.table(paste(IDD,'GUIDE_ANNOTATIONS.txt',sep='/'),sep='\t',header=TRUE,row.names=1,
                                 stringsAsFactors = FALSE)
cat('\n+ Done!\n')

cat('\n- Creating raw guide counts object...')
LS_raw_guide_counts<-read.table(paste(IDD,'GUIDE_COUNTS.txt',sep='/'),sep='\t',header=TRUE,row.names=1,
                                stringsAsFactors = FALSE)
cat('\n+ Done!\n')

cat('\n- Creating CRISPR scores object...')
LS_crispr_scores<-read.table(paste(IDD,'CRISPR_SCORES.txt',sep='/'),sep='\t',header=TRUE,row.names=1,
                             stringsAsFactors = FALSE)
cat('\n+ Done!\n')

cat('\n- Creating average across replicates guide counts object...')
LS_averaged_guide_counts<-averageReplicates_and_finalize(LS_raw_guide_counts)
cat('\n+ Done!\n')

cat('\n- Creating geneWise guide sets object...\n')
LS_geneWise_guides<-buildGeneWise_guideSets(LS_guide_annotations)
cat('\n+ Done!\n')

cat('\n- Saving...')
save(LS_guide_annotations,file=paste(RDD,'LS_guide_annotations.rdata',sep='/'))
save(LS_raw_guide_counts,file=paste(RDD,'LS_raw_guide_counts.rdata',sep='/'))
save(LS_crispr_scores,file=paste(RDD,'LS_crispr_scores.rdata',sep='/'))
save(LS_averaged_guide_counts,file=paste(RDD,'LS_averaged_guide_counts.rdata',sep='/'))
save(LS_geneWise_guides,file=paste(RDD,'LS_geneWise_guide.rdata',sep='/'))
cat('\n+ Done!')







