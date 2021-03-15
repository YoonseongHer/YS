library(GenomicRanges)
library(Biobase)
gr <- GRanges(seqnames = c('chr1','chr2','chr3')
              ,ranges = IRanges(start = c(1,1,1)
                                ,end = c(10,20,30))
              ,strand = c('-','+','+'))

gr

gr$score = 1:3
gr$name = letters[1:3]
gr

seqnames(gr)
start(gr)
end(gr)
strand(gr)

length(gr)
width(gr)

class(mcols(gr))
gr$score

data(eset)

# featureData to GRanges
head(fData(eset))

fdata.gr<-GRanges(seqnames = fData(eset)$seqnames
        , ranges = IRanges(start = fData(eset)$start
                           ,end = fData(eset)$end)
        ,strand = fData(eset)$strand)

fdata.gr$gene_name <- fData(eset)$gene_name

#chr17:7582285-7773654
#chr13 :47501983-49283632

gr = GRanges(seqnames = c("chr17",'chr13')
        ,ranges = IRanges(start=c(7582285,47501983)
                          ,end = c(7773654,49283632))
        )

overlaps <- findOverlaps(query = gr, subject = fdata.gr)

queryHits(overlaps)
subjectHits(overlaps)

fdata.gr$gene_name[subjectHits(overlaps)][queryHits(overlaps)==1]
