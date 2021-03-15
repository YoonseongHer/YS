library('Biobase')

getwd()

# read cufflinks outputs
cuff.mut1 <-read.delim('result/cufflinks/EPHB6_MUT1/genes.fpkm_tracking')
cuff.mut2 <- read.delim('result/cufflinks/EPHB6_MUT2/genes.fpkm_tracking')
cuff.wt1 <- read.delim('result/cufflinks/EPHB6_WT1//genes.fpkm_tracking')
cuff.wt2 <- read.delim('result/cufflinks/EPHB6_WT2//genes.fpkm_tracking')

#check data type
class(cuff.mut1)

#check header
head(cuff.mut1)

#check order geneID

cuff.mut1$tracking_id==cuff.mut2$tracking_id

length(cuff.mut2$tracking_id)
# make unique key gene ID

kid <- unique(c(cuff.mut1$tracking_id
,cuff.mut2$tracking_id
,cuff.wt1$tracking_id
,cuff.wt2$tracking_id))

length(kid)
class(kid)
##################

# 1. Make empty matrix

assay <- matrix(NA, nrow=length(kid), ncol=4)

rownames(assay) <- kid
colnames(assay) <- c('EPHB6_MUT1','EPHB6_MUT2','EPHB6_WT1','EPHB6_WT2')
head(assay) # initially all values = NA

# 2. Insert FPKM values to matrix for each sample
# find matched index with kid, and ordered 
# FPKM into 'assay' matrix
assay[,1] <- cuff.mut1$FPKM[match(kid, cuff.mut1$tracking_id)]
head(assay)
assay[,2] <- cuff.mut2$FPKM[match(kid, cuff.mut2$tracking_id)]
head(assay)
assay[,3] <- cuff.wt1$FPKM[match(kid, cuff.wt1$tracking_id)]
head(assay)
assay[,4] <- cuff.wt2$FPKM[match(kid, cuff.wt2$tracking_id)]
head(assay)

head(assay)
sum(is.na(assay))
####### Create PhenoData (sample infomation)
sample.id <- colnames(assay)

group <- rep(c('MUT','WT'),each=2)
group

phenoData <- data.frame(sample.id, group,row.names = sample.id)
phenoData

rownames(phenoData) <- sample.id
phenoData

####### FratureData(Gene infomation)
# For R project, data() function can load 'RData' in 'data' folder
data("gencode.v22.df")
# same with load('data/gencode.v22.df.rda)
head(gencode.v22.df)

# match the order with kid

matched.index<-match(kid, gencode.v22.df$gene_id)

featureData <- gencode.v22.df[matched.index,]
head(featureData)

rownames(featureData) <- kid

############ Create ExpresssionSet

##### check

colnames(assay) == rownames(phenoData)
sum(!(rownames(assay) == rownames(featureData)))

# ExpresstionSet(assayData = assay, phnoData=phenoData,featureData=featureData) # error

eset <- ExpressionSet(assayData = assay, 
              phnoData=AnnotatedDataFrame(phenoData),
              featureData=AnnotatedDataFrame(featureData))

save(eset, file='data/eset.rda')

library(Biostrings)

library(rtracklayer)
gtf = import('data/test.gtf.gz')

head(gtf)
class(gtf)
class(gencode.v22.df)
head(gencode.v22.df)
colnames(featureData)
colnames(mcols(gtf))

match(kid,gtf$gene_id)

lat = letters[1:10]

tes = c('b','d','e','z')

match(tes,lat)
gtf[1:10,]
