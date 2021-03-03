# read.delim(), read.table

ex_tab<-read.delim(file='extdata/example_tab.txt',header=T,sep='\t')

head(ex_tab)

ex_tab2<-read.table(file='extdata/example_tab.txt',header = T)

head(read.delim('extdata/example_tab.txt'))

# read.csv()

ex_csv <- read.csv(file='extdata/example.csv')
ex_csv

# read.xlsx()

read.xlsx()

ex_xlsx2 <- read.xlsx(file='extdata/example.xlsx',sheetIndex = 1)
head(ex_xlsx)

# write.table()
head(ex_csv)
class(ex_csv)
write.table(x=ex_csv,file='output/ex.txt')

write.table(x=ex_csv,file='output/ex2.txt',sep='\t',quote = F,row.names=F)

write.csv(x=ex_csv,file='output/ex2.txt',quote = F,row.names=F)

write.xlsx(x=ex_csv,file='output/ex2.xlsx',row.names=F)

ex_tab2<-read.table(file='extdata/example_tab.txt',header = F,skip=10)

ex_tab2

# save()

save(ex_csv,file='data/ex.rda')
load(file='data/ex.rda')

data(ex)