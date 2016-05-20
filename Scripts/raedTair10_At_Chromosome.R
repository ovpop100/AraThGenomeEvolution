### read Arabidopsis reference genome (TAIR10)

source("https://bioconductor.org/biocLite.R")

#biocLite("genbankr")
suppressPackageStartupMessages(library(genbankr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(seqinr))

### read genBank files (all 5 chromosomes)

#chr 1
gb1 = readGenBank("/Users/ovip/desktop/HHU_Projects/git_projects/Arabidopsis/input_files/Tair10_At_genome/NC_003070_chr1.gb")
#chr 2
gb2 = readGenBank("/Users/ovip/desktop/HHU_Projects/git_projects/Arabidopsis/input_files/Tair10_At_genome/NC_003071_chr2.gb")
#chr 3
gb3 = readGenBank("/Users/ovip/desktop/HHU_Projects/git_projects/Arabidopsis/input_files/Tair10_At_genome/NC_003074_chr3.gb")
#chr 4
gb4 = readGenBank("/Users/ovip/desktop/HHU_Projects/git_projects/Arabidopsis/input_files/Tair10_At_genome/gbkfile-4.gb")
#chr 5
gb5 = readGenBank("/Users/ovip/desktop/HHU_Projects/git_projects/Arabidopsis/input_files/Tair10_At_genome/NC_003076_chr5.gb")

#### convert data into dataframe 
## gene INFO
df1=GenomicRanges::as.data.frame(genes(gb1))
df2=GenomicRanges::as.data.frame(genes(gb2))
df3=GenomicRanges::as.data.frame(genes(gb3))
df4=GenomicRanges::as.data.frame(genes(gb4))
df5=GenomicRanges::as.data.frame(genes(gb5))

## CDS INFO
cd1=GenomicRanges::as.data.frame(cds(gb1))
cd2=GenomicRanges::as.data.frame(cds(gb2))
cd3=GenomicRanges::as.data.frame(cds(gb3))
cd4=GenomicRanges::as.data.frame(cds(gb4))
cd5=GenomicRanges::as.data.frame(cds(gb5))

## intergenic regions
# extract intergenic regions and convert them to factors (list)

ig1=GenomicRanges::as.factor(intergenic(gb1))
ig2=GenomicRanges::as.factor(intergenic(gb2))
ig3=GenomicRanges::as.factor(intergenic(gb3))
ig4=GenomicRanges::as.factor(intergenic(gb4))
ig5=GenomicRanges::as.factor(intergenic(gb5))

# substitute : with a -
ig1=strsplit(gsub(":",'-',ig1),'-')
ig2=strsplit(gsub(":",'-',ig2),'-')
ig3=strsplit(gsub(":",'-',ig3),'-')
ig4=strsplit(gsub(":",'-',ig4),'-')
ig5=strsplit(gsub(":",'-',ig5),'-')

# create a matrix with start and end for each intergenic region for all chromosomes
igg1=t(matrix(as.numeric(unlist(ig1)),3,length(ig1)))
igg2=t(matrix(as.numeric(unlist(ig2)),3,length(ig2)))
igg3=t(matrix(as.numeric(unlist(ig3)),3,length(ig3)))
igg4=t(matrix(as.numeric(unlist(ig4)),3,length(ig4)))
igg5=t(matrix(as.numeric(unlist(ig5)),3,length(ig5)))

## combine the matrices to one genome matrix with all chromosomes start and end of intergenic region
igg=rbind(igg1,igg2,igg3,igg4,igg5)
colnames(igg)=c('chr','start','end')



#### EXTRACT GENOMIC SEQUENCE for each chromosome
nucSeq1=as.character(as.matrix(getSeq(gb1)))
nucSeq2=as.character(as.matrix(getSeq(gb2)))
nucSeq3=as.character(as.matrix(getSeq(gb3)))
nucSeq4=as.character(as.matrix(getSeq(gb4)))
nucSeq5=as.character(as.matrix(getSeq(gb5)))
## convert DNA strings to integers (a=1,c=2,g=3,t=4) all other characters like (N,M,Y,D etc...) are omitted and
## represented as NaNs. ACTG accounts for 99.5% of the genome. 
numSeq1=s2n(nucSeq1,base4 = F)
numSeq2=s2n(nucSeq2,base4 = F)
numSeq3=s2n(nucSeq3,base4 = F)
numSeq4=s2n(nucSeq4,base4 = F)
numSeq5=s2n(nucSeq5,base4 = F)


## create a big LIST with all infos and save it

refGenome=list(
    README=c("LIST containing information on ORF data: geneInfo, exon data: cdsInfo, genomic Sequence as ACTG XYDMW
    and numeric version of acgt coded as 1,2,3,4. The other characters are NA. Intergenic positions start, stop are 
    available as matrix in intergenicPos"),
  
    chr1=list(
              geneInfo=df1,
              cdsInfo=cd1,
              genomicSeq=nucSeq1,
              numericSeq=numSeq1
    ),
    chr2=list(
              geneInfo=df2,
              cdsInfo=cd2,
              genomicSeq=nucSeq2,
              numericSeq=numSeq2
    ),
    chr3=list(
              geneInfo=df3,
              cdsInfo=cd3,
              genomicSeq=nucSeq3,
              numericSeq=numSeq3
    ),
    chr4=list(
              geneInfo=df4,
              cdsInfo=cd4,
              genomicSeq=nucSeq4,
              numericSeq=numSeq4
    ),
    chr5=list(
              geneInfo=df5,
              cdsInfo=cd5,
              genomicSeq=nucSeq5,
              numericSeq=numSeq5
    ),
    
    intergenicPos=igg
    
)

save(refGenome,file="/Users/ovip/desktop/HHU_Projects/git_projects/Arabidopsis/input_files/refGenome_structureR")



### some modified function in package genbankr:::.seqTypeFromLocus
#new version 
#function(locus) {
#   gsub("^[[:space:]]*LOCUS[[:space:]]+[^[:space:]]+[[:space:]]+[[:digit:]]+[[:space:]]+([^[:space:]]+).*",
#        "\\1",
#        locus)
# }
# 
# ###old version
# function (locus) 
# {
#   gsub("^[[:space:]]*LOCUS[[:space:]]+[[:alnum:]]+[[:space:]]+[[:digit:]]+[[:space:]]+([^[:space:]]+).*", 
#        "\\1", locus)
# }
# 
# test string for the modified function
# k="LOCUS       NC_000932             154478 bp    DNA     circular PLN 26-MAR-2010"
# 
# test gbk file
# gb = readGenBank("/Users/ovip/desktop/HHU_Projects/git_projects/Arabidopsis/input_files/Tair10_At_genome/test.gb")
