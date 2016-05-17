library(vegan)
library(ape)
require(phangorn)
library(MonoPhy)
library(phytools)

#### identification of SNPS in genes that are reated to envrionmental conditions.
## in order to distinct between vertical inherited SNPS one has to identify groups of ECOtypes that evolve from same common ancestor
## hypothesis here is: after divergence from A. lyrata, ecotyoes evolving from same ancenstor should show similar pattern in 
## the amount of accumulated non synonymous substitutions (dN) compared to the others. Analysis of ther amount of dN per gene followed
## by cluster analysis should group ecotypes that have a common ancestor together


dnds_mat <- read.csv("~/Desktop/HHU_Projects/git_projects/Arabidopsis/input_files/dn_mat.txt", header=FALSE)

country <- read.table("~/Desktop/HHU_Projects/git_projects/Arabidopsis/input_files/country.txt", quote="\"", comment.char="")

fun <- read.table("~/Desktop/HHU_Projects/git_projects/Arabidopsis/input_files/functions.txt", quote="\"", comment.char="")

colnames(dnds_mat)=country$V1

## remove all rows with rowsum = 0 

xi=which(rowSums(dnds_mat)>0)
dnds_mat=dnds_mat[xi,]
fun=fun$V1[xi]

# category abundance
tt=table(fun)
# dnds mean per main cat
avg_fun=aggregate(dnds_mat,by=list(fun),FUN="mean",simplify = T)

# add names to rows
rownames(avg_fun)=names(tt)

avg_fun=avg_fun[,-1]

# normalize dnds table by max dnds per gene
xm=apply(dnds_mat, 1, max)
xx=dnds_mat/xm

#distance between ecotypes
vd=vegdist(t(xx),method = "bray")
# hierarchical clustering
#hc=hclust(vd)
#plot(hc)
# phylogenetic tree - neighbor joining
t=nj(vd)
#root the tree
rt=midpoint(t)

# extract clusters with minmal dN distance between the eco types.
# use k-mean cluster algorithm
dd=as.matrix(vd)

## in order to which number of clusters is represents best
kk_sm=rep(0,79)
kk_mn=rep(0,79)
for (i in 1:79){
km=kmeans(dd,i)
kk_sm[i]=sum(km$withinss)
kk_mn[i]=mean(km$withinss)
}



