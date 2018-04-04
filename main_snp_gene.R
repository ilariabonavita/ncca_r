########################
# 
# Main script for running ncca for testying association with phenotype and genotype, association 
# is done genewise along the genome.
# 
#########################

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(CCP)
library(DBI)

# base directory with data
base.dir <- args[1]                       # base.dir <- '/Users/ilaria_bonavita/statgen/dyslexia/ags'
# name and extention of the genotype 
snp.file <- args[2]                       # snp.file <- 'ags_snp.common.order.RDS'
# name and extention of the phenotype file 
pheno.file <- args[3]                     # pheno.file <- 'ags_traits.phe.common.order.RDS'
# portion of training examples
train.frac <- as.numeric(args[4])         # train.frac <- 0.60
# nearest neighbours for distance matrix in ncca
NNs <- as.numeric(args[5])                # NNs <- 20
# number of the desired canonical variates 
dim.out <- as.numeric(args[6])            # dim.out <- 2

setwd(base.dir)
source('/Users/ilaria_bonavita/statgen/dyslexia/ncca.rsvd.R')
source('/Users/ilaria_bonavita/statgen/dyslexia/snp4gene_db.R')

# create experiment directory to save results
exp.dir <- paste0('snp_gene_train',train.frac,'_dimout',dim.out,'_NNs',NNs,'onegene')
dir.create(file.path(base.dir,exp.dir))
# Load data
snp <- readRDS(file.path(base.dir,'00_data',snp.file))
pheno <- readRDS(file.path(base.dir,'00_data',pheno.file))


#Create map db
# the db is used to map the snps to the gene they belong to
print('Creating db....')
rds.file <- file.path(base.dir,'00_data','rs.pos.holland.txt')
map <- init.map(rds = rds.file)
mapdb <- setup.db(dbname = 'mapfin2')

# or use a db already exsistent 
#mapdb <- dbConnect(RSQLite::SQLite(),"/Users/ilaria_bonavita/statgen/dyslexia/ags/00_data/map.sqlite")

# Get train and test indeces - randomly selected
set.seed(123)
N <- nrow(snp)
print(N)
N_paired <- floor(train.frac*N)
print(N_paired)
PairedIndices <- sample(1:N, N_paired)
UnpairedIndices <- setdiff(1:N,PairedIndices)


Y <- pheno[PairedIndices,]
YV <- pheno[UnpairedIndices,]

XX <- snp[PairedIndices,]
XXV <- snp[UnpairedIndices,]

### uncomment these two lines if you want to perform PERMUTATION test
# XX <- XX[sample(nrow(XX)),]
# XXV <- XXV[sample(nrow(XXV)),]
######

# get gene list, every genes is mapped to a set of snps
# we select genes that include at least 4 snps include - this number can be changed
gene.list <- dbGetQuery(mapdb, 'SELECT gene FROM sgmap GROUP BY gene HAVING count(gene) > 3 ORDER BY count(gene) DESC')
gene.list <- gene.list$gene[-1]
saveRDS(gene.list,file.path(base.dir,'00_data','gene.list.RDS'))

# crate empty matrices for saving canonical correlations and p-values for train and validation sets
pv <- matrix(0,ncol = length(gene.list), nrow = dim.out)
cc <- matrix(0,ncol = length(gene.list), nrow = dim.out)
pvV <- matrix(0,ncol = length(gene.list), nrow = dim.out)
ccV <- matrix(0,ncol = length(gene.list), nrow = dim.out)
colnames(pv) <- gene.list
colnames(cc) <- gene.list
colnames(pvV) <- gene.list
colnames(ccV) <- gene.list

# perform ncca for every gene-snps vs phenotype
for(j in c(1:length(gene.list))){ 
  gene <- gene.list[j]
  # print(gene)
  # select the snps included in gene j
  snp.gene <- snps4gene(mapdb,gene)
  X <- XX[,na.omit(pmatch(snp.gene$snp,colnames(XX),duplicates.ok = T))]
  XV <- XXV[,na.omit(pmatch(snp.gene$snp,colnames(XX)))]
  
  nout <- ncca.rsvd(X,Y,XV=XV,YV=YV,d=dim.out,nx=NNs,ny=NNs)
  
  pv[,j] <- nout$p_XY$p.value
  pvV[,j] <- nout$p_XVYV$p.value
  cc[,j] <- nout$cor_XY
  ccV[,j] <- nout$cor_XVYV
  
  # uncomment if you want the ncca-transformaed data matrices 
  #X_new <- nout$X_new
  #Y_new <- nout$Y_new
  #XV_new <- nout$XV_new
  #YV_new <- nout$YV_new
  
  #print(j)
  #print('done')
}

# write the ordered list of p-values for the validtion set
pvVord <- data.frame(pval=pvV[1,order(pvV[1,])]) 

# save results 
saveRDS(as.data.frame(pv),file.path(base.dir,exp.dir,'pval.RDS'))
saveRDS(as.data.frame(cc),file.path(base.dir,exp.dir,'cc.RDS'))
saveRDS(as.data.frame(pvV),file.path(base.dir,exp.dir,'pvalV.RDS'))
saveRDS(as.data.frame(ccV),file.path(base.dir,exp.dir,'ccV.RDS'))
saveRDS(pvVord,file.path(base.dir,exp.dir,'pvVord.RDS'))


# read results
#exp.dir <- '/Users/ilaria_bonavita/statgen/dyslexia/finland/snp_gene_train0.6_dimout2_NNs18'
# cc <- readRDS(file.path(exp.dir,'cc.RDS'))
# ccV <-readRDS(file.path(exp.dir,'ccV.RDS'))
# pval <- readRDS(file.path(exp.dir,'pval.RDS'))
# pvalV <- readRDS(file.path(exp.dir,'pvalV.RDS'))
#         
# plot(-log(pvalV[1,],10),type = 'l')
