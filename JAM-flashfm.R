

library(flashfm); library(data.table); library(R2BGLiMS); library(parallel); library(MASS)

source("./SW.R")

# reg indicator: is a data frame/matrix containing the following information of the genome region to finemap.
reg.indicator <- read.table('./LoR_FEB21.txt',sep = '\t',fill=T, header=T)
# > head(reg.indicator,4)
# n.regions    p.value chr. lower.bound upper.bound Length Tot.snps ntraits number.trait.1
# 1         1 pval=1e-06    1   109811755   109825591  13836      106       2             17
# 2         2 pval=1e-06    1    55517883    55674945 157062     1057       2             17
# 3         3 pval=1e-06    1    17700017    17873294 173277     1705       2             14
# 4         4 pval=1e-06    1   180122985   180192543  69558      606       3              5
# number.trait.2 number.trait.3 number.trait.4 number.trait.5 number.trait.6
# 1             14             NA             NA             NA             NA
# 2             14             NA             NA             NA             NA
# 3             17             NA             NA             NA             NA
# 4              6              3             NA             NA             NA

# reg is the region id and it goes from 1:nrow(reg.indicator).
# set a value or set it as an iteration variable.
reg <- 1 

##### CHR ID #####
k <- reg.indicator$chr.[reg]

##### TRAITS #####

qt <- reg.indicator[reg,c(9:ncol(reg.indicator[reg,]))][!is.na(reg.indicator[reg,c(9:ncol(reg.indicator[reg,]))])]
traits <- read.table("./trait-codes.txt",header=FALSE,as.is=TRUE,fill=TRUE)
ts <- sapply(qt,function(x) traits[which(traits[,1]==x),2])
M <- length(qt)

rel <- read.table("./related_samples.txt",header=FALSE)
relids <- as.numeric(substring(rel[,1],first=2))
y <- read.table("./pooled_ug_gwas_wgs_with_raw_and_tfmd_phe_20170622.txt",as.is=TRUE,sep="\t",header=TRUE)
y <- y[,-(1:42)]
y <- y[-relids,qt]
colnames(y)<-ts

covY <- cov(y,use="pairwise.complete.obs")
diag(covY) <- 1 ## Uganda traits are transformed to N(0,1) and mean and var estimates 
## based on unrelated not as good since know truth
Vy <- diag(covY)
ybar <- rep(0,M)
names(ybar) <- ts
names(Vy) <- ts


##### FOLDERS AND PATHS #####
regname <- paste0("chr",k,"-",reg.indicator$lower.bound[reg],":",reg.indicator$upper.bound[reg],"/")
path_JAM<-'/JAM-flashfm/'
dir.create(paste0(path_JAM,regname))
dir.create(paste0(path_JAM,regname,'tmpJAM'))
fstub <- paste0(path_JAM,regname,'tmpJAM/',"chr",k,"-",reg.indicator$lower.bound[reg],":",reg.indicator$upper.bound[reg])
savetub <- paste0(path_JAM,regname,'tmpJAM/')
DIRout <- paste0(path_JAM,regname)

##### INTERVAL BOUNDARIES #####
lb<-reg.indicator$lower.bound[reg]
ub<-reg.indicator$upper.bound[reg]

##### Loading GWAS data #####
aux<-matrix(c(1:length(ts),qt),nrow=length(ts),ncol=2)
betas<-GWAS<-vector("list",M)
names(betas) <- ts

for (i in 1:M){
  GWAS[[i]] <-read.table(file=paste0("./chr.",k,".qt.",aux[i,2],".assoc.txt"),header=TRUE,as.is=TRUE)
  GWAS[[i]] <- GWAS[[i]][!duplicated(GWAS[[i]]$ps),]
  ind.reg <- which((GWAS[[i]]$ps<=ub)&(GWAS[[i]]$ps>=lb))
  GWAS[[i]] <- GWAS[[i]][ind.reg,]
  GWAS[[i]]$rs <- paste0("chr",k,"_",unlist(lapply(strsplit(GWAS[[i]]$rs, split="\\:"), function(x) x[2])))
  betas[[i]] <-GWAS[[i]][,"beta_uganda"]
  names(betas[[i]]) <-GWAS[[i]]$rs
  rownames(GWAS[[i]]) <- GWAS[[i]]$rs
}

##### EFFECTIVE SAMPLE SIZE #####
#  effective sample size for each trait was computed using the flashfm::Neff function 
# with the vectors of reference (or minor) allele frequencies
# and the standard errors of SNP effect estimates for all SNPs in the trait GWAS.
# The effective sample size results are in Supp Data SD1.18.

Ne <- numeric(M)
for(i in 1:M) Ne[i] <- scan(paste0("./Neff.qt.",qt[i],".txt"))[2]


##### GENOTYPE MATRIX #####
Gmat<-fread(paste0("./chr",k,'-',reg.indicator$lower.bound[reg],':',reg.indicator$upper.bound[reg],
                   '/Gmat_file.txt'), sep=' ', header = FALSE)
Gmat <- Gmat[,-1]
rnames <- paste0("chr",gsub(":","_",Gmat$V2))
aux <- Gmat[,(1:3)]
Gmat <- Gmat[,-(1:3)]
Gmat <- as.matrix(Gmat)
X <- t(Gmat) # snps columns
colnames(X)<-rnames

rs.names<-beta.names<-list()
rs.names[[1]]=beta.names[[1]]<-colnames(X)
for(i in 2:(M+1)){
  rs.names[[i]]<-GWAS[[i-1]]$rs
  beta.names[[i]]<-names(betas[[i-1]])
}

ind.rs.names<-Reduce(intersect, rs.names) 
ind.beta.names<-Reduce(intersect, beta.names) 
X<-X[,ind.rs.names]

dim.GWAS<-length.betas<-c()
# check intersection between GWAS and GMAT
for (i in 1:M){
  GWAS[[i]]<-GWAS[[i]][ind.rs.names,]
  betas[[i]]<-betas[[i]][ind.beta.names]
  dim.GWAS[i]<-dim(GWAS[[i]])[1]
  length.betas[i]<-length(betas[[i]])
}

##### SNP INFORMATION MATRIX #####
snpinfo <- GWAS[[1]][,c("rs","ps","allele1","allele0")]
rownames(snpinfo) <- snpinfo$rs
snpinfo <- snpinfo[colnames(X),]
colnames(snpinfo)<-c('rs','BP','A1', 'A2')
Nqq<-makeNlist.rel(Ne,y)$Nqq

print('pre-processing done')

####### JAM EXPANDED #######
main.input <- JAMexpanded.multi(betas,X,snpinfo,ybar,Vy,N=Ne,chr=k,fstub=fstub,mafthr=0,
                                path2plink="./plink",r2=0.9,
                                save.path=savetub,related=TRUE,y=y)
names(main.input$SM) <- ts
save(main.input, file=paste0(DIRout,'main.input.Rdata'))
print('flashfm.input done')


###### SUMMARY STATS #######
ss.stats <- flashfm::summaryStats(Xmat=TRUE,ybar.all=ybar,main.input=main.input)
save(ss.stats, file=paste0(DIRout,'ss.stats.Rdata'))
print('Summary stats done')


###### FLASH-FM #######
start.time <- Sys.time()
if (reg.indicator[reg,8]==4 | reg.indicator[reg,8]==5){
  fm.multi <- flashfm::flashfm(main.input, TOdds=1,covY=covY,ss.stats=ss.stats,cpp=0.95,maxmod=NULL,fastapprox=FALSE, NCORES=M)
} else {
  fm.multi <- flashfm::flashfm(main.input, TOdds=1,covY=covY,ss.stats=ss.stats,cpp=0.99,maxmod=NULL,fastapprox=FALSE, NCORES=M)
}
end.time <- Sys.time()
time <- end.time - start.time
save(fm.multi, file=paste0(DIRout,'fm.multi.Rdata'))
print('flashfm done')

print(as.numeric(time, units="secs"))
print(as.numeric(time, units="mins"))
print(as.numeric(time, units="hours"))
write.table(as.numeric(time, units="secs"), file=paste0(DIRout,'time.txt'), col.names = F, row.names = F)

###### SNPGROUPS AND PP SUMMARIZE #######
snpGroups <- makeSNPgroups2(main.input,fm.multi,is.snpmat=T,r2.minmerge=0.6)
save(snpGroups, file=paste0(DIRout,'snpgroups.Rdata'))

mpp.pp <- PPsummarise(fm.multi, snpGroups, minPP=0.05) 
save(mpp.pp, file=paste0(DIRout,'mpp.pp.Rdata'))
print('snpGroups and pp summarize')

###### StepWise method ######
sw <- vector("list",M)
bmod <- vector("list",M)
for(i in 1:M) {
  names(GWAS[[i]])[c(9,6,7)]<-c("af","beta","se")
  bmod[[i]]<-stepwise(GWAS[[i]],X,Nqq[i,i],fstub)
  sw[[i]]<-data.frame(snps=bmod[[i]]$snp,groupIDs.fn(snpGroups,Msnps=bmod[[i]]$snp),pvalue=bmod[[i]]$pvalue)
}
save(sw, file=paste0(DIRout,'sw.Rdata'))
print('stepwise done')


###################################################################################################
############################################ THE END! #############################################
###################################################################################################
