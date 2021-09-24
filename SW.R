
stepwise <- function(GWAS,Gmat,N,fstub) {
 tsnps <- colnames(Gmat) 
 rownames(GWAS) <- GWAS$rs 
 zout <- data.frame(GWAS[tsnps,c("rs","chr","ps","allele1","allele0","af","beta","se")])
 names(zout) <- c("rsid","chromosome","position","allele1","allele2","maf","beta","se")
 raf <- zout$maf
 zout$maf <- raf*(raf<=0.5) + (1-raf)*(raf>0.5)
 #zout$flip <- 1*(raf>0.5)
 zfile <- paste0(fstub,".z")
 write.table(zout,file=zfile,quote=FALSE,row.names=FALSE,col.names=TRUE)
 message("z file written to ",zfile)
 
ldout <- cor(Gmat,Gmat)
#ldout <- cor(refG,refG) 
diag(ldout) <- 1
ldfile <- paste0(fstub,".ld")
write.table(ldout,file=ldfile,quote=FALSE,row.names=FALSE,col.names=FALSE)

 nsnps <- ncol(Gmat)
 
 
 # output files from FINEMAP
 snpfile <- paste0(fstub,".snp")
 confile <- paste0(fstub,".config")
 logfile <- paste0(fstub,".log")
 crfile <- paste0(fstub,".cred") # /net/bliss-04-nfs/scratch/jennifer/Uganda/FM.cred
 
 mfile <- paste0(fstub,".master")
 write.table("z;ld;snp;config;cred;log;n_samples",file=mfile,quote=FALSE,col.names=FALSE,row.names=FALSE,append=FALSE)
write.table(paste(zfile,ldfile,snpfile,confile,crfile,logfile,N,sep=";"),file=mfile,quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE)
 
 FMpath <- "./finemap_v1.4_x86_64"
 #finemap_v1.3.1_x86_64/finemap_v1.3.1_x86_64"
 #finemap_v1.1_x86_64/finemap_v1.1_x86_64"
 
 fline <- paste0(FMpath," --cond --log --corr-config 0.8 --cond-pvalue 0.000001  --in-files ",mfile)
 system(fline)

fmresults <- read.table(confile,header=TRUE,as.is=TRUE,sep=" ")
msnps <- fmresults$config
if(length(msnps) > 1 & msnps[1]=="null")  fmresults <- fmresults[-which(msnps=="null"),]
pv <- fmresults$pvalue 
msnps <- fmresults$config
snpPP <- data.frame(snp=msnps, pvalue=pv,stringsAsFactors=FALSE)

return(snpPP)
}


 gmodK.fn <- function(modsnps,snpgroups) {
 tmp <- as.character(modsnps)
 msnps <- unlist(strsplit(tmp,"%"))

 gsnps <- groupIDs.fn(snpgroups,msnps)
 out <- paste(gsnps[order(gsnps)],collapse="%")

 return(out)

 }

ingroup.fn <- function(snp,group) {
# used by groupIDs.fn
# psnp <- unlist(strsplit(snp,"[.]"))[2]
# pos <- unlist(strsplit(psnp,"_"))[1]
# ind <- grep(pos,group)
 ind <- grep(snp,group)
 1*(length(ind)>0)
 }

groupIDs.fn <- function(snpgroups,Msnps) {
# outputs snp group for each snp in Msnps; if not in a group, output rsid
 ng <- length(snpgroups)
 ns <- length(Msnps)
 names(snpgroups) <- LETTERS[1:ng]
 G <- character(ns)
 for(i in 1:ns) {
  if(Msnps[i]=="1") {
    G[i] <- "null"
  } else {
  check <- unlist(lapply(snpgroups,ingroup.fn,snp=Msnps[i]))
  if(sum(check)>0) { G[i] <- names(which(unlist(check)>0))
  } else {
          	G[i] <- Msnps[i]
                }
   }
  }
 return(G)
}

 

PPmodGroups <- function(bmod,snpgroups,minPP=0.05) {
# bmod <- best.models(SM)
 snpmods <- bmod$snps
 pp <- bmod$PP
 Gmods <- apply(as.matrix(snpmods),1,gmodK.fn,snpgroups)
 modpp <- data.frame(model=Gmods,PP=pp)
 mods <- unique(Gmods)
 GPP <- apply(as.matrix(mods),1,function(x) sum(modpp$PP[modpp$model==x]))
 out <- data.frame(model=mods,GPP=GPP)
 out <- out[order(out$GPP,decreasing=TRUE),]
 Gout <- out[out$GPP>minPP,]
 Sout <- data.frame(SNPmod=snpmods,Gmod=Gmods,PP=pp)
 return(list(snp=Sout,group=Gout))
 }




