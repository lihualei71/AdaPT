#####################################################################################
## This script implements AdaPT and Methods on GEOquery Data:
## The code is from Ang Li and Rina Foygel Barber (2016)
##   "Multiple testing with the structure adaptive Benjamini-Hochberg algorithm"
## We modify the source code from 
##   "http://www.stat.uchicago.edu/~rina/sabha/gene_drug_data_example.R"
## The modification includes 
##   1) Parameter in Adaptive Seqstep: Lei & Fithian (2016) sets s = alpha instead 
##      fixing it to be 0.05;
##   2) Adding IHW (Ignatiadis et al. 2016) and AdaPT;
##   3) The function to plot Figure 1.
#####################################################################################

setwd("your path")

#####################################################################################
## download & process the gene/drug data
#####################################################################################

gene_drug_get_data = function(data_downloaded = FALSE){
  if(data_downloaded){Data = read.table('gene_drug_data.txt')}else{
    if(!require(GEOquery)){
      install.packages("GEOquery")
    }
    library("GEOquery")
    gds2324 = getGEO('GDS2324')
    eset = GDS2eSet(gds2324, do.log2=TRUE)
    Data = exprs(eset)
    write.table(Data,'gene_drug_data.txt')
  }
  Data
}

gene_drug_get_pvals = function(Data){
  
  # code for processing this data set & for accumulation tests
  # the code in this function is from the paper:
  # Ang Li and Rina Foygel Barber, "Accumulation tests for FDR control in ordered hypothesis testing" (2015). Available from http://arxiv.org/abs/1505.07352
  
  # the following code is copied (with some edits) from:
  # http://www.stat.uchicago.edu/~rina/accumulationtests/gene_dosage_experiment.R
  
  ##### Two-sample t-tests on a data matrix
  # This function inputs a data matrix X with n columns. The columns indexed by g_1 and g_2 belong to sample 1 and sample 2, respectively (e.g. cases and controls). For each row of X, the function runs a two-sample t-test for that row.
  ttest_mat=function(X,g1,g2){
    # g1 & g2 give columns for groups 1 & 2
    n1=length(g1);n2=length(g2)
    means1=rowMeans(X[,g1]);means2=rowMeans(X[,g2])
    vars1=rowSums((X[,g1]-means1%*%t(rep(1,n1)))^2)/(n1-1);vars2=rowSums((X[,g2]-means2%*%t(rep(1,n2)))^2)/(n2-1)
    sds_diff = sqrt(vars1/n1 + vars2/n2)
    tstats = (means1-means2)/sds_diff
    dfs = (vars1/n1+vars2/n2)^2 / ((vars1/n1)^2/(n1-1) + (vars2/n2)^2/(n2-1))
    pvals=2*(1-pt(abs(tstats),dfs))
    output=list()
    output$tstats=tstats
    output$dfs= dfs
    output$pvals= pvals
    output
  }
  
  # The next function is the same, except that instead of performing two-sided t-tests, we perform a one-sided t-test for each row of X. The signs s_i specify, for the i-th row of X, whether we are testing for a positive effect or a negative effect.
  signed_ttest_mat=function(X,g1,g2,s){
    n1=length(g1);n2=length(g2)
    means1=rowMeans(X[,g1]);means2=rowMeans(X[,g2])
    vars1=rowSums((X[,g1]-means1%*%t(rep(1,n1)))^2)/(n1-1);vars2=rowSums((X[,g2]-means2%*%t(rep(1,n2)))^2)/(n2-1)
    sds_diff = sqrt(vars1/n1 + vars2/n2)
    tstats = s*(means1 - means2)/sds_diff
    dfs = (vars1/n1+vars2/n2)^2 / ((vars1/n1)^2/(n1-1) + (vars2/n2)^2/(n2-1))
    pvals=(1-pt(tstats,dfs))
    output=list()
    output$tstats=tstats
    output$dfs= dfs
    output$pvals= pvals
    output
  }
  
  # Each row of "Data" is a gene (total: n=22283 genes). 
  # The 25 columns of "Data" correspond to 5 trials each at 5 different dosages: columns 1-5 are zero dose (control group), followed by columns 6-10 at the lowest dose, etc. The entries of "Data" are log-transformed gene expression levels.
  n=dim(Data)[1]
  highdose=21:25;lowdose=6:10;control=1:5
  
  ### Computing p-values
  # Next, we will use the highest-dosage data to produce an ordering on the low-dosage vs control-group p-values, resulting in an ordered sequence of p-values which we will use for the accumulation tests.
  # First, for each gene, produce a pval for highest dose compared to mean of lowest dose + control. Record also the sign of the effect (increased or reduced gene expression at the highest dose, compared to pooled low-dose and control-group data).
  ttest_highdose=ttest_mat(Data,highdose,c(lowdose, control))
  pvals_highdose=ttest_highdose$pvals
  test_signs=sign(ttest_highdose$tstats)
  
  pvals_highdose_small=ttest_mat(Data,highdose[1:2],c(lowdose,control))$pvals
  set.seed(1);pvals_random=1:n
  
  # Next, for each gene we will perform a one-sided t-test (using the sign of the estimated high-dose effect), and then use a permutation test to get a final p-value for this gene. These p-values will then be reordered according to the high-dose results above.
  signed_pvals_lowdose_ttest=signed_ttest_mat(Data,lowdose,control,test_signs)$pvals
  
  signed_pvals_lowdose_ttest_permuted=matrix(0,choose(10,5),n)
  nchoosek=combn(1:10,5)
  for(i in 1:choose(10,5)){
    permord=c(nchoosek[,i],setdiff(1:10,nchoosek[,i]))
    signed_pvals_lowdose_ttest_permuted[i,]=signed_ttest_mat(Data[,permord],lowdose,control,test_signs)$pvals
  }
  signed_pvals_permutation_test=colSums(abs(signed_pvals_lowdose_ttest_permuted)-rep(1,choose(10,5))%*%t(abs(signed_pvals_lowdose_ttest))<=0)/choose(10,5)
  
  output = list()
  output$pvals = signed_pvals_permutation_test*(1-1/(1+choose(10,5))) # multiplying to avoid pvalues exactly equal to 1 due to discretization
  output$ord = order(abs(pvals_highdose))
  output$high_pvals = abs(pvals_highdose)
  output$ord_small = order(abs(pvals_highdose_small))
  output$mod_pvals = abs(pvals_highdose_small)
  output$ord_random = order(abs(pvals_random))
  output
  
}


#####################################################################################
## set up plotting functions for gene/drug data example
#####################################################################################
gene_drug_plot_results <- function(NumRej,filename){
  ### Plot results for each method
  cols <- c('black', 'black', 'black', 'red', 'blue', 'blue','orange', 'orange', 'red')
  ltys <- c(1,5,3,4,1,5,1,5,1)
  pchs <- c(20,15,17,19,20,15,20,15,1)
  methods <- c('SeqStep', 'Accum. Test', 'ForwardStop', 'Ada. SeqStep', 'BH', 'Storey', 'SABHA', 'IHW', 'AdaPT')
  
  pdf(filename,width=12.5,height=4)
  par(mfrow=c(1,3),oma=c(0,0,3,0))
  titles <- c('Highly informative ordering', 'Moderately informative ordering', 'Original ordering')
  for(j in 3:1){
    plot(0:1,0:1,type='n',xlim=range(alphalist),ylim=c(0,4200),xlab=expression(paste('Target FDR level ',alpha)),ylab='# of discoveries', main = titles[j], axes=FALSE, font.main=1, cex.main = 2, cex.lab = 2)
    axis(side=1,at=c(0,0.1,0.2,0.3), cex.axis = 2)
    axis(side=2, cex.axis = 2)
    alpha_pt=c(10,20,30)
    for(i in length(methods):1){
      points(alphalist,NumRej[i,,j],type='l',col=cols[i],lty=ltys[i])
      points(alphalist[alpha_pt],NumRej[i,alpha_pt,j],col=cols[i],pch=pchs[i])
    }
    if(j==3){
      legend(-0.005,4000,methods,col=cols,lty=ltys,pch=pchs,seg.len=3,cex=1.5)
    }
  }   
  mtext("Gene/drug response data: results",outer=TRUE,cex=2,font=2)
  dev.off()
}

#####################################################################################
## run gene/drug data example
#####################################################################################

alphalist <- seq(0.01,0.3,by=0.01) # target FDR level
tau <- 0.5; eps <- 0.1 # parameters for SABHA
thr <- 0.5 # parameter for Storey-BH
thr1 <- 0.1; thr2 <- 0.5 # parameters for adaptive SeqStep

Data = gene_drug_get_data()
# if data already downloaded & saved:
# Data <- gene_drug_get_data(data_downloaded = TRUE)
output <- gene_drug_get_pvals(Data)

### set up all methods
## Accumulation Tests (Li & Barber, 2016a)
source('https://www.stat.uchicago.edu/~rina/accumulationtests/accumulation_test_functions.R')

## Auxiliary functions of SABHA (Li & Barber, 2016b)
source('http://www.stat.uchicago.edu/~rina/sabha/All_q_est_functions.R')

## AdaPT 
source('AdaPT.R')

## Independent Hypothesis Weighting (Ignatiadis et al. 2016)
source("https://bioconductor.org/biocLite.R")
biocLite("IHW")
library("IHW")


## BH procedure (Benjamini & Hochberg, 1995)
BH_method <- function(pvals,alpha){
  khat <- max(c(0,which(sort(pvals)<=alpha*(1:length(pvals))/length(pvals))))
  which(pvals<=alpha*khat/length(pvals))
}

## Storey's BH Procedure (Storey et al. 2004)
Storey_method <- function(pvals,alpha,thr){
  est_proportion_nulls <- min(1,mean(pvals>thr)/(1-thr))
  pvals[pvals>thr] <-  Inf
  khat <- max(c(0,which(sort(pvals)<=alpha/est_proportion_nulls*(1:length(pvals))/length(pvals))))
  which(pvals<=alpha/est_proportion_nulls*khat/length(pvals))
}

## SABHA (Li & Barber, 2016b)
SABHA_method <- function(pvals, qhat, alpha, tau){
  # Use the original, or estimated q as input
  pvals[pvals>tau] <- Inf
  khat <- max(c(0,which(sort(qhat*pvals)<=alpha*(1:length(pvals))/length(pvals))))
  which(qhat*pvals<=alpha*khat/length(pvals))
}

## Adaptive Seqstep (Lei & Fithian, 2016)
Adaptive_SeqStep_method <- function(pvals, alpha, thr1, thr2){ # Lei & Fithian 2016's method
  # thr1 & thr2 correspond to s & lambda (with s<=lambda) in their paper
  fdphat <- thr1 / (1-thr2) * (1 + cumsum(pvals>thr2)) / pmax(1, cumsum(pvals<=thr1))
  if(any(fdphat<=alpha)){
    khat <- max(which(fdphat<=alpha))
    return(which(pvals[1:khat]<=thr1))
  }else{
    return(NULL)
  }
  
}

num_alpha <- length(alphalist)
max_alpha <- max(alphalist)
# gather results
NumRej <- array(0,c(9,num_alpha,3))
res.AdaPT <- list()
n <- length(pvals)
set.seed(1)
# methods: 1 SeqStep, 2 HingeExp, 3 ForwardStop, 4 Adaptive SeqStep, 5 BH, 6 Storey-BH, 7 SABHA, 8 IHW, 9 AdaPT
for(j in 1:3){
  if(j==1){pvals <- output$pvals[output$ord]} 
  else if (j==2){pvals <- output$pvals[output$ord_small]}
  else if (j==3){pvals <- output$pvals[output$ord_random]}
  qhat <- Solve_q_step(pvals,tau,eps)
  for(i in 1:num_alpha){
    NumRej[1,i,j] <- SeqStep(pvals,alpha=alphalist[i],C=2)
    NumRej[2,i,j] <- HingeExp(pvals*(1-1/(1+choose(10,5))),alpha=alphalist[i])
    NumRej[3,i,j] <- ForwardStop(pvals*(1-1/(1+choose(10,5))),alpha=alphalist[i])
    NumRej[4,i,j] <- length(Adaptive_SeqStep_method(pvals, alphalist[i], alphalist[i], tau))
    NumRej[5,i,j] <- length(BH_method(pvals, alphalist[i]))
    NumRej[6,i,j] <- length(Storey_method(pvals, alphalist[i], thr))
    NumRej[7,i,j] <- length(SABHA_method(pvals, qhat, alphalist[i], tau))
    NumRej[8,i,j] <- rejections(ihw(pvals, 1:length(pvals), alphalist[i]))
  }
  if (j < 3) {
      s0 <- c(rep(0.45, floor(n / 2)), rep(0.05, n - floor(n / 2)))
  } else {
      s0 <- rep(0.45, n)
  }
  res.AdaPT[[j]] <- AdaPT(1:length(pvals), pvals, s0 = s0, q.list = alphalist)
  NumRej[9,,j] <- res.AdaPT[[j]]$num_rej
}

## Plot Figure 2.
gene_drug_plot_results(NumRej,'~/study/Courses/adaptive-fdr/adapt_proc/figs/gene_drug_results.pdf')

## Function to plot s(x) and signal strength (Figure 3-5)
illustrate <- function(s, pi, mu, pmax = 0.2, length.out = 3, title){
  par(mar=c(4.1,4.1,2,0.15), mfrow = c(2, 1))
  what.type <- ifelse(pvals < s, 1, ifelse(pvals > 1-s, 2, 3))
  plot((1:n*pmax)/n, type="n", pch=".",xaxs="i",yaxs="i",ylab="p-value",xlab="",
       col=c("red","blue","black")[what.type],yaxt="n",main=title)
  axis(2,at=seq(0, pmax, length.out = length.out), labels = seq(0, pmax, length.out = length.out))
  hhi <- hlo <- .5
  vlo <- -.01; vhi=1.01
  polygon(x=c(1:n,n:1),y=c(s,rep(vlo,n)),col="#FFDDDD",border="red")
  polygon(x=c(1:n,n:1),y=c(1-s,rep(vhi,n)),col="light blue",border="blue")
  points(pvals,pch=".",col=c("red","blue","black")[what.type])
  box()
  signal <- pi / (1 + mu) + (1 - pi) / 2
  plot(signal,ylim=c(0,0.5), type = 'l', lwd = 2,xaxs="i",yaxs="i",ylab="estimated signal strength",xlab="x (index)")
}

## Figure 3 (Left)
pdf("original_02.pdf")
illustrate(res.AdaPT[[3]]$s[,20], res.AdaPT[[3]]$pi[,20], res.AdaPT[[3]]$mu[,20], 0.1, title = "alpha = 0.2")
dev.off()

## Figure 3 (Right)
pdf("original_03.pdf")
illustrate(res.AdaPT[[3]]$s[,30], res.AdaPT[[3]]$pi[,30], res.AdaPT[[3]]$mu[,30], 0.1, title = "alpha = 0.3")
dev.off()

## Figure 4 (Left)
pdf("high_005.pdf")
illustrate(res.AdaPT[[1]]$s[,5], res.AdaPT[[1]]$pi[,5], res.AdaPT[[1]]$mu[,5], 0.3, title = "alpha = 0.05")
dev.off()

## Figure 4 (Right)
pdf("high_01.pdf")
illustrate(res.AdaPT[[1]]$s[,10], res.AdaPT[[1]]$pi[,10], res.AdaPT[[1]]$mu[,10], 0.3, title = "alpha = 0.1")
dev.off()

## Figure 5 (Left)
pdf("moderate_005.pdf")
illustrate(res.AdaPT[[2]]$s[,5], res.AdaPT[[2]]$pi[,5], res.AdaPT[[2]]$mu[,5], 0.1, title = "alpha = 0.05")
dev.off()

## Figure 6 (Left)
pdf("moderate_01.pdf")
illustrate(res.AdaPT[[2]]$s[,10], res.AdaPT[[2]]$pi[,10], res.AdaPT[[2]]$mu[,10], 0.1, title = "alpha = 0.1")
dev.off()


#### Appendix C
## Sensitivity to the choice of x.
x.list <- list(matrix(0, n, 3), matrix(0, n, 3))
x.list[[1]][, 1] <- 1:n
x.list[[1]][, 2] <- output$high_pvals[output$ord]
x.list[[1]][, 3] <- log(x.list[[1]][, 2] / (1 - x.list[[1]][, 2]))
x.list[[2]][, 1] <- 1:n
x.list[[2]][, 2] <- output$mod_pvals[output$ord_small]
x.list[[2]][, 3] <- log(x.list[[2]][, 2] / (1 - x.list[[2]][, 2]))
NumRej.x <- array(0,c(3,num_alpha,2))
n <- length(pvals)
s0 <- c(rep(0.45, floor(n / 2)), rep(0.05, n - floor(n / 2)))
set.seed(1)
for(j in 1:2){
  if(j==1){
    pvals <- output$pvals[output$ord]
  } else if (j==2){
    pvals <- output$pvals[output$ord_small]
  }
  for (i in 1:3){
    x <- x.list[[j]][, i]
    res.AdaPT <- AdaPT(x, pvals, s0 = s0, q.list = alphalist)
    NumRej.x[i,,j] <- res.AdaPT$num_rej
  }
}

cols=1:3
ltys=1:3
pchs=1:3
methods=c("x=order", "x=pvals", "x=logit(pvals)")

pdf('~/study/Courses/adaptive-fdr/adapt_proc/figs/gene_drug_results_x.pdf',width=8.5,height=4)
par(mfrow=c(1,2),oma=c(0,0,3,0))
titles=c('Highly informative ordering', 'Moderately informative ordering')
for(j in 2:1){
  plot(0:1,0:1,type='n',xlim=range(alphalist),ylim=c(0,4200),xlab=expression(paste('Target FDR level ',alpha)),ylab='# of discoveries', main = titles[j], axes=FALSE, font.main=1, cex.main = 1.2, cex.lab = 1.2)
  axis(side=1,at=c(0,0.1,0.2,0.3), cex.axis = 1.2)
  axis(side=2, cex.axis = 1.3)
  alpha_pt=c(10,20,30)
  for(i in length(methods):1){
    points(alphalist,NumRej.x[i,,j],type='l',col=cols[i],lty=ltys[i])
    points(alphalist[alpha_pt],NumRej.x[i,alpha_pt,j],col=cols[i],pch=pchs[i])
  }
  if(j==2){
    legend(-0.005,4000,methods,col=cols,lty=ltys,pch=pchs,seg.len=3,cex=1)
  }
}   
mtext("Gene/drug response data: results (with different x)",outer=TRUE,cex=1.3,font=2)
dev.off()


## Sensitivity to the number of knots
df.list <- c(5, 10, 20, 30)
NumRej.knots <- array(0,c(length(df.list),num_alpha,3))
n <- length(pvals)
set.seed(1)
for(j in 1:3){
  if(j==1){
    pvals <- output$pvals[output$ord]
  } else if (j==2){
    pvals <- output$pvals[output$ord_small]
  } else if (j==3){
    pvals <- output$pvals[output$ord_random]
  }
  if (j < 3) {
    s0 <- c(rep(0.45, floor(n / 2)), rep(0.05, n - floor(n / 2)))
  } else {
    s0 <- rep(0.45, n)
  }
  for (i in 1:length(df.list)){
      res.AdaPT <- AdaPT(1:length(pvals), pvals, s0 = s0, q.list = alphalist, df = df.list[i])
      NumRej.knots[i,,j] <- res.AdaPT$num_rej
  }
}

cols=1:length(df.list)
ltys=1:length(df.list)
pchs=1:length(df.list)
methods=paste0("df=", df.list)

pdf('~/study/Courses/adaptive-fdr/adapt_proc/figs/gene_drug_results_knots.pdf',width=12.5,height=4)
par(mfrow=c(1,3),oma=c(0,0,3,0))
titles=c('Highly informative ordering', 'Moderately informative ordering', 'Original ordering')
for(j in 3:1){
  plot(0:1,0:1,type='n',xlim=range(alphalist),ylim=c(0,4200),xlab=expression(paste('Target FDR level ',alpha)),ylab='# of discoveries', main = titles[j], axes=FALSE, font.main=1, cex.main = 2, cex.lab = 2)
  axis(side=1,at=c(0,0.1,0.2,0.3), cex.axis = 2)
  axis(side=2, cex.axis = 2)
  alpha_pt=c(10,20,30)
  for(i in length(methods):1){
    points(alphalist,NumRej.knots[i,,j],type='l',col=cols[i],lty=ltys[i])
    points(alphalist[alpha_pt],NumRej.knots[i,alpha_pt,j],col=cols[i],pch=pchs[i])
  }
  if(j==3){
    legend(-0.005,4000,methods,col=cols,lty=ltys,pch=pchs,seg.len=3,cex=1.5)
  }
}   
mtext("Gene/drug response data: results (with different number of knots)",outer=TRUE,cex=2,font=2)
dev.off()


## Sensitivity to delta
delta.list <- c(0.025, 0.05, 0.1, 0.2)
NumRej.delta <- array(0,c(length(delta.list),num_alpha,3))
NumRej.step <- array(NA,c(length(delta.list),num_alpha,3))
steps <- matrix(0, length(delta.list), 3)
n <- length(pvals)
set.seed(1)
for(j in 1:3){
  if(j==1){
    pvals <- output$pvals[output$ord]
  } else if (j==2){
    pvals <- output$pvals[output$ord_small]
  } else if (j==3){
    pvals <- output$pvals[output$ord_random]
  }
  if (j < 3) {
    s0 <- c(rep(0.45, floor(n / 2)), rep(0.05, n - floor(n / 2)))
  } else {
    s0 <- rep(0.45, n)
  }
  for (i in 1:length(delta.list)){
    res.AdaPT <- AdaPT(1:length(pvals), pvals, s0 = s0, delta.high = delta.list[i], delta.low = delta.list[i],
                       q.list = alphalist, df = 10)
    NumRej.delta[i,,j] <- res.AdaPT$num_rej
    NumRej.step[i,,j] <- res.AdaPT$step
  }
}

cols=1:length(delta.list)
ltys=1:length(delta.list)
pchs=1:length(delta.list)
methods=paste0("delta=", delta.list)

pdf('~/study/Courses/adaptive-fdr/adapt_proc/figs/gene_drug_results_delta.pdf',width=12.5,height=8.5)
par(mfrow=c(2,3),oma=c(0,0,3,0))
titles=c('Highly informative ordering', 'Moderately informative ordering', 'Original ordering')
for(j in 3:1){
  plot(0:1,0:1,type='n',xlim=range(alphalist),ylim=c(0,4200),xlab=expression(paste('Target FDR level ',alpha)),ylab='# of discoveries', main = titles[j], axes=FALSE, font.main=1, cex.main = 2, cex.lab = 2)
  axis(side=1,at=c(0,0.1,0.2,0.3), cex.axis = 2)
  axis(side=2, cex.axis = 1.5)
  alpha_pt=c(10,20,30)
  for(i in length(methods):1){
    points(alphalist,NumRej.delta[i,,j],type='l',col=cols[i],lty=ltys[i])
    points(alphalist[alpha_pt],NumRej.delta[i,alpha_pt,j],col=cols[i],pch=pchs[i])
  }
  if(j==3){
    legend(-0.005,4000,methods,col=cols,lty=ltys,pch=pchs,seg.len=3,cex=1.5)
  }
} 
for(j in 3:1){
  plot(0:1,0:1,type='n',xlim=range(alphalist),ylim=c(0,1000),xlab=expression(paste('Target FDR level ',alpha)),ylab='# of Iterations', main = titles[j], axes=FALSE, font.main=1, cex.main = 2, cex.lab = 2)
  axis(side=1,at=c(0,0.1,0.2,0.3), cex.axis = 2)
  axis(side=2, cex.axis = 2)
  alpha_pt=c(10,20,30)
  for(i in length(methods):1){
    points(alphalist,NumRej.step[i,,j],type='l',col=cols[i],lty=ltys[i])
    points(alphalist[alpha_pt],NumRej.step[i,alpha_pt,j],col=cols[i],pch=pchs[i])
  }
  if(j==3){
    legend(-0.005,4000,methods,col=cols,lty=ltys,pch=pchs,seg.len=3,cex=1.5)
  }
}   
mtext("Gene/drug response data: results (with different delta)",outer=TRUE,cex=2,font=2)
dev.off()

