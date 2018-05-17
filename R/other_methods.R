### set up all methods
## Accumulation Tests (Li & Barber, 2016a)
source('https://www.stat.uchicago.edu/~rina/accumulationtests/accumulation_test_functions.R')

## Auxiliary functions of SABHA (Li & Barber, 2016b)
source('http://www.stat.uchicago.edu/~rina/sabha/All_q_est_functions.R')

## Independent Hypothesis Weighting (Ignatiadis et al. 2016)
tmp <- try(library("IHW"))
if (class(tmp) == "try-error"){
    source("https://bioconductor.org/biocLite.R")
    biocLite("IHW")
    library("IHW")
}


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

## Barber-Candes procedure (Barber and Candes, 2015)
BC_method <- function(pvals,alpha){
    sorted.mask.pvals <- sort(pmin(pvals,1-pvals))
    fdphat <- sapply(sorted.mask.pvals, function(thresh){
        (1+sum(pvals>=1-thresh))/max(1,sum(pvals<=thresh))
    })
    khat <- which(fdphat<=alpha)
    if (length(khat) == 0){
        return(khat)
    } else {
        khat <- max(khat)
        phat <- sorted.mask.pvals[khat]
        return(which(pvals<=phat))
    }
}

## SABHA (Li & Barber, 2016b)
source('All_q_est_functions.R')
Solve_q_ordered_simple <- function(pvals, tau, eps, ADMM_params){
    target.num <- 5000
    n <- length(pvals)
    if (n <= target.num){
        qhat <- Solve_q_ordered(pvals, tau, eps, ADMM_params) 
        return(qhat)
    }
    num.reps <- ceiling(n / target.num)
    new.pvals <- sapply(1:target.num, function(i){
        ind <- pmin((i - 1) * num.reps + (1:num.reps), n)
        pvals[ind[1]]
    })
    qhat <- Solve_q_ordered(new.pvals, tau, eps, ADMM_params)
    qhat <- rep(qhat, each = num.reps)[1:n]
    return(qhat)
}

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

## Independent Filtering (with oracle threshold)
tmp <- try(library("genefilter"))
if (class(tmp) == "try-error"){
    source("https://bioconductor.org/biocLite.R")
    biocLite("genefilter")
    library("genefilter")
}

IF.oracle <- function(pvals, x, alpha){
    theta.list <- seq(0, 1, 0.05)
    R.list <- sapply(theta.list, function(theta){
        R1 <- filtered_R(filter = x, test = pvals,
                         theta = theta, method = "BH",
                         alpha = alpha)
        R2 <- filtered_R(filter = rev(x), test = rev(pvals),
                         theta = theta, method = "BH",
                         alpha = alpha)
    })
    return(max(R.list))
}

## IHW (with oracle nbins)
ihw.oracle <- function(pvals, x, alpha){
    nrejs <- sapply(1:15, function(nbins){
        rejections(ihw(pvals, x, alpha, nbins = nbins))
    })
    max(nrejs)
}

