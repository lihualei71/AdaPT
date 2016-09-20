##############################################################################
## These functions implement AdaPT from the paper:
## Lihua Lei & William Fithian,
##   "AdaPT: An interactive procedure for multiple testing with side information"
## Available from http://arxiv.org/abs/
##############################################################################

## Function to calculate FDPhat (1 + At) / (Rt v 1).
## Inputs:
##    pvals: p-values
##    s: s(x)
fdp.hat <- function(pvals, s) {
    A <- sum(pvals > 1 - s)
    R <- sum(pvals <= s)
    return((1 + A) / max(R, 1))
}

## Function to calculate the number of rejections (p-values below s(x)).
## Inputs:
##    pvals: p-values
##    s: s(x)
num.rej <- function(pvals, s) {
    sum(pvals <= s)
}


## Function to initialize Algorithm 2 (discussed in Appendix A.2)
## Inputs: 
##   x: covariates
##   p: p-values
##   pmin: s(x)
##   pmax: 1 - s(x)
##   df: the degree of freedom used for natural spline. 
init.em <- function(x, p, pmin, pmax, 
                    df = 10){
    n <- length(p)
    z <- ifelse(p < pmin | p > pmax, 1, (pmin + 1 - pmax) / (pmin - pmax))
    mod <- lm(z ~ ns(x, df = df))
    pix <- pmax(predict(mod, type = "response"), 0)
    imputed.p <- sapply(1:n, function(i){
        if (p[i] < pmin[i] | p[i] > pmax[i]){
            temp <- runif(1) < (1 + pix[i]) / 2
            return(ifelse(temp, min(p[i], 1 - p[i]), max(p[i], 1 - p[i])))
        } else {
            return(p[i])
        }        
    })
    imputed.p[imputed.p == 1] <- 0.9999
    mod2 <- glm(-log(imputed.p) ~ ns(x, df = df), family = Gamma())
    mux <- predict(mod2, type = "response")
    return(list(mux = mux, pix = pix))
}

## Function to update pi and mu in Algorithm 2.
## Inputs: 
##   x: covariates
##   p: p-values
##   pmin: s(x)
##   pmax: 1 - s(x)
##   pix0: pi(x) at step 0. Initialized as in Appendix A.2 if not set.
##   mux0: mu(x) at step 0. Initialized as in Appendix A.2 if not set.
##   df: the degree of freedom used for natural spline. 
glm.mixem.censor <- function(x, p, pmin, pmax,
                             pix0 = NULL, mux0 = NULL,
                             df = 10,
                             num.steps = 5){
    n <- length(p)
    if (is.null(pix0) || is.null(mux0)){
        init <- init.em(x, p, pmin, pmax, df)
        if (is.null(pix0)){
            pix0 <- pmin(init$pix, 1)
        }
        if (is.null(mux0)){
            mux0 <- pmax(init$mux, 1)
        }
    } 
    pix <- pix0
    mux <- mux0
    converge <- FALSE
    for (i in 1: num.steps){
        z <- ifelse(p < pmin | p > pmax, 
                    1 / (1 + 2 * (1 - pix) / pix * mux / (p^(1 / mux - 1) + (1 - p)^(1 / mux - 1))),
                    1 / (1 + (1 - pix) / pix * mux * p^(1 - 1 / mux)))
        z <- pmax(z, 0.00001)
        p.new <- ifelse(p < pmin | p > pmax,
                        exp((log(p) * p^(1 / mux - 1) + log(1 - p) * (1 - p)^(1 / mux - 1)) / (p^(1 / mux - 1) + (1 - p)^(1 / mux - 1))),
                        p)
        p.new[p.new == 0] <- 0.00001
        mod <- suppressWarnings(glm(z ~ ns(x, df = df), family = binomial()))
        pix.new <- predict(mod, type = "response")
        mod2 <- glm(-log(p.new) ~ ns(x, df = df), family = Gamma(), weights = z)
        mux.new <- predict(mod2, type = "response")
        if (max(abs(pix.new - pix)) < 0.01) {
            converge <- TRUE
            break
        }
        pix <- pmax(pmin(pix.new, 1), 0)
        mux <- pmax(mux.new, 1)
    }
    pix <- pmax(pmin(pix.new, 1), 0)
    mux <- pmax(mux.new, 1)
    return(list(mux = mux, pix = pix, converge = converge))
}

## Function to update s(x) (discussed in section 4.2)
## Inputs: 
##   x: covariates
##   pvals: p-values
##   s: s(x)
##   delta: tolerance parameter (discussed in section 4.2)
##   pix0: pi(x) at step 0. Initialized as in Appendix A.2 if not set.
##   mux0: mu(x) at step 0. Initialized as in Appendix A.2 if not set.
##   ...: other args passed to glm.mixem.censor
find.snew.mix <- function(x, pvals, s, delta,
                          pix0 = NULL, mux0 = NULL, ...){
    n <- length(s)
    temp <- glm.mixem.censor(x, pvals, s, 1-s, pix0, mux0, ...)
    pix <- temp$pix
    mux <- temp$mux
    lower.c <- 0
    upper.c <- 1
    level.curve.frac <- function(c, mux, s, pix){
        if (c == 0){return(0)}
        pix[pix == 0] <- 0.00001
        s.new <- (1/c + (1-pix)/pix*mux*(1-c)/c)^(mux / (1 - mux))
        s.new[mux <= 1] <- 0
        s.new <- pmin(s.new, s)
        sum(s.new[s.new > 0] * (1 - pix[s.new > 0])) /
            sum(s[s.new > 0] * (1 - pix[s.new > 0]))
    }
    lower.frac <- level.curve.frac(lower.c, mux, s, pix)
    upper.frac <- level.curve.frac(upper.c, mux, s, pix)
    while (lower.frac < 1 - delta - 0.000001 & upper.frac > 1 - delta + 0.000001) {
        mid.c <- (lower.c + upper.c) / 2
        mid.frac <- level.curve.frac(mid.c, mux, s, pix)
        if (mid.frac < 1 - delta) {
            lower.c <- mid.c
            lower.frac <- mid.frac
        } else {
            upper.c <- mid.c
            upper.frac <- mid.frac
        }
    }
    c <- ifelse(lower.frac > 1 - delta - 0.000001, lower.c, upper.c)
    s.new <- pmin((1/c + (1-pix)/pix*mux*(1-c)/c)^(mux / (1 - mux)), s)
    s.new[pix <= 0.05] <- 0
    s.new <- ifelse(s.new < (1 - delta) * s, (1 - delta) * s, s.new)    
    return(list(s.new = s.new, pix = pix, mux = mux))
}

## Main Function (AdaPT)
## Inputs: 
##   x: covariates
##   pvals: p-values
##   s0: s(x) at step 0
##   delta.high: tolerance parameter when fdphat > q + .05 (discussed in section 4.2 and 4.3)
##   delta.low: tolerance parameter when fdphat <= q + .05 (discussed in section 4.2 and 4.3)
##   B: buffer size (discussed in section 4.2)
##   q.list: a list of FDR levels (discussed in section 4.3)
##   quiet: display the information of intermediate steps, including 
##          FDPhat, number of rejections, the plot of s(x), \hat{\pi}(x) 
##          and \hat{\mu}(x).
##   ...: other args passed to find.snew.fun
library('splines')
AdaPT <- function(x, pvals, s0 = rep(0.45, length(pvals)), 
                  delta.high = 0.05, delta.low = 0.05, 
                  B = 5,
                  q.list = seq(0.05, 0.3, 0.01), 
                  quiet = FALSE,
                  find.snew.fun = find.snew.mix,
                  ...){
    n <- length(pvals)
    m <- length(q.list)
    s <- s0
    fdp <- fdp.hat(pvals, s)    
    qind <- m
    q <- q.list[qind]
    step <- 0
    s.return <- matrix(0, n, m)
    pi.return <- matrix(0, n, m)
    mu.return <- matrix(0, n, m)
    fdp.return <- rep(0, m)
    step.return <- rep(0, m)
    num_rej.return <- rep(0, m)
    Rej.buffer <- rep(n, B) 
    par(mfrow = c(3, 1))
    while (fdp > q | qind > 0){
        step <- step + 1
        if (step == 1){
            pix <- NULL
            mux <- NULL
        }
        delta <- ifelse(fdp > q + 0.05, delta.high, delta.low)
        temp <- find.snew.fun(x, pvals, s, delta, pix, mux, ...)
        s.new <- temp$s.new
        pix <- temp$pix
        mux <- temp$mux
        if (!quiet){
            plot(x, s, type = 'l')
            lines(x, s.new, col = "red")
            plot(x, pix, type = 'l')
            plot(x, 1 / (1 + mux), type = 'l')
            abline(h = 1, col = "blue")
        }
        s <- s.new
        fdp <- fdp.hat(pvals, s)
        R <- sum(pvals <= s)
        Rej.buffer[2:B] <- Rej.buffer[1:(B-1)]
        Rej.buffer[1] <- R
        if (R < 10) {
            break
        }
        if (fdp <= q) {
            s.return[, qind] <- s
            pi.return[, qind] <- pix
            mu.return[, qind] <- mux
            num_rej.return[qind] <- R
            fdp.return[qind] <- fdp
            step.return[qind] <- step
            if (qind == 1) break
            qind <- max(which(q.list < fdp))
            q <- q.list[qind]
        }
        if (!quiet){
            print(paste0("Step ", step, ": FDP ", fdp, ", Number of Rej. ", R))
        }
        if (all(Rej.buffer == Rej.buffer[1])){
            s <- s * (1 - delta.low)
        }
    }
    return(list(num_rej = cummax(num_rej.return), fdp = fdp.return, s = s.return,
                pi = pi.return, mu = mu.return, step = step.return))
}
