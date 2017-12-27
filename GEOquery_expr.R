source("other_methods.R")
source("AdaPT.R")

## Load Data: with variable name "output"

if (file.exists("GEOquery.RData")){
    load("GEOquery.RData")
} else {
    source("GEOquery_get_pvals.R")
    load("GEOquery.RData")    
}

#####################################################################################
## Compare methods with highly informative/moderately informative/
## original ordering
#####################################################################################

alphalist <- seq(0.01, 0.3, by=0.01) # target FDR level
tau <- 0.5; eps <- 0.1 # parameters for SABHA
thr <- 0.5 # parameter for Storey-BH
thr1 <- 0.1; thr2 <- 0.5 # parameters for adaptive SeqStep
ADMM_params <- c(10^2, 10^3, 2, 5000, 1e-3) # alpha_ADMM,beta,eta,max_iters,converge_thr
set.seed(1)
n <- length(output$pvals)

num_alpha <- length(alphalist)
max_alpha <- max(alphalist)
# gather results
NumRej <- array(0, c(13, num_alpha, 2))
res.AdaPT <- list()
x <- data.frame(x = 1:n)

# methods: 1 SeqStep, 2 HingeExp, 3 ForwardStop, 4 Adaptive SeqStep, 5 BH, 6 Storey-BH, 7 Barber-Candes, 8 SABHA (step), 9 SABHA (ordered), 10 IHW, 11, IHW (oracle), 12 IF (oracle), 13 AdaPT
for(j in 1:2){
    if (j==1){
        pvals <- output$pvals[output$ord]
    } else if (j==2){
        pvals <- output$pvals[output$ord_small]
    } 

    qhat_step <- Solve_q_step(pvals, tau, eps)
    qhat_ordered <- Solve_q_ordered_simple(pvals, tau, eps, ADMM_params)
    for(i in 1:num_alpha){
        NumRej[1,i,j] <- SeqStep(pvals, alpha = alphalist[i], C = 2)
        NumRej[2,i,j] <- HingeExp(pvals, alpha = alphalist[i])
        NumRej[3,i,j] <- ForwardStop(pvals, alpha = alphalist[i])
        NumRej[4,i,j] <- length(Adaptive_SeqStep_method(pvals, alphalist[i], alphalist[i], tau))
        NumRej[5,i,j] <- length(BH_method(pvals, alphalist[i]))
        NumRej[6,i,j] <- length(Storey_method(pvals, alphalist[i], thr))
        NumRej[7,i,j] <- length(BC_method(pvals, alphalist[i]))
        NumRej[8,i,j] <- length(SABHA_method(pvals, qhat_step, alphalist[i], tau))
        NumRej[9,i,j] <- length(SABHA_method(pvals, qhat_ordered, alphalist[i], tau))
        NumRej[10,i,j] <- rejections(ihw(pvals, 1:n, alphalist[i]))
        NumRej[11,i,j] <- ihw.oracle(pvals, 1:n, alphalist[i])
        NumRej[12,i,j] <- IF.oracle(pvals, 1:n, alphalist[i])
    }
    
    pi.formulas <- paste0("ns(x, df = ", 6:10, ")")
    mu.formulas <- paste0("ns(x, df = ", 6:10, ")")
    formulas <- expand.grid(pi.formulas, mu.formulas)
    
    res.AdaPT[[j]] <-
        tmp3 <- AdaPT.glm(x, pvals,
                  dist = beta.family(),
                  pi.formulas = formulas[, 1],
                  mu.formulas = formulas[, 2])

    NumRej[13,,j] <- res.AdaPT[[j]]$num_rej[1:num_alpha]

    print(paste0(j, "-th step finished!"))
    save(file = "data/GEOquery_res.RData", NumRej)
}

result <- list(NumRej = NumRej,
               res.AdaPT = res.AdaPT)

save(file = "data/GEOquery_res.RData", result)    

#####################################################################################
## Compare methods with random ordering
#####################################################################################

## repeat.times <- 100

## # gather results
## NumRej <- array(0, c(13, num_alpha, repeat.times))

## # methods: 1 SeqStep, 2 HingeExp, 3 ForwardStop, 4 Adaptive SeqStep, 5 BH, 6 Storey-BH, 7 Barber-Candes, 8 SABHA (step), 9 SABHA (ordered), 10 IHW, 11, IHW (oracle), 12 IF (oracle), 13 AdaPT
## for(j in 1:repeat.times){
##     if (j %% 5 == 1){
##         set.seed(floor((j + 1) / 5))
##     }
##     pvals <- output$pvals[sample(n,n)]

##     qhat_step <- Solve_q_step(pvals, tau, eps)
##     qhat_ordered <- Solve_q_ordered(pvals, tau, eps, ADMM_params)
##     for(i in 1:num_alpha){
##         NumRej[1,i,j] <- SeqStep(pvals, alpha = alphalist[i], C = 2)
##         NumRej[2,i,j] <- HingeExp(pvals, alpha = alphalist[i])
##         NumRej[3,i,j] <- ForwardStop(pvals, alpha = alphalist[i])
##         NumRej[4,i,j] <- length(Adaptive_SeqStep_method(pvals, alphalist[i], alphalist[i], tau))
##         NumRej[5,i,j] <- length(BH_method(pvals, alphalist[i]))
##         NumRej[6,i,j] <- length(Storey_method(pvals, alphalist[i], thr))
##         NumRej[7,i,j] <- length(BC_method(pvals, alphalist[i]))
##         NumRej[8,i,j] <- length(SABHA_method(pvals, qhat_step, alphalist[i], tau))
##         NumRej[9,i,j] <- length(SABHA_method(pvals, qhat_ordered, alphalist[i], tau))
##         NumRej[10,i,j] <- rejections(ihw(pvals, 1:n, alphalist[i]))
##         NumRej[11,i,j] <- ihw.oracle(pvals, 1:n, alphalist[i])
##         NumRej[12,i,j] <- IF.oracle(pvals, 1:n, alphalist[i])
##     }
##    
##     s0 <- rep(0.45, n)
##     pi.formulas <- paste0("ns(x, df = ", 6:15, ")")
##     mu.formulas <- paste0("ns(x, df = ", 6:15, ")")
##     formulas <- expand.grid(pi.formulas, mu.formulas)
##
##     res.AdaPT[[j]] <-
##         AdaPT.glm(x, pvals,
##                   dist = beta.family(),
##                   pi.formulas = formulas[, 1],
##                   mu.formulas = formulas[, 2],
##                   s0 = s0, delta = 1 / n,
##                   alphas = seq(0.01, 1, 0.01),
##                   period.up = 10)
##     NumRej[13,,j] <- res.AdaPT[[j]]$num_rej[1:num_alpha]
##     print(paste0(j, "-th step finished!"))
##     save(file = "data/GEOquery_random.RData", NumRej)        
## }
