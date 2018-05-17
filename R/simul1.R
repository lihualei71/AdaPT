source("AdaPT.R")
source('summarize_methods.R')
source("useful_functions.R")
source("All_q_est_functions.R")

repeat.times <- as.numeric(Sys.getenv("times"))
seed <- as.numeric(Sys.getenv("seed"))

simul1 <- function(x, mu, H0, 
                   pi.formula, mu.formula, 
                   alpha.list,
                   repeat.times,
                   ...){
    n <- length(mu)
    n1 <- n2 <- floor(sqrt(n))
    m <- length(alpha.list)
    summary.FDP <- list()
    summary.power <- list()
    for (j in 1:7){
        summary.FDP[[j]] <- matrix(rep(0, repeat.times * m),
                                   ncol = repeat.times)
        summary.power[[j]] <- matrix(rep(0, repeat.times * m),
                                     ncol = repeat.times)
    }

    ## Settings for SABHA
    rhoG <- 0.9037629 # max_k ||(D_G+)_k||_2
    TV_bd <- 10
    tau <- 0.5
    eps <- 0.1
    ADMM_params <- c(10^2, 10^3, 2, 5000, 1e-3) 
    
    for (i in 1:repeat.times){
        z <- rnorm(n) + mu
        pvals <- 1 - pnorm(z)
        
        BH.result <- summary.BH(pvals, H0, alpha.list)
        summary.FDP[[1]][, i] <- BH.result[, 2]
        summary.power[[1]][, i] <- BH.result[, 3]

        storey.result <- summary.storey(pvals, H0,
                                        alpha.list = alpha.list)
        summary.FDP[[2]][, i] <- storey.result[, 2]
        summary.power[[2]][, i] <- storey.result[, 3]

        BC.result <- summary.BC(pvals, H0, alpha.list)
        summary.FDP[[3]][, i] <- BC.result[, 2]
        summary.power[[3]][, i] <- BC.result[, 3]

        IHW.result <- summary.IHW(pvals, H0, alpha.list)
        summary.FDP[[4]][, i] <- IHW.result[, 2]
        summary.power[[4]][, i] <- IHW.result[, 3]

        qhat <- Solve_q_TV_2dim(matrix(pvals, n1, n2), tau, eps, TV_bd, ADMM_params)
        SABHA.result <- summary.SABHA(pvals, H0, qhat, alpha.list)
        summary.FDP[[5]][, i] <- SABHA.result[, 2]
        summary.power[[5]][, i] <- SABHA.result[, 3]

        SABHA.factor <- 1 + 1 / eps / (1 - tau) / sqrt(n) +
            2 * rhoG * TV_bd * sqrt(log(n)) / eps^2 / n
        SABHA.result2 <- summary.SABHA(pvals, H0, qhat,
                                       alpha.list / SABHA.factor)
        summary.FDP[[6]][, i] <- SABHA.result2[, 2]
        summary.power[[6]][, i] <- SABHA.result2[, 3]
        
        res.AdaPT <- try(AdaPT.gam(x, pvals,
                                   dist = beta.family(),
                                   pi.formula = pi.formula,
                                   mu.formula = mu.formula,
                                   alphas = alpha.list))
        if (class(res.AdaPT) != "try-error"){
            AdaPT.result <- summary.AdaPT(res.AdaPT, H0, pvals)
            summary.FDP[[7]][, i] <- AdaPT.result[, 2]
            summary.power[[7]][, i] <- AdaPT.result[, 3]
        }

        print(paste0(i, "-th step finishes!"))
    }
    avg.FDP <- lapply(summary.FDP, function(FDP){
        apply(FDP, 1, function(x){mean(x, na.rm = TRUE)})
    })
    avg.power <- lapply(summary.power, function(power){
        apply(power, 1, function(x){mean(x, na.rm = TRUE)})
    })
    return(list(FDP = avg.FDP, power = avg.power))
}

####### Generate x
set.seed(seed)
n <- 2500
x1 <- seq(-100, 100, length.out = 50)
x2 <- seq(-100, 100, length.out = 50)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")
pi.formula <- mu.formula <- "s(x1, x2)"
alpha.list <- seq(0.01, 0.3, 0.01)
output.filename <- paste0("../data/simul1_seed_", seed, ".RData")

## Case 1: a circle in the center
H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
mu <- ifelse(H0, 2, 0)

result1 <- simul1(x, mu, H0, rho, type,
                  pi.formula = pi.formula,
                  mu.formula = mu.formula,
                  alpha.list = alpha.list,
                  repeat.times = repeat.times)
result <- list(result1)
save(file = output.filename, result)

## Case 2: a circle in the corner
H0 <- apply(x, 1, function(coord){sum((coord - 65)^2) < 900})
mu <- ifelse(H0, 2, 0)

result2 <- simul1(x, mu, H0, rho, type,
                  pi.formula = pi.formula,
                  mu.formula = mu.formula,
                  alpha.list = alpha.list,
                  repeat.times = repeat.times)
result <- list(result1, result2)
save(file = output.filename, result)


## Case 3: a thin ellipsoid
shape.fun <- function(coord){
    transform.coord <- c(coord[1] + coord[2], coord[2] - coord[1])/sqrt(2)
    transform.coord[1]^2 / 100^2 + transform.coord[2]^2 / 15^2 < 1
}
H0 <- apply(x, 1, shape.fun)
mu <- ifelse(H0, 2, 0)

result3 <- simul1(x, mu, H0, rho, type,
                  pi.formula = pi.formula,
                  mu.formula = mu.formula,
                  alpha.list = alpha.list,
                  repeat.times = repeat.times)

## Save data
result <- list(result1, result2, result3)
save(file = output.filename, result)
