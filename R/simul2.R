library("dplyr")
source('summarize_methods.R')
source("useful_functions.R")
source("AdaPT.R")

repeat.times <- as.numeric(Sys.getenv("times"))
seed <- as.numeric(Sys.getenv("seed"))

simul2 <- function(x, mu, pi,
                   alpha.list,
                   repeat.times){
    n <- length(mu)
    m <- length(alpha.list)
    summary.FDP <- list()
    summary.power <- list()
    for (j in 1:5){
        summary.FDP[[j]] <- matrix(rep(0, repeat.times * m),
                                   ncol = repeat.times)
        summary.power[[j]] <- matrix(rep(0, repeat.times * m),
                                     ncol = repeat.times)
    }
    xx <- data.frame(x[, 1:2])
    names(xx) <- c("x1", "x2")
    
    for (i in 1:repeat.times){
        H0 <- as.logical(ifelse(runif(n) < pi, 1, 0))
        y <- ifelse(H0, rexp(n, 1/mu), rexp(n, 1))
        pvals <- exp(-y)
        
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

        res.AdaPT <- try(
            AdaPT.glmnet(x, pvals, beta.family(),
                         glmnet.args = list(nlambda = 30),
                         alphas = alpha.list))
        if (class(res.AdaPT) != "try-error"){
            AdaPT.result <- summary.AdaPT(res.AdaPT, H0, pvals)
            summary.FDP[[4]][, i] <- AdaPT.result[, 2]
            summary.power[[4]][, i] <- AdaPT.result[, 3]
        }

        res.AdaPT.oracle <- try(
            AdaPT.glm(xx, pvals, beta.family(),
                      pi.formulas = "x1 + x2",
                      mu.formulas = "x1 + x2",
                      alphas = alpha.list))
        if (class(res.AdaPT.oracle) != "try-error"){
            AdaPT.oracle.result <- summary.AdaPT(res.AdaPT.oracle, H0, pvals)
            summary.FDP[[5]][, i] <- AdaPT.oracle.result[, 2]
            summary.power[[5]][, i] <- AdaPT.oracle.result[, 3]
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

#### Simulation 2
m <- 100
n <- 2000
alpha.list <- seq(0.01, 0.3, 0.01)
output.filename <- paste0("../data/simul2_seed_", seed, ".RData")

x <- matrix(runif(n * m), n, m)
pi1 <- 0.3

beta.pi <- c(3, 3, rep(0, m-2))
beta0.pi <- uniroot(function(b){
    mean(inv.logit(x %*% beta.pi + b)) - pi1
}, c(-100, 100))$root
pi <- inv.logit(x %*% beta.pi + beta0.pi)

beta.mu <- c(2, 2, rep(0, m-2))
beta0.mu <- 0
mu <- pmax(1, x %*% beta.mu + beta0.mu)

result <- simul2(x, mu, pi, alpha.list, repeat.times)

save(file = output.filename, result)
