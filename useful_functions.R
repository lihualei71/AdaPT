logit <- function(x){
    log(x / (1 - x))
}

inv.logit <- function(x){
    exp(x) / (1 + exp(x))
}

fdp.hat <- function(pvals, s) {
    ## Function to calculate FDPhat (1 + At) / (Rt v 1).
    ## Inputs:
    ##    pvals: p-values
    ##    s: s(x)
    ## Output:
    ##    fdphat
    A <- sum(pvals > 1 - s)
    R <- sum(pvals <= s)
    return((1 + A) / max(R, 1))
}

find.newname <- function(names.vec){
    name <- "aaa"
    while (name %in% names.vec){
        name <- paste0(name, "a")
    }
    return(name)
}

complete.formulas <- function(params, response.name){
    if (!is.null(params[["formula"]])){
        completed.formula <- as.formula(
            paste(response.name, "~", params[["formula"]]))
        params[["formula"]] <- completed.formula
    }
    if (!is.null(params[["alter.formulas"]])){
        completed.alter.formulas <-
            sapply(params[["alter.formulas"]], function(formula){
                as.formula(paste(response.name, "~", formula))
            })
        params[["alter.formulas"]] <- completed.alter.formulas
    }
    return(params)
}

## Generate p-values
library("MASS")
pvals.gen <- function(n, mu, rho, type, x = NULL){
    if (rho >= 0){
        z1 <- rnorm(n)
        z2 <- rnorm(1)
        z <- z1 * sqrt(1 - rho^2) + z2 * rho + mu
    } else {
        if (type == 1){
            new.rho <- rho / n
            Sigma <- diag(rep(1 - new.rho, n)) + new.rho
            z <- mvrnorm(1, mu, Sigma)
        } else if (type == 2){
            distance <- as.matrix(dist(x))
            min.dist <- min(distance[distance > 0])
            gamma <- -log(rho) / min.dist
            Sigma <- exp(-gamma * mat)
            z <- mvrnorm(1, mu, Sigma)
        } else if (type == 3){
            distance <- as.matrix(dist(x))
            min.dist <- min(distance[distance > 0])
            gamma <- -log(rho) / min.dist^2
            Sigma <- exp(-gamma * mat^2)
            z <- mvrnorm(1, mu, Sigma)
        }

    }
    pvals <- 1 - pnorm(z)
    return(pvals)
}
