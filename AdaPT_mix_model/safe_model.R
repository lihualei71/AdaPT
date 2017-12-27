################################################################
## Safe versions of model fitting.
## Maximally circumvent unexpected errors during the process.
##
## Required Input:
##     data: a data-frame
## Output:
##     mod: model
##     fitv: fitted value
## 
## Main modification:
## Add an argument "alter.formulas" to incorporate a list of
## alternative models when the fitting fails. For instance,
## the default model to use is a GLM with natural spline basis
## ns(x, df = 4) and the alternative models are GLMs with natural
## spline basis ns(x, df = 3), ns(x, df = 2), ns(x, df = 1). This
## effectively avoids the potential colinearity problem during the
## process.
################################################################

library("mgcv")
library("HDtweedie")

safe.fit <- function(algo, data.args,
                     algo.args, alter.args = NULL,
                     ...){
    ## Safe version of arbitrary model-fitting algorithm
    oldw <- getOption("warn")
    options(warn = -1)

    extra.args <- list(...)
    args <- c(data.args, algo.args, extra.args)
    mod <- try(do.call(algo, args))
    
    if (class(mod)[1] != "try-error" || is.null(alter.args)){
        options(warn = oldw)        
        return(list(mod = mod, args = args))
    } 

    rm(args)
    gc()
    for (algo.args in alter.args){
        args <- c(data.args, algo.args, extra.args)
        mod <- try(do.call(algo, args))

        if (class(mod)[1] != "try-error"){
            options(warn = oldw)            
            return(list(mod = mod, args = args))
        }

        rm(args)
        gc()
    }

    options(warn = oldw)
    stop("model fitting fails!")
}

safe.glm <- function(formula, data,
                     alter.formulas = NULL,
                     family = gaussian(),
                     ...){
    ## Safe GLM
    glm <- stats::glm
    algo <- function(formula, data, family, ...){
        if (family$link %in% c("inverse", "log")){
            mod <- try(glm(formula = formula, data = data,
                           family = family, ...),
                       silent = TRUE)
            if (class(mod)[1] != "try-error"){
                return(mod)
            }
            tmp.mat <- model.matrix(formula, data = data)
            p <- ncol(tmp.mat) - 1
            start <- c(1, rep(0, p))
            mod <- glm(formula = formula, data = data,
                       family = family, start = start, ...)
        } else {
            mod <- glm(formula = formula, data = data,
                       family = family, ...)
        }
        return(mod)
    }

    data.args <- c(list(data = data))
    algo.args <- c(list(formula = formula, family = family))
    alter.args <- lapply(alter.formulas, function(formula){
        c(list(formula = formula, family = family))
    })

    result <- safe.fit(algo, data.args, algo.args, alter.args,
                       ...)
    mod <- result$mod
    fitv <- predict(mod, type = "response")
    
    return(list(mod = mod, fitv = fitv))
}

safe.gam <- function(formula, data,
                     alter.formulas = NULL,
                     family = gaussian(),
                     ...){
    ## Safe GAM
    gam <- mgcv::gam

    algo <- function(formula, data, family, ...){
        if (family$link %in% c("inverse", "log")){
            mod <- try(gam(formula = formula, data = data,
                           family = family, ...),
                       silent = TRUE)
            if (class(mod)[1] != "try-error"){
                return(mod)
            }
            tmp.mat <- model.matrix(formula, data = data)
            p <- ncol(tmp.mat) - 1
            start <- c(1, rep(0, p))
            mod <- gam(formula = formula, data = data,
                       family = family, start = start, ...)
        } else {
            mod <- gam(formula = formula, data = data,
                       family = family, ...)
        }
        return(mod)
    }

    data.args <- c(list(data = data))    
    algo.args <- c(list(formula = formula, family = family))
    alter.args <- lapply(alter.formulas, function(formula){
        c(list(formula = formula, family = family))
    })

    result <- safe.fit(algo, data.args, algo.args, alter.args,
                       ...)
    mod <- result$mod
    fitv <- predict(mod, type = "response")
    
    return(list(mod = mod, fitv = fitv))
}

safe.glmnet <- function(x, y,
                        family = gaussian(),
                        ...){
    ## Safe GLMnet
    if (class(family)[1] == "family"){
        family <- family$family
    }

    if (family %in% c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian")){
        algo <- function(...)
            glmnetUtils::cv.glmnet(..., family = family)
    } else if (family == "Gamma"){
        algo <- function(...)
            HDtweedie::cv.HDtweedie(..., p = 2)
    }

    data.args <- list(x = x, y = y)
    algo.args <- list()
    alter.args <- NULL
    result <- safe.fit(algo, data.args, algo.args, alter.args, ...)    
    mod <- result$mod
    fitv <- as.numeric(predict(mod, newx = x, s = "lambda.min", type = "response"))
    return(list(mod = mod, fitv = fitv))
}

safe.logistic.glm <- function(formula, data,
                              alter.formulas = NULL,
                              ...){
    safe.glm(formula, data, alter.formulas, 
             family = binomial(),
             ...)
}

safe.logistic.gam <- function(formula, data,
                              alter.formulas = NULL,
                              ...){
    safe.gam(formula, data, alter.formulas, 
             family = binomial(),
             ...)
}

safe.logistic.glmnet <- function(x, y,
                                 ...){
    safe.glmnet(x, y,
                family = "binomial",
                ...)
}

safe.gaussian.glm <- function(formula, data,
                              alter.formulas = NULL,
                              ...){
    safe.glm(formula, data, alter.formulas, 
             family = gaussian(),
             ...)
}

safe.gaussian.gam <- function(formula, data,
                              alter.formulas = NULL,
                              ...){
    safe.gam(formula, data, alter.formulas, 
             family = gaussian(),
             ...)
}

safe.gaussian.glmnet <- function(x, y,
                                 ...){
    safe.glmnet(x, y,
                family = "gaussian",
                ...)
}

safe.Gamma.glm <- function(formula, data,
                           alter.formulas = NULL,
                           ...){
    safe.glm(formula, data, alter.formulas, 
             family = Gamma(),
             ...)
}

safe.Gamma.gam <- function(formula, data,
                           alter.formulas = NULL,
                           ...){
    safe.gam(formula, data, alter.formulas, 
             family = Gamma(),
             ...)
}

safe.Gamma.glmnet <- function(x, y,
                              ...){
    safe.glmnet(x, y,
                family = "Gamma",
                ...)
}
