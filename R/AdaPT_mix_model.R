################################################################
## EM algorithm and Local FDR estimation for the mixture model.
## 
## Components:
## (0) useful_functions.R
##     
## (1) exp_family.R:
##     A class for distributions of p-values:
##     log f(p; mu) = g(p)(eta(mu) - eta(mu*)) - (A(mu) - A(mu*)).
##     See Section 4.
##     Used for deriving the E-step and estimating local FDR.
##
##     Attributes:
##        g: transformation of p-value.
##        g.inv: inverse function of g.
##        eta: name of link function. Support all link functions in glm. See ?family for details.
##        mu.star: the parameter value for U([0, 1]).
##        A: log-partition function.
##        name (optional): name of the instance.
##        family (optional): the exponential family used for GLM/GAM. The family for which eta is the canonical link is recommended.
## 
## (2) safe_model.R
##     Safe versions of model fitting.
##     Maximally circumvent unexpected errors during the process. 
##
##     Required Input:
##         data: a data-frame
##     Output:
##         mod: model
##         fitv: fitted value
##     
##     Main modification:
##     Add an argument "alter.formulas" to incorporate a list of
##     alternative models when the fitting fails. For instance,
##     the default model to use is a GLM with natural spline basis
##     ns(x, df = 4) and the alternative models are GLMs with natural
##     spline basis ns(x, df = 3), ns(x, df = 2), ns(x, df = 1). This
##     effectively avoids the potential colinearity problem during the
##     process.
##
## (3) init_EM_mix.R
##     Initialize pix and mux. See Appendix A.2
##
##     Required Input:
##        x: covariate.
##        p: p-values.
##        pmin: lower threshold, i.e. s0(x) in AdaPT.
##        pmax: upper threshold, i.e. 1 - s0(x) in AdaPT.
##        dist: an "exp_family" object.
##        pi.fit.fun: a function to fit initial pix using \tilde{J_i}.
##        pi.fit.args: other arguments passed into pi.fit.fun.
##        mu.fit.fun: a model to fit initial mux using (x,min(p,1-p)).
##        mu.fit.args: other arguments passed into mu.fit.fun.
##     Output:
##        pix, mux: initial guess of pi(x) and mu(x).
##
## (4) Estep_mix.R
##     Compute E-step. See Appendix A.1.
##
##     Required Input:
##        x: covariate.
##        p: p-values.
##        pmin: lower threshold, i.e. s(x) in AdaPT.
##        pmax: upper threshold, i.e. 1 - s(x) in AdaPT.
##        pix: pi(x) in last step.
##        mux: mu(x) in last step.
##        dist: distribution family for p-values in "exp_family" class.
##     Output:
##        Hhat: conditional expectation of labels.
##        yhat: imputed y-values.
##        phat: g^{-1}(yhat).
##
## (5) Mstep_mix.R
##     Compute M-step. See Appendix A.1.
##
##     Required Input: 
##        x: covariate.
##        Hhat, yhat, phat: output of the E-step.
##        pi.fun: a function to fit pix using (x,Hhat).
##        mu.fun: a function to fit mux using (x,Hhat,yhat,phat).
##     Output:
##        pix: updated pi(x).
##        mux: updated mu(x).
## 
## (6) EM_mix.R
##     EM algorithm to fit a mixture model.
##     Combine init_EM_mix.R, Estep.R and Mstep.R together.
##
##     Required Input:
##        x: covariate
##        p: p-values.
##        pmin: lower threshold, i.e. s(x) in AdaPT.
##        pmax: upper threshold, i.e. 1 - s(x) in AdaPT.
##        dist: distribution family for p-values in "exp_family" class.
##        init.fun: a function for initialization when pix0 or mux0 is NULL.
##        init.params: extra parameters for init.fun (besided x, p, pmin, pmax, dist.)
##        Mstep.fun: a function for computing the M-step.
##        Mstep.params: extra parameters for init.fun (besided x and dist.)
##        params0: initial pi(x) and mu(x). Use initilization methods in init_EM_mix.R if pix=NULL or mux=NULL.
##        num.steps: maximal number of steps.
##        tol: tolerance for early stopping. Stop the procedure if ||mux.new - mux.old||_{\infty} < tol.
##     Output:
##        params: fitted parameters including pix and mux.
## 
## (7) lfdr_mix.R
##     Compute (over-estimated) Local fdr:
##     lfdr_i = f(1|x)/f(p_i|x).
##     
##     Required Input:
##        p: p-values.
##        dist: distribution family for p-values in "exp_family" class.
##        params: parameters including pix and mux.
##    
##     Output:
##        lfdr: local fdr for each p-value.
##
## (8) compute_threshold_mix.R
##     Compute the threshold curve s(x) given a level curve of local fdr. See Section 4.2.
##
##     Required Input:
##        dist: an "exp_family" object.
##        pix, mux: parameters.
##        lfdr.lev: target local-fdr level.
##     Output:
##        s: the threshold curve for a given local-fdr level.
##
## (9) EM_mix_select_model.R
##     Model selection based on partially masked data using information criteria.
##
##     Required Input:
##        x: covariate
##        p: p-values.
##        pmin: lower threshold, i.e. s(x) in AdaPT.
##        pmax: upper threshold, i.e. 1 - s(x) in AdaPT.
##        dist: distribution family for p-values in "exp_family" class.
##        cv.args: a list of arguments for cross-validation.
##        criterion: AIC or BIC.
##        ...: other arguments.
##     Output:
##        model: selected model. an element from cv.args.
################################################################

source("useful_functions.R")
source("AdaPT_mix_model/gen_EM_algo.R")
source("AdaPT_mix_model/exp_family.R")
source("AdaPT_mix_model/safe_model.R")
source("AdaPT_mix_model/init_EM_mix.R")
source("AdaPT_mix_model/Estep_mix.R")
source("AdaPT_mix_model/Mstep_mix.R")
source("AdaPT_mix_model/EM_mix.R")
source("AdaPT_mix_model/lfdr_mix.R")
source("AdaPT_mix_model/compute_threshold_mix.R")
source("AdaPT_mix_model/EM_mix_model_select.R")
