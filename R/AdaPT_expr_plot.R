################################################################
## Plots in Section 5
################################################################
source("AdaPT.R")
source("AdaPT_mix_summary.R")
source("expr_plot_fun.R")
source("useful_functions.R")
library("latex2exp")

#### Figure 1
n <- 50
set.seed(1)
pvals <- pnorm(rnorm(n,mean=c(seq(-3,0,length.out=n/5), rep(0,4*n/5))))
pvals <- pmax(pvals,0.005)

s <-.15
lambda <- .5
k <- 30
hlo <- .5
hhi <- .5
vlo <- 0
vhi <- 1

pvals <- pmax(pvals,0.01)
w <- 0.42 * exp(-0.043 * (0:n))
s <- 1

#### Figure 1
n <- 50
set.seed(1)
pvals <- pnorm(rnorm(n,mean=c(seq(-3,0,length.out=n/5), rep(0,4*n/5))))
pvals <- pmax(pvals,0.005)

s <-.15
lambda <- .5
k <- 30
hlo <- .5
hhi <- .5
vlo <- 0
vhi <- 1

pvals <- pmax(pvals,0.01)
w <- 0.42 * exp(-0.043 * (0:n))
s <- 1

pdf("../figs/Fig1_left.pdf", height = 3, width = 4)
par(mar = c(3, 2.5, 3.1, 0.1))
plot(pvals, ylim = 0:1, xlim = c(.5, n + hhi), xaxt = "n",
     yaxt = "n", xaxs = "i", yaxs = "i", type = "n",
     xlab = "", ylab = "")
title("AdaPT (Intermediate Stage)")
title(ylab = expression(paste("p-value ", p[i])), line = 1.5)
title(xlab = expression(paste("predictor ", x[i])), line = 1.5)
axis(2, at = c(0, 1), labels = c("0", "1"))
abline(h = 0:1)
polygon(x = 0.5 + c(0, n, n:0), y = c(vlo, vlo, s * rev(w)),
        col = "#FFDDDD", border = "red")
polygon(x = 0.5 + c(0, n, n:0), y = c(vhi, vhi, rev(1 - w)),
        col = "light blue", border = "blue")
which.cens <- which((pvals < s * w[-1]) | (pvals > 1 - w[-1]))
points((1:n)[-which.cens], pvals[-which.cens], pch = 16, cex = 0.8)
points(which.cens, pvals[which.cens], pch = 16, cex = 0.8,
       col = ifelse(pvals[which.cens] < 0.5, "red", "blue"))
legend("topright", bg = "white",lty = 1,col = c("blue", "red"),
       legend = c(expression(1 - s[t](x)),expression(s[t](x))))
box()
dev.off()

pdf("../figs/Fig1_right.pdf",height = 3, width = 4)
par(mar = c(3, 2.5, 3.1, 0.1))
plot(pvals, ylim = 0:1, xlim = c(.5, n+hhi), xaxt = "n",
     yaxt = "n", xaxs = "i", yaxs = "i", type = "n",
     xlab = "", ylab = "")
title("AdaPT (Analyst's View)")
title(ylab = expression(paste("p-value ", p[i])), line = 1.5)
title(xlab = expression(paste("predictor ", x[i])), line = 1.5)
axis(2, at = c(0,1), labels = c("0", "1"))
abline(h = 0:1)
polygon(x = 0.5 + c(0, n, n:0), y = c(vlo, vlo, s * rev(w)),
        col = "#FFDDDD", border = "red")
polygon(x = 0.5 + c(0, n, n:0), y = c(vhi, vhi, rev(1 - w)),
        col = "light blue", border = "blue")
which.cens <- which((pvals < s * w[-1]) | (pvals > 1 - w[-1]))
points((1:n)[-which.cens], pvals[-which.cens], pch = 16, cex = 0.8)
points(which.cens, pvals[which.cens], pch = 21,
       col = "purple", bg = "white", cex = 0.8)
points(which.cens, 1 - pvals[which.cens], pch = 21, col = "purple",
       bg = "white", cex = 0.8)
legend("topright",bg = "white",lty = 1,col = c("blue", "red"),
       legend = c(expression(1 - s[t](x)), expression(s[t](x))))
box()
dev.off()

#### Real data examples
methods <- c('SeqStep', 'HingeExp', 'ForwardStop', 'Ada. SeqStep', 'BH', 'Storey-BH', 'Barber-Candes', 'SABHA (step)', 'SABHA (ordered)', 'IHW', 'IHW (oracle)', 'IF (oracle)', 'AdaPT')
## Green for non-adaptive ordered testing procedures: SeqStep, Accumulation Test, Forward Stop and Adaptive SeqStep; Blue for non-adaptive BH-type methods: BH, Storey, BC; Orange for IHW and IF; Black for SABHA; Red for AdaPT
cols <- c('green', 'green', 'green', 'green', 'blue', 'blue','blue', 'black', 'black', 'orange', 'orange', 'orange', 'red')
ltys <- c(3, 5, 4, 1, 6, 2, 5, 1, 3, 2, 4, 6, 1)
pchs <- c(13, 12, 11, 10, 6, 5, 4, 2, 1, 6, 5, 7, 1)

## gene dosage experiment
load("../data/GEOQuery_res.RData")

NumRej1 <- result$NumRej[,,2]
obj1 <- result$res.AdaPT[[2]]
NumRej2 <- result$NumRej[,,1]
obj2 <- result$res.AdaPT[[1]]

load("../data/GEOquery_random.RData")
NumRej0 <- NumRej

pdf("../figs/Fig2.pdf", width = 10, height = 4)
par(mfrow = c(1, 3), oma = c(0, 0, 3, 0), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot.results(NumRej0, methods, "Random ordering", cols, ltys, pchs, c(0, 1500), "# of Rejections", TRUE, 1.2)
plot.results(NumRej1, methods, "Moderately informative ordering", cols, ltys, pchs, c(0, 4000), "# of Rejections", FALSE)
plot.results(NumRej2, methods, "Highly informative Ordering", cols, ltys, pchs, c(0, 4000), "", FALSE)
mtext("Number of Rejections (Gene/Drug Response)",outer=TRUE,cex=1.8,font=2)
dev.off()

pdf("../figs/Fig3_left.pdf")
info.plot(obj1, 0.05, "alpha = 0.05",
          xlab = "x = rank", disp.pmax = 0.1,
          legend.position = "topright")
dev.off()

pdf("../figs/Fig3_right.pdf")
info.plot(obj1, 0.1, "alpha = 0.1",
          xlab = "x = rank", disp.pmax = 0.1,
          legend.position = "topright")
dev.off()

pdf("../figs/Fig4_left.pdf")
info.plot(obj2, 0.05, "alpha = 0.05",
          xlab = "x = rank", disp.pmax = 0.3,
          legend.position = "topright")
dev.off()

pdf("../figs/Fig4_right.pdf")
info.plot(obj2, 0.1, "alpha = 0.1",
          xlab = "x = rank", disp.pmax = 0.3,
          legend.position = "topright")
dev.off()

corr.GEOquery.moderate <- corr.lfdr(obj1, num.steps = 100)

corr.GEOquery.high <- corr.lfdr(obj2, num.steps = 100)

corr.GEOquery <- list(corr.GEOquery.moderate,
                      corr.GEOquery.high)

pdf("../figs/Fig5.pdf", width = 7, height = 4.5)
plot.corr(corr.GEOquery, "Correlation of Estimated Local FDR", ylim = c(0.85, 1), cols = c("red", "blue"), ltys = 1:2, legend = c("moderately informative ordering", "highly informative ordering"), cex.lab = 1.3, cex.legend = 1.3, cex.main = 1.6)
dev.off()

#### Bottomly experiment
load("../data/Bottomly_res.RData")
NumRej <- result$NumRej
obj <- result$res.AdaPT

pdf("../figs/Fig11_left.pdf")
par(mfrow = c(1, 1), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot.results(NumRej, methods, "Number of Rejections (Bottomly)", cols, ltys, pchs, c(0, 5300), "# of Rejections")
dev.off()

pdf("../figs/Fig11_right.pdf")
info.plot(obj, 0.1, "Bottomly, alpha = 0.1",
          xlab = "x = log(baseMean)",
          disp.pmax = 0.1, legend.position = "topleft")
dev.off()

corr.Bottomly <- corr.lfdr(obj, num.steps = 100)

pdf("../figs/Fig11_middle.pdf")
plot.corr(corr.Bottomly, "Correlation of Estimated Local FDR (Bottomly)", ylim = c(0.98, 1), cex.main = 1.7, cex.lab = 1.7)
dev.off()

#### airway experiment
load("../data/airway_res.RData")
NumRej <- result$NumRej
obj <- result$res.AdaPT

pdf("../figs/Fig12_left.pdf")
par(mfrow = c(1, 1), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot.results(NumRej, methods, "Number of Rejections (airway)", cols, ltys, pchs, c(0, 12000), "# of Rejections")
dev.off()

pdf("../figs/Fig12_right.pdf")
info.plot(obj, 0.1, "airway, alpha = 0.1",
          xlab = "x = log(baseMean)",
          disp.pmax = 0.1, legend.position = "topleft")
dev.off()

corr.airway <- corr.lfdr(obj, num.steps = 100)

pdf("../figs/Fig12_middle.pdf")
plot.corr(corr.airway, "Correlation of Estimated Local FDR (airway)", ylim = c(0.98, 1), cex.main = 1.7, cex.lab = 1.7)
dev.off()

#### pasilla experiment
load("../data/pasilla_res.RData")
NumRej <- result$NumRej
obj <- result$res.AdaPT

pdf("../figs/Fig13_left.pdf")
par(mfrow = c(1, 1), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot.results(NumRej, methods, "Number of Rejections (pasilla)", cols, ltys, pchs, c(0, 2000), "# of Rejections")
dev.off()

pdf("../figs/Fig13_right.pdf")
info.plot(obj, 0.1, "pasilla, alpha = 0.1",
          xlab = "x = log(baseMean)",
          disp.pmax = 0.1, legend.position = "topleft")
dev.off()

corr.pasilla <- corr.lfdr(obj, num.steps = 100)

pdf("../figs/Fig13_middle.pdf")
plot.corr(corr.pasilla, "Correlation of Estimated Local FDR (pasilla)", ylim = c(0.98, 1), cex.main = 1.7, cex.lab = 1.7)
dev.off()

#### proteomics experiment
load("../data/proteomics_res.RData")
NumRej <- result$NumRej
obj <- result$res.AdaPT

pdf("../figs/Fig14_left.pdf")
par(mfrow = c(1, 1), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot.results(NumRej, methods, "Number of Rejections (SILAC)", cols, ltys, pchs, c(0, 1200), "# of Rejections")
dev.off()

pdf("../figs/Fig14_right.pdf")
info.plot(obj, 0.1, "SILAC, alpha = 0.1",
          xlab = "x = log(number of peptides)",
          disp.pmax = 0.1, legend.position = "topleft")
dev.off()

corr.proteomics <- corr.lfdr(obj, num.steps = 100)

pdf("../figs/Fig14_middle.pdf")
plot.corr(corr.proteomics, "Correlation of Estimated Local FDR (SILAC)", ylim = c(0.98, 1), cex.main = 1.7, cex.lab = 1.7)
dev.off()

#### Simulation 1 meta info
set.seed(1)
n <- 2500
x1 <- seq(-100, 100, length.out = 50)
x2 <- seq(-100, 100, length.out = 50)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")
pi.formula <- mu.formula <- "s(x1, x2)"
alpha.list <- seq(0.01, 0.3, 0.01)

## Truth
pdf("../figs/Fig6.pdf", width = 5, height = 1.75)
par(mfrow = c(1, 3), mar = c(2, 2, 2, 2))
n <- 2500
x1 <- seq(-100, 100, length.out = 50)
x2 <- seq(-100, 100, length.out = 50)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")
mains <- c("Circle in the middle",
           "Circle in the corner",
           "Thin ellipse")
for (i in 1:3){
    if (i == 1){
        H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
    } else if (i == 2){
        H0 <- apply(x, 1, function(coord){sum((coord - 65)^2) < 900})
    } else {
        shape.fun <- function(coord){
            transform.coord <- c(coord[1] + coord[2], coord[2] - coord[1])/sqrt(2)
            transform.coord[1]^2 / 100^2 + transform.coord[2]^2 / 15^2 < 1
        }
        H0 <- apply(x, 1, shape.fun)
    }
    plot.mask.2d(x, H0, cex = 0.5,
                 col.bg = "#A9A9A9", col.fg = "#000000",
                 main = mains[i], xaxt = "n", yaxt = "n")
    axis(side = 1, at = c(-100, 0, 100))
    axis(side = 2, at = c(-100, 0, 100))
}
dev.off()

#### Simulation 1
## Plots for simulation 1
methods <- c('BH', 'Storey-BH', 'Barber-Candes', 'SABHA', 'IHW', 'AdaPT')
## Blue for non-adaptive BH-type methods: BH, Storey, BC; Orange for IHW; Black for SABHA; Red for AdaPT
cols <- c('blue', 'blue','blue', 'black', 'orange', 'red')
## ltys <- c(1,2,3,2,3,1)
ltys <- 6:1
## pchs <- c(2,3,1,1,2,1)
pchs <- 6:1
inds <- c(1:5,7)

load("../data/simul1.RData")
titles <- c("Circle in the middle", "Circle in the corner",
            "Thin ellipse")
FDR.filename <- "../figs/Fig7_upper.pdf"
ylim <- c(0, 0.35)
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
    FDP <- result[[k]]$FDP[inds,]
    legend <- (k == 1)
    plot.results(FDP, methods, titles[k], cols, ltys, pchs,
                 ylim = ylim, ylab = "FDR",
                 legend = legend, cex.legend = 1.1)
}
dev.off()

power.filename <- "../figs/Fig7_lower.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
    power <- result[[k]]$power[inds,]
    legend <- FALSE
    plot.results(power, methods, titles[k], cols, ltys, pchs,
                 ylim = c(0, 1.05), ylab = "power",
                 legend = legend)
}
dev.off()

## Estimated local FDR in Case 1 (for illustration)
H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
mu <- ifelse(H0, 2, 0)
pvals <- pvals.gen(n, mu, 0, 0)

res.AdaPT <- AdaPT.gam(x, pvals, beta.family(),
                       pi.formula, mu.formula)

pdf("../figs/Fig8.pdf", width = 10, height = 3.7)
par(mfrow = c(1, 3))
plot.lfdr.2d(x, 1-res.AdaPT$lfdr[, 80],
             main = "Estimated local FDR (alpha = 0.5)",
             cex = 1.5, cex.main = 1.7,
             xaxt = "n", yaxt = "n", pch = 20)
plot.lfdr.2d(x, 1-res.AdaPT$lfdr[, 30],
             main = "Estimated local FDR (alpha = 0.3)",
             cex = 1.5, cex.main = 1.7,
             xaxt = "n", yaxt = "n", pch = 20)
plot.lfdr.2d(x, 1-res.AdaPT$lfdr[, 10],
             main = "Estimated local FDR (alpha = 0.1)",
             cex = 1.5, cex.main = 1.7,
             xaxt = "n", yaxt = "n", pch = 20)
dev.off() 


#### Simulation 2 meta info
m <- 100
n <- 2000

x <- matrix(runif(n * m), n, m)
pi1 <- 0.3

beta.pi <- c(3, 3, rep(0, m-2))
beta0.pi <- uniroot(function(b){
    mean(inv.logit(x %*% beta.pi + b)) - pi1
}, c(-100, 100))$root
pi <- inv.logit(x %*% beta.pi + beta0.pi)

beta.mu <- c(2, 2, rep(0, m-2))
beta0.mu <- 0
mu <- 1 / pmax(1, x %*% beta.mu + beta0.mu)

pdf("../figs/Fig9.pdf", width = 10, height = 3.7)
par(mfrow = c(1, 2), cex = 1.5, mar = c(4, 5, 1, 3))
hist(pi, 20, xlab = TeX("$\\pi_{1i}$"), main = "")
hist(1 / mu, 20, xlab = TeX("$\\mu_{i}$"), main = "")
dev.off()


## Simulation 2
## Plots for simulation 2
methods <- c('BH', 'Storey-BH', 'Barber-Candes', 'AdaPT', 'AdaPT (oracle)')
## Blue for non-adaptive BH-type methods: BH, Storey, BC; Orange for IHW; Black for SABHA; Red for AdaPT
cols <- c('blue', 'blue','blue', 'red', 'red')
ltys <- c(1,2,3,2,1)
pchs <- c(2,3,1,2,1)   

load("../data/simul2.RData")

filename <- "../figs/Fig10.pdf"
pdf(filename, width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5, 1, 3), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
FDP <- result$FDP
plot.results(FDP, methods, '', cols, ltys, pchs,
             ylim = c(0, 0.35), ylab = "FDR",
             cex.legend = 1.2)
power <- result$power
plot.results(power, methods, '', cols, ltys, pchs,
             ylim = c(0, 0.45), ylab = "power",
             legend = FALSE)
dev.off()
