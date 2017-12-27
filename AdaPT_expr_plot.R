################################################################
## Plots in Section 5
################################################################
source("AdaPT_mix.R")
source("AdaPT_mix_summary.R")
source("expr_plot_fun.R")

#### Real data examples
methods <- c('SeqStep', 'HingeExp', 'ForwardStop', 'Ada. SeqStep', 'BH', 'Storey-BH', 'Barber-Candes', 'SABHA (step)', 'SABHA (ordered)', 'IHW', 'IHW (oracle)', 'IF (oracle)', 'AdaPT')
## Black for non-adaptive ordered testing procedures: SeqStep, Accumulation Test, Forward Stop and Adaptive SeqStep; Blue for non-adaptive BH-type methods: BH, Storey, BC; Orange for IHW and IF; Green for SABHA; Red for AdaPT (old and new)
cols <- c('black', 'black', 'black', 'black', 'blue', 'blue','blue', 'green', 'green', 'orange', 'orange', 'orange', 'red')
ltys <- c(4,3,2,1,1,2,3,2,1,3,2,1,1)
pchs <- c(4,3,2,1,2,3,1,1,4,2,3,5,1)


## gene dosage experiment
load("data/GEOQuery_res.RData")

inds <- c(1:11, 12, 14)
NumRej1 <- result$NumRej[inds,,3]
obj1 <- result$res.AdaPT.beta[[3]]
NumRej2 <- result$NumRej[inds,,2]
obj2 <- result$res.AdaPT.beta[[2]]
NumRej3 <- result$NumRej[inds,,1]
obj3 <- result$res.AdaPT.beta[[1]]

#### Preprocessing
## summary.NumRej <- array(0,c(17,30,100))
## for (seed in 0:19){
##     filename <- paste0("data/GEOquery_random_", seed, ".RData")
##     load(filename)
##     tmp.inds <- (seed * 5 + 1):((seed + 1) * 5)
##     summary.NumRej[,,tmp.inds] <- NumRej
## }
## NumRej <- apply(summary.NumRej, MARGIN = c(1, 2), mean)
## save(file = "data/GEOquery_random.RData", NumRej)

load("data/GEOquery_random.RData")
NumRej0 <- NumRej[inds,]

pdf("figs/GEOquery_rejs.pdf", width = 10, height = 8)
par(mfrow = c(2, 2), oma = c(0, 0, 3, 0), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot.results(NumRej0, methods, "Random ordering", cols, ltys, pchs, c(0, 1500), "# of Rejections", TRUE, 0.95)
plot.results(NumRej1, methods, "Original ordering", cols, ltys, pchs, c(0, 1500), "", FALSE)
plot.results(NumRej2, methods, "Moderately informative ordering", cols, ltys, pchs, c(0, 4000), "# of Rejections", FALSE)
plot.results(NumRej3, methods, "Highly informative Ordering", cols, ltys, pchs, c(0, 4000), "", FALSE)
mtext("Number of Rejections (Gene/Drug Response)",outer=TRUE,cex=1.8,font=2)
dev.off()

pdf("figs/GEOquery_original_info_20.pdf")
info.plot(obj1, 0.2, "alpha = 0.2", "x = rank", disp.pmax = 0.1, legend.position = "topright")
dev.off()

pdf("figs/GEOquery_original_info_25.pdf")
info.plot(obj1, 0.25, "alpha = 0.25", "x = rank", disp.pmax = 0.1, legend.position = "topright")
dev.off()

pdf("figs/GEOquery_moderate_info_05.pdf")
info.plot(obj2, 0.05, "alpha = 0.05", "x = rank", disp.pmax = 0.1, legend.position = "topright")
dev.off()

pdf("figs/GEOquery_moderate_info_10.pdf")
info.plot(obj2, 0.1, "alpha = 0.1", "x = rank", disp.pmax = 0.1, legend.position = "topright")
dev.off()

pdf("figs/GEOquery_high_info_05.pdf")
info.plot(obj3, 0.05, "alpha = 0.05", "x = rank", disp.pmax = 0.3, legend.position = "topright")
dev.off()

pdf("figs/GEOquery_high_info_10.pdf")
info.plot(obj3, 0.1, "alpha = 0.1", "x = rank", disp.pmax = 0.3, legend.position = "topright")
dev.off()

corr.GEOquery.original <- corr.lfdr(obj1, num.steps = 100,
                                    dist = beta.family())

corr.GEOquery.moderate <- corr.lfdr(obj2, num.steps = 100,
                                    dist = beta.family())

corr.GEOquery.high <- corr.lfdr(obj3, num.steps = 100,
                                dist = beta.family())

corr.GEOquery <- list(corr.GEOquery.original,
                      corr.GEOquery.moderate,
                      corr.GEOquery.high)

pdf("figs/GEOquery_corr.pdf")
plot.corr(corr.GEOquery, "Correlation of Estimated Local FDR", ylim = c(0.85, 1), cols = c("black", "red", "blue"), ltys = 1:3, legend = c("original ordering", "moderately informative ordering", "highly informative ordering"), cex.lab = 1.3, cex.legend = 1.3, cex.main = 1.6)
dev.off()

#### Bottomly experiment
load("data/Bottomly_res.RData")
inds <- 1:13
NumRej <- result$NumRej[inds,]
obj <- result$res.AdaPT.beta

pdf("figs/Bottomly_rejs.pdf")
par(mfrow = c(1, 1), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot.results(NumRej, methods, "Number of Rejections (Bottomly)", cols, ltys, pchs, c(0, 5300), "# of Rejections")
dev.off()

pdf("figs/Bottomly_info_10.pdf")
info.plot(obj, 0.1, "Bottomly, alpha = 0.1", "x = log(baseMean)", disp.pmax = 0.1, legend.position = "topleft")
dev.off()

pdf("figs/Bottomly_info_05.pdf")
info.plot(obj, 0.05, "Bottomly, alpha = 0.05", "x = log(baseMean)", disp.pmax = 0.1, legend.position = "topleft")
dev.off()

corr.Bottomly <- corr.lfdr(obj, num.steps = 100,
                           dist = beta.family())

pdf("figs/Bottomly_corr.pdf")
plot.corr(corr.Bottomly, "Correlation of Estimated Local FDR (Bottomly)", ylim = c(0.98, 1), cex.main = 1.7, cex.lab = 1.7)
dev.off()

#### airway experiment
load("data/airway_res.RData")
inds <- 1:13
NumRej <- result$NumRej[inds,]
obj <- result$res.AdaPT.beta

pdf("figs/airway_rejs.pdf")
par(mfrow = c(1, 1), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot.results(NumRej, methods, "Number of Rejections (airway)", cols, ltys, pchs, c(0, 12000), "# of Rejections")
dev.off()

pdf("figs/airway_info_10.pdf")
info.plot(obj, 0.1, "airway, alpha = 0.1", "x = log(baseMean)", disp.pmax = 0.1, legend.position = "topleft")
dev.off()

pdf("figs/airway_info_05.pdf")
info.plot(obj, 0.05, "airway, alpha = 0.05", "x = log(baseMean)", disp.pmax = 0.1, legend.position = "topleft")
dev.off()

corr.airway <- corr.lfdr(obj, num.steps = 100,
                         dist = beta.family())

pdf("figs/airway_corr.pdf")
plot.corr(corr.airway, "Correlation of Estimated Local FDR (airway)", ylim = c(0.98, 1), cex.main = 1.7, cex.lab = 1.7)
dev.off()

#### proteomics experiment
load("data/proteomics_res.RData")
inds <- 1:13
NumRej <- result$NumRej[inds,]
obj <- result$res.AdaPT.beta

pdf("figs/proteomics_rejs.pdf")
par(mfrow = c(1, 1), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot.results(NumRej, methods, "Number of Rejections (SILAC)", cols, ltys, pchs, c(0, 1200), "# of Rejections")
dev.off()

pdf("figs/proteomics_info_10.pdf")
info.plot(obj, 0.1, "proteomics, alpha = 0.1", "x = log(number of peptides)", disp.pmax = 0.1, legend.position = "topleft")
dev.off()

pdf("figs/proteomics_info_05.pdf")
info.plot(obj, 0.05, "proteomics, alpha = 0.05", "x = log(number of peptides)", disp.pmax = 0.1, legend.position = "topleft")
dev.off()

corr.proteomics <- corr.lfdr(obj, num.steps = 100,
                             dist = beta.family())

pdf("figs/proteomics_corr.pdf")
plot.corr(corr.proteomics, "Correlation of Estimated Local FDR (SILAC)", ylim = c(0.98, 1), cex.main = 1.7, cex.lab = 1.7)
dev.off()

#### pasilla experiment
load("data/pasilla_res.RData")
inds <- 1:13
NumRej <- result$NumRej[inds,]
obj <- result$res.AdaPT.beta

pdf("figs/pasilla_rejs.pdf")
par(mfrow = c(1, 1), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
plot.results(NumRej, methods, "Number of Rejections (pasilla)", cols, ltys, pchs, c(0, 2000), "# of Rejections")
dev.off()

pdf("figs/pasilla_info_10.pdf")
info.plot(obj, 0.1, "pasilla, alpha = 0.1", "x = log(baseMean)", disp.pmax = 0.1, legend.position = "topleft")
dev.off()

pdf("figs/pasilla_info_05.pdf")
info.plot(obj, 0.05, "pasilla, alpha = 0.05", "x = log(baseMean)", disp.pmax = 0.1, legend.position = "topleft")
dev.off()

corr.pasilla <- corr.lfdr(obj, num.steps = 100,
                          dist = beta.family())

pdf("figs/pasilla_corr.pdf")
plot.corr(corr.pasilla, "Correlation of Estimated Local FDR (pasilla)", ylim = c(0.98, 1), cex.main = 1.7, cex.lab = 1.7)
dev.off()

#### Simulation results
exprs <- data.frame(rho = c(0, 5, -5, 8, 8),
                    type = c(0, 0, 1, 2, 3))
    
## Read Data on FDR and Power
for (i in 1:nrow(exprs)){
    rho <- as.numeric(exprs$rho[i])
    type <- as.numeric(exprs$type[i])
    data.convex.sum <- list()
    for (k in 1:3){
        data.convex.sum[[k]] <- list(FDP = list(), power = list())
    }
    for (k in 1:3){
        for (j in 1:6){
            data.convex.sum[[k]]$FDP[[j]] <- rep(0, 30)
            data.convex.sum[[k]]$power[[j]] <- rep(0, 30)    
        }
    }

    for (seed in 1:10){
        filename <- paste0("data/data_convex_rho_", rho, "_type_", type, "_seed_", seed, ".RData")
        load(filename)
        for (k in 1:3){
            for (j in 1:6){
                data.convex.sum[[k]]$FDP[[j]] <-
                    data.convex.sum[[k]]$FDP[[j]] +
                        data.convex[[k]]$FDP[[j]]
                data.convex.sum[[k]]$power[[j]] <-
                    data.convex.sum[[k]]$power[[j]] +
                        data.convex[[k]]$power[[j]]
            }
        }
    }
    for (k in 1:3){
        for (j in 1:6){
            data.convex.sum[[k]]$FDP[[j]] <-
                data.convex.sum[[k]]$FDP[[j]] / 10
            data.convex.sum[[k]]$power[[j]] <-
                data.convex.sum[[k]]$power[[j]] / 10
        }
    }

    result <- list()
    for (k in 1:3){
        tmp <- list()
        tmp$FDP <- Reduce(rbind, data.convex.sum[[k]]$FDP)
        rownames(tmp$FDP) <- NULL
        tmp$power <- Reduce(rbind, data.convex.sum[[k]]$power)
        rownames(tmp$power) <- NULL
        result[[k]] <- tmp
    }

    if (rho != 8){
        filename <- paste0("data/simul_", rho, ".RData")
    } else {
        filename <- paste0("data/simul_", rho, "_type_",
                           type, ".RData")
    }

    save(file = filename, result)
}

methods <- c("Storey-BH", "Barber-Candes", "IHW", "AdaPT")
cols <- c("blue", "blue", "orange", "red")
ltys <- c(2, 3, 1, 1)
pchs <- c(2, 3, 1, 1)


for (i in 1:nrow(exprs)){
    rho <- as.numeric(exprs$rho[i])
    type <- as.numeric(exprs$type[i])
    if (rho != 8){
        filename <- paste0("data/simul_", rho, ".RData")
    } else {
        filename <- paste0("data/simul_", rho, "_type_",
                           type, ".RData")
    }
    load(filename)
    titles <- c("Circle in the middle",
        "Circle in the corner",
        "Thin ellipse")    

    if (rho == 0){
        FDR.filename <- "figs/simul_low_FDR.pdf"
        ylim <- c(0, 0.32)
    } else if (rho == 5){
        FDR.filename <- "figs/simul_low_FDR_5.pdf"
        ylim <- c(0, 0.48)
    } else if (rho == -5){
        FDR.filename <- "figs/simul_low_FDR_-5.pdf"
        ylim <- c(0, 0.32)
    } else if (rho == 8){
        FDR.filename <- paste0("figs/simul_low_FDR_8_", type, ".pdf")
        ylim <- c(0, 0.32)
    }
    pdf(FDR.filename, width = 9, height = 2.8)
    par(mfrow = c(1, 3), mar = c(4, 4, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
    for (k in 1:3){
        FDP <- result[[k]]$FDP[2:5, ]
        legend <- (k == 1)
        plot.results(FDP, methods, titles[k], cols, ltys, pchs,
                     ylim = ylim, ylab = "FDR",
                     legend = legend)
    }
    dev.off()
    
    power.filename <-
        switch(as.character(rho),
               "0" = "figs/simul_low_power.pdf",
               "5" = "figs/simul_low_power_5.pdf",
               "-5" = "figs/simul_low_power_-5.pdf",
               "8" = paste0("figs/simul_low_power_8_", type, ".pdf"))
    pdf(power.filename, width = 9, height = 2.8)
    par(mfrow = c(1, 3), mar = c(4, 4, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
    for (k in 1:3){
        power <- result[[k]]$power[2:5, ]
        legend <- (k == 1)
        plot.results(power, methods, titles[k], cols, ltys, pchs,
                     ylim = c(0, 1.05), ylab = "power",
                     legend = legend)
    }
    dev.off()
}
