plot.results <- function(vals, methods, title,
                         cols, ltys, pchs,
                         ylim, ylab,
                         legend = TRUE,
                         cex.legend = 1){
### Plot results for each method
    alphalist <- seq(0.01, 0.3, 0.01)
    plot(0:1, 0:1, type = 'n',
         xlim = range(alphalist), ylim = ylim,
         xlab = expression(paste('Target FDR level ',alpha)),
         ylab = ylab,
         main = title, axes = FALSE)
    axis(side = 1, at = c(0, 0.1, 0.2, 0.3))
    axis(side = 2)
    alpha_pt = c(10, 20, 30)
    for (i in 1:length(methods)){
        points(alphalist, vals[i, ],
               type = 'l', col = cols[i], lty = ltys[i])
        points(alphalist[alpha_pt], vals[i, alpha_pt],
               col = cols[i], pch = pchs[i])
    }
    if (legend){
        legend("topleft", methods,
               col = cols, lty = ltys, pch = pchs,
               seg.len = 3, cex = cex.legend, bty = "n")
    }
}

plot.corr <- function(corr.dfs, title,
                      ylim = c(0.7, 1),
                      cols = 1,
                      ltys = 1,
                      legend = NULL,
                      legend.position = "bottomright",
                      cex.legend = 1.3,
                      ...){
    par(mfrow = c(1, 1), ...)
    plot(0:1, 0:1, xlim = c(0.5, 0.01), ylim = ylim, type = 'n',
         xlab = expression(paste("Target FDR ", alpha, " (large --> small)")),
         ylab = 'Correlation', main = title, axes=FALSE)
    axis(side = 1, at = c(0.5, 0.4, 0.3, 0.2, 0.1, 0))
    axis(side = 2)
    if (class(corr.dfs)[1] == "list"){
        for (i in 1:length(corr.dfs)){
            corr.df <- corr.dfs[[i]]
            corr.df <- corr.df[corr.df$alpha <= 0.5, ]
            points(corr.df$alpha, corr.df$corr, type = 'l',
                   lwd = 2, col = cols[i], lty = ltys[i])
        }
        if (!is.null(legend)){
            legend(legend.position, legend,
                   col = cols, lty = ltys,
                   bty = "n", cex = cex.legend)
        }
    } else {
        corr.dfs <- corr.dfs[corr.dfs$alpha <= 0.5, ]
        points(corr.dfs$alpha, corr.dfs$corr, type = 'l',
               lwd = 2, col = cols, lty = ltys)
    }
}
