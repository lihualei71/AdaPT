################################################################
## Plot functions for AdaPT objects.
################################################################

plot.mix.1d <- function(x, s, s.new, params, lfdr){
    ## Plot results during the process.
    par(mfrow = c(4, 1))
    xx <- sort(as.numeric(x[, 1]), index.return = TRUE)
    ind <- xx$ix
    xx <- xx$x
    plot(xx, s[ind], type = 'l',
         ylim = c(min(c(s, s.new)), max(c(s, s.new))),
         main = "s(x)")
    lines(xx, s.new[ind], col = "red")
    plot(xx, lfdr[ind], ylim = c(0, 1),
         main = "lfdr(x)")
    plot(xx, params$pix[ind], type = 'l', ylim = c(0, 1),
         main = "pi(x)")
    plot(xx, params$mux[ind], type = 'l',
         main = "mu(x)")

}

plot.mask.2d <- function(x, mask, main, cex,
                         col.bg = "#FFB6C1",
                         col.fg = "#800000",
                         ...){
    par(...)
    color <- ifelse(mask, col.fg, col.bg)
    plot(x[, 1], x[, 2], type = "n", xlab = "", ylab = "",
         main = main)
    points(x[, 1], x[, 2], col = color, cex = cex)
}

plot.lfdr.2d <- function(x, lfdr, main, cex,
                         col.fg = "#000080",
                         col.bg = "#ADD8E6",
                         ...){
    par(...)    
    rbPal <- colorRampPalette(c(col.bg, col.fg))
    color <- rbPal(5)[as.numeric(cut(lfdr, breaks=5))]
    plot(x[, 1], x[, 2], type = "n", xlab = "", ylab = "",
         main = main)
    points(x[, 1], x[, 2], col = color, cex = cex)
}
