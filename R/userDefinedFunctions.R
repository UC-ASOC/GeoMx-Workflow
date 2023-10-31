##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Heewon Seo (Heewon.Seo@UCalgary.ca)
# Updated on Oct 02, 2023
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# User defined functions
pasteListValues <- function(x, comment = TRUE){
        desc <- c()
        listNames = names(x)
        for (i in seq_along(x) ) {
                if (comment) {
                        desc <- rbind(desc, paste0("# ", listNames[i], "=", x[[i]]))
                } else {
                        desc <- rbind(desc, paste0(listNames[i], "=", x[[i]]))
                }
        }
        return(desc)
}

histQC <- function(df = NULL,
                   annotation = NULL,
                   colBy = NULL,
                   rowBy = NULL,
                   threshold = NULL,
                   scaleTrans = NULL,
                   xLabel = NULL,
                   cols = NULL) {
        if (is.null(xLabel)) { xLabel <- annotation }
        gg <- ggplot(df, aes_string(x = paste0("unlist(`", annotation, "`)"), fill = colBy)) +
                geom_histogram(bins = 50) +
                scale_fill_manual(values = cols) +
                facet_grid(as.formula(paste(rowBy, "~", colBy))) +
                geom_vline(xintercept = threshold, lty = "dashed", color = "black") +
                theme_bw() +
                guides(fill = "none") +
                labs(
                        title = "",
                        x = xLabel,
                        y = "Count"
                )
        theme(text = element_text(size = 12))


        if (!is.null(scaleTrans)) {
                gg <- gg + scale_x_continuous(trans = scaleTrans)
        }
        return(gg)
}

boxQC <- function(df = NULL,
                  x = NULL,
                  y = NULL,
                  colBy = NULL,
                  scaleTrans = NULL,
                  xLabel = NULL,
                  yLabel = NULL,
                  cols = NULL) {
        if (is.null(xLabel)) { xLabel <- x }
        if (is.null(yLabel)) { yLabel <- y }
        dodge <- position_dodge(width = 0.5)
        violin <- ggplot(data = df, aes_string(x = paste0("unlist(`", x, "`)"), y = paste0("unlist(`", y, "`)"), fill = paste0("unlist(`", x, "`)"))) +
                geom_violin(position = dodge, size = 0) +
                geom_boxplot(width = 0.1, position = dodge, fill="white") +
                facet_grid(as.formula(paste("~", colBy))) +
                scale_fill_manual(values = cols) +
                labs(
                        title = "",
                        x = xLabel, 
                        y = yLabel
                ) +
                theme_bw() +
                theme(
                        axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                        legend.position = "none", 
                        text = element_text(size = 12)
                )
        if (!is.null(scaleTrans)) {
                violin <- violin + scale_y_continuous(trans = scaleTrans)
        }

        return(violin)
}

panel.reg <- function(x, y, col) abline(lm(y~x), col=col) 
panel.lm =  function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...)  {
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) panel.reg(x[ok], y[ok], col.smooth)
}
panel.cor <- function(x, y, prefix = "", ...) {
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
        pearson <- cor.test(x, y, method="pearson", use="complete.obs")
        spearman <- cor.test(x, y, method="spearman", use="complete.obs")

        pearson_coef <- format(round(pearson$estimate, 4), nsmall = 4)
        pearson_p <- formatC(pearson$p.value, format="e", digits=2)
        spearman_coef <- format(round(spearman$estimate, 4), nsmall = 4)
        spearman_p <- formatC(spearman$p.value, format="e", digits=2)
	
        text(0.5, 0.5, paste("Pearson= ", pearson_coef, "\nP-value=", pearson_p, "\n\nSpearman= ", spearman_coef, "\nP-value=", spearman_p, sep=""), cex = 1.2, font = 6)
}

getHistoryFile <- function() {
  c_args <- commandArgs()
  r_file <- c_args[grepl("\\.R$", c_args, ignore.case = TRUE)]
  r_file <- gsub("--file=", "", r_file)
  r_file <- normalizePath(r_file)
  return(r_file)
}
message("\nUser defined functions loaded.")