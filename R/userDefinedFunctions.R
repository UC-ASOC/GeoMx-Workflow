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

getHistoryFile <- function() {
  c_args <- commandArgs()
  r_file <- c_args[grepl("\\.R$", c_args, ignore.case = TRUE)]
  r_file <- gsub("--file=", "", r_file)
  r_file <- normalizePath(r_file)
  return(r_file)
}
message("\nUser defined functions loaded.")