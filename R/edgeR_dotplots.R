#' Dot plots based on edgeR results for microbiome data
#'
#' Plots results of \code{get_edgeR_results_all_levels()} using the \code{ggplot2} engine.
#'
#' @param formla formula. specifies the design matrix used by \code{edgeR::glmFit}.
#' @param ttle character. Title of plot.
#' @param varname character. Label for independent variable, for plotting.
#' @param pseq object of class \code{phyloseq}
#' @param coef integer. Specifies which linear model coefficient to test (default 2).
#' @param printIfEmpty logical. If no results pass alpha, whether or not to plot an empty plot.
#' @param invisbl logical. Whether or not to forgoe plotting. Either way, invisibly returns the \code{topTags} object returned by \code{get_edgeR_results}.
#' @param color character. Variable  to color by, from either the \code{tax_table} or the \code{topTags$table}.
#' @param sortby character. Variable to sort the plot by, from either the \code{tax_table} or the \code{topTags$table}.
#' @param ... further arguments passed to \code{get_edgeR_results}.
#' @export
plot_edgeR <- function(formla, ttle=NULL, varname=NULL,
                       pseq=NYC_HANES, coef=2, printIfEmpty=FALSE,
                       invisbl=FALSE, color="FDR", sortby=NULL, ...) {
  #########################################################################
  # Get a dot plot of log-fold-change of all significant OTUs against one #
  # factor level specified by 'coef'. Also invisibly returns a data frame #
  # of edgeR results equivalent to get_edgeR_results()$table.             #
  #########################################################################

  # formla and ... are passed to get_edgeR_results
  suppressPackageStartupMessages({

    require(dplyr)
  })

  DGELRT = get_edgeR_results(formla, pseq, coef = coef, ...)

  if(!printIfEmpty & nrow(DGELRT$table) == 0) {
    return("No results")
  }

  #extract the main results table
  dat <- data.frame(DGELRT$table)

  #drop NA Genus factor levels
  dat <- dat[!is.na(dat$Genus), ]


  #drop NA

  #sort genus factor levels by log fold change if sortby is missing
  if(missing(sortby)) {
    ord <- order(dat$logFC, decreasing = TRUE)
  } else {
    ord <- order(dat[[sortby]], dat$logFC, decreasing=TRUE)
  }
  dat$Genus <- as.character(dat$Genus)
  dat$Genus <- factor(dat$Genus, levels = unique(dat$Genus[ord]))


  #create title
  if(is.null(ttle)) {
    metadata <- data.frame(sample_data(pseq))
    #extract the variable name from the formula
    colname <- strsplit(as.character(as.formula(formla)), "\\ ")[[2]][1]
    if(is.null(varname)) varname <- colname
    if(class(metadata[,colname]) %in% c("character","factor")) {
      #if it's a categorical variable
      ttle=paste0(varname, ": ", levels(metadata[,colname])[coef], #comparing the level specified by 'coef'
                  " vs. ", levels(metadata[,colname])[1])        #to the reference level
    } else {
      #if it's a quantitative variable
      ttle=paste0(varname, ": dy/dx")
    }
  }

  phyla_colors <- c('#00B2FFFF','#FF9900FF','#FF0000FF','#0DFF00FF','#00FF9FFF','#A600FFFF','#B9FF00FF','#FF00ACFF','#0006FFFF')

  #plot
  p <-
    ggplot2::ggplot(dat, ggplot2::aes(x=logFC, y=Genus, color=eval(as.name(color)))) +
    ggplot2::geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
    ggplot2::geom_point(size=3) +
    ggplot2::guides(color=ggplot2::guide_legend(title=color)) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -90, hjust = 0, vjust=0.5),
          axis.text.y = ggplot2::element_text(face = "italic"),
          axis.text = ggplot2::element_text(size=12),
          legend.text = ggplot2::element_text(size=12),
          legend.title = ggplot2::element_text(size=14),
          title = ggplot2::element_text(size=15)) +
    ggplot2::scale_y_discrete(position="right") +
    ggplot2::labs(x="Log2-fold change") +
    ggplot2::ggtitle(ttle)
  #add the colors set permanently to phyla (for multiple plots)
  if(color=="Phylum") p <- p + ggplot2::scale_color_manual(values=phyla_colors)
  if(!invisbl) print(p)
  invisible(dat)
}


#functions to drop nulls from list of models
nullify_nondataframes <- function(lst)
  lapply(lst, function(l) lapply(l, function(k) if(is.data.frame(k)) {k} else {NULL}) )
drop_nulls_from_list <- function(lst){
  lst <- lst[ !sapply( lst, is.null ) ]
  if( !is.data.frame(lst) ){
    lst <- lapply( lst, drop_nulls_from_list)
  }
  return(lst[sapply(lst, function(i) length(i) > 0)  ])
}







