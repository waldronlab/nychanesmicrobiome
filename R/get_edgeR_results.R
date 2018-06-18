#' User friendly interface to edgeR for microbiome data
#'
#' Returns an object of class \code{topTags} containing results filtered by FDR < alpha for one comparison. \code{get_edgeR_results_all_levels()} returns a list of \code{topTags} objects returned by \code{get_edgeR_results()} for all levels of the first variable in the formula.
#'
#' @param formla formula. specifies the design matrix used by \code{edgeR::glmFit}.
#' @param pseq object of class \code{phyloseq}
#' @param method character. Specifies the method of multiple testing correction to apply, either "BH" (Benjamini-Hochberg) or "IHW" (independent hypothesis weighting).
#' @param coef integer. Specifies which linear model coefficient to test (default 2).
#' @param alph numeric. Specifies what FDR level is considered alpha, and only keeps results with FDR less than this number.
#' @param filtering logical. Whether or not to apply pre-filtering.
#' @param countMinimum integer. If \code{filtering==TRUE}, the minimum count required for an OTU to be retained. Requires either \code{percentMinimumHaveCount} or \code{nMinimumHaveCount}.
#' @param percentMinimumHaveCount numeric. If \code{filtering==TRUE}, the minimum percentage of samples that must have \code{countMinimum} counts for an OTU to be retained.
#' @param nMinimumHaveCount numeric. If \code{filtering==TRUE}, the minimum number of samples that must have \code{countMinimum} counts for an OTU to be retained.
#' @param NA.RM logical. Whether or not to remove samples with NA values in the \code{sample_data()}.
#' @param vars character vector, containing names of each independent variable for which you would like a separate model.
#' @param varlabels character vector, containing variable labels for printing
#' @param adjusted_for character vector, containing names of the additional adjustment variables to include in each model
#' @param to.data.frame logical. Whether or not to return a combined data.frame of all model results, rather than the default list.
#' @param ... further arguments passed to \code{get_edgeR_results}
#'
#' @export
get_edgeR_results <- function(formla, pseq=NYC_HANES, method=c("BH","IHW"),
                              coef=2, alph=0.01, filtering=method[1] == "BH",
                              countMinimum = 8,
                              percentMinimumHaveCount = NULL,
                              nMinimumHaveCount = 3, NA.RM=TRUE) {
  ###################################################################
  # Get a GLM-like object filtered by FDR p values less than 'alph' #
  # for one comparison, specified by 'coef'                         #
  ###################################################################

  #coef=2 sets default to test first coefficient
  # If filtering = TRUE, user must provide either percentMinimumHaveCount or
  # nMinimumHaveCount, indicating how many samples must meet the countMinimum criteria
  # If method="IHW", no filtering is performed and and filtering arguments are ignored.


  #create model matrix
  dsg.mtrx <- model.matrix(formla, data=data.frame(sample_data(pseq)))
  #drop na's
  if(NA.RM) pseq <- prune_samples(samples = rownames(dsg.mtrx), x = pseq)

  #extract otu table to use as a count matrix
  otus <- as(otu_table(pseq), "matrix")

  #initialize DGEList object
  dge <- edgeR::DGEList(counts = otus)

  #filtering: keep OTU if minimum number of samples has minimum count
  if(filtering & !method[1] == "IHW") {
    if(!is.null(percentMinimumHaveCount)) {
      nMinimumHaveCount = ncol(otus) * percentMinimumHaveCount
    } else if(is.null(nMinimumHaveCount)) {
      stop("User must enter either nMinimumHaveCount or percentMinimumHaveCount.")
    }
    keep <- rowSums(dge$counts >= countMinimum ) >= nMinimumHaveCount
    dge <- dge[keep, , FALSE]
  }


  #set up the glm fit
  disp <- estimateDisp(dge, design = dsg.mtrx)
  fit <- glmFit(disp, dsg.mtrx)
  lrt <- glmLRT(fit, coef = coef)
  results <- topTags(lrt, n=nrow(lrt$table))

  if(method[1]=="IHW") {
    lrt$table$mean_abundance <- rowSums(otus) / ncol(otus)
    ihwResult <- IHW::ihw(PValue ~ mean_abundance, data=lrt$table, alpha=alph)
    results$table$ihw_pvalue <- ihwResult@df$weighted_pvalue
    #filter out results that didn't pass alpha
    results <- results[results$table$ihw_pvalue < alph, ]
  } else if(method[1]=="BH") {
    #filter out results that didn't pass alpha
    results <- results[results$table$FDR < alph, ]
  } else stop("Argument 'method' must be either 'BH' or 'IHW'.")

  # attach taxonomy names
  if(nrow(results$table) > 0) { #if any results pass the filter
    results$table <-
      cbind(results$table, as.matrix(tax_table(pseq)[rownames(results$table),]))
  }
  results
}

#' @export
#' @rdname get_edgeR_results
get_edger_results_all_levels <- function(formla, ...) {
  #the first variable in formula should be the one you want all levels of
  var1 <- strsplit(as.character(formla)[2], " \\+ ")[[1]][1]
  lapply(2:length(levels(metadata[[var1]])),
         function(i) get_edgeR_results(formla = formla,coef=i, ...))
}

#' @export
#' @rdname get_edgeR_results
edger_list_to_data.frame <- function(list_models) {
  #@param edger_list A list of objects each returned by
  #                          lapply(2:length(levels(x)), function(i) get_edgeR_results(~x,coef=i))
  var_df <- function(edger_list) {
    list_edger_results <- lapply(seq_along(edger_list), function(i)
      data.frame(logFC = edger_list[[i]]@.Data[[1]]$logFC,
                 OTU = rownames(edger_list[[i]]@.Data[[1]]),
                 coef = edger_list[[i]]@.Data[[3]][1],
                 FDR = edger_list[[i]]@.Data[[1]]$FDR)
    )
    df <- do.call(rbind, list_edger_results)
    # df$OTU <- factor(df$OTU, labels = paste0("Oligotype ", levels(as.factor(df$OTU))))
    df <- aggregate(logFC~OTU+coef+FDR,df,function(x) x[which.max(abs(x))]) #when more than one OTU for a given genus was significant, select the max logFC to display
    df
  }
  for(i in seq_along(list_models)) list_models[[i]][sapply(list_models[[i]], nrow)==0] <- NULL
  list_models[sapply(list_models, length)==0] <- NULL
  
  df_models <-  do.call(rbind, lapply(seq_along(list_models), function(i) data.frame(var_df(list_models[[i]]), variable=names(list_models)[i])) )
  df_models
}

#' @export
#' @rdname get_edgeR_results
get_all_edgeR_models <-  function(vars, varlabels=vars, adjusted_for, to.data.frame=TRUE, ...) {
  as.form <- function(vars) as.formula(paste0("~", paste(vars, collapse="+")))
  
  if(missing(adjusted_for)) {
    models = lapply(seq_along(vars),
                    function(i) get_edger_results_all_levels(as.form(vars[i]), ...))
  } else {
    models = lapply(seq_along(vars),
                    function(i) get_edger_results_all_levels(as.form(c(vars[i], adjusted_for)), ...))
  }
  names(models) <- varlabels
  if(to.data.frame) edger_list_to_data.frame(models) else models
}



edger_get_p <- function(formla, pseq=NYC_HANES, alph=1, ...)
{

  #######################################################
  # Get a vector of raw p-values for regression of all  #
  # OTUs against all levels of one metadata variable    #
  #######################################################
  vect <- sample_data(pseq)[,as.character(formla)[2]][[1]] # exctract the vector to get information about it

  n_coefs <- ifelse(class(vect) %in% c("character","factor"),
                    length(levels(factor(vect))) - 1, #the number of coeficients to test will be n levels - 1
                    1) # or 1 if it's a numeric etc.

  list_fits <- lapply(seq_len(n_coefs) + 1, # adding one because the coefficients start at 2 (skip the intercept)
                      function(coef)
                        get_edgeR_results(formla, pseq, coef=coef, alph=alph)  )

  p_vect <- unlist(sapply(seq_along(list_fits), function(i) list_fits[[i]]$table[,"PValue"]))
  p_vect
}
