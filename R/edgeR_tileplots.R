
edger_summary_plot <- function(list_models, ttle="Significant genera by sociodemographic characteristics",
                               top_n_genera=NULL) {
  
  
  require(magrittr)
  require(dplyr)
  list_models <- nullify_nondataframes(list_models)
  list_models <- drop_nulls_from_list(list_models)
  
  #create a data frame with a column of all genera appearing in any of the models
  all_sig_OTUs <- data.frame(genus =
                               na.omit(unique(unlist(lapply(list_models, function(model)
                                 unlist(lapply(model, function(dataframe)
                                   as.data.frame(dataframe)[, "Genus"])))))))
  
  
  #function to count how many of the models (including all factor levels) contain a particular OTU
  count_of_variables_with_OTU <- function(otu_vect)  {
    sapply(otu_vect, function(otu)
      sum(unlist(lapply(list_models, function(model)
        any( unlist(lapply(model, function(df)
          otu %in% df[,"Genus"]))))  ))
    )
  }
  
  #create a column in the data frame with the count for each genus
  all_sig_OTUs$count <- count_of_variables_with_OTU(all_sig_OTUs$genus)
  
  #return only top genera if top_n_genera is a number
  if(!is.null(top_n_genera)) {
    stopifnot(is.numeric(top_n_genera))
    all_sig_OTUs <- all_sig_OTUs[order(all_sig_OTUs$count, decreasing = TRUE), ][1:top_n_genera,]
  }
  
  #sort models by number of significant genera
  n_genera <- sapply(names(list_models),
                     function(i) length(unique(unlist( lapply(list_models[[i]], function(j) j$Genus)))))
  n_genera <- n_genera[order(n_genera, decreasing = TRUE)]
  list_models <- list_models[ names(n_genera) ]
  
  #update list names to include number of significant genera (not doing is if
  #top_n_genera because numbers don't make sense)
  if(is.null(top_n_genera)) names(list_models) <- paste0(names(n_genera), " (", n_genera, ")")
  
  #generate a matrix indicating which OTU appeared in which model
  matrix_OTU_appears_in_model <-
    t(
      sapply( all_sig_OTUs$genus, function(otu)
        unlist(
          lapply(list_models, function(model)
            as.numeric(
              any(
                unlist(
                  lapply(model, function(df)
                    otu %in% df[,"Genus"]
                  )
                )
              )
            )
          )
        )
      )
    )
  
  #if these two procedures produced the same result,
  if( assertthat::are_equal(all_sig_OTUs$count, rowSums(matrix_OTU_appears_in_model))) {
    #add the rownames
    rownames(matrix_OTU_appears_in_model) <- all_sig_OTUs$genus
    
    #bind it to the data frame
    all_sig_OTUs <- cbind(all_sig_OTUs, as.data.frame(matrix_OTU_appears_in_model))
    all_sig_OTUs <- dplyr::arrange(all_sig_OTUs, desc(count)) # sort by count of appearances
    all_sig_OTUs <- dplyr::select(all_sig_OTUs, -count) # remove the count variable, not necessary
  }
  
  melted_appearances <- reshape2::melt(matrix_OTU_appears_in_model)
  melted_appearances[,1] <- factor(melted_appearances[,1],
                                  levels = all_sig_OTUs$genus)
  
  ggplot2::ggplot(melted_appearances, ggplot2::aes(x=Var1, y=Var2, fill=value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(name="Appeared") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust=1, hjust = 1),
                   legend.position = "none") +
    ggplot2::labs(y="Sociodemographic Variable", x="Genus", title=ttle)
  
}

edger_logfc_tileplot <- function(list_list_models) {
  
  dfs <- lapply(seq_along(list_list_models), function(i) edger_logfc_tileplot_df(list_list_models[[i]], oligo_name = names(list_list_models)[i]))
  df_tile <- do.call(rbind, dfs)
  
  ggplot2::ggplot(df_tile, ggplot2::aes(x=OTU, y=coef, fill=logFC)) +
    ggplot2::geom_tile() +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust=1, hjust = 1)) +
    ggplot2::scale_fill_gradient2(limits=c(-3,3.5),low = "seagreen4", mid = "white", high = "purple3") +
    ggplot2::facet_wrap(~oligo_name,ncol=1, dir="h", scales="free")
}

edger_logfc_tileplot_df <- function(list_models, oligo_name) {
  
  var_df <- function(edger_list) {
    list_edger_results <- lapply(seq_along(edger_list), function(i)
      data.frame(logFC = edger_list[[i]]@.Data[[1]]$logFC,
                 OTU = rownames(edger_list[[i]]@.Data[[1]]),
                 coef = edger_list[[i]]@.Data[[3]][1])
    )
    df <- do.call(rbind, list_edger_results)
    # df$OTU <- factor(df$OTU, labels = paste0("Oligotype ", levels(as.factor(df$OTU))))
    df <- aggregate(logFC~OTU+coef,df,function(x) x[which.max(abs(x))]) #when more than one OTU for a given genus was significant, select the max logFC to display
    df
  }
  
  
  for(i in seq_along(list_models)) list_models[[i]][sapply(list_models[[i]], nrow)==0] <- NULL
  list_models[sapply(list_models, length)==0] <- NULL
  
  df_tile <- lapply(list_models, var_df)
  for(i in seq_along(df_tile)) df_tile[[i]]$variable <- names(df_tile)[i]
  df_tile <- do.call(rbind, df_tile)
  df_tile$OTU <- factor(df_tile$OTU, labels = paste0("O", 1:nlevels(df_tile$OTU)))
  df_tile$oligo_name <- oligo_name
  
  df_tile
}


