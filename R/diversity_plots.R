#' Alpha and beta diversity plots stratified by sample variables
#'
#' @param pseq An object of class \code{phyloseq}
#' @param vars Character vector. Names of variables in \code{sample_data} to stratify by.
#' @param metadata Data.frame. Defaults to \code{sample_data(pseq)[,vars]}, but this allows the use of say, a calculated variable without having to modify the \code{pseq}.
#' @param notch Logical. Whether or not to return notched boxplots.
#' @param dist_mtrx Matrix. Matrix of between-group distances, e.g. as returned by \code{phyloseq::distance()}.
#' @export
plot_alpha_by <- function(pseq, vars, metadata=sample_data(pseq)[,vars], notch=TRUE, adjust_formula, alpha_p) {
  require(dplyr)
  require(magrittr)
  require(reshape)
  require(phyloseq)
  
  df_alpha <- cbind(metadata, 
                    estimate_richness(pseq, measures = "Chao1")["Chao1"])
  
  if(!missing(adjust_formula)) {
    #this allows residuals instead of raw values
    mod <- lm(as.formula(paste0("Chao1 ~ ", adjust_formula)), data=df_alpha)
    df_alpha$Chao1 <-  resid(mod) + coef(mod)[1]
    df_alpha <- df_alpha[! sapply(names(df_alpha), function(i) length(grep(i,adjust_formula))>0)]
  }
  
  measure.vars <-names(df_alpha)[-which(names(df_alpha)=="Chao1")]
  
  
  
  if(missing(alpha_p)) alpha_p <-  sapply(measure.vars, function(i) 
    kruskal.test(Chao1 ~ eval(as.name(i)), data=df_alpha)$p.value) %>% round(3)
  
  df_alpha$id <- seq_len(nrow(df_alpha))
  
  
  melt_alpha <- na.omit(melt(df_alpha, measure.vars = measure.vars, id.vars="id"))
  melt_alpha %<>% left_join(dplyr::select(df_alpha, id, Chao1), by="id")
  #sort levels of each variable in melted form
  melt_alpha$value <- factor(melt_alpha$value, 
                             levels=sapply(measure.vars, function(VAR) levels(df_alpha[,VAR])) %>% unlist %>% unique
  )
  melt_alpha %>% 
    ggplot(aes(x = value, y = Chao1)) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(aes(fill=variable), alpha=0.2, notch=notch) +
    geom_jitter(aes(color=variable), width = 0.3, shape=19, size=1) +
    facet_grid(. ~ variable, scales="free_x", space="free_x") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust=1),
          axis.text = element_text(size=14),
          axis.ticks = element_blank(),
          axis.title = element_text(size=15),
          strip.text=element_text(size=15),
          panel.grid.minor = element_blank(), 
          
          legend.position="none") +
    labs(x=NULL, y="Chao1 Alpha Diversity")
}

#' @export
#' @rdname plot_alpha_by
plot_beta_by <- function(pseq, vars, dist_mtrx) {
  
  suppressPackageStartupMessages({
    require(magrittr)
    require(ggplot2)
  })
  
  distances_within <- function(dist_mtrx, pseq, fact) {
    metadata <- sample_data(pseq)
    fact <- na.omit(eval(as.name(fact), envir = metadata))
    group_ids <- lapply(levels(fact), function(i) rownames(metadata[fact == i, ] ))
    names(group_ids) <- levels(fact)
    within_distances <- lapply(
      levels(fact), 
      function(i) as.matrix(dist_mtrx)[group_ids[[i]], group_ids[[i]]])
    result <- unlist(within_distances)
    result[result>0]
  }
  distances_between <- function(dist_mtrx, pseq, fact) {
    metadata <- sample_data(pseq)
    fact <- na.omit(eval(as.name(fact), envir = metadata))
    dist_mtrx <- as.matrix(dist_mtrx)
    group_ids <- lapply(levels(fact), function(i) rownames(metadata[fact == i, ] ))
    names(group_ids) <- levels(fact)
    between_distances <- lapply(
      levels(fact), 
      function(i) dist_mtrx[group_ids[[i]], -which(colnames(dist_mtrx) %in% group_ids[[i]])])
    result <- unlist(between_distances)
    result[result>0]
  }
  
  within = do.call(rbind, 
                   lapply(vars, 
                          function(i) data.frame(dist=distances_within(dist_mtrx, pseq, i),
                                                 var=i,
                                                 loc="within")))
  between = do.call(rbind,
                    lapply(vars,
                           function(i) data.frame(dist=distances_between(dist_mtrx, pseq, i),
                                                  var=i,
                                                  loc="between") ))
  betadiv = rbind(within, between)
  
  betadiv %>% 
    ggplot(aes(x = var, y = dist, fill=loc)) +
    geom_boxplot() +
    facet_grid(. ~ var, scales="free_x", space="free_x") +
    theme_minimal() 
}