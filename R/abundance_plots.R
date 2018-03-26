#' Area and bar plots of sample-level relative abundance at specified taxonomic levels
#'
#' @param pseq An object of class \code{phyloseq}
#' @param top_n Integer. Number of taxonomic groups to display, sorted by relative abundance
#' @param taxrank Character. Taxonomic rank to display. Should be a column name of the \code{taxa_table} in \code{pseq}.
#' @param prop_of Character. Either "top_n" or "total". Indicates whether heights of area/bars should represent the proportion out of the displayed taxonomic groups, or the out of the total.
#' @param orderby_OTU Character (optional). The OTU name (i.e. rowname in \code{otu_table}) whose relative abundance you would like to sort by. This is what creates the visual gradient from left to right of, say, Firmicutes.
#' @param PALETTE Character of same length \code{top_n} of palette for area or bar chart. 
#' @param type Character. Either area or bar, specifying which type of plot. Area displays all samples left to right, whereas bars give mean abundances.
#' @param brewerpal Character. Name of RColorBrewer palette to use.
#' @param by Character. Variable in \code{sample_data} to stratify plots by.
#' @param ttle Character. Plot title.
#' @export
plot_abundance <- function(pseq, top_n=8, taxrank="Phylum", prop_of="top_n", type="area",
                           orderby_OTU, PALETTE)
{
  suppressPackageStartupMessages({
    require(dplyr)
    require(magrittr)
    require(phyloseq)
    require(ggplot2)
  })
  #glom the pseq by specified taxrank
  glom = tax_glom(pseq, taxrank = taxrank) 
  
  #if applicable, calculate abundance as a proportion of total n
  if(prop_of=="total") glom <- transform_sample_counts(glom, function(i) i / sum(i))
  
  #if top_n is larger than the number of taxa at this rank, replace it with the number of taxa
  if(is.infinite(top_n)) {
    top_n <- nrow(tax_table(glom)) 
  } else {
    top_n = min(c(top_n, nrow(tax_table(glom))))
  }
  
  
  #subset the pseq object to only top taxa
  top_taxa = names(sort(taxa_sums(glom), TRUE)[1:top_n])
  top_subset <- prune_taxa(top_taxa, glom)
  
  #sort by taxa abundance
  taxa_to_sort <- data.frame(
    id=1:top_n, 
    taxa = as.character(tax_table(top_subset)[,taxrank]),
    otu = as.character(taxa_names(top_subset)),
    row.names = as.character(taxa_names(top_subset))
  )
  taxa_order <- top_subset %>% otu_table %>% rowSums %>% sort(T) %>% names
  taxa_to_sort <- taxa_to_sort[taxa_order, ]
  
  #if applicable, calculate abundance as a proportion of top n
  if(prop_of=="top_n") top_subset <- transform_sample_counts(top_subset, function(i) i / sum(i))
  
  #melt for use in plotting
  areadat <- psmelt(top_subset) %>% dplyr::select(OTU, Sample, Abundance) %>%  
    mutate(Sample = as.numeric(factor(Sample)),    # geom_area only seems to work with a numeric x
           OTU = factor(OTU,                         # sort by taxa
                        levels=taxa_to_sort$otu,   
                        labels=taxa_to_sort$taxa))
  
  #if missing, orderby_OTU is the most abundant taxa
  if(missing(orderby_OTU)) orderby_OTU <- taxa_to_sort$taxa[1] 
  
  # get order by particular OTU abundance
  OTU_order <- areadat %>%     
    filter(OTU==orderby_OTU) %>% 
    arrange(Abundance) %>% 
    dplyr::select(Sample)
  
  # apply OTU order to data.frame
  areadat %<>% mutate(Sample = as.numeric( 
    factor(Sample, levels=factor(OTU_order$Sample)))) 
  
  #get percentages for label
  mean_abundance <- apply(otu_table(top_subset)[taxa_order,], 1, function(i) round(100*mean(i), 1))
  sd_abundance <- apply(otu_table(top_subset)[taxa_order,], 1, function(i) round(100*sd(i), 1))
  levels(areadat$OTU) <- paste0(levels(areadat$OTU), " (", mean_abundance, " ± ", sd_abundance, "%)")
  
  #create plot
  if(type=="area") {
    plt <- ggplot(arrange(areadat, desc(OTU), Abundance),
                  aes(x=Sample, y=Abundance, fill=OTU))+
      geom_area() + theme_minimal() +
      labs(fill=taxrank)
    
    #apply palette if applicable
    if(!missing(PALETTE)) plt <- plt + scale_fill_manual(values=PALETTE)
    
  } else if(type=="bar") {
    plt <- ggplot(arrange(areadat, desc(OTU), Abundance),
                  aes(x=OTU, y=Abundance, fill=OTU)) +
      geom_col() + theme_minimal() +
      theme(axis.text.x = element_text(vjust=1, hjust=1, angle=50),
            legend.position = "none")
    
    #apply palette if applicable
    if(!missing(PALETTE)) plt <- plt + scale_fill_manual(values=PALETTE)
  }
  
  plt + ggtitle(paste(taxrank, "abundance"))
  
}

#' @export
#' @rdname plot_abundance
plot_abundance_by <- function(pseq, brewerpal="Set1", by, type=c("area", "bar"), 
                              ttle=paste0("Relative phylum abundances by ", by)){
  suppressPackageStartupMessages({
    library(RColorBrewer)
    library(dplyr)
    library(ggplot2)
  })
  
  set.seed(14)
  palet=sample(brewer.pal(9, brewerpal), size = 9, replace = FALSE)
  
  #extract by vector and its levels
  #by_vect <- eval(quote(by), sample_data(pseq))
  by_vect <- eval(parse(text=paste0("sample_data(pseq)$",by)))
  by_levels <- levels(factor(by_vect))
  
  #glom by phylum
  glom <- tax_glom(pseq, taxrank = "Phylum")
  
  #subset most abundant 8 phyla for each level of factor
  phyla8_subsets <- vector("list", length(by_levels))
  for(j in seq_along(by_levels)) { 
    samp <- prune_samples(by_vect==by_levels[[j]], glom)
    phyla8_subsets[[j]] <-
      prune_taxa(names(sort(taxa_sums(samp), TRUE)[1:8]), 
                 samp)
    rm(samp)
  }
  
  
  #sort by phylum abundance
  phyla_to_sort <- vector("list", length(by_levels))
  for(i in seq_along(by_levels)) {
    phyla_to_sort[[i]] <- data.frame(id=1:8, phyla = as.character(tax_table(phyla8_subsets[[i]])[,"Phylum"]), 
                                     otu = as.character(taxa_names(phyla8_subsets[[i]])))
    rownames(phyla_to_sort[[i]]) <- phyla_to_sort[[i]]$otu
    phylum_ranks <- phyla8_subsets[[i]] %>% otu_table %>% rowSums %>% sort(TRUE) %>% names
    phyla_to_sort[[i]] <- phyla_to_sort[[i]][phylum_ranks, ]
  }
  
  #calculate abundance as a proportion of top 8
  prop <- vector("list", length(by_levels))
  for(i in seq_along(by_levels)) {
    prop[[i]] <- transform_sample_counts(phyla8_subsets[[i]], function(j) j / sum(j))
  }
  
  
  #melt for use in plotting
  areadat <- vector("list", length(by_levels))
  for(i in seq_along(by_levels)) {
    areadat[[i]] <- psmelt(prop[[i]]) %>% dplyr::select(OTU, Sample, Abundance) %>%  
      mutate(Sample = as.numeric(factor(Sample)))    # geom_area only seems to work with a numeric x
    areadat[[i]]$OTU = factor(areadat[[i]]$OTU,                         # sort by phyla
                              levels=phyla_to_sort[[i]]$otu,   
                              labels=as.character( phyla_to_sort[[i]]$phyla ) ) 
    areadat[[i]]$level <- by_levels[i]
  }
  
  
  if(type[1]=="bar") {
    
    #summarize mean abundance for bars
    bardat <- do.call(rbind, areadat) %>%
      group_by(level, OTU) %>% 
      summarise(mean_abundance = mean(Abundance)) 
    bardat$level <- factor(bardat$level, levels=by_levels)
    
    bar_abundance <- ggplot(data=arrange(bardat, level), aes(x=level, y=mean_abundance, fill=OTU)) +
      geom_bar(position="fill", stat="identity") +
      labs(fill="Phylum", x="", y="Relative Abundance") +
      scale_fill_manual(values=palet) +
      theme(legend.position="right") +
      ggtitle(ttle)
    bar_abundance
    
  } else if (type[1]=="area") {
    
    firmicutes_orders <- vector("list", length(by_levels))
    for(i in seq_along(by_levels)) {
      firmicutes_orders[[i]] <- areadat[[i]] %>%     # get order by firmicutes abundance
        filter(OTU=="Firmicutes") %>% 
        arrange(Abundance) %>% 
        dplyr::select(Sample)
      
      #apply order to 'sample' variable in melted abundances
      areadat[[i]]$Sample = as.numeric(factor(areadat[[i]]$Sample, levels=factor(firmicutes_orders[[i]]$Sample)))
    }
    
    
    areadat_bind <- do.call(rbind, areadat)
    
    #area plot
    area_abundance <- ggplot(arrange(areadat_bind, desc(OTU), Abundance), aes(x=Sample, y=Abundance, fill=OTU))+
      geom_area()+
      labs(fill="Phylum")+
      scale_fill_manual(values=palet) +
      facet_grid(. ~ level, scales = "free_x")+
      theme(axis.text.x = element_blank(),
            legend.position="bottom") +
      ggtitle(ttle)
    area_abundance
  }
  
  
  
}

