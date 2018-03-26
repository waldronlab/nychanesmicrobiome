#' Import NYC HANES-II microbiome data as a \code{phyloseq} object.
#' 
#' @export
loadQiimeData <- function(){
  
  otu_table <- system.file("extdata","otu_table_mc10_w_tax.biom", package="nychanesmicrobiome", mustWork = TRUE)
  tree_otu <- system.file("extdata","rep_set.tre", package="nychanesmicrobiome", mustWork = TRUE)
  rep_set <- system.file("extdata","rep_set.fna", package="nychanesmicrobiome", mustWork = TRUE)
  
  phylo <- phyloseq::import_biom(BIOMfilename = otu_table, 
                                 treefilename = phyloseq::read_tree(tree_otu), 
                                 refseqfilename = rep_set, 
                                 refseqFunction = phyloseq::parse_taxonomy_default)
  colnames(phyloseq::tax_table(phylo)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  #Remove controls
  phylo <- prune_samples(!(sample_names(phylo) %in% c("20151013CZTME1","20151020CZE3","20151020CZE4","20151020TME1","NC1","NC2")), phylo)
  #Remove replicates
  phylo <- prune_samples(!(sample_names(phylo) %in% c("NYDH0036","NYDH0051","NYDH0060","NYDH0152","NYDH0213","NYDH0487","NYDH0492","NYDH0522","NYDH0527","NYDH0545R","NYDH0649c", "NYDH0661", "NYDH0691", "NYDH0893","NYDH0931","NYDH0988","NYDH1042","NYDH1460","NYDH1353"))
                         , phylo)
  
  metadata <- sas7bdat::read.sas7bdat(system.file("extdata","public_v2_010518.sas7bdat", package="nychanesmicrobiome", mustWork = TRUE))
  link <- sas7bdat::read.sas7bdat(system.file("extdata","linked_for_cuny_031618.sas7bdat", package="nychanesmicrobiome", mustWork = TRUE))
  
  
  new_metadata <- dplyr::left_join(link, metadata, by='KEY')
  
  sample_names(phylo) <- stringr::str_replace_all(sample_names(phylo), 'c|R', '')
  
  
  f_sample_selection <- system.file("extdata", "smokingsampleselection.csv", package="nychanesmicrobiome", mustWork = TRUE)
  sample_selection <- read.csv(f_sample_selection)
  
  altsmokers <- sample_selection[sample_selection$final.altsmokers,'key']
  never <- sample_selection[sample_selection$final.nevers,'key']
  former <- sample_selection[sample_selection$final.formers,'key']
  secondhand <- sample_selection[sample_selection$final.secondhand,'key']
  cigarettes <- sample_selection[sample_selection$final.cigarettes,'key']
  
  sample_selection[sample_selection$key %in% altsmokers, 'smokingstatus'] <- 'alternativeonly'
  sample_selection[sample_selection$key %in% never, 'smokingstatus'] <- 'never'
  sample_selection[sample_selection$key %in% former, 'smokingstatus'] <- 'former'
  sample_selection[sample_selection$key %in% secondhand, 'smokingstatus'] <- 'secondhand'
  sample_selection[sample_selection$key %in% cigarettes, 'smokingstatus'] <- 'cigarette'
  
  new_metadata_smokingstatus <- dplyr::full_join(new_metadata, sample_selection, by=c('KEY' = 'key'))
  rownames(new_metadata_smokingstatus) <- new_metadata_smokingstatus$Burklab_ID
  sample_data(phylo) <- new_metadata_smokingstatus
  
  
  #Remove samples with less than 1000 reads
  phylo <- phyloseq::prune_samples(phyloseq::sample_sums(phylo)>1000, phylo)
  
  
  #Remove OTU classified as chloroplasts and mitochondria
  phylo <- phyloseq::subset_taxa(phylo, !Class %in% c("D_2__Chloroplast") & !Family %in% c("D_4__Mitochondria"))
  
  #Merge splitted genera
  splitted_genera <- sapply(data.frame(phyloseq::tax_table(phylo)[,"Genus"]), function(x)  grep(" [1-9]",x))
  phyloseq::tax_table(phylo)[splitted_genera,"Genus"] <- sapply(phyloseq::tax_table(phylo)[splitted_genera,"Genus"], function(x) gsub(" [1-9]","", x))
  drop_prefixes <- function(x) {
    x <- as.character(x)
    x <- strsplit(x,"__")
    x <- sapply(x, `[`, 2)
    x
  }
  phyloseq::tax_table(phylo)[,"Phylum"] <- drop_prefixes(phyloseq::tax_table(phylo)[,"Phylum"])
  phyloseq::tax_table(phylo)[,"Genus"] <- drop_prefixes(phyloseq::tax_table(phylo)[,"Genus"] )
  
  phylo
}

#' @export
loadQiimeITSData <- function(){
  
  otu_table <- "../input/qiime_data/ITS/otu_table_mc2_w_tax_w_metadata.txt"
  mapping_file <- "../input/qiime_data/ITS/mappingFile_ITS.tsv_corrected.txt"
  rep_set <- "../input/qiime_data/ITS/rep_set.fna"
  
  phyloseq::import_biom(BIOMfilename = otu_table,  
                        refseqfilename = rep_set, 
                        refseqFunction = phyloseq::parse_taxonomy_default)
}

