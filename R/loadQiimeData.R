#' Import NYC HANES-II microbiome data as a \code{phyloseq} object.
#' 
#' @export
loadQiimeData <- function(metadata = sas7bdat::read.sas7bdat("http://nychanes.org/wp-content/uploads/sites/6/2018/01/public_v2_010518.sas7bdat")){
  
  otu_table <- system.file("extdata","otu_table_mc10_w_tax.biom", package="nychanesmicrobiome", mustWork = TRUE)
  tree_otu <- system.file("extdata","rep_set.tre", package="nychanesmicrobiome", mustWork = TRUE)
  rep_set <- system.file("extdata","rep_set.fna", package="nychanesmicrobiome", mustWork = TRUE)
  
  phylo <- import_biom(BIOMfilename = otu_table, 
                                 treefilename = read_tree(tree_otu), 
                                 refseqfilename = rep_set, 
                                 refseqFunction = parse_taxonomy_default)
  colnames(tax_table(phylo)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  #Remove controls
  phylo <- prune_samples(!(sample_names(phylo) %in% c("20151013CZTME1","20151020CZE3","20151020CZE4","20151020TME1","NC1","NC2")), phylo)
  #Remove replicates
  phylo <- prune_samples(!(sample_names(phylo) %in% c("NYDH0036","NYDH0051","NYDH0060","NYDH0152","NYDH0213","NYDH0487","NYDH0492","NYDH0522","NYDH0527","NYDH0545R","NYDH0649c", "NYDH0661", "NYDH0691", "NYDH0893","NYDH0931","NYDH0988","NYDH1042","NYDH1460","NYDH1353"))
                         , phylo)
  
  link <- read.delim(system.file("extdata","original_map.tsv", package="nychanesmicrobiome", mustWork = TRUE))
  
  #drop the smoking status variable that is coded by annotateFullDataset()
  metadata$smokingstatus <- NULL
  
  new_metadata <- dplyr::left_join(link, metadata, by='KEY')
  
  sample_names(phylo) <- stringr::str_replace_all(sample_names(phylo), 'c|R', '')
  
  
  f_sample_selection <- system.file("extdata", "smokingsampleselection.tsv", package="nychanesmicrobiome", mustWork = TRUE)
  sample_selection <- read.delim(f_sample_selection)
  
  sample_selection$smokingstatus <- factor(as.matrix(sample_selection[,-1]) %*% 1:5, 
                                            levels=1:5, 
                                            labels=c("alternativeonly","never","former","secondhand","cigarette"))
  
  new_metadata_smokingstatus <- dplyr::full_join(new_metadata, sample_selection, by=c('KEY' = 'key'))
  rownames(new_metadata_smokingstatus) <- new_metadata_smokingstatus$Burklab_ID
  sample_data(phylo) <- new_metadata_smokingstatus
  
  
  #Remove samples with less than 1000 reads
  phylo <- prune_samples(sample_sums(phylo)>1000, phylo)
  
  
  #Remove OTU classified as chloroplasts and mitochondria
  phylo <- subset_taxa(phylo, !Class %in% c("D_2__Chloroplast") & !Family %in% c("D_4__Mitochondria"))
  
  #Merge splitted genera
  splitted_genera <- sapply(data.frame(tax_table(phylo)[,"Genus"]), function(x)  grep(" [1-9]",x))
  tax_table(phylo)[splitted_genera,"Genus"] <- sapply(tax_table(phylo)[splitted_genera,"Genus"], function(x) gsub(" [1-9]","", x))
  drop_prefixes <- function(x) {
    x <- as.character(x)
    x <- strsplit(x,"__")
    x <- sapply(x, `[`, 2)
    x
  }
  tax_table(phylo)[,"Phylum"] <- drop_prefixes(tax_table(phylo)[,"Phylum"])
  tax_table(phylo)[,"Genus"] <- drop_prefixes(tax_table(phylo)[,"Genus"] )
  
  phylo
}

#' @export
loadQiimeITSData <- function(){
  
  otu_table <- "../input/qiime_data/ITS/otu_table_mc2_w_tax_w_metadata.txt"
  mapping_file <- "../input/qiime_data/ITS/mappingFile_ITS.tsv_corrected.txt"
  rep_set <- "../input/qiime_data/ITS/rep_set.fna"
  
  import_biom(BIOMfilename = otu_table,  
                        refseqfilename = rep_set, 
                        refseqFunction = parse_taxonomy_default)
}

