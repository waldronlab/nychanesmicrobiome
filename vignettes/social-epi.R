## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  echo = FALSE,
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>",
  fig.wide=TRUE
  )

suppressPackageStartupMessages({
  library(phyloseq)
  library(tableone)
  library(mosaic)
  library(ape)
  library(vegan)
  library(fpc)
  library(cluster)
  library(reshape2)
  library(dplyr)
  library(magrittr)
  library(ggplot2)
  library(gtable)
  library(grid)
  library(knitr)
  library(htmlTable)
  library(kableExtra)
  library(gridExtra)
  library(sas7bdat)
})

## ----load_data-----------------------------------------------------------
#pilot microbiome data
suppressPackageStartupMessages({
  NYC_HANES_raw <- loadQiimeData()
  NYC_HANES <- annotateFactors(NYC_HANES_raw)
})
metadata <- data.frame(sample_data(NYC_HANES))

#full NYC-HANES2 dataset 
nychanes_sas <- read.sas7bdat("http://nychanes.org/wp-content/uploads/sites/6/2015/11/NYCHANES_2013_14_V2_09212016.sas7bdat")
nychanes_full <- annotateFullDataset(nychanes_sas)

## ----table1, results='asis'----------------------------------------------

variablesToLook <- c("SPAGE", "AGEGRP5C", "GENDER","EDU4CAT","INC3C","DMQ_2","RACE","US_BORN",
                     "DMQ_7YEAR","A1C","GLUCOSE", "AGEGRP3C", "EDU3CAT",
                     "OHQ_3", "OHQ_5_3CAT","DBQ_10_3CAT")



table1_data <- rbind(transform(formatMetadata(NYC_HANES), dataset="Oral Microbiome Subsample"),
                     transform(formatMetadata(nychanes_full), dataset="Full NYC HANES Sample"))
names(table1_data) <- c(names(formatMetadata(NYC_HANES)), "dataset")


table1 <- CreateTableOne(vars=names(table1_data)[-which(names(table1_data)=="dataset")], strata="dataset",data=table1_data, includeNA = TRUE, test=FALSE)

table1 <- print(table1, smd=TRUE, printToggle=FALSE)


rownames(table1) <- rownames(prettytable_autoindent(vars = names(table1_data)[-which(names(table1_data)=="dataset")],data = table1_data, 
                         includeNA = TRUE))


rownames(table1)[1] <- "Total"
rownames(table1) <- gsub("NA", "<i>Missing</i>", rownames(table1))
rownames(table1)[2] <- "Age in years -- median [range]"
table1[2,1:2 ] <- c(paste0(median(table1_data$`Age (yrs)`[table1_data$dataset=="Oral Microbiome Subsample"]), " [", min(table1_data$`Age (yrs)`[table1_data$dataset=="Oral Microbiome Subsample"]),
                      " to ", max(table1_data$`Age (yrs)`[table1_data$dataset=="Oral Microbiome Subsample"]), "]"),
                 paste0(median(table1_data$`Age (yrs)`[table1_data$dataset=="Full NYC HANES Sample"],na.rm=TRUE), " [", 
                        min(table1_data$`Age (yrs)`[table1_data$dataset=="Full NYC HANES Sample"],na.rm=TRUE),
                      " to ", max(table1_data$`Age (yrs)`[table1_data$dataset=="Full NYC HANES Sample"],na.rm=TRUE), "]"))


knitr::kable(table1, format="markdown", caption="Table 1. Sample demographics and oral health behavioral characteristics.")


## ----cramersV, fig.cap="Examining collinearity among sociodemographic variables. Data are absolute value of pairwise Cramer’s V correlation coefficient between sociodemographic factor levels. Data are from the full sample (n=1,527) of the New York City Health and Nutrition Examination Survey, 2013-2014. Abbreviations: cat=categories; US=United States.."----
cramer_vars <- c('AGEGRP5C','GENDER','EDU4CAT','INC3C','DMQ_2','RACE','US_BORN','OHQ_3','OHQ_5_3CAT','SMOKER4CAT','DBQ_10_3CAT')
cramer_varlabs <-c("Age (5 cat)","Gender", "Education (4 cat)","Family Income (3 cat)",
                "Marital Status","Race/Ethnicity","US- vs. Foreign-Born","Gum Disease","Mouthwash Use",
                "Smoking Status (4 cat)","Sugar-sweetened beverages (3 cat)")
df_cramer <- nychanes_full[,cramer_vars]
names(df_cramer) <- cramer_varlabs


cramersV_matrix <- function (data) 
{
    mtrx <- sapply(1:ncol(data), function(var1) sapply(1:ncol(data), 
        function(var2) lsr::cramersV(data[, var1], data[, var2])))
    diag(mtrx) <- 1
    rownames(mtrx) <- colnames(mtrx) <- colnames(data)
    mtrx
}

mtrx_cramersV <- abs(cramersV_matrix(data=df_cramer))

color_args <- list(name="Cramer's V",low = "blue", 
  high = "darkorange")


# Make a triangle
mtrx_cramersV[lower.tri(mtrx_cramersV)] <- NA
diag(mtrx_cramersV) <- NA
mtrx_cramersV <- as.data.frame(mtrx_cramersV)

# Convert to a data frame, and add labels
df_cramer <- as.data.frame(mtrx_cramersV)
df_cramer$var2 <- seq_along(rownames(df_cramer))

# Reshape to suit ggplot, remove NAs, and sort the labels
df_cramer <- na.omit(reshape2::melt(df_cramer, id.vars="var2"))
df_cramer$variable <- factor(df_cramer$variable, levels=rev(levels(df_cramer$variable)))

p1 <- ggplot(df_cramer, aes(x=variable, y=var2, fill=value)) +
    geom_tile() +
    do.call(scale_fill_gradient, color_args) +
    scale_y_continuous(breaks=1:11, labels=rev(levels(df_cramer$variable))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust = 1),
          legend.position = "none") +
    labs(x=NULL, y=NULL)

#create boxplot
mtrx_cramersV$var2 <- rownames(mtrx_cramersV)
data_boxplot <-  reshape2::melt(mtrx_cramersV, "var2") %>% na.omit %>%
  mutate(grp=NA, alph=-1*value + 1,
         Name = paste0(var2, "*", variable))
p2 <- ggplot(data_boxplot, aes(y=value, x=grp, col=value)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(size=2) +
  do.call(scale_color_gradient, color_args) +
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(0, 0.1, 4, 0.5), "cm")) 

gt1 <- ggplot_gtable(ggplot_build(p1))
gt2 <- ggplot_gtable(ggplot_build(p2))
gt <- gtable(widths = unit(c(5.5, 2.5), "null"), height = unit(8, "null"))


# Instert gt1, gt2 and gt3 into the new gtable
gt <- gtable_add_grob(gt, gt1, 1, 1)
gt <- gtable_add_grob(gt, gt2, 1, 2)
grid.newpage()
grid.draw(gt)

grid.rect(x = 0.5, y = 0.5, height = 0.995, width = 0.995, default.units = "npc", 
    gp = gpar(col="transparent", fill = NA, lwd = 1))


## ----phylum_abundance, fig.cap="Genus- and phylum-level relative abundances. Data are percent of overall communities within samples, summarized as mean ± standard deviation of percent across samples. Data are from the oral microbiome subsample (n=296) of the New York City Health and Nutrition Examination Survey, 2013-2014."----
get_phyla_colors <- function(v,s) {
  set.seed(48)
  rev(sample(rainbow(n = 8,  start = 0.1, 
                         end = 1, s = s, v=v),
                 size=8, replace=FALSE))[c(2,1,3:8)]
  
}
phyla_colors <- c(rbind(get_phyla_colors(v=.7,s=.5)[c(2,4,6,8)], 
                        get_phyla_colors(v=.9,s=.7)[c(1,3,5,7)]))


g <- plot_abundance(NYC_HANES, PALETTE = phyla_colors, top_n = 8, prop = "total") +
  labs(x = "└────────────── 296 Samples ───────────────┘", y = "", title = "") + 
  scale_x_continuous(labels = NULL) + 
  scale_y_continuous(labels = seq(0,100,25))

p <- plot_abundance(NYC_HANES, top_n = 8, PALETTE =  phyla_colors, taxrank="Genus", 
               prop = "total", type="area") 
p <- p + labs(x="", y="", title="") + scale_x_continuous(labels=NULL) + scale_y_continuous(labels=seq(0,100,25))

grid.arrange(p + theme(plot.margin=unit(c(1,0.95,-0.7,0.1), "cm"),
                      panel.grid = element_blank()), 
             g + theme(plot.margin=unit(c(-0.7,0.7,1,-0.1), "cm"),
                      panel.grid = element_blank()), nrow=2, left="Abundance (%)")

## ----alpha_diversity, fig.cap="Alpha diversity by Sociodemographic Characteristics. Chao1 alpha diversity of 16S rRNA oral microbiome samples. Measures were compared using a null hypothesis of no difference between groups (Kruskal-Wallis test, p > 0.1 for all tests). Data are from the oral microbiome subsample (n=296) of the New York City Health and Nutrition Examination Survey, 2013-2014. Abbreviations: GED=General equivalency diploma; PR=Puerto Rico; US=United States.", fig.height=7.2, fig.width=12----

df_alpha <- cbind(table1_data[table1_data$dataset=="Oral Microbiome Subsample",
                              -which(names(table1_data) %in% c("Years in United States", "Age (yrs)",
                                                               "Gum disease (self-reported)",     
                                                               "Mouthwash use (times per week)",  
                                                               "Smoking status",
                                                               "Sugar-sweetened beverages (per week)",
                                                               "dataset"))  ])

names(df_alpha)[c(3,4,7)] <- c("Educational \nachievement", "Family \nincome","Place of \nbirth")

#calculate p-values
df_alpha <- cbind(df_alpha, estimate_richness(NYC_HANES, measures = "Chao1")["Chao1"])

 alpha_p <- "\n(P=" %>% paste0(sapply(names(df_alpha), function(i) 
    kruskal.test(Chao1 ~ eval(as.name(i)), data=df_alpha)$p.value) %>% formatC(digits=2, format = "f")   ) %>%
   paste0(")")

 names(df_alpha)[-8] %<>% paste0(alpha_p)


plot_alpha_by(NYC_HANES, metadata=df_alpha[-8], notch=FALSE)


## ----Import oligotyping output-------------------------------------------

f <- c("Neisseria","Prevotella","Streptococcus")

phylo.otut <- list()

for (x in f) {
  otumap <- read.table(paste0("../inst/extdata/", x,"_MATRIX-COUNT.txt"),  header = T,row.names = 1) %>% t
  otumap <- otumap[,rownames(sample_data(NYC_HANES))]
  taxmat <- matrix(sample(letters, 7*nrow(otumap), replace = TRUE), nrow = nrow(otumap), ncol = 7)
  rownames(taxmat) <- rownames(otumap)
  colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  phylo.otut[x] <- phyloseq(otu_table(otumap, taxa_are_rows = TRUE), 
                            tax_table(taxmat),
                            sample_data(NYC_HANES)) 
}


## ----tileplot_combine, fig.width=10, fig.height=8, fig.cap="Differential abundance by sociodemographic characteristics. OTUs (A) and oligotypes (B) meeting unadjusted FDR <0.01 in negative binomial log-linear GLMs using edgeR. Filled tiles in (A) indicate the genus had at least one OTU differentially abundant by at least one coefficient contrast within the sociodemographic factor. Most commonly differential genera in (A) included Prevotella (n=8) and Lactobacillus (n=7)."----

grid.arrange(plot_logfc, plot_oligo, heights=c(4,3.3))


## ----edger_crude_vs_adjusted_calc, cache=TRUE----------------------------


vars <- c("GENDER","AGEGRP3C","EDU3CAT", "INC3C", "DMQ_2", "RACE", "US_BORN" )
varlabels <- c("Gender", "Age", "Education", "Family income", "Marital status", 
               "Race/ethnicity", "US- vs. foreign-born")

df_crude <- df_tile
df_crude_nofilter <- get_all_edgeR_models(vars=vars, varlabels=varlabels, alph=1)
df_ohq <- get_all_edgeR_models(vars=vars, varlabels = varlabels, adjusted_for="OHQ_5_3CAT", alph=1) %>% dplyr::rename(logFC_ohq = logFC)

df_smoke <- get_all_edgeR_models(vars=vars, varlabels = varlabels, adjusted_for="smokingstatus", alph=1) %>% dplyr::rename(logFC_smoke = logFC)

df_age_gender <- get_all_edgeR_models(vars=vars, varlabels = varlabels, adjusted_for=c("AGEGRP3C","GENDER"), alph=1) %>% dplyr::rename(logFC_age_gender = logFC)

df_sugar <- get_all_edgeR_models(vars=vars, varlabels = varlabels, adjusted_for=c("DBQ_10_3CAT"), alph=1) %>% dplyr::rename(logFC_sugar = logFC)


## ----edger_crude_vs_adjusted_print, fig.width=7, fig.height=4.5, fig.cap="Distribution of absolute values of log-fold change (logFC). Values are from crude and adjusted edgeR models for all OTUs meeting FDR<0.01 in the crude model, for each sociodemographic variable."----

df_crude[,c(1,2,4)] %<>% lapply(as.character)
df_crude$variable <-  gsub("[0-9]| \\(|\\)","",df_crude$variable)
df_crude_nofilter[,c(1,2,4)] %<>% lapply(as.character)
df_ohq[,c(1,2,4)] %<>% lapply(as.character)
df_smoke[,c(1,2,4)] %<>% lapply(as.character)
df_age_gender[,c(1,2,4)] %<>% lapply(as.character)
df_sugar[,c(1,2,4)] %<>% lapply(as.character)

create_crude_vs_adjusted_boxplot_df <- function(df_crude) {
  df_crude %<>% inner_join(df_ohq, by=c("OTU","coef","variable")) %>%
    inner_join(df_smoke, by=c("OTU","coef","variable")) %>%
    inner_join(df_age_gender, by=c("OTU","coef","variable"))%>%
    inner_join(df_sugar, by=c("OTU","coef","variable"))
  
  df_crude[,grep("logFC", names(df_crude))] %<>% lapply(function(i) abs(as.numeric(i)))
  df_crude$smoke_diff <- (df_crude$logFC_smoke - df_crude$logFC) / df_crude$logFC
  df_crude$ohq_diff <- (df_crude$logFC_ohq - df_crude$logFC) / df_crude$logFC
  df_crude$agegender_diff <- (df_crude$logFC_age_gender - df_crude$logFC) / df_crude$logFC
  df_crude$sugar_diff <- (df_crude$logFC_sugar - df_crude$logFC) / df_crude$logFC
  df_crude
}

df_fdr <- create_crude_vs_adjusted_boxplot_df(df_crude)
df_nofdr <- create_crude_vs_adjusted_boxplot_df(df_crude_nofilter)

#for boxplot comparing absolute values (FDR)
df_fdr_melt <- melt(df_fdr, id.vars = "variable",measure.vars = c("logFC","logFC_ohq","logFC_smoke","logFC_age_gender","logFC_sugar"),
                    variable.name = "logFC")
names(df_fdr_melt) <- c("variable", "logFC","value" )

dodge_amount = 1

figure5 <- ggplot(df_fdr_melt, aes(y=value, x=variable, fill=logFC)) +
  stat_boxplot(geom="errorbar", position=position_dodge(dodge_amount)) +
  geom_boxplot(outlier.shape=NA, alpha=.3, position=position_dodge(dodge_amount)) +
  geom_point(aes(color=logFC), size=.5, position=position_jitterdodge(jitter.width=.3, dodge.width = dodge_amount)) +
  scale_fill_discrete(labels=c("Crude","Adjusted for Mouthwash Use",
                               "Adjusted for Smoking","Adjusted for Age and Gender","Adjusted for Sugar-sweetened Beverage Consumption"), name="") +
  scale_color_discrete(guide="none")+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle=40,hjust=1,vjust=1), legend.position = "top", panel.spacing = unit(.6, "lines"),
        strip.text = element_blank()) +
  labs(y="Absolute Value of logFC", x=NULL) +
  facet_grid(.~variable, scales="free_x", space="free_x")

figure5

## ----results='asis'------------------------------------------------------


df_table <- df_crude %>% mutate(direction=ifelse(as.numeric(logFC)>0, "Greater abundance in:","Lower abundance in:"))

commalist <- function(x) {
  sentence = x[1]
  if(length(x) > 1) {
    for(i in 2:(length(x))) {
      sentence <- paste0(sentence, "<br> ", x[i])
    }
  }
  sentence
}


selected_genera = c("Streptococcus", "Lactobacillus", "Lactococcus", "Prevotella", "Porphyromonas", "Fusobacterium")

reference_groups = sapply(vars, function(i) levels(metadata[[i]])[1])
names(reference_groups) <- varlabels




df_table %>% filter(Genus %in% selected_genera) %>% group_by(Genus, direction) %>% summarize(lst = commalist(unique(as.character(coef_lab))) ) %>% reshape2::dcast(Genus ~ direction) %>% kable(format="markdown", caption = "Table 2. Differential abundance findings for OTUs selected based on clinical relevance. Greater or lower abundance indicates false discovery rate (FDR) <0.01.") %>%
  kable_styling(bootstrap_options = "condensed")


tfooter = paste0("Reference groups for sociodemographic variables are as follows: ",
                 paste(paste0(names(reference_groups), ": ", reference_groups), collapse=", "))
tfooter


## ----permanova_one_way, cache=TRUE---------------------------------------
dist_perm <- phyloseq::distance(NYC_HANES, method="wunifrac")

permanova_var_names <- c("AGEGRP5C","GENDER","EDU4CAT","INC3C","DMQ_2","RACE","US_BORN","smokingstatus")

df_permanova <- metadata[, permanova_var_names]
#creating an NA factor level for INC3C because permanova doesn't allow NAs
df_permanova$INC3C %<>% as.character
df_permanova$INC3C[is.na(df_permanova$INC3C)] <- "NA"


list_permanovas <- list(
  `Age (5 cat)` = adonis2(dist_perm ~ AGEGRP5C, data=df_permanova),
  Gender = adonis2(dist_perm ~ GENDER, data=df_permanova),
  `Education (4 cat)` = adonis2(dist_perm ~ EDU4CAT, data=df_permanova),
  `Income (tertiles)` = adonis2(dist_perm ~ INC3C, data=df_permanova),
  `Marital Status` = adonis2(dist_perm ~ DMQ_2, data=df_permanova),
  `Race/ethnicity` = adonis2(dist_perm ~ RACE, data=df_permanova),
  `US-born` = adonis2(dist_perm ~ US_BORN, data=df_permanova),
  `Smoking status` = adonis2(dist_perm ~ smokingstatus, data=df_permanova)
)

SumOfSqs <- unlist(lapply(list_permanovas, function(i) i[["SumOfSqs"]][1]))
r_squared <- sapply(list_permanovas, function(i) 1 - (i$SumOfSqs[2]/sum(i$SumOfSqs)))
perm_results <- data.frame(
  SumOfSqs,
  F = unlist(lapply(list_permanovas, function(i) i[["F"]][1])),
  p_value = unlist(lapply(list_permanovas, function(i) i[["Pr(>F)"]][1]))
)

perm_results <- perm_results[names(sort(SumOfSqs, decreasing = TRUE)),]

## ----permanova_oneway_print----------------------------------------------
knitr::kable( perm_results)


## ----compare_beta_diversity, cache=TRUE, fig.width=5, fig.height=5, fig.cap="Within and between group beta diversity estimate distributions."----

pb <- plot_beta_by(pseq=NYC_HANES, vars = permanova_var_names, dist_mtrx = dist_perm)

pb <- pb+ theme(legend.title = element_blank(),
                strip.text.x = element_blank(),
                axis.text.x = element_text(angle=60,hjust=1,vjust=1),
                legend.position = "bottom") +
  scale_x_discrete(name=NULL, breaks=permanova_var_names,
                   labels=c("Age (5 cat)", "Gender","Education (4 cat)","Income (3 cat)",
                            "Marital Status","Race/ethnicity","Nativity","Smoking Status")) +
  scale_fill_brewer(palette = "Set2", labels=c("Within group","Between groups")) +
  scale_color_brewer(palette = "Set2", guide="none") +
  labs(y="Weighted UniFrac Distance")


pb

## ------------------------------------------------------------------------
sessionInfo()

