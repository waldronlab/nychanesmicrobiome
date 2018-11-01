code_smoking_selection_variable <- function(nychanes_sas) {
  smq <- as.matrix(nychanes_sas[, grep("SMQ_13", colnames(nychanes_sas)) ])
  smoke_df <- data.frame(key= nychanes_sas$KEY,
                         hookah=apply(smq, 1, function(x) any(x==6)),
                         ecig=apply(smq, 1, function(x) any(x==7)),
                         cigar=apply(smq, 1, function(x) any(x==3)),
                         never=(nychanes_sas$SMQ_1==2 & nychanes_sas$SMQ_12==2),
                         cigarette=apply(smq, 1, function(x) any(x==1)) & nychanes_sas$SMQ_1==1,  #cigarettes in last 5 days & smoked 100 cigarettes in lifetime
                         otherprod=apply(smq[, grep("SMQ_13_2", colnames(smq))], 1, function(x) any(x %in% 2:8)),  #other nicotine products
                         stringsAsFactors = FALSE)
  smoke_df$former <- (nychanes_sas$SMQ_1==1 & nychanes_sas$SMQ_3==3 & nychanes_sas$SMQ_12==2)  #lifetime 100+ & currently not at all & no hookah/ecig/cigar
  smoke_df[is.na(smoke_df)] <- FALSE
  smoke_df$cotinine <- nychanes_sas$COTININE
  
  altsmokers <- as.character(smoke_df$key)[with(smoke_df, (hookah | ecig | cigar))] 
  nevers <- as.character(smoke_df$key)[!is.na(smoke_df$cotinine) & smoke_df$never & (smoke_df$cotinine < 0.05)] 
  
  formers <- as.character(smoke_df$key)[!is.na(smoke_df$cotinine) & smoke_df$former & (smoke_df$cotinine < 0.05)] 
  secondhand <- smoke_df[!is.na(smoke_df$cotinine) & (smoke_df$former | smoke_df$never) & !(smoke_df$key %in% altsmokers) & !smoke_df$otherprod & smoke_df$cotinine > 1 & smoke_df$cotinine < 14, ]
  secondhand <- as.character(secondhand$key) #all
  cigarettes <- smoke_df[!is.na(smoke_df$cotinine) & smoke_df$cigarette & !(smoke_df$key %in% altsmokers) & !smoke_df$otherprod, ]
  cigarettes <- cigarettes[rank(-cigarettes$cotinine) <= 90, ]
  cigarettes <- as.character(cigarettes$key)
  
  smoke_df$final.altsmokers <- smoke_df$key %in% altsmokers
  smoke_df$final.nevers <- smoke_df$key %in% nevers
  smoke_df$final.formers <- smoke_df$key %in% formers
  smoke_df$final.secondhand <- smoke_df$key %in% secondhand
  smoke_df$final.cigarettes <- smoke_df$key %in% cigarettes
  
  smoke_df$smokingstatus <- "alternativeonly"
  smoke_df$smokingstatus[smoke_df$never] <- "never"
  smoke_df$smokingstatus[smoke_df$former] <- "former"
  smoke_df$smokingstatus[smoke_df$final.secondhand] <- "secondhand"
  smoke_df$smokingstatus[smoke_df$cigarette & !smoke_df$final.altsmokers] <- "cigarette"
  
  final <- cbind(nychanes_full, smokingstatus = smoke_df$smokingstatus)
  
  final
  
}

#### ORIGINAL FILE ####
# library(sas7bdat)
# ##unzip ../input/NYC-HANES.zip to get this file:
# jesica=sas7bdat::read.sas7bdat("../input/microbiome.sas7bdat")
# jesica2=sas7bdat::read.sas7bdat("../input/NYCHANE2013-4_03212016.sas7bdat") 
# 
# ##Remove any subjects not consenting to future use:
# jesica <- jesica[!is.na(jesica$Micro_Test) & jesica$Micro_Test==1, ]
# jesica <- jesica[!is.na(jesica$key), ]
# stopifnot(identical(all(jesica$Micro_Test == 1), TRUE))
# 
# ## samples missing the saliva sample:
# saliva.missing <- c("N-22102", "N-20192", "N-22098")
# summary(jesica$SALIVA_ID %in% saliva.missing)
# jesica <- jesica[!jesica$SALIVA_ID %in% saliva.missing, ]
# 
# dim(jesica)  #1347

# levi <- data.frame(key= nychanes_full$KEY,
#  hookah=apply(smq, 1, function(x) any(x==6)),
#  ecig=apply(smq, 1, function(x) any(x==7)),
#  cigar=apply(smq, 1, function(x) any(x==3)),
#  never=(nychanes_full$SMQ_1==2 & nychanes_full$SMQ_12==2),
#  cigarette=apply(smq, 1, function(x) any(x==1)) & nychanes_full$SMQ_1==1,  #cigarettes in last 5 days & smoked 100 cigarettes in lifetime
#  otherprod=apply(smq[, grep("SMQ_13_2", colnames(smq))], 1, function(x) any(x %in% 2:8)),  #other nicotine products
#  stringsAsFactors = FALSE)
# levi$former <- (nychanes_full$SMQ_1==1 & nychanes_full$SMQ_3==3 & nychanes_full$SMQ_12==2)  #lifetime 100+ & currently not at all & no hookah/ecig/cigar
# levi[is.na(levi)] <- FALSE
# levi$cotinine <- nychanes_full$COTININE
#levi$never[levi$cigarette | levi$cigar | levi$ecig | levi$hookah] <- FALSE

#table(jesica$Hookah_pipe, levi$hookah)

# hist(levi$cotinine[levi$never==TRUE], xlim=c(0, .5), breaks="FD")
# 
# ## Random number generation
# set.seed(1)
# randnums <- matrix(sample(1:20000), ncol=10)
# samplefun <- function(obj, col, n){
#     ##obj = randnums, col is which column to use, output is similar to sample(1:n)
#     order(obj[1:n, col])
# }

# altsmokers <- as.character(levi$key)[with(levi, (hookah | ecig | cigar))] #all
# nevers <- as.character(levi$key)[!is.na(levi$cotinine) & levi$never & (levi$cotinine < 0.05)] #45
#ids <- samplefun(randnums, 1, length(nevers))
#nevers <- nevers[ids <= 45]

# formers <- as.character(levi$key)[!is.na(levi$cotinine) & levi$former & (levi$cotinine < 0.05)] #45
# #ids <- samplefun(randnums, 2, length(formers))
# #formers <- formers[ids <= 45]
# secondhand <- levi[!is.na(levi$cotinine) & (levi$former | levi$never) & !(levi$key %in% altsmokers) & !levi$otherprod & levi$cotinine > 1 & levi$cotinine < 14, ]
# secondhand <- as.character(secondhand$key) #all
# cigarettes <- levi[!is.na(levi$cotinine) & levi$cigarette & !(levi$key %in% altsmokers) & !levi$otherprod, ]
# cigarettes <- cigarettes[rank(-cigarettes$cotinine) <= 90, ]
# cigarettes <- as.character(cigarettes$key)
# 
# levi$final.altsmokers <- levi$key %in% altsmokers
# levi$final.nevers <- levi$key %in% nevers
# levi$final.formers <- levi$key %in% formers
# levi$final.secondhand <- levi$key %in% secondhand
# levi$final.cigarettes <- levi$key %in% cigarettes
# 
# levi$smokingstatus <- "alternativeonly"
# levi$smokingstatus[levi$never] <- "never"
# levi$smokingstatus[levi$former] <- "former"
# levi$smokingstatus[levi$final.secondhand] <- "secondhand"
# levi$smokingstatus[levi$cigarette] <- "cigarette"
# 
# final <- cbind(nychanes_full, levi)
# selected.log <- apply(levi[, grep("final", colnames(levi))], 1, any)
# summary(selected.log)
# #final <- final[selected.log, ]

# tmp <- final[, colnames(final) %in% colnames(levi)][, -1:-2]
# for (i in 1:ncol(tmp)) tmp[, i] <- as.numeric(tmp[, i])
# tmp <- tmp[order(tmp$final.altsmokers, tmp$final.cigarettes, tmp$final.secondhand, tmp$final.formers, tmp$final.nevers), ]
# tmp$cotinine <- scale(tmp$cotinine)
# tmp <- t(as.matrix(tmp))
# pdf("heatmap.pdf")
# library(gplots)
# heatmap.2(tmp, trace="none", scale="none", dendrogram="none", margins=c(2, 10), labCol=rep("", ncol(tmp)), Rowv=FALSE, Colv=FALSE)
# dev.off()
# system("open heatmap.pdf")
# 
# table(final$hookah, apply(final[, c("ecig", "cigar", "cigarette")], 1, any))
# table(final$ecig, apply(final[, c("hookah", "cigar", "cigarette")], 1, any))
# table(final$cigar, apply(final[, c("hookah", "ecig", "cigarette")], 1, any))
# table(apply(final[, c("cigar", "hookah", "ecig", "cigarette")], 1, sum))
# 
# table(final$hookah, final$ecig)
# table(final$hookah, final$cigar)
# 
# table(final$ecig, final$cigarette)
# table(final$ecig, final$cigar)
# 
# final$order <- samplefun(randnums, 10, nrow(final))
# write.csv(final, "../output/pilot_10-29-2015.csv")
# 
# jesica.match <- jesica2[match(final$key, jesica2$KEY), ]
# stopifnot(all.equal(as.character(jesica.match$KEY), as.character(final$key)))
# 
# summary(colnames(final) %in% colnames(jesica.match))
# summary(toupper(colnames(final)) %in% toupper(colnames(jesica.match)))
# missingcols <- colnames(final)[!toupper(colnames(final)) %in% toupper(colnames(jesica.match))]
# 
# jesica.match <- cbind(jesica.match, final[, missingcols])
# 
# write.csv(jesica.match, file="../output/pilot_05-11-2016.csv")
# 
# ##finalorig <- read.csv("../output/pilot_10-26-2015.csv", stringsAsFactors=FALSE)
# ##summary(finalorig$key %in% final$key)


