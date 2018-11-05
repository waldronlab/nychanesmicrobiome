#' Format sample data files as needed for sociodemographic and smoking analyses.
#' 
#' @param NYC_HANES An object of class \code{pseq}.
#' @param dat Data.frame.
#' 
#' @details These are convenience functions that accomplish all recoding and formatting of sample data for consistency across sociodemographics and smoking analyses. 
#' \code{annotateFactors} formats the \code{phyloseq} object of the n=297 pilot study, and \code{annotateFullDataset} formats the full NYC HANES data.frame. \code{formatMetadata} produces a dataset where factor levels represent very nice-looking labels for sociodemographics variables used in plotting and tables. 
#' 
#' @export
annotateFactors <- function(NYC_HANES){
	colnames(tax_table(NYC_HANES)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")


	sample_data(NYC_HANES)$SMOKER3CAT <- factor(sample_data(NYC_HANES)$SMOKER3CAT, labels  = c("Never smoker",
																							  "Current smoker",
																							  "Former smoker"))

	sample_data(NYC_HANES)$COTININE <- as.numeric(sample_data(NYC_HANES)$COTININE)

	sample_data(NYC_HANES)$COTININECAT <- factor(sample_data(NYC_HANES)$COTININECAT, labels = c("Active smoker",
																								"Non-smoker + secondhand smoke",
																								"Non-smoker no secondhand smoke",
																								"Skipped"))

	sample_data(NYC_HANES)$RACE <- factor(sample_data(NYC_HANES)$RACE, labels = c("Non-Hispanic White",
																				  "Non-Hispanic Black",
																				  "Hispanic",
																				  "Asian",
																				  "Other"))

	sample_data(NYC_HANES)$GENDER <- factor(sample_data(NYC_HANES)$GENDER, labels = c("Male","Female"))

	sample_data(NYC_HANES)$AGEGRP5C <- factor(sample_data(NYC_HANES)$AGEGRP5C, labels = c("20-29",
																						  "30-39",
																						  "40-49",
																						  "50-59",
																						  "60 AND OVER"))

	sample_data(NYC_HANES)$SR_ACTIVE <- factor(sample_data(NYC_HANES)$SR_ACTIVE, labels = c("Very active",
																										"Somewhat active",
																										"Not very active/not active at all"))

	sample_data(NYC_HANES)$EDU4CAT <- relevel(factor(sample_data(NYC_HANES)$EDU4CAT, labels = c( "Less than High school diploma",
															   "High school graduate/GED",
															   "Some College or associate's degree",
															   "College graduate or more")), "College graduate or more")

	sample_data(NYC_HANES)$OHQ_3 <- factor(sample_data(NYC_HANES)$OHQ_3, labels = c("Yes","No",NA))

	sample_data(NYC_HANES)$OHQ_5 <- as.numeric(sample_data(NYC_HANES)$OHQ_5)

	sample_data(NYC_HANES)$OHQ_5_3CAT <- cut(sample_data(NYC_HANES)$OHQ_5, breaks=c(-1,0,5,7),
	                                         labels = c("0","1-5","6-7"), include.lowest = FALSE)

	sample_data(NYC_HANES)$DBTS_NEW <- factor(sample_data(NYC_HANES)$DBTS_NEW, labels = c("Yes",
																						  "No",NA))

	sample_data(NYC_HANES)$smokingstatus <- factor(sample_data(NYC_HANES)$smokingstatus, levels = c("cigarette",
																									           "never",
																									           "former",
																									           "alternativeonly",
																									           "secondhand"),
																									           labels = c("Cigarette","Never smoker","Former smoker","Alternative smoker","Secondhand"))

	sample_data(NYC_HANES)$CURRENT_SMKER <- factor(sample_data(NYC_HANES)$CURRENT_SMKER, labels = c("Yes","No"))

	sample_data(NYC_HANES)$CIGARETTES <- factor(sample_data(NYC_HANES)$CIGARETTES, labels = c("Yes","No",NA))

	sample_data(NYC_HANES)$BMI <- as.double(sample_data(NYC_HANES)$BMI)

	sample_data(NYC_HANES)$SPAGE <- as.integer(sample_data(NYC_HANES)$SPAGE)

	sample_data(NYC_HANES)$MERCURYU <- as.double(sample_data(NYC_HANES)$MERCURYU )

	sample_data(NYC_HANES)$SMOKER <- factor(sample_data(NYC_HANES)$SMOKER, labels = c("Smoker",
	                                                                                  "Non smoker"))

	sample_data(NYC_HANES)$SMQ_13_1_1 <- factor(sample_data(NYC_HANES)$SMQ_13_1_1, labels = c("Cigarettes",
	                                                                                      "Cigars/Cigarillos",
	                                                                                      #"Chewing tobacco",
	                                                                                      #"Snuff",
	                                                                                      "Hookah Pipe",
	                                                                                      "E-Cigsarettes",
	                                                                                      #"Nicotine patches/Gum/Other",
	                                                                                      "NA"))

	sample_data(NYC_HANES)$SMQ_7 <- as.integer(sample_data(NYC_HANES)$SMQ_7)

	sample_data(NYC_HANES)$INC25KMOD <- factor(sample_data(NYC_HANES)$INC25KMOD, labels = c("Less Than $20,000",
	                                                                                  "$20,000-$49,999",
	                                                                                  "$50,000-$74,999",
	                                                                                  "$75,000-$99,999",
	                                                                                  "$100,000 or More",
	                                                                                  "NA"))

	sample_data(NYC_HANES)$POVGROUP6_0812CT <- factor(sample_data(NYC_HANES)$POVGROUP6_0812CT, labels = c("0 to 5%",
	                                                                                                      "5 to 10%",
	                                                                                                      "10 to 20%",
	                                                                                                      "20 to 30%",
	                                                                                                      "30 to 40%",
	                                                                                                      "40 to 100%"))
	sample_data(NYC_HANES)$SMQ_4UNIT <- factor(sample_data(NYC_HANES)$SMQ_4UNIT, labels = c(1,4,52,NA))

	sample_data(NYC_HANES)$US_BORN <- factor(sample_data(NYC_HANES)$US_BORN,
	                                         levels=c("1","2","NA"),
	                                         labels = c("US-Born, 50 States, DC, PR and Territories",
	                                                    "Other", NA))

	sample_data(NYC_HANES)$DMQ_7YEAR <- as.integer(sample_data(NYC_HANES)$DMQ_7YEAR)

	sample_data(NYC_HANES)$DMQ_2 <- factor(sample_data(NYC_HANES)$DMQ_2,
	                                       levels=1:6,
	                                       labels=c("Married","Widowed","Divorced",
	                                                "Separated","Never married",
	                                                "Living with partner"))

	sample_data(NYC_HANES)$AGEGRP3C <- factor(sample_data(NYC_HANES)$AGEGRP3C, levels=1:3,
	                                          labels=c("20-34", "35-64", "65 and over"))

	sample_data(NYC_HANES)$EDU3CAT <- relevel(factor(sample_data(NYC_HANES)$EDU3CAT, levels=1:3,
	                                         labels=c("High School Diploma or Less", "Some College or Associate's Degree",
	                                                  "College Graduate or More")),"College Graduate or More")


	sample_data(NYC_HANES)$INC10K <- as.numeric(as.character(sample_data(NYC_HANES)$INC10K))
	sample_data(NYC_HANES)$INC3C <- relevel(cut(sample_data(NYC_HANES)$INC10K,
	                                      breaks = stats::quantile(sample_data(NYC_HANES)$INC10K, probs=c(0, 1/3, 2/3, 1), na.rm=TRUE),
	                                      include.lowest = TRUE,
	                                      labels=c("Less Than $30,000", "$30,000 - $60,000", "$60,000 or more")), "$60,000 or more")


	sample_data(NYC_HANES)$A1C <- as.integer(sample_data(NYC_HANES)$A1C)

	sample_data(NYC_HANES)$GLUCOSE <- as.integer(sample_data(NYC_HANES)$GLUCOSE)

	sample_data(NYC_HANES)$DBQ_10 <- as.numeric(sample_data(NYC_HANES)$DBQ_10)
	sample_data(NYC_HANES)$DBQ_10 <- ifelse(sample_data(NYC_HANES)$DBQ_10UNIT=="1",
	                                        sample_data(NYC_HANES)$DBQ_10*7,
	                                        ifelse(sample_data(NYC_HANES)$DBQ_10UNIT=="3",
	                                               sample_data(NYC_HANES)$DBQ_10/4.33333,
	                                               sample_data(NYC_HANES)$DBQ_10))

	sample_data(NYC_HANES)$DBQ_10_3CAT <- cut(sample_data(NYC_HANES)$DBQ_10, breaks=c(-1,.999,5,999),
	                             labels = c("0-<1","1-5","6 or more"))

	NYC_HANES
}


#' @export
#' @rdname annotateFactors
annotateFullDataset <- function(dat) {
  # Similar to annotateFactors but for use with the complete
  # NYCHANES data as data.frame

  #converting NaN to NA
  dat[ is.na(dat) ] <- NA

  dat$RACE <- factor(dat$RACE, labels = c("Non-Hispanic White",
                                                      "Non-Hispanic Black",
                                                      "Hispanic",
                                                      "Asian",
                                                      "Other"))

  dat$GENDER <- factor(dat$GENDER, labels = c("Male","Female"))

  dat$AGEGRP5C <- factor(dat$AGEGRP5C, labels = c("20-29",
                                                              "30-39",
                                                              "40-49",
                                                              "50-59",
                                                              "60 AND OVER"))

  dat$AGEGRP4C <- factor(dat$AGEGRP4C, labels = c("20-34",
                                                              "35-49",
                                                              "50-64",
                                                              "65 AND OVER"))

  dat$AGEGRP3C <- factor(dat$AGEGRP3C, labels = c("20-34",
                                                              "35-64",
                                                              "65 AND OVER"))

  dat$SR_ACTIVE <- factor(dat$SR_ACTIVE, labels = c("Very active",
                                                                "Somewhat active",
                                                                "Not very active/not active at all"))

  dat$EDU3CAT <- factor(dat$EDU3CAT, labels = c("High school diploma or less",
                                                            "Some college or associate's degree",
                                                            "College graduate or more"))

  dat$EDU4CAT <- factor(dat$EDU4CAT, labels = c( "Less than High school diploma",
                                                             "High school graduate/GED",
                                                             "Some College or associate's degree",
                                                             "College graduate or more"))


  dat$OHQ_3 <- factor(dat$OHQ_3, labels = c("Yes",
                                                        "No"))
  dat$OHQ_5 <- as.numeric(sample_data(dat)$OHQ_5)

  dat$OHQ_5_3CAT <- cut(dat$OHQ_5, breaks=c(-1,0,5,7),
                                           labels = c("0","1-5","6-7"), include.lowest = FALSE)


  dat$DBTS_NEW <- factor(dat$DBTS_NEW, labels = c("Yes","No"))


  dat$CURRENT_SMKER <- factor(dat$CURRENT_SMKER, labels = c("Yes","No"))

  dat$BMI <- as.double(dat$BMI)

  dat$SPAGE <- as.integer(dat$SPAGE)

  dat$MERCURYU <- as.double(dat$MERCURYU )

  dat$SMOKER <- factor(dat$SMOKER, labels = c("Smoker",
                                                          "Non smoker"))

  dat$SMQ_13_1_1 <- factor(dat$SMQ_13_1_1, labels = c("Cigarettes",
                                                                  "Cigars/Cigarillos",
                                                                  "Chewing tobacco",
                                                                  "Snuff",
                                                                  "Hookah Pipe",
                                                                  "E-Cigsarettes",
                                                                  "Nicotine patches/Gum/Other"))

  dat$SMOKER3CAT <- factor(dat$SMOKER3CAT, labels  = c("Never smoker",
                                                                   "Current smoker",
                                                                   "Former smoker"))

  dat$SMQ_7 <- as.integer(dat$SMQ_7)

  dat$INC25KMOD <- factor(dat$INC25KMOD, labels = c("Less Than $20,000",
                                                                "$20,000-$49,999",
                                                                "$50,000-$74,999",
                                                                "$75,000-$99,999",
                                                                "$100,000 or More"))


  dat$INC10K <- as.numeric(as.character(dat$INC10K))
  dat$INC3C <- cut(dat$INC10K,
                                      breaks = stats::quantile(dat$INC10K, probs=c(0, 1/3, 2/3, 1), na.rm=TRUE),
                                      include.lowest = TRUE,
                                      labels=c("Less Than $30,000", "$30,000 - $60,000", "$60,000 or more"))


  dat$POVGROUP6_0812CT <- factor(dat$POVGROUP6_0812CT, labels = c("0 to 5%",
                                                                              "5 to 10%",
                                                                              "10 to 20%",
                                                                              "20 to 30%",
                                                                              "30 to 40%",
                                                                              "40 to 100%"))
  dat$US_BORN <- factor(dat$US_BORN,
                              levels=c("1","2"),
                              labels = c("US-Born, 50 States, DC, PR and Territories",
                                         "Other"))

  dat$DMQ_7YEAR <- as.integer(dat$DMQ_7YEAR)

  dat$DMQ_2 <- factor(dat$DMQ_2,
                            levels=1:6,
                            labels=c("Married","Widowed","Divorced",
                                     "Separated","Never married",
                                     "Living with partner"))

  dat$A1C <- as.integer(dat$A1C)

  dat$GLUCOSE <- as.integer(dat$GLUCOSE)

  dat$GLUCOSESI <- as.numeric(dat$GLUCOSESI)

  dat$DBQ_10 <- as.numeric(dat$DBQ_10)
  dat$DBQ_10UNIT <- as.character(dat$DBQ_10UNIT)
  dat$DBQ_10 <- ifelse(dat$DBQ_10 == 0, 0,
                       ifelse(dat$DBQ_10UNIT=="1",
                              dat$DBQ_10*7,
                              ifelse(dat$DBQ_10UNIT=="3",
                                     dat$DBQ_10/4.33333,
                                     dat$DBQ_10)))

  dat$DBQ_10_3CAT <- cut(dat$DBQ_10, breaks=c(-1,.999,5,999),
                                            labels = c("0-<1","1-5","6 or more"))

  dat
}

#' @export
#' @rdname annotateFactors
formatMetadata <- function(dat) {
  #returns a data frame of the sample_data with variable names formatted for printing
  #publication-ready tables



  rawdata <- if("phyloseq" %in% class(dat)) phyloseq::sample_data(dat) else dat

  newdata <- data.frame(
    `Age (yrs)` = rawdata$SPAGE,
    `Age group` = rawdata$AGEGRP5C,
    Gender = rawdata$GENDER,
    `Educational achievement` = rawdata$EDU4CAT,
    `Annual family income` = rawdata$INC3C,
    `Marital Status` = rawdata$DMQ_2,
    `Race/ethnicity` = rawdata$RACE,
    `Place of birth` = rawdata$US_BORN,
    `Gum disease (self-reported)` = rawdata$OHQ_3,
    `Mouthwash use (times per week)` = rawdata$OHQ_5_3CAT,
    `Sugar-sweetened beverages (per week)` = rawdata$DBQ_10_3CAT,
    `Smoking status` = rawdata$smokingstatus,
    check.names=FALSE
  )



  levels(newdata$`Age group`) <- c('20-29', '30-39', '40-49', '50-59', '60 and over')
  levels(newdata$`Mouthwash use (times per week)`) <- c("None", "1 to 5", "6 to 7")


  newdata$`Place of birth` <- factor(as.character(newdata$`Place of birth`),
                                     labels=c("US, PR and Territories","Other"))

  newdata

}


