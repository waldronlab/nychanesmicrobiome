#' Create an auto-indented table based on tableone::CreateTableOne.
#' 
#' @param vars Character vector. Names of variables to summarize
#' @param data Data.frame. Where to search for \code{vars}.
#' @param nspaces Integer. Number of spaces to add to the beginning of indented rows.
#' @param includeNA Logical. Whether or not to treat NAs as a factor level.
#' @param nonnormal Character vector. Names of variables in \code{vars} to treat as non-normally distributed for summary and test functions.
#' @param ... Arguments passed to \code{CreateTableOne}
#' 
#' @export
prettytable_autoindent <- function(vars, data, nspaces=6, includeNA = TRUE, nonnormal=NULL, ...) { 
  #... args get passed to CreateTableOne
  
  table1 <- print(tableone::CreateTableOne(vars=vars, data=data, includeNA = includeNA, ...),
                  printToggle = FALSE, nonnormal=nonnormal)
  
  spaces = paste( rep("&nbsp;", nspaces), collapse="")
 
  varTypes <- data.frame(type = unlist(sapply(data[,vars], class)))
  varTypes$varname <- rownames(varTypes)
  
  
  varTypes$nrow <- sapply(seq_along(vars) , function(i) {
    if (varTypes$type[i] %in% c("factor","character")) {
      ifelse(length(levels(factor(data[, varTypes$varname[i]], 
                                  exclude=ifelse(includeNA, 9999, NA)))) == 2, 1,
             length(levels(factor(data[, varTypes$varname[i]], 
                                  exclude=ifelse(includeNA, 9999, NA)))) #using 9999 to represent excluding nothing
      )
    } else 1
  })
  varTypes$end <- cumsum(ifelse(varTypes$nrow == 1, 1, varTypes$nrow + 1))
  varTypes$start <- sapply(seq_along(varTypes$nrow), 
                           function(i) ifelse(i==1, 1, varTypes$end[i-1] + 1)) + 1
  varTypes$end <- varTypes$end + 1
  varTypes$start_indenting <- ifelse(varTypes$nrow == 1, NA, varTypes$start + 1)
  
  
  index <- na.omit(unlist(
    sapply(seq_len(length(vars)), 
           function(i) 
             if (is.na(varTypes$start_indenting[i])) NA 
           else seq(varTypes$start_indenting[i], varTypes$end[i]))
  ))
  
  
  rownames(table1)[index] <- paste0(spaces, rownames(table1)[index])
  table1
  
}