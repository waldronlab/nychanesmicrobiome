# Microbiome analysis of the NYC-HANES study

To reproduce this analysis, you need two packages from CRAN (www.cran.r-project.org). To install these from the R command line:
```
install.packages(c("BiocManager", "devtools"))
```

Then you can install this package and build its vignettes as follows:
```
BiocManager::install(“waldronlab/nychanesmicrobiome”)
```

On some machines, the above command may not by default install all the required dependencies. If this happens to you, you can install the `nychanesmicrobiome` package and its dependencies as follows:

```
BiocManager::install(c("waldronlab/nychanesmicrobiome", "phyloseq", "lsr", "sas7bdat", "kableExtra", 
          "ggplot2", "IHW", "dplyr", "epitools", "statmod", "DESeq2", "edgeR", "GSVA", 
          "EnrichmentBrowser", "magrittr", "dplyr", "vegan", "reshape2", "BiocStyle"))
```


After installing the `nychanesmicrobiome` package, load it by doing:
```
library(nychanesmicrobiome)
```

Then see the code and results of previous analyses by doing:
```
browseVignettes("nychanesmicrobiome")
```

With the `nychanesmicrobiome` package loaded, you should be able to reproduce any of the code found in these "vignette" analyses. The `nychanesmicrobiome` package provides a number of documented functions for analysing this dataset. Some functions are defined in the package to keep the analysis scripts cleaner and separate out components used by multiple analysts:
```
help(package="nychanesmicrobiome")
```

In particular, `loadQiimeData()` will import the NYC HANES microbiome data as a phyloseq object, including cleaning and merging the participant data with the microbiome data.
