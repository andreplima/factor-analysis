#-----------------------------------------------------------------------------------------------------------
# This code performs confirmatory factor analysis on models induced from applying EFA to the
# Serendipity'2018 dataset, which has been generously made available by the authors:
# Kotkov, D., Konstan, J. A., Zhao, Q., & Veijalainen, J. (2018, April). Investigating serendipity in
#   recommender systems based on real  user feedback. In Proceedings of the 33rd Annual ACM Symposium on
#   Applied Computing (pp. 1341-1350). ACM.
# The complete dataset is available at: https://grouplens.org/datasets/serendipity-2018/
#
# We would also like to acknowledge the following authors for the R code for confirmaroty factor analysis,
# from which we benefited greatly:
# Yves Rosseel (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical
#   Software, 48(2), 1-36. URL http://www.jstatsoft.org/v48/i02/
#-----------------------------------------------------------------------------------------------------------

# removes the current environment variables
rm(list = ls())

# sets up some environment configuration options
ECO_CUT_LEVEL  = 0.3
ECO_CUT_2CHK   = 0.1
ECO_PRECISION  = 3
ECO_OPTIMETHOD = 'minres'

# sets the working directory
setwd("C:/Users/andre/OneDrive/Documentos/gitrepos/factor-analysis")
imageDirectory<-"C:/Users/andre/OneDrive/Documentos/gitrepos/factor-analysis/datasets/serendipity2018"
#filename = "answers-all-all-median.dat"
#filename = "answers-1399-all-median.dat"
filename = "answers-606-all-median.dat"


# installs necessary packages
#install.packages("lavaan")
#install.packages("MVN")

# loads libraries
library(svDialogs)
library(psych)
library(lavaan)
library(semPlot)
#library(MVN)
library(likert)
library(gsubfn)
library(dplyr)
library(tidyr)
library(knitr)
options(knitr.kable.NA = '') # this will hide missing values in the kable table


#---------------------------------------------------------------------------------------------------------
# General definitions
#---------------------------------------------------------------------------------------------------------

typeSer2018Data <- function(odata) {
  data <- odata
  
  #data$s3 <- NULL # removes the item that is becomes isolated 
  #idxs = c("s1", "s2", "s4", "s5", "s6", "s7", "s8")
  idxs = paste("s", 1:8, sep="")
  
  data[idxs] <- lapply(data[idxs], factor, levels = 1:5)
  data["rating"] <- lapply(data["rating"],  factor, levels = seq(.5,5,by=.5))
  datalk <- likert(data[idxs])
  return(list(odata, data, datalk))
}

loadSer2018Data <- function(datasetDirectory, filename){
  data <- read.delim(file.path(datasetDirectory, filename), header = TRUE)
  data$q <- NULL  # removes the item that is related to recency of event
  return(typeSer2018Data(data))
}

induceModel <- function(modelspec, trainingData) {
  postulated <- cfa(modelspec, data=trainingData, std.lv=TRUE, ordered=serVars)
  warnings()
  semPaths(postulated, what="std", edge.label.cex = 1, curvePivot = TRUE, whatLabels = "par", layout = "tree2", title=TRUE, XKCD=FALSE)
  return(postulated)
}

assessModel <- function(postulated) {
  
  # Applies some quality tests to the model that was fitted to the data
  # Some reference values were obtained from this document: https://www.cscu.cornell.edu/news/Handouts/SEM_fit.pdf
  v   <- round(fitMeasures(postulated, "cfi"), 3)
  MQ2 <- v > 0.9
  cat("\n")
  cat("   Comparative Fit Index (CFI) measures whether the model fits the data better than the baseline model.\n")
  cat("   .. CFI:", if (MQ2) "Passed" else "Failed", "   ... (value =", v, "> 0.9)\n")
  
  v   <- round(fitMeasures(postulated, "tli"), 3)
  MQ1 <- v > 0.95
  cat("\n")
  cat("   Tucker Lewis Index (TLI) is similar to CFI, but is non-normed and, thus, sensitive to sample size.\n")
  cat("   .. TLI:", if (MQ1) "Passed" else "Failed", "   ... (value =", v, "> 0.95)\n")
  
  v1  <- round(fitMeasures(postulated, "logl"), 3)
  v2  <- round(fitMeasures(postulated, "unrestricted.logl"), 3)
  v1 
  MQ3 <- v1 < v2
  if(!is.numeric(MQ3)) {
    MQ3 <- FALSE
  }
  cat("\n")
  cat("   Log likelihood from postulated model should be lower than that of the unrestricted model.\n")
  cat("   .. Test:", if (MQ3) "Passed" else "Failed", "   ... (postulated =", v1, ", unrestricted =", v2, ")\n")
  
  v1  <- round(fitMeasures(postulated, "rmsea"), 3)
  v2  <- round(fitMeasures(postulated, "rmsea.pvalue"), 3)
  MQ4 <- v1 < 0.08 && v2 > 0.05
  cat("\n")
  cat("   Root mean square error of approximation (RMSEA).\n")
  cat("   .. Test:", if (MQ4) "Passed" else "Failed", "   ... (RMSEA =", v1, "< 0.08), p-value =", v2, "> 0.05 )\n")
  
  cat("\n")
  cat("   For later reference:\n")
  cat("   AIC:", round(fitMeasures(postulated, "aic"), 3), "\n")
  cat("   BIC:", round(fitMeasures(postulated, "bic"), 3), "\n")

  MQ <- MQ1 && MQ2 && MQ3 && MQ4
  cat("\n")
  if (MQ) {
    cat("   All model quality tests were ok.\n")
  } else {
    cat("** At least one of the model quality tests has failed.\n")
  }

}

showEstimates <- function(postulated) {
  cat("-- Model parameters:\n")
  print(
    parameterEstimates(postulated, standardized=TRUE) %>%
      filter(op == "=~") %>%
      select('Latent Factor'=lhs, Indicator=rhs, B=est, SE=se, Z=z, 'p-value'=pvalue, Beta=std.all) %>%
      kable(digits = 3, format="pandoc", caption="Factor Loadings")
  )
  cat("\n")
  print(warnings())
  
}

showResiduals <- function(postulated) {
  lavTables(postulated)
}

compareModels <- function(postulated, baseline) {
  res <- anova(postulated, baseline)
  print(res)
  MQ5 <- res$`Chisq diff`[2] > 0 && res$`Pr(>Chisq)`[2] < 0.05
  cat("\n")
  cat("-- Result from model comparison:", if (MQ5) "Passed (postulated is more specific)" else "Failed (baseline is mode specific)", "\n")
  cat("\n")
  print(warnings())
}

Ser.model.sf <- ' singlefactor   =~ s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8 + rating'

# specification of the postulated model
Ser.model.p1 <- ' relevance      =~ s7 + rating
                  novelty        =~ s1 + s2
                  unexpectedness =~ s4 + s5 + s6 + s8 + s3'


# specification of the postulated model with s3 excluded
Ser.model.p2 <- ' relevance      =~ s7 + rating
                  novelty        =~ s1 + s2
                  unexpectedness =~ s4 + s5 + s6 + s8'


# specification of the alternative model
Ser.model.p3 <- ' relevance      =~ s7 + rating
                  unexpectedness =~ s4 + s5 + s6 + s8'


# clears the current console content
cat("\014")

# loads the preprocessed data
# (assumes that any measures needed to tackle outliers or missing/invalid data have already been taken
cat("-- Loading the dataset.\n")

list[oserData, tserData, lserData] <- loadSer2018Data(imageDirectory, filename)
nv = ncol(oserData) # nv stands for the number of variables in the dataset
nf = ncol(oserData) # nf stands for the number of factors to be extracted
ss = nrow(oserData) # ss stands for the sample size (number of cases)
serData <- tserData
serVars <- names(tserData) #append(paste("s", 1:8, sep=""), "rating")
cat("   The datafile ", filename, " has", ss, "samples and ", nf, "variables.\n")
cat("-- You are ready to go on with your analyses in interative mode.\n")

# The remainder of this script has been organised as a "buffett": general definitions and data 
# have been loaded, and you can interactively perform analyses and comparisons on your own:

# -- to fit a confirmatory model to the data:
#postulated.nf3 <- induceModel(Ser.model, serData)

# -- to apply some quality tests to the model that was fitted to the data:
# assessModelfit(postulated)

# -- to inspect the model parameters:
# showEstimates(postulated)

# -- to inspect the model residuals:
# showResiduals(postulated)

# -- to inspect suggestes modifications:
# modificationIndices(postulated, sort.=TRUE, minimum.value=3)

# -- to compare models:
#baseline <- cfa(Ser.model, data=serData, std.lv=TRUE, ordered=serVars, orthogonal = TRUE)
#compareModels(postulated, baseline)  

model.sf <- induceModel(Ser.model.sf, serData)
model.p1 <- induceModel(Ser.model.p1, serData)
model.p2 <- induceModel(Ser.model.p2, serData)
model.p3 <- induceModel(Ser.model.p3, serData)

#compareModels(model.p2, model.p1)
#compareModels(model.p3, model.p1)
#compareModels(model.p3, model.p2)

