#-----------------------------------------------------------------------------------------------------------
# This script has been adapted from
# Hartman, Rose. (2019). Understanding Data. 
#   Online http://www.understandingdata.net/2017/03/22/cfa-in-lavaan/#refs
#-----------------------------------------------------------------------------------------------------------

# removes the current environment variables
rm(list = ls())
plot.new()

# environment configuration options
ECO_CUT_LEVEL  = 0.3
ECO_CUT_2CHK   = 0.1
ECO_PRECISION  = 3
ECO_OPTIMETHOD = 'minres'

# installs necessary packages
#install.packages("lavaan")
#install.packages("MVN")

# loads libraries
library(lavaan)
library(MVN)

#---------------------------------------------------------------------------------------------------------
# General definitions
#---------------------------------------------------------------------------------------------------------

Serendip.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(Serendip.model, data=HolzingerSwineford1939, std.lv=TRUE, missing="fiml")              

fitMeasures(fit, "cfi") > 0.9
fitMeasures(fit, "tli") > 0.9

fitMeasures(fit, "logl") < fitMeasures(fit, "unrestricted.logl")

fitMeasures(fit, "aic")
fitMeasures(fit, "bic")

# tests if the model has close fit according to RMSEA
fitMeasures(fit, "rmsea.pvalue") > 0.05




library(dplyr)
library(tidyr)
library(knitr)
options(knitr.kable.NA = '') # this will hide missing values in the kable table

parameterEstimates(fit, standardized=TRUE) %>%
  filter(op == "=~") %>%
  select('Latent Factor'=lhs, Indicator=rhs, B=est, SE=se, Z=z, 'p-value'=pvalue, Beta=std.all) %>%
  kable(digits = 3, format="pandoc", caption="Factor Loadings")




cov_table <- residuals(fit, type = "cor")$cov
cov_table[upper.tri(cov_table)] <- NA # erase the upper triangle
diag(cov_table) <- NA # erase the diagonal 0's
kable(cov_table, format="pandoc", digits=2) # makes a nice table and rounds everyhing to 2 digits
# do we have categorical variables? replace previous snippet with lavTables(fit)


modificationIndices(fit, sort.=TRUE, minimum.value=3)


# comparing models

# 1. postulated model vs postulated model without covariance among the factors

fit_orth <- cfa(Serendip.model, data=HolzingerSwineford1939, std.lv=TRUE,  missing="fiml", orthogonal = TRUE)
fit_orth
fit
anova(fit, fit_orth)

# 2. postulated model vs postulated model with full covariance among the factors

Serendip.model.one <- ' ability  =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 '
fit_one <- cfa(Serendip.model.one, data=HolzingerSwineford1939, std.lv=TRUE,  missing="fiml")
fit_one
fit
anova(fit, fit_one)


# 3. postulated model vs alternative model

Descriptives
         Skew    Kurtosis
x1 -0.2543455  0.30753382
x5 -0.3497961 -0.55253689
x9  0.2038709  0.28990791

x6  0.8579486  0.81655717
x8  0.5252580  1.17155564
x2  0.4700766  0.33239397

x3  0.3834294 -0.90752645
x4  0.2674867  0.08012676
x7  0.2490881 -0.30740386


Serendip.model.alt <- ' visual  =~ x1 + x5
                  textual =~ x2 + x3 + x6 + x8
                  speed   =~ x4 + x7 + x9 '

Serendip.model.alt <- ' visual  =~ x1 + x5 + x9
                  textual =~ x2 + x6 + x8
                  speed   =~ x3 + x4 + x7 '

fit_alt <- cfa(Serendip.model.alt, data=HolzingerSwineford1939, std.lv=TRUE,  missing="fiml")
fit_alt
fit
anova(fit, fit_alt)



# clears the current console content
cat("\014")

cat("---------------------------------------------------------------------------------------------------------\n")
cat("Stage 1 - loads the preprocessed data and applies data quality tests\n")
cat("---------------------------------------------------------------------------------------------------------\n")

# loads the preprocessed data
# (assumes that any measures required to tackling missing and invalid data, as well as outliers, have already been
# taken)
cat("-- Loading the dataset and computing its correlation matrix.\n")
raqData <- read.delim(file.path(imageDirectory, filename), header = TRUE)
nv = ncol(raqData) # nv stands for the number of variables in the dataset
nf = ncol(raqData) # nf stands for the number of factors to be extracted
ss = nrow(raqData) # ss stands for the sample size (number of cases)

# creates the correlation matrix
raqMatrix <- cor(raqData)
round(raqMatrix, 2) #xxx rounding may be hurtful?
cat("   The dataset ", filename, " has", ss, "samples and ", nf, "variables.\n")

cat("\n")
cat("-- Applying data quality tests.\n")

# applies Bartlett's sphericity test
res = cortest.bartlett(raqMatrix, n=ss)
DQ1 = res$p.value < 5e-2
cat("   Bartlett test checks if the correlation matrix is significantly different from the identity.\n")
cat("   .. Batlett test:", if (DQ1) "Passed" else "Failed", "   ... (p-value =", res$p.value, "< 0.05)\n")

# applies a test based on the Keiser-Meyer-Olkin (KMO) measure for sampling adequacy
res = kmo(raqData)
DQ2 = res$overall > .5
cat("   KMO measure assesses the relation between partial and total correlations.\n")
cat("   ..", res$report, "\n")
cat("   .. KMO test:", if (DQ2) "Passed" else "Failed", "       ... (KMO =", res$overall, "> 0.5)\n")

# applies the determinant test to assure correlations are not too high
DQ3 = det(raqMatrix) > 1e-5
cat("   Multicollinearity test checks for issues with too much communality among variables.\n")
cat("   .. Multicollinearity test: ", if (DQ3) "Passed" else "Failed", '\n')

# displays histograms for each variable (indicator) and
res <- displayDataQuality(raqData)
cat("\n")
cat("   Checking the normality for each variable")
cat("\n")
print(res$univariateNormality)
cat("\n")

# checks the test results
if (DQ1 && DQ2 && DQ3) {
  cat("\n")
  cat("   All data quality tests were ok.\n")
  cat("\n")
  readline(prompt="Press [enter] to continue")
} else {
  stop("   At least one of the data quality tests failed.\n")
}
cat("\n")

iter = 0
cat("\n")
cat("---------------------------------------------------------------------------------------------------------\n")
cat("Stage 2 - Performing factor extraction (iteration ", iter, ", #factors =", nf, ")\n")
cat("---------------------------------------------------------------------------------------------------------\n")

# applies PCA to identify eingenvectors
pc1 <- principal(raqMatrix, nfactors = nf, rotate = "none")
#pc1 <- fac(raqMatrix, nfactors = nf, n.obs = ss, fm = ECO_OPTIMETHOD, rotate = "none")
print(pc1)

# displays the scree plot
dev.off()
par(mfrow=c(1,3))
displayScreePlot(pc1$values)

KC1 = FALSE
KC2 = FALSE
RC1 = FALSE
RC2 = FALSE
UC  = FALSE
while (!UC){

  iter = iter + 1

  # IMPORTANT: after consulting the scree plot and applying the Kaiser criteria (eigenvalue > 1),
  #            the number of factors to be extracted is selected and extracted
  #            the Kaiser's criterion was adopted (instead of the Jollife's) because in our initial analysis:
  #            (1) the number of variables was smaller than 30 and communalities were consistently larger than .7
  #            (2) the sample size is larger than 250  and communalities were larger than .6 on average

  user.input <- as.integer(dlgInput("After considering the scree plot and eigenvalues from PCA, what is the number of factors to extract? Type '0' to proceed with the current number of factors.", nf)$res)
  if(user.input == 0) {
    if((!KC1 && !KC2) || (!RC1 || !RC2)) {
      cat("** WARNING: proceeding with the analysis even though the quality criteria have not been fully met.\n\n")
      UC = TRUE
    }
  } else {

    nf = user.input

    cat("\n")
    cat("---------------------------------------------------------------------------------------------------------\n")
    cat("Stage 2 - Performing factor extraction (iteration ", iter, ", #factors =", nf, ")\n")
    cat("---------------------------------------------------------------------------------------------------------\n")

    # applies FA to extract a reduced number of factors
    #pc2 <- principal(raqMatrix, nfactors = nf, rotate = "none")
    pc2 <- fac(raqMatrix, nfactors = nf, n.obs = ss, fm = ECO_OPTIMETHOD, rotate = "none")
    print(pc2)

    # assesses the communality of the selected factors to check if Kaiser's criteria still apply
    # (1) the number of variables is smaller than 30 and extracted communalities are consistently larger than .7
    # (2) the sample size is larger than 250 and and extracted communalities are larger than .6 on average

    cat("\n\n")
    cat("-- Assessing Kaiser's criteria\n")
    cat("   number of variables is ", nv, "\n")
    cat("   number of variables with communality larger than 0.7: ", sum(pc2$communality > .7), "\n")
    cat("   average communality per variable ...................: ", round(mean(pc2$communality), ECO_PRECISION), "(preferably larger than 0.6)\n")
    KC1 = nv < 30  && sum(pc2$communality   > .7) == nv
    KC2 = ss > 250 && mean(pc2$communality) > .6
    cat("   .. Kaiser's criterion #1: ", if (KC1) "Passed" else "Failed", "\n")
    cat("   .. Kaiser's criterion #2: ", if (KC2) "Passed" else "Failed", "\n")
    if (KC1 || KC2) {
      cat("   On Kaiser's criteria, the analysis may proceeed.\n")
    } else {
      cat("   On Kaiser criteria, the analysis should not proceed.\n")
    }
    cat("\n")

    # based on the obtained factor scores, reconstructs the correlation matrix, and compares it to the the original correlation matrix
    # (a.k.a computes the residuals of the original/reconstructed model)
    cat("-- Assessing a relative measure of fit between original and reconstructed correlation matrix\n")
    RC1 = pc2$fit.off > 0.9
    cat("   .. Fit test (residuals) .:", if (RC1) "Passed" else "Failed", "   ... (fit =", round(pc2$fit.off, ECO_PRECISION), "> 0.9)\n\n")

    cat("-- Assessing an absolute measure of fit between original and reconstructed correlation matrix\n")
    residmatrix <- factor.residuals(raqMatrix, pc2$loadings)
    residuals <- as.matrix(residmatrix[upper.tri(residmatrix)])
    propLargeResid = residual.stats(residuals, nf)
    RC2 = propLargeResid < 0.5
    cat("   .. Fit test (residuals) .:", if (RC2) "Passed" else "Failed", "   ... (fit =", propLargeResid, "< 0.5)\n")

    cat("\n")
    cat("-- Number of factors extracted:", nf, "\n")
    cat("\n")

    # updates the graphical display
    dev.off()
    par(mfrow=c(1,3))
    displayScreePlot(pc1$values)
    displayResidualsHist(residuals, nf)

    readline(prompt="Press [enter] to continue")
  }
}

iter = 0
UC = FALSE
while(!UC) {

  iter = iter + 1
  cat("---------------------------------------------------------------------------------------------------------\n")
  cat("Stage 3 - Performing factor rotation (iteration", iter, ", #factors =", nf, ")\n")
  cat("---------------------------------------------------------------------------------------------------------\n")


  if(iter == 1) {
    user.input <- dlgInput("Which rotation method should be applied? Type 'accept' to proceed to the next stage.", "promax")$res
  } else {
    user.input <- dlgInput("Which rotation method should be applied? Type 'accept' to proceed to the next stage.", "accept")$res
  }
  if(user.input == "accept") {
    UC = TRUE
    cat("-- Rotation method selected:", rotationMethod, "\n")
  } else {
    rotationMethod = user.input
    #pc4 <- principal(raqMatrix, nfactors = nf, rotate = rotationMethod)
    pc4 <- fac(raqMatrix, nfactors = nf, n.obs = ss, fm = ECO_OPTIMETHOD, rotate = rotationMethod)
    cat("-- Rotation method selected:", rotationMethod, "\n")
    print.psych(pc4, cut = ECO_CUT_LEVEL, sort = TRUE)

    cat("\n\n")
    cat("-- Factor structure (for visual inspection -- double check if this structure matrix is similar to the pattern matrix)\n")
    print(factor.structure(pc4, cut = ECO_CUT_LEVEL))
    readline(prompt="Press [enter] to continue")
    #xxx what should be done if the factor loadings from the pattern and the structure matrices sharply disagree?
  }
}

# displays the factor-variable interactions graph, with loadings from the pattern matrix
displayFactorGraph(pc4)


cat("\n")
cat("---------------------------------------------------------------------------------------------------------\n")
cat("Stage 4 - Collecting factor scores (#factors =", nf, ")\n")
cat("---------------------------------------------------------------------------------------------------------\n")

# must use original data instead of the correlation matrix
#pc5 <- principal(raqData, nfactors = nf, rotate = rotationMethod, scores = TRUE)
pc5 <- fac(raqData, nfactors = nf, n.obs = ss, fm = ECO_OPTIMETHOD, rotate = rotationMethod, scores = "regression")
newraqData <- cbind(raqData, pc5$scores)
cat("Factor scores have been appended to the original data and is available in the 'newraqData' object\n")

cat("\n")
cat("---------------------------------------------------------------------------------------------------------\n")
cat("Stage 5 - Performing reliability analysis for each subscale (#factors =", nf, ")\n")
cat("---------------------------------------------------------------------------------------------------------\n")
# assumes that the interpretative stage has not been completed by the analyst at this point
# each cluster of highly correlated items has been associated to a concept in the domain of analysis (subscale)

clusters = as.data.frame(factor2cluster(pc4$loadings))
for(factorName in names(clusters)) {
  cat("-- Performing reliability analysis for ", factorName, "\n")
  factorComp <- c()
  factorKeys <- c()
  for(itemKeys in clusters[factorName]) {
    itemID = 0
    for(itemKey in itemKeys) {
      itemID = itemID + 1
      if(itemKey != 0) {
        factorComp <- c(factorComp, itemID)
        factorKeys <- c(factorKeys, itemKey)
      }
    }
  }
  #cat("   Factor composition:", factorComp, "\n")
  #cat("   Factor keys ......:", factorKeys, "\n")
  cat("   Factor composition:", rownames(clusters)[factorComp], "\n")
  cat("   Factor keys ......:", factorKeys, "\n")
  factorData <- raqData[, factorComp]
  if(length(factorComp) > 1) {
    print(psych::alpha(factorData, keys = factorKeys))
  } else {
    cat("   Factor has a single item; no reliability will be computed.\n")
  }
  readline(prompt="Press [enter] to continue")
}


cat("\n")
cat("---------------------------------------------------------------------------------------------------------\n")
cat("The analysis has completed.\n")
cat("---------------------------------------------------------------------------------------------------------\n")

# displays the loadings per variable on each factor
displayLoadings(pc4)

