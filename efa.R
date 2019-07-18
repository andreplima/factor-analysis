#-----------------------------------------------------------------------------------------------------------
# This code performs factor analysis on the Serendipity'2018 dataset, which has been generously made 
# available by the authors:
# Kotkov, D., Konstan, J. A., Zhao, Q., & Veijalainen, J. (2018, April). Investigating serendipity in
#   recommender systems based on real user feedback. In Proceedings of the 33rd Annual ACM Symposium on
#   Applied Computing (pp. 1341-1350). ACM.
# The complete dataset is available at: https://grouplens.org/datasets/serendipity-2018/
#
# We would also like to acknowledge the following authors for the R code for exploratory factor analysis, 
# from which we benefited greatly:
# Revelle, W. (2018) psych: Procedures for Personality and Psychological Research, Northwestern University,
#    Evanston, Illinois, USA, https://CRAN.R-project.org/package=psych Version = 1.8.12.
# Field, A. P., Miles, J. N. V., & Field, Z. C. (2012). Discovering Statistics Using R. (pp. 749-811). London, Sage.
#-----------------------------------------------------------------------------------------------------------

# removes the current environment variables
rm(list = ls())
plot.new()

# sets up some environment configuration options
ECO_CUT_LEVEL  = 0.3
ECO_CUT_2CHK   = 0.1
ECO_PRECISION  = 3
#ECO_OPTIMETHOD = 'wls'
ECO_OPTIMETHOD = 'minres'


# sets the working directory
setwd("C:/Users/andre/OneDrive/Documentos/gitrepos/factor-analysis")
datasetDirectory <- "C:/Users/andre/OneDrive/Documentos/gitrepos/factor-analysis/datasets/serendipity2018"
outputDirectory  <- "C:/Users/andre/OneDrive/Documentos/gitrepos/factor-analysis/outputs"
#filename = "answers-all-all-median.dat"
#filename = "answers-1399-all-median.dat"
filename = "answers-606-all-median.dat"


# installs necessary packages
#install.packages("corpcor")
#install.packages("GPArotation")
#install.packages("psych")
#install.packages("pastecs")
#install.packages("svDialogs")
#install.packages("reshape2")
#install.packages("ggplot2")

# loads libraries
library(corpcor)
library(GPArotation)
library(psych)
library(MASS) # used by kmo function
library(svDialogs)
library(reshape2)
library(ggplot2)
library(grid)
library(MVN)
library(likert)
library(polycor)
library(gsubfn)

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

# KMO Kaiser-Meyer-Olkin Measure of Sampling Adequacy
# Function by G. Jay Kerns, Ph.D., Youngstown State University (http://tolstoy.newcastle.edu.au/R/e2/help/07/08/22816.html)
# Adapted by Andre P Lima for using with a polychoric correlation matrix
kmo = function(data, corMatrix){
  #X   <- cor(as.matrix(data));
  X   <- corMatrix$correlations;
  iX  <- ginv(X);
  S2  <- diag(diag((iX^-1)));
  AIS <- S2%*%iX%*%S2;                      # anti-image covariance matrix
  IS  <- X+AIS-2*S2;                        # image covariance matrix
  Dai <- sqrt(diag(diag(AIS)));
  IR  <- ginv(Dai)%*%IS%*%ginv(Dai);        # image correlation matrix
  AIR <- ginv(Dai)%*%AIS%*%ginv(Dai);       # anti-image correlation matrix
  a   <- apply((AIR - diag(diag(AIR)))^2, 2, sum);
  AA  <- sum(a);
  b   <- apply((X - diag(nrow(X)))^2, 2, sum);
  BB  <- sum(b);
  MSA <- b/(b+a);                        # indiv. measures of sampling adequacy
  AIR <- AIR-diag(nrow(AIR))+diag(MSA);  # Examine the anti-image of the correlation matrix, i.e., the negative of the partial
                                         # correlations, contitioned on all other variables.
  kmo <- BB/(AA+BB);                     # overall KMO statistic

  # reports the result
  if      (kmo >= 0.00 && kmo < 0.50) {test <- 'The KMO test yields a degree of common variance unacceptable for FA.'}
  else if (kmo >= 0.50 && kmo < 0.60) {test <- 'The KMO test yields a degree of common variance miserable.'}
  else if (kmo >= 0.60 && kmo < 0.70) {test <- 'The KMO test yields a degree of common variance mediocre.'}
  else if (kmo >= 0.70 && kmo < 0.80) {test <- 'The KMO test yields a degree of common variance middling.' }
  else if (kmo >= 0.80 && kmo < 0.90) {test <- 'The KMO test yields a degree of common variance meritorious.' }
  else                                {test <- 'The KMO test yields a degree of common variance marvelous.' }

  ans <- list(overall    = round(kmo, ECO_PRECISION),
              report     = test,
              individual = MSA) #,
              #AIS        = AIS,
              #AIR        = AIR )
  return(ans)
}

residual.stats <- function(residuals, nf){

  large.resid       <- abs(residuals) > 0.05
  numberLargeResids <- sum(large.resid)
  propLargeResid    <- numberLargeResids/nrow(residuals)
  rmsr              <- sqrt(mean(residuals^2))

  cat("   .. Number of residuals = ", nrow(residuals), "\n")
  cat("   .. Number     of absolute residuals > 0.05 = ", numberLargeResids, "\n")
  cat("   .. Proportion of absolute residuals > 0.05 = ", propLargeResid, "\n")
  cat("   .. Root mean squared residual = ", rmsr, "\n")
  cat("   .. Normality of residuals:", if (shapiro.test(residuals)$p.value > .05) "Passed" else "Failed", "(Shapiro-Wilk test)\n")

  return(propLargeResid)
}

factor.structure <- function(fa, cut = ECO_CUT_2CHK, decimals = 2){
  structure.matrix <- fa.sort(fa$loadings %*% fa$Phi)
  structure.matrix <- data.frame(ifelse(abs(structure.matrix) < cut, "", round(structure.matrix, decimals)))

  return(structure.matrix)
}

displayLikert <- function(data_likert) {
  plot(data_likert, include.histogram = TRUE)
  Sys.sleep(.1)
}

displayScreePlot <- function(eigenvals) {
  # draws the scree plot (integral, non-rotated PCA)
  par(mfg=c(1,3))
  plot(eigenvals, type = "b", main = "Scree Plot, non-rotated PCA", xlab = "component", ylab = "eigenvalue", cex.main=1.5, cex.axis=1.2, cex.lab=1.5)
  abline(h=1, col="red", lty=2, lwd=1)
  text(1.2, 1.1, "Kaiser", col="red")
  abline(h=.7, col="red", lty=2, lwd=1)
  text(1.2, 0.8, "Jollife", col="red")
  Sys.sleep(0.1)
}

displayResidualsHist <- function(residuals, nf) {
  # plots the residuals histogram
  par(mfg=c(1,2))
  hist(residuals, main = paste("2. Histogram of residuals for ", nf, "factors"), col="gray", freq=FALSE)
  x <- seq(min(residuals), max(residuals), by=10^-ECO_PRECISION)
  curve(dnorm(x, mean=0, sd=sd(residuals)), add=TRUE, col="red")
  Sys.sleep(0.1)
}

displayFactorGraph <- function(results) {
  # plots the factor-variable graph
  par(mfg=c(1,3))
  fa.diagram(results, main = "3. Factors and standardised loadings\n(from the pattern matrix)", marg = c(.5,.5, 5,.5), cex.main=1.5, cex.axis=1.2, cex.lab=1.5)
  Sys.sleep(0.1)
}

displayLoadings <- function(results) {

  plot.new()
  
  # The results$loading element is an S3 object and cannot be directly coerced into a dataframe.
  # Removing the class attribute gives you a matrix which can then be passed to 'melt'
  # Tip from https://stackoverflow.com/questions/15585870/psych-getting-factor-loadings-as-data-frame-for-latex-export
  loadings.m <- melt(unclass(fa.sort(results$loadings)), varnames=c("Item", "Factor"), value.name="Loading")

  # For each item from the survey questionnaire, plots the loading as length and fill color of a bar
  # note that the length will be the absolute value of the loading but the fill color will be the signed value,
  # more on this below
  
  # And then you ask me: why "print" a ggplot object? That's why:
  # https://stackoverflow.com/questions/38068774/rstudio-suddenly-stopped-showing-plots-in-in-plot-pane
  
  print(
    ggplot(loadings.m, aes(Item, abs(Loading), fill=Loading)) +
    coord_flip()                 + # flips  the axes so the items appear in the (common) y axis
    ylab("Loading Strength")     + # improves y-axis label
    facet_wrap(~ Factor, nrow=1) + # places the factors in separate facets
    geom_bar(stat="identity")    + # makes  the bars
    # defines the fill color gradient: blue=positive, red=negative
    scale_fill_gradient2(name = "Loading", high = "red", mid = "white", low = "blue", midpoint=0, guide=F) +
    theme_bw(base_size=12)         # uses a black-and-white theme with set font size
  )
  
  Sys.sleep(.1)
}

# clears the current console content
cat("\014")

cat("---------------------------------------------------------------------------------------------------------\n")
cat("Stage 1 - loads the preprocessed data and applies data quality tests\n")
cat("---------------------------------------------------------------------------------------------------------\n")

# loads the preprocessed data
# (assumes that any measures needed to tackle outliers or missing/invalid data have already been taken
cat("-- Loading the dataset and computing its correlation matrix.\n")

list[oserData, serData, serDataLk] <- loadSer2018Data(datasetDirectory, filename)
nv = ncol(serData) # nv stands for the number of variables in the dataset
nf = ncol(serData) # nf stands for the number of factors to be extracted
ss = nrow(serData) # ss stands for the sample size (number of cases)
cat("   The datafile ", filename, " has", ss, "samples and ", nf, "variables.\n")

# creates the correlation matrix
serMatrix <- hetcor(serData)

cat("\n")
cat("-- Applying data quality tests.\n")

# applies Bartlett's sphericity test
res = cortest.bartlett(serMatrix$correlations, n=ss)
DQ1 = res$p.value < 5e-2
cat("   Bartlett test checks if the correlation matrix is significantly different from the identity.\n")
cat("   .. Batlett test:", if (DQ1) "Passed" else "Failed", "   ... (p-value =", res$p.value, "< 0.05)\n")

# applies a test based on the Keiser-Meyer-Olkin (KMO) measure for sampling adequacy
res = kmo(serData, serMatrix)
DQ2a = res$overall > .5
DQ2b = sum(res$individual < .5) == 0
DQ2 = DQ2a && DQ2b
cat("   KMO measure assesses the relation between partial and total correlations.\n")
cat("   ..", res$report, "\n")
cat("   .. KMO test ...:", if (DQ2a) "Passed" else "Failed", "   ... (KMO =", res$overall, "> 0.5)\n")
cat("   .. KMO MSA ....:", if (DQ2b) "Passed" else "Failed", "   ... (", sum(res$individual < .5), " variables with KMO under 0.5)\n")

# applies the determinant test to assure correlations are not too high
DQ3 = det(serMatrix$correlations) > 1e-5
cat("   Multicollinearity test checks for issues with too much communality among variables.\n")
cat("   .. Multicollinearity test: ", if (DQ3) "Passed" else "Failed", '\n')

# displays histograms for each variable (indicator) and
displayLikert(serDataLk)

# checks the test results
if (DQ1 && DQ2 && DQ3) {
  cat("\n")
  cat("   All data quality tests were ok.\n")
  cat("\n")
  readline(prompt="Press [enter] to continue")
} else {
  cat("\n")
  cat("** At least one of the data quality tests has failed.\n")
  user.input <- as.integer(dlgInput("At least one of the data quality tests has failed. Type '0' to proceed with the analysis.", 0)$res)
  if(user.input == 0) {
    cat("** WARNING: proceeding with the analysis even though the quality criteria have not been fully met.\n")
  } else {
    stop("** The analysis was interrupted because violations of distributional assumptions were found.\n")
  }
}
cat("\n")

iter = 0
cat("\n")
cat("---------------------------------------------------------------------------------------------------------\n")
cat("Stage 2 - Performing factor extraction (iteration ", iter, ", #factors =", nf, ")\n")
cat("---------------------------------------------------------------------------------------------------------\n")

# applies PCA to identify eingenvectors
pc1 <- principal(serMatrix$correlations, nfactors = nf, rotate = "none")
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
    UC = TRUE
    if((!KC1 && !KC2) || (!RC1 || !RC2)) {
      cat("** WARNING: proceeding with the analysis even though the quality criteria have not been fully met.\n\n")
    }
  } else {

    nf = user.input

    cat("\n")
    cat("---------------------------------------------------------------------------------------------------------\n")
    cat("Stage 2 - Performing factor extraction (iteration ", iter, ", #factors =", nf, ")\n")
    cat("---------------------------------------------------------------------------------------------------------\n")

    # applies FA to extract a reduced number of factors
    pc2 <- fac(serMatrix$correlations, nfactors = nf, n.obs = ss, fm = ECO_OPTIMETHOD, rotate = "none")
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
      cat("   On Kaiser's criteria, the analysis should not proceed.\n")
    }
    cat("\n")

    # based on the obtained factor scores, reconstructs the correlation matrix, and compares it to the the original correlation matrix
    # (a.k.a computes the residuals of the original/reconstructed model)
    cat("-- Assessing a relative measure of fit between original and reconstructed correlation matrix\n")
    RC1 = pc2$fit.off > 0.9
    cat("   .. Fit test (residuals) .:", if (RC1) "Passed" else "Failed", "   ... (fit =", round(pc2$fit.off, ECO_PRECISION), "> 0.9)\n\n")

    cat("-- Assessing an absolute measure of fit between original and reconstructed correlation matrix\n")
    residmatrix <- factor.residuals(serMatrix$correlations, pc2$loadings)
    residuals <- as.matrix(residmatrix[upper.tri(residmatrix)])
    propLargeResid = residual.stats(residuals, nf)
    RC2 = propLargeResid < 0.5
    cat("   .. Fit test (residuals) .:", if (RC2) "Passed" else "Failed", "   ... (prop =", propLargeResid, "< 0.5)\n")

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
    pc4 <- fac(serMatrix$correlations, nfactors = nf, n.obs = ss, fm = ECO_OPTIMETHOD, rotate = rotationMethod)
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
cat("Stage 4 - Estimating factor scores (#factors =", nf, ")\n")
cat("---------------------------------------------------------------------------------------------------------\n")

# must use original data instead of the correlation matrix
pc5 <- fac(oserData, nfactors = nf, n.obs = ss, fm = ECO_OPTIMETHOD, rotate = rotationMethod, scores = "regression", cor="poly")
newserData <- cbind(serData, pc5$scores)
cat("Factor scores have been appended to the original data and is available in the 'newserData' object\n")

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
  cat("   Factor composition:", rownames(clusters)[factorComp], "\n")
  cat("   Factor keys ......:", factorKeys, "\n")
  factorData <- oserData[, factorComp]
  if(length(factorComp) > 1) {
    print(psych::alpha(factorData, keys = factorKeys))
  } else {
    cat("   Factor has a single item; no reliability will be computed.\n")
  }
  readline(prompt="Press [enter] to continue")
}

# displays the loadings per variable on each factor
displayLoadings(pc4)

cat("\n")
cat("---------------------------------------------------------------------------------------------------------\n")
cat("The analysis has completed.\n")
cat("---------------------------------------------------------------------------------------------------------\n")
