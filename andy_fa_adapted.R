#-----------------------------------------------------------------------------------------------------------
# This code performs factor analysis on the Serendipity'2018 dataset, which has been gently made available
# by the authors:
# Kotkov, D., Konstan, J. A., Zhao, Q., & Veijalainen, J. (2018, April). Investigating serendipity in
#   recommender systems based on real  user feedback. In Proceedings of the 33rd Annual ACM Symposium on
#   Applied Computing (pp. 1341-1350). ACM.
# The whole dataset is available at: https://grouplens.org/datasets/serendipity-2018/
#
# We would also like to acknolwdge the following authors for the R code for exploratory factor analysis, from
# which we greatly benefitted:
# Field, A. P., Miles, J. N. V., & Field, Z. C. (2012). Discovering Statistics Using R: and Sex and Drugs and
#    Rock 'N' Roll. (pp. ). #London Sage
#-----------------------------------------------------------------------------------------------------------

# removes the current environment variables
rm(list = ls())

# sets the working directory
setwd("D:/Users/Andre/Google Drive/Doutorado/SCC55951 - IHC Practice/factor-analysis")
imageDirectory<-"D:/Users/Andre/Google Drive/Doutorado/SCC55951 - IHC Practice/factor-analysis/datasets/serendipity2018"
filename = "answers.dat"

# installs necessary packages
#install.packages("corpcor")
#install.packages("GPArotation")
#install.packages("psych")
#install.packages("pastecs")
#install.packages("svDialogs")

# loads libraries
library(corpcor)
library(GPArotation)
library(psych)
library(MASS) # used by kmo function
library(svDialogs)

#---------------------------------------------------------------------------------------------------------
# General definitions
#---------------------------------------------------------------------------------------------------------

# KMO Kaiser-Meyer-Olkin Measure of Sampling Adequacy
# Function by G. Jay Kerns, Ph.D., Youngstown State University (http://tolstoy.newcastle.edu.au/R/e2/help/07/08/22816.html)
kmo = function(data){
  X   <- cor(as.matrix(data));
  iX  <- ginv(X);
  S2  <- diag(diag((iX^-1)));
  AIS <- S2%*%iX%*%S2;                      # anti-image covariance matrix
  IS  <- X+AIS-2*S2;                         # image covariance matrix
  Dai <- sqrt(diag(diag(AIS)));
  IR  <- ginv(Dai)%*%IS%*%ginv(Dai);         # image correlation matrix
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

  ans <- list(overall    = kmo,
              report     = test,
              individual = MSA) #,
              #AIS        = AIS,
              #AIR        = AIR )
  return(ans)
}

residual.stats <- function(matrix, nf){
  residuals         <- as.matrix(matrix[upper.tri(matrix)])
  large.resid       <- abs(residuals) > 0.05
  numberLargeResids <- sum(large.resid)
  propLargeResid    <- numberLargeResids/nrow(residuals)
  rmsr              <- sqrt(mean(residuals^2))

  cat("   .. Number of residuals = ", nrow(residuals), "\n")
  cat("   .. Number     of absolute residuals > 0.05 = ", numberLargeResids, "\n")
  cat("   .. Proportion of absolute residuals > 0.05 = ", propLargeResid, "\n")
  cat("   .. Root mean squared residual = ", rmsr, "\n")
  cat("   .. Normality of residuals (Shapiro-Wilk test): ", if (shapiro.test(residuals)$p.value > .05) "Passed" else "Failed", "\n")
  hist(residuals, main = paste("2. Histogram of residuals for ", nf, "factors"))
  Sys.sleep(.1)
  
  return(propLargeResid)
}

factor.structure <- function(fa, cut = 0.2, decimals = 2){
  structure.matrix <- fa.sort(fa$loadings %*% fa$Phi)
  structure.matrix <- data.frame(ifelse(abs(structure.matrix) < cut, "", round(structure.matrix, decimals)))

  return(structure.matrix)
}

# clears the current console content
cat("\014")

#---------------------------------------------------------------------------------------------------------
# Stage 1 - loads the dataset and applies data quality tests
#---------------------------------------------------------------------------------------------------------

cat("---------------------------------------------------------------------------------------------------------\n")
cat("Stage 1 - loads the preprocessed data and applies data quality tests\n")
cat("---------------------------------------------------------------------------------------------------------\n")

# loads the preprocessed data 
# (assumes that any measures required to tackling missing and invalid data, as well as outliers, have already been 
# taken)
cat("-- Loading the dataset and computing its correlation matrix.\n")
raqData <- read.delim(file.path(imageDirectory, filename), header = TRUE)
nv = ncol(raqData) # nd stands for the number of variables in the dataset
nf = ncol(raqData) # nf stands for the number of factors to be extracted
ss = nrow(raqData) # ss stands for the sample size (number of cases)

# creates the correlation matrix
raqMatrix <- cor(raqData)
round(raqMatrix, 2) #xxx rounding may be hurtful?
cat("   The dataset has", ss, "samples and ", nf, "variables.\n")

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

if (DQ1 && DQ2 && DQ3) {
  cat("   All data quality tests were ok.\n")
	readline(prompt="Press [enter] to continue")
} else {
  stop("   At least one of the data quality tests failed.\n")
}

#---------------------------------------------------------------------------------------------------------
# Stage 2 - factor extraction
#---------------------------------------------------------------------------------------------------------
cat("\n")

iter = 0
cat("---------------------------------------------------------------------------------------------------------\n")
cat("Stage 2 - Performing factor extraction (iteration ", iter, ", #factors =", nf, ")\n")
cat("---------------------------------------------------------------------------------------------------------\n")

# applies PCA to identify eingenvectors that satisfies a cascaded optimised variance
pc1 <- principal(raqMatrix, nfactors = nf, rotate = "none")
print(pc1)

KC1 = FALSE
KC2 = FALSE
RC1 = FALSE
RC2 = FALSE
while ((!KC1 && !KC2) || (!RC1 || !RC2)){
  
  iter = iter + 1
  
  # shows the scree plot (integral, non-rotated PCA)
  par(mfrow=c(1,2))
  plot(pc1$values, type = "b", main = "1. Scree Plot from integral, non-rotated PCA", xlab = "component", ylab = "eigenvalue")
  abline(h=1, col="red", lty=2, lwd=1)
  text(1.2, 1.1, "Kaiser", col="red")
  abline(h=.7, col="red", lty=2, lwd=1)
  text(1.2, 0.8, "Jollife", col="red")
  Sys.sleep(.1)

  # IMPORTANT: after consulting the scree plot and applying the Kaiser criteria (eigenvalue > 1),
  #            the number of factors to be extracted is selected and extracted
  #            the Kaiser's criterion was adopted (instead of the Jollife's) because in our initial analysis:
  #            (1) the number of variables is smaller than 30 and extracted communalities are consistently larger than .7
  #            (2) the sample size is larger than 250 and and extracted communalities are larger than .6 on average
  user.input <- as.integer(dlgInput("After considering the scree plot and eigenvalues from PCA, what is the number of factors to extract? Type '0' to proceed with the current number of factors.", 1)$res)
  if(user.input == 0) {
    cat("** WARNING: proceeding with the analysis even though the quality criteria have not been fully met.\n\n")
    break
  } else {
    nf = user.input
  }
  
  
  #cat("\014")
  cat("\n")
  cat("---------------------------------------------------------------------------------------------------------\n")
  cat("Stage 2 - Performing factor extraction (iteration ", iter, ", #factors =", nf, ")\n")
  cat("---------------------------------------------------------------------------------------------------------\n")
  
  # applies PCA to extract a reduced number of eigenvectors
  pc2 <- principal(raqMatrix, nfactors = nf, rotate = "none")
  print(pc2)
  
  # assesses the communality of the selected factors to check if Kaiser's criteria still apply
  # (1) the number of variables is smaller than 30 and extracted communalities are consistently larger than .7
  # (2) the sample size is larger than 250 and and extracted communalities are larger than .6 on average
  
  cat("\n\n")
  cat("-- Assessing Kaiser's criteria\n")
  cat("   number of variables is ", nv, "\n")
  cat("   extracted communalities consistently larger than 0.7: ", sum(pc2$communality > .7), "\n")
  cat("   extracted communalities on average: ", mean(pc2$communality), "(preferably larger than 0.6)\n")
  KC1 = nv < 30  && sum(pc2$communality > .2) == nv
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
  cat("   .. Fit test (residuals):", if (RC1) "Passed" else "Failed", "   ... (fit =", pc2$fit.off, "> 0.9)\n\n")
  
  cat("-- Assessing an absolute measure of fit between original and reconstructed correlation matrix\n")
  propLargeResid = residual.stats(factor.residuals(raqMatrix, pc2$loadings), nf)
  RC2 = propLargeResid < 0.5
  cat("   .. Fit test (residuals):", if (RC2) "Passed" else "Failed", "   ... (fit =", propLargeResid, "< 0.5)\n")
  
  cat("\n")
  readline(prompt="Press [enter] to continue")
  
}

iter = 0
UC1 = FALSE
while(!UC1) {

  iter = iter + 1
  cat("---------------------------------------------------------------------------------------------------------\n")
  cat("Stage 3 - Performing factor rotation (iteration", iter, ", #factors =", nf, ")\n")
  cat("---------------------------------------------------------------------------------------------------------\n")
  
  user.input <- dlgInput("Which rotation method should be applied? Type 'accept' to proceed to the next stage.", "promax")$res
  if(user.input != "accept") {

    rotationMethod = user.input
    pc4 <- principal(raqMatrix, nfactors = nf, rotate = rotationMethod)
    cat("-- Rotation method selected:", rotationMethod, "\n")
    print.psych(pc4, cut = 0.3, sort = TRUE) #xxx should we change the cut level?
    
    cat("\n\n")
    cat("-- Factor structure (for visual inspection -- double check if this structure matrix is similar to the pattern matrix)\n")
    print(factor.structure(pc4, cut = 0.3))
    readline(prompt="Press [enter] to continue")
  } else {
    UC1 = TRUE
    cat("-- Rotation method selected:", rotationMethod, "\n")
  }
}
 
#xxx what should be done if factor loadings appear in different ranks across pattern and factor matrices? 

cat("\n")
cat("---------------------------------------------------------------------------------------------------------\n")
cat("Stage 4 - Collecting factor scores (#factors =", nf, ")\n")
cat("---------------------------------------------------------------------------------------------------------\n")

# must use original data instead of the correlation matrix
pc5 <- principal(raqData, nfactors = nf, rotate = rotationMethod, scores = TRUE)
newraqData <- cbind(raqData, pc5$scores)
cat("Factor scores have been appended to the original data and is available in the 'newraqData' variable.\n")

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
  cat("   Factor composition:", factorComp, "\n")
  cat("   Factor keys ......:", factorKeys, "\n")
  factorData <- raqData[, factorComp]
  print(alpha(factorData, keys = factorKeys))
  readline(prompt="Press [enter] to continue")
}

cat("\n")
cat("---------------------------------------------------------------------------------------------------------\n")
cat("The analysis has completed.\n")
cat("---------------------------------------------------------------------------------------------------------\n")
