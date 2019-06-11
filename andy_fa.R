#-----------------------------------------------------------------------------------------------------------
# This script has been adapted from
# Field, A. P., Miles, J. N. V., & Field, Z. C. (2012). Discovering Statistics Using R: and Sex and Drugs and
#    Rock 'N' Roll. (pp. ). #London Sage
#-----------------------------------------------------------------------------------------------------------

# sets the working directory
setwd("D:/Users/Andre/Google Drive/Doutorado/SCC55951 - IHC Practice/Workarea")
imageDirectory<-"D:/Users/Andre/Google Drive/Doutorado/SCC55951 - IHC Practice/Workarea/datasets/dsur"
filename = "raq.dat"

# installs packages
#install.packages("corpcor")
#install.packages("GPArotation")
#install.packages("psych")
#install.packages("pastecs")

# loads libraries
library(corpcor)
library(GPArotation)
library(psych)
library(MASS) # used by kmo function

#---------------------------------------------------------------------------------------------------------
# Problem general definitions
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
  AIR <- AIR-diag(nrow(AIR))+diag(MSA);  # Examine the anti-image of the correlation matrix, i.e., the negative of the partial correlations,
                                         # contitioned on all other variables.
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

residual.stats <- function(matrix){
  residuals         <- as.matrix(matrix[upper.tri(matrix)])
  large.resid       <- abs(residuals) > 0.05
  numberLargeResids <- sum(large.resid)
  propLargeResid    <- numberLargeResids/nrow(residuals)
  rmsr              <- sqrt(mean(residuals^2))

  cat("Root means squared residual = ", rmsr, "\n")
  cat("Number of absolute residuals > 0.05 = ", numberLargeResids, "\n")
  cat("Proportion of absolute residuals > 0.05 = ", propLargeResid, "\n")
  hist(residuals)
}

factor.structure <- function(fa, cut = 0.2, decimals = 2){
  structure.matrix <- fa.sort(fa$loadings %*% fa$Phi)
  structure.matrix <- data.frame(ifelse(abs(structure.matrix) < cut, "", round(structure.matrix, decimals)))

  return(structure.matrix)
}

cat("\014")

#---------------------------------------------------------------------------------------------------------
# Stage 1 - initial preparation
#---------------------------------------------------------------------------------------------------------

# loads the preprocessed data (assumes that any measures required to tackling missing and invalid data, as well as outliers, have been taken already)
raqData <- read.delim(file.path(imageDirectory, filename), header = TRUE)
nf = ncol(raqData) # number of factors to be extracted
ss = nrow(raqData) # sample size (number of cases)

# creates the correlation matrix
raqMatrix <- cor(raqData)
round(raqMatrix, 2) #xxx rounding may be hurtful?

# applies Bartlett's sphericity test
res = cortest.bartlett(raqMatrix, n=ss)
cat("Bartlett test checks if the R matrix is significantly different from the identity. p-value =", res$p.value)
if (res$p.value < 5e-2) "Batlett: Passed" else "Bartlet: Failed"

# applies KMO test for sampling adequacy
res = kmo(raqData)
res$report

# computes the determinant
if (det(raqMatrix) > 1e-5) "Multicollinearity: not a problem" else "Multicollinearity: there is an issue with the data"

#---------------------------------------------------------------------------------------------------------
# Stage 2 - factor extraction
#---------------------------------------------------------------------------------------------------------

# applies PCA to identify eingenvectors that satisfies a cascaded optimised variance
pc1 <- principal(raqMatrix, nfactors = nf, rotate = "none")

# shows the scree plot
plot(pc1$values, type = "b")

# IMPORTANT: after consulting the scree plot and applying the Kaiser criteria,
#            the number of factors to be extracted is selected and extracted
nf = 4
pc2 <- principal(raqMatrix, nfactors = nf, rotate = "none")

# explores the residuals
#factor.model(pc2$loadings)
#factor.residuals(raqMatrix, pc2$loadings)
pc2$fit.off

#residuals<-factor.residuals(raqMatrix, pc2$loadings)
#residuals<-as.matrix(residuals[upper.tri(residuals)])
#large.resid<-abs(residuals) > 0.05
#sum(large.resid)
#sum(large.resid)/nrow(residuals)
#sqrt(mean(residuals^2))
#hist(residuals)

#resids <- factor.residuals(raqMatrix, pc2$loadings)
#residual.stats(resids)
residual.stats(factor.residuals(raqMatrix, pc2$loadings))

#---------------------------------------------------------------------------------------------------------
# Stage 3 - factor rotation
#---------------------------------------------------------------------------------------------------------

pc4 <- principal(raqMatrix, nfactors = nf, rotate = "oblimin") #xxx trocar para promax
print.psych(pc4, cut = 0.3, sort = TRUE) #readaptar o cut
#pc4$loadings%*%pc4$Phi

factor.structure(pc4, cut = 0.3)

#---------------------------------------------------------------------------------------------------------
# Stage 4 - factor scores
#---------------------------------------------------------------------------------------------------------
# must use original data instead of the correlation matrix
pc5 <- principal(raqData, nfactors = 4, rotate = "oblimin", scores = TRUE)
#pc5$scores
head(pc5$scores, 10)
newraqData <- cbind(raqData, pc5$scores)

# rebuilds the correlation between extracted factors
# --- check with results obtained from print.psych
cor(pc5$scores);
round(cor(pc5$scores), 2)

# rebuilds a row in the structure matrix (an entry for a specific question/item)
# --- check with results from factor.structure
round(cor(pc5$scores, raqData$Q01),2)
round(cor(pc5$scores, raqData$Q06),2)
round(cor(pc5$scores, raqData$Q18),2)

#---------------------------------------------------------------------------------------------------------
# Stage 5 - reliability analysis
#---------------------------------------------------------------------------------------------------------
# assumes the interpretative stage has been successfully completed by the analyst
# each cluster of highly correlated items has been associated to a concept in the domain of analysis (subscale)
# note that items are monotonically ordered in each subscale
computerFear<-raqData[,c(6, 7, 10, 13, 14, 15, 18)]
statisticsFear <- raqData[, c(1, 3, 4, 5, 12, 16, 20, 21)]
mathFear <- raqData[, c(8, 11, 17)]
peerEvaluation <- raqData[, c(2, 9, 19, 22, 23)]

# computes the Cronbach's alpha for each subscale
# --- check for items with raw alpha higher than the overall subscale alpha
alpha(computerFear)
alpha(statisticsFear, keys = c(1, -1, 1, 1, 1, 1, 1, 1)) # keys must be informed when there are items with negative correlation
alpha(mathFear)
alpha(peerEvaluation)
alpha(statisticsFear) #for illustrative pruposes
