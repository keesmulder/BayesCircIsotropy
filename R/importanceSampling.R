
# Function to find the mode of a vm kernel.
findVMPostMode <- function(vmPostKernel, th_bar = 0, startingValues = c(th_bar, 1)) {
  optim(startingValues, function(x) -vmPostKernel(x))$par
}





# Function to generate a Bivariate normal distribution based on a conjugate von Mises posterior.
generateBinormImportanceDensityFromConjugate <- function(th, covShrink = 1, R0 = 0, c0 = 0, Qmc = 1000, independent = TRUE) {

  # Get a sample from the posterior.
  postSample <- vmConjugateSample(th, R0 = R0, c0 = c0, Q = Qmc)

  # Estimated covariance matrix
  postCov <- cov(postSample)/covShrink
  if (independent) postCov[1, 2] <- postCov[2, 1] <- 0

  # The von Mises log-posterior function
  vMPost <- vmConjugatePostKernel(th, R0 = R0, c0 = c0, log = TRUE)

  # Find MAPs, the posterior modes
  postMode <- optim(c(meanDir(th), 1), function(x) -vMPost(x))$par

  list(
    pdf         = function(x) dmvnorm(x, mean = postMode, sigma = postCov),
    samplingFun = function(n) rmvnorm(n, mean = postMode, sigma = postCov),
    mean        = postMode,
    sigma       = postCov
  )
}



# General function to generate a Bivariate normal distribution, given a posterior kernel.
generateBinormImportanceDensityPostKernel <- function(postKernel, covShrink = 1, Qmc = 1000, independent = TRUE) {

  # Get a sample from the posterior.
  postSample <- vmConjugateSample(th, R0 = R0, c0 = c0, Q = Qmc)

  # Estimated covariance matrix
  postCov <- cov(postSample)/covShrink
  if (independent) postCov[1, 2] <- postCov[2, 1] <- 0

  # The von Mises log-posterior function
  vMPost <- vmConjugatePostKernel(th, R0 = R0, c0 = c0, log = TRUE)

  # Find MAPs, the posterior modes
  postMode <- optim(c(meanDir(th), 1), function(x) -vMPost(x))$par

  list(
    pdf         = function(x) dmvnorm(x, mean = postMode, sigma = postCov),
    samplingFun = function(n) rmvnorm(n, mean = postMode, sigma = postCov),
    mean        = postMode,
    sigma       = postCov
  )
}


# Estimate marginal likelihood from importance sample
estimateMargLikImpSamp <- function(posteriorKernel, importancePdf, importanceSamplingFun,
                                   Qis = 10000) {

  # Obtain a sample from the importance density
  importanceSample <- importanceSamplingFun(Qis)

  # Evaluate the posterior kernel and importance density at the values from the importance sample.
  posteriorProb  <- apply(importanceSample, 1, posteriorKernel)
  importanceProb <- apply(importanceSample, 1, importancePdf)

  mean(posteriorProb / importanceProb)
}
