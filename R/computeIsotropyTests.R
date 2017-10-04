
# Compute a set of frequentist tests for Isotropy.
frequentistIsotropyTests <- function(th) {

  # Prevent the warnings thrown by the circular package.
  oldw <- getOption("warn")
  options(warn = -1)

  ray  <- circular::rayleigh.test(th)
  kui  <- circular::kuiper.test(th)
  rao  <- circular::rao.spacing.test(th)
  wat  <- circular::watson.test(th)

  options(warn = oldw)

  c(rayStat = ray$statistic, rayp = ray$p.value,
    kuiStat = kui$statistic,
    raoStat = rao$statistic,
    watStat = wat$statistic)

}

computeIsotropyMarginalLikelihood    <- function(th) (2*pi) ^ (-length(th))
computeLogIsotropyMarginalLikelihood <- function(th) -length(th) * log(2*pi)

computeIsotropyDIC    <- function(th)  2 * length(th) * log( 2*pi )



# Computes the bayes factor in favor of the alternative hypothesis of von Mises
# distributed data.
# This version uses a conjugate prior.
# conj_prior has mu0, R0, and c0, in that order.
computeIsotropyBFHaVMConj <- function(th,
                                      conj_prior,
                                      kappaMax = 20) {

  # Sufficient statistics for the von Mises posterior, with conjugate prior values
  n  <- length(th)
  R  <- getR(th)
  R0 <- conj_prior[2]
  c0 <- conj_prior[3]


  # First integral to be computed numerically, the inverse of the normalizing
  # constant of the prior.
  NCPrior <- integrate(function(kp) {
    exp(logBesselI(R0 * kp, nu = 0) - c0 * logBesselI(kp, nu = 0))
  }, lower = 0, upper = kappaMax)$value

  # Second integral to be computed numerically.
  uniF <- function(kp) {
    exp(logBesselI(R0 * kp, nu = 0) + logBesselI(R * kp, nu = 0) -
          (n + c0) * logBesselI(kp, nu = 0))
  }

  partialIntegral <- integrate(uniF, lower = 0, upper = kappaMax)$value

  postLogMargLik <-  log(partialIntegral) - log(NCPrior) - n * log(2 * pi)

  # Compute the marginal likelihood of the null hypothesis.
  isoLogMargLik <- computeLogIsotropyMarginalLikelihood(th)

  # Compute the log and regular Bayes Factor
  logBF <- postLogMargLik - isoLogMargLik
  BF    <- exp(logBF)

  # The Bayes Factor is the ratio of marginal likelihoods.
  tibble::lst(BF, logBF,
              pHa                     = BF / (1 + BF),
              pH0                     = 1 / (1 + BF),
              bfChoseHa               = BF > 1,
              HaLogMargLik = postLogMargLik,
              H0LogMargLik = isoLogMargLik)
}



# Computes the bayes factor in favor of the alternative hypothesis of von Mises
# distributed data.
# This version is used for the Jeffreys prior.
computeIsotropyBFHaVMJeffreysPrior <- function(th,
                                               kappaMin = 0,
                                               kappaMax = 20) {

  # Generate a log Jeffreys Prior.
  logPriorFunction <- createVMJeffreysLogPrior(truncLower = kappaMin,
                                               truncUpper = kappaMax)

  # Sufficient statistics for the concentration in a von Mises likelihood
  n  <- length(th)
  R  <- getR(th)


  # The integral to be computed numerically.
  uniF <- function(kp) {
    exp(logPriorFunction(kp) +
          logBesselI(R * kp, nu = 0) -
          n * logBesselI(kp, nu = 0))
  }

  # The 2 * pi in the BF comes from the prior which is circular uniform over the
  # mean direction.
  BF <- integrate(uniF, lower = kappaMin, upper = kappaMax)$value * (2 * pi)
  logBF <- log(BF)
  postLogMargLik <-  logBF - n * log(2 * pi)

  # Compute the marginal likelihood of the null hypothesis.
  isoLogMargLik <- computeLogIsotropyMarginalLikelihood(th)

  # Return a list of results.
  tibble::lst(BF, logBF,
              pHa                     = BF / (1 + BF),
              pH0                     = 1 / (1 + BF),
              bfChoseHa               = BF > 1,
              HaLogMargLik = postLogMargLik,
              H0LogMargLik = isoLogMargLik)
}



# Computes the bayes factor in favor of the alternative hypothesis of von Mises
# distributed data.
# This version uses a general prior.
computeIsotropyBFHaVMGeneral <- function(th,
                                         likelihoodGenFunc,
                                         priorFunction,
                                         kappaMax = 20) {

  # Obtain the posterior.
  postKernel <- createPosteriorFunction(data = th,
                                        likelihoodGenFunc = likelihoodGenFunc,
                                        priorFunction = priorFunction)

  # Compute the marginal likelihood of the alternative hypothesis with bivariate integration.
  postMargLik <- cubature::adaptIntegrate(postKernel,
                                          lower = c(0, 0),
                                          upper = c(2*pi, kappaMax))$integral

  # Compute the marginal likelihood of the null hypothesis.
  isoLogMargLik <- computeLogIsotropyMarginalLikelihood(th)

  # Compute the log and regular Bayes Factor
  logBF <- log(postMargLik) - isoLogMargLik
  BF    <- exp(logBF)

  # The Bayes Factor is the ratio of marginal likelihoods.
  tibble::lst(BF, logBF,
              pHa                     = BF / (1 + BF),
              pH0                     = 1 / (1 + BF),
              bfChoseHa               = BF > 1,
              HaLogMargLik = log(postMargLik),
              H0LogMargLik = isoLogMargLik)
}


# Compute the Bayes Factor in favor of the kernel density alternative.
# kappaMax MUST MATCH THE PRIOR FUNCTION MAXIMUM or the computed posterior mode might not work
computeIsotropyBFHaVMKDE <- function(th,
                                     priorFunction,
                                     computePosteriorMode = TRUE,
                                     kappaMax = 100) {

  # Create a likelihood function, unnormalized for numerical stability.
  VMKDELogLik    <- createKDEkpLogLik(th, normalize = TRUE)

  # The posterior function, still off by a factor (2 pi)^n, due to lack of
  # normalization.
  VMKDEPost   <- function(kp) exp(log(priorFunction(kp)) + VMKDELogLik(kp))

  # Marginal likelihood for the null hypothesis of isotropy.
  isoLogMargLik <- computeLogIsotropyMarginalLikelihood(th)

  # Marginal likelhood of the kernel density alternative.
  postLogMargLik <- log(integrate(VMKDEPost, 0, upper = kappaMax)$value)

  # Compute the log Bayes Factor
  logBF <- postLogMargLik - isoLogMargLik
  BF    <- exp(logBF)

  # Compute Posterior mode if desired.
  if (computePosteriorMode) {
    kpMode <- optimize(VMKDEPost,
                       lower = 0, upper = kappaMax, maximum = TRUE)$maximum
  } else {
    kpMode <- NULL
  }

  # Return the results/
  tibble::lst(BF, logBF,
              pHa                     = BF / (1 + BF),
              pH0                     = 1 / (1 + BF),
              bfChoseHa               = BF > 1,
              kpMode,
              HaLogMargLik = postLogMargLik,
              H0LogMargLik = isoLogMargLik)
}







