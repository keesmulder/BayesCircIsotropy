

logBesselI <- function(x, nu) log(besselI(x, nu, expon.scaled = TRUE)) + x

A    <- function(kp) besselI(kp, 1) / besselI(kp, 0)
logA <- function(kp) logBesselI(kp, 1) - logBesselI(kp, 0)

# Evaluate the (possibly truncated) Jeffreys prior of the von Mises distribtuion.
vmJeffreysPrior <- function(kp, n0 = 1, log = FALSE, truncateAt = 0) {
  Akp <- A(kp)
  if (log) {
    log(kp >= truncateAt) +
      0.5 * (2 * log(n0) + log(kp) + log(Akp) + log(0.5 + besselI(kp, 2)/(2 * besselI(kp, 0)) - Akp ^ 2))
  } else {
    (kp >= truncateAt) *
      sqrt(n0 ^ 2 * kp * Akp * (0.5 + besselI(kp, 2)/(2 * besselI(kp, 0)) - Akp ^ 2))
  }
}



LogVMJeffreysPrior <- function(kp) {
  logAkp <- logA(kp)
  0.5 * (log(kp) +
           logAkp +
           log(0.5 + exp(logBesselI(kp, 2) -
                           log(2) -
                           logBesselI(kp, 0)) -
                           exp(logAkp) ^ 2)
         )
}

# Create a function that returns the (possibly truncated) Jeffreys prior of the
# von Mises distribution. inputLength dictates the length of the input vector,
# which might be either c(mu, kp) or simply kp.
createVMJeffreysPrior <- function(truncLower = 0, truncUpper = 20,
                                  normalize = TRUE, inputLength = 2) {
  unnormalizedFun <- function(kp) {
    logAkp <- logBesselI(kp, 1) - logBesselI(kp, 0)
    return(
      (kp >= truncLower) *
        (kp <= truncUpper) *
      exp(0.5 * (log(kp) +
               logAkp +
               log(0.5 + exp(logBesselI(kp, 2) -
                               log(2) -
                               logBesselI(kp, 0)) -
                     exp(logAkp) ^ 2)))
      )
  }

  if (normalize) {
    nc <- integrate(unnormalizedFun, truncLower, truncUpper)$value
  } else {
    nc <- 1
  }

  if (inputLength == 1) {
    return(function(kp) unnormalizedFun(kp) / (nc * 2 * pi))
  } else if (inputLength == 2) {
    return(function(x) unnormalizedFun(x[2]) / (nc * 2 * pi))
  } else {stop("Impossible number of arguments required.")}
}



# Create a function that returns  the (possibly truncated) log Jeffreys prior of
# the von Mises distribution.
createVMJeffreysLogPrior <- function(truncLower = 0, truncUpper = Inf,
                                     normalize = TRUE) {
  unnormalizedFun <- function(kp) {
    logAkp <- logBesselI(kp, 1) - logBesselI(kp, 0)
    return(
      (kp >= truncLower) *
        (kp <= truncUpper) *
        exp(0.5 * (log(kp) +
                     logAkp +
                     log(0.5 + exp(logBesselI(kp, 2) -
                                     log(2) -
                                     logBesselI(kp, 0)) -
                           exp(logAkp) ^ 2)))
    )
  }

  if (normalize) {
    nc <- integrate(unnormalizedFun, truncLower, truncUpper)$value
  } else {
    nc <- 1
  }

  # Return the normalized log function
  function(kp) {
    logAkp <- logBesselI(kp, 1) - logBesselI(kp, 0)
    0.5 * (log(kp) +
             logAkp +
             log(0.5 + exp(logBesselI(kp, 2) -
                             log(2) -
                             logBesselI(kp, 0)) -
                   exp(logAkp) ^ 2)) - log(nc * 2 * pi) +
      log(kp >= truncLower) + log(kp <= truncUpper)
  }
}





# n0 is the size of the dataset for which the Fisher's information is used.
# Usually, this is 1. It does not matter too much, as it is just a
# multiplicative constant, although it still determines the area under the curve
# of the posterior, which means this value also factors in the Bayes Factor if
# the improper prior is mistakenly used.
createVMJeffreysPostKernel <- function(th, n0 = 1, log = FALSE) {

  # Obtain the likelihood part of the posterior.
  vmLikelihood <- createVMLikelihoodFunction(th, log = log)

  # Get the posterior kernel
  if (log) {
    vmJeffPostKernel <- function(x) vmJeffreysPrior(x[2], n0 = n0, log = TRUE) + vmLikelihood(x)
  } else {
    vmJeffPostKernel <- function(x) vmJeffreysPrior(x[2], n0 = n0, log = FALSE) * vmLikelihood(x)
  }
  vmJeffPostKernel
}


# Previously not the true kernel, as the prior did not integrate to 1. Setting normalizePrior to FALSE here just fudges things up.
vmConjKappaOnlyPostKernel <- function(th, R0 = sqrt(2), c0 = 2, log = FALSE, normalizePrior = TRUE) {

  # Obtain the likelihood part of the posterior.
  vmLikelihood         <- createVMLikelihoodFunction(th, log = log)
  vmConjKappaOnlyPrior <- createVMKappaOnlyConjPrior(R0 = R0, c0 = c0, normalize = normalizePrior)

  # Get the posterior kernel
  if (log) {
    vmConjPostKernel <- function(x) vmConjKappaOnlyPrior(x) + vmLikelihood(x)
  } else {
    vmConjPostKernel <- function(x) vmConjKappaOnlyPrior(x) * vmLikelihood(x)
  }
  vmConjPostKernel
}


# General function that takes data, a function to generate a likelihood from the
# data, and a prior to return a posterior function.
createLogPosteriorFunction <- function(data, likelihoodGenFunc, priorFunction, log = FALSE) {

  # Get the correct likelihood or loglikelihood
  hasLogOption <- "log" %in% names(formals(likelihoodGenFunc))

  if (hasLogOption) {
    likelihoodFunction <- likelihoodGenFunc(data, log = TRUE)

    posteriorFunction <- function(parameter) {
      priorFunction(parameter, log = TRUE) +  likelihoodFunction(parameter, log = TRUE)
    }
  }  else {
    likelihoodFunction <- likelihoodGenFunc(data)

    posteriorFunction <- function(parameter) {
      priorFunction(parameter, log = TRUE) +  log(likelihoodFunction(parameter, log = TRUE))
    }
  }

  posteriorFunction

}


# General function that takes data, a function to generate a likelihood from the
# data, and a prior to return a posterior function.
createPosteriorFunction <- function(data, likelihoodGenFunc, priorFunction) {

  likelihoodFunction <- likelihoodGenFunc(data)

  function(parameter) priorFunction(parameter) * likelihoodFunction(parameter)
}


# Create my own kappa prior function that corresponds with two angles at 90 degrees.
createVMKappaOnlyConjPrior <- function(R0 = sqrt(2), c0 = 2, kappaMax = 100,
                                       normalize = TRUE, truncateAt = 0) {

  if (truncateAt == 0) {
    unnormalizedKpFun <- function(x) besselI(R0 * x, 0) / (besselI(x, 0) ^ c0)
  } else {
    unnormalizedKpFun <- function(x) besselI(R0 * x, 0) / (besselI(x, 0) ^ c0) * (x > truncateAt)
  }


  if (normalize) {
    nc <- integrate(unnormalizedKpFun, 0, kappaMax)$value
    return(function(x) Vectorize(unnormalizedKpFun(x[2])/(nc*2*pi)))
  } else {
    return(function(x) Vectorize(unnormalizedKpFun(x[2])))
  }
}


# Create my own kappa conditional posterior, where we take mu_cur = mu_n so that the cosine term drops
createKappaCondPosterior <- function(Rn = sqrt(2), eta = 2, kappaMax = 100, normalize = TRUE) {
  unnormalizedKpFun <- function(x) exp(Rn * x) / (besselI(x, 0)^eta)

  if (normalize) {
    nc <- integrate(unnormalizedKpFun, 0, kappaMax)$value
    return(function(x) Vectorize(unnormalizedKpFun(x)/nc))
  } else {
    return(function(x) Vectorize(unnormalizedKpFun(x)))
  }
}



# Create a posterior function for the parameters of a von Mises distribution with conjugate prior.
vmConjugatePostKernel <- function(th, R0 = 0, c0 = 0, mu0 = meanDir(th), log = FALSE) {
  n <- length(th)

  # Sufficient tatistics for the von Mises posterior
  Cn  <- R0 * cos(mu0) + sum(cos(th))
  Sn  <- R0 * sin(mu0) + sum(sin(th))
  Rn  <- sqrt(Cn ^ 2 + Sn ^ 2)
  mun <- atan2(Sn, Cn)

  # The kernel of the von Mises log-posterior
  logFun <- function(x) {
    -(n + c0) * log(2 * pi) -
      (n + c0) * logBesselI(nu = 0, x = x[2]) +
      Rn*x[2]*cos(mun - x[1])
  }

  # Construct outputfunction.
  if (log) {
    vmPostKernel <- function(x) {
      if (x[2] < 0) {
        return(0)
      } else {
        return(logFun(x))
      }
    }
  } else {
    vmPostKernel <- function(x) {
      if (x[2] < 0) {
        return(0)
      } else {
        return(exp(logFun(x)))
      }
    }
  }
  vmPostKernel
}


# Create a von Mises likelihood function.
createVMLikelihoodFunction <- function(th, log = FALSE) {
  vmConjugatePostKernel(th, R0 = 0, c0 = 0, log = log)
}
