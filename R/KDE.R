# In this file, the following terms are used:
#  x:  a new angle to be evaluated.
#  th: the existing dataset.
#  kp: bandwith / kappa.


# Unnormalized log pdf of the von Mises distribution
logVM <- function(th, mu, kp) kp * cos(th - mu)

# Evaluate kernel density estimate given parameters.
evalVMKDE <- function(x, th, kp) {
  sum(exp(logVM(x, th, kp))) / length(th)
}

# Create a pdf function for a specific dataset and bandwith.
createVMKDEpdf <- function(th, kp) {
  Vectorize(function(x) evalVMKDE(x, th, kp))
}

# Create a likelihood for the bandwidth with LOO-CV.
# If the likelihood is not normalized, it is off by a factor (2*pi)^n.
createKDEkpLik <- function(th, normalize = TRUE) {
  if (normalize) { nc <- 2 * pi } else {nc <- 1}

  Vectorize(
    function(kp) {
      (nc * besselI(kp, 0)) ^ (-length(th)) *
        prod(sapply(1:length(th),
                    function(i) evalVMKDE(th[i], th[-i], kp)))
    }
  )
}


# Create a log-likelihood for the bandwidth with LOO-CV.
# If the likelihood is not normalized, it is off by a factor (2*pi)^n.
createKDEkpLogLik <- function(th, normalize = TRUE) {
  if (normalize) { lnc <- log(2 * pi) } else {lnc <- 0 }

  Vectorize(
    function(kp) {
      -length(th) * (lnc + logBesselI(kp, 0)) +
        sum(sapply(1:length(th),
                   function(i) log(evalVMKDE(th[i], th[-i], kp))))
    }
  )
}






