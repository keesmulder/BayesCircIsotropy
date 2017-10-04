
# Generate an MCMC sample from the von Mises posterior from a conjugate prior.
vmConjugateSample <- function(th, Q = 1000, fixedMu = NA, fixedKp = NA, R0 = sqrt(2), c0 = 2) {

  mus <- kps <- numeric(Q)
  mus[1] <- mu_cur <- 0
  kps[1] <- kp_cur <- 0

  eta <- length(th) + c0

  ct <- sum(cos(th))
  st <- sum(sin(th))
  rt <- sqrt(ct^2 + st^2) + R0
  mt <- atan2(st, ct)

  muFixed <- !is.na(fixedMu)
  kpFixed <- !is.na(fixedKp)
  if(muFixed) mus[1] <- mu_cur <- fixedMu
  if(kpFixed) kps[1] <- kp_cur <- fixedKp

  for (i in 2:Q) {
    if(!muFixed) mu_cur <- rvmc(1, mt, rt*kp_cur)

    if (!kpFixed) {
      etag   <- -rt * cos(mu_cur - mt)
      kp_cur <- sampleKappa(etag, eta)[1]
    }

    mus[i] <- mu_cur
    kps[i] <- kp_cur
  }

  cbind(mu = mus, kp = kps)
}


# MH
generateNextKappaWithMHStep <- function(mu_cur, kp_cur, logPostKernel, kappaProposalDistribution) {

  # Draw new kappa candidate
  kp_can <- kappaProposalDistribution$samplingFun(1, kp_cur)

  # Due to numerical issues, the rounded kappa candidate might be zero, which
  # causes issues in the MH-ratio. This is fixed by returning the current value.
  # This step is only reasonable if p(kappa = 0) = 0 under the prior (and
  # therefore, the posterior). This is reasonable in the current applications,
  # but this does mean that this can not always be applied. If this condidition
  # is not met, we might retry any zero proposal.
  if (kp_can == 0) return(kp_cur)

  log_prop_cur_to_can <- kappaProposalDistribution$logpdf(kp_can, kp_cur)
  log_prop_can_to_cur <- kappaProposalDistribution$logpdf(kp_cur, kp_can)

  log_prop_cur <- logPostKernel(c(mu_cur, kp_cur))
  log_prop_can <- logPostKernel(c(mu_cur, kp_can))

  log_MH_ratio <- log_prop_can + log_prop_can_to_cur -
    log_prop_cur - log_prop_can_to_cur

  if (log_MH_ratio > log(runif(1))) {
    return(kp_can)
  } else {
    return(kp_cur)
  }
}


# Generate an MCMC sample from the von Mises posterior assuming the distribution
# of the mean direction is still von Mises. Then the non-conjugate prior is on
# kappa, while the prior on mu is either circular uniform or conjugate. This is
# useful for the jeffreys posterior.
# The posterior resultant lenght Rn is required because it is used to sample mu easily.
vmMuDistrVMSample <- function(th,
                              logPostKernel, kappaProposalDistribution,
                              Q = 1000,
                              fixedMu = NA, fixedKp = NA,
                              mu_start = 0, kp_start = 1) {

  mus <- kps <- numeric(Q)
  mus[1] <- mu_cur <- mu_start
  kps[1] <- kp_cur <- kp_start

  mun <- meanDir(th)
  Rn  <- getR(th)

  muFixed <- !is.na(fixedMu)
  kpFixed <- !is.na(fixedKp)
  if (muFixed) mus[1] <- mu_cur <- fixedMu
  if (kpFixed) kps[1] <- kp_cur <- fixedKp
  for (i in 1:Q) {
    if (!muFixed) mu_cur <- rvmc(1, mun, Rn * kp_cur)

    if (!kpFixed) {
      kp_cur <- generateNextKappaWithMHStep(mu_cur, kp_cur,
                                            logPostKernel, kappaProposalDistribution)
    }

    mus[i] <- mu_cur
    kps[i] <- kp_cur
  }

  cbind(mu = mus, kp = kps)
}





# General rejection sampler, used here to obtain samples from the prior.
rejectionSampler <- function(targetFunction, n = 1, xMin = 0, xMax = 15) {

  optimum <- optimize(targetFunction, c(xMin, xMax), maximum = TRUE)

  maxHeight <- optimum$objective
  xLength   <- xMax - xMin

  sample <- numeric(n)

  for (i in 1:n) {

    done <- FALSE

    while (!done) {

      # The candidate s uniformly drawn between bounds.
      candidate <- runif(1, xMin, xMax)


      # Accept or reject the candidate.
      if(runif(1, 0, maxHeight) < targetFunction(candidate)) {
        sample[i] <- candidate
        done <- TRUE
      }
    }
  }
  sample
}
