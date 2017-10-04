# Generate von Mises data with a Dirichlet process on the means.
generateVMDPData <- function(n, kp = 5) {
  th <- mus <- numeric(n)

  for (i in 1:n) {

    # The component to belong to.
    choice <- sample(1:i, 1)

    # Create a new component with probability 1/i
    if (choice == i) {
      mus[i] <- runif(1, 0, 2*pi)

      # Otherwise, use the as the mean a previously used mean.
    } else {
      mus[i] <- mus[choice]
    }
    th[i] <- rvmc(1, mus[i], kp)
  }
  th %% (2*pi)
}


# Generate Data with von Mises distribution centered with equal probability on
# each previously sampled data point or a new one.
generateVMPseudoDPData <- function(n, kp = 5) {
  th <- numeric(n)
  for (i in 1:n) {

    # The component to belong to.
    choice <- sample(1:i, 1)

    # Create a new component with probability 1/i
    if (choice == i) {
      th[i] <- runif(1, 0, 2*pi)

    # Otherwise,, use the as the mean a previously drawn value.
    } else {
      th[i] <- rvmc(1, th[choice], kp)
    }
  }
  th %% (2*pi)
}


