#### CP 2019/06/26
### This function computes iaaft and ebisuzaki surrogate for time series. It is a copy-paste of kahaaga/tstools. Hower, tstools package requires the pacman library which can only be installed on version >= 3.5.0 of R, which I currently don't have. 

iaaft_surrogate <- function(series, n.max.iter = 150) {
#  if (!is_valid_input(series)) {
#    rlang::abort("Surrogate generation failed. Input series is not valid!")
#  }

  n <- length(series)

  # Fourier transform of the original series
  original_fft <- stats::fft(series)

  # Amplitudes of the original Fourier transform
  original_fft_amplitudes <- Mod(original_fft)

  # Sample a Gaussian vector and sort it
  gaussian <- sort(stats::rnorm(n, 0, 2), index.return = T)$ix

  # Sort the original series according to the sorted Gaussian. After a number
  # of iterations, this will be our final surrogate series.
  series.randsorted <- series[gaussian] ###CP. If I've understood the algorithm well enough, we should have series=series[gaussian]. If we don't, we forget this distribution after a few steps, as we change series.randsorted to surrogate AND use the sorting of series instead of the gaussian sorting for surrogate

  # Iterate until autocorrelation functions match sufficiently.
  iteration <- 0
  convergence.achieved <- FALSE

  tolerance <- 0.00001

  # Compute root mean square difference between original and randomly sorted
  # series
  acf.diff.old <- sqrt(sum((stats::acf(series, n - 1, plot = F)$acf -
                  stats::acf(series.randsorted, n - 1, plot = F)$acf) ^ 2))


  while (!convergence.achieved && iteration <= n.max.iter) {
    series.randsorted.fft <- stats::fft(series.randsorted)

    # Extract the phases from the Fourier transform of the randomly
    # sorted series.
    randomised.phases <- Arg(series.randsorted.fft)

    # New spectrum preserving original amplitudes but using the randomised
    # phases.
    new_fft <- original_fft_amplitudes * exp(randomised.phases * 1i)

    # fft with inverse = T returns unnormalised values, so normalise by n
    surrogate <- Re(stats::fft(new_fft, inverse = T)) / n # Sample the real part

    surrogate[sort(surrogate, index.return = T)$ix] <- sort(series)

    series.randsorted <- surrogate

    # Check for convergence
    acf.diff.new <- sqrt(sum((stats::acf(series, n - 1, plot = F)$acf -
                    stats::acf(surrogate, n - 1, plot = F)$acf)^2))
	#print(acf.diff.new)
	#print("old")
	#print(acf.diff.old)

    if (abs(acf.diff.old - acf.diff.new) < tolerance) {
        #cat("\nConvergence achieved after ", iteration, " iterations.\n")
        convergence.achieved <- T
    } else {
        acf.diff.old <- acf.diff.new
    }
    iteration <- iteration + 1
  }
  return(surrogate)
}

ebisuzaki_surrogate <- function(series) {

	require(rEDM)

  surrogate <- make_surrogate_data(ts = series,
                            method = "ebisuzaki",
                            num_surr = 1)
  return(as.numeric(surrogate))
}
