### This script details the procedure of kernel function selection for weighting policies over time after adoption 
### This intends to imitate the effect of policy "positivity" delay before or decay after reaching the maximum intended impact (standardized to 1)

kernel.matix.fun <- function(window.size, n.years, kernel.type, plot.kernel = TRUE, ...) {
  message("\nNOTE:\nThis function computes and plots all available kernels, all which will be returned in a list.\nHowever, the output will only include the full matrix of weights based on the chosen kernel type and its parameters.\n")
  set.seed(35)
  par(mfrow = c(3,2))#, mar = c(1.5, 0.5, 4, 4))
  # window.size <- window.size   # Size of the search window
  # n.years <- n.years           # Number of years in the data
  
  # Gamma distribution properties and parameters 
  # shape (k) * scale (theta) = mean
  # (k - 1) * theta = mode, when k >= 1
  k <- 2
  theta <- 2
  gamma.kernel <- dgamma(x = 0:(window.size - 1), shape = k, scale = theta)
  gamma.kernel.std <- gamma.kernel/max(gamma.kernel)
  plot(x = 0:(window.size - 1), y = gamma.kernel.std, xaxp = c(0, window.size, window.size), xlab = "Time since policy introduction", ylab = "Weight", pch = 16, , ylim = c(0,1))
  abline(v = ifelse(k >= 1, (k - 1) * theta, 0), lty = 2)
  text(x = length(gamma.kernel.std) - 4, y = max(gamma.kernel.std), pos = 1, label = paste0("Gamma\nk = ", k, ", \u03B8", " = ", theta), font = c(4, 2))
  
  # Poisson distribution with lambda parameter
  lambda <- 2
  pois.kernel <- dpois(x = 0:(window.size - 1), lambda = lambda)
  pois.kernel.std <- pois.kernel/max(pois.kernel)
  plot(x = 0:(window.size - 1), y = pois.kernel.std, xaxp = c(0, window.size, window.size), xlab = "Time since policy introduction", ylab = "Weight", pch = 16, , ylim = c(0,1))
  abline(v = which.max(pois.kernel.std) - 1, lty = 2)
  text(x = length(pois.kernel.std) - 4, y = max(pois.kernel.std), pos = 1, label = paste0("Poisson\n\u03BB", " = ", lambda), font = c(4, 2))
  
  
  # Negative binomial distribution with r and p parameters
  r <- 3    # Number of successful trials (may be non-integer)
  p <- 0.4  # Probability between 0 and 1
  nb.kernel <- dnbinom(x = 0:(window.size - 1), size = r, prob = p)
  nb.kernel.std <- nb.kernel/max(nb.kernel)
  plot(x = 0:(window.size - 1), y = nb.kernel.std, xaxp = c(0, window.size, window.size), xlab = "Time since policy introduction", ylab = "Weight", pch = 16, ylim = c(0,1))
  abline(v = which.max(nb.kernel.std) - 1, lty = 2)
  text(x = length(nb.kernel.std) - 6, y = max(nb.kernel.std), pos = 1, label = paste0("Negative Binomial\nr = ", r, ", p = ", p), font = c(4, 2))
  
  # Cauchy distribution
  # x_0 location parameter is the mode. Mean is undefined
  x <- 3            # Location parameter
  gamma.par <- 15   # Scale parameter
  cauchy.kernel <- dcauchy(x = 0:(window.size - 1), location = x, scale = gamma.par, log = FALSE)
  cauchy.kernel.std <- cauchy.kernel/max(cauchy.kernel)
  plot(x = 0:(window.size - 1), y = cauchy.kernel.std, xaxp = c(0, window.size, window.size), xlab = "Time since policy introduction", ylab = "Weight", pch = 16, ylim = c(0,1))
  abline(v = which.max(cauchy.kernel.std) - 1, lty = 2)
  text(x = length(cauchy.kernel.std) - 5, y = max(cauchy.kernel.std), pos = 1, label = paste0("Cauchy\nx0 = ", x, ", \u03B3 = ", gamma.par), font = c(4, 2))
  
  # Symmetric sigmoid kernel
  sigmoid.scale.par <- 10
  a <- 1
  b <- 0
  scaling.factor <- ifelse(window.size >= sigmoid.scale.par, window.size/sigmoid.scale.par, 1)
  x.par <- (((-window.size/2) + 0.5):(window.size/2) - 0.5)/scaling.factor
  sigmoid.kernel <- 1/(1 + exp(-a * (x.par - b))) 
  sigmoid.kernel.std <- sigmoid.kernel/max(abs(sigmoid.kernel))
  plot(x = 0:(window.size - 1), y = sigmoid.kernel.std, xaxp = c(0, window.size, window.size), xlab = "Time since policy introduction", ylab = "Weight", pch = 16, ylim = c(0,1))
  abline(v = which.max(sigmoid.kernel.std) - 1, lty = 2)
  text(x = length(sigmoid.kernel.std) - 5, y = max(sigmoid.kernel.std) - 0.2, pos = 1, label = "Symmetric\nsigmoid\na = 1, b = 0", font = c(4, 2))
  
  
  # Creating a matrix of kernel weights
  if (kernel.type == "nb") {
    kern.type <- "nb.kernel.std"
  } else if (kernel.type == "gamma") {
    kern.type <- "gamma.kernel.std"
  } else if (kernel.type == "pois") {
    kern.type <- "pois.kernel.std" 
  } else if (kernel.type == "cauchy") {
    kern.type <- "cauchy.kernel.std"
  } else if (kernel.type == "sigmoid") {
    kern.type <- "sigmoid.kernel.std"
  } else {
    stop("This kernel is not defined! Choose one of: [gamma, pois, nb, cauchy, sigmoid]")
  }
  
  kern.mat <- matrix(0, nrow = n.years, ncol = n.years + window.size) # Expanding the matrix beyond the actual number of years in the data to fit the kernel (truncated in the end)
  for (i in 1:(n.years)) {
      kern.mat[i, i:(i - 1 + window.size)] <- get(kern.type)
  }
  kern.mat <- kern.mat[1:n.years, 1:n.years]
  
  return(list(gamma_kernel = gamma.kernel.std, 
              poisson_kernel = pois.kernel.std, 
              negative_binomial_kernel = nb.kernel.std, 
              cauchy_kernel = cauchy.kernel.std,
              sigmoid_kernel = sigmoid.kernel.std,
              chosen_kernel_matrix = kern.mat))
}


