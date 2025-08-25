apply.positivity.weights.fun <- function(data, policy.cols, kernel.matrix) {
  require(tidyverse)
  split.df <- split(data, list(data$state_name, data$policy_category, data$policy_class))
  out <- lapply(split.df, function(x) {
      x1 <- x[, policy.cols] 
      x2 <- t(as.matrix(x1))
      
      x3 <- matrix(0, nrow = nrow(x2), ncol = ncol(x2))
      for (i in 1:nrow(x2)) {
        # Finding the first instance of the policy in effect (assuming the data is chronologically sorted, this is the year of policy adoption)
        kernel.mat.index <- min(which(!is.na(x2[i,])))
        
        if (!is.na(kernel.mat.index)) {
          x3[i,] <- x2[i,] * kernel.mat[kernel.mat.index,]
        }
      }
      colnames(x3) <- colnames(x2)
      rownames(x3) <- rownames(x2)
      x3 <- as.tibble(t(x3))
  })
 
  return(bind_cols(data[, 1:4], bind_rows(out)))
}


