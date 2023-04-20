calculate_2p_norm <- function(X, p) {
  x_norm <- sum(apply(X^2, 1, sum)^(p/2))
  return(x_norm)
}

SRFD <- function(X){
  
  K <- nrow(X)
  N <- ncol(X)
  C <- param@healthy_pattern_num + (param@class_num - 1) * param@cancer_pattern_num
  
  # Initialization
  W <- matrix(runif(K * C), K, C)
  R <- matrix(runif(C * N), C, N)
  R <- R / rep(colSums(R), each = nrow(R))

  # Structural mask
  healthy_top_mask <- matrix(0, nrow = param@healthy_pattern_num, ncol = N)
  healthy_left_mask <- matrix(1, nrow = C - param@healthy_pattern_num, ncol = param@train_sample_num[1])
  cancer_mask <- pracma::blkdiag(matrix(1, nrow = param@cancer_pattern_num, ncol = param@train_sample_num[2]))
  if(length(param@train_sample_num)>2) {
    for (i in 2:(length(param@train_sample_num)-1)) { 
        cancer_mask <- pracma::blkdiag(cancer_mask, matrix(1, nrow = param@cancer_pattern_num, ncol = param@train_sample_num[i+1]))
    }
  }
  cancer_mask <- 1 - cancer_mask
  train_mask <- cbind(healthy_left_mask, cancer_mask)
  mask <- rbind(healthy_top_mask, train_mask)
  # clean up
  rm(healthy_top_mask, healthy_left_mask, cancer_mask, train_mask)
  gc()
  
  # objective function
  Omega <- 1/2*param@eta*norm(mask*R, "F")^2
  # calculate_2p_norm(X-W*H,p) 
  err <- NULL
  err[1] <- 1/2*calculate_2p_norm(X - W %*% R, param@p) + Omega
  
  iter <- 1
  while (iter < param@maximum_iterations) {
    cat(".")
    if((iter/80)%%1==0) cat("[", iter, "]\n")
    
    Z <- X - W %*% R
    L2_p <- param@p / (2 * (colSums(Z^2)) ^ (1 - param@p/2)) 
    D <- diag(L2_p)
    
    # update W
    factorW_numerator <- X %*% D %*% t(R) 
    factorW_denominator <- W %*% (R %*% D %*% t(R)) 
    W <- W * (factorW_numerator / factorW_denominator)
    W[W > 1] <- 1.0 - 1e-4
    
    # update H
    factorH_numerator <- t(W) %*% X %*% D
    factorH_denominator <- t(W) %*% W %*% R %*% D + param@eta * mask * R
    R <- R * (factorH_numerator / factorH_denominator)
    R <- R / rep(colSums(R), each = nrow(R)) # normalization

    iter <- iter + 1
    Omega <- 1/2 * param@eta * norm(mask * R, "F")^2
    err[iter] <- 1/2 * calculate_2p_norm(X - W %*% R, param@p) + Omega
    
    if (abs(err[iter-1] - err[iter]) < param@convergence_threshold) {
      break
    }
  }
  rm(Z, D, L2_p, X, K, N, C, Omega, factorW_numerator, factorW_denominator, factorH_numerator, factorH_denominator, mask)
  gc()
  return(list(W=W, R=R, err=err))
}  