#library(pracma)

#train_W <- x$W
#predicted_label <- rep(0, ncol(train_data))

deconvolution_with_reference <- function(test_data, train_W, predicted_label, n.tries=3) {
  
  set.seed(param@seed)
  R <- param@healthy_pattern_num + (param@class_num-1)*param@cancer_pattern_num
  
  test_H <- lapply(1:ncol(test_data), function(i)
    {
      if(predicted_label[i] == 0) { # no structural constraints
        mask <- rep(0,R)
      } else if(k == 1) { # predicted healthy controls
        mask <- c(rep(0, param@healthy_pattern_num), rep(1, param@cancer_pattern_num*(param@class_num-1)))
      } else { # predicted tumour types
        mask <- c(rep(0, param@healthy_pattern_num), rep(1, param@cancer_pattern_num*(param@class_num-1)))
        mask[(param@healthy_pattern_num + param@cancer_pattern_num*(predicted_label[i]-2)+1):(param@healthy_pattern_num + param@cancer_pattern_num*(predicted_label[i]-1))] <- 0
      }
    
      cat(".")
      if((i/80)%%1==0) cat("[", i, "]\n")
    
      x.fmincon <- NULL  
      attempts <- 0
      fun <- function(y) { norm(test_data[,i]-train_W %*% y, type = "F") }
      while(is.null(x.fmincon) && attempts < n.tries) {
        attempts <- attempts + 1
        tryCatch( {
          x.fmincon <- pracma::fmincon(x0 = runif(R, min = 0, max = 1), fn = fun, A = NULL, b = NULL, Aeq = rbind(rep(1, R), mask), beq = c(1, 0), lb = rep(0, R), ub = rep(1, R))$par 
          }, error = function(cond) {
            message(paste("\nError acting on data[,", i, "], attempt ", attempts, " of ", n.tries, sep=""))
            message(cond, "\n")
            return(NULL)
          } )
      }
      return(x.fmincon)
  })
  tryCatch( { test_H <- do.call("rbind", test_H) }, error = function(cond) { message("Failed to generate matrix, try increaseing n.tries"); return(NULL) })
  return(t(test_H))
}



