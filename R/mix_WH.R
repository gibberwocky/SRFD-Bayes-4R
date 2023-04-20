#train_W <- x$W
#extend_H <- os_R
#extend_labels <- os_labels

mix_WH <- function(train_W, extend_H, extend_labels, reconstruction_err, param) {
  library(MASS)
  extend_WH <- train_W %*% extend_H
  extend_id <- unique(extend_labels)
  extend_sample_num <- sapply(extend_id, function(x) sum(extend_labels == x))
  extend_sample_num <- c(0, extend_sample_num)
  train_sample_num <- c(0, param@train_sample_num)
  extend_data <- numeric(0)
  for (i in 1:length(extend_id)){
    extend_class_WH <- extend_WH[,sum(extend_sample_num[1:i])+1:sum(extend_sample_num[1:i+1])]
    rec_class_err <- reconstruction_err[,sum(train_sample_num[1:i])+1:sum(train_sample_num[1:i+1])]
    rec_class_err <- as.vector(rec_class_err)
    fit <- fitdistr(rec_class_err, "normal")
    mu <- fit$estimate[1]
    sigma <- fit$estimate[2]
    extend_err <- matrix(nrow = nrow(extend_class_WH), ncol = ncol(extend_class_WH))
    for (j in 1:ncol(extend_class_WH)) {
      extend_err[,j] <- rnorm(nrow(extend_class_WH),mu,sigma)
    }
    extend_WH_err <- extend_class_WH + extend_err
    extend_WH_err[extend_WH_err > 1] <- 1
    extend_WH_err[extend_WH_err < 0] <- 0
    extend_data <- cbind(extend_data, extend_WH_err)
  }
  return(extend_data)
}
