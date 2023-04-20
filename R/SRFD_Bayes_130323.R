#train_data <- train.data
#test_data <- test.data

SRFD_Bayes <- function(train_data, test_data, method="SVM") {

  if (param@dataset_dir=="simulation_dataset") {
      
    train_labels <- train_data[nrow(train_data),]
    sorted_train_index <- order(train_labels)
    train_data <- train_data[1:(nrow(train_data)-1),sorted_train_index]
    
    test_labels <- test_data[nrow(test_data),]
    sorted_test_index <- order(test_labels)
    test_data <- test_data[1:(nrow(test_data)-1),sorted_test_index]
    
    train_class_index <- unique(train_labels)
    test_class_index <- unique(test_labels)
    
    param@class_num <<- length(train_class_index)
    param@train_sample_num <<- unlist(lapply(1:length(train_class_index), function(i) sum(train_labels == train_class_index[i])))
    param@test_sample_num <<- unlist(lapply(1:length(test_class_index), function(i) sum(test_labels == test_class_index[i])))
    #param@class_num <- length(train_class_index)
    #param@train_sample_num <- unlist(lapply(1:length(train_class_index), function(i) sum(train_labels == train_class_index[i])))
    #param@test_sample_num <- unlist(lapply(1:length(test_class_index), function(i) sum(test_labels == test_class_index[i])))

  } else if(param@dataset_name == "Xu_data") {
    param@class_num <<- 2
    param@train_sample_num <<- c(418,704)
    param@test_sample_num <<- c(417,346)
  } else if(param@dataset_name == 'Chen_data') {
    param@class_num <<- 6
    param@train_sample_num <<- c(207,5,20,39,45,35)
    param@test_sample_num <<- c(207,2,3,29,11,34)
  }
    
  param@train_gt_label <<- unlist(lapply(1:length(param@train_sample_num), function(i) rep(i, param@train_sample_num[i])))
  param@test_gt_label <<- unlist(lapply(1:length(param@test_sample_num), function(i) rep(i, param@test_sample_num[i])))
  #param@train_gt_label <- unlist(lapply(1:length(param@train_sample_num), function(i) rep(i, param@train_sample_num[i])))
  #param@test_gt_label <- unlist(lapply(1:length(param@test_sample_num), function(i) rep(i, param@test_sample_num[i])))
  
  # Semi-reference-free deconvolution (SRFD)
  cat("\n1. Performing semi-reference-free deconvolution (SRFD)\n")
  x <- SRFD(train_data)

#  if(param@oversample){ # not yet tested
#    os_R <- BorderlineSMOTE(train_R, param$train_gt_label, 40, fraction = param)
#    reconstruction_err <- train_data - train_W %*% train_R
#    os_data <- MixWH(train_W, os_R, os_labels, reconstruction_err, param)
#    train_data <- cbind(train_data, os_data)
#    param$train_gt_label <- c(param$train_gt_label, os_labels)
#    param$train_gt_label <- sort(param$train_gt_label)
#    train_data <- train_data[,order(param$train_gt_label)]
#  }
  # Compute source fraction (deconcolution without structural constraints)
  cat("\n2. Estimating source fractions\n")
  train_sf <- deconvolution_with_reference(train_data, x$W, rep(0, ncol(train_data)), n.tries=3)
  test_sf <- deconvolution_with_reference(test_data, x$W, rep(0, ncol(test_data)), n.tries=3)
  
  # Bayesian diagnosis
  cat("\n3. Performing Bayesian inference\n")
  bayes_pred_prob_all <- list()
  for (j in 1:10) {
    cat(".")
    bd <- Bayes_diagnosis(x = train_data, y = param@train_gt_label, newx = test_data, train_sf, test_sf, method=method, param=param)
    bayes_pred_prob_all[[j]] <- bd$posterior
  }
  mean_bayes_pre_prob <- Reduce("+", bayes_pred_prob_all) / length(bayes_pred_prob_all)
  SRFD_Bayes_prob <- apply(mean_bayes_pre_prob, 1, max)
  SRFD_Bayes_prediction <- apply(mean_bayes_pre_prob, 1, function(x) which.max(x))
  
  # For cancer detection
  SRFD_Bayes_prob[which(SRFD_Bayes_prediction==1)] <- 1-SRFD_Bayes_prob[which(SRFD_Bayes_prediction==1)]
  
  # Final results
  tumor_fraction <- colSums(test_sf[(param@healthy_pattern_num+1):nrow(test_sf),])
  rm(train_data, test_data, train_sf, test_sf, train_labels, test_labels, train_class_index, test_class_index, sorted_train_index, sorted_test_index, bayes_pred_prob_all, bd)
  gc()
  return(list(train_reference = x$W, SRFD_Bayes_prob = SRFD_Bayes_prob, SRFD_Bayes_prediction = SRFD_Bayes_prediction, tumor_fraction = tumor_fraction))
}
