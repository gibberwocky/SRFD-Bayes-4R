#diagnosis_results <- SRFD_Bayes_results
#predict_label <- diagnosis_results$SRFD_Bayes_prediction
#num_in_class <- param@test_sample_num

evaluate_diagnosis <- function(diagnosis_results){
  cat("\nconfusion_matrix: \n")
  confusion_matrix <- compute_confusion_matrix(diagnosis_results$SRFD_Bayes_prediction, param@test_sample_num)
  print(confusion_matrix)
  ACC <- diag(confusion_matrix)/param@test_sample_num
  return(mean(ACC))
}

## compute confusion matrix
compute_confusion_matrix <- function(predict_label,num_in_class) {
  num_class <- length(num_in_class)
  num_in_class <- c(0,num_in_class)
  confusion_matrix <- matrix(nrow=num_class,ncol=num_class)
  for (ci in 1:num_class) {
    for (cj in 1:num_class) {
      summer <- 0
      c_start <- sum(num_in_class[1:ci]) + 1
      c_end <- sum(num_in_class[1:ci+1])
      summer <- length(which(predict_label[c_start:c_end]==cj))
      confusion_matrix[ci, cj] <- summer
    }
  }
  return(confusion_matrix)
}

